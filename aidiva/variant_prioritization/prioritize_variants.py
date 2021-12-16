import argparse
import networkx as nx
import numpy as np
import os
import pandas as pd
import multiprocessing as mp
import pickle
import re
from functools import partial
from itertools import combinations
import logging
import pysam

if not __name__=="__main__":
    from . import get_HPO_similarity_score as gs


CODING_VARIANTS = ["splice_acceptor_variant",
                   "splice_donor_variant",
                   "stop_gained",
                   "frameshift_variant",
                   "stop_lost",
                   "start_lost",
                   "inframe_insertion",
                   "inframe_deletion",
                   "missense_variant",
                   "protein_altering_variant",
                   "splice_region_variant",
                   "incomplete_terminal_codon_variant",
                   "start_retained_variant",
                   "stop_retained_variant",
                   "synonymous_variant",
                   "coding_sequence_variant"]
                   #"5_prime_UTR_variant",
                   #"3_prime_UTR_variant"]


logger = logging.getLogger(__name__)


cadd_identifier = "CADD_PHRED"
duplication_identifier = "segmentDuplication"
repeat_identifier = "simpleRepeat"


def prioritize_variants(variant_data, hpo_resources_folder, reference, num_cores, family_file=None, family_type="SINGLE", hpo_list=None, gene_exclusion_list=None):
    # load HPO resources
    gene_2_HPO_f = hpo_resources_folder + "gene2hpo_test.pkl"
    hgnc_2_gene_f = hpo_resources_folder + "hgnc2gene_test.pkl"
    gene_2_interacting_f = hpo_resources_folder + "gene2interacting_test.pkl"
    HPO_graph_file = hpo_resources_folder + "hpo_graph.gpkl"
    hpo_list_file = hpo_list
    gene_exclusion_file = gene_exclusion_list

    hgnc_2_gene = pickle.load(open(hgnc_2_gene_f, "rb"))
    gene_2_interacting = pickle.load(open(gene_2_interacting_f, "rb"))
    gene_2_HPO = pickle.load(open(gene_2_HPO_f, "rb"))

    HPO_graph = nx.read_gpickle(HPO_graph_file)
    information_content_per_node = nx.get_node_attributes(HPO_graph, "IC")
    node_ancestor_mapping = {hpo_term: nx.ancestors(HPO_graph, hpo_term) for hpo_term in HPO_graph}

    genes2exclude = set()
    if gene_exclusion_file is not None:
        if os.path.isfile(gene_exclusion_file):
            with open(gene_exclusion_file, "r") as exclusion_file:
                for line in exclusion_file:
                    if line.startswith("#"):
                        continue
                    if line == "\n":
                        continue
                    gene = line.rstrip().split("\t")[0].upper()
                    genes2exclude.add(gene)
        else:
            logger.error("The specified gene exclusion list %s is not a valid file" % (gene_exclusion_file))
            logger.warn("No genes will be excluded during filtering!")

    HPO_query = set()
    HPO_query_distances = 0
    if hpo_list_file is not None:
        if os.path.isfile(hpo_list_file):
            with open(hpo_list_file, "r") as hpo_file:
                for line in hpo_file:
                    HPO_term = line.rstrip()
                    HPO_query.add(HPO_term)
            HPO_query = list(HPO_query)
            HPO_query.sort() # makes sure that the gene symbols are ordered (could lead to problems otherwise)
        else:
            logger.error("The specified HPO list %s is not a valid file" % (hpo_list_file))
            logger.warn("HPO score finalization will be skipped!")

    # read family relationship (from PED file)
    family = dict()
    if family_file is not None:
        if os.path.isfile(family_file):
            with open(family_file, "r") as fam_file:
                for line in fam_file:
                    line = line.rstrip()
                    splitline = line.split("\t")

                    if splitline[5] == "2":
                        family[splitline[1]] = 1
                    elif splitline[5] == "1":
                        family[splitline[1]] = 0
                    else:
                        logger.error("There was a problem with the given PED file describing the family relations.")
        else:
            logger.error("The specified family file %s is not a valid file" % (family_file))
            logger.warn("Inheritance assessment will be skipped!")

    variant_data = parallelize_dataframe_processing(variant_data, partial(parallelized_variant_processing, family, family_type, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, HPO_query_distances, information_content_per_node, node_ancestor_mapping, reference), num_cores)
    # TODO: Add flag to choose the final score for sorting (if DB annotation from ClinVar/HGMD should be used)
    variant_data = variant_data.sort_values(["FINAL_AIDIVA_SCORE"], ascending=[False])
    variant_data = variant_data.reset_index(drop=True)

    return variant_data


def parallelize_dataframe_processing(variant_data, function, num_cores):
    num_partitions = num_cores * 2

    if len(variant_data) <= num_partitions:
        dataframe_splitted = np.array_split(variant_data, 1)
    else:
        dataframe_splitted = np.array_split(variant_data, num_partitions)

    try:
        pool = mp.Pool(num_cores)
        variant_data = pd.concat(pool.map(function, dataframe_splitted))
    finally:
        pool.close()
        pool.join()

    return variant_data


def parallelized_variant_processing(family, family_type, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, HPO_query_distances, ic_per_nodes, node_ancestor_mapping, reference, variant_data):
    variant_data = check_inheritance(variant_data, family_type, family)
    # TODO: implement handling of annotations from pathogenicity (and disease) databases -> ClinVar, HGMD, (OMIM)
    variant_data[["HPO_RELATEDNESS", "HPO_RELATEDNESS_INTERACTING", "FINAL_AIDIVA_SCORE"]] = variant_data.apply(lambda variant: pd.Series(compute_hpo_relatedness_and_final_score(variant, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, HPO_query_distances, ic_per_nodes, node_ancestor_mapping)), axis=1)
    variant_data[["FINAL_AIDIVA_SCORE_DB"]] = variant_data.apply(lambda variant: pd.Series(check_databases_for_pathogenicity_classification(variant)), axis=1)
    variant_data[["FILTER_PASSED", "FILTER_COMMENT"]] = variant_data.apply(lambda variant: pd.Series(check_filters(variant, genes2exclude, HPO_query, reference)), axis=1)

    return variant_data


def check_databases_for_pathogenicity_classification(variant):
    improved_score = np.nan

    clinvar_classification = str(variant["CLINVAR_DETAILS"]).split("%")[0]
    if "(replaced" in clinvar_classification:
        clinvar_classification = clinvar_classification.split("(replaced")[0]
    elif "/" in clinvar_classification:
        if clinvar_classification == "benign/likely_benign" or clinvar_classification == "likely_benign/benign":
            clinvar_classification  = "likely_benign"
        elif clinvar_classification == "pathogenic/likely_pathogenic" or clinvar_classification == "likely_pathogenic/pathogenic":
            clinvar_classification = "likely_pathogenic"
        else:
            logger.warn(f"Found unknown ClinVar classification {clinvar_classification}")

    hgmd_classification = str(variant["HGMD_CLASS"])

    if clinvar_classification == "benign":
        clinvar_score = 0.0
    elif clinvar_classification == "likely_benign":
        clinvar_score = 0.25
    elif clinvar_classification == "uncertain_significance":
        clinvar_score = 0.5
    elif clinvar_classification == "likely_pathogenic":
        clinvar_score = 0.75
    elif clinvar_classification == "pathogenic":
        clinvar_score = 1.0
    else:
        # handle missing database scores as if they have uncertain significance
        clinvar_score = 0.5
    
    if hgmd_classification == "DM":
        hgmd_score = 1.0
        #if float(variant["HGMD_RANKSCORE"]) > 0.5:
        #    hgmd_score = 1.0
    elif hgmd_classification == "DM?":
        hgmd_score = 0.75
    #elif hgmd_classification == "DP":
    #    hgmd_score = 0.5
    #elif hgmd_classification == "FP":
    #    hgmd_score = 0.5
    #elif hgmd_classification == "DFP":
    #    hgmd_score = 0.5
    #elif hgmd_classification == "R":
    #    hgmd_score = np.nan
    else:
        hgmd_score = np.nan


    # TODO: what happens if ClinVar and HGMD have contradicting classifications???



    if np.isnan(clinvar_score) and not np.isnan(hgmd_score):
        database_score = hgmd_score
    elif np.isnan(hgmd_score) and not np.isnan(clinvar_score):
        database_score = clinvar_score
    elif not np.isnan(hgmd_score) and not np.isnan(clinvar_score):
        database_score = (hgmd_score + clinvar_score) / 2
    else:
        database_score = np.nan
    
    if not np.isnan(float(variant["FINAL_AIDIVA_SCORE"])):
        if not np.isnan(database_score):
            improved_score = 0.8 * float(variant["FINAL_AIDIVA_SCORE"]) + 0.2 * database_score
        else:
            improved_score = float(variant["FINAL_AIDIVA_SCORE"])

    return improved_score


def compute_hpo_relatedness_and_final_score(variant, genes2exclude, gene_2_HPO, hgnc_2_gene, gene_2_interacting, HPO_graph, HPO_query, HPO_query_distances, ic_per_nodes, node_ancestor_mapping):
    if HPO_query:
        if np.isnan(variant["AIDIVA_SCORE"]) and ((str(variant["rf_score"]) == "nan") or (str(variant["rf_score"]) == "")) and ((str(variant["ada_score"]) == "nan") or (str(variant["ada_score"]) == "")):
            final_score = np.nan
            hpo_relatedness = np.nan
            hpo_relatedness_interacting = np.nan

        else:
            if np.isnan(variant["AIDIVA_SCORE"]):
                # if both scores are present use maximum to benefit of the strengths of each model
                if ((str(variant["rf_score"]) != "nan") and (str(variant["rf_score"]) != "")) and ((str(variant["ada_score"]) != "nan") and (str(variant["ada_score"]) != "")):
                    pathogenictiy_prediction = max([float(variant["rf_score"]), float(variant["ada_score"])], default=np.nan)
                elif ((str(variant["rf_score"]) != "") and (str(variant["rf_score"]) != "nan")) and ((str(variant["ada_score"]) == "") or (str(variant["ada_score"]) == "nan")):
                    pathogenictiy_prediction = float(variant["rf_score"])
                elif ((str(variant["ada_score"]) != "") and (str(variant["ada_score"]) != "nan")) and ((str(variant["rf_score"]) == "") or (str(variant["rf_score"]) == "nan")):
                    pathogenictiy_prediction = float(variant["ada_score"])
                else:
                    pathogenictiy_prediction = np.nan
                
                logger.info("Using dbscSNV prediction instead of AIdiva prediction for splicing variant!")

            else:
                pathogenictiy_prediction = float(variant["AIDIVA_SCORE"])

            variant_gene = str(variant["SYMBOL"]).upper()
            hgnc_id = str(variant["HGNC_ID"])

            ## TODO: check here if the gene symbol is the approved one (use hgnc id)

            gene_distances = []
            gene_distances_interacting = []

            if variant_gene in gene_2_interacting.keys():
                interacting_genes = [gene_interaction["interacting_gene"] for gene_interaction in gene_2_interacting[variant_gene]]
            else:
                interacting_genes = []

            # we use the hgnc ID to prevent problems if a given gene symbol isn't used anymore   
            # let variants with a high AIDIVA_SCORE (>=0.8) pass to prevent missed hits
            if (variant_gene not in genes2exclude) or (pathogenictiy_prediction >= 0.8):
                if variant_gene in gene_2_HPO.keys():
                    gene_HPO_list = gene_2_HPO.get(variant_gene, [])

                else:
                    if (str(hgnc_id) != "nan") and (hgnc_id in hgnc_2_gene.keys()):
                        gene_symbol = hgnc_2_gene[hgnc_id]
                        gene_HPO_list = gene_2_HPO.get(gene_symbol, [])

                    else:
                        logger.warn("The processed gene %s (HGNC:%s) is not found in the current HPO resources" % (variant_gene, hgnc_id))
                        gene_HPO_list = []

                if gene_HPO_list:
                    g_dist = gs.calculate_hpo_set_similarity(HPO_graph, HPO_query, gene_HPO_list, ic_per_nodes, node_ancestor_mapping)
                    gene_distances.append(g_dist)

                # TODO: check with HGNC ID to prevent use of previous gene symbols
                for interacting_gene in interacting_genes:
                    if (interacting_gene in genes2exclude) and (pathogenictiy_prediction < 0.8):
                        continue

                    if interacting_gene in gene_2_HPO.keys():
                        gene_HPO_list = gene_2_HPO.get(interacting_gene, [])

                    else:
                        logger.warn("The processed interacting gene %s is not found in the current HPO resources" % (interacting_gene))
                        gene_HPO_list = []

                    if gene_HPO_list:
                        g_dist = gs.calculate_hpo_set_similarity(HPO_graph, HPO_query, gene_HPO_list, ic_per_nodes, node_ancestor_mapping)
                        gene_distances_interacting.append(g_dist)
                    

                #if gene_distances or gene_distances_interacting:
                hpo_relatedness = max(gene_distances, default=0.0)
                hpo_relatedness_interacting = max(gene_distances_interacting, default=0.0)

                ## TODO: try different weighting of AIDIVA_SCORE and HPO_RELATEDNESS and HPO_RELATEDNESS_INTERACTING (eg 0.6 and 0.3 and 0.1)
                # predicted pathogenicity has a higher weight than the HPO relatedness
                final_score = (pathogenictiy_prediction * 0.6 + float(hpo_relatedness) * 0.3 + float(hpo_relatedness_interacting) * 0.1)

            else:
                final_score = np.nan
                hpo_relatedness = np.nan
                hpo_relatedness_interacting = np.nan

    else:
        if np.isnan(variant["AIDIVA_SCORE"]):
            if ((str(variant["rf_score"]) == "") or (str(variant["rf_score"]) == "nan")) and ((str(variant["ada_score"]) == "") or (str(variant["ada_score"]) == "nan")):
                pathogenictiy_prediction = np.nan

            else:
                if (str(variant["rf_score"]) != "") and (str(variant["rf_score"]) != "nan"):                        
                    rf_score = float(variant["rf_score"])

                else:
                    rf_score = np.nan
                
                if (str(variant["ada_score"]) != "") and (str(variant["ada_score"]) != "nan"):
                    ada_score = float(variant["ada_score"])

                else:
                    ada_score = np.nan

                pathogenictiy_prediction = max([rf_score, ada_score], default=np.nan)
                logger.info("Using dbscSNV prediction instead of AIdiva prediction for splicing variant!")

        else:
            pathogenictiy_prediction = float(variant["AIDIVA_SCORE"])

        final_score = pathogenictiy_prediction
        hpo_relatedness = np.nan
        hpo_relatedness_interacting = np.nan

    return [hpo_relatedness, hpo_relatedness_interacting, final_score]


def check_inheritance(variant_data, family_type="SINGLE", family=None):
    variant_columns = variant_data.columns
    if family_type == "SINGLE":
        variant_data["COMPOUND"] = 0
        variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive_single(variant, variant_columns), axis=1)

        variant_data_grouped = [group for key, group in variant_data.groupby("SYMBOL")]

        for group in variant_data_grouped:
            check_compound_single(group, variant_columns)
        variant_data = pd.concat(variant_data_grouped)

    elif family_type == "TRIO":
        if not family is None:
            variant_data["COMPOUND"] = 0
            variant_data["DOMINANT_DENOVO"] = variant_data.apply(lambda variant: check_denovo(variant, family), axis=1)

            # TODO: do we need a check for affected family members?
            variant_data["DOMINANT"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
            variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
            variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)

            variant_data_grouped = [group for key, group in variant_data.groupby("SYMBOL")]

            affected_child = ""
            parent_1 = ""
            parent_2 = ""

            # TODO: handle the case if affected status is given as unknown
            # TODO: change to be less strict and to not rely on unaffected family members
            for name in family.keys():
                if family[name] == 1:
                    affected_child = name
                elif family[name] == 0:
                    if not parent_1:
                        parent_1 = name
                        continue
                    elif not parent_2:
                        parent_2 = name
                        continue
                    else:
                        logger.error("There was a problem processing the loaded PED file")

            if  affected_child and parent_1 and parent_2:
                for group in variant_data_grouped:
                    check_compound(group, affected_child, parent_1, parent_2)

            variant_data = pd.concat(variant_data_grouped)

        else:
            logger.error("If family type (TRIO) is used a proper PED file defining the family relations is required!")
            logger.warn("Inheritance assessment will be skipped!")
        
    elif (family_type == "MULTI") and (family is not None):
        if family is not None:
            variant_data["DOMINANT"] = variant_data.apply(lambda variant: check_dominant(variant, family), axis=1)
            variant_data["XLINKED"] = variant_data.apply(lambda variant: check_xlinked(variant, family), axis=1)
            variant_data["RECESSIVE"] = variant_data.apply(lambda variant: check_recessive(variant, family, family_type), axis=1)

        else:
            logger.error("If family type (MULTI) is used a proper PED file defining the family relations (if known) is required!")
            logger.warn("Inheritance assessment will be skipped!")
        
    else:
        logger.error("Unsupported family type (%s) was given!" % family_type)
        logger.warn("Inheritance assessment will be skipped!")

    variant_columns = variant_data.columns
    variant_data["INHERITANCE"] = variant_data.apply(lambda variant: add_inheritance_mode(variant, variant_columns), axis=1)

    return variant_data


def add_inheritance_mode(variant, variant_columns):
    inheritance_list = []

    if "DOMINANT" in variant_columns:
        if variant["DOMINANT"] == 1:
            inheritance_list.append("DOMINANT")

    if "DOMINANT_DENOVO" in variant_columns:
        if variant["DOMINANT_DENOVO"] == 1:
            inheritance_list.append("DOMINANT_DENOVO")
    
    if "RECESSIVE" in variant_columns:
        if variant["RECESSIVE"] == 1:
            inheritance_list.append("RECESSIVE")
    
    if "COMPOUND" in variant_columns:
        if variant["COMPOUND"] == 1:
            inheritance_list.append("COMPOUND")

    if "XLINKED" in variant_columns:
        if variant["XLINKED"] == 1:
            inheritance_list.append("XLINKED")
    
    inheritance_mode = "&".join(inheritance_list)

    return inheritance_mode


# Returns longest homopolymer
def longest_homopolymer(sequence):
    if len(sequence) == 0:
        return 0

    runs = ''.join('*' if x == y else ' ' for x, y in zip(sequence, sequence[1:]))
    star_strings = runs.split()
    if len(star_strings) == 0:
        return 1

    return 1 + max(len(stars) for stars in star_strings)


# Returns most frequent base in a string
def frequent_base(sequence):
    # Get the most common base and its percentage of length of the sequence
    nucleotide_list = list(sequence)
    most_common_base = max([nucleotide_list.count(base) for base in set(nucleotide_list)])
    most_common_base_percentage = round(float(most_common_base) / len(sequence), 2)

    return most_common_base_percentage


# Homopolymer filter
def homopolymer_filter(sequence):
    base_type_zero_count = 0  # number of base types (A, C, G, T) not represented in the sequence
    low_complexity_flag = 0  # low complexity = (less than 3 base types represented or homopolymers of length >= 4)
    homopolymer_flag = 0  # homopolymer of length >= 5

    if sequence != '.':

        # Get the longest k-mer
        max_homopolymer_length = longest_homopolymer(sequence)

        # Get the frequency of the homopolymer base
        max_base_freq = frequent_base(sequence)

        # Count base types
        nucleotide_list = list(sequence)
        for base in ['A', 'C', 'G', 'T']:
            if nucleotide_list.count(base) == 0:
                base_type_zero_count += 1

        # Set homopolymer flag
        if max_homopolymer_length >= 5 or max_base_freq >= 0.8:
            homopolymer_flag = 1

        # Set low complexity flag
        if base_type_zero_count > 1 or max_homopolymer_length >= 4:
            low_complexity_flag = 1

    return homopolymer_flag, low_complexity_flag


def check_filters(variant, genes2exclude, HPO_query, reference):
    variant_genes = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(variant_genes.split(";"))

    # TODO: change to use biopython instead of pysam to prevent the need for both libraries
    in_fasta = pysam.FastaFile(reference)
    #reference_sequence = SeqIO.index(reference, "fasta")

    consequences = str(variant["Consequence"])
    found_consequences = consequences.split("&")
    seg_dup = float(variant[duplication_identifier])
    repeat = str(variant[repeat_identifier])
    filter_comment = ""

    try:
        maf = float(variant["MAX_AF"])
    except Exception as e:
        logger.warn("Allele frequency entry could not be identified, use 0.0 instead")
        maf = 0.0

    # filter for low confidence regions (these are regions where often false positives were observed)
    if "low_conf_region" in variant["FILTER"]:
        filter_passed = 0 # low confidence region
        filter_comment = "low confidence region"

        return filter_passed, filter_comment

    elif "off-target" in variant["FILTER"]:
        filter_passed = 0 # off target
        filter_comment = "off target"

        return filter_passed, filter_comment

    # high confidence: 0 - 0.15
    # mid confidence: 0.15 - 0.75
    # see https://onlinelibrary.wiley.com/doi/full/10.1002/humu.23674 for more detailed information about the abb score
    # if high confidence (> 0.15) filtering is to strict mid confidence (> 0.75) filtering could be applied
    #if float(variant["ABB_SCORE"]) > 0.75:
    #    filter_passed = 0 # abb score to high -> crap region (low confidence variant)
    #    filter_comment = "ABB_SCORE"
    #
    #    return filter_passed, filter_comment

    # exclude gene, if it is in the exclusion list
    if len(genes2exclude & genenames) > 0:
        for gene in genenames:
            # let variants with a high AIDIVA_SCORE (>=0.8) pass to be more sensitive
            if (gene.upper() in genes2exclude) and (float(variant["AIDIVA_SCORE"]) < 0.8):
                filter_passed = 0 # gene in exclusion list
                filter_comment = "gene exclusion"

                return filter_passed, filter_comment

    # let variants with a high AIDIVA_SCORE (>=0.8) pass to be more sensitive
    if (repeat != "NA") and (repeat != "") and (repeat != "nan") and (float(variant["AIDIVA_SCORE"]) < 0.8):
        filter_passed = 0 # tandem repeat
        filter_comment = "tandem repeat"

        return filter_passed, filter_comment
    
    # let variants with a high AIDIVA_SCORE (>=0.8) pass to be more sensitive
    if (seg_dup > 0) and (float(variant["AIDIVA_SCORE"]) < 0.8):
        filter_passed = 0 # segmental duplication
        filter_comment = "segmental duplication"

        return filter_passed, filter_comment
    
    # TODO: is this filter too strict? should we also filter inframe indels?
    if (len(variant["REF"]) > 1 or len(variant["ALT"]) > 1):
        # Get sequence context (vicinity) of a variant for homopolymer check (5 bases up- and down-stream)
        # Get fewer bases when variant is at the start or end of the sequence
        num_bases = 5
        pos_start = max(int(variant["POS"]) - (num_bases + 1), 1)
        pos_end = min(int(variant["POS"]) + num_bases, in_fasta.get_reference_length(variant["CHROM"]))

        try:
            sequence_context = in_fasta.fetch(variant["CHROM"], pos_start, pos_end)
            sequence_context = sequence_context.upper()
        except FileNotFoundError:
            sequence_context = '.'
        
        homopolymer_flag, low_complexity_flag = homopolymer_filter(sequence_context)

        if low_complexity_flag and homopolymer_flag:
            filter_passed = 0
            filter_comment = "homopolymer and low complexity region"

            return filter_passed, filter_comment

        if homopolymer_flag:
            filter_passed = 0
            filter_comment = "homopolymer"

            return filter_passed, filter_comment

        if low_complexity_flag:
            filter_passed = 0
            filter_comment = "low complexity region"

            return filter_passed, filter_comment
    
    # TODO: check position of frameshift in CDS region
    #
    #
    #

    if not str(variant["REPEATMASKER"]).strip():
        filter_passed = 0 # masked repeat region
        filter_comment = "masked repeat region"

        return filter_passed, filter_comment

    ## TODO: change frequency based on inheritance mode (hom/het)
    if maf <= 0.01:
        if any(term for term in CODING_VARIANTS if term in found_consequences):
            if not np.isnan(variant["FINAL_AIDIVA_SCORE"]):
                if len(HPO_query) >= 1:
                    if float(variant["HPO_RELATEDNESS"]) > 0.0:
                        filter_passed = 1
                        filter_comment = "passed all"

                    elif float(variant["HPO_RELATEDNESS_INTERACTING"]) > 0.0:
                        filter_passed = 1
                        filter_comment = "HPO related to interacting genes"
                    
                    else:
                        filter_passed = 0 # no relation to reported HPO terms
                        filter_comment = "no HPO relation"
                
                else:
                    filter_passed = 1 # skip hpo filter if no terms are present
                    filter_comment = "no HPO terms given"
            
            else:
                filter_passed = 0 # no prediction present (eg. variant type not covered by the used ML models)
                filter_comment = "missing AIDIVA_SCORE"

        else:
            filter_passed = 0 # non coding
            filter_comment = "non coding"
            
    else:
        filter_passed = 0 # allele frequency to high
        filter_comment = "high MAF"

    return filter_passed, filter_comment


def check_compound(gene_variants, affected_child, parent_1, parent_2):
    num_variant_candidates = gene_variants.shape[0]

    if num_variant_candidates >= 2:
        candidate_indices = [x for x in combinations(gene_variants.index.tolist(), 2)]
        for index_pair in candidate_indices:
            affected_child_zygosity_a = gene_variants.loc[index_pair[0], "GT." + affected_child].replace("|", "/")
            affected_child_zygosity_b = gene_variants.loc[index_pair[1], "GT." + affected_child].replace("|", "/")

            parent_1_zygosity_a = gene_variants.loc[index_pair[0], "GT." + parent_1].replace("|", "/")
            parent_1_zygosity_b = gene_variants.loc[index_pair[1], "GT." + parent_1].replace("|", "/")

            parent_2_zygosity_a = gene_variants.loc[index_pair[0], "GT." + parent_2].replace("|", "/")
            parent_2_zygosity_b = gene_variants.loc[index_pair[1], "GT." + parent_2].replace("|", "/")

            if (affected_child_zygosity_a == "0/1") and (affected_child_zygosity_b == "0/1"):
                if ((parent_1_zygosity_a == "0/0") and (parent_2_zygosity_a == "0/1")) and ((parent_1_zygosity_b == "0/1") and (parent_2_zygosity_b == "0/0")):
                    gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                    gene_variants.loc[index_pair[1], "COMPOUND"] = 1

                elif ((parent_1_zygosity_a == "0/1") and (parent_2_zygosity_a == "0/0")) and ((parent_1_zygosity_b == "0/0") and (parent_2_zygosity_b == "0/1")):
                    gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                    gene_variants.loc[index_pair[1], "COMPOUND"] = 1

            if ("." in affected_child_zygosity_a) or ("." in affected_child_zygosity_b):
                logger.warn("Skip variant pair, uncalled genotype in affected sample!")


def check_compound_single(gene_variants, variant_columns):
    genotype_column = [column for column in variant_columns if column.startswith("GT.")][0]
    num_variant_candidates = gene_variants.shape[0]

    if num_variant_candidates >= 2:
        candidate_indices = [x for x in combinations(gene_variants.index.tolist(), 2)]
        for index_pair in candidate_indices:
            genotype_variant_a = gene_variants.loc[index_pair[0], genotype_column].replace("|", "/")
            genotype_variant_b = gene_variants.loc[index_pair[1], genotype_column].replace("|", "/")

            if (((genotype_variant_a == "0/1") and (genotype_variant_b == "0/1"))):
                gene_variants.loc[index_pair[0], "COMPOUND"] = 1
                gene_variants.loc[index_pair[1], "COMPOUND"] = 1


def check_denovo(variant, family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        # check if sample is found in pedigree
        # sample info complete?
        ## TODO: do we need this check???
        if name in check_samples:
            check_samples[name] = 1

        # heterozygous in affected individual - good
        if zygosity == "0/1" and family[name] == 1:
            judgement = 1
            continue

        # hom ref, not affected - good
        elif zygosity == "0/0" and family[name] == 0 :
            judgement = 1
            continue

        # heterozygous in non-affected - bad
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 0
            break

        # hom ref in affected - bad
        elif zygosity == "0/0" and family[name] == 1:
            judgement = 0
            break

        # homozygous can"t be denovo
        elif zygosity == "1/1":
            judgement = 0
            break

        else:
            # reject missing genotype here or beforehand
            pass
        
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement


def check_dominant(variant, family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        # affected family members should have the mutation (hom ref not allowed)
        if zygosity == "0/0" and family[name] == 1:
            judgement = 0
            break

        # affected family members might be het
        elif zygosity == "0/1" and family[name] == 1:
            judgement = 1
            continue

        # affected family members might be hom alt
        # that"s the major difference to de novo...
        elif zygosity == "1/1" and family[name] == 1:
            judgement = 1
            continue

        # non-affected family members must not have the mutation - hom ref is OK
        elif zygosity == "0/0" and family[name] == 0:
            judgement = 1
            continue

        # non-affected family members must not have the mutation - het is bad
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 0
            break

        # non-affected family members must not have the mutation - hom alt is worst
        elif zygosity == "1/1" and family[name] == 0:
            judgement = 0
            break
        
        else:
            # reject missing genotype here or beforehand
            pass

    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement


def check_dominant_single(variant, variant_columns):
    genotype_column = [column for column in variant_columns if column.startswith("GT.")][0]
    judgement = 0

    if (variant[genotype_column] == "0/1"):
                judgement = 1

    return judgement


def check_recessive(variant, family, family_type):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        # affected individuals have to be homozygous
        if zygosity == "1/1" and family[name] == 1:
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == "0/0" or zygosity == "0/1" ) and family[name] == 1:
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 1
            continue

        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == "0/0" and family[name] == 0 and family_type == "MULTI":
            judgement = 1
            continue

        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == "0/0" and family[name] == 0 and family_type == "TRIO":
            judgement = 0
            break

        # non-affected individuals must not be hom alt
        elif zygosity == "1/1" and family[name] == 0:
            judgement = 0
            break

        else:
            # reject missing genotype here or beforehand
            pass
        
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return judgement


def check_recessive_single(variant, variant_columns):
    genotype_column = [column for column in variant_columns if column.startswith("GT.")][0]
    variant_genotype = variant[genotype_column].replace("|", "/")
    is_recessive = 0

    # recessive variants have to be homozygous in affected samples
    if (variant_genotype == "1/1"):
        is_recessive = 1
    
    return is_recessive


def check_xlinked(variant, family):
    judgement = 0
    check_samples = dict()
    inheritance_logic = dict()

    if not ((variant["CHROM"] == "X") or (variant["CHROM"] == "x") or (variant["CHROM"] == "chrX") or (variant["CHROM"] == "chrx") or (variant["CHROM"] == "23")):
        return 0

    # create data structure for completeness check
    for name in family.keys():
        check_samples[name] = 0

        zygosity = variant["GT." + name].replace("|", "/")

        if name in check_samples:
            check_samples[name] = 1

        if family[name] == 0:
            inheritance_logic[name] = zygosity

        # affected individuals have to be homozygous
        if zygosity == "1/1" and family[name] == 1:
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == "0/0" or zygosity == "0/1" ) and family[name] == 1:
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == "0/1" and family[name] == 0:
            judgement = 1
            continue

        # non-affected individuals might be hom ref
        elif zygosity == "0/0" and family[name] == 0:
            judgement = 1
            continue

        # non-affected individuals must not be hom alt
        elif zygosity == "1/1" and family[name] == 0:
            judgement = 0
            break

        else:
            # reject missing genotype here or beforehand
            pass
        
    # sanity check
    het_checker = 0
    hom_checker = 0

    for values in inheritance_logic.values():
        # mother
        if values == "0/1":
            het_checker = 1
        
        # father
        if values == "0/0":
            hom_checker = 1

    if het_checker == 1 and hom_checker == 1:
        judgement = 1
    else:
        judgement = 0

    # another sanity check
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return judgement


if __name__=="__main__":
    import get_HPO_similarity_score as gs

    parser = argparse.ArgumentParser(description = "Filter variants and finalize the AIDIVA_SCORE based on the given HPO terms (if this information is present)")
    parser.add_argument("--in_file", type=str, dest="in_file", required=True, help="Tab separated input annotated and scored file [required]")
    parser.add_argument("--out_file", type=str, dest="out_filename", required=True, help="Name to save the results [required]")
    parser.add_argument("--family", type=str, dest="family", required=False, help="Tab separated list of samples annotated with affection status. [required]")
    parser.add_argument("--family_type", type=str, choices=["TRIO", "FAMILY", "SINGLE"], dest="family_type", required=False, help="Choose if the data you provide is a trio or a larger family [required]")
    parser.add_argument("--gene_exclusion", type=str, dest="gene_exclusion_list", required=False, help="List of genes that should be excluded in the prioritization")
    parser.add_argument("--hpo_list", type=str, dest="hpo_list", default=None, required=False, help="List of HPO terms that are observed in the patient. These terms are used to adjust the AIDIVA_SCORE\n")
    parser.add_argument("--hpo_resources", type=str, dest="hpo_resources", default="../../data/", required=True, help="Folder where the HPO resources (HPO_graph,...) are found\n")
    args = parser.parse_args()

    input_data = pd.read_csv(args.in_file, sep="\t", low_memory=False)

    prioritized_variants = prioritize_variants(input_data, args.hpo_resources, args.family, args.family_type, args.hpo_list, args.gene_exclusion_list, args.low_conf_region)
    prioritized_variants.to_csv(args.out_filename, sep="\t", index=False)
