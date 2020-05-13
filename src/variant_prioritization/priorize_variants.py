import argparse
import csv
import re
from scipy.stats import poisson
import os
import shutil
import pickle
import networkx
import sys
import pandas as pd
from operator import itemgetter
import get_HPO_similarity_score as gs


variant_consequences = {'transcript_ablation': 'non_exonic',
                        'splice_acceptor_variant': 'exonic;splicing',
                        'splice_donor_variant': 'exonic;splicing',
                        'stop_gained': 'exonic',
                        'frameshift_variant': 'exonic',
                        'stop_lost': 'exonic',
                        'start_lost': 'exonic',
                        'transcript_amplification': 'non_exonic',
                        'inframe_insertion': 'exonic',
                        'inframe_deletion': 'exonic',
                        'missense_variant': 'exonic',
                        'protein_altering_variant': 'exonic',
                        'splice_region_variant': 'splicing', # maybe not
                        'incomplete_terminal_codon_variant': 'exonic',
                        'start_retained_variant': 'exonic', # maybe not
                        'stop_retained_variant': 'exonic', # maybe not
                        'synonymous_variant': 'exonic',
                        'coding_sequence_variant': 'exonic',
                        'mature_miRNA_variant': 'non_exonic',
                        '5_prime_UTR_variant': 'UTR_exonic',
                        '3_prime_UTR_variant': 'UTR_exonic',
                        'non_coding_transcript_exon_variant': 'non_exonic',
                        'intron_variant': 'intronic',
                        'NMD_transcript_variant': 'non_exonic',
                        'non_coding_transcript_variant': 'non_exonic',
                        'upstream_gene_variant': 'non_exonic',
                        'downstream_gene_variant': 'non_exonic',
                        'TFBS_ablation': 'non_exonic',
                        'TFBS_amplification': 'non_exonic',
                        'TF_binding_site_variant': 'non_exonic',
                        'regulatory_region_ablation': 'non_exonic',
                        'regulatory_region_amplification': 'non_exonic',
                        'feature_elongation': 'non_exonic',
                        'regulatory_region_variant': 'non_exonic',
                        'feature_truncation': 'non_exonic',
                        'intergenic_variant': 'non_exonic'}


# TODO use pandas dataframes instead of the csv library
def prioritize_variants(in_file, out_file, filtered_file, fam_file, inheritance, family_type, white_list, gene_exclusion):
    variant_data = pd.read_csv(in_file, sep="\t", low_memory=False)

    # read the gene exclusion list
    # an empty set will not filter out anything, if gene exclusion list is not provided
    genes2exclude = set()
    genes_known = set()

    ### HPO files load
    script_path = os.path.dirname(os.path.abspath(__file__))

    # TODO add these filepaths to the config file
    gene_2_HPO_f = os.path.join(script_path, '../../res/gene_2_HPO.p')
    HPO_graph_file = os.path.join(script_path, '../../res/v2_ready_graph.pk')
    hpo_dict_file = os.path.join(script_path, '../../res/HPO_gene_assiciation.p')

    gene_2_HPO = pickle.load(open(gene_2_HPO_f, 'rb'))
    HPO_graph_nodes, HPO_graph_edges = pickle.load(open(HPO_graph_file, 'rb'))

    HPO_graph = networkx.Graph()
    HPO_graph.add_nodes_from(HPO_graph_nodes)
    HPO_graph.add_edges_from(HPO_graph_edges)

    if gene_exclusion:
        gene_exclusion = open(gene_exclusion, "r")
        for gene in gene_exclusion:
            gene = gene.rstrip()
            genes2exclude.add(gene)
        gene_exclusion.close()

    #fill HPO query terms
    if white_list == None:
        white_list = 'None'
        HPO_query = list()
    else:
        HPO_query = set()
        if os.path.isfile(white_list):
            with open(white_list, 'r') as w:
                HPO_dict = pickle.load(open(hpo_dict_file,'rb'))
                HPO_query = list()
                for line in w:
                    HPO_term = line.rstrip('\n')
                    try:
                        HPO_query.append(HPO_term)
                    except:
                        print('%s not found in database' % (HPO_term))
            HPO_query= list(set(HPO_query))
        else:
            print('The specified HPO list %s is not a valid file' % (white_list))
            print('eDiVA will proceed as without any HPO list')

    # read family relationships
    # TODO change to ped file
    family = dict()
    with open(fam_file, "r") as fam_file:
        for line in fam_file:
            if line.startswith('sample'):
                continue
            line = line.rstrip('\n')
            splitline = line.split('\t')
            print(splitline)
            family[splitline[0]] = splitline[1]

    print(family)



    variant_data["HPO_RELATEDNESS", "FINAL_RANK"] = variant_data.apply(lambda variant: pd.Series(compute_hpo_relatedness_and_final_rank(variant, genes2exclude, HPO_graph, gene_2_HPO, HPO_query)), axis=1)
    variant_data["RECESSIVE", "DOMINANT_DENOVO", "DOMINANT_INHERITED", "XLINKED", "COMPOUND", "FILTER_PASSED"] = variant_data.apply(lambda variant: pd.Series(check_inheritance_and_filters(variant, genes2exclude, HPO_list, family, family_type, names)), axis=1)

    ## TODO: Chek recessive variants for possible compounds
    ## Compoundizer method applied on each gene set???
    #variant_data["COMPOUND"] = variant_data.apply(lambda variant: pd.Series(check_compounds(variant)))

    variant_data.to_csv(out_file, sep='\t', encoding='utf-8', index=False)
    variant_data[variant_data["FILTER_PASSED"] == 1].to_csv(filtered_file, sep='\t', encoding='utf-8', index=False)


def compute_hpo_relatedness_and_final_rank(variant, genes2exclude, HPO_graph, gene_2_HPO, HPO_query):
    genecolumn = re.sub("\(.*?\)", "", str(variant["SYMBOL"]))
    genenames = set(genecolumn.split(";"))

    query_dist = 0
    gene_distances = []
    processed_HPO_genes = dict()
    for gene_id in genenames:
        #if not gene_id in genes2exclude: # TODO: decide where to handle the geneexclusion list
        if gene_id in processed_HPO_genes.keys():
            gene_distances.append(processed_HPO_genes[gene_id])
        else:
            #process ex novo
            #get HPOs related to the gene
            gene_HPO_list = gs.extract_HPO_related_to_gene(gene_2_HPO, gene_id)
            ## TODO: calculate whole dict (query_dist) in the beginning and pass as parameter
            (g_dist,query_dist) = gs.list_distance(HPO_graph, HPO_query, gene_HPO_list, query_dist)
            gene_distances.append(g_dist)
            processed_HPO_genes[gene_id] = g_dist
    hpo_relatedness = str(max(gene_distances))

    final_rank = str((float(variant["Rank"]) + float(hpo_relatedness)) / 2)

    return [hpo_relatedness, final_rank]


def check_inheritance_and_filters(variant, genes2exclude, HPO_list, family, familytype):
    genecolumn = re.sub("\(.*?\)", "", variant["SYMBOL"])
    genenames = set(genecolumn.split(";"))

    consequences = variant["Consequence"]
    seg_dup = float(variant["SegDupMax"])
    tandem = variant["SimpleTandemRepeatRegion"]
    cadd = variant["CADD_PHRED"]

    try:
        maf = max(float(variant["AA_AF"]), float(variant["AFR_AF"]), float(variant["AMR_AF"]), float(variant["EA_AF"]), float(variant["EAS_AF"]), float(variant["EUR_AF"]), float(variant["SAS_AF"]), float(variant["gnomAD_AFR_AF"]), float(variant["gnomAD_AMR_AF"]), float(variant["gnomAD_ASJ_AF"]), float(variant["gnomAD_EAS_AF"]), float(variant["gnomAD_FIN_AF"]), float(variant["gnomAD_NFE_AF"]), float(variant["gnomAD_OTH_AF"]), float(variant["gnomAD_SAS_AF"]))
    except Exception as e:
        print("Allele frequency could not be identified, use 0.0 instead")
        maf = 0.0

    dominant_denovo = denovo(family)

    ## TODO: do we need a check for affected family members?
    dominant_inherited = dominant(family)
    if variant["Chr"] == "X" or variant["Chr"] == "x" or variant["Chr"] == "23":
        xlinked = xlinked(family)
    else:
        xlinked = 0
    recessive = recessive(family, familytype)

    found_consequences = [variant_consequences[consequence] for consequence in consequences.split("&")]

    #conditions = dict()
    #conditions['RECESSIVE'] = (0.03, -1)
    #conditions['DOMINANT_DENOVO'] = (0.01, -1)
    #conditions['DOMINANT_INHERITED'] = (0.01, -1)
    #conditions['COMPOUND'] = (0.02, -1)
    #conditions['XLINKED'] = (0.01, -1)
    #conditions['COMPOUND_SINGLE_SAMPLE'] = conditions['DOMINANT_DENOVO']
    #conditions['UNKNOWN'] = conditions['DOMINANT_INHERITED']

    #(MAF_threshold,CADD_threshold) = conditions["DOMINANT_INHERITED"]

    # exclude gene, if it is on the exclusion list
    if len(genes2exclude & genenames) > 0:
        filter_passed = 0 # gene in exclusion list
        return [recessive, dominant_denovo, dominant_inherited, xlinked, 0, filter_passed]

    elif not (tandem == 'NA' or tandem == '' or tandem == 'nan'):
        filter_passed = 0 # tandem repeat
        return [recessive, dominant_denovo, dominant_inherited, xlinked, 0, filter_passed]

    ## TODO: remove
    #elif cadd > 0 and cadd < CADD_threshold :
    #     filter_passed = 0

    ## TODO: filter later compound only less than 0.01
    elif maf <= 0.02:
        if ('exonic' in found_consequences or 'splicing' in found_consequences or 'exonic;splicing' in found_consequences):
            if (not "synonymous_variant" in consequences.split("&")) and ("unknown" != consequences) and ("UNKNOWN" != consequences):
                if (seg_dup == 0):
                    filter_passed = 1
                    if len(HPO_list) > 1 and 'NONE' not in HPO_list:
                        if float(variant["HPO_RELATEDNESS"]) > 0:
                            filter_passed = 1
                        else:
                            filter_passed = 0 # no relation to reported HPO terms
                # e.g. intronic variants fitting the criteria
                else:
                    filter_passed = 0 # segment duplication
            else:
                filter_passed = 0 # synonymous variant  or unknown effect
        else:
            filter_passed = 0 # not exonic
    else:
        filter_passed = 0 # allele frequency to high

    return [recessive, dominant_denovo, dominant_inherited, xlinked, 0, filter_passed]


def denovo(family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[sam] = 0

        if "REF." + name in variant.columns() and "ALT." + name in variant.columns():
            zygosity = variant[name]
            refcoverage = variant["REF." + str(name)].replace('.', '0') # could be numeric or .
            altcoverage = variant["ALT." + str(name)].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = '.'
            altcoverage = '.'

        # check if sample is found in pedigree
        # sample info complete?
        ## TODO: do we need this check???
        if name in check_samples:
            check_samples[name] = 1

        # heterozygous in affected individual - good
        if zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue

        # hom ref, not affected - good
        elif zygosity == '0/0' and family[name] == '0' :
            if int(altcoverage) <= max(3, (float(altcoverage) + float(refcoverage)) / 10) :
                judgement = 1
            else :
                judgement =0
            continue

        # heterozygous in non-affected - bad
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 0
            break

        # hom ref in affected - bad
        elif zygosity == '0/0' and family[name] == '1':
            judgement = 0
            break

        # homozygous can't be denovo
        elif zygosity == '1/1':
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':

            # if vcf file was not supplemented by pileup data
            # reject it variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 0
                continue

            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref, non called genotype
            # poisson for low coverage and percentage for high coverage
            # poisson 10 reads (poisson average rate of success = 5) and alt reads = 0 - should get still accepted
            elif (poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007) and (altcoverage / coverage <= 0.05):
                judgement = 1
                continue

            # not necessary to check for hom alt
            # coverage too low?
            else:
                judgement = 0
                break

        # do not accept missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':

            # except if vcf file was not supplemented by pileup data
            # reject variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.':
                judgement = 0
                continue

            # do not be that grateful, if there is coverage data
            # if the SNP caller could not call a genotype in an area, where there was coverage
            # we rather trust the SNP caller, than starting to call SNP based on pileup coverage
            else:
                judgement = 0
                break

    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement


def dominant(family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[sam] = 0

        if "REF." + name in variant.columns() and "ALT." + name in variant.columns():
            zygosity = variant[name]
            refcoverage = variant["REF." + str(name)].replace('.', '0') # could be numeric or .
            altcoverage = variant["ALT." + str(name)].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = '.'
            altcoverage = '.'

        if name in check_samples:
            check_samples[name] = 1

        # affected family members should have the mutation (hom ref not allowed)
        if zygosity == '0/0' and family[name] == '1':
            judgement = 0
            break

        # affected family members might be het
        elif zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue

        # affected family members might be hom alt
        # that's the major difference to de novo...
        elif zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue

        # non-affected family members must not have the mutation - hom ref is OK
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue

        # non-affected family members must not have the mutation - het is bad
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 0
            break

        # non-affected family members must not have the mutation - hom alt is worst
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':

            # if vcf file was not supplemented
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue

            # do not do any other judgements, i.e.
            # being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            # also tolerate low coverage
            judgement = 1
            continue

        # accept some missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':

            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue

            ### do not do any other judgements, i.e.
            ### being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            ### also tolerate low coverage
            ##judgement = 1
            ##continue

            # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                break

            # hom alt
            elif poisson.cdf(float(refcoverage), float(coverage / 2)) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 1
                continue

            # het
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                continue

        else:
            pass

    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return judgement

def recessive(family, familytype):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[sam] = 0

        if "REF." + name in variant.columns() and "ALT." + name in variant.columns():
            zygosity = variant[name]
            refcoverage = variant["REF." + str(name)].replace('.', '0') # could be numeric or .
            altcoverage = variant["ALT." + str(name)].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = '.'
            altcoverage = '.'

        if name in check_samples:
            check_samples[name] = 1

        # affected individuals have to be homozygous
        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 1
            continue

        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'FAMILY':
            judgement = 1
            continue

        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'TRIO':
            judgement = 0
            break

        # non-affected individuals must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx

            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                break

            # hom alt
            elif poisson.cdf(float(refcoverage), float(coverage / 2)) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 0
                break

            # het, which is OK
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                continue

            # coverage too low?
            else:
                # accept missing values in family interrogations
                if familytype == 'FAMILY':
                    judgement = 1
                    continue
                # do not accept missing values in trio setups
                elif familytype == 'TRIO':
                    judgement = 0
                    break

                # for security reasons
                judgement = 0
                break

        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == './.' and family[name] == '1':
            judgement = 0
            break

    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return judgement


def xlinked(family):
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for name in family.keys():
        check_samples[sam] = 0

        if "REF." + name in variant.columns() and "ALT." + name in variant.columns():
            zygosity = variant[name]
            refcoverage = variant["REF." + str(name)].replace('.', '0') # could be numeric or .
            altcoverage = variant["ALT." + str(name)].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = variant[name]
            refcoverage = '.'
            altcoverage = '.'

        if name in check_samples:
            check_samples[name] = 1

        if family[name] == '0':
            inheritance_logic[name] = zygosity

        # affected individuals have to be homozygous
        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue

        # affected individuals should not be hom ref or het
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            break

        # non-affected individuals might be het
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 1
            continue

        # non-affected individuals might be hom ref
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue

        # non-affected individuals must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break

        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx

            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                inheritance_logic[name] = '0/0'
                judgement = 1
                continue

            # hom alt
            if poisson.cdf(float(refcoverage), float(coverage) / 2) <= 0.007 and refcoverage / coverage <= 0.05:
                inheritance_logic[name] = '1/1'
                judgement = 0
                break

            # het, which is OK
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                inheritance_logic[name] = '0/1'
                judgement = 1
                continue

            # coverage too low?
            else:
                judgement = 0
                break

        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == './.' and family[name] == '1':
            judgement = 0
            break

    # sanity check
    het_checker = 0
    hom_checker = 0

    for values in inheritance_logic.values():
        if values == '0/1':
            het_checker = 1
        if values == '0/0':
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


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = 'filter SNPs according to their family distribution')
    parser.add_argument('--infile', type=str, dest='infile', required=True, help='Input annotated ranked file [required]')
    parser.add_argument('--outfile', type=str, dest='outfile', required=True, help='Output file wit all variants. [required]')
    parser.add_argument('--filteredoutfile', type=str, dest='filteredfile', required=True, help='filtered comma separated list of SNPs annotated with inheritance pattern, only reporting the requested variants. [required]')
    parser.add_argument('--family', type=str, dest='famfile', required=True, help='tab separated list of samples annotated with affection status. [required]')
    parser.add_argument('--inheritance', type=str, choices=['DOMINANT_DENOVO', 'DOMINANT_INHERITED', 'RECESSIVE', 'XLINKED', 'COMPOUND', 'UNKNOWN'], dest='inheritance', required=True, help="""choose a inheritance model [required]
    DOMINANT_INHERITED: used for families
    DOMINANT_DENOVO: apply to novel variants seen in the affected individuals
    RECESSIVE: detect recessive, homozygous variants (if trio is specified the script will require that all non-affected are heterozygous)
    XLINKED: used for X linked recessive variants in trios only
    COMPOUND: detect compound heterozygous recessive variants
    UNKNOWN: skip inheritance filtering if not known
    """)
    parser.add_argument('--familytype', type=str, choices=['TRIO', 'FAMILY', 'SINGLE'], dest='familytype', required=True, help="choose if the data you provide is a trio or a larger family")
    parser.add_argument('--geneexclusion', type=str, dest='geneexclusion', required=False, help='[Analysis of DNA sequence variants detected by high-throughput sequencing; DOI: 10.1002/humu.22035].')
    parser.add_argument('--HPO_list', '--white_list', type=str, dest='white_list', default=None, required=False, help='--HPO_list \t a .txt file with the list of HPO terms describing the disease. It will be used to flag all genes related to the HPO terms. It works with Refseq and UCSC naming.\n')

    args = parser.parse_args()
    prioritize_variants(args.infile, args.outfile, args.filteredfile, args.famfile, args.inheritance, args.familytype, args.white_list, args.geneexclusion)
