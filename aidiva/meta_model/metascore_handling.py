import argparse
import gzip
import json
import logging
import os
import pickle
import numpy as np
import pandas as pd
import random
import re


RANDOM_SEED = 14038

VARIANT_CONSEQUENCES = {"transcript_ablation": 1,
                        "splice_acceptor_variant": 2,
                        "splice_donor_variant": 3,
                        "stop_gained": 4,
                        "frameshift_variant": 5,
                        "stop_lost": 6,
                        "start_lost": 7,
                        "transcript_amplification": 8,
                        "inframe_insertion": 9,
                        "inframe_deletion": 10,
                        "missense_variant": 11,
                        "protein_altering_variant": 12,
                        "splice_region_variant": 13,
                        "splice_donor_5th_base_variant": 14,
                        "splice_donor_region_variant": 15,
                        "splice_polypyrimidine_tract_variant": 16,
                        "incomplete_terminal_codon_variant": 17,
                        "start_retained_variant": 18,
                        "stop_retained_variant": 19,
                        "synonymous_variant": 20,
                        "coding_sequence_variant": 21,
                        "mature_miRNA_variant": 22,
                        "5_prime_UTR_variant": 23,
                        "3_prime_UTR_variant": 24,
                        "non_coding_transcript_exon_variant": 25,
                        "intron_variant": 26,
                        "NMD_transcript_variant": 27,
                        "non_coding_transcript_variant": 28,
                        "upstream_gene_variant": 29,
                        "downstream_gene_variant": 30,
                        "TFBS_ablation": 31,
                        "TFBS_amplification": 32,
                        "TF_binding_site_variant": 33,
                        "regulatory_region_ablation": 34,
                        "regulatory_region_amplification": 35,
                        "feature_elongation": 36,
                        "regulatory_region_variant": 37,
                        "feature_truncation": 38,
                        "intergenic_variant": 39,
                        # use "unknown" consequence as default if new consequence terms are added to the database that are not yet implemented (this prevents the program from exiting with an error)
                        "unknown": 40}

logger = logging.getLogger(__name__)


def extract_rf_rank_and_score(gene_rank_score_information, no_variant):
    gene_entries = gene_rank_score_information.split(";")
    gene_info_dict_rf = {}

    for gene_entry in gene_entries:
        gene = gene_entry.replace(")", "").split("(")[0].upper()
        gene_info = gene_entry.replace(")", "").split("(")[1]

        splitted_gene_info = gene_info.split(", ")
        gene_consequence = splitted_gene_info[0].split(": ")[1]
        gene_rank = splitted_gene_info[2].split(": ")[1]
        gene_score = splitted_gene_info[3].split(": ")[1]

        if not no_variant:
            gene_variant = splitted_gene_info[5].split(": ")[1]

        else:
            gene_variant = "Unknown"

        if gene in gene_info_dict_rf.keys():
            gene_info_dict_rf[gene].append({"rank": gene_rank, "score": gene_score, "variant_type": gene_consequence, "variant": gene_variant})

        else:
            gene_info_dict_rf[gene] = [{"rank": gene_rank, "score": gene_score, "variant_type": gene_consequence, "variant": gene_variant}]

    return gene_info_dict_rf


def choose_desired_transcript(gene_transcripts):
    if len(gene_transcripts) > 2:
        logger.warning("Check sample. Reason: more than two transcripts for one gene.")
        print(f"WARNING: Check sample. Reason: more than two transcripts for one gene.")
        consequence_value = 50

        for transcript in gene_transcripts:
            transcript_gene = gene_transcripts[0].replace(")", "").split("(")[0].upper()
            transcript_consequence = choose_desired_variant_type(gene_transcripts[0].replace(")", "").split("(")[1])

            if VARIANT_CONSEQUENCES[transcript_consequence] < consequence_value:
                gene = transcript_gene
                variant_type = transcript_consequence

    else:
        transcript_one_gene = gene_transcripts[0].replace(")", "").split("(")[0].upper()
        transcript_one_consequence = choose_desired_variant_type(gene_transcripts[0].replace(")", "").split("(")[1])

        transcript_two_gene = gene_transcripts[1].replace(")", "").split("(")[0].upper()
        transcript_two_consequence = choose_desired_variant_type(gene_transcripts[1].replace(")", "").split("(")[1])

        gene = transcript_one_gene
        variant_type = transcript_one_consequence

        if VARIANT_CONSEQUENCES[transcript_two_consequence] < VARIANT_CONSEQUENCES[transcript_one_consequence]:
            gene = transcript_two_gene
            variant_type = transcript_two_consequence

    return gene, variant_type


def choose_desired_variant_type(variant_consequences):
    found_consequences = variant_consequences.split("+")

    # use most severe consequence for filtering if overlapping consequences are present
    if len(found_consequences) > 1:
        most_severe_consequence = found_consequences[0]

        for consequence in found_consequences:
            if VARIANT_CONSEQUENCES[consequence] < VARIANT_CONSEQUENCES[most_severe_consequence]:
                most_severe_consequence = consequence

    else:
        most_severe_consequence = found_consequences[0]

    return most_severe_consequence


def extract_eb_rank_and_score(gene_rank_score_information, no_variant):
    gene_entries = gene_rank_score_information.split(";")
    gene_info_dict_eb = {}

    for gene_entry in gene_entries:
        #transcript_info = gene_entry.split("=")[0]
        #rank_score_info = gene_entry.split("=")[1]

        #rank = rank_score_info.replace(")", "").split("(")[0]
        #score = rank_score_info.replace(")", "").split("(")[1].split("/")[0]
        #genotype = rank_score_info.replace(")", "").split("(")[1].split("/")[1]

        gene_name = gene_entry.replace(")", "").split("(")[0].upper()
        gene_info = gene_entry.replace(")", "").split("(")[1]

        splitted_gene_info = gene_info.split(", ")
        gene_consequence = splitted_gene_info[0].split(": ")[1]
        gene_rank = splitted_gene_info[2].split(": ")[1]
        gene_score = splitted_gene_info[3].split(": ")[1]

        if not no_variant:
            gene_variant = splitted_gene_info[5].split(": ")[1]

        else:
            gene_variant = "Unknown"

        # check for multiple transcripts
        splitted_gene_names = gene_name.split("/")

        if len(splitted_gene_names) == 1:
            current_gene = splitted_gene_names[0].upper()
            variant_type_info = splitted_gene_info[1]
            #print(splitted_gene_info)

            variant_type = choose_desired_variant_type(variant_type_info)

        elif len(splitted_gene_names) > 1:
            current_gene, variant_type = choose_desired_transcript(gene_entry)

        if current_gene in gene_info_dict_eb.keys():
            gene_info_dict_eb[current_gene].append({"rank": gene_rank, "score": gene_score, "variant_type": gene_consequence, "variant": gene_variant})

        else:
            gene_info_dict_eb[current_gene] = [{"rank": gene_rank, "score": gene_score, "variant_type": gene_consequence, "variant": gene_variant}]

    return gene_info_dict_eb


def get_gene_list(data, model_type):
    if model_type == "RF":
        gene_list = subdata["RF_top10_aiDIVA_v101"].values[0]

    elif model_type == "EB dom":
        gene_list = subdata["EB_top10_GSvar_v2_dominant"].values[0]

    elif model_type == "EB rec":
        gene_list = subdata["EB_top10_GSvar_v2_recessive"].values[0]

    else:
        print("ERROR!")

    return gene_list


def get_llm_gene_candidates(data):
    first_gene = str(data["1st ranked Gene"].values[0]).upper()
    second_gene = str(data["2nd ranked Gene"].values[0]).upper()
    third_gene = str(data["3rd ranked Gene"].values[0]).upper()

    return first_gene, second_gene, third_gene


def create_table_rf_based(in_data_rf_based, rf_gene_list, no_variant):
    meta_table_dict_list = []

    #rf_gene_list = get_gene_list(in_data_rf_based, "RF")
    rf_gene_info_dict = extract_rf_rank_and_score(rf_gene_list, no_variant)

    rf_genes = set(list(rf_gene_info_dict.keys()))

    current_gene_set = set()
    current_gene_set |= rf_genes

    rf_gene_one_llm, rf_gene_two_llm, rf_gene_three_llm = get_llm_gene_candidates(in_data_rf_based)

    for gene in current_gene_set:
        if not no_variant:
            gene_variant = rf_gene_info_dict[gene][0]["variant"]

        else:
            gene_variant = "Unknown"

        if gene in rf_gene_info_dict.keys():
            if len(rf_gene_info_dict[gene]) > 1:
                consequence_value = 50

                for transcript in rf_gene_info_dict[gene]:
                    if VARIANT_CONSEQUENCES[transcript["variant_type"]] < consequence_value:
                        rf_rank = int(transcript["rank"])
                        rf_score = float(transcript["score"])

            else:
                rf_rank = int(rf_gene_info_dict[gene][0]["rank"])
                rf_score = float(rf_gene_info_dict[gene][0]["score"])

        else:
            rf_rank = 11
            rf_score = 0.0

        if gene in rf_gene_info_dict.keys():
            if gene == rf_gene_one_llm:
                rf_rank_llm = 1

            elif gene == rf_gene_two_llm:
                rf_rank_llm = 2

            elif gene == rf_gene_three_llm:
                rf_rank_llm = 3

            else:
                rf_rank_llm = 4

        else:
            rf_rank_llm = 4

        #print(f"{gene}\t{rf_rank}\t{rf_score}\t{rf_rank_llm}\n")

        meta_table_dict_list.append({"gene_name": gene, "variant": gene_variant, "rf_rank": rf_rank, "rf_score": rf_score, "rf_rank_llm": rf_rank_llm})

    meta_table = pd.DataFrame.from_dict(meta_table_dict_list)

    return meta_table


def create_table_rf_and_evidence_based(in_data_rf_based, rf_gene_list, in_data_eb_dom, eb_dom_gene_list, in_data_eb_rec, eb_rec_gene_list, no_variant):
    meta_table_dict_list = []

    #rf_gene_list = get_gene_list(in_data_rf_based, sample_id, "RF")
    rf_gene_info_dict = extract_rf_rank_and_score(rf_gene_list, no_variant)
    #print("RF:", rf_gene_info_dict)

    #eb_dom_gene_list = get_gene_list(in_data_eb_dom, sample_id, "EB dom")
    eb_dom_gene_info_dict = extract_eb_rank_and_score(eb_dom_gene_list, no_variant)
    #print("EB dom:", eb_dom_gene_info_dict)

    #eb_rec_gene_list = get_gene_list(in_data_eb_rec, sample_id, "EB rec")
    eb_rec_gene_info_dict = extract_eb_rank_and_score(eb_rec_gene_list, no_variant)
    #print("EB rec:", eb_rec_gene_info_dict)

    rf_genes = set(list(rf_gene_info_dict.keys()))
    eb_dom_genes = set(list(eb_dom_gene_info_dict.keys()))
    eb_rec_genes = set(list(eb_rec_gene_info_dict.keys()))

    current_gene_set = set()
    current_gene_set |= rf_genes
    current_gene_set |= eb_dom_genes
    current_gene_set |= eb_rec_genes

    rf_gene_one_llm, rf_gene_two_llm, rf_gene_three_llm = get_llm_gene_candidates(in_data_rf_based)
    eb_dom_gene_one_llm, eb_dom_gene_two_llm, eb_dom_gene_three_llm = get_llm_gene_candidates(in_data_eb_dom)
    eb_rec_gene_one_llm, eb_rec_gene_two_llm, eb_rec_gene_three_llm = get_llm_gene_candidates(in_data_eb_rec)

    for gene in current_gene_set:
        gene_variant = ""

        if gene in rf_gene_info_dict.keys():
            if gene_variant == "":
                gene_variant = rf_gene_info_dict[gene][0]["variant"]

            elif gene_variant != "" and gene_variant != rf_gene_info_dict[gene][0]["variant"]:
                gene_variant = gene_variant + ";" + rf_gene_info_dict[gene][0]["variant"]

            #gene_variant = rf_gene_info_dict[gene][0]["variant"]
            if len(rf_gene_info_dict[gene]) > 1:
                consequence_value = 50

                for transcript in rf_gene_info_dict[gene]:
                    if transcript["variant_type"] == "":
                        transcript["variant_type"] = "unknown"

                    if VARIANT_CONSEQUENCES[transcript["variant_type"]] < consequence_value:
                        rf_rank = int(transcript["rank"])
                        rf_score = float(transcript["score"])

            else:
                rf_rank = int(rf_gene_info_dict[gene][0]["rank"])
                rf_score = float(rf_gene_info_dict[gene][0]["score"])

        else:
            rf_rank = 11
            rf_score = 0.0

        if gene in eb_dom_gene_info_dict.keys():
            if gene_variant == "":
                gene_variant = eb_dom_gene_info_dict[gene][0]["variant"]

            elif gene_variant != "" and gene_variant != eb_dom_gene_info_dict[gene][0]["variant"]:
                gene_variant = gene_variant + ";" + eb_dom_gene_info_dict[gene][0]["variant"]

            #gene_variant = eb_dom_gene_info_dict[gene][0]["variant"]
            if len(eb_dom_gene_info_dict[gene]) > 1:
                consequence_value = 50

                for transcript in eb_dom_gene_info_dict[gene]:
                    if transcript["variant_type"] == "":
                        transcript["variant_type"] = "unknown"

                    if VARIANT_CONSEQUENCES[transcript["variant_type"]] < consequence_value:
                        eb_dom_rank = int(transcript["rank"])
                        eb_dom_score = float(transcript["score"])

            else:
                eb_dom_rank = int(eb_dom_gene_info_dict[gene][0]["rank"])
                eb_dom_score = float(eb_dom_gene_info_dict[gene][0]["score"])

        else:
            eb_dom_rank = 11
            eb_dom_score = 0.0

        if gene in eb_rec_gene_info_dict.keys():
            if gene_variant == "":
                gene_variant = eb_rec_gene_info_dict[gene][0]["variant"]
            
            elif gene_variant != "" and gene_variant != eb_rec_gene_info_dict[gene][0]["variant"]:
                gene_variant = gene_variant + ";" + eb_rec_gene_info_dict[gene][0]["variant"]
            
            #gene_variant = eb_rec_gene_info_dict[gene][0]["variant"]
            if len(eb_rec_gene_info_dict[gene]) > 1:
                consequence_value = 50

                for transcript in eb_rec_gene_info_dict[gene]:
                    if transcript["variant_type"] == "":
                        transcript["variant_type"] = "unknown"

                    if VARIANT_CONSEQUENCES[transcript["variant_type"]] < consequence_value:
                        eb_rec_rank = int(transcript["rank"])
                        eb_rec_score = float(transcript["score"])

            else:
                eb_rec_rank = int(eb_rec_gene_info_dict[gene][0]["rank"])
                eb_rec_score = float(eb_rec_gene_info_dict[gene][0]["score"])

        else:
            eb_rec_rank = 11
            eb_rec_score = 0.0

        if eb_dom_rank < eb_rec_rank:
            eb_model = "dom"
            eb_rank = eb_dom_rank
            eb_score = eb_dom_score

        elif eb_dom_rank > eb_rec_rank:
            eb_model = "rec"
            eb_rank = eb_rec_rank
            eb_score = eb_rec_score

        # rank of causal variant in both models equals
        else:
            if eb_dom_score > eb_rec_score:
                eb_model = "dom"
                eb_rank = eb_dom_rank
                eb_score = eb_dom_score

            elif eb_dom_score < eb_rec_score:
                eb_model = "rec"
                eb_rank = eb_rec_rank
                eb_score = eb_rec_score

            # equal scores use dominant model as default
            else:
                eb_model = "dom"
                eb_rank = eb_dom_rank
                eb_score = eb_dom_score

        if gene in rf_gene_info_dict.keys():
            if gene == rf_gene_one_llm:
                rf_rank_llm = 1

            elif gene == rf_gene_two_llm:
                rf_rank_llm = 2

            elif gene == rf_gene_three_llm:
                rf_rank_llm = 3

            else:
                rf_rank_llm = 4

        else:
            rf_rank_llm = 4

        if (gene in eb_dom_gene_info_dict.keys()) or (gene in eb_rec_gene_info_dict.keys()):
            if gene == eb_dom_gene_one_llm or gene == eb_rec_gene_one_llm:
                eb_rank_llm = 1

            elif gene == eb_dom_gene_two_llm or gene == eb_rec_gene_two_llm:
                eb_rank_llm = 2

            elif gene == eb_dom_gene_three_llm or gene == eb_rec_gene_three_llm:
                eb_rank_llm = 3

            else:
                eb_rank_llm = 4

        else:
            eb_rank_llm = 4

        #print(f"{gene}\t{rf_rank}\t{rf_score}\t{rf_rank_llm}\t{eb_rank}\t{eb_score}\t{eb_rank_llm}\n")
        meta_table_dict_list.append({"gene_name": gene, "variant": gene_variant, "rf_rank": rf_rank, "rf_score": rf_score, "rf_rank_llm": rf_rank_llm, "eb_model": eb_model, "eb_rank": eb_rank, "eb_score": eb_score, "eb_rank_llm": eb_rank_llm})

    meta_table = pd.DataFrame.from_dict(meta_table_dict_list)

    return meta_table


def import_model(model_file):
    if model_file.endswith(".gz"):
        model_to_import = gzip.open(model_file, "rb")

    else:
        model_to_import = open(model_file, "rb")

    model = pickle.load(model_to_import)
    model_to_import.close()

    return model


def meta_scoring(metascore_table, meta_model, rf_only):
    #metascore_table = pd.read_csv(in_file, sep="\t", low_memory=False)
    metascore_model = import_model(meta_model)

    label_column = "label"

    if rf_only:
        feature_list = ["rf_rank", "rf_score", "rf_rank_llm"]

    else:
        feature_list = ["rf_rank", "rf_score", "rf_rank_llm", "eb_rank", "eb_score", "eb_rank_llm", "dominant", "recessive"]

        metascore_table["dominant"] = metascore_table["eb_model"]
        metascore_table["recessive"] = metascore_table["eb_model"]
        metascore_table.loc[metascore_table['eb_model'] == 'dom', 'dominant'] = 1
        metascore_table.loc[metascore_table['eb_model'] == 'rec', 'dominant'] = 0
        metascore_table.loc[metascore_table['eb_model'] == 'dom', 'recessive'] = 0
        metascore_table.loc[metascore_table['eb_model'] == 'rec', 'recessive'] = 1

    metascore_table_X = np.asarray(metascore_table[feature_list])

    score_prediction = pd.DataFrame(metascore_model.predict_proba(metascore_table_X), columns=["Probability_Benign", "Probability_Pathogenic"])
    metascore_table["metascore_prediction"] = score_prediction["Probability_Pathogenic"]
    metascore_table["meta_rank"] = metascore_table["metascore_prediction"].rank(method="min", ascending=False).astype(int)
    metascore_table.sort_values("meta_rank", inplace=True)
    metascore_table.reset_index(inplace=True, drop=True)

    #metascore_table.to_csv(out_file, sep="\t", index=False)

    return metascore_table


def main(in_rf_llm, in_eb_dom_llm, in_eb_rec_llm, rf_file, eb_dom_file, eb_rec_file, sample_id, out_file, no_variant):
    rf_llm_res = pd.read_csv(in_rf_llm, sep="\t", low_memory=False)
    eb_dom_llm_res = pd.read_csv(in_eb_dom_llm, sep="\t", low_memory=False)
    eb_rec_llm_res = pd.read_csv(in_eb_rec_llm, sep="\t", low_memory=False)

    in_data_rf = pd.read_csv(rf_file, sep="\t", low_memory=False)

    top_ranking_genes_random_forest = top_ranking.extract_top_ranking_entries_random_forest_based(sample_id, in_data_rf, 10)
    top_ranking_genes_evidence_dominant = top_ranking.extract_top_ranking_entries_evidence_based(sample_id, eb_dom_file, 10)
    top_ranking_genes_evidence_recessive = top_ranking.extract_top_ranking_entries_evidence_based(sample_id, eb_rec_file, 10)

    metascore_table = create_table_rf_and_evidence_based(rf_llm_res, top_ranking_genes_random_forest, eb_dom_llm_res, top_ranking_genes_evidence_dominant, eb_rec_llm_res, top_ranking_genes_evidence_recessive, no_variant)
    metascore_table.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    import get_top_ranking_genes as top_ranking

    parser = argparse.ArgumentParser(description = "Create table with prompts for chatGPT")
    parser.add_argument("--in_rf_llm", type=str, dest="in_rf_llm", required=True, help="List of sample IDs")
    parser.add_argument("--in_eb_dom_llm", type=str, dest="in_eb_dom_llm", required=False, help="File with the RF ranks and scores")
    parser.add_argument("--in_eb_rec_llm", type=str, dest="in_eb_rec_llm", required=False, help="File with the RF ranks and scores")
    parser.add_argument("--rf_file", type=str, dest="rf_file", required=False, help="File with the RF ranks and scores")
    parser.add_argument("--eb_dom_file", type=str, dest="eb_dom_file", required=False, help="File with the EB dom ranks and scores")
    parser.add_argument("--eb_rec_file", type=str, dest="eb_rec_file", required=False, help="File with the EB rec ranks and scores")
    parser.add_argument("--sample_id", type=str, dest="sample_id", required=False, help="File with the EB rec ranks and scores")
    parser.add_argument("--out_file", type=str, dest="out_file", required=True, help="Output file")
    parser.add_argument("--no_variant", action="store_true", dest="no_variant", required=False, help="Do not store variant")
    args = parser.parse_args()

    main(args.in_rf_llm, args.in_eb_dom_llm, args.in_eb_rec_llm, args.rf_file, args.eb_dom_file, args.eb_rec_file, args.sample_id, args.out_file, args.no_variant)
