###########################################################
### Small script to easily recreate all used resources! ###
###########################################################

import argparse
import gzip
import json
import logging
import math
import pandas as pd
import sys
import networkx as nx


logger = logging.getLogger("HPO resource generator")
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


# get mapping gene -> HPOs
# download from HPO charite phenotype to gene
# wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt
def generate_gene2hpo_dict(gene2phenotype_list, gene2hpo_dict):
    logger.info("Generate gene to HPO mapping...")
    gene_2_HPO = dict()

    with open(gene2phenotype_list) as rd:
        for line in rd:
            if line.startswith("#"):
                pass

            else:
                splitted_line = line.strip().split("\t")
                # Format: HPO-id<tab>HPO label<tab>entrez-gene-id<tab>entrez-gene-symbol<tab>Additional Info from G-D source<tab>G-D source<tab>disease-ID for link
                gene_symbol = splitted_line[3].upper()
                hpo_term = splitted_line[0]
                to_add = gene_2_HPO.get(gene_symbol,[])
                to_add.append(hpo_term)
                to_add = list(set(to_add))
                gene_2_HPO[gene_symbol] = to_add

    with open(gene2hpo_dict, "w") as gene2hpo_f:
        json.dump(gene_2_HPO, gene2hpo_f)

    logger.info(f"Gene to HPO mapping successfully generated and saved as {gene2hpo_dict}")


# get mapping HPO -> translation
# download data:
# wget http://purl.obolibrary.org/obo/hp.obo -> ontology file
def generate_hpo2name_dict(hpo_ontology, hpo2name_dict):
    hpo_mapping = dict()
    token = False

    with open(hpo_ontology, "r") as ontology:
        for line in ontology:
            if line.startswith("[Term]"):
                if token:
                    hpo_mapping[hpo_term] = dict()
                    hpo_mapping[hpo_term]["name"] = name
                    hpo_mapping[hpo_term]["synonyms"] = synonyms

                token = True
                synonyms = []
                name = ""
            
            elif line.startswith("id:"):
                hpo_term = line.strip().split("id: ")[1]
            
            elif line.startswith("name:"):
                name = line.strip().split("name: ")[1]
            
            elif line.startswith('synonym:'):
                synonyms.append(line.strip().split('"')[1])
            
            elif line.startswith("[Typedef]"):
                break

            else:
                continue

        # add information of the last entry of the file to the HPO edges file
        if token:
            hpo_mapping[hpo_term] = dict()
            hpo_mapping[hpo_term]["name"] = name
            hpo_mapping[hpo_term]["synonyms"] = synonyms
            
    with open(hpo2name_dict, "w") as hpo2name_f:
        json.dump(hpo_mapping, hpo2name_f)

    logger.info(f"HPO to name mapping successfully generated and saved as {hpo2name_dict}")


# download data:
# wget http://purl.obolibrary.org/obo/hp.obo -> ontology file
# wget http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa
def create_hpo_graph(hpo_ontology, phenotype_hpoa, hpo_graph_file, hpo_replacement_information_file):
    logger.info("Extract HPO edges...")
    hpo_edges = dict()
    replacements = dict()
    considerations = dict()
    alternatives = dict()
    token = False
    obsolete = False

    with open(hpo_ontology, "r") as ontology:
        for line in ontology:
            if line.startswith("id: HP"):
                if token and not obsolete:
                    hpo_edges[hpo_term] = parents
                    if replacement != "":
                        replacements[hpo_term] = replacement

                token = True
                hpo_term = line.strip().split("id: ")[1]
                parents = []
                replacement = ""
                obsolete =False

            elif line.startswith("is_a:"):
                parents.append(line.strip().split("is_a: ")[1].split(" !")[0])

            elif line.startswith('is_obsolete:'):
                obsolete =True

            elif line.startswith('replaced_by:'):
                replacement = line.strip().split('replaced_by: ')[1]
                obsolete =False

            elif line.startswith('alt_id:'):
                alt = line.strip().split('alt_id: ')[1]
                if hpo_term in alternatives.keys():
                    current_alternatives = alternatives[hpo_term]
                    current_alternatives.append(alt)
                    alternatives[hpo_term] = current_alternatives

                else:
                    alternatives[hpo_term] = [alt]

            elif line.startswith('consider:'):
                consider = line.strip().split('consider: ')[1]
                if hpo_term in considerations.keys():
                    current_considerations = considerations[hpo_term]
                    current_considerations.append(consider)
                    considerations[hpo_term] = current_considerations

                else:
                    considerations[hpo_term] = [consider]

        # add information of the last entry of the file to the HPO edges file
        if token and not obsolete:
            hpo_edges[hpo_term] = parents
            if replacement != "":
                replacements[hpo_term] = replacement

        hpo_edges["replacements"] = replacements

    logger.info("Generate HPO graph...")
    # idea is a graph with attribute the IC value per node calculated

    # generate graph with counts:
    counts_dict = dict()
    total_counts = 0

    phenotype_info_table = pd.read_csv(phenotype_hpoa, sep="\t", comment="#", low_memory=False)
    hpo_counts = phenotype_info_table["hpo_id"].value_counts()

    for hpo_term, count_value in hpo_counts.items():
        counts_dict[hpo_term] = count_value
        total_counts += count_value

    # get replacements of obsolete nodes
    replacements = dict(hpo_edges.get("replacements", []))

    # let"s build a graph
    hpo_graph = nx.DiGraph()
    hpo_graph.graph["alternatives"] = alternatives
    hpo_graph.graph["considerations"] = considerations
    hpo_graph.graph["replacements"] = replacements

    hpo_replacement_information = {"alternatives": alternatives, "considerations": considerations, "replacements": replacements}

    # TODO: remove support for networkx v1.x and get rid of the duplicated if-else conditions (networkx v2.x and v3.x use the same syntax for the graphs)
    for node in hpo_edges.keys():
        hpo_graph.add_node(node)
        if str(nx.__version__).startswith("1."):
            hpo_graph.node[node]["count"] = 0.0

        elif str(nx.__version__).startswith("2."):
            hpo_graph.nodes[node]["count"] = 0.0

        elif str(nx.__version__).startswith("3."):
            hpo_graph.nodes[node]["count"] = 0.0

        else:
            logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

        ancestors = [(x,node) for x in hpo_edges[node]]
        hpo_graph.add_edges_from(ancestors)

        if node in replacements.keys():
            if str(nx.__version__).startswith("1."):
                hpo_graph.node[node]['replaced_by'] = replacements[node]
                hpo_graph.node[node]['IC'] = -math.log(1.0 / total_counts)

            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]['replaced_by'] = replacements[node]
                hpo_graph.nodes[node]['IC'] = -math.log(1.0 / total_counts)

            elif str(nx.__version__).startswith("3."):
                hpo_graph.nodes[node]["replaced_by"] = replacements[node]
                hpo_graph.nodes[node]["IC"] = -math.log(1.0 /total_counts)

            else:
                logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

    for node in hpo_graph.nodes():
        if str(nx.__version__).startswith("1."):
            hpo_graph.node[node]["count"] = 0.0

        elif str(nx.__version__).startswith("2."):
            hpo_graph.nodes[node]["count"] = 0.0

        elif str(nx.__version__).startswith("3."):
            hpo_graph.nodes[node]["count"] = 0.0

        else:
            logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

    # populate with raw counts
    for node in counts_dict.keys():
        if str(nx.__version__).startswith("1."):
            hpo_graph.node[node]["count"] = counts_dict[node]

        elif str(nx.__version__).startswith("2."):
            hpo_graph.nodes[node]["count"] = counts_dict[node]

        elif str(nx.__version__).startswith("3."):
            hpo_graph.nodes[node]["count"] = counts_dict[node]

        else:
            logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

    hpo_graph.nodes(data="count")

    # now fill it with the actual value.
    for node in hpo_edges.keys():
        descendants = nx.descendants(hpo_graph,node)
        if str(nx.__version__).startswith("1."):
            count = hpo_graph.node[node]["count"]

        elif str(nx.__version__).startswith("2."):
            count = hpo_graph.nodes[node]["count"]

        elif str(nx.__version__).startswith("3."):
            count = hpo_graph.nodes[node]["count"]

        else:
            logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

        for descendant in descendants:
            if str(nx.__version__).startswith("1."):
                count += hpo_graph.node[descendant]["count"]

            elif str(nx.__version__).startswith("2."):
                count += hpo_graph.nodes[descendant]["count"]

            elif str(nx.__version__).startswith("3."):
                count += hpo_graph.nodes[descendant]["count"]

            else:
                logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

        if count > 0 :
            if str(nx.__version__).startswith("1."):
                hpo_graph.node[node]["IC"] =  -math.log(float(count) / total_counts)

            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]["IC"] =  -math.log(float(count) / total_counts)

            elif str(nx.__version__).startswith("3."):
                hpo_graph.nodes[node]["IC"] =  -math.log(float(count) / total_counts)

            else:
                logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

        else :
            if str(nx.__version__).startswith("1."):
                hpo_graph.node[node]["IC"] = -math.log(1.0 / total_counts) # missing nodes, set as rare as possible

            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]["IC"] = -math.log(1.0 / total_counts) # missing nodes, set as rare as possible

            elif str(nx.__version__).startswith("3."):
                hpo_graph.nodes[node]["IC"] = -math.log(1.0 / total_counts) # missing nodes, set as rare as possible

            else:
                logger.error("There seems to be a problem with your installation of NetworkX, make sure that you have either v1, v2 or v3 installed!")

    nx.write_gexf(hpo_graph, hpo_graph_file)
    with open(hpo_replacement_information_file, "w") as hpo_replacement_information_f:
        json.dump(hpo_replacement_information, hpo_replacement_information_f)

    logger.info(f"HPO graph successfully generated and saved as {hpo_graph_file}, HPO replacement information saved to {hpo_replacement_information_file}")


#  wget https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
def create_gene2hgnc_mapping(hgnc_symbol_file, hgnc_2_gene):
    with open(hgnc_symbol_file, "r") as file:
        gene_dict = dict()

        for line in file:
            if line.startswith("#") or line.startswith("hgnc_id"):
                continue

            splitted_line = line.split("\t")
            hgnc_id = splitted_line[0].replace("HGNC:", "")
            gene_dict[hgnc_id] = splitted_line[1].upper()

    with open(hgnc_2_gene, "w") as hgnc_2_gene_f:
        json.dump(gene_dict, hgnc_2_gene_f)
    logger.info(f"Gene symbol to HGNC mapping successfully generated and saved as {hgnc_2_gene}")


# StringDB v11.0b
# wget https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
# wget https://version-11-0b.string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
def create_gene2interacting_mapping(string_mapping, string_db_links, string_interactions):
    with gzip.open(string_mapping, "rt") as string_mapping_file:
        string2name = dict()

        for line in string_mapping_file:
            if line.startswith("#") or (line == "\n"):
                continue

            splitted_line = line.replace("\n", "").split("\t")

            string_id = str(splitted_line[2])
            gene_name = str(splitted_line[1])
            print(gene_name, string_id)

            string2name[string_id] = gene_name.upper()

    with gzip.open(string_db_links, "rt") as string_links_file: 
        string_interaction_mapping = dict()

        for line in string_links_file:
            if line.startswith("protein") or (line == "\n"):
                continue

            splitted_line = line.replace("\n", "").split(" ")
            protein_a = splitted_line[0]
            protein_b = splitted_line[1]
            confidence_experimental = int(splitted_line[6])
            confidence_database = int(splitted_line[7])

            if ((protein_a in string2name.keys()) and (protein_b in string2name.keys())) and ((confidence_experimental >= 900) or (confidence_database >= 900)):
                gene_name = string2name[protein_a]
                interacting_gene = string2name[protein_b]

                if gene_name in string_interaction_mapping.keys():
                    string_interaction_mapping[gene_name].append({"interacting_gene": interacting_gene, "confidence_experimental": confidence_experimental, "confidence_database": confidence_database})

                else:
                    string_interaction_mapping[gene_name] = [{"interacting_gene": interacting_gene, "confidence_experimental": confidence_experimental, "confidence_database": confidence_database}]
    
    with open(string_interactions, "w") as string_interactions_f:
        json.dump(string_interaction_mapping, string_interactions_f)

    logger.info(f"Gene to interacting gene mapping successfully generated and saved as {string_interactions}")

# this file can be obtained through ensembl biomart
# wget -O grch38_ensembl_transcript_length_and_strand.tsv 'https://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
# wget -O grch37_ensembl_transcript_length_and_strand.tsv 'https://grch37.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "cds_length" /><Attribute name = "transcript_length" /><Attribute name = "strand" /></Dataset></Query>'
def create_transcript_length_mapping(transcript_lengths, transcript_information):
    transcript_mapping = {}

    with open(transcript_information, "r") as transc:
        for line in transc:
            if line.startswith("Gene") or line.startswith("#"):
                continue

            else:
                # header:
                # #Gene stable ID [TAB] Transcript stable ID [TAB] CDS Length [TAB] Transcript length (including UTRs and CDS) [TAB] Strand
                splitted = line.split("\t")
                transcript_id = splitted[1]
                cds_length = splitted[2]

                if cds_length != "":
                    transcript_mapping[transcript_id] = cds_length

    with open(transcript_lengths, "w") as transcript_lengths_f:
        json.dump(transcript_mapping, transcript_lengths_f)

    logger.info(f"Transcript to CDS length mapping successfully generated and saved as {transcript_lengths}")


if __name__=="__main__":
    parser = argparse.ArgumentParser("Script to generate the HPO resources needed in the prioritization step of aiDIVA")
    parser.add_argument("--hpo_ontology", type=str, dest="hpo_ontology", metavar="hp.obo", required=False, help="File containing the HPO ontology\n")
    parser.add_argument("--gene_phenotype", type=str, dest="gene_phenotype", metavar="phenotype_to_genes.txt", required=False, help="File that contains information about the genes and the associated phenotypes\n")
    parser.add_argument("--gene_hpo", type=str, dest="gene_hpo", metavar="gene2hpo.json", required=False, help="File to save the generated gene2hpo_dict\n")
    parser.add_argument("--hpo_name", type=str, dest="hpo_name", metavar="hpo2name.json", required=False, help="File to save the generated hpo2name_dict\n")
    parser.add_argument("--phenotype_hpoa", type=str, dest="phenotype_hpoa", metavar="phenotype.hpoa", required=False, help="File containing inforamtion to all hpo ids, is used to compute the hpo counts\n")
    parser.add_argument("--hpo_graph", type=str, dest="hpo_graph", metavar="hpo_graph.gexf", required=False, help="File to save the generated hpo graph\n")
    parser.add_argument("--hpo_replacements", type=str, dest="hpo_replacements", metavar="hpo2replacement.json", required=False, help="File to save the hpo replacement information\n")
    parser.add_argument("--hgnc_symbols", type=str, dest="hgnc_symbols", metavar="hgnc_complete_set.txt", required=False, help="File containing the approved hgnc genes and their previous gene symbols if there are any\n")
    parser.add_argument("--hgnc_gene", type=str, dest="hgnc_gene", metavar="hgnc2gene.json", required=False, help="File to save the generated hgnc to gene mapping\n")
    parser.add_argument("--string_links", type=str, dest="string_links", metavar="9606.protein.links.detailed.v11.0.txt.gz", required=False, help="File with the protein relations found in the STRING database\n")
    parser.add_argument("--string_mapping", type=str, dest="string_mapping", metavar="human.name_2_string.tsv.gz", required=False, help="File with mapping of string ids to gene symbols\n")
    parser.add_argument("--gene_interacting", type=str, dest="gene_interacting", metavar="gene2interacting.json", required=False, help="File to save the generated gene to interacting genes mapping\n")
    parser.add_argument("--transcript_length_mapping", type=str, dest="transcript_length_mapping", metavar="transcript2length.json", required=False, help="File to save the mapping of transcript ids to CDS length\n")
    parser.add_argument("--transcript_information", type=str, dest="transcript_information", metavar="transcript_length_information.tsv", required=False, help="File containing the CDS length information per transcript\n")
    args = parser.parse_args()

    if args.gene_phenotype and args.gene_hpo:
        generate_gene2hpo_dict(args.gene_phenotype, args.gene_hpo)

    if args.hpo_ontology and args.hpo_name:
        generate_hpo2name_dict(args.hpo_ontology, args.hpo_name):

    if args.hpo_ontology and args.phenotype_hpoa and args.hpo_graph and args.hpo_replacements:
        create_hpo_graph(args.hpo_ontology, args.phenotype_hpoa, args.hpo_graph, args.hpo_replacements)

    if args.hgnc_symbols and args.hgnc_gene:
        create_gene2hgnc_mapping(args.hgnc_symbols, args.hgnc_gene)

    if args.string_mapping and args.string_links and args.gene_interacting:
        create_gene2interacting_mapping(args.string_mapping, args.string_links, args.gene_interacting)

    if args.transcript_information and args.transcript_length_mapping:
        create_transcript_length_mapping(args.transcript_length_mapping, args.transcript_information)
