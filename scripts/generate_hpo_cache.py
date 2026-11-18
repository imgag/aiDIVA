import argparse
import gzip
import json
import logging
import math
import pandas as pd
import networkx as nx


# generate mapping gene -> HPOs
def generate_gene2hpo_dict(gene2phenotype_list, cache_path):
    print("Generate gene to HPO mapping...")
    gene_2_HPO = dict()

    with open(gene2phenotype_list) as rd:
        for line in rd:
            if line.startswith("#"):
                continue

            else:
                splitted_line = line.strip().split("\t")
                # Format: HPO-id<tab>HPO label<tab>entrez-gene-id<tab>entrez-gene-symbol<tab>Additional Info from G-D source<tab>G-D source<tab>disease-ID for link
                gene_symbol = splitted_line[3].upper()
                hpo_term = splitted_line[0]
                to_add = gene_2_HPO.get(gene_symbol,[])
                to_add.append(hpo_term)
                to_add = sorted(set(to_add))
                gene_2_HPO[gene_symbol] = to_add

    write_cache(gene_2_HPO, "hpo2gene.json", "JSON", cache_path)
    print(f"Gene to HPO mapping successfully loaded and cached")


# build hpo graph with networkx to easily access the information
def create_hpo_graph(hpo_ontology, phenotype_hpoa, cache_path):
    print("Extract HPO edges...")
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

    print("Generate HPO graph...")
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

    # let's build a graph
    hpo_graph = nx.DiGraph()
    hpo_graph.graph["alternatives"] = alternatives
    hpo_graph.graph["considerations"] = considerations
    hpo_graph.graph["replacements"] = replacements

    hpo_replacement_information = {"alternatives": alternatives, "considerations": considerations, "replacements": replacements}

    # TODO: simplify the duplicated if-else conditions
    for node in hpo_edges.keys():
        hpo_graph.add_node(node)
        hpo_graph.nodes[node]["count"] = 0.0

        ancestors = [(x,node) for x in hpo_edges[node]]
        hpo_graph.add_edges_from(ancestors)

        if node in replacements.keys():
            hpo_graph.nodes[node]["replaced_by"] = replacements[node]
            hpo_graph.nodes[node]["IC"] = -math.log(1.0 /total_counts)

    for node in hpo_graph.nodes():
        hpo_graph.nodes[node]["count"] = 0.0

    # populate with raw counts
    for node in counts_dict.keys():
        hpo_graph.nodes[node]["count"] = counts_dict[node]

    hpo_graph.nodes(data="count")

    # now fill it with the actual value.
    for node in hpo_edges.keys():
        descendants = nx.descendants(hpo_graph,node)
        count = hpo_graph.nodes[node]["count"]

        for descendant in descendants:
            count += hpo_graph.nodes[descendant]["count"]

        if count > 0 :
            hpo_graph.nodes[node]["IC"] = -math.log(float(count) / total_counts)

        else :
            hpo_graph.nodes[node]["IC"] = -math.log(1.0 / total_counts) # missing nodes, set as rare as possible

    write_cache(hpo_replacement_information, "hpo2replacement.json", "JSON", cache_path)
    write_cache(hpo_graph, "hpo_graph.gexf", "GEXF", cache_path)
    print(f"HPO graph successfully generated and HPO replacement information prepared.")


# prepare internal mapping of hgnc ids to the gene symbol
def create_gene2hgnc_mapping(hgnc_symbol_file, cache_path):
    with open(hgnc_symbol_file, "r") as file:
        gene_dict = dict()

        for line in file:
            if line.startswith("#") or line.startswith("hgnc_id"):
                continue

            splitted_line = line.split("\t")
            hgnc_id = splitted_line[0].replace("HGNC:", "")
            gene_dict[hgnc_id] = splitted_line[1].upper()

    write_cache(gene_dict, "hgnc2gene.json", "JSON", cache_path)
    print(f"Gene symbol to HGNC mapping successfully prepared.")


# prepare internal mapping of interacting genes based on information from the StringDB
def create_gene2interacting_mapping(string_mapping, string_db_links, cache_path):
    string2name = dict()

    string_mapping_aliases = pd.read_csv(string_mapping, sep="\t", low_memory=False)
    # aiDIVA uses the HGNC gene symbols so we filter the alias list to only include those for the mapping
    string_mapping_aliases = string_mapping_aliases[string_mapping_aliases["source"] == "Ensembl_HGNC_symbol"]

    for index, row in string_mapping_aliases.iterrows():
        string_id = str(row["#string_protein_id"])
        gene_name = str(row["alias"])

        if string_id not in string2name.keys():
            string2name[string_id] = gene_name.upper()
        
        else:
            print(f"string_id ({string_id}) already in mapping dictionary!")

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

    write_cache(string_interaction_mapping, "gene2interacting.json", "JSON", cache_path)
    print(f"Gene to interacting gene mapping successfully prepared.")


# create internal mapping of the transcript ids to their respective length based on ensembl data
def create_transcript_length_mapping(transcript_information, cache_path):
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

    write_cache(transcript_mapping, "transcript2length.json", "JSON", cache_path)
    print(f"Successfully loaded and cached transcript to CDS length mapping.")


# generate mapping HPO -> translation
def generate_hpo2name_dict(hpo_ontology, cache_path):
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

    write_cache(hpo_mapping, "hpo2translation.json", "JSON", cache_path)
    print(f"HPO to name mapping successfully loaded and prepared")


def write_cache(object_to_cache, cached_name, cache_mode, cache_path):
    if cache_mode == "JSON":
        with open(str(cache_path + cached_name), "w") as cache_file:
            json.dump(object_to_cache, cache_file)

    elif cache_mode == "GEXF":
        nx.write_gexf(object_to_cache, str(cache_path + cached_name))

    else:
        print("ERROR: Unsupported cache type given!")


if __name__=="__main__":
    parser = argparse.ArgumentParser("Script to generate the HPO resources needed in the prioritization step of aiDIVA")
    parser.add_argument("--cache_path", type=str, dest="cache_path", metavar="cache_path/", required=True, help="Path to store cached files")
    parser.add_argument("--hpo_ontology", type=str, dest="hpo_ontology", metavar="hp.obo", required=False, help="File containing the HPO ontology")
    parser.add_argument("--gene_phenotype", type=str, dest="gene_phenotype", metavar="phenotype_to_genes.txt", required=False, help="File that contains information about the genes and the associated phenotypes")
    parser.add_argument("--phenotype_hpoa", type=str, dest="phenotype_hpoa", metavar="phenotype.hpoa", required=False, help="File containing inforamtion to all hpo ids, is used to compute the hpo counts")
    parser.add_argument("--hgnc_symbols", type=str, dest="hgnc_symbols", metavar="hgnc_complete_set.txt", required=False, help="File containing the approved hgnc genes and their previous gene symbols if there are any")
    parser.add_argument("--string_links", type=str, dest="string_links", metavar="9606.protein.links.detailed.v12.0.txt.gz", required=False, help="File with the protein relations found in the STRING database")
    parser.add_argument("--string_aliases", type=str, dest="string_aliases", metavar="9606.protein.aliases.v12.0.txt.gz", required=False, help="File with mapping of string ids to gene symbols")
    parser.add_argument("--transcript_information", type=str, dest="transcript_information", metavar="transcript_length_information.tsv", required=False, help="File containing the CDS length information per transcript")
    args = parser.parse_args()

    if args.gene_phenotype:
        generate_gene2hpo_dict(args.gene_phenotype, args.cache_path)

    if args.hpo_ontology and args.phenotype_hpoa:
        create_hpo_graph(args.hpo_ontology, args.phenotype_hpoa, args.cache_path)
        generate_hpo2name_dict(args.hpo_ontology, args.cache_path)

    if args.hgnc_symbols:
        create_gene2hgnc_mapping(args.hgnc_symbols, args.cache_path)

    if args.string_aliases and args.string_links:
        create_gene2interacting_mapping(args.string_aliases, args.string_links, args.cache_path)

    if args.transcript_information:
        create_transcript_length_mapping(args.transcript_information, args.cache_path)