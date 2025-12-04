import gzip
import logging
import math
import pandas as pd
import sys
import networkx as nx


logger = logging.getLogger(__name__)


# generate mapping gene -> HPOs
def generate_gene2hpo_dict(gene2phenotype_list):
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

    logger.info(f"Gene to HPO mapping successfully loaded and prepared")
    return gene_2_HPO


# build hpo graph with networkx to easily access the information
def create_hpo_graph(hpo_ontology, phenotype_hpoa):
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

    logger.info(f"HPO graph successfully generated and HPO replacement information prepared.")
    return hpo_graph, hpo_replacement_information


# prepare internal mapping of hgnc ids to the gene symbol
def create_gene2hgnc_mapping(hgnc_symbol_file):
    with open(hgnc_symbol_file, "r") as file:
        gene_dict = dict()

        for line in file:
            if line.startswith("#") or line.startswith("hgnc_id"):
                continue

            splitted_line = line.split("\t")
            hgnc_id = splitted_line[0].replace("HGNC:", "")
            gene_dict[hgnc_id] = splitted_line[1].upper()

    logger.info(f"Gene symbol to HGNC mapping successfully prepared.")
    return gene_dict


# prepare internal mapping of interacting genes based on information from the StringDB
def create_gene2interacting_mapping(string_mapping, string_db_links):
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
            print(f"WARNING: string_id ({string_id}) already in mapping dictionary!")

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
    
    logger.info(f"Gene to interacting gene mapping successfully prepared.")
    return string_interaction_mapping


# create internal mapping of the transcript ids to their respective length based on ensembl data
def create_transcript_length_mapping(transcript_information):
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

    logger.info(f"Successfully loaded and prepared transcript to CDS length mapping.")
    return transcript_mapping
