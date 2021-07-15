import argparse
import gzip
import math
import networkx as nx
import pickle


# get mapping gene -> HPOs
# download from HPO charite phenotype to gene
# wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/phenotype_to_genes.txt
def generate_gene2hpo_dict(gene2phenotype_list, gene2hpo_dict):
    print("Generate gene to HPO mapping...")
    gene_2_HPO = dict()

    with open(gene2phenotype_list) as rd:
        for line in rd:
            if line.startswith("#"):
                pass
            else:
                ff = line.strip().split("\t")
                # Format: HPO-id<tab>HPO label<tab>entrez-gene-id<tab>entrez-gene-symbol<tab>Additional Info from G-D source<tab>G-D source<tab>disease-ID for link
                key = ff[3].upper()
                HPO = ff[0]
                to_add = gene_2_HPO.get(key,[])
                to_add.append(HPO)
                to_add = list(set(to_add))
                gene_2_HPO[key] = to_add

    pickle.dump(gene_2_HPO, open(gene2hpo_dict,"wb"))
    print("Gene to HPO mapping successfully generated and saved as %s" % (gene2hpo_dict))


# download data
# wget https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
def extract_hpo_graph_edges(hpo_ontology, hpo_edges_file):
    print("Extract HPO edges...")
    out_HPO = dict()
    replacements = []
    token = False
    obsolete = False

    ontology = open(hpo_ontology, "r")
    for line in ontology:
        if line.startswith("id: HP:"):
            if token and not obsolete:
                out_HPO[name] = parents
                if replacement != "":
                    replacements.append((name,replacement))

            token = True
            name = line.strip().split("id: ")[1]
            parents = []
            replacement = ""
            obsolete =False

        elif line.startswith("is_a:"):
            parents.append(line.strip().split("is_a: ")[1].split(" !")[0])
        elif line.startswith('is_obsolete:'):
            obsolete =True
        elif line.startswith('replaced_by:'):
            # add a field to say it's replaced
            replacement = line.strip().split('replaced_by: ')[1]
            obsolete =False #means we can backtrack it
        #elif line.startswith('alt_id:'):
            #add alternative nodes, will be later added with
            # replacement field for the most common one
            #alt = line.strip().split('alt_id: ')[1]
            #alternatives.append((name,alt))
        #elif line.startswith('consider:'):
            #add alternative nodes, will be later added with
            # replacement field for the most common one
            #alt = line.strip().split('consider: ')[1]
            #alternatives.append((alt,name))
            #obsolete =False # means we can backtrack it

    # add information of the last entry of the file to the HPO edges file
    if token and not obsolete:
        out_HPO[name] = parents
        if replacement != "":
            replacements.append((name,replacement))
    
    out_HPO['replacements'] = replacements
    #out_HPO['alternatives'] = alternatives

    ontology.close()
    pickle.dump(out_HPO, open(hpo_edges_file, "wb"))
    print("HPO edges successfully extracted and saved as %s" % (hpo_edges_file))


# counts as wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
# awk -F '\t'  '{print $5}' < phenotype_annotation.tab | sort  | uniq -c | awk '{print $2 "\t" $1}' > HPO_counts.txt
def generate_hpo_graph(hpo_counts, hpo_edges_file, hpo_graph_file):
    print("Generate HPO graph...")
    offset =1000 # penalization for the distance
    # idea is a graph with attribute the IC value per node calculated
    output = hpo_graph_file # contains the nodes and edges to import in a DiGraph

    # generate graph with counts:
    counts_dict = dict()
    tot = 0
    with open(hpo_counts) as count_file:
        for line in count_file:
            splitted_line = line.strip().split("\t")
            counts_dict[splitted_line[0]] = int(splitted_line[1])
            tot += int(splitted_line[1])

    # load dict with edges
    edges = pickle.load(open(hpo_edges_file, "rb"))

    #get replacements of obsolete nodes
    replacements = dict(edges.get('replacements', []))
    temp_value = edges.pop('replacements', None)

    # let"s build a graph
    hpo_graph = nx.DiGraph()

    #populate with alternatives
    #mark alternatives as replaced, it's the same for us.
    #alternatives = edges.get('alternatives',[])
    #tmpval = edges.pop('alternatives',None)

    for node in edges.keys():
        hpo_graph.add_node(node)
        if str(nx.__version__).startswith("1."):
            hpo_graph.node[node]["count"] = 0.0
        elif str(nx.__version__).startswith("2."):
            hpo_graph.nodes[node]["count"] = 0.0
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        ancestors = [(x,node) for x in edges[node]]
        hpo_graph.add_edges_from(ancestors)
        if node in replacements.keys():
            if str(nx.__version__).startswith("1."):
                hpo_graph.node[node]['replaced_by'] = replacements[node]
                hpo_graph.node[node]['IC'] = -math.log(1.0 / tot)
            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]['replaced_by'] = replacements[node]
                hpo_graph.nodes[node]['IC'] = -math.log(1.0 / tot)
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    for node in hpo_graph.nodes():
        if str(nx.__version__).startswith("1."):
            hpo_graph.node[node]["count"] = 0.0
        elif str(nx.__version__).startswith("2."):
            hpo_graph.nodes[node]["count"] = 0.0
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    # populate with raw counts
    for node in counts_dict.keys():
        if str(nx.__version__).startswith("1."):
            hpo_graph.node[node]["count"] = counts_dict[node]
        elif str(nx.__version__).startswith("2."):
            hpo_graph.nodes[node]["count"] = counts_dict[node]
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    hpo_graph.nodes(data="count")
    # now fill it with the actual value.

    for node in edges.keys():
        descendants = nx.descendants(hpo_graph,node)
        if str(nx.__version__).startswith("1."):
            count = hpo_graph.node[node]["count"]
        elif str(nx.__version__).startswith("2."):
            count = hpo_graph.nodes[node]["count"]
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        for descendant in descendants:
            if str(nx.__version__).startswith("1."):
                count += hpo_graph.node[descendant]["count"]
            elif str(nx.__version__).startswith("2."):
                count += hpo_graph.nodes[descendant]["count"]
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        if count > 0 :
            if str(nx.__version__).startswith("1."):
                hpo_graph.node[node]["IC"] =  -math.log(float(count) / tot)
            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]["IC"] =  -math.log(float(count) / tot)
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")
        else :
            if str(nx.__version__).startswith("1."):
                hpo_graph.node[node]["IC"] = -math.log(1.0 / tot) # missing nodes, set as rare as possible
            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]["IC"] = -math.log(1.0 / tot) # missing nodes, set as rare as possible
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    # add edges weight
    for node_a, node_b in hpo_graph.edges():
        if str(nx.__version__).startswith("1."):
            hpo_graph[node_a][node_b]['dist'] = offset + abs(hpo_graph.node[node_a]['IC'] - hpo_graph.node[node_b]['IC'])
        elif str(nx.__version__).startswith("2."):
            hpo_graph[node_a][node_b]['dist'] = offset + abs(hpo_graph.nodes[node_a]['IC'] - hpo_graph.nodes[node_b]['IC'])
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")
    # count is the IC of the node then

    # convert directed graph to an undirected graph
    hpo_graph = hpo_graph.to_undirected()

    if str(nx.__version__).startswith("1."):
        pickle.dump([hpo_graph.nodes(data=True), hpo_graph.edges(data=True)], open(output, "wb"))
        print("NOTE: You pickled the graph using NetworkX v1 the pickled graph it is also upwards compatible with v2!")
    elif str(nx.__version__).startswith("2."):
        pickle.dump([hpo_graph.nodes(data=True), hpo_graph.edges(data=True)], open(output, "wb"))
        print("NOTE: You pickled the graph using NetworkX v2 it is not backwards compatible!")
    else:
        print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    print("HPO graph successfully generated and saved as %s" % (hpo_graph_file))


#  wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt
def create_gene2hgnc_mapping(hgnc_symbol_file, hgnc_2_gene):
    file = open(hgnc_symbol_file, "r")
    gene_dict = dict()

    for line in file:
        if line.startswith("#") or line.startswith("hgnc_id"):
            continue

        splitted_line = line.split("\t")
        hgnc_id = splitted_line[0].replace("HGNC:", "")
        gene_dict[hgnc_id] = splitted_line[1].upper()

    file.close()
    pickle.dump(gene_dict, open(hgnc_2_gene, "wb"))

    print("Gene symbol to HGNC mapping successfully generated and saved as %s" % (hgnc_2_gene))


# wget https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz
# wget https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz
def create_gene2interacting_mapping(string_mapping, string_db_links, string_interactions):
    string_mapping_file = gzip.open(string_mapping, "rt")
    string2name = dict()

    for line in string_mapping_file:
        if line.startswith("#") or (line == "\n"):
            continue
        splitted_line = line.replace("\n", "").split("\t")
        string2name[splitted_line[2]] = splitted_line[1].upper()

    string_mapping_file.close()

    string_links_file = gzip.open(string_db_links, "rt")
    string_interaction_mapping = dict()

    for line in string_links_file:
        if line.startswith("protein") or (line == "\n"):
            continue
        splitted_line = line.replace("\n", "").split(" ")
        prot1 = splitted_line[0]
        prot2 = splitted_line[1]
        conf_exp = int(splitted_line[6])
        conf_db = int(splitted_line[7])

        if ((prot1 in string2name.keys()) and (prot2 in string2name.keys())) and ((conf_exp >= 900) or (conf_db >= 900)):
            gene_name = string2name[prot1]
            interacting_gene = string2name[prot2]
            if gene_name in string_interaction_mapping.keys():
                string_interaction_mapping[gene_name].append({"interacting_gene": interacting_gene, "confidence_experimental": conf_exp, "confidence_database": conf_db})
            else:
                string_interaction_mapping[gene_name] = [{"interacting_gene": interacting_gene, "confidence_experimental": conf_exp, "confidence_database": conf_db}]
    
    string_links_file.close()
    pickle.dump(string_interaction_mapping, open(string_interactions, "wb"))

    print("Gene to interacting gene mapping successfully generated and saved as %s" % (string_interactions))


if __name__=="__main__":
    parser = argparse.ArgumentParser("Script to generate the HPO resources needed in the prioritization step of AIdiva")
    parser.add_argument("--hpo_ontology", type=str, dest="hpo_ontology", metavar="hp.obo", required=True, help="File containing the HPO ontology\n")
    parser.add_argument("--gene_phenotype", type=str, dest="gene_phenotype", metavar="phenotype_to_genes.txt", required=True, help="File that contains information about the genes and the associated phenotypes\n")
    parser.add_argument("--gene_hpo", type=str, dest="gene_hpo", metavar="gene2hpo.pkl", required=True, help="File to save the generated gene2hpo_dict\n")
    parser.add_argument("--hpo_edges", type=str, dest="hpo_edges", metavar="hpo_edges.pkl", required=True, help="File where the extracted hpo edges are stored\n")
    parser.add_argument("--hpo_counts", type=str, dest="hpo_counts", metavar="HPO_counts.txt", required=True, help="File containing the hpo counts needed for the hpo graph construction\n")
    parser.add_argument("--hpo_graph", type=str, dest="hpo_graph", metavar="hpo_graph.pkl", required=True, help="File to save the generated hpo graph\n")
    parser.add_argument("--hgnc_symbols", type=str, dest="hgnc_symbols", metavar="hgnc_approved_symbols.txt", required=True, help="File containing the approved hgnc genes and their previous gene symbols if there are any\n")
    parser.add_argument("--hgnc_gene", type=str, dest="hgnc_gene", metavar="hgnc2gene.pkl", required=True, help="File to save the generated hgnc to gene mapping\n")
    parser.add_argument("--string_links", type=str, dest="string_links", metavar="9606.protein.links.detailed.v11.0.txt.gz", required=True, help="File with the protein relations found in the STRING database\n")
    parser.add_argument("--string_mapping", type=str, dest="string_mapping", metavar="human.name_2_string.tsv.gz", required=True, help="File with mapping of string ids to gene symbols\n")
    parser.add_argument("--gene_interacting", type=str, dest="gene_interacting", metavar="gene2interacting.pkl", required=True, help="File to save the generated gene to interacting genes mapping\n")
    args = parser.parse_args()

    generate_gene2hpo_dict(args.gene_phenotype, args.gene_hpo)
    extract_hpo_graph_edges(args.hpo_ontology, args.hpo_edges)
    generate_hpo_graph(args.hpo_counts, args.hpo_edges, args.hpo_graph)
    create_gene2hgnc_mapping(args.hgnc_symbols, args.hgnc_gene)
    create_gene2interacting_mapping(args.string_mapping, args.string_links, args.gene_interacting)
