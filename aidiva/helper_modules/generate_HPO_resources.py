import argparse
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
                key = ff[3]
                HPO = ff[0]
                to_add = gene_2_HPO.get(key,[])
                to_add.append(HPO)
                to_add = list(set(to_add))
                gene_2_HPO[key] = to_add

    pickle.dump(gene_2_HPO, open(gene2hpo_dict,"wb"))
    print("Gene to HPO mapping successfully generated and saved as %s" % (gene2hpo_dict))


# get mapping HPO -> genes
# download from HPO charite phenotype to gene
# wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/phenotype_to_genes.txt
## TODO: remove this method (is obsolete)
def generate_hpo2gene_dict(gene2phenotype_list, hpo2gene_dict):
    print("Generate HPO to gene mapping...")
    HPO_association = dict()
    with open(gene2phenotype_list) as rd:
        for line in rd:
            if line.startswith("#"):
                pass
            else:
                ff = line.strip().split("\t")
                key = ff[0]
                value = ff[3]
                try:
                    HPO_association[key].append(value)
                except:
                    HPO_association[key] = [value]

    print("HPO to gene mapping successfully generated and saved as %s" % (hpo2gene_dict))
    pickle.dump(HPO_association, open(hpo2gene_dict, "wb"))


# download data
# wget https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
def extract_hpo_graph_edges(hpo_ontology, hpo_edges_file):
    print("Extract HPO edges...")
    out_HPO = dict()
    replacements = []
    token = False
    obsolete = False

    with open(hpo_ontology) as ontology:
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

    out_HPO[name] = parents
    out_HPO['replacements'] = replacements
    #out_HPO['alternatives'] = alternatives
    pickle.dump(out_HPO, open(hpo_edges_file, "wb"))
    print("HPO edges successfully extracted and saved as %s" % (hpo_edges_file))


def check_graph_qualtiy(DG):
    # find if all ancestors have IC <= sons
    # if not, why:
    for node in DG:
        ancestors = nx.ancestors(DG, node)
        if str(nx.__version__).startswith("1."):
            ancestors_val = [DG.node[x]["count"]  - DG.node[node]["count"] for x in ancestors]
        elif str(nx.__version__).startswith("2."):
            ancestors_val = [DG.nodes[x]["count"]  - DG.nodes[node]["count"] for x in ancestors]
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        problematic = [i for i, e in enumerate(ancestors_val) if e > 0]
        for i in problematic:
            print(node)
            print(list(ancestors)[i])
            print(ancestors_val[i])


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

    print(tot)

    # load dict with edges
    edges = pickle.load(open(hpo_edges_file, "rb"))
    print(len(edges.keys()))

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

    print("edges")
    print(hpo_graph.number_of_edges())
    print("nodes")
    print(hpo_graph.number_of_nodes())

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
                hpo_graph.node[node]["IC"] = -math.log(1.0 / tot) #missing nodes, set as rare as possible
            elif str(nx.__version__).startswith("2."):
                hpo_graph.nodes[node]["IC"] = -math.log(1.0 / tot) #missing nodes, set as rare as possible
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")
            #print(node)
            #print(hpo_graph.node[node])

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

    #pickle.dump(hpo_graph,open(output,"wb"))
    if str(nx.__version__).startswith("1."):
        pickle.dump([hpo_graph.nodes(data=True), hpo_graph.edges(data=True)], open(output, "wb"))
        print("NOTE: You pickled the graph with NetworkX v1 the pickled graph it is also upwards compatible with v2!")
    elif str(nx.__version__).startswith("2."):
        pickle.dump([hpo_graph.nodes(data=True), hpo_graph.edges(data=True)], open(output, "wb"))
        print("NOTE: You pickled the graph with NetworkX v2 it is not backwards compatible!")
    else:
        print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    print("HPO graph successfully generated and saved as %s" % (hpo_graph_file))


if __name__=="__main__":
    parser = argparse.ArgumentParser("Script to generate the HPO resources needed in the prioritization step of AIdiva")
    parser.add_argument("--hpo_ontology", type=str, dest="hpo_ontology", metavar="hp.obo", required=True, help="File containing the HPO ontology\n")
    parser.add_argument("--gene_phenotype", type=str, dest="gene_phenotype", metavar="phenotype_to_genes.txt", required=True, help="File that contains information about the genes and the associated phenotypes\n")
    parser.add_argument("--gene_hpo", type=str, dest="gene_hpo", metavar="gene2hpo.pkl", required=True, help="File to save the generated gene2hpo_dict\n")
    parser.add_argument("--hpo_edges", type=str, dest="hpo_edges", metavar="hpo_edges.pkl", required=True, help="File where the extracted hpo edges are stored\n")
    parser.add_argument("--hpo_counts", type=str, dest="hpo_counts", metavar="HPO_counts.txt", required=True, help="File containing the hpo counts needed for the hpo graph construction\n")
    parser.add_argument("--hpo_graph", type=str, dest="hpo_graph", metavar="hpo_graph.pkl", required=True, help="File to save the generated hpo graph\n")
    args = parser.parse_args()

    generate_gene2hpo_dict(args.gene_phenotype, args.gene_hpo)
    extract_hpo_graph_edges(args.hpo_ontology, args.hpo_edges)
    generate_hpo_graph(args.hpo_counts, args.hpo_edges, args.hpo_graph)
