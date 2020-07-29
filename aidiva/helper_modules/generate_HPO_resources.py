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
# then call this python ... hp.obo myHPO_edges.pk
def extract_hpo_graph_edges(hpo_ontology, hpo_edges_file):
    print("Extract HPO edges...")
    out_HPO = dict()
    token = False

    with open(hpo_ontology) as rd:
        for line in rd:
            if line.startswith("id: HP:"):
                if token:
                    out_HPO[name] = parents
                token = True
                name = line.strip().split("id: ")[1]
                parents = []
            elif line.startswith("is_a:"):
                parents.append(line.strip().split("is_a: ")[1].split(" !")[0])

    out_HPO[name] = parents
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


# usage: python me edges.pk ontology.txt graph.pk
# counts as wget wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
# awk -F '\t'  '{print $5}' < phenotype_annotation.tab | sort  |uniq -c | awk '{print $2 "\t" $1}' > HPO_counts.txt
def generate_hpo_graph(hpo_counts, hpo_edges_file, hpo_graph_file):
    print("Generate HPO graph...")
    # idea is a graph with attribute the IC value per node
    # calculated
    counts = hpo_counts
    output = hpo_graph_file

    # generate graph with counts:
    counts_d = dict()
    tot = 0
    with open(counts) as rd:
        for line in rd:
            ff = line.strip().split("\t")
            counts_d[ff[0]] = int(ff[1])
            tot += int(ff[1])

    print(tot)

    # load dict with edges
    edges = pickle.load(open(hpo_edges_file,"rb"))
    print(len(edges.keys()))

    # let"s build a graph
    DG = nx.DiGraph()
    # DG.add_edges_from([(1,2)])
    for k in edges.keys():
        DG.add_node(k)
        if str(nx.__version__).startswith("1."):
            DG.node[k]["count"] = 0.0
        elif str(nx.__version__).startswith("2."):
            DG.nodes[k]["count"] = 0.0
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        ancestors = [(x,k) for x in edges[k]]
        DG.add_edges_from(ancestors)

    # nx.set_node_attributes(DG, 0,"count",)

    print("edges")
    print(DG.number_of_edges())
    print("nodes")
    print(DG.number_of_nodes())

    for k in DG.nodes():
        if str(nx.__version__).startswith("1."):
            DG.node[k]["count"] = 0.0
        elif str(nx.__version__).startswith("2."):
            DG.nodes[k]["count"] = 0.0
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    # populate with raw counts
    for k in counts_d.keys():
        if str(nx.__version__).startswith("1."):
            DG.node[k]["count"] = counts_d[k]
        elif str(nx.__version__).startswith("2."):
            DG.nodes[k]["count"] = counts_d[k]
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

    DG.nodes(data="count")
    # now fill it with the actual value.
    for k in edges.keys():
        desc = nx.descendants(DG,k)
        if str(nx.__version__).startswith("1."):
            count = DG.node[k]["count"]
        elif str(nx.__version__).startswith("2."):
            count = DG.nodes[k]["count"]
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        for i in desc:
            if str(nx.__version__).startswith("1."):
                count += DG.node[i]["count"]
            elif str(nx.__version__).startswith("2."):
                count += DG.nodes[i]["count"]
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        if count > 0 :
            if str(nx.__version__).startswith("1."):
                DG.node[k]["count"] =  -math.log(float(count) / tot)
            elif str(nx.__version__).startswith("2."):
                DG.nodes[k]["count"] =  -math.log(float(count) / tot)
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")
        else :
            if str(nx.__version__).startswith("1."):
                DG.node[k]["count"] = -math.log(1.0 / tot) #missing nodes, set as rare as possible
            elif str(nx.__version__).startswith("2."):
                DG.nodes[k]["count"] = -math.log(1.0 / tot) #missing nodes, set as rare as possible
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")
            # print k
            # print DG.node[k]

    # count is the IC of the node then
    #pickle.dump(DG,open(output,"wb"))
    if str(nx.__version__).startswith("1."):
        pickle.dump([DG.nodes(data=True), DG.edges(data=True)], open(output, "wb"))
        print("NOTE: You pickled the graph with NetworkX v1 the pickled graph is also importable with v2!")
    elif str(nx.__version__).startswith("2."):
        pickle.dump([DG.nodes(data=True), DG.edges(data=True)], open(output, "wb"))
        print("NOTE: You pickled the graph with NetworkX v2 it is not backwards compatible!")
    else:
        print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")
    print("HPO graph successfully generated and saved as %s" % (hpo_graph_file))


if __name__=="__main__":
    parser = argparse.ArgumentParser("Script to generate the HPO resources needed in the prioritization step of AIdiva")
    parser.add_argument("--hpo_ontology", type=str, dest="hpo_ontology", metavar="hp.obo", required=True, help="File containing the HPO ontology\n")
    parser.add_argument("--gene_phenotype", type=str, dest="gene_phenotype", metavar="phenotype_to_genes.txt", required=True, help="File that contains information about the genes and the associated phenotypes\n")
    parser.add_argument("--gene_hpo", type=str, dest="gene_hpo", metavar="gene2hpo.pkl", required=True, help="File to save the generated gene2hpo_dict\n")
    parser.add_argument("--hpo_gene", type=str, dest="hpo_gene", metavar="hpo2gene.pkl", required=True, help="File to save the generated hpo2gene_dict\n")
    parser.add_argument("--hpo_edges", type=str, dest="hpo_edges", metavar="hpo_edges.pkl", required=True, help="File where the extracted hpo edges are stored\n")
    parser.add_argument("--hpo_counts", type=str, dest="hpo_counts", metavar="HPO_counts.txt", required=True, help="File containing the hpo counts needed for the hpo graph construction\n")
    parser.add_argument("--hpo_graph", type=str, dest="hpo_graph", metavar="hpo_graph.pkl", required=True, help="File to save the generated hpo graph\n")
    args = parser.parse_args()

    generate_gene2hpo_dict(args.gene_phenotype, args.gene_hpo)
    generate_hpo2gene_dict(args.gene_phenotype, args.hpo_gene)
    extract_hpo_graph_edges(args.hpo_ontology, args.hpo_edges)
    generate_hpo_graph(args.hpo_counts, args.hpo_edges, args.hpo_graph)
