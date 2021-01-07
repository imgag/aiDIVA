## how we measure the similarity between two lists w/ IC per each node
## we have a DAG strucutre
## goal is for each Gene !! output a "semantic distance"
#  based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2756558/ [but different]
#  with this two equal nodes will have distance "0"
#  maximum distance is -2log(1/tot) ~~ 25
import networkx as nx
import pickle
import numpy as np


def list_distance(DG, Q, G, Query_distances):
    #idea is :
    # for each query HPO calculate all distances
    # store them in a dict with HPOs as keys
    # value is the minimum value of distance on the query HPOs
    # So than for the list of genes it"s enough to
    # collect the values at columns names
    # and if missing set "1"
    #cover cases where no HPO from Query
    # or no HPO provided, or no HPO
    # associated with the gene
    if "NONE" in Q or "NONE" in G:
        return (0, Query_distances)
    if len(Q) < 1 or len(G) < 1:
        return (0, Query_distances)
    offset = 1000
    if Query_distances == 0:
        for k_q in Q:
            if k_q not in list(DG.nodes()):
                # missing node (obsolete not updated or just wrong value)
                continue

            if str(nx.__version__).startswith("1."):
                k_q = DG.node[k_q].get("replaced_by", k_q)
            elif str(nx.__version__).startswith("2."):
                k_q = DG.nodes[k_q].get("replaced_by", k_q)
            else:
                print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

            distance =  nx.shortest_path_length(DG, k_q, weight="dist")
            if Query_distances == 0:
                Query_distances = {key: float(value) % offset for (key, value) in distance.items()}
                print("calc whole dist")
            else:
                for k in Query_distances.keys():
                    try:
                        Query_distances[k] = min([Query_distances[k] , float(distance[k]) % offset])
                    except:
                        Query_distances[k] = float(Query_distances[k]) % offset

        if Query_distances == 0:
            # can happen when the original list has no updated HPO or wrong values
            return (0, 0)

        # IC stored as count
        Query_distances["maxval"] = 2 * (max([d["IC"] for n, d in DG.nodes(data=True)]))

    # now I have the query distances value
    # map the genes HPO and extract values.
    # missing one : print it and add it to the db
    maxval = Query_distances["maxval"]
    results = [Query_distances.get(q_g,maxval) for q_g in G]
    final_value = np.mean(results) / maxval

    if final_value > 1:
        final_value = 1 # borderline cases where go up an down to get to the other node

    return (1 - final_value, Query_distances)


## TODO: the following method is obsolete, "list_distance" can be used instead
def precompute_query_distances(DG, Q, Query_distances):
    offset = 1000
    for k_q in Q:
        if k_q not in list(DG.nodes()):
            # missing node (obsolete not updated or just wrong value)
            continue

        if str(nx.__version__).startswith("1."):
            k_q = DG.node[k_q].get("replaced_by", k_q)
        elif str(nx.__version__).startswith("2."):
            k_q = DG.nodes[k_q].get("replaced_by", k_q)
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        ## TODO: compute shortest path lengths for all nodes not only for k_q
        distance =  nx.shortest_path_length(DG, k_q, weight="dist")
        if Query_distances == 0:
            Query_distances = {key: float(value) % offset for (key, value) in distance.items()}
            print("calc whole dist")
        else:
            for k in Query_distances.keys():
                try:
                    Query_distances[k] = min([Query_distances[k] , float(distance[k]) % offset])
                except:
                    Query_distances[k] = float(Query_distances[k]) % offset

    if Query_distances == 0:
        # can happen when the original list has no updated HPO or wrong values
        return 0

    # IC stored as count
    Query_distances["maxval"] = 2 * (max([d["IC"] for n, d in DG.nodes(data=True)]))

    return Query_distances


def extract_HPO_related_to_gene(gene_2_HPO, gene):
    # gene_2_HPO : dict with [gene] --- HPO_list
    if type(gene_2_HPO) is dict:
        gene_2_HPO_dict = gene_2_HPO
    else:
        gene_2_HPO_dict = pickle.load(open(gene_2_HPO, "rb"))
    outlist = gene_2_HPO_dict.get(gene, [])

    return outlist
