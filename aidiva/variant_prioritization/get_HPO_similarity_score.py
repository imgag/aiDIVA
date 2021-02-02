## how we measure the similarity between two lists w/ IC per each node
## we have a DAG strucutre
## goal is for each Gene !! output a "semantic distance"
#  based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2756558/ [but different]
#  with this two equal nodes will have distance "0"
#  maximum distance is -2log(1/tot) ~~ 25
import networkx as nx
import pickle
import numpy as np


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
def list_distance(HPO_graph, HPO_query, gene_HPO_list, HPO_query_distances):
    if "NONE" in HPO_query or "NONE" in gene_HPO_list:
        return 0
    
    if len(HPO_query) < 1 or len(gene_HPO_list) < 1:
        return 0
    
    # map the genes HPO and extract values.
    maxval = HPO_query_distances["maxval"]
    results = [HPO_query_distances.get(gene_HPO, maxval) for gene_HPO in gene_HPO_list]
    final_value = np.mean(results) / maxval

    if final_value > 1:
        final_value = 1 # borderline cases where go up an down to get to the other node

    return 1 - final_value


def precompute_query_distances(HPO_graph, HPO_query):
    HPO_query_distances = dict()
    offset = 1000
    for hpo_term in HPO_query:
        if hpo_term not in list(HPO_graph.nodes()):
            # missing node (obsolete not updated or just wrong value)
            print("INFO: %s not in HPO graph!" % hpo_term)
            continue

        if str(nx.__version__).startswith("1."):
            hpo_term = HPO_graph.node[hpo_term].get("replaced_by", hpo_term)
        elif str(nx.__version__).startswith("2."):
            hpo_term = HPO_graph.nodes[hpo_term].get("replaced_by", hpo_term)
        else:
            print("ERROR: There seems to be a problem with your installation of NetworkX, make sure that you have either v1 or v2 installed!")

        ## TODO: compute shortest path lengths for all nodes not only for hpo_term
        computed_distances =  nx.shortest_path_length(HPO_graph, hpo_term, weight="dist")
        
        if not HPO_query_distances:
                HPO_query_distances = {hpo_id: float(hpo_dist) % offset for (hpo_id, hpo_dist) in computed_distances.items()}
                print("calc whole dist")
        else:
            for hpo_id in computed_distances.keys():
                if hpo_id in HPO_query_distances.keys():
                    HPO_query_distances[hpo_id] = min([HPO_query_distances[hpo_id] , float(computed_distances[hpo_id]) % offset])
                else:
                    HPO_query_distances[hpo_id] = float(computed_distances[hpo_id]) % offset

    
    # IC stored as count
    HPO_query_distances["maxval"] = 2 * (max([node_data["IC"] for node, node_data in HPO_graph.nodes(data=True)]))

    # now I have the query distances value
    return HPO_query_distances
