
import numpy as np
import logging


logger = logging.getLogger(__name__)


def find_mica_IC(HPO_term_a, HPO_term_b, ic_per_nodes, node_ancestor_mapping):
    hpo_nodes_shared_ancestors = node_ancestor_mapping[HPO_term_a].intersection(node_ancestor_mapping[HPO_term_b])
    mica_ic = max([float(ic_per_nodes[node]) for node in hpo_nodes_shared_ancestors], default=0.0)
    
    return mica_ic


def compute_similarity_between_nodes(hpo_term_a, hpo_term_b, ic_per_nodes, node_ancestor_mapping, similarity_measure="Lin"):
    node_a_ic = float(ic_per_nodes[hpo_term_a])
    node_b_ic = float(ic_per_nodes[hpo_term_b])
    mica_ic = find_mica_IC(hpo_term_a, hpo_term_b, ic_per_nodes, node_ancestor_mapping)
    
    if similarity_measure == "Lin":
        if mica_ic == 0.0 or (node_a_ic == 0.0 and node_b_ic == 0.0):
            nodes_similarity = 0.0
        else:
            nodes_similarity = (2.0 * mica_ic) / (node_a_ic + node_b_ic)
        
    elif similarity_measure == "Resnik":
        nodes_similarity = mica_ic
        
    elif similarity_measure == "Jiang-Conrath":
        nodes_similarity = (1 - (node_a_ic + node_b_ic - 2.0 * mica_ic))
        
    elif similarity_measure == "Relevance":
        #nodes_similarity = ((2.0 * mica_ic) / (node_a_ic + node_b_ic)) * (1 - p(mica))
        pass

    elif similarity_measure == "Information-Coefficient":
        nodes_similarity = ((2.0 * mica_ic) / (node_a_ic + node_b_ic)) * (1 - (1 / (1 + mica_ic)))

    elif similarity_measure == "Graph-IC":
        pass

    elif similarity_measure == "Wang":
        pass
        
    else:
        logger.error("An error occured while determining the similarity measure!")
        
    return nodes_similarity


def calculate_hpo_set_similarity(hpo_graph, hpo_term_set_a, hpo_term_set_b, ic_per_nodes, node_ancestor_mapping, hpo_replacement_information):
    checked_term_set_a = []
    checked_term_set_b = []
    
    # alternatives and/or considerations are ignored for now
    alternatives = hpo_replacement_information["alternatives"]
    considerations = hpo_replacement_information["considerations"]
    replacements = hpo_replacement_information["replacements"]
    
    for term_a in hpo_term_set_a:
        if term_a not in hpo_graph:
            if term_a in replacements.keys():
                checked_term_set_a.append(replacements[term_a])
                logger.debug(f"{term_a} (sample) not in HPO graph! Replacement ({replacements[term_a]}) found will use this term instead!")
            elif term_a in alternatives.keys():
                #checked_term_set_a.extend(alternatives[term_a])
                logger.debug(f"{term_a} (sample) not in HPO graph! Alternatives ({alternatives[term_a]}) found! HPO term will be skipped!")
            elif term_a in considerations.keys():
                #checked_term_set_a.extend(considerations[term_a])
                logger.debug(f"{term_a} (sample) not in HPO graph! Considerations ({considerations[term_a]}) found! HPO term will be skipped!")
            else:
                logger.debug(f"{term_a} (sample) not in HPO graph! HPO term will be skipped!")
        else:
            checked_term_set_a.append(term_a)
    
    for term_b in hpo_term_set_b:
        if term_b not in hpo_graph:
            if term_b in replacements.keys():
                checked_term_set_b.append(replacements[term_b])
                logger.debug(f"{term_b} (gene) not in HPO graph! Replacement ({replacements[term_b]}) found will use this term instead!")
            elif term_b in considerations.keys():
                #checked_term_set_b.extend(considerations[term_b])
                logger.debug(f"{term_b} (gene) not in HPO graph! Considerations ({considerations[term_b]}) found! HPO term will be skipped!")
            elif term_b in alternatives.keys():
                #checked_term_set_b.extend(alternatives[term_b])
                logger.debug(f"{term_b} (gene) not in HPO graph! Alternatives ({alternatives[term_b]}) found! HPO term will be skipped!")
            else:
                logger.debug(f"{term_b} (gene) not in HPO graph! HPO term will be skipped!")
        else:
            checked_term_set_b.append(term_b)
    
    if checked_term_set_a and checked_term_set_b:
        similarities_a_to_b = [max([compute_similarity_between_nodes(term_a, term_b, ic_per_nodes, node_ancestor_mapping) for term_b in checked_term_set_b], default=0.0) for term_a in checked_term_set_a]
        #similarities_b_to_a = [max([compute_similarity_between_nodes(term_b, term_a, ic_per_nodes, node_ancestor_mapping) for term_a in checked_term_set_a], default=0.0) for term_b in checked_term_set_b]
        
        set_a_to_b_similarity = np.median(similarities_a_to_b)
        #set_b_to_a_similarity = np.median(similarities_b_to_a)

        #if set_a_to_b_similarity != 0.0 or set_b_to_a_similarity != 0.0:
        #    hpo_set_similarity = (set_a_to_b_similarity + set_b_to_a_similarity) / 2
        #else:
        #    hpo_set_similarity = 0.0
        
        return set_a_to_b_similarity

    else:
        logger.debug(f"Sample HPO set ({hpo_term_set_a}) and/or Gene HPO set ({hpo_term_set_b}) was empty after checking if the terms are part of the used HPO graph! Maybe no supported term was in the set! Similarity set to 0.0!")
        
        return 0.0
