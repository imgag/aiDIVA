## how we measure the similarity between two lists w/ IC per each node
## we have a DAG strucutre
## goal is for each Gene !! output a 'semantic distance'
#  based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2756558/ [but different]
#  with this two equal nodes will have distance '0'
#  maximum distance is -2log(1/tot) ~~ 25
import networkx as nx
import pickle
import numpy as np
import math
import random

def calc_me(DG, a, b, PW =False):
    #actual calculation of IC distance
    #return IC(a) + IC(b) -2*IC(MICA)
    # MICA = Max IC Ancestor
        
        
    if  any(x not in DG.nodes()for x in [a,b]):
        #means one key is not in the DG nodes,
        # it can happen so we need to be safe
        #return max possible value
        return 2*(max([d['IC'] for n,d in DG.nodes_iter(data=True)]))
    #check for obsolete nodes
    #substitute by the replacement if obsolete
    a = DG.nodes[a].get('replaced_by',a)
    b = DG.nodes[b].get('replaced_by',b)
    
    
    if  any(x not in DG.nodes()for x in [a,b]):
        #means one key is not in the DG nodes,
        # it can happen so we need to be safe
        #return max possible value
        return 2*(max([d['IC'] for n,d in DG.nodes_iter(data=True)]))

    
    if a==b :
        return 0.0 
        
    # 
    # IC_a = DG.node[a]['IC']
    # IC_b = DG.node[b]['IC']
    #     
    # ancestors_a = list(nx.ancestors(DG,a))
    # ancestors_b = list(nx.ancestors(DG,b))
    # 
    # ancestors_a.append(a)
    # ancestors_b.append(b)
    # 
    # common_ancestors = list(set(ancestors_a) & set(ancestors_b))
    # ancestors_val = [DG.node[x]['IC'] for x in common_ancestors]
    # 
    # distance = IC_a + IC_b -2.0*max(ancestors_val)
    offset =1000
    distance = nx.shortest_path_length(DG,a,b,weight='dist')%offset
    print(distance)
    return distance

def list_distance(DG,Q,G,Query_distances):
    #idea is :
    # for each query HPO calculate all distances
    # store them in a dict with HPOs as keys
    #   value is the minimum value of distance on the query HPOs
    # So than for the list of genes it's enough to
    # collect the values at columns names
    # and if missing set '1'
    #cover cases where no HPO from Query
    # or no HPO provided, or no HPO
    # associated with the gene
    if 'NONE' in Q or 'NONE' in G:
        return (0,Query_distances)
    if len(Q) <1 or len(G) < 1:
        return (0,Query_distances)
    offset =1000
    if Query_distances == 0:        
        # #build it        
        for k_q in Q:
            if  k_q not in list(DG.nodes()):
                #missing node (obsolete not updated or just wrong value)
               continue
            k_q = DG.nodes[k_q].get('replaced_by',k_q)
            distance =  nx.shortest_path_length(DG,k_q,weight='dist')
            if Query_distances ==0:
                  Query_distances = {key: float(value)%offset for (key, value) in distance.items()}
                  print('calc whole dist')
            else:
                for k in Query_distances.keys():
                 try:
                    Query_distances[k] = min([Query_distances[k] , float(distance[k])%offset] )
                 except:
                     Query_distances[k] = float(Query_distances[k])%offset
        if Query_distances == 0:
          #can happen when the original list has no updated HPO or wrong values
          return (0,0)
        Query_distances['maxval']=2*(max([d['IC'] for n,d in DG.nodes(data=True)]))
        # Query_distances['maxval']=2*(max([d['IC'] for n,d in DG.nodes_iter(data=True)]))
    #now I have the query distances value
    # map the genes HPO and extract values.
    # missing one : print it and add it to the db
    #results = []
    maxval = Query_distances['maxval']
    results = [Query_distances.get(q_g,maxval) for q_g in G] 
    #for q_g in G:
    #    q_g = DG.node[q_g].get('replaced_by',q_g)
    #    results.append(Query_distances.get(q_g,2*(max([d['IC'] for n,d in DG.nodes_iter(data=True)]))))
    final_value = np.mean(results)/maxval
    if final_value > 1:
     final_value = 1 #borderline cases whiere go up an down to get to the other node
    return (1-final_value,Query_distances)   

def calc_distance(DG,query,gene,Query_distances=0):
    ### DEPRECATED
    ## Distance (Query, Gene)
    ##
    ## Query = HPO list from user
    ## Gene  = HPO associated to each gene
    
    #asymmetric one
    if len(query)*len(gene) ==0:
        #one of the lists is empty at least
        return 0
    
    #avg [ sum_{t_i \in Q}  min_{t_2 \in G} ( IC(t_1) + IC(t_2) - 2*IC(MICA(t_1,t_2)  ) )       ]  
    #graph contains IC
    
    distances = []
    distances =[ float(min([calc_me(DG,qg,x) for x in gene])) for qg in query]
        
    
    final_value = np.mean(distances)/(2*(max([d['IC'] for n,d in DG.nodes_iter(data=True)])))
    #print distances
    #the division is to ensure a maximum to 1
    #print final_value
    return (1-final_value)


def check_qualtiy(DG):
    #find if all ancestors have IC <= sons
    # if not, why :
    for node in DG:
        ancestors = nx.ancestors(DG,node)
        ancestors_val = [DG.node[x]['IC']  - DG.node[node]['IC'] for x in ancestors]
        problematic = [i for i, e in enumerate(ancestors_val) if e > 0]
        for i in problematic:
            print(node)
            print(list(ancestors)[i])
            print(ancestors_val[i])
    return None


def get_DG_edges(HPO, outfile):
    #This one generates a dict file to generate edges of the HPO graph
    #download data
    #wget https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
    #then call this python ... hp.obo myHPO_edges.pk
    import sys
    import pickle
    listfile = HPO
    
    out_HPO = dict()
    replacements =[]
    alternatives =[]
    token = False
    obsolete =False
    
    with open(listfile) as rd:
        
        for line in rd:
            if line.startswith('id: HP:'):
                if token and not obsolete:
                    out_HPO[name]=parents
                    if repl !='':
                        replacements.append((name,repl))
                        
                token=True
                name = line.strip().split('id: ')[1]
                parents = []
                repl =''
                obsolete =False
                
            elif line.startswith('is_a:'):
                parents.append(line.strip().split('is_a: ')[1].split(' !')[0])
                
            elif line.startswith('replaced_by:'):
                #add a field to say it's replaced
                repl = line.strip().split('replaced_by: ')[1]
                obsolete =False #means we can backtrack it 
            elif line.startswith('is_obsolete:'):
                obsolete =True
            elif line.startswith('alt_id:'):
                #add alternative nodes, will be later added with
                # replacement field for the most common one
                alt = line.strip().split('alt_id: ')[1]
                alternatives.append((name,alt))
            elif line.startswith('consider:'):
                #add alternative nodes, will be later added with
                # replacement field for the most common one
                alt = line.strip().split('consider: ')[1]
                alternatives.append((alt,name))
                obsolete =False #means we can backtrack it
                
                
    out_HPO[name] = parents
    out_HPO['replacements'] = replacements
    out_HPO['alternatives'] = alternatives

    pickle.dump( out_HPO, open( outfile,'wb'))
    

def generate_HPO_graph(edges_file,counts,output):
 offset =1000 #penalization for the distance       
 #usage: python me edges.pk ontology.txt graph.pk
 # counts as wget wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
 # awk -F '\t'  '{print $5}' < phenotype_annotation.tab | sort  |uniq -c | awk '{print $2 "\t" $1}' > HPO_counts.txt
 #idea is a graph with attribute the IC value per node
 #  calculated 
  
 # generate graph with counts:
 counts_d=dict()
 tot = 0
 with open(counts) as rd:
     for line in rd:
         ff=line.strip().split('\t')
         counts_d[ff[0]] = int(ff[1])
         tot += int(ff[1])

 print(tot)

 # load dict with edges
 edges =pickle.load(open(edges_file,'rb'))
 print(len(edges.keys()))
 
 #get replacements of obsolete nodes
 replacements = dict(edges.get('replacements',[]))
 tmpval = edges.pop('replacements',None)

 #let's build a graph
 DG = nx.DiGraph()
 
 #populate with alternatives
 #mark alternatives as replaced, it's the same for us.
 alternatives = edges.get('alternatives',[])
 tmpval = edges.pop('alternatives',None)
 
 # DG.add_edges_from([(1,2)])
 for k in edges.keys():
     DG.add_node(k)
     DG.nodes[k]['count']=0.0
     ancestors = [(x,k) for x in edges[k]]
     DG.add_edges_from(ancestors)
     if k in replacements.keys():
        DG.nodes[k]['replaced_by']=replacements[k]
        DG.nodes[k]['IC'] = -math.log(1.0/tot)
 
 #nx.set_node_attributes(DG, 0,'count',)

 print('edges')
 print(DG.number_of_edges())
 print('nodes')
 print(DG.number_of_nodes())

 for k in DG.nodes():
     DG.nodes[k]['count']=0.0
 
 #populate with raw counts
 for k in counts_d.keys():
     DG.nodes[k]['count'] = counts_d[k]
  
 DG.nodes(data='count')
 #now fill it with the actual value.
 for k in edges.keys():
     desc = nx.descendants(DG,k)
     count = DG.nodes[k]['count']
     for i in desc:
         count += DG.nodes[i]['count']
         
     if count >0 :
         DG.nodes[k]['IC'] =  -math.log(float(count)/tot)
     else :
         DG.nodes[k]['IC'] = -math.log(1.0/tot) #missing nodes, set as rare as possible
         #print k
         #print DG.node[k]
 

 
 # add edges weight
 for a,b in DG.edges():
    DG[a][b]['dist']=offset+abs(DG.nodes[a]['IC'] - DG.nodes[b]['IC'])
 
 
 #alternatives fill in IC and count
 for node,k in alternatives:
    DG.add_node(k)
    DG.nodes[k]['count']=0.0
    DG.nodes[k]['replaced_by']=node
    DG.nodes[k]['IC'] = DG.nodes[node]['IC']
 #count is the IC of the node then : IC = information content
 
 G = DG.to_undirected()
 DG= G
 pickle.dump(DG,open(output,'wb'))
 return None

def generate_gene_2_HPO_dict(HPO_info,outfile):
    #get mapping gene -> HPOs
    #download from HPO charite ALL_FREQ gene to phenotype
    #wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
    gene_2_HPO = dict()

    with open(HPO_info) as rd:
        for line in rd:
            if line.startswith('#'):pass
            else:
                ff = line.strip().split('\t')
                #format #Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-Name<tab>HPO-Term-ID
                key = ff[1]
                HPO =  ff[-1]
                to_add = gene_2_HPO.get(key,[])
                to_add.append(HPO)
                to_add= list(set(to_add))
                gene_2_HPO[key] = to_add
    
    
    pickle.dump(gene_2_HPO, open(outfile,'wb'))
    return None
     
    
def extract_HPO_related_to_gene(gene_2_HPO,gene):
    #gene_2_HPO : dict with [gene] --- HPO_list
    if type(gene_2_HPO) is dict:
            gene_2_HPO_dict = gene_2_HPO
    else:
        gene_2_HPO_dict = pickle.load(open(gene_2_HPO,'rb'))
    outlist =  gene_2_HPO_dict.get(gene,[])  
    return outlist
    
    
def alter_HPO_list(DG,HPO):
    #way to get a list of HPO
    # for each one of these you can
    #   - keep it
    #   - choose an ancestor
    #   - choose a descendant
    #   - remove it
    #   - choose a HPO unrelated
    # all with same  priority
    out_list =[]
    toadd =''
    for hpo in HPO:
        if 'NONE' == hpo :
            out_list = []
            break
        #check replacement
        hpo = DG.nodes[hpo].get('replaced_by',hpo)
        p_val = random.uniform(0,4)
        if p_val < 1:
            #keep it
            out_list.append(hpo)
            continue
        elif p_val < 2:
            #ancestor
            ancestors = list(nx.ancestors(DG,hpo))
            if len(ancestors) >0:
                toadd=random.choice(ancestors)
                out_list.append(toadd)
            continue
        elif p_val < 3:
            #descendants, if none, nothing
            desc = list(nx.descendants(DG,hpo))
            if len(desc) >0:
                toadd=random.choice(desc)
                out_list.append(toadd)
            continue
            #remove it
            
        else:
            ancestors = nx.ancestors(DG,hpo)
            desc = nx.descendants(DG,hpo)
            remaining = list(set(DG.node.keys()) - (ancestors |desc |set(hpo)))
            if len(remaining) >0:
                toadd=random.choice(remaining)
                out_list.append(toadd)

    if len(out_list) <1:
        out_list=['NONE']
    return out_list


def attempt_graph_populate_dist(DG,offset=1000):
    #attempt to use the directed graph to build a new graph with
    # node1 --> ancestor [dist = length]
    # node1 --> descendant [dist =length]
    # so ideally then we can use that to search for shortest path length
    # and get same value as calc_me
    
    GG=nx.Graph()

    for root_id in DG.nodes():
        from_root = (nx.shortest_path_length(DG,root_id,weight='dist'))
        for k,v in from_root.items():
            GG.add_node(k)
            GG.nodes[k]['IC '] = DG.nodes[k]['IC']
            links  = [(root_id ,k)]
            GG.add_edges_from(links)
            GG[root_id][k]['dist']=offset+abs(DG.nodes[root_id]['IC'] - DG.nodes[k]['IC'])
            
            
    ## add replaced!
    replaced_nodes =[x  for x in GG.nodes(data=True)  if 'replaced_by'in x[1].keys()]
    for node_info in replaced_nodes:
        k = node_info[0]
        node_dict = node_info[1]
        GG.add_node(k)
        GG.nodes[k]['IC '] = node_dict['IC']
        GG.nodes[k]['replaced_by '] = node_dict['replaced_by']
    
    return GG
            
def calc_pairwise(DG,outfile):
    #generates a file too big for the moment
    count =0
    #select all keys
    all_dists = dict()
    #remove all replaced_by
    offset=1000
    GG = attempt_graph_populate_dist(DG,offset)
    kk = [x  for x in DG.node.keys() if 'replaced_by' not in DG.node[x].keys()]
    DG = GG    
    kk_y = kk
    #for all keys
    for key_x in kk:
        print("%s %s" % (key_x, str(count)))
        count +=1
        #calc distance
        dists = nx.shortest_path_length(GG,key_x,weight='dist')
        #pop k from key_y
        #Store them in a dict key_x,key_y = [val]
        kk_y.pop(0)
        tmp_keys  =[':'.join([key_x,y]) for y in kk_y]
        tmp_vals  =[dists[y]%offset  for y in kk_y]
        tmp_dict  = dict(zip(tmp_keys, tmp_vals))
        all_dists.update(tmp_dict)
        
        #if keyx == keyy dont store it: dist =0
        #another time
    pickle.dump(all_dists,open(outfile,'wb'))      
        
        
    return None






