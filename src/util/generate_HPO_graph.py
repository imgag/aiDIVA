#https://networkx.github.io/documentation/stable/tutorial.html#directed-graphs
import networkx as nx
import cPickle as pickle
import sys
import math



def check_qualtiy(DG):
    #find if all ancestors have IC <= sons
    # if not, why :
    for node in DG:
        ancestors = nx.ancestors(DG,node)
        ancestors_val = [DG.node[x]['count']  - DG.node[node]['count'] for x in ancestors]
        problematic = [i for i, e in enumerate(ancestors_val) if e > 0]
        for i in problematic:
            print node
            print list(ancestors)[i]
            print ancestors_val[i]




#usage: python me edges.pk ontology.txt graph.pk
# counts as wget wget http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab
# awk -F '\t'  '{print $5}' < phenotype_annotation.tab | sort  |uniq -c | awk '{print $2 "\t" $1}' > HPO_counts.txt
#idea is a graph with attribute the IC value per node
#  calculated
counts = sys.argv[2]
output = sys.argv[3]

# generate graph with counts:
counts_d=dict()
tot = 0
with open(counts) as rd:
    for line in rd:
        ff=line.strip().split('\t')
        counts_d[ff[0]] = int(ff[1])
        tot += int(ff[1])

print tot

# load dict with edges
edges =pickle.load(open(sys.argv[1],'rb'))
print( len(edges.keys()))

#let's build a graph
DG = nx.DiGraph()
# DG.add_edges_from([(1,2)])
for k in edges.keys():
    DG.add_node(k)
    DG.node[k]['count']=0.0
    ancestors = [(x,k) for x in edges[k]]
    DG.add_edges_from(ancestors)

#nx.set_node_attributes(DG, 0,'count',)

print 'edges'
print DG.number_of_edges()
print 'nodes'
print DG.number_of_nodes()

for k in DG.nodes():
    DG.node[k]['count']=0.0

#populate with raw counts
for k in counts_d.keys():
    DG.node[k]['count'] = counts_d[k]



DG.nodes(data='count')
#now fill it with the actual value.
for k in edges.keys():
    desc = nx.descendants(DG,k)
    count = DG.node[k]['count']
    for i in desc:
        count += DG.node[i]['count']

    if count >0 :
        DG.node[k]['count'] =  -math.log(float(count)/tot)
    else :
        DG.node[k]['count'] = -math.log(1.0/tot) #missing nodes, set as rare as possible
        #print k
        #print DG.node[k]

#count is the IC of the node then
pickle.dump(DG,open('DG_temp.pk','wb'))
     
