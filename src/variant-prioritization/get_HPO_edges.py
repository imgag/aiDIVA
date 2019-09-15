#download data
#wget https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
#then call this python ... hp.obo myHPO_edges.pk
import sys
import pickle
listfile = sys.argv[1]

out_HPO = dict()
token = False
with open(listfile) as rd:
    
    for line in rd:
        if line.startswith('id: HP:'):
            if token:
                out_HPO[name]=parents
            token=True
            name = line.strip().split('id: ')[1]
            parents = []
            
        elif line.startswith('is_a:'):
            parents.append(line.strip().split('is_a: ')[1].split(' !')[0])
            
            
out_HPO[name]=parents
pickle.dump( out_HPO, open( sys.argv[2],'wb'))