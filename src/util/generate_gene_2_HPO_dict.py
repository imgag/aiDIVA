#get mapping gene -> HPOs
#download from HPO charite ALL_FREQ gene to phenotype
#wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
import cPickle as pickle
import sys

gene_2_HPO = dict()

with open(sys.argv[1]) as rd:
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

fname=sys.argv[2]
pickle.dump(gene_2_HPO, open(fname,'wb'))
