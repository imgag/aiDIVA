import argparse
import csv
import pprint
import re
from scipy.stats import poisson
import os
#import MySQLdb
import shutil
import pickle
import sys
from operator import itemgetter

try:
    import xlsxwriter
    import xlrd
    writeXLS = True
except:
    writeXLS = False
    print('No XlS writer')

# Note header:
#[   'Chr',
#    'Position',
#    'Reference',
#    'Alteration',
#    'Function(Refseq)',
#    'Gene(Refseq)',
#    'ExonicFunction(Refseq)',
#    'AminoAcidChange(Refseq)',
#    'Function(Ensembl)',
#    'Gene(Ensembl)',
#    'ExonicFunction(Ensembl)',
#    'AminoAcidChange(Ensembl)',
#    'Function(Known)',
#    'Gene(Known)',
#    'ExonicFunction(Known)',
#    'AminoAcidChange(Known)',
#    'dbsnpIdentifier',
#    'dbSNPfrequency',
#    'EurEVSFrequecy',
#    'AfrEVSFrequecy',
#    'TotalEVSFrequecy',
#    'Eur1000GenomesFrequency',
#    'Afr1000GenomesFrequency',
#    'Amr1000GenomesFrequency',
#    'Asia1000GenomesFrequency',
#    'Total1000GenomesFrequency',
#    'SegMentDup',
#    'PlacentalMammalPhyloP',
#    'PrimatesPhyloP',
#    'VertebratesPhyloP',
#    'PlacentalMammalPhastCons',
#    'PrimatesPhastCons',
#    'VertebratesPhastCons',
#    'Score1GERP++',
#    'Score2GERP++',
#    'SIFTScore',
#    'polyphen2',
#    'MutAss',
#    'Condel',
#    'pos_samples(sampleid>zygosity>Cov>AF)',
#    'neg_samples(sampleid>zygosity>Cov>AF)',
#    'rank']


parser = argparse.ArgumentParser(description = 'filter SNPs according to their family distribution')
    
parser.add_argument('--infile', type=argparse.FileType('r'), dest='infile', required=True, help='comma separated list of SNPs annotated with mutation impact data. [required]')
parser.add_argument('--outfile', type=argparse.FileType('w'), dest='outfile', required=True, help='comma separated list of SNPs annotated with inheritance pattern. [required]')
parser.add_argument('--filteredoutfile', type=argparse.FileType('w'), dest='filteredfile', required=True, help='filtered comma separated list of SNPs annotated with inheritance pattern, only reporting the requested variants. [required]')
parser.add_argument('--family',  type=argparse.FileType('r'), dest='famfile', required=True, help='tab separated list of samples annotated with affection status. [required]')
parser.add_argument('--inheritance', choices=['dominant_denovo', 'dominant_inherited', 'recessive', 'Xlinked', 'compound'], dest='inheritance', required=True, help="""choose a inheritance model [required]
dominant_inherited: used for families
dominant_denovo: apply to novel variants seen in the affected individuals

recessive: detect recessive, homozygous variants (if trio is specified the script will require that all non-affected are heterozygous)
Xlinked: used for X linked recessive variants in trios only
compound: detect compound heterozygous recessive variants
""")
parser.add_argument('--familytype', choices=['trio', 'family'], dest='familytype', required=True, help="choose if the data you provide is a trio or a larger family")
parser.add_argument('--geneexclusion',  type=argparse.FileType('r'), dest='geneexclusion', required=False, help='[Analysis of DNA sequence variants detected by high-throughput sequencing; DOI: 10.1002/humu.22035]. [required]')
parser.add_argument('--HPO_list','--white_list',type=str,dest='white_list',default=None,required=False,help='--HPO_list \t a .txt file with the list of HPO terms describing the disease. It will be used to flag all genes related to the HPO terms. It works with Refseq and UCSC naming.\n')
parser.add_argument('--csvfile', dest='csvfile', required=False, help='csv file with username and user email address. [optional]')

args = parser.parse_args()

def main (args):
    #mailer_path='/home/rrahman/soft/python-mailer/pymailer.py'
    #pp = pprint.PrettyPrinter( indent=4) # DEBUG
    names = list()
    # read the gene exclusion list
    # an empty set will not filter out anything, if gene exclusion list is not provided
    genes2exclude = set()
    genes_known   = set()
    if args.geneexclusion:
        for gene in args.geneexclusion:
            gene = gene.rstrip()
            genes2exclude.add(gene)

    if args.white_list == None:
        args.white_list = 'None'
        related_genes = set()
    else:
        related_genes = set()
        if os.path.isfile(args.white_list):
            script_dir = os.path.dirname(sys.argv[0])
            with open(args.white_list,'r') as w:
                HPO_dict  = pickle.load(open('HPO_gene_assiciation.p','rb'))
                related_genes = list()
                for line in w:
                    HPO_term = line.rstrip('\n')
                    try:
                        related_genes.extend(HPO_dict[HPO_term])
                    except:
                        print('%s not found in database'%(HPO_term))
                    
            genes2exclude = set(genes2exclude) - set(related_genes)
            related_genes= set(related_genes)
        else:
            print('The specified HPO list %s is not a valid file'%(args.white_list))
            print('eDiVA will proceed as without any HPO list')

        
    # read family relationships
    family = dict()
    for line in args.famfile:
        if line.startswith('sample'):
            continue
        line = line.rstrip('\n')
        splitline = line.split('\t')
        family[splitline[0]] = splitline[1]
        #names.append(splitline[0])


    ####
    # example
    # sample	affected
    # VH017	1
    # VH018	0
    # VH019	0
    ####

    # read all data
    alldata = list(csv.reader(args.infile))
   
    header = alldata.pop(0)
    if args.inheritance == 'compound' :
      alldata = sorted(alldata, key=itemgetter(0,1))
    header =[x.replace('#','',1) for x in header]
    #print header
    #raise
    out         = csv.writer(args.outfile)
    outfiltered = csv.writer(args.filteredfile)
    
    ############
    # Filter out SNPs from filtered file where : 
    # kick when variant function : synonymous,unknown => line[ExonicFunction(Refseq)]
    # keep when function should be exonic,splicing => line[Function(Refseq)]
    # kick when segmentdup should be 0 => line[SegMentDup]
    # only consider Refseq annotation for performing the first criteria
    #############

    index_sample      = min([identifycolumns(header,x) for x in family.keys()    ])#+identifycolumns(header, 'ExAC_SAS') #samples(sampleid>zygosity>DPRef>DPAlt>AF)
    index_MAF1k       = identifycolumns(header, 'Total1000GenomesFrequency')
    index_MAFevs      = identifycolumns(header, 'TotalEVSFrequency')
    index_MAF_exac    = identifycolumns(header, 'ExAC_AF')
    index_function    = identifycolumns(header, 'Function(Refseq)')
    index_varfunction = identifycolumns(header, 'ExonicFunction(Refseq)')
    index_segdup      = identifycolumns(header, 'SegMentDup')
    index_gene        = identifycolumns(header, 'Gene(Refseq)')
    index_str         = identifycolumns(header, 'SimpleTandemRepeatLength')
    index_CADD        = identifycolumns(header, 'Cadd2')
    #print index_sample
    sample_annot_size = 6 # sampleid - DP - REF - ALT - AF - GQ
    print(family.keys())
    for i in range(index_sample,index_sample+sample_annot_size*len(family.keys()),sample_annot_size):
        names.append(header[i])

    #print header[index_sample]
    
    ############
    # both the output file headers should be consistent
    #############
    
    header.append('inheritance')
    header.append('filter')
    header.append('Related to HPO')
    #header.extend(['OMIM_name','OMIM_ID','clinical_significance', 'disease_name', 'clinical_review',' access_number'])
    outfiltered.writerow(header)
    out.writerow(header)

    
    # init for compound
    initializer = 0
    #if args.inheritance == 'compound':
    #    compound_gene_storage = []
    #    pp.pprint(line)
    #    pp.pprint([line[index_gene], index_gene])
    #    old_gene   = re.sub('\(.*?\)','',line[index_gene]) # PKHD1L1(NM_177531:exon75:c.12330+1G>A) transformed to PKHD1L1. Also works for ";" separated multiple annotations
    #    pp.pprint(old_gene)
    #    old_gene_set = set(old_gene.split(';')) # try to remove double occurences of the same gene names
    #    pp.pprint(old_gene_set)
    #    new_gene   = re.sub('\(.*?\)','',line[index_gene])
    #    new_gene_set = set(new_gene.split(';'))
    
    # start reading data
    for line in alldata:
         # init for compound
        if args.inheritance == 'compound' and initializer == 0:
            compound_gene_storage = []
            old_gene   = re.sub('\(.*?\)','',line[index_gene]) # PKHD1L1(NM_177531:exon75:c.12330+1G>A) transformed to PKHD1L1. Also works for ";" separated multiple annotations
            old_gene_set = set(old_gene.split(';')) # try to remove double occurences of the same gene names
            new_gene   = re.sub('\(.*?\)','',line[index_gene])
            new_gene_set = set(new_gene.split(';'))
            initializer = 1
        
        # read minor allele frequencies
        try:
            MAF1k      = line[index_MAF1k].replace('NA','0')
            MAFevs     = line[index_MAFevs].replace('NA','0')
            if 'NA' == line[index_MAF_exac]:
                MAFexac ='0'
            else:
                MAFexac= line[index_MAF_exac]
        except:
            #pp.pprint(line)
            pass
        try:
            # avoiding problems with NAs
            MAF    = max(float(MAF1k), float(MAFevs),float(MAFexac))

        except:
            print(('Freq error 1k %s EVS %s ExAC %s')%(MAF1k,MAFevs,MAFexac))
            MAF    = 0
        try :
            CADD = float(line[index_CADD])
        except:
            CADD  = 0 
        
        # read sample names and according zygosity NOW IT's a LIST
        sampledata = line[index_sample:index_sample+sample_annot_size*len(family.keys())]

        if len(sampledata)==0:
            print(line)
            print(len(line))
            raise
        elif sampledata[0] =='0':
            print(len(line))
            print(line)
            raise
        # read simple tandem repeat information
        tandem = line[index_str]
        
        # check, if gene is on the gene exclusion list.
        genenames = set()
        genecolumn   = re.sub('\(.*?\)','',line[index_gene])
        genenames = set(genecolumn.split(';'))
        
    
        if len(related_genes & genenames) >0:
            known= 'yes'
        else:
            known = 'no'
        judgement = int()
        ###
        # look for de novo variants
        ###
        if args.inheritance == 'dominant_denovo':
            
            judgement = denovo(sampledata, family,names)
            # exclude gene, if it is on the exclusion list
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                line.append(known)
                out.writerow(line)
                continue

            # check before all others, if variant locates to simple tandem repeat region
            elif judgement == 1 and not tandem == 'NA':
                line.append('denovo')
                line.append('tandem')
                out.writerow(line)
                continue
            
            # Exclude non damaging CADD variants            
            elif judgement == 1 and CADD<=19 and CADD>0 :
                line.append('dominant')
                line.append('cadd')
                line.append(known)
                out.writerow(line)
                continue
            
            elif judgement == 1 and MAF <= 0.01:
                
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                    if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                        if (line[index_segdup] == '0'):
                            line.append('denovo')
                            line.append('pass')
                            line.append(known)
                            out.writerow(line)
                            outfiltered.writerow(line)
                # e.g. intronic variants fitting the criteria
                        else:
                            line.append(args.inheritance)
                            line.append('filtered')
                            line.append(known)
                            out.writerow(line)
                    else:
                        line.append(args.inheritance)
                        line.append('filtered')
                        line.append(known)
                        line.append(known)
                        out.writerow(line)
                else:
                    line.append(args.inheritance)
                    line.append('filtered')
                    line.append(known)
                    out.writerow(line)
                continue
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.01:
                line.append('denovo')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
            
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
        
        ###
        # look for familial dominant variants. (being tolerant for missing values)
        ###
        elif args.inheritance == 'dominant_inherited':

            judgement = dominant(sampledata, family,names)
            sample_annot_size = len(sampledata)/len(names)
            total_affected=0
            for i in range(0,sample_annot_size*len(names),sample_annot_size):                
                sam = sampledata[i]
                features    = sampledata[i:i+sample_annot_size]#sam.split(':')
                name        = names[i/sample_annot_size]
                total_affected += int(family[name] )
            if total_affected==1:judgement *=0 #denovo mutation excluded       

            # exclude gene, if it is on the exclusion list
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                line.append(known)
                out.writerow(line)
                continue
            
            # check before all others, if variant locates to simple tandem repeat region
            elif judgement == 1 and not tandem == 'NA':
                line.append('dominant')
                line.append('tandem')
                line.append(known)
                out.writerow(line)
                continue
             #Exclude non damaging CADD variants            
            elif judgement == 1 and CADD<=20 and CADD>0 :
                line.append('dominant')
                line.append('cadd')
                line.append(known)
                out.writerow(line)
                continue
            
            
            elif judgement == 1 and MAF <= 0.01:#0.05
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                    if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                        if (line[index_segdup] == '0'):
                            line.append('dominant')
                            line.append('pass')
                            line.append(known)
                            out.writerow(line)
                            outfiltered.writerow(line)
                # e.g. intronic variants fitting the criteria
                        else:
                            line.append(args.inheritance)
                            line.append('filtered')
                            line.append(known)
                            out.writerow(line)
                    else:
                        line.append(args.inheritance)
                        line.append('filtered')
                        line.append(known)
                        out.writerow(line)
                else:
                    line.append(args.inheritance)
                    line.append('filtered')
                    line.append(known)
                    out.writerow(line)
                continue
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.01:#0.05
                line.append('dominant')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
            
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
        
        ###
        # look for recessive variants (be aware of trio and family inheritance)
        ###
        elif args.inheritance == 'recessive':
            judgement = recessive(sampledata, family, args.familytype,names)
            
            # exclude gene, if it is on the exclusion list
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                line.append(known)
                out.writerow(line)
                continue
            
            # check before all others, if variant locates to simple tandem repeat region
            elif judgement == 1 and not tandem == 'NA':
                line.append('recessive')
                line.append('tandem')
                line.append(known)
                out.writerow(line)
                continue
            
            elif judgement == 1 and MAF <= 0.03:
                
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                    if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                        if (line[index_segdup] == '0'):
                            line.append('recessive')
                            line.append('pass')
                            line.append(known)
                            out.writerow(line)
                            outfiltered.writerow(line)
                # e.g. intronic variants fitting the criteria
                        else:
                            line.append(args.inheritance)
                            line.append('filtered')
                            line.append(known)
                            out.writerow(line)
                    else:
                        line.append(args.inheritance)
                        line.append('filtered')
                        line.append(known)
                        out.writerow(line)
                else:
                    line.append(args.inheritance)
                    line.append(known)
                    line.append('filtered_index')
                    out.writerow(line)
                continue
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.03:
                line.append('recessive')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
            
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
        
        ###
        # look for X linked recessive variants in trios
        ###
        elif args.inheritance == 'Xlinked':
            
            index_chromosome = identifycolumns(header, 'Chr')
            
            # skip all variants not located 
            if line[index_chromosome].lower() == 'x' or line[index_chromosome] == '23':
                pass
            else:
                continue
            
            judgement = xlinked(sampledata, family,names)
            
            # X linked should only work on trios
            if not args.familytype == 'trio':
                line.append('Trio_only')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
            
            # exclude gene, if it is on the exclusion list
            elif len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                line.append(known)
                out.writerow(line)
                continue
            
            # check before all others, if variant locates to simple tandem repeat region
            elif judgement == 1 and not tandem == 'NA':
                line.append('Xlinked')
                line.append('tandem')
                line.append(known)
                out.writerow(line)
                continue
            
            elif judgement == 1 and MAF <= 0.01:
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                    if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                        if (line[index_segdup] == '0'):
                            line.append('Xlinked')
                            line.append('pass')
                            line.append(known)
                            out.writerow(line)
                            outfiltered.writerow(line)
                # e.g. intronic variants fitting the criteria
                        else:
                            line.append(args.inheritance)
                            line.append('filtered')
                            line.append(known)
                            out.writerow(line)
                    else:
                        line.append(args.inheritance)
                        line.append('filtered')
                        line.append(known)
                        out.writerow(line)
                else:
                    line.append(args.inheritance)
                    line.append('filtered')
                    line.append(known)
                    out.writerow(line)
                continue
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.01:
                line.append('Xlinked')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
            
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                line.append(known)
                out.writerow(line)
                continue
        
        ###
        # look for compound heterozygous variants
        ###
        elif args.inheritance == 'compound':
            """
            valid combinations:
            child	parent	parent
            1/1 invalid
            0/1	0/0	0/0 (denovo)
            0/1	0/1	0/0
            0/1	0/0	0/1
            0/1	0/1	0/1
            
            invalid combinations:
            2&4 (parent 50% chance sick)
            3&4 (parent 50% chance sick)
            """
            
            new_gene = re.sub('\(.*?\)','',line[index_gene])
            new_gene_set = set(new_gene.split(';'))
            
            # if not old_gene == new_gene:
            # check if the names are the same
            # sometimes gene looks like 'PKHD1L1;PKHD1L1'
            # print line
            # print old_gene_set
            # print new_gene_set
            # print len(old_gene_set - new_gene_set)
            # print '2'
            if len(old_gene_set - new_gene_set) > 0:
                
                #pp.pprint(['old: ', old_gene, 'new: ',new_gene, 'orig: ', line[index_gene]])

                comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
                extension = []

                pass_ = 0
                if len(compound_gene_storage) == 1:
                    extension.extend(['NOT_compound','filtered'])
                else:
                    extension.append('compound')
                    if comp_judgement==1:
                        extension.append('pass')
                        pass_ = 1
                    else:
                        extension.append('filtered')                

                    for row in compound_gene_storage:
                        genecolumn2   = re.sub('\(.*?\)','',row[index_gene])
                        genenames2 = set(genecolumn2.split(';'))
                        
                        if len(genes_known & genenames2) >0:
                            known= 'yes'
                        else:
                            known = 'no'
                        row.extend(extension)
                        row.append(known)
                        out.writerow(row)
                        if pass_>0:outfiltered.writerow(row)

                # reset values
                compound_gene_storage = []
                old_gene     = new_gene
                old_gene_set = new_gene_set
                pass
            
            

            judgement = compound(sampledata, family,names)
            
            # top SNP
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                line.append(known)
                out.writerow(line)
            
            # check before all others, if variant locates to simple tandem repeat region
            elif judgement == 1 and not tandem == 'NA':
                line.append('compound')
                line.append('tandem')
                line.append(known)
                out.writerow(line)
           
            ## Exclude non damaging CADD variants            
            #elif judgement == 1 and CADD<=19 and CADD>0 :
            #    line.append('dominant')
            #    line.append('cadd')
            #    line.append(known)
            #    out.writerow(line)
            #    continue
           
           
            elif judgement == 1 and MAF <= 0.02:#0.03
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                    if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                        if (line[index_segdup] == '0'):
                            compound_gene_storage.append(line)
                            
                        else:
                            line.append(args.inheritance)
                            line.append('filteredsd')#
                            line.append(known)
                            out.writerow(line)
                    else:
                        line.append(args.inheritance)
                        line.append('filteredvf')#
                        line.append(known)
                        out.writerow(line)
                else:
                    line.append(args.inheritance)
                    line.append('filteredex')#
                    line.append(known)
                    out.writerow(line)
                
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.02:#0.03
                line.append('compound')
                line.append('filtered')#
                line.append(known)
                out.writerow(line)
            
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                line.append('filtered')
                line.append(known)
                out.writerow(line)

    else:
        # clean up for last gene
        if args.inheritance == 'compound':
            comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
            genecolumn   = re.sub('\(.*?\)','',line[index_gene])
            genenames = set(genecolumn.split(';'))

            if len(old_gene_set - new_gene_set) >= 0:
                #pp.pprint(['old: ', old_gene, 'new: ',new_gene, 'orig: ', line[index_gene]])
                comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
                extension = []
                pass_ = 0
                if len(compound_gene_storage) == 1: extension.extend(['NOT_compound','filtered'])
                else:
                    extension.append('compound')
                    if comp_judgement==1:
                        extension.append('pass')
                        pass_ = 1
                    else:
                        extension.append('filteredjg')                
                    for row in compound_gene_storage:
                        genecolumn2   = re.sub('\(.*?\)','',row[index_gene])
                        genenames2 = set(genecolumn2.split(';'))
                        
                        if len(genes_known & genenames2) >0:
                            known= 'yes'
                        else:
                            known = 'no'
                        row.extend(extension)
                        row.append(known)
                        out.writerow(row)
                        if pass_>0:outfiltered.writerow(row)
        
     ### write an xls output
    if writeXLS == True:
        ## open output file for re-reading
        args.filteredfile.close()
        fh = open(args.filteredfile.name, 'r+')
        # open output file for re-reading
        
        excel_name = 'variant_prioritization_report.xlsx'#args.filteredfile.name + ".xlsx"
        inheritance_file=args.filteredfile.name + ".xlsx"
        

        excel_path  =  os.path.dirname(args.outfile.name).split('/')
        excel_path  = '/'.join(excel_path[:-1])
    if excel_path=='':
        pass
    else:
        excel_path+='/'
    
    tmp_name = excel_path+'tmp.xlsx'
    excel_name = excel_path+excel_name
    
    # open xls file for writing
    print("Printing the Excel file:")
    xls = xlsxwriter.Workbook(tmp_name)
    sheet_name = 'ediva_filtered' + args.inheritance
    print(sheet_name)
    sheet_name=sheet_name[6:30]
    worksheet = xls.add_worksheet(sheet_name)
    row_xls=0
    fh.seek(0)
    ## DB parameters
    #username = 'edivapublic';
    #database = 'eDiVa_annotation';
    #dbhost = 'mysqlsrv-ediva.linux.crg.es';
    #passw = 'x86d2k1B';
        
    #db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    #user=username, # your username
    #passwd=passw)#, # your password
    #db=database) # name of the data base        
    #cur = db.cursor()
    
    
    # read line by line and transform to xls
    for line in fh:
        #line.rstrip('\n')
        data = line.strip().split(',')
        if row_xls>0:
            ## Adding Omim and Clinvar annotation 18-02-2015
            #gene_name = data[8]
            #sql = ("SELECT gene_name , title_mim_number ,details_mim_number "+
            #       "FROM eDiVa_annotation.Table_gene2mim, eDiVa_annotation.Table_omim "+
            #       "where eDiVa_annotation.Table_gene2mim.mim_number = eDiVa_annotation.Table_omim.mim_number "+
            #       "and eDiVa_annotation.Table_gene2mim.gene_name ='%s';"%gene_name)
                #sql = "select chr,pos,lengthofrepeat,copyNum,region from ediva_public_omics.Table_simpleRepeat;"
            #cur.execute(sql)

            #omim_disease = "."
            #omim_web="."
            
            #for row in cur:
            #    omim_disease+=str(row[1])+" | "
            #    omim_web+=str(row[2])+" | "
                
                
            #sql_clinvar =("SELECT clinical_significance, disease_name, clinical_review, access_number "+
            #              "FROM eDiVa_annotation.Table_clinvar "+
            #              "WHERE chr='%s' and pos='%s' and ref='%s'and alt='%s';"
            #              %(data[0],data[1],data[2],data[3])
            #)
            
            #cur.execute(sql_clinvar)         
            #clinvar_clinical_significance = "."
            #clinvar_disease_name = "."
            #clinvar_clinical_review = "."
            #clinvar_access_number = "."
            #for row in cur:
            #    clinvar_clinical_significance +=str(row[0])+"  "
            #    clinvar_disease_name +=str(row[1])+"  "
            #    clinvar_clinical_review +=str(row[2])+"  "
            #    clinvar_access_number +=str(row[3])+"  "
            #added_annotation= [omim_disease,omim_web,clinvar_disease_name,clinvar_access_number,clinvar_clinical_review,clinvar_clinical_significance]
            #added_annotation = [x[1:] if len(x)>1  else x for x in added_annotation]
            #data.extend(added_annotation)
            ######## -till here
            worksheet.write_row(row_xls, 0, data)
            #print row_xls
        else:
            #data.extend(['OMIM_name','OMIM_ID','clinical_significance', 'disease_name', 'clinical_review',' access_number'])
            worksheet.write_row(row_xls, 0, data)
        row_xls += 1
    #cur.close()
    #db.close()
    try:
        print("Hello")
        xls.close()
    except:
        print('error line 700')
        pass
        #xls = xlsxwriter.Workbook(tmp_name)
        print(inheritance_file)
        #shutil.copyfile(tmp_name,inheritance_file)
        os.system('mv %s %s'%(tmp_name,inheritance_file))
        xls = xlsxwriter.Workbook(tmp_name)
        #check if already exist
        
        try:
            print("Updating the existing file with the current inheritance sheet")
            workbook_rd = xlrd.open_workbook(inheritance_file)
            worksheets = workbook_rd.sheet_names()
            for worksheet_name in worksheets:
                try:
                    worksheet_rd = workbook_rd.sheet_by_name(worksheet_name)
                    worksheet_wr = xls.add_worksheet(worksheet_name)
                    
                    num_rows = worksheet_rd.nrows - 1
                    curr_row = -1
                    while curr_row < num_rows:
                            curr_row += 1
                            row_content = worksheet_rd.row(curr_row)
                            row_content_values = [x.value for x in row_content]
                            worksheet_wr.write_row(curr_row, 0, row_content_values)
                            #print row
             
                except:
                    print("There was a problem in processing %s sheet. \nIt may be because the sheet was already there before"%(worksheet_name))
                    #raise
            #os.remove(excel_name)
            #workbook_rd.close()
        
        except:
            
            pass
        
        try:
            print("Updating the existing file with the already existing sheets")
            workbook_rd = xlrd.open_workbook(excel_name)
            worksheets = workbook_rd.sheet_names()
            for worksheet_name in worksheets:
                try:
                    worksheet_rd = workbook_rd.sheet_by_name(worksheet_name)
                    worksheet_wr = xls.add_worksheet(worksheet_name)
                    
                    num_rows = worksheet_rd.nrows - 1
                    curr_row = -1
                    while curr_row < num_rows:
                            curr_row += 1
                            row_content = worksheet_rd.row(curr_row)
                            row_content_values = [x.value for x in row_content]
                            worksheet_wr.write_row(curr_row, 0, row_content_values)
                            #print row
                except:
                    print("There was a problem in processing %s sheet. \nIt may be because the sheet was already there before"%(worksheet_name))
                    #raise
            #os.remove(excel_name)
        
        except:
            pass

        
        fh.close()
        xls.close()
        #print excel_name
        #print inheritance_file
        #print tmp_name
        cmd='mv %s %s'%(tmp_name,excel_name)
        os.system(cmd)
        #os.rename(tmp_name, excel_name)   
    if args.csvfile != None and os.path.isfile(args.csvfile):
        #mailCmd = 'python '+ mailer_path +' -s /home/rrahman/soft/python-mailer/family.html '+ str(args.csvfile) +' Variant Prioritization'
        mailCmd = 'python ' + ' -s /home/rrahman/soft/python-mailer/family.html ' + str(args.csvfile) + ' Variant Prioritization'
        #print mailCmd
        os.system(mailCmd)
    exit(0)


###################################################
# sub routines
###################################################

def compound(sampledata, family,names,debug=False):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata#.split(';')

    check_samples = dict()
    sample_annot_size = int(len(sampledata)/len(names))

    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0

    for i in range(0,int(sample_annot_size*len(names)),sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[int(i/sample_annot_size)]
        # error catching because of wrong splitting,
        #e.g. 40ACVi>0/1>99>0.333;0.167,40ACVm>0/1>99>0.333;0.167,40ACVp>0/2>99>0.333;0.167
        if len(features) > 1 and name in family:
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
        
            #sub_pp.pprint([refcoverage, altcoverage])
            
            if name in check_samples:
                check_samples[name] = 1
    
            # homo alt is not expected in compound
            if zygosity == '1/1' and family[name] == '1':
                #sub_pp.pprint("dropped out in 1/1 1")
                judgement = 0
                if debug ==True:
                    print("763")
                    print(name)
                    print(family[name])
                    print(zygosity)
                break
            
            # het is good (could be an inherited variant or de novo)
            if zygosity == '0/1' and family[name] == '1':
                #sub_pp.pprint("dropped out in 0/1 1")
                judgement = 1
                continue
            
            # het or hom ref for parents is good
            elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '0':
                #sub_pp.pprint("dropped out in 0/0 0")
                judgement = 1
                continue
            
            # parents must not be hom alt
            elif zygosity == '1/1' and family[name] == '0':
                #sub_pp.pprint("dropped out in 1/1 0")
                judgement = 0
                if debug ==True:
                    print('781')
                break
            
            # offspring should have the variant
            elif zygosity == '0/0' and family[name] == '1':
                #sub_pp.pprint("dropped out in 0/0 1")
                judgement = 0
                if debug ==True:
                    print('789')
                break
            
            # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
            elif zygosity == './.' and family[name] == '0':
                # which chance has the current read distribution to miss out on an alt read
                # e.g. 10ref, 0alt
                # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
                # http://stattrek.com/online-calculator/poisson.aspx
                # use poisson distribution
                # poisson_miss = poisson.cdf(0.0, float(ND_coverage)/2)
                
                # if vcf file was not supplemented by pileup data
                # accept variants which could not be called in the parents
                if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                    judgement = 1
                    continue
                
                try:
                    int(refcoverage)
                except:
                    #sub_pp.pprint(sampledata)
                    exit(0)
                
                # hom ref
                #if int(refcoverage) >= 8 and int(altcoverage) == 0:
                #    #sub_pp.pprint("dropped out in ./. 0 8 0")
                #    judgement = 1
                #    continue
                refcoverage = float(refcoverage)
                altcoverage = float(altcoverage)
                
                coverage = refcoverage + altcoverage
                
                if coverage == 0:
                    judgement = 0
                    if debug==True :
                        print('824')
                    break
                
                # hom ref
                # poisson for low coverage and percentage for high coverage
                elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                    judgement = 1
                    continue
                
                # hom alt
                #elif int(altcoverage) >=8 and int(refcoverage) == 0:
                #    #sub_pp.pprint("dropped out in ./. 0 0 8")
                #    judgement = 0
                #    break
                
                # hom alt
                elif poisson.cdf( float(refcoverage), float(coverage/2) ) <= 0.007  and refcoverage / coverage <= 0.05:
                    judgement = 0
                    if debug ==True:
                        print('824')
                    break
                
                # coverage too low?
                else:
                    #sub_pp.pprint("dropped out in ./. 0 else")
                    judgement = 0
                    if debug ==True:
                        print("850")
                    break
            
            # do not accept missing values for affected individuals
            elif zygosity == './.' and family[name] == '1':
                #sub_pp.pprint("dropped out in ./. 1")
                judgement = 0
                if debug ==True:
                    print('857')
                break
    

    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            if debug ==True:
                print('863')
                print(check_samples)
            break
    if debug:
        print(check_samples.values())
        print(judgement)
    return(judgement)

def compoundizer(variantlist, family, index_sample,names):
    
    
    sub_pp = pprint.PrettyPrinter(indent = 5)
    
    #sub_pp.pprint(variantlist)
    
    # initialize
    ticker_dict = dict()
    for member in family.keys():
        if family[member] == '0':
            ticker_dict[member] = list()
    #Parents names    
    name1 = list(ticker_dict.keys())[0]
    name2 = list(ticker_dict.keys())[1]

    judgement = 0
    
    # check line by line, if this variant could support a compound het

    for variantline in variantlist:
        
        # produce a list with all the sample data
        sampledata = variantline[index_sample:len(variantline)-1]#.split(';')
        
        sample_annot_size = int(len(sampledata)/len(names))
        # produce a dictionary to save the zygosities
        zygosities = dict()
        
        #if not len(sampledata) == 3:
        #    print 'Something went wrong. The samples provided did not contain a trio case?'
        #    judgement = 0
        #    break
        #

        for i in range(0,int(sample_annot_size*len(names)),sample_annot_size):
            sam = sampledata[i]
            features    = sampledata[i:i+sample_annot_size]#sam.split(':')
            #print sampledata
            #print names
            name        = names[int(i/sample_annot_size)]
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
            # check if, we are looking at the offspring
            
            if not list(ticker_dict.keys())[0] == name and not list(ticker_dict.keys())[1] == name:
                continue
            
            zygosities[name] = zygosity

        # check the entered values for possible compound supporters
        # denovo
        if zygosities[name1]   == '0/0' and zygosities[name2] == '0/0':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(0)
            pass
    
        elif zygosities[name1]   == '0/1' and zygosities[name2] == '0/0':
            ticker_dict[name1].append(1)
            ticker_dict[name2].append(0)
            pass
        
        elif zygosities[name1]   == '0/0' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(1)
            pass
        
        elif zygosities[name1]   == '0/1' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(0)
            pass
        
        elif zygosities[name1]   == '0/1' and zygosities[name2] == './.':
            ticker_dict[name1].append(1)
            ticker_dict[name2].append(0)
            pass
        
        elif zygosities[name1]   == './.' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(1)
            pass
        else:
            pass
    
    #sub_pp.pprint(ticker_dict)
    
    # check if there are enough supporters
    name1_sum = sum(ticker_dict[name1])
    name2_sum = sum(ticker_dict[name2])

    
    if (name1_sum + name2_sum) >= 2 and name1_sum >= 1 and name2_sum >= 1:
        judgement = 1
    else:
        judgement = 0
        

        
    return(judgement)

def denovo(sampledata, family,names):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    samples = sampledata#.split(';')
    judgement = 0
    
    check_samples = dict()
    
    # create data structure for completeness check
    for sam in family.keys():
         check_samples[sam] = 0

    # go into the variant data
    sample_annot_size = len(sampledata)/len(names)
    for i in range(0,int(sample_annot_size*len(names)),int(sample_annot_size)):
        sam = samples[i]
        features    = samples[i:int(i+sample_annot_size)]#sam.split(':')
        name        = names[int(i/sample_annot_size)]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            
            continue
        
        if len(features)>=3:
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity    = features[0]
            refcoverage = '.'
            altcoverage = '.'
        # check if sample is found in pedigree
        
        # sample info complete?
        if name in check_samples:
            check_samples[name] = 1
        
        #sub_pp.pprint(family)
        
        # heterozygous in affected individual - good
        if zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue
        
        # hom ref, not affected - good
        elif zygosity == '0/0' and family[name] == '0' :
            if int(altcoverage)<3 : 
                judgement = 1
            else :
                judgement =0
            continue
        
        # heterozygous in non-affected - bad
        elif zygosity == '0/1' and family[name] == '0':
            #sub_pp.pprint("heterozygous in non-affected")# DEBUG
            judgement = 0
            break
        
        # hom ref in affected - bad
        elif zygosity == '0/0' and family[name] == '1':
            #sub_pp.pprint("hom ref in affected")# DEBUG
            judgement = 0
            break
        
        # homozygous can't be denovo
        elif zygosity == '1/1':
            #sub_pp.pprint("homozygous can't be denovo")# DEBUG
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            
            # if vcf file was not supplemented by pileup data
            # reject it variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 0
                continue
            
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                #sub_pp.pprint("coverage")# DEBUG
                judgement = 0
                break
            
            # hom ref, non called genotype
            # poisson for low coverage and percentage for high coverage
            # poisson 10 reads (poisson average rate of success = 5) and alt reads = 0 - should get still accepted
            elif (poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007) and (altcoverage / coverage <= 0.05):
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
                judgement = 1
                continue
            
            ## het, non called genotype
            #elif int(altcoverage) >=1:
            #    judgement = 0
            #    break
            
            # not necessary to check for hom alt
            
            # coverage too low?
            else:
                #sub_pp.pprint("else")# DEBUG
                #sub_pp.pprint(poisson.cdf( float(altcoverage), float(coverage)/2 ))# DEBUG
                #sub_pp.pprint(altcoverage / coverage)# DEBUG
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':
            
            # except if vcf file was not supplemented by pileup data
            # reject variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.':
                judgement = 0
                continue
            
            # do not be that grateful, if there is coverage data
            # if the SNP caller could not call a genotype in an area, where there was coverage
            # we rather trust the SNP caller, than starting to call SNP based on pileup coverage
            else:
                judgement = 0
                break
    
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break
    
    return(judgement)

def dominant(sampledata, family,names):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata#.split(';')
    check_samples = dict()
    
    for samp in family.keys():
         check_samples[samp] = 0
    
    judgement = 0
    
    sample_annot_size = len(sampledata)/len(names)
    for i in range(0,sample_annot_size*len(names),sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[i/sample_annot_size]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            
            continue
        
        if len(features)>=3:
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity    = features[0]
            refcoverage = '.'
            altcoverage = '.'

        if check_samples.has_key(name):
            check_samples[name] = 1
    
        # affected family members should have the mutation (hom ref not allowed)
        if zygosity == '0/0' and family[name] == '1':
            judgement = 0
            break
        
        # affected family members might be het
        elif zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue
        
        # affected family members might be hom alt
        # that's the major difference to de novo...
        elif zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue
        
        # non-affected family members must not have the mutation - hom ref is OK
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue
        
        # non-affected family members must not have the mutation - het is bad
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 0
            break
        
        # non-affected family members must not have the mutation - hom alt is worst
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            
            # if vcf file was not supplemented
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue
            
            # do not do any other judgements, i.e.
            # being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            # also tolerate low coverage
            judgement = 1
            continue
        
        # accept some missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue
            
            ### do not do any other judgements, i.e.
            ### being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            ### also tolerate low coverage
            ##judgement = 1
            ##continue
            
            # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
            
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                judgement = 0
                break
            
            # hom ref
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                #sub_pp.pprint(['./. 0 8 0', name])
                break
            
            # hom alt
            #elif int(altcoverage) >=8 and int(refcoverage) == 0:
            elif poisson.cdf( float(refcoverage), float(coverage/2) ) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 1
                #sub_pp.pprint(['./. 0 0 8', name])
                continue
            
            # het
            #elif int(refcoverage) >= 8 and not int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                #sub_pp.pprint(['./. 0 8 not 0', name])
                continue
        
        else:
            pass
    
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break
    #if judgement==1:
    #    print samples
    #    print family
        
    return(judgement)

def identifycolumns (header, question):
    try:
        index = header.index(question)
    except:
        exit("ERROR: %s column could not be identified in the annotated data" % question)
    return( int(index) )

def recessive(sampledata, family, familytype,names):
    #sub_pp = pprint.PrettyPrinter(indent = 8)
    
    # get samples as data structure
    all_samples = dict()
    samples = sampledata#.split(';')

    check_samples = dict()
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    
    sample_annot_size = len(sampledata)/len(names)
    for i in range(0,sample_annot_size*len(names),sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[i/sample_annot_size]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            
            continue
        
        if len(features)>=3:
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity    = features[0]
            refcoverage = '.'
            altcoverage = '.'
            
        
        if check_samples.has_key(name):
            check_samples[name] = 1
        
        # affected individuals have to be homozygous
        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            #sub_pp.pprint(['1/1 1', name])
            continue
        
        # affected individuals should not be hom ref or het
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            #sub_pp.pprint(['0/0 0/1 1 ', name])
            break
        
        # non-affected individuals might be het
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 1
            #sub_pp.pprint(['0/1 0', name])
            continue
        
        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'family':
            judgement = 1
            #sub_pp.pprint(['0/0 0 family', name])
            continue
        
        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'trio':
            judgement = 0
            #sub_pp.pprint(['0/0 0 trio', name])
            break
        
        # non-affected individuals must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            #sub_pp.pprint(['1/1 0', name])
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                #sub_pp.pprint(['./. 0 . .', name])
                continue
            
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                judgement = 0
                break
            
            # hom ref
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                #sub_pp.pprint(['./. 0 8 0', name])
                break
            
            # hom alt
            #elif int(altcoverage) >=8 and int(refcoverage) == 0:
            elif poisson.cdf( float(refcoverage), float(coverage/2) ) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 0
                #sub_pp.pprint(['./. 0 0 8', name])
                break
            
            # het, which is OK
            #elif int(refcoverage) >= 8 and not int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                #sub_pp.pprint(['./. 0 8 not 0', name])
                continue
            
            # coverage too low?
            else:
                
                # accept missing values in family interrogations
                if familytype == 'family':
                    judgement = 1
                    #sub_pp.pprint(['family', name])
                    continue
                # do not accept missing values in trio setups
                elif familytype == 'trio':
                    judgement = 0
                    #sub_pp.pprint(['trio 0', name])
                    break
                
                # for security reasons
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == './.' and family[name] == '1':
            
            judgement = 0
            break
    
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return(judgement)

def xlinked(sampledata, family,names):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    
    # get samples as data structure
    all_samples = dict()
    samples = sampledata#.split(';')

    check_samples = dict()
    
    # collect the two parents - one should be het (mother), one whould be hom ref (father)
    # finally should contain two samples, one het & one hom ref
    inheritance_logic = dict()
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    
    sample_annot_size = len(sampledata)/len(names)
    for i in range(0,sample_annot_size*len(names),sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[i/sample_annot_size]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            
            continue
        
        if len(features)>=3:
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
            
        else:
            #stick with genotype and the others are empty
            zygosity    = features[0]
            refcoverage = '.'
            altcoverage = '.'
        
        if check_samples.has_key(name):
            check_samples[name] = 1
        
        if family[name] == '0':
            inheritance_logic[name] = zygosity
        
        # affected individuals have to be homozygous
        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue
        
        # affected individuals should not be hom ref or het
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            break
        
        # non-affected individuals might be het
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 1
            continue
        
        # non-affected individuals might be hom ref
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue
        
        # non-affected individuals must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue
            
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                judgement = 0
                break
            
            # hom ref
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                inheritance_logic[name] = '0/0'
                judgement = 1
                continue
            
            # hom alt
            #elif int(altcoverage) >=8 and int(refcoverage) == 0:
            if poisson.cdf( float(refcoverage), float(coverage)/2 ) <= 0.007 and refcoverage / coverage <= 0.05:
                inheritance_logic[name] = '1/1'
                judgement = 0
                break
            
            # het, which is OK
            #elif int(refcoverage) >= 8 and not int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) >= 0.007 or altcoverage / coverage >= 0.05:
                inheritance_logic[name] = '0/1'
                judgement = 1
                continue
            
            # coverage too low?
            else:
                
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == './.' and family[name] == '1':
            
            judgement = 0
            break
    
    # sanity check
    het_checker = 0
    hom_checker = 0
    
    for values in inheritance_logic.values():
        if values == '0/1':
            het_checker = 1
        if values == '0/0':
            hom_checker = 1
    if het_checker == 1 and hom_checker == 1:
        judgement = 1
    else:
        judgement = 0
    
    # another sanity check
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return(judgement)




main(args)
