#!/usr/bin/env python


###################################################
# sub routines
###################################################

def compound(sampledata, family,names,debug=False):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata#.split(';')
    check_samples = dict()
    sample_annot_size = len(sampledata)/len(names)
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    #define type of  compount variant
    # if less than two samples have 0/1 we are in front of a de novo
    if ''.join(samples).count('0/1') <=1: #de novo
        type_check = 'compound_denovo'
    else:
        type_check ='compound'
    for i in range(0,sample_annot_size*len(names),sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[i/sample_annot_size]
        # error catching because of wrong splitting,
        #e.g. 40ACVi>0/1>99>0.333;0.167,40ACVm>0/1>99>0.333;0.167,40ACVp>0/2>99>0.333;0.167
        if len(features) > 1 and family.has_key(name):
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
          
            if check_samples.has_key(name):
              check_samples[name] = 1
  
          # homo alt is not expected in compound
            if zygosity== '0/0':
                if family[name]== '0':  judgement = 1                
                else:                   judgement = 0
                    
            elif zygosity == '1/1':
                if family[name] == '0': judgement = 0
                else:                   judgement = 1 #parents can be hom alt
                    
            elif zygosity=='0/1':       judgement = 1
            else: #'./.' or weirdos
                if family[name] =='0':  judgement = 0
                else:                   judgement = 0
                    
            judgement *= evaluate_genotype_quality(zygosity,refcoverage,altcoverage,type_check)
            ## escape element
        if judgement==0:
            break

    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            if debug ==True:
                print '863'
                print check_samples
            break
    if debug:
        print check_samples.values()
        print judgement
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
    name1 = ticker_dict.keys()[0]
    name2 = ticker_dict.keys()[1]
    
    judgement = 0
    
    # check line by line, if this variant could support a compound het
    
    for variantline in variantlist:
        # produce a list with all the sample data
        sampledata = variantline[index_sample:len(variantline)-1]#.split(';')
        sample_annot_size = len(sampledata)/len(names)
        # produce a dictionary to save the zygosities
        zygosities = dict()
        
        #if not len(sampledata) == 3:
        #    print 'Something went wrong. The samples provided did not contain a trio case?'
        #    judgement = 0
        #    break
        #
        for i in range(0,sample_annot_size*len(names),sample_annot_size):
            sam = sampledata[i]
            features    = sampledata[i:i+sample_annot_size]#sam.split(':')
            #print sampledata
            #print names
            name        = names[i/sample_annot_size]
            zygosity    = features[0]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
            # check if, we are looking at the offspring
            
            if not ticker_dict.keys()[0] == name and not ticker_dict.keys()[1] == name:
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
    for i in range(0,sample_annot_size*len(names),sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[i/sample_annot_size]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print family
            
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
        if check_samples.has_key(name):
                  check_samples[name] = 1
        
        #sub_pp.pprint(family)
                ## check for the genotypes and 
        if zygosity== '0/0':
            if family[name]== '0':  judgement=1                
            else:                   judgement = 0
                
        elif zygosity == '1/1':     judgement = 0
                
        elif zygosity=='0/1':
            if family[name] == '0': judgement = 0
            else:                   judgement = 1
        else: #'./.' or weirdos
            if family[name] =='0':  judgement = 0
            else:                   judgement = 0
                
        judgement *= evaluate_genotype_quality(zygosity,refcoverage,altcoverage,'dominant_denovo')
        ## escape element
        if judgement==0:
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
    total_affected =0
    #print 'start'
    for i in range(0,sample_annot_size*len(names),sample_annot_size):
        
        sam = samples[i]
        
        features    = samples[i:i+sample_annot_size]#sam.split(':')
        name        = names[i/sample_annot_size]
        total_affected += int(family[name])
        
        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print family            
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
   
        if zygosity== '0/0':
            if family[name]== '0':  judgement = 1                 
            else:                   judgement = 0
                
        elif zygosity == '1/1':
            if family[name] == '0': judgement = 0
            else:                    judgement= 1
                
        elif zygosity=='0/1':
            if family[name] == '0': judgement = 0
            else:                   judgement = 1
        else: #'./.' or weirdos
            if family[name] =='0':  judgement = 0
            else:                   judgement = 0
        #print zygosity
        #print family[name]
        #print name
        
        judgement *= evaluate_genotype_quality(zygosity,refcoverage,altcoverage,'dominant_inherited')
        ## escape element
        if judgement==0:
            break
    
    #print len(samples)
    #print total_affected
    #print sample_annot_size
    #print family
    #print samples
    if total_affected ==1: judgement = 0 #because it's a denovo
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
            print family
            
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
            
        ## check for the genotypes and 
        if zygosity== '0/0':
            if family[name]== '0':
                if familytype=='trio':  judgement = 0
                else :                  judgement = 1
            else:                       judgement = 0
                
        elif zygosity == '1/1':
            if family[name] == '0': judgement = 0
            else:                   judgement = 1
                
        elif zygosity=='0/1':
            if family[name] == '0': judgement = 1
            else:                   judgement = 0
        else: #'./.' or weirdos
            if family[name] =='0':  judgement = 0
            else:                   judgement =0
        
        
        judgement *= evaluate_genotype_quality(zygosity,refcoverage,altcoverage,'recessive')
        ## escape element
        if judgement==0:
            break
    
    for vals in check_samples.values():
        # it's a double check in case we didn't process all samples
        #it means judgement is 0'  so we exit
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
            print family
            
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
         ## check for the genotypes and 
        if zygosity== '0/0':
            if family[name]== '0':judgement = 1
            else:                 judgement = 0 #0/0 1
                
        elif zygosity == '1/1':
            if family[name] == '0':judgement = 0
            else:                  judgement = 1 #1/1 1
                
        elif zygosity=='0/1':
            if family[name] == '0':judgement = 1
            else:                  judgement = 0
        else: #'./.' or weirdos
            if family[name] =='0': judgement = 0
            else:                  judgement = 0
                
        judgement *= evaluate_genotype_quality(zygosity,refcoverage,altcoverage,'Xlinked')
        ## escape element
        if judgement==0:
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
        judgement *= 1
    else:
        judgement *= 0
    
    # another sanity check
    for vals in check_samples.values():
        if vals == 0:
            judgement *= 0
            break

    return(judgement)

def check_thresholds(args,line,genes2exclude,genenames,indexes,MAF_threshold,judgement,out,outfiltered,known,compound_gene_storage=None):
    (index_MAF1k,index_MAFevs,index_MAF_exac,index_function,index_varfunction,index_segdup,index_gene,index_str,index_qual)=indexes
    
    ##MAF as maximum 1000G EVS and EXAC    
    try:
        MAF1k      = line[index_MAF1k].replace('NA','0')
        MAFevs     = line[index_MAFevs].replace('NA','0')
        if 'NA' == line[index_MAF_exac]:
            MAFexac ='0'
        else:
            MAFexac= line[index_MAF_exac]
    except:
        pp.pprint(line)
    try:
        # avoiding problems with NAs
        MAF    = max(float(MAF1k), float(MAFevs),float(MAFexac))

    except:
        print ('Freq error 1k %s EVS %s ExAC %s')%(MAF1k,MAFevs,MAFexac)
        MAF    = 0
    tandem = line[index_str]
    if line[index_qual]=='.':
       line[index_qual] = '10'
    #Exclusion list check    
    if len(genes2exclude & genenames) > 0:
        line.append('NOT_' + args.inheritance)
        line.append('exclusionlist')
        line.append(known)
        out.writerow(line)
        return 
    if judgement ==1 and float(line[index_qual] ) < 10:
        line.append(args.inheritance)
        line.append('filtered QUAL')
        line.append(known)
        out.writerow(line)
    # check before all others, if variant locates to simple tandem repeat region
    elif judgement == 1 and not tandem == 'NA':
        line.append(args.inheritance)
        line.append('filtered tandem')
        line.append(known)
        out.writerow(line)
        return 
    
    elif judgement == 1 and MAF <= MAF_threshold:
        
        if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
            if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                if (line[index_segdup] == '0'):
                    if args.inheritance != 'compound':
                        line.append(args.inheritance)
                        line.append('pass')
                        line.append(known)
                        outfiltered.writerow(line)
                        out.writerow(line)
                    else:
                        compound_gene_storage.append(line)
                        out.writerow(line)
                    
                    
                    return
        
                else:
                    # e.g. variant in seg dup area
                    line.append(args.inheritance)
                    line.append('filtered seg_dup')
                    line.append(known)
                    out.writerow(line)
                    return
            else:
                # e.g. var function synonimous 
                line.append(args.inheritance)
                line.append('filtered varfunction')
                line.append(known)
                out.writerow(line)
                return
        else:
            # e.g. intronic variants fitting the criteria
            line.append(args.inheritance)
            line.append('filtered non coding')
            line.append(known)
            out.writerow(line)
            return
        return 
    
    # fits inheritance, but is too frequent in the population
    elif judgement == 1 and MAF > MAF_threshold:
        line.append(args.inheritance)
        line.append('filtered MAF too high')
        line.append(known)
        out.writerow(line)
        return 

    # does not fit anything
    else:
        line.append('NOT_' + args.inheritance)
        #line.append('bad_inheritance')
        line.append('filtered')
        line.append(known)
        out.writerow(line)
        return        
    
    return None

def evaluate_genotype_quality(zygosity,refcoverage,altcoverage,inheritance):

    if refcoverage == '.' or refcoverage == '' or altcoverage==',' or altcoverage== '':
        return 1
    # if we don't have information about coverage we trust the genotyping.
    
    try:
        ref =int(refcoverage)

    except:
        print 'error' 
        print refcoverage
        print zygosity
        print altcoverage
        print inheritance
        raise
    alt =int(altcoverage)
    
    ### inheritnace based:
    
    if inheritance == 'recessive':
        if zygosity == '0/0': return 1
        elif zygosity == '0/1':
            if ref+alt >=6 and ref >= 2 and alt >= 2: return 1
            else :  return 0
        elif zygosity == '1/1':
            if alt >=5 and float(alt)/(ref+alt)> 0.9: return 1
            else: return 0
        else :
            return 1

    elif inheritance == 'dominant_denovo':
        if zygosity == '0/0':
            if ref+alt >=8 and alt <=1 : return 1 # 8,1
            else: return 0
        elif zygosity == '0/1':
            if alt >=3 and float(alt)/(ref+alt)>= 0.2: return 1 # 5, 0.2
            else :  return 0
        elif zygosity == '1/1':
            if alt >=3 and float(alt)/(ref+alt)> 0.6: return 1 # 5, 0.9
            else: return 0
        else :
            return 1
    
    elif inheritance == 'dominant_inherited':
        if zygosity == '0/0':
            if ref+alt >=3 and alt <=2 : return 1 #8 ,1
            else: return 0
        elif zygosity == '0/1':
            if alt >=5 and float(alt)/(ref+alt)>= 0.2: return 1 # 5 0.2
            else :  return 0
        elif zygosity == '1/1':
            if alt >=5 and float(alt)/(ref+alt)> 0.3: return 1 # 5 0.9
            else: return 0
        else :
            return 1
    
    elif inheritance == 'compound':
        if zygosity == '0/0':
            if ref+alt >=8 and alt <=1 : return 1  # 8, 1
            else: return 0
        elif zygosity == '0/1':
            if ref+alt >=6 and ref >= 2 and alt >= 2: return 1 # 6 2,2
            else :  return 0
        elif zygosity == '1/1':
            if alt >=5 and float(alt)/(ref+alt)> 0.9: return 1 # 5,0.9
            else: return 0
        else :
            return 1
    
    elif inheritance == 'compound_denovo':
        #0/1 need for denovo stringency
        if zygosity == '0/0':
            if ref+alt >=8 and alt <=1 : return 1  # 8, 1
            else: return 0
        elif zygosity == '0/1':
            if alt >=2 and float(alt)/(ref+alt)>= 0.20: return 1 #2 and 0.20
            else :  return 0
        elif zygosity == '1/1':
            if alt >=5 and float(alt)/(ref+alt)> 0.9: return 1
            else: return 0
        else : return 1
        
    elif inheritance =='Xlinked':
        if zygosity == '0/0':
            if ref+alt >=8 and alt <=1 : return 1
            else: return 0
        elif zygosity == '0/1':
            if alt >=5 and float(alt)/(ref+alt)>= 0.2: return 1
            else :  return 0
        elif zygosity == '1/1':
            if alt >=5 and float(alt)/(ref+alt)> 0.9: return 1
            else: return 0
        else :
            return 1
        
        
        
        
if __name__=='__main__':
    import argparse
    import csv
    import pprint
    import re
    import cPickle as pickle
    from scipy.stats import poisson
    import os
    import MySQLdb
    import shutil
    import sys
    from operator import itemgetter

    
    try:
        import xlsxwriter
        import xlrd
        writeXLS = True
    except:
        writeXLS = False
        print 'No XlS writer'
    
    sample_annot_size = 6 # sampleid - DP - REF - ALT - AF - GQ
    
    
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
    mailer_path='/home/rrahman/soft/python-mailer/pymailer.py'
    pp = pprint.PrettyPrinter( indent=4) # DEBUG
    names = list()
    # read the gene exclusion list
    # an empty set will not filter out anything, if gene exclusion list is not provided
    genes2exclude = set()
    if args.geneexclusion:
        for gene in args.geneexclusion:
            gene = gene.rstrip()
            genes2exclude.add(gene)
    if args.white_list == None:
        args.white_list = 'None'
        related_genes=set()
    else:
        related_genes=set()
        if os.path.isfile(args.white_list):
            script_dir = os.path.dirname(sys.argv[0])
            with open(args.white_list,'r') as w:
                HPO_dict  = pickle.load(open(script_dir +'/HPO_gene_assiciation.p','rb'))
                related_genes = list()
                for line in w:
                    HPO_term = line.rstrip('\n')
                    try:
                        related_genes.extend(HPO_dict[HPO_term])
                    except:
                        print '%s not found in database'%(HPO_term)
        else:
            print 'The specified HPO list %s is not a valid file'%(args.white_list)
            print 'eDiVA will proceed as without any HPO list'
                
        genes2exclude = set(genes2exclude) - set(related_genes)
        related_genes= set(related_genes)
        
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
    index_qual         = identifycolumns(header, 'QUAL')
    
    indexes=[index_MAF1k,index_MAFevs,index_MAF_exac,index_function,index_varfunction,index_segdup,index_gene,index_str,index_qual]
    
    #print index_sample
    print family.keys()
    for i in range(index_sample,index_sample+sample_annot_size*len(family.keys()),sample_annot_size):
        names.append(header[i])

    print header[index_sample]
    
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
       
        # read sample names and according zygosity NOW IT's a LIST
        sampledata = line[index_sample:index_sample+sample_annot_size*len(family.keys())]

        if len(sampledata)==0:
            print line
            print len(line)
            raise
        elif sampledata[0] =='0':
            print line
            raise
        
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
            MAF_threshold=0.01
            judgement = denovo(sampledata, family,names)
            check_thresholds(args,line,genes2exclude,genenames,indexes,MAF_threshold,judgement,out,outfiltered,known)      
        ###
        # look for familial dominant variants. (being tolerant for missing values)
        ###
        elif args.inheritance == 'dominant_inherited':
            MAF_threshold=0.01#05
            judgement = dominant(sampledata, family,names)
            check_thresholds(args,line,genes2exclude,genenames,indexes,MAF_threshold,judgement,out,outfiltered,known)
        ###
        # look for recessive variants (be aware of trio and family inheritance)
        ###
        elif args.inheritance == 'recessive':
            MAF_threshold=0.03
            judgement = recessive(sampledata, family, args.familytype,names)
            check_thresholds(args,line,genes2exclude,genenames,indexes,MAF_threshold,judgement,out,outfiltered,known)
            
        ###
        # look for X linked recessive variants in trios
        ###
        elif args.inheritance == 'Xlinked':            
            index_chromosome = identifycolumns(header, 'Chr')
            # skip all variants not located 
            if line[index_chromosome].lower() == 'x' or line[index_chromosome] == '23':
                MAF_threshold=0.01
                judgement = xlinked(sampledata, family,names)
                
            # X linked should only work on trios
                if not args.familytype == 'trio':
                    line.append('Trio_only')
                    line.append('filtered')
                    line.append(known)
                    out.writerow(line)
                else:
                    check_thresholds(args,line,genes2exclude,genenames,indexes,MAF_threshold,judgement,out,outfiltered,known)

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
            if len(old_gene_set - new_gene_set) > 0:
                
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
                        extension.append('filtered')                

                    for row in compound_gene_storage:
                        genecolumn2   = re.sub('\(.*?\)','',row[index_gene])
                        genenames2 = set(genecolumn2.split(';'))
                        
                        if len(related_genes & genenames2) >0:
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
            
            MAF_threshold=0.02
            judgement = compound(sampledata, family,names)
            check_thresholds(args,line,genes2exclude,genenames,indexes,MAF_threshold,judgement,out,outfiltered,known,compound_gene_storage)
               
    else:
        # clean up for last gene
        if args.inheritance == 'compound':
            comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
            genecolumn   = re.sub('\(.*?\)','',line[index_gene])
            genenames = set(genecolumn.split(';'))
        
            if len(old_gene_set - new_gene_set) > 0:
                
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
                        extension.append('filtered')                
                    for row in compound_gene_storage:
                        genecolumn2   = re.sub('\(.*?\)','',row[index_gene])
                        genenames2 = set(genecolumn2.split(';'))
                        
                        if len(related_genes & genenames2) >0:
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
        print "Printing the Excel file:"
        xls = xlsxwriter.Workbook(tmp_name)
        sheet_name = 'ediva_filtered' + args.inheritance
	print sheet_name
        sheet_name=sheet_name[6:30]
        worksheet = xls.add_worksheet(sheet_name)
        row_xls=0
        fh.seek(0)
        ## DB parameters
        username = 'edivapublic';
        database = 'eDiVa_annotation';
        dbhost = 'mysqlsrv-ediva.linux.crg.es';
        passw = 'x86d2k1B';
         
        db = MySQLdb.connect(host=dbhost, # your host, usually localhost
        user=username, # your username
        passwd=passw)#, # your password
        #db=database) # name of the data base        
        cur = db.cursor()
        
        
        # read line by line and transform to xls
        for line in fh:
            #line.rstrip('\n')
            data = line.strip().split(',')
            if row_xls>0:
                ## Adding Omim and Clinvar annotation 18-02-2015
                gene_name = data[8]
                sql = ("SELECT gene_name , title_mim_number ,details_mim_number "+
                       "FROM eDiVa_annotation.Table_gene2mim, eDiVa_annotation.Table_omim "+
                       "where eDiVa_annotation.Table_gene2mim.mim_number = eDiVa_annotation.Table_omim.mim_number "+
                       "and eDiVa_annotation.Table_gene2mim.gene_name ='%s';"%gene_name)
                    #sql = "select chr,pos,lengthofrepeat,copyNum,region from ediva_public_omics.Table_simpleRepeat;"
                cur.execute(sql)
    
                omim_disease = "."
                omim_web="."
                
                for row in cur:
                    omim_disease+=str(row[1])+"  "
                    omim_web+=str(row[2])+"  "
                    
                    
                sql_clinvar =("SELECT clinical_significance, disease_name, clinical_review, access_number "+
                              "FROM eDiVa_annotation.Table_clinvar "+
                              "WHERE chr='%s' and pos='%s' and ref='%s'and alt='%s';"
                              %(data[0],data[1],data[2],data[3])
                )
                
                cur.execute(sql_clinvar)         
                clinvar_clinical_significance = "."
                clinvar_disease_name = "."
                clinvar_clinical_review = "."
                clinvar_access_number = "."
                for row in cur:
                    clinvar_clinical_significance +=str(row[0])+"  "
                    clinvar_disease_name +=str(row[1])+"  "
                    clinvar_clinical_review +=str(row[2])+"  "
                    clinvar_access_number +=str(row[3])+"  "
                added_annotation= [omim_disease,omim_web,clinvar_disease_name,clinvar_access_number,clinvar_clinical_review,clinvar_clinical_significance]
                added_annotation = [x[1:] if len(x)>1  else x for x in added_annotation]
                data.extend(added_annotation)
                ######## -till here
                worksheet.write_row(row_xls, 0, data)
                #print row_xls
            else:
                data.extend(['OMIM_name','OMIM_ID','clinical_significance', 'disease_name', 'clinical_review',' access_number'])
                worksheet.write_row(row_xls, 0, data)
            row_xls += 1
        cur.close()
        db.close()
	try:
         xls.close()
	except:
                print 'error line 700'
		pass
        #xls = xlsxwriter.Workbook(tmp_name)
        print inheritance_file
        #shutil.copyfile(tmp_name,inheritance_file)
        os.system('mv %s %s'%(tmp_name,inheritance_file))
        xls = xlsxwriter.Workbook(tmp_name)
        #check if already exist
        
        try:
            print "Updating the existing file with the current inheritance sheet"
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
                    print "There was a problem in processing %s sheet. \nIt may be because the sheet was already there before"%(worksheet_name)
                    #raise
            #os.remove(excel_name)
            #workbook_rd.close()
        
        except:
            
            pass
            
        
        try:
            print "Updating the existing file with the already existing sheets"
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
                    print "There was a problem in processing %s sheet. \nIt may be because the sheet was already there before"%(worksheet_name)
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
    
        if args.csvfile!=None and os.path.isfile(args.csvfile):
            mailCmd = 'python '+ mailer_path +' -s /home/rrahman/soft/python-mailer/family.html '+ str(args.csvfile) +' Variant Prioritization'
            #print mailCmd
            os.system(mailCmd)
