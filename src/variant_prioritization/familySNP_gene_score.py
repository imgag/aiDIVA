import argparse
import csv
import re
from scipy.stats import poisson
import os
import shutil
import pickle
import networkx
import sys
import pandas as pd
from operator import itemgetter
import get_HPO_similarity_score as gs


variant_consequences = {'transcript_ablation': 'non_exonic',
                        'splice_acceptor_variant': 'exonic;splicing',
                        'splice_donor_variant': 'exonic;splicing',
                        'stop_gained': 'exonic',
                        'frameshift_variant': 'exonic',
                        'stop_lost': 'exonic',
                        'start_lost': 'exonic',
                        'transcript_amplification': 'non_exonic',
                        'inframe_insertion': 'exonic',
                        'inframe_deletion': 'exonic',
                        'missense_variant': 'exonic',
                        'protein_altering_variant': 'exonic',
                        'splice_region_variant': 'splicing', # maybe not
                        'incomplete_terminal_codon_variant': 'exonic',
                        'start_retained_variant': 'exonic', # maybe not
                        'stop_retained_variant': 'exonic', # maybe not
                        'synonymous_variant': 'exonic',
                        'coding_sequence_variant': 'exonic',
                        'mature_miRNA_variant': 'non_exonic',
                        '5_prime_UTR_variant': 'UTR_exonic',
                        '3_prime_UTR_variant': 'UTR_exonic',
                        'non_coding_transcript_exon_variant': 'non_exonic',
                        'intron_variant': 'intronic',
                        'NMD_transcript_variant': 'non_exonic',
                        'non_coding_transcript_variant': 'non_exonic',
                        'upstream_gene_variant': 'non_exonic',
                        'downstream_gene_variant': 'non_exonic',
                        'TFBS_ablation': 'non_exonic',
                        'TFBS_amplification': 'non_exonic',
                        'TF_binding_site_variant': 'non_exonic',
                        'regulatory_region_ablation': 'non_exonic',
                        'regulatory_region_amplification': 'non_exonic',
                        'feature_elongation': 'non_exonic',
                        'regulatory_region_variant': 'non_exonic',
                        'feature_truncation': 'non_exonic',
                        'intergenic_variant': 'non_exonic'}


# TODO use pandas dataframes instead of the csv library
def main_program(infile, outfile, filteredfile, famfile, inheritance, familytype, white_list, geneexclusion):
    infile = open(infile, "r")
    outfile = open(outfile, "w")
    filteredfile = open(filteredfile, "w")
    famfile = open(famfile, "r")
    
    names = list()

    # read the gene exclusion list
    # an empty set will not filter out anything, if gene exclusion list is not provided
    genes2exclude = set()
    genes_known = set()

    ### HPO files load
    script_path = os.path.dirname(os.path.abspath(__file__))
    
    # TODO add these filepaths to the config file
    gene_2_HPO_f = os.path.join(script_path, '../../res/gene_2_HPO.p')
    dg_f = os.path.join(script_path, '../../res/v2_ready_graph.pk')
    hpo_dict_file = os.path.join(script_path, '../../res/HPO_gene_assiciation.p')

    gene_2_HPO = pickle.load(open(gene_2_HPO_f, 'rb'))
    graph_nodes, graph_edges = pickle.load(open(dg_f, 'rb'))

    DG = networkx.Graph()
    DG.add_nodes_from(graph_nodes)
    DG.add_edges_from(graph_edges)

    query_dist = 0
    if geneexclusion:
        geneexclusion = open(geneexclusion, "r")
        for gene in geneexclusion:
            gene = gene.rstrip()
            genes2exclude.add(gene)
        geneexclusion.close()

    #fill HPO query terms
    if white_list == None:
        white_list = 'None'
        HPO_query = list()
    else:
        HPO_query = set()
        if os.path.isfile(white_list):
            with open(white_list, 'r') as w:
                HPO_dict = pickle.load(open(hpo_dict_file,'rb'))
                HPO_query = list()
                for line in w:
                    HPO_term = line.rstrip('\n')
                    try:
                        HPO_query.append(HPO_term)
                    except:
                        print('%s not found in database' % (HPO_term))
            HPO_query= list(set(HPO_query))
        else:
            print('The specified HPO list %s is not a valid file' % (white_list))
            print('eDiVA will proceed as without any HPO list')

    # read family relationships
    family = dict()
    for line in famfile:
        if line.startswith('sample'):
            continue
        line = line.rstrip('\n')
        splitline = line.split('\t')
        print(splitline)
        family[splitline[0]] = splitline[1]

    # read all data
    #alldata = list(csv.reader(infile, delimiter="\t"))
    #header = alldata.pop(0)
    
    if inheritance == 'COMPOUND':
        alldata = sorted(alldata, key=itemgetter(0,1)) # sorts alldata according to the first and second element of each row
        alldata = pd.DataFrame(alldata, columns=header)
    else:
        alldata = pd.read_csv(infile, sep="\t", low_memory=False)
        
    #header =[x.replace('#','',1) for x in header]

    out = csv.writer(outfile, delimiter="\t")
    outfiltered = csv.writer(filteredfile, delimiter="\t")

    ############
    # Filter out variants from filtered file where :
    # in case of Xlinked inheritance only variants on the X chromosome are present
    # kick when variant function : synonymous,unknown => line[ExonicFunction(ENSEMBL)]
    # keep when function should be exonic,splicing => line[Function(ENSEMBL)]
    # kick when segmentdup should be 0 => line[SegMentDup]
    # only consider ENSEMBL annotation for performing the first criteria
    #############

    print(family)


    # TODO if the pandas integration works the indices are obsolete 
    index_sample = min([identifycolumns(header,x) for x in family.keys()])#+identifycolumns(header, 'ExAC_SAS') #samples(sampleid>zygosity>DPRef>DPAlt>AF)

    index_MAF1k_AA = identifycolumns(header, 'AA_AF')
    index_MAF1k_AFR = identifycolumns(header, 'AFR_AF')
    index_MAF1k_AMR = identifycolumns(header, 'AMR_AF')
    index_MAF1k_EA = identifycolumns(header, 'EA_AF')
    index_MAF1k_EAS = identifycolumns(header, 'EAS_AF')
    index_MAF1k_EUR = identifycolumns(header, 'EUR_AF')
    index_MAF1k_SAS = identifycolumns(header, 'SAS_AF')
    
    index_MAF_gnomAD_AFR = identifycolumns(header, 'gnomAD_AFR_AF')
    index_MAF_gnomAD_AMR = identifycolumns(header, 'gnomAD_AMR_AF')
    index_MAF_gnomAD_ASJ = identifycolumns(header, 'gnomAD_ASJ_AF')
    index_MAF_gnomAD_EAS = identifycolumns(header, 'gnomAD_EAS_AF')
    index_MAF_gnomAD_FIN = identifycolumns(header, 'gnomAD_FIN_AF')
    index_MAF_gnomAD_NFE = identifycolumns(header, 'gnomAD_NFE_AF')
    index_MAF_gnomAD_OTH = identifycolumns(header, 'gnomAD_OTH_AF')
    index_MAF_gnomAD_SAS = identifycolumns(header, 'gnomAD_SAS_AF')

    index_function = identifycolumns(header, 'Consequence')
    index_varfunction = identifycolumns(header, 'Consequence') # this column doesn't exist anymore TODO remove
    index_segdup = identifycolumns(header, 'SegDupMax') # either SegmentDuplication or SegDupMax, if the first one is chosen we need to choose the max. value if there are multiple
    index_gene = identifycolumns(header, 'SYMBOL')
    index_str = identifycolumns(header, 'SimpleTandemRepeatRegion')
    index_CADD = identifycolumns(header, 'CADD_PHRED')
    index_rank = identifycolumns(header, 'Rank')

    sample_annot_size = 6 # sampleid - DP - REF - ALT - AF - GQ
    print(family.keys())
    for i in range(index_sample,index_sample+sample_annot_size*len(family.keys()),sample_annot_size):
        names.append(header[i])

    ############
    # both the output file headers should be consistent
    #############

    header.append('inheritance')
    header.append('filter')
    header.append('HPO_relatedness')
    header.append('Final_rank')
    index_HPO = identifycolumns(header, 'HPO_relatedness')
    index_final = identifycolumns(header, 'Final_rank')

    outfiltered.writerow(header)
    out.writerow(header)

    # init for compound
    initializer = 0
    processed_HPO_genes=dict()
    compound_gene_storage = []
    
    alldata.apply()
    
    # start reading data
    for line in alldata:
        # init for compound
        if inheritance == 'COMPOUND' and initializer == 0:
            compound_gene_storage = []
            old_gene = re.sub('\(.*?\)','',line[index_gene]) # PKHD1L1(NM_177531:exon75:c.12330+1G>A) transformed to PKHD1L1. Also works for ";" separated multiple annotations
            old_gene_set = set(old_gene.split(';')) # try to remove double occurences of the same gene names
            new_gene = re.sub('\(.*?\)','',line[index_gene])
            new_gene_set = set(new_gene.split(';'))
            initializer = 1

        # read minor allele frequencies
        try:
            # avoiding problems with NAs (NAs should have been filled in the input data)
            MAF1k_AA = [index_MAF1k_AA]
            MAF1k_AFR = line[index_MAF1k_AFR]
            MAF1k_AMR = line[index_MAF1k_AMR]
            MAF1k_EA = line[index_MAF1k_EA]
            MAF1k_EAS = line[index_MAF1k_EAS]
            MAF1k_EUR = line[index_MAF1k_EUR]
            MAF1k_SAS = line[index_MAF1k_SAS]
            
            MAFgnomAD_AFR = line[index_MAF_gnomAD_AFR]
            MAFgnomAD_AMR = line[index_MAF_gnomAD_AMR]
            MAFgnomAD_ASJ = line[index_MAF_gnomAD_ASJ]
            MAFgnomAD_EAS = line[index_MAF_gnomAD_EAS]
            MAFgnomAD_FIN = line[index_MAF_gnomAD_FIN]
            MAFgnomAD_NFE = line[index_MAF_gnomAD_NFE]
            MAFgnomAD_OTH = line[index_MAF_gnomAD_OTH]
            MAFgnomAD_SAS = line[index_MAF_gnomAD_SAS]
            
            MAF = max(float(MAF1k_AA), float(MAF1k_AFR), float(MAF1k_AMR), float(MAF1k_EA), float(MAF1k_EAS), float(MAF1k_EUR), float(MAF1k_SAS), float(MAFgnomAD_AFR), float(MAFgnomAD_AMR), float(MAFgnomAD_ASJ), float(MAFgnomAD_EAS), float(MAFgnomAD_FIN), float(MAFgnomAD_NFE), float(MAFgnomAD_OTH), float(MAFgnomAD_SAS))
        except:
            print('Freq error 1k_AA %s 1k_AFR %s 1k_AMR %s 1k_EA %s 1k_EA %s 1k_EAS %s 1k_EUR %s 1k_SAS %s gnomAD_AFR %s gnomAD_AMR %s gnomAD_ASJ %s gnomAD_EAS %s gnomAD_FIN %s gnomAD_FIN %s gnomAD_NFE %s gnomAD_OTH %s gnomAD_SAS %s')%(MAF1k_AA, MAF1k_AFR, MAF1k_AMR, MAF1k_EA, MAF1k_EAS, MAF1k_EUR, MAF1k_SAS, MAFgnomAD_AFR, MAFgnomAD_AMR, MAFgnomAD_ASJ, MAFgnomAD_EAS, MAFgnomAD_FIN, MAFgnomAD_NFE, MAFgnomAD_OTH, MAFgnomAD_SAS)
            MAF = 0

        try :
            CADD = float(line[index_CADD])
        except:
            CADD = 0

        # read sample names and according zygosity NOW IT's a LIST
        sampledata = line[index_sample:index_sample+sample_annot_size*len(family.keys())]

        if len(sampledata) == 0:
            print(line)
            print(len(line))
            raise
        
        elif sampledata[0] == '0':
            print(len(line))
            print(line)
            raise

        # read simple tandem repeat information
        tandem = line[index_str]

        # check, if gene is on the gene exclusion list.
        genenames = set()
        genecolumn = re.sub('\(.*?\)','',line[index_gene])
        genenames = set(genecolumn.split(';'))

        #part of the distance of the query vs gene - HPO terms
        gene_distances = []
        for gene_id in genenames:
            if gene_id in processed_HPO_genes.keys():
                gene_distances.append(processed_HPO_genes[gene_id])
            else:
                #process ex novo
                #get HPOs related to the gene
                gene_HPO_list = gs.extract_HPO_related_to_gene(gene_2_HPO, gene_id)
                (g_dist,query_dist) = gs.list_distance(DG, HPO_query, gene_HPO_list, query_dist)
                gene_distances.append(g_dist)
                processed_HPO_genes[gene_id] = g_dist
        final_distance = str(max(gene_distances))
        known = final_distance
        judgement = int()

        ###
        # look for de novo variants
        ###
        if inheritance == 'DOMINANT_DENOVO':
            judgement = denovo(sampledata, family,names)
            filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage)
            continue

        ###
        # look for familial dominant variants. (being tolerant for missing values)
        ###
        elif inheritance == 'DOMINANT_INHERITED':
            judgement = dominant(sampledata, family, names)
            sample_annot_size = int(len(sampledata) / len(names))
            total_affected = 0
            for i in range(0, int(sample_annot_size * len(names)), sample_annot_size):
                sam = sampledata[i]
                features = sampledata[i:i+sample_annot_size]
                name = names[int(i / sample_annot_size)]
                total_affected += int(family[name])

            if total_affected == 1 and not familytype == 'SINGLE':
                judgement *= 0 #denovo mutation excluded if not single sample case
            filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage)
            continue

        ###
        # look for recessive variants (be aware of trio and family inheritance)
        ###
        elif inheritance == 'RECESSIVE':
            # TODO single sample case?
            judgement = recessive(sampledata, family, familytype, names)
            filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage)
            continue

        ###
        # look for X linked recessive variants in trios
        ###
        elif inheritance == 'XLINKED':
            index_chromosome = identifycolumns(header, 'Chr')
            
            # skip all variants not located on the X chromosome (both out files will only contain variants located on the X chromosome)
            if line[index_chromosome] == 'X' or line[index_chromosome] == 'x' or line[index_chromosome] == '23':
                pass
            else:
                continue

            judgement = xlinked(sampledata, family, names)
            filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage, familytype)
        
        ###
        # look for compound heterozygous variants
        ###
        elif inheritance == 'COMPOUND' :
            new_gene = re.sub('\(.*?\)', '', line[index_gene])
            new_gene_set = set(new_gene.split(';'))

            if len(old_gene_set - new_gene_set) > 0:

                if  len(names) == 3:
                    comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
                    extension = []
                    pass_ = 0
                    
                    for row in compound_gene_storage:
                        rrow = ','.join(row)

                        if len(compound_gene_storage) == 1:
                            rrow=rrow.replace(',COMPOUND,PASS', 'NOT_compound,filtered')
                            out.writerow(rrow.split(','))
                        else:
                            if comp_judgement==1:
                               outfiltered.writerow(rrow.split(','))
                            else:
                                rrow=rrow.replace(',COMPOUND,PASS', 'NOT_compound,filtered')
                                out.writerow(rrow.split(','))

                elif len(names) == 1:
                    #just check there's more than one and print it
                    if len(compound_gene_storage) == 1:
                        row = compound_gene_storage[0]
                        rrow = ','.join(row).replace(',COMPOUND_SINGLE_SAMPLE,PASS', 'NOT_compound,filtered')
                        out.writerow(rrow.split(','))
                        
                    elif len(compound_gene_storage) > 1:
                        for row in compound_gene_storage:
                            rrow = ','.join(row)
                            outfiltered.writerow(rrow.split(','))
                            out.writerow(rrow.split(','))

                # reset values
                compound_gene_storage = []
                old_gene = new_gene
                old_gene_set = new_gene_set

            if len(names) == 3:
                judgement = compound(sampledata, family, names)
                compound_gene_storage = filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage)
                
            elif len(names) == 1:
                #then use the denovo filter strategy
                judgement = denovo(sampledata, family, names)
                compound_gene_storage = filter_line(judgement, line, MAF, CADD, tandem, 'COMPOUND_SINGLE_SAMPLE', index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage)
                
            continue
        
        ###
        # in case of single sample data with unknown inheritance, skip the inheritance filtering and perform filtering only based on the other criteria
        ###
        elif inheritance == 'UNKNOWN':
            judgement = 1
            filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, HPO_query, compound_gene_storage)
            continue

        else:
            # clean up for last gene
            if inheritance == 'COMPOUND':
                if  len(names) ==3:
                    comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
                    genecolumn = re.sub('\(.*?\)', '', line[index_gene])
                    genenames = set(genecolumn.split(';'))
                    
                    if len(old_gene_set - new_gene_set) > 0:
                        comp_judgement = compoundizer(compound_gene_storage, family, index_sample,names)
                        extension = []
                        pass_ = 0
                        
                        for row in compound_gene_storage:
                            rrow = ','.join(row)
                            if len(compound_gene_storage) == 1:
                                rrow = rrow.replace(',COMPOUND,PASS', ',NOT_compound,filtered')
                                out.writerow(rrow.split(','))
                            else:
                                if comp_judgement==1:
                                    outfiltered.writerow(rrow.split(','))
                                    out.writerow(rrow.split(',') )
                                else:
                                    rrow = rrow.replace(',COMPOUND,PASS', ',NOT_compound,filtered')
                                    out.writerow(rrow.split(',') )
                                    
                elif len(names) == 1:
                    #just check there's more than one and print it
                    if len(compound_gene_storage) == 1:
                        row = compound_gene_storage[0]
                        rrow = ','.join(row).replace(',COMPOUND_SINGLE_SAMPLE,PASS', 'NOT_compound,filtered')
                        out.writerow(rrow.split(','))
                        
                    elif len(compound_gene_storage) > 1:
                        for row in compound_gene_storage:
                            rrow = ','.join(row)
                            outfiltered.writerow(rrow.split(','))
                            out.writerow(rrow.split(','))

    exit(0)


###################################################
# sub routines
###################################################
def compound(sampledata, family,names,debug=False):
    # get samples as data structure
    all_samples = dict()
    samples = sampledata

    check_samples = dict()
    sample_annot_size = int(len(sampledata) / len(names))

    for samp in family.keys():
        check_samples[samp] = 0

    judgement = 0

    for i in range(0,int(sample_annot_size * len(names)), sample_annot_size):
        sam = samples[i]
        features    = samples[i:i+sample_annot_size]
        name        = names[int(i/sample_annot_size)]
        # error catching because of wrong splitting,
        #e.g. 40ACVi>0/1>99>0.333;0.167,40ACVm>0/1>99>0.333;0.167,40ACVp>0/2>99>0.333;0.167
        if len(features) > 1 and name in family:
            zygosity    = features[0]
            refcoverage = features[2].replace('.','0') # could be numeric or .
            altcoverage = features[3].replace('.','0') # could be numeric or .

            if name in check_samples:
                check_samples[name] = 1

            # homo alt is not expected in compound
            if zygosity == '1/1' and family[name] == '1':
                judgement = 0
                break

            # het is good (could be an inherited variant or de novo)
            if zygosity == '0/1' and family[name] == '1':
                judgement = 1
                continue

            # het or hom ref for parents is good
            elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '0':
                judgement = 1
                continue

            # parents must not be hom alt
            elif zygosity == '1/1' and family[name] == '0':
                judgement = 0
                break

            # offspring should have the variant
            elif zygosity == '0/0' and family[name] == '1':
                judgement = 0
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

                # TODO do we need the following try block or can it be removed !?!
                try:
                    int(refcoverage)
                except:
                    exit(0)

                refcoverage = float(refcoverage)
                altcoverage = float(altcoverage)
                coverage = refcoverage + altcoverage

                if coverage == 0:
                    judgement = 0
                    break

                # hom ref
                # poisson for low coverage and percentage for high coverage
                elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                    judgement = 1
                    continue

                # hom alt
                elif poisson.cdf( float(refcoverage), float(coverage/2) ) <= 0.007  and refcoverage / coverage <= 0.05:
                    judgement = 0
                    break

                # coverage too low?
                else:
                    judgement = 0
                    break

            # do not accept missing values for affected individuals
            elif zygosity == './.' and family[name] == '1':
                judgement = 0
                break

    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break
        
    return(judgement)


def compoundizer(variantlist, family, index_sample, names):
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
        sample_annot_size = 6
        sampledata = variantline[index_sample:index_sample+sample_annot_size*len(family.keys())]

        # produce a dictionary to save the zygosities
        zygosities = dict()

        for i in range(0,sample_annot_size * len(names), sample_annot_size):
            sam = sampledata[i]
            features = sampledata[i:i+sample_annot_size]
            name = names[int(i / sample_annot_size)]
            zygosity = features[0]
            refcoverage = features[2].replace('.', '0') # could be numeric or .
            altcoverage = features[3].replace('.', '0') # could be numeric or .
            
            # check if, we are looking at the offspring
            if not list(ticker_dict.keys())[0] == name and not list(ticker_dict.keys())[1] == name:
                continue

            zygosities[name] = zygosity

        # check the entered values for possible compound supporters
        # denovo
        if zygosities[name1] == '0/0' and zygosities[name2] == '0/0':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(0)
            pass

        elif zygosities[name1] == '0/1' and zygosities[name2] == '0/0':
            ticker_dict[name1].append(1)
            ticker_dict[name2].append(0)
            pass

        elif zygosities[name1] == '0/0' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(1)
            pass

        elif zygosities[name1] == '0/1' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(0)
            pass

        elif zygosities[name1] == '0/1' and zygosities[name2] == './.':
            ticker_dict[name1].append(1)
            ticker_dict[name2].append(0)
            pass

        elif zygosities[name1] == './.' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(1)
            pass
        else:
            pass

    # check if there are enough supporters
    name1_sum = sum(ticker_dict[name1])
    name2_sum = sum(ticker_dict[name2])


    if (name1_sum + name2_sum) >= 2 and name1_sum >= 1 and name2_sum >= 1:
        judgement = 1
    else:
        judgement = 0

    return(judgement)


def denovo(sampledata, family, names):
    # get samples as data structure
    samples = sampledata
    judgement = 0
    check_samples = dict()

    # create data structure for completeness check
    for sam in family.keys():
         check_samples[sam] = 0

    # go into the variant data
    sample_annot_size = int(len(sampledata) / len(names))
    for i in range(0,(sample_annot_size * len(names)), sample_annot_size):
        sam = samples[i]
        features = samples[i:i+sample_annot_size]
        name = names[int(i / sample_annot_size)]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            continue

        if len(features) >= 3:
            zygosity = features[0]
            refcoverage = features[2].replace('.', '0') # could be numeric or .
            altcoverage = features[3].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = features[0]
            refcoverage = '.'
            altcoverage = '.'
            
        # check if sample is found in pedigree
        # sample info complete?
        if name in check_samples:
            check_samples[name] = 1

        # heterozygous in affected individual - good
        if zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue

        # hom ref, not affected - good
        elif zygosity == '0/0' and family[name] == '0' :
            if int(altcoverage) <= max(3, (float(altcoverage) + float(refcoverage)) / 10) :
                judgement = 1
            else :
                judgement =0
            continue

        # heterozygous in non-affected - bad
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 0
            break

        # hom ref in affected - bad
        elif zygosity == '0/0' and family[name] == '1':
            judgement = 0
            break

        # homozygous can't be denovo
        elif zygosity == '1/1':
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
                judgement = 0
                break

            # hom ref, non called genotype
            # poisson for low coverage and percentage for high coverage
            # poisson 10 reads (poisson average rate of success = 5) and alt reads = 0 - should get still accepted
            elif (poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007) and (altcoverage / coverage <= 0.05):
                judgement = 1
                continue

            # not necessary to check for hom alt
            # coverage too low?
            else:
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


def dominant(sampledata, family, names):
    # get samples as data structure
    all_samples = dict()
    samples = sampledata
    check_samples = dict()

    for samp in family.keys():
         check_samples[samp] = 0

    judgement = 0

    sample_annot_size = int(len(sampledata) / len(names))
    for i in range(0, int(sample_annot_size * len(names)), sample_annot_size):
        sam = samples[i]
        features = samples[i:i+sample_annot_size]
        name = names[int(i / sample_annot_size)]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            continue

        if len(features) >= 3:
            zygosity = features[0]
            refcoverage = features[2].replace('.', '0') # could be numeric or .
            altcoverage = features[3].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = features[0]
            refcoverage = '.'
            altcoverage = '.'

        if name in check_samples:
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
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                break

            # hom alt
            elif poisson.cdf(float(refcoverage), float(coverage / 2)) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 1
                continue

            # het
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                continue

        else:
            pass

    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return(judgement)


def identifycolumns (header, question):
    try:
        index = header.index(question)
    except:
        exit("ERROR: %s column could not be identified in the annotated data" % question)
    return int(index)


def recessive(sampledata, family, familytype, names):
    # get samples as data structure
    all_samples = dict()
    samples = sampledata
    check_samples = dict()

    for samp in family.keys():
        check_samples[samp] = 0

    judgement = 0

    sample_annot_size = int(len(sampledata) / len(names))
    for i in range(0, int(sample_annot_size * len(names)), sample_annot_size):
        sam = samples[i]
        features = samples[i:i+sample_annot_size]
        name = names[int(i / sample_annot_size)]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            continue

        if len(features) >= 3:
            zygosity = features[0]
            refcoverage = features[2].replace('.', '0') # could be numeric or .
            altcoverage = features[3].replace('.', '0') # could be numeric or .
        else:
            #stick with genotype and the others are empty
            zygosity = features[0]
            refcoverage = '.'
            altcoverage = '.'

        if name in check_samples:
            check_samples[name] = 1

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

        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'FAMILY':
            judgement = 1
            continue

        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'TRIO':
            judgement = 0
            break

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
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue

            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            coverage = refcoverage + altcoverage

            if coverage == 0:
                judgement = 0
                break

            # hom ref
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 0
                break

            # hom alt
            elif poisson.cdf(float(refcoverage), float(coverage / 2)) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 0
                break

            # het, which is OK
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                continue

            # coverage too low?
            else:
                # accept missing values in family interrogations
                if familytype == 'FAMILY':
                    judgement = 1
                    continue
                # do not accept missing values in trio setups
                elif familytype == 'TRIO':
                    judgement = 0
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
    # get samples as data structure
    all_samples = dict()
    samples = sampledata
    check_samples = dict()

    # collect the two parents - one should be het (mother), one whould be hom ref (father)
    # finally should contain two samples, one het & one hom ref
    inheritance_logic = dict()

    for samp in family.keys():
        check_samples[samp] = 0

    judgement = 0

    sample_annot_size = int(len(sampledata) / len(names))
    for i in range(0, int(sample_annot_size * len(names)), sample_annot_size):
        sam = samples[i]
        features = samples[i:i+sample_annot_size]
        name = names[int(i / sample_annot_size)]

        # check if sample is found in pedigree
        try:
            family[name]
        except:
            # if not found, go on to next sample
            print(family)
            continue

        if len(features) >= 3:
            zygosity = features[0]
            refcoverage = features[2].replace('.', '0') # could be numeric or .
            altcoverage = features[3].replace('.', '0') # could be numeric or .

        else:
            #stick with genotype and the others are empty
            zygosity = features[0]
            refcoverage = '.'
            altcoverage = '.'

        if name in check_samples:
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
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) <= 0.007 and altcoverage / coverage <= 0.05:
                inheritance_logic[name] = '0/0'
                judgement = 1
                continue

            # hom alt
            if poisson.cdf(float(refcoverage), float(coverage) / 2) <= 0.007 and refcoverage / coverage <= 0.05:
                inheritance_logic[name] = '1/1'
                judgement = 0
                break

            # het, which is OK
            elif poisson.cdf(float(altcoverage), float(coverage) / 2) >= 0.007 or altcoverage / coverage >= 0.05:
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


def filter_line(judgement, line, MAF, CADD, tandem, inheritance, index_function, index_varfunction, index_segdup, out, outfiltered, genes2exclude, genenames, known, index_rank, hpolist, compound_gene_storage, familytype='TRIO'):
    #deal with each line before writing to output
    #remember continue after
    conditions = dict()
    conditions['RECESSIVE'] = (0.03, -1)
    conditions['DOMINANT_DENOVO'] = (0.01, -1)
    conditions['DOMINANT_INHERITED'] = (0.01, -1)
    conditions['COMPOUND'] = (0.02, -1)
    conditions['XLINKED'] = (0.01, -1)
    conditions['COMPOUND_SINGLE_SAMPLE'] = conditions['DOMINANT_DENOVO']
    conditions['UNKNOWN'] = conditions['DOMINANT_INHERITED']
    
    (MAF_threshold,CADD_threshold) = conditions[inheritance]
    if not familytype == 'TRIO' and inheritance=='XLINKED' :
        cause, result = ('trio_only', 'filtered')

    # exclude gene, if it is on the exclusion list
    elif len(genes2exclude & genenames) > 0:
        cause, result= ('NOT_' + inheritance, 'exclusionlist')

    elif judgement == 1 and not (tandem == 'NA' or tandem == '' or tandem == 'nan'):
        cause, result = (inheritance, 'tandem')
        
    elif judgement == 1 and CADD > 0 and CADD < CADD_threshold :
         cause, result = (inheritance, 'low CADD')
         
    elif judgement == 1 and MAF <= MAF_threshold:
        found_consequences = [variant_consequences[consequence] for consequence in line[index_function].split("&")]
        if ('exonic' in found_consequences or 'splicing' in found_consequences or 'exonic;splicing' in found_consequences):
            if (line[index_function] != 'synonymous_variant' and line[index_function] != 'unknown' and line[index_function] != 'UNKNOWN'):
                if (float(line[index_segdup]) == 0)  :
                    cause, result = (inheritance, 'PASS')
                    if len(hpolist) > 1 and 'NONE' not in hpolist:
                        if float(known) > 0: ######## quota to include it in the output : '0'
                            cause, result = (inheritance, 'PASS')
                        else:
                            cause, result = (inheritance, 'filterd_semantic_distance')
                # e.g. intronic variants fitting the criteria
                else:
                    cause,result = (inheritance, 'filtered seg. dup region')
            else:
                cause, result = (inheritance,'filtered effect')
        else:
            cause, result = (inheritance, 'filtered not exonic')
    
    # fits inheritance, but is too frequent in the population
    elif judgement == 1 and MAF > MAF_threshold:
        cause, result = (inheritance, 'filtered MAF')
    else:
        cause, result = ('NOT_' + inheritance, 'filtered')

    line.append(cause)
    line.append(result)
    line.append(known)
    
    #calc overall score
    final_rank = str((float(line[index_rank]) + float(known)) / 2)
    line.append(final_rank)
    
    if result != 'PASS' :
        out.writerow(line)
    if result == 'PASS':
        if inheritance !='COMPOUND' and inheritance != 'COMPOUND_SINGLE_SAMPLE':
            out.writerow(line)
            outfiltered.writerow(line)
        else:
            compound_gene_storage.append(line)

    return compound_gene_storage


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = 'filter SNPs according to their family distribution')
    parser.add_argument('--infile', type=str, dest='infile', required=True, help='Input annotated ranked file [required]')
    parser.add_argument('--outfile', type=str, dest='outfile', required=True, help='Output file wit all variants. [required]')
    parser.add_argument('--filteredoutfile', type=str, dest='filteredfile', required=True, help='filtered comma separated list of SNPs annotated with inheritance pattern, only reporting the requested variants. [required]')
    parser.add_argument('--family', type=str, dest='famfile', required=True, help='tab separated list of samples annotated with affection status. [required]')
    parser.add_argument('--inheritance', type=str, choices=['DOMINANT_DENOVO', 'DOMINANT_INHERITED', 'RECESSIVE', 'XLINKED', 'COMPOUND', 'UNKNOWN'], dest='inheritance', required=True, help="""choose a inheritance model [required]
    DOMINANT_INHERITED: used for families
    DOMINANT_DENOVO: apply to novel variants seen in the affected individuals
    RECESSIVE: detect recessive, homozygous variants (if trio is specified the script will require that all non-affected are heterozygous)
    XLINKED: used for X linked recessive variants in trios only
    COMPOUND: detect compound heterozygous recessive variants
    UNKNOWN: skip inheritance filtering if not known
    """)
    parser.add_argument('--familytype', type=str, choices=['TRIO', 'FAMILY', 'SINGLE'], dest='familytype', required=True, help="choose if the data you provide is a trio or a larger family")
    parser.add_argument('--geneexclusion', type=str, dest='geneexclusion', required=False, help='[Analysis of DNA sequence variants detected by high-throughput sequencing; DOI: 10.1002/humu.22035]. [required]')
    parser.add_argument('--HPO_list', '--white_list', type=str, dest='white_list', default=None, required=False, help='--HPO_list \t a .txt file with the list of HPO terms describing the disease. It will be used to flag all genes related to the HPO terms. It works with Refseq and UCSC naming.\n')

    args = parser.parse_args()
    main_program(args.infile, args.outfile, args.filteredfile, args.famfile, args.inheritance, args.familytype, args.white_list, args.geneexclusion)
