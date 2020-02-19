#!/usr/bin/env python3
import sys
import os
import subprocess
import argparse
import random
import numpy


def calculate_variant(single_variant, number_mat_reads, number_pat_reads, number_affected_reads, compound_status):
    single_variant = single_variant.split('\t')
    
    ### Calculating a binomial distribution for 0/1 or beta distribution for 1/1 and 0/0
    if number_mat_reads > 0:
        ### 0/1 heteroygous
        mat_01_ref = int(numpy.random.binomial(number_mat_reads, 0.5, 1))
        mat_01_alt = number_mat_reads - mat_01_ref
        ### 1/1 Homozygous alt
        mat_11_alt = int(round(numpy.random.beta(number_mat_reads, 1) * number_mat_reads))
        mat_11_ref = number_mat_reads - mat_11_alt
        ### 0/0 Homozygus ref
        mat_00_alt = int(round(numpy.random.beta(1, number_mat_reads) * number_mat_reads))
        mat_00_ref = number_mat_reads - mat_00_alt
    else:
        mat_01_ref, mat_01_alt, mat_11_alt, mat_11_ref, mat_00_alt, mat_00_ref = (0, 0, 0, 0, 0, 0)
        
    if number_pat_reads > 0:
        pat_01_ref = int(numpy.random.binomial(number_pat_reads, 0.5, 1))
        pat_01_alt = number_pat_reads - pat_01_ref
        pat_11_alt = int(round(numpy.random.beta(number_pat_reads, 1) * number_pat_reads))
        pat_11_ref = number_pat_reads - pat_11_alt
        pat_00_alt = int(round(numpy.random.beta(1, number_pat_reads) * number_pat_reads))
        pat_00_ref = number_pat_reads - pat_00_alt
    else:
        pat_01_ref, pat_01_alt, pat_11_alt, pat_11_ref, pat_00_alt, pat_00_ref = (0, 0, 0, 0, 0, 0)
    
    if number_affected_reads > 0:
        affected_01_ref = int(numpy.random.binomial(number_affected_reads, 0.5, 1))
        affected_01_alt = number_affected_reads - affected_01_ref
        affected_11_alt = int(round(numpy.random.beta(number_affected_reads, 1) * number_affected_reads))
        affected_11_ref = number_affected_reads - affected_11_alt
        affected_00_alt = int(round(numpy.random.beta(1, number_affected_reads) * number_affected_reads))
        affected_00_ref = number_affected_reads - affected_00_alt
    else:
        affected_01_ref, affected_01_alt, affected_11_alt, affected_11_ref, affected_00_alt, affected_00_ref = (0, 0, 0, 0, 0, 0)
    
    ### Inheritance modes
    # M 0/1; P 0/0; C0/1
    if single_variant[6] == 'Dominant':
        mat_info = '0/1' + ':' + str(mat_01_ref) + ',' + str(mat_01_alt) + ':' + str(number_mat_reads) + ':0:.'
        pat_info = '0/0' + ':' + str(pat_00_ref) + ',' + str(pat_00_alt) + ':' + str(number_pat_reads) + ':0:.'
        affected_info = '0/1' + ':' + str(affected_01_ref) + ',' + str(affected_01_alt) + ':' + str(number_affected_reads) + ':0:.'
    
    # M 0/1; P 0/1; C1/1
    elif single_variant[6] == 'Recessive':
        mat_info = '0/1' + ':' + str(mat_01_ref) + ',' + str(mat_01_alt) + ':' + str(number_mat_reads) + ':0:.'
        pat_info = '0/1' + ':' + str(pat_01_ref) + ',' + str(pat_01_alt) + ':' + str(number_pat_reads) + ':0:.'
        affected_info = '1/1' + ':' + str(affected_11_ref) + ',' + str(affected_11_alt) + ':' + str(number_affected_reads) + ':0:.'
    
    # M 0/0; P 0/0; C1/1 or C0/1
    elif single_variant[6] == 'DeNovo':
        mat_info = '0/0' + ':' + str(mat_00_ref) + ',' + str(mat_00_alt) + ':' + str(number_mat_reads) + ':0:.'
        pat_info = '0/0' + ':' + str(pat_00_ref) + ',' + str(pat_00_alt) + ':' + str(number_pat_reads) + ':0:.'
        affected_info = '0/1' + ':' + str(affected_01_ref) + ',' + str(affected_01_alt) + ':' + str(number_affected_reads) + ':0:.'
    
    # M-P 0/1; P-M 0/0; C1/1
    elif single_variant[6] == 'Compound':
        if compound_status == -1:
            mat_info = '0/1' + ':' + str(mat_01_ref) + ',' + str(mat_01_alt) + ':' + str(number_mat_reads) + ':0:.'
            pat_info = '0/0' + ':' + str(pat_00_ref) + ',' + str(pat_00_alt) + ':' + str(number_pat_reads) + ':0:.'
            affected_info = '0/1' + ':' + str(affected_01_ref) + ',' + str(affected_01_alt) + ':' + str(number_affected_reads) + ':0:.'
        else:
            mat_info = '0/0' + ':' + str(mat_00_ref) + ',' + str(mat_00_alt) + ':' + str(number_mat_reads) + ':0:.'
            pat_info = '0/1' + ':' + str(pat_01_ref) + ',' + str(pat_01_alt) + ':' + str(number_pat_reads) + ':0:.'
            affected_info = '0/1' + ':' + str(affected_01_ref) + ',' + str(affected_01_alt) + ':' + str(number_affected_reads) + ':0:.'
    
    else:
        print(single_variant)
    
    ### Writing the variant line for the vcf
    CHR = single_variant[0]
    POS = single_variant[1]
    ID = single_variant[2]
    REF = single_variant[3]
    ALT = single_variant[4]
    INFO = 'SIMULATED=true' 
    INFO = INFO + ';HPO=' + single_variant[-1]
    FORMAT = 'GT:AD:DP:GQ:PL'
    line_for_vcf = CHR + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + '.' + '\t' + '.' + '\t' + INFO + '\t' + FORMAT + '\t' + affected_info + '\t' + pat_info + '\t' + mat_info + '\n'

    return line_for_vcf


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mat_bam', type=str, dest='mat_bam', metavar='alignment.bam', required=True, help='Maternal bam file\n')
    parser.add_argument('--pat_bam', type=str, dest='pat_bam', metavar='alignment.bam', required=True, help='Paternal bam file\n')
    parser.add_argument('--affected_bam', type=str, dest='affected_bam', metavar='alignment.bam', required=True, help='Affected bam file\n')
    parser.add_argument('--fam_vcf', type=str, dest='fam_vcf', metavar='fam.vcf', required=True, help='Combined vcf of trio/family\n')
    parser.add_argument('--varfile', type=str, dest='varfile', metavar='varfile', required=True, help='File containing the simulated variants and inheritance type')
    parser.add_argument('--is_compound', action='store_true', dest='is_compound', required=False, help='Flag to indicate whether the given variants are compound variants. If set the simulator assumes that always two continuous variants belong to the same disease.')
    parser.add_argument('--samtools_path', type=str, dest='samtools_path', metavar='samtools.path', required=False, help='Path to the executable of samtools if not present in the PATH variable')
    parser.add_argument('--bcftools_path', type=str, dest='bcftools_path', metavar='bcftools.path', required=False, help='Path to the executable of bcftools if not present in the PATH variable')
    args = parser.parse_args()
    
    compound_status = -1
    variant_list = open(args.varfile).read().splitlines()
    print(variant_list[0])
    
    tmp = [x  for x in variant_list if not(x.startswith('#'))]
    variant_list = tmp
    variant_position = 'variant_position.bed'
    open_script = open(variant_position, 'w')
    
    for var in variant_list:
        if not(var.startswith('#')):
            # gather info from the variant list and store them in a bed file called variant_posi..
            single_variant = var.split()
            pos_1 = str(int(single_variant[1]) - 1)
            chr_pos = single_variant[0] + '\t' + pos_1 + '\t' + single_variant[1] + '\n'
            open_script.write(chr_pos)
    open_script.close()
    
    simulation_base_name = args.fam_vcf[0:-4]
    ## open reads from mother and father and son
    mat_reads = 'mat_reads'
    pat_reads = 'pat_reads'
    affected_reads = 'affected_reads'
    
    if args.samtools_path:
        subprocess.call('%ssamtools mpileup -l %s %s > mat_reads'%(args.samtools_path, variant_position, args.mat_bam), stdout=subprocess.PIPE, shell=True)
        subprocess.call('%ssamtools mpileup -l %s %s > pat_reads'%(args.samtools_path, variant_position, args.pat_bam), stdout=subprocess.PIPE, shell=True)
        subprocess.call('%ssamtools mpileup -l %s %s > affected_reads'%(args.samtools_path, variant_position, args.affected_bam), stdout=subprocess.PIPE, shell=True)
    else:
        subprocess.call('samtools mpileup -l %s %s > mat_reads'%(variant_position, args.mat_bam), stdout=subprocess.PIPE, shell=True)
        subprocess.call('samtools mpileup -l %s %s > pat_reads'%(variant_position, args.pat_bam), stdout=subprocess.PIPE, shell=True)
        subprocess.call('samtools mpileup -l %s %s > affected_reads'%(variant_position, args.affected_bam), stdout=subprocess.PIPE, shell=True)
    
    vcf_lines = list()
    mat_reads_dict = dict()
    pat_reads_dict = dict()
    affected_reads_dict = dict()
    
    with open(mat_reads, 'r') as mat_reads, open(pat_reads, 'r') as pat_reads, open(affected_reads, 'r') as affected_reads:
        for line in mat_reads:
            fields = line.split('\t')
            mat_reads_dict[';'.join([fields[0], fields[1]])] = max(0, int(fields[3]))
        for line in pat_reads:
            fields = line.split('\t')
            pat_reads_dict[';'.join([fields[0], fields[1]])] = max(0, int(fields[3]))
        for line in affected_reads:
            fields = line.split('\t')
            affected_reads_dict[';'.join([fields[0], fields[1]])] = max(0, int(fields[3]))
        pass
    
    for i in range(len(variant_list)):
        var = variant_list[i]
        fields = var.split('\t')
        number_mat_reads =  mat_reads_dict.get(';'.join([fields[0], fields[1]]) ,0)
        number_pat_reads = pat_reads_dict.get(';'.join([fields[0], fields[1]]) ,0)
        number_affected_reads = affected_reads_dict.get(';'.join([fields[0], fields[1]]) ,0)
        
        try:
            line_for_vcf = calculate_variant(var, number_mat_reads, number_pat_reads, number_affected_reads, compound_status)
        except:
            print(number_pat_reads)
            print(number_affected_reads)
            print(number_mat_reads)
            raise
        
        compound_status *= -1
        vcf_lines.append(line_for_vcf)
    

    
    #######
    ##  Results generation as output:
    ######
    #simulated_vcf = 'simulated.vcf'
    #simulated_vcf = open(simulated_vcf, 'w')
    #fam_vcf = args.fam_vcf
    #fam_vcf = open(fam_vcf, 'r')
    #for line in fam_vcf:
    #    simulated_vcf.write(line)
        
    #Appending variant lines at the end of the VCF file
    count = 1
    
    if args.is_compound:
        for first_entry, second_entry in zip(vcf_lines[0::2], vcf_lines[1::2]):
            print(first_entry)
            print(second_entry)
            fam_vcf = args.fam_vcf
            simulated_vcf = 'simulated_' + str(count) + '.vcf'
            simulated_vcf = open(simulated_vcf, 'w')
            fam_vcf = open(fam_vcf, 'r')
            
            for line in fam_vcf:
                simulated_vcf.write(line)
            
            simulated_vcf.write(first_entry)
            simulated_vcf.write(second_entry)
            
            simulated_vcf.close()
            fam_vcf.close()
        
            if args.bcftools_path:
                subprocess.call(['%sbcftools sort %s > %s'%(args.bcftools_path, 'simulated_' + str(count) + '.vcf', simulation_base_name.split('/')[-1] + '_simulated_variant_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
            else:
                subprocess.call(['bcftools sort %s > %s'%('simulated_' + str(count) + '.vcf', simulation_base_name.split('/')[-1] + '_simulated_variant_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
        
            subprocess.call(['gzip -9 %s'%(simulation_base_name.split('/')[-1] + '_simulated_variant_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
            subprocess.call(['rm %s'%('simulated_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
        
    else:
        for line_for_vcf in vcf_lines:
            print(line_for_vcf)
            fam_vcf = args.fam_vcf
            simulated_vcf = 'simulated_' + str(count) + '.vcf'
            simulated_vcf = open(simulated_vcf, 'w')
            fam_vcf = open(fam_vcf, 'r')
            
            for line in fam_vcf:
                simulated_vcf.write(line)
            
            simulated_vcf.write(line_for_vcf)
            
            simulated_vcf.close()
            fam_vcf.close()
        
            if args.bcftools_path:
                subprocess.call(['%sbcftools sort %s > %s'%(args.bcftools_path, 'simulated_' + str(count) + '.vcf', simulation_base_name.split('/')[-1] + '_simulated_variant_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
            else:
                subprocess.call(['bcftools sort %s > %s'%('simulated_' + str(count) + '.vcf', simulation_base_name.split('/')[-1] + '_simulated_variant_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
        
            subprocess.call(['gzip -9 %s'%(simulation_base_name.split('/')[-1] + '_simulated_variant_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
            subprocess.call(['rm %s'%('simulated_' + str(count) + '.vcf')], stdout=subprocess.PIPE, shell=True)#!/usr/bin/env python
