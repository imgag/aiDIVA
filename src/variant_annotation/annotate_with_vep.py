import subprocess
import argparse


def call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, additional_annotation_dict, only_basic=False):
    # the path to the executable
    vep_command = vep_annotation_dict["vep"] + " "
    
    # add essential parameters
    vep_command = vep_command + "--offline" + " "
    vep_command = vep_command + "--cache" + " "
    vep_command = vep_command + "--dir_cache " + vep_annotation_dict["cache-path"] + " "
    vep_command = vep_command + "--sift s" + " "
    vep_command = vep_command + "--polyphen s" + " "
    vep_command = vep_command + "--symbol" + " "
    vep_command = vep_command + "--biotype" + " "
    vep_command = vep_command + "--variant_class" + " "
    
    # allele frequencies to include
    vep_command = vep_command + "--af" + " "
    vep_command = vep_command + "--af_1kg" + " "
    vep_command = vep_command + "--af_esp" + " "
    vep_command = vep_command + "--af_gnomad" + " "
    
    # vep plugins to use
    if not only_basic:
        vep_command = vep_command + "--dir_plugin " + vep_annotation_dict["plugin-path"] + " "
        vep_command = vep_command + "--plugin Condel," + vep_annotation_dict["condel"] + ",s" + " "
        vep_command = vep_command + "--plugin CADD," + vep_annotation_dict["cadd-snps"] + "," + vep_annotation_dict["cadd-indel"] + " "
        vep_command = vep_command + "--plugin REVEL," + vep_annotation_dict["revel"] + " "
        vep_command = vep_command + "--plugin dbNSFP," + vep_annotation_dict["dbNSFP"] + ",MutationAssessor_score,Eigen-raw,Eigen-phred" + " "
    
        # custom annotations that should be used during the annotation
        # TODO does not work at the moment need to do it manually
        #print(additional_annotation_dict)
        #for key in additional_annotation_dict:
        #    if additional_annotation_dict[key].endswith(".bw"):
        #        vep_command = vep_command + "--custom " + additional_annotation_dict[key] + "," + key + ",bigwig,exact,0" + " "
        #    else:
        #        print("ERROR: Given file seems to be no bigwig file")
    
    vep_command = vep_command + "-i " + input_vcf_file + " "
    vep_command = vep_command + "-o " + output_vcf_file + " "
    vep_command = vep_command + "--fork " + str(vep_annotation_dict["num-threads"]) + " "
    vep_command = vep_command + "--vcf"
    
    print("Starting VEP annotation ...")
    
    subprocess.run(vep_command, shell=True, check=True)
    
    print("Finished VEP annotation!")
    print("The annotated VCF is saved as %s" % (output_vcf_file))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description = 'AIdiva -- Annotation with VEP')
    parser.add_argument('--in_data', type=str, dest='in_data', metavar='data.vcf', required=True, help='VCF file containing the data, you want to annotate with VEP\n')
    parser.add_argument('--out_data', type=str, dest='out_data', metavar='out.vcf', required=True, help='Specifies the extended output file\n')
    args = parser.parse_args()
    
    input_vcf_file = args.in_data
    output_vcf_file = args.out_data
    
    vep_annotation_dict = {"vep": "/mnt/users/ahbranl1/data_vep/ensembl-vep-release-98.3/vep",
                           "cache-path": "/mnt/users/ahbranl1/data_vep/cache",
                           "plugin-path": "/mnt/users/ahboced1/data_vep/plugins",
                           "condel": "/mnt/users/ahboced1/Tools/vep_data/plugins/config/Condel/config",
                           "cadd-snps": "/mnt/share/data/dbs/CADD/whole_genome_SNVs.tsv.gz",
                           "cadd-indel": "/mnt/share/data/dbs/CADD/InDels.tsv.gz",
                           "revel": "/mnt/share/data/dbs/REVEl/revel_all_chromosomes.tsv.gz",
                           "dbNSFP": "/mnt/users/ahbranl1/data_vep/dbNSFP/dbNSFP_hg19_3.5.gz"}
    
    additional_annotation_dict = {"phyloP46mammal": "/home/dominic/Masterarbeit/databases/phyloP46way.placentalMammal.complete.onlyChrom.bw",
                                  "phyloP46primate": "/home/dominic/Masterarbeit/databases/phyloP46way.primate.complete.onlyChrom.bw",
                                  "phyloP46vertebrate": "/home/dominic/Masterarbeit/databases/phyloP46way.vertebrate.complete.onlyChrom.bw",
                                  "pastCons46mammal": "/home/dominic/Masterarbeit/databases/phastCons46way.placentalMammal.complete.onlyChrom.bw",
                                  "pastCons46primate": "/home/dominic/Masterarbeit/databases/phastCons46way.primate.complete.onlyChrom.bw",
                                  "pastCons46vertebrate": "/home/dominic/Masterarbeit/databases/phastCons46way.vertebrate.complete.onlyChrom.bw"}
    
    call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, additional_annotation_dict)



