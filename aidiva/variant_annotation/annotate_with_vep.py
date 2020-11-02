import subprocess
import argparse
import tempfile
import time
import os


def call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, basic=False, num_cores=1):
    # the path to the executable
    database_path = "/mnt/storage1/share/data/dbs/"
    # database_path = os.path.dirname(__file__) + "/../../annotation_resources/

    tool_path = "/mnt/storage1/share/opt/"
    # tool_path = os.path.dirname(__file__) + "/../../tools/"

    vep_command = tool_path + vep_annotation_dict["vep"] + " "

    bed_annotation = vep_annotation_dict["custom"]["bed-files"]
    bigwig_annotation = vep_annotation_dict["custom"]["bigwig-files"]
    #print("ENV:", os.environ["PERL5LIB"])

    # set the correct paths to the needed perl modules
    os.environ["PERL5LIB"] = tool_path + "/ensembl-vep-release-100.3/Bio/:" + tool_path + "/ensembl-vep-release-100.3/cpan/lib/perl5/:" + os.environ["PERL5LIB"]

    #print("ENV_mod:", os.environ["PERL5LIB"])

    # add essential parameters
    vep_command = vep_command + "--offline" + " "
    vep_command = vep_command + "--cache" + " "
    vep_command = vep_command + "--dir_cache " + database_path + vep_annotation_dict["cache-path"] + " "
    vep_command = vep_command + "--sift s" + " "
    vep_command = vep_command + "--polyphen s" + " "
    vep_command = vep_command + "--symbol" + " "
    vep_command = vep_command + "--biotype" + " "
    vep_command = vep_command + "--variant_class" + " "

    # allele frequencies to include
    vep_command = vep_command + "--max_af" + " "
    # the following AF annotations could be dropped since we only need the max AF
    vep_command = vep_command + "--af" + " "
    vep_command = vep_command + "--af_1kg" + " "
    vep_command = vep_command + "--af_esp" + " "
    vep_command = vep_command + "--af_gnomad" + " "

    # vep plugins to use
    if not basic:
        #vep_command = vep_command + "--plugin Condel," + vep_annotation_dict["condel"] + ",s" + " "
        vep_command = vep_command + "--plugin CADD," + database_path + vep_annotation_dict["cadd-snps"] + "," + database_path + vep_annotation_dict["cadd-indel"] + " "
        vep_command = vep_command + "--plugin REVEL," + database_path + vep_annotation_dict["revel"] + " "
        #vep_command = vep_command + "--plugin dbNSFP," + vep_annotation_dict["dbNSFP"] + ",MutationAssessor_score,Eigen-raw,Eigen-phred" + " "

        #vcf_annotation = vep_annotation_dict["custom"]["vcf-files"]

        if bed_annotation:
            for key in bed_annotation:
                vep_command = vep_command + "--custom " + database_path + bed_annotation[key]["file"]  + "," + key + ",bed," + bed_annotation[key]["method"] + ",0" + " "

        # does not work (seems to be a problem with VEP)
        #if vcf_annotation:
        #    for key in vcf_annotation:
        #        vep_command = vep_command + "--custom " + vcf_annotation[key]["file"] + "," + vcf_annotation[key]["prefix"] + ",vcf," + vcf_annotation[key]["method"] + ",0," + key + " "

        if bigwig_annotation:
            for key in bigwig_annotation:
                vep_command = vep_command + "--custom " + database_path + bigwig_annotation[key]["file"] + "," + key + ",bigwig," + bigwig_annotation[key]["method"] + ",0" + " "
    else:
        vep_command = vep_command + "--custom " + database_path + bed_annotation["simpleRepeat"]["file"]  + "," + "simpleRepeat" + ",bed," + bed_annotation["simpleRepeat"]["method"] + ",0" + " "

    vep_command = vep_command + "-i " + input_vcf_file + " "
    vep_command = vep_command + "-o " + output_vcf_file + " "
    vep_command = vep_command + "--fork " + str(vep_annotation_dict["num-threads"]) + " "
    vep_command = vep_command + "--vcf" + " "
    vep_command = vep_command + "--no_stats" + " "
    vep_command = vep_command + "--force_overwrite"

    subprocess.run(vep_command, shell=True, check=True)
    print("The annotated VCF is saved as %s" % (output_vcf_file))


def annotate_from_vcf(input_vcf_file, output_vcf_file, num_cores=1):
    tmp = tempfile.NamedTemporaryFile(mode="w+b", suffix=".config", delete=False)

    database_path = "/mnt/storage1/share/data/dbs/"
    # database_path = os.path.dirname(__file__) + "/../../annotation_resources/

    tool_path = "/mnt/storage1/share/opt/"
    # tool_path = os.path.dirname(__file__) + "/../../tools/"

    command = tool_path + "ngs-bits-current/VcfAnnotateFromVCF -config_file " + tmp.name + " -in " + input_vcf_file + " -out " + output_vcf_file + " -threads " + str(num_cores)

    try:
        tmp.write(str(database_path + "Condel/hg19_precomputed_Condel.vcf.gz\t\tCONDEL\t\ttrue\n").encode())
        tmp.write(str(database_path + "Eigen/hg19_Eigen-phred_coding_chrom1-22.vcf.gz\t\tEIGEN_PHRED\t\ttrue\n").encode())
        tmp.write(str(database_path + "fathmm-XF/hg19_fathmm_xf_coding.vcf.gz\t\tFATHMM_XF\t\ttrue\n").encode())
        tmp.write(str(database_path + "MutationAssessor/hg19_precomputed_MutationAssessor.vcf.gz\t\tMutationAssessor\t\ttrue\n").encode())
        tmp.close()

        subprocess.run(command, shell=True, check=True)
    finally:
        os.remove(tmp.name)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Annotation with VEP")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="data.vcf", required=True, help="VCF file containing the data, you want to annotate with VEP\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.vcf", required=True, help="Specifies the extended output file\n")
    args = parser.parse_args()

    input_vcf_file = args.in_data
    output_vcf_file = args.out_data

    ## TODO: change to readfrom yaml file
    vep_annotation_dict = {"vep": "/mnt/users/ahbranl1/data_vep/ensembl-vep-release-98.3/vep",
                           "cache-path": "/mnt/users/ahbranl1/data_vep/ensembl-vep-cache/cache",
                           "num-threads": 10,
                           "plugin-path": "/mnt/users/ahboced1/data_vep/plugins",
                           "condel": "/mnt/users/ahboced1/Tools/vep_data/plugins/config/Condel/config",
                           "cadd-snps": "/mnt/share/data/dbs/CADD/whole_genome_SNVs.tsv.gz",
                           "cadd-indel": "/mnt/share/data/dbs/CADD/InDels.tsv.gz",
                           "revel": "/mnt/share/data/dbs/REVEl/revel_all_chromosomes.tsv.gz",
                           "dbNSFP": "/mnt/users/ahbranl1/data_vep/dbNSFP/dbNSFP_hg19_3.5.gz",
                           "custom": {"bed-files": {"simpleRepeat": {"file": "/mnt/users/ahboced1/databases/hg19/simpleRepeats.bedGraph.gz", "method": "overlap"},
                           "segmentDuplication": {"file": "/mnt/users/ahboced1/databases/hg19/segmentDuplication.bedGraph.gz", "method": "overlap"},
                           "ABB_SCORE": {"file": "/mnt/users/ahboced1/databases/hg19/abb_score.bedGraph.gz", "method": "exact"}}}}

    call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, True)
