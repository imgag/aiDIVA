import subprocess
import argparse
import tempfile
import time
import os


def call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, basic=False, expanded=False, num_cores=1):
    # the path to the executable
    vep_command = vep_annotation_dict["vep"] + "/" + "vep "

    bed_annotation = vep_annotation_dict["custom"]["bed-files"]
    bigwig_annotation = vep_annotation_dict["custom"]["bigwig-files"]
    #vcf_annotation = vep_annotation_dict["custom"]["vcf-files"]

    # set the correct paths to the needed perl modules
    if "PERL5LIB" in os.environ:
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/:" + os.environ["PERL5LIB"]
    else:
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/"

    # add essential parameters
    vep_command = vep_command + "--species homo_sapiens --assembly GRCh37" + " "
    vep_command = vep_command + "--offline" + " "
    vep_command = vep_command + "--cache" + " "
    vep_command = vep_command + "--dir_cache " + vep_annotation_dict["db"] + vep_annotation_dict["cache-path"] + " "
    vep_command = vep_command + "--gencode_basic" + " "
    vep_command = vep_command + "--symbol" + " "
    vep_command = vep_command + "--biotype" + " "
    vep_command = vep_command + "--variant_class" + " "

    if not expanded:
        # allele frequencies to include
        vep_command = vep_command + "--max_af" + " "
        # the following AF annotations could be dropped since we only need the max AF
        #vep_command = vep_command + "--af" + " "
        #vep_command = vep_command + "--af_1kg" + " "
        #vep_command = vep_command + "--af_esp" + " "
        #vep_command = vep_command + "--af_gnomad" + " "
        vep_command = vep_command + "--custom " + vep_annotation_dict["db"] + bed_annotation["simpleRepeat"]["file"]  + "," + "simpleRepeat" + ",bed," + bed_annotation["simpleRepeat"]["method"] + ",0" + " "
        vep_command = vep_command + "--custom " + vep_annotation_dict["db"] + bed_annotation["oe_lof"]["file"]  + "," + "oe_lof" + ",bed," + bed_annotation["oe_lof"]["method"] + ",0" + " "

    # vep plugins to use
    if not basic:
        vep_command = vep_command + "--sift s" + " "
        vep_command = vep_command + "--polyphen s" + " "
        #vep_command = vep_command + "--plugin Condel," + vep_annotation_dict["condel"] + ",s" + " "
        vep_command = vep_command + "--plugin CADD," + vep_annotation_dict["db"] + vep_annotation_dict["cadd-snps"] + " " #"," + vep_annotation_dict["db"] + vep_annotation_dict["cadd-indel"] + " "
        vep_command = vep_command + "--plugin REVEL," + vep_annotation_dict["db"] + vep_annotation_dict["revel"] + " "
        #vep_command = vep_command + "--plugin dbNSFP," + vep_annotation_dict["dbNSFP"] + ",MutationAssessor_score,Eigen-raw,Eigen-phred" + " "

        if bed_annotation:
            for key in bed_annotation:
                if (key != "simpleRepeat") and (key != "oe_lof"):
                    vep_command = vep_command + "--custom " + vep_annotation_dict["db"] + bed_annotation[key]["file"]  + "," + key + ",bed," + bed_annotation[key]["method"] + ",0" + " "

        # does not work (seems to be a problem with VEP)
        #if vcf_annotation:
        #    for key in vcf_annotation:
        #        vep_command = vep_command + "--custom " + vcf_annotation[key]["file"] + "," + vcf_annotation[key]["prefix"] + ",vcf," + vcf_annotation[key]["method"] + ",0," + key + " "

        if bigwig_annotation:
            for key in bigwig_annotation:
                vep_command = vep_command + "--custom " + vep_annotation_dict["db"] + bigwig_annotation[key]["file"] + "," + key + ",bigwig," + bigwig_annotation[key]["method"] + ",0" + " "

    vep_command = vep_command + "-i " + input_vcf_file + " "
    vep_command = vep_command + "-o " + output_vcf_file + " "
    vep_command = vep_command + "--fork " + str(num_cores) + " "
    vep_command = vep_command + "--vcf" + " "
    vep_command = vep_command + "--no_stats" + " "
    vep_command = vep_command + "--force_overwrite"

    subprocess.run(vep_command, shell=True, check=True)
    print("The annotated VCF is saved as %s" % (output_vcf_file))


def annotate_consequence_information(input_vcf_file, output_vcf_file, vep_annotation_dict, num_cores=1):
	# the path to the executable
    vep_command = vep_annotation_dict["vep"] + "/" + "vep "

    # set the correct paths to the needed perl modules
    if "PERL5LIB" in os.environ:
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/:" + os.environ["PERL5LIB"]
    else:
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/"

    # add essential parameters
    vep_command = vep_command + "--species homo_sapiens --assembly GRCh37" + " "
    vep_command = vep_command + "--offline" + " "
    vep_command = vep_command + "--cache" + " "
    vep_command = vep_command + "--dir_cache " + vep_annotation_dict["db"] + vep_annotation_dict["cache-path"] + " "
    vep_command = vep_command + "--gencode_basic" + " "


    vep_command = vep_command + "-i " + input_vcf_file + " "
    vep_command = vep_command + "-o " + output_vcf_file + " "
    vep_command = vep_command + "--fields Consequence "
    vep_command = vep_command + "--vcf_info_field CONS "
    vep_command = vep_command + "--fork " + str(num_cores) + " "
    vep_command = vep_command + "--vcf" + " "
    vep_command = vep_command + "--no_stats" + " "
    vep_command = vep_command + "--force_overwrite"

    subprocess.run(vep_command, shell=True, check=True)
    print("The annotated VCF is saved as %s" % (output_vcf_file))


def annotate_from_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, num_cores):
    tmp = tempfile.NamedTemporaryFile(mode="w+b", suffix=".config", delete=False)

    vcf_annotation = vep_annotation_dict["custom"]["vcf-files"]
    command = vep_annotation_dict["ngs-bits"] + "/" + "VcfAnnotateFromVcf -config_file " + tmp.name + " -in " + input_vcf_file + " -out " + output_vcf_file + " -threads " + str(num_cores)

    try:
        tmp.write(str(vep_annotation_dict["db"] + vcf_annotation["CONDEL"]["file"] + "\t\tCONDEL\t\ttrue\n").encode())
        tmp.write(str(vep_annotation_dict["db"] + vcf_annotation["EIGEN_PHRED"]["file"] + "\t\tEIGEN_PHRED\t\ttrue\n").encode())
        tmp.write(str(vep_annotation_dict["db"] + vcf_annotation["FATHMM_XF"]["file"] + "\t\tFATHMM_XF\t\ttrue\n").encode())
        tmp.write(str(vep_annotation_dict["db"] + vcf_annotation["MutationAssessor"]["file"] + "\t\tMutationAssessor\t\ttrue\n").encode())
        tmp.write(str(vep_annotation_dict["db"] + vcf_annotation["gnomAD"]["file"] + "\tgnomAD\tAN,Hom\t\ttrue\n").encode())
        #tmp.write(str(vep_annotation_dict["db"] + vcf_annotation["CAPICE"]["file"] + "\t\tCAPICE\t\ttrue\n").encode())
        tmp.close()

        subprocess.run(command, shell=True, check=True)
    finally:
        os.remove(tmp.name)


def sort_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict):
    command = vep_annotation_dict["ngs-bits"] + "/" + "VcfSort -in " + input_vcf_file + " -out " + output_vcf_file
    subprocess.run(command, shell=True, check=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Annotation with VEP")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="data.vcf", required=True, help="VCF file containing the data, you want to annotate with VEP\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.vcf", required=True, help="Specifies the extended output file\n")
    args = parser.parse_args()

    input_vcf_file = args.in_data
    output_vcf_file = args.out_data

    ## TODO: change to readfrom yaml file
    vep_annotation_dict2 = {"vep": "/mnt/users/ahbranl1/data_vep/ensembl-vep-release-98.3/vep",
                           "cache-path": "/mnt/users/ahbranl1/data_vep/ensembl-vep-cache/cache",
                           "num-threads": 10,
                           "plugin-path": "/mnt/users/ahboced1/data_vep/plugins",
                           "condel": "/mnt/users/ahboced1/Tools/vep_data/plugins/config/Condel/config",
                           "cadd-snps": "/mnt/share/data/dbs/CADD/whole_genome_SNVs.tsv.gz",
                           "cadd-indel": "/mnt/share/data/dbs/CADD/InDels.tsv.gz",
                           "revel": "/mnt/share/data/dbs/REVEL/revel_all_chromosomes.tsv.gz",
                           "dbNSFP": "/mnt/users/ahbranl1/data_vep/dbNSFP/dbNSFP_hg19_3.5.gz",
                           "custom": {"bed-files": {"simpleRepeat": {"file": "/mnt/users/ahboced1/databases/hg19/simpleRepeats.bedGraph.gz", "method": "overlap"},
                           "segmentDuplication": {"file": "/mnt/users/ahboced1/databases/hg19/segmentDuplication.bedGraph.gz", "method": "overlap"},
                           "ABB_SCORE": {"file": "/mnt/users/ahboced1/databases/hg19/abb_score.bedGraph.gz", "method": "exact"}}}}

    call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict2, True)
