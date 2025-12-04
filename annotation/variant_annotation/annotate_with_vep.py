import subprocess
import argparse
import tempfile
import os
import logging
import yaml
import sys


logger = logging.getLogger(__name__)


def call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_dict, annotation_dict, assembly_build="GRCh38", basic=False, expanded=False, num_cores=1):
    vep_mode = vep_dict["vep-mode"]

    if vep_mode == "LOCAL":
        logger.info("Using local VEP installation.")

        # the path to the executable
        vep_command = f"{vep_dict['vep-local']['vep']}/vep"

        # set the correct paths to the needed perl modules
        if "PERL5LIB" in os.environ:
            os.environ["PERL5LIB"] = f"{os.environ['PERL5LIB']}:{vep_dict['vep-local']['vep']}/Bio/:{vep_dict['vep-local']['vep-plugin-path']}:{vep_dict['vep-local']['vep-cpan']}"

        else:
            os.environ["PERL5LIB"] = f"{vep_dict['vep-local']['vep']}/Bio/:{vep_dict['vep-local']['vep-plugin-path']}:{vep_dict['vep-local']['vep-cpan']}"

        plugin_path = vep_dict['vep-local']['vep-plugin-path'] + "/"

    elif vep_mode == "APPTAINER":
        logger.info("Using containerized VEP with Apptainer.")
        vep_command = f"apptainer exec --bind {annotation_dict['resources-base-path']}:{annotation_dict['resources-base-path']}:ro,{input_vcf_file}:{input_vcf_file}:ro"
        vep_command = f"{vep_command} {vep_dict['vep-container']['vep-image']} vep"

    elif vep_mode == "SINGULARITY":
        logger.info("Using containerized VEP with Singularity.")
        vep_command = f"singularity exec --bind {annotation_dict['resources-base-path']}:{annotation_dict['resources-base-path']}:ro,{input_vcf_file}:{input_vcf_file}:ro"
        vep_command = f"{vep_command} {vep_dict['vep-container']['vep-image']} vep"

    #elif vep_mode == "DOCKER":
    #    logger.info("Using containerized VEP with Docker.")
    #    vep_command = f"docker run {vep_dict['vep-container']['vep-docker']} vep"

    else:
        sys.exit()

    cache_path = annotation_dict['vep-cache'] + "/"

    # add essential parameters
    vep_command = f"{vep_command} --species homo_sapiens --assembly {assembly_build} "
    vep_command = f"{vep_command} --offline"
    vep_command = f"{vep_command} --cache"
    vep_command = f"{vep_command} --dir_cache {cache_path}"

    if vep_mode == "LOCAL":
        vep_command = f"{vep_command} --dir_plugins {plugin_path}"

    vep_command = f"{vep_command} --gencode_basic"
    vep_command = f"{vep_command} --symbol"
    vep_command = f"{vep_command} --biotype"
    vep_command = f"{vep_command} --variant_class"
    vep_command = f"{vep_command} --af_gnomadg"

    # vep plugins to use
    if not basic:
        vep_command = f"{vep_command} --sift s"
        vep_command = f"{vep_command} --polyphen s"

        # specify all the Plugins that VEP should use
        vep_command = f"{vep_command} --plugin AlphaMissense,file={annotation_dict['vep-plugin-files']['AlphaMissense']}"
        vep_command = f"{vep_command} --plugin dbNSFP,{annotation_dict['vep-plugin-files']['dbNSFP']},transcript_match=1,MutationAssessor_score"

    if not expanded:
        vep_command = f"{vep_command} --plugin SpliceAI,snv={annotation_dict['vcf-files']['SpliceAI-SNV']},indel={annotation_dict['vcf-files']['SpliceAI-InDel']}"

    # all vcf annotations
    if not basic:
        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['CONDEL']},short_name=VEP_CONDEL,format=vcf,type=exact,fields=CONDEL"
        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['EIGEN_PHRED']},short_name=VEP_EIGEN,format=vcf,type=exact,fields=EIGEN_PHRED"
        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['FATHMM_XF']},short_name=VEP_FATHMM_XF,format=vcf,type=exact,fields=FATHMM_XF"

        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['CAPICE']},short_name=VEP_CAPICE,format=vcf,type=exact,fields=CAPICE"
        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['CADD']},short_name=VEP_CADD,format=vcf,type=exact,fields=CADD"
        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['REVEL']},short_name=VEP_REVEL,format=vcf,type=exact,fields=REVEL"

        vep_command = f"{vep_command} --custom file={annotation_dict['bigwig-files']['phyloP_primate']},short_name=phyloP_primate,format=bigwig,type=exact,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bigwig-files']['phyloP_mammal']},short_name=phyloP_mammal,format=bigwig,type=exact,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bigwig-files']['phyloP_vertebrate']},short_name=phyloP_vertebrate,format=bigwig,type=exact,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bigwig-files']['phastCons_primate']},short_name=phastCons_primate,format=bigwig,type=exact,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bigwig-files']['phastCons_mammal']},short_name=phastCons_mammal,format=bigwig,type=exact,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bigwig-files']['phastCons_vertebrate']},short_name=phastCons_vertebrate,format=bigwig,type=exact,coords=0"

    if not expanded:
        vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['gnomAD']},short_name=gnomAD,format=vcf,type=exact,fields=AN%Hom"

        # HGMD needs a valid license, therefore we check if the file exists otherwise this annotation is skipped
        if os.path.isfile(f"{annotation_dict['vcf-files']['hgmd']}"):
            vep_command = f"{vep_command} --custom file={annotation_dict['vcf-files']['hgmd']},short_name=HGMD,format=vcf,type=exact,fields=CLASS%RANKSCORE"

        else:
            logger.warning("HGMD file is not found! Skip HGMD annotation!")

        vep_command = f"{vep_command} --custom file={annotation_dict['bed-files']['segmentDuplication']},short_name=SegDup,format=bed,type=overlap,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bed-files']['simpleRepeat']},short_name=SimpleRepeats,format=bed,type=overlap,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bed-files']['oe_lof']},short_name=oe_lof,format=bed,type=overlap,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['bed-files']['repeatMasker']},short_name=REPEATMASKER,format=bed,type=overlap,coords=0"
        vep_command = f"{vep_command} --custom file={annotation_dict['low-confidence']},short_name=low_conf_region,format=bed,type=overlap,coords=0"

        if os.path.isfile(f"{annotation_dict['bed-files']['omim']}"):
            vep_command = f"{vep_command} --custom file={annotation_dict['bed-files']['omim']},short_name=OMIM,format=bed,type=overlap,coords=0"

        else:
            logger.warning("OMIM file is not found! Skip OMIM annotation!")

    vep_command = f"{vep_command} -i " + input_vcf_file + " "
    vep_command = f"{vep_command} -o " + output_vcf_file + " "

    # the developers of VEP recommend to use 4 threads
    if num_cores >= 4:
        vep_command = f"{vep_command} --fork 4 "

    elif num_cores == 2:
        vep_command = f"{vep_command} --fork 2 "

    vep_command = f"{vep_command} --format vcf" + " " # we need this to prevent vep from not working if the VCF file has no variant entries
    vep_command = f"{vep_command} --vcf" + " "
    vep_command = f"{vep_command} --no_stats" + " "
    vep_command = f"{vep_command} --force_overwrite"

    print(vep_command)
    subprocess.run(vep_command, shell=True, check=True)
    logger.debug("The VEP annotated VCF is saved as %s" % (output_vcf_file))


## TODO do we need to add the filtered_folder as additional bind path?
def call_vep_and_annotate_consequence_information(input_vcf_file, output_vcf_file, vep_dict, annotation_dict, assembly_build="GRCh38", num_cores=1):
    vep_mode = vep_dict["vep-mode"]

    if vep_mode == "LOCAL":
        logger.info("Using local VEP installation for consequence annotation.")
        # the path to the executable
        vep_command = f"{vep_dict['vep-local']['vep']}/vep"

        # set the correct paths to the needed perl modules
        if "PERL5LIB" in os.environ:
            os.environ["PERL5LIB"] = f"{os.environ['PERL5LIB']}:{vep_dict['vep-local']['vep']}/Bio/:{vep_dict['vep-local']['vep-cpan']}"

        else:
            os.environ["PERL5LIB"] = f"{vep_dict['vep-local']['vep']}/Bio/:{vep_dict['vep-local']['vep-plugin-path']}:{vep_dict['vep-local']['vep-cpan']}"

    elif vep_mode == "APPTAINER":
        logger.info("Using containerized VEP with Apptainer for consequence annotation.")
        vep_command = f"apptainer exec --bind {annotation_dict['resources-base-path']}:{annotation_dict['resources-base-path']}:ro,{input_vcf_file}:{input_vcf_file}:ro"
        vep_command = f"{vep_command} {vep_dict['vep-container']['vep-image']} vep"

    elif vep_mode == "SINGULARITY":
        logger.info("Using containerized VEP with Singularity for consequence annotation.")
        vep_command = f"singularity exec --bind {annotation_dict['resources-base-path']}:{annotation_dict['resources-base-path']}:ro,{input_vcf_file}:{input_vcf_file}:ro"
        vep_command = f"{vep_command} {vep_dict['vep-container']['vep-image']} vep"

    #elif vep_mode == "DOCKER":
    #    logger.info("Using containerized VEP with Docker for consequence annotation.")
    #    vep_command = f"docker run -v {annotation_dict['resources-base-path']}:{annotation_dict['resources-base-path']}:ro -v {input_vcf_file}:{input_vcf_file}:ro
    #    vep_command = f"{vep_command} {vep_dict['vep-container']['vep-docker']} vep"

    else:
        sys.exit()

    cache_path = annotation_dict['vep-cache'] + "/"

    # add essential parameters
    vep_command = f"{vep_command} --species homo_sapiens --assembly {assembly_build}"
    vep_command = f"{vep_command} --offline"
    vep_command = f"{vep_command} --cache"
    vep_command = f"{vep_command} --dir_cache {cache_path}"
    vep_command = f"{vep_command} --gencode_basic"

    vep_command = f"{vep_command} -i {input_vcf_file}"
    vep_command = f"{vep_command} -o {output_vcf_file}"
    vep_command = f"{vep_command} --fields Consequence"
    vep_command = f"{vep_command} --vcf_info_field CONS"

    # the developers of VEP recommend to use 4 threads
    if num_cores >= 4:
        vep_command = f"{vep_command} --fork 4 "

    elif num_cores == 2:
        vep_command = f"{vep_command} --fork 2 "

    vep_command = f"{vep_command} --vcf"
    vep_command = f"{vep_command} --no_stats"
    vep_command = f"{vep_command} --force_overwrite"

    print(vep_command)
    subprocess.run(vep_command, shell=True, check=True)
    logger.debug(f"The consequence annotated VCF is saved as {output_vcf_file}")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Annotation with VEP")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="data.vcf", required=True, help="VCF file containing the data, you want to annotate with VEP\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.vcf", required=True, help="Specifies the annotated output file\n")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the annotation parameters\n")
    parser.add_argument("--reference", type=str, dest="reference", metavar="GRCh38.fa", required=True, help="Reference Genome used during variant calling\n")
    parser.add_argument("--basic", dest="basic", action="store_true", required=False, help="Flag to perform basic annotation on InDels\n")
    parser.add_argument("--expanded", dest="expanded", action="store_true", required=False, help="Flag to perform annotation on expanded InDels\n")
    parser.add_argument("--inhouse_sample", dest="inhouse_sample", action="store_true", required=False, help="Flag to indicate that we are annotating an inhouse sample (skips leftNormalize since it is already performed)\n")
    parser.add_argument("--threads", dest="threads", metavar=1, required=False, help="Number of threads to use during annotation\n")
    args = parser.parse_args()

    input_vcf_file = args.in_data
    output_vcf_file = args.out_data

    # parse configuration file
    with open(args.config, "r") as config_file:
        configuration = yaml.load(config_file, Loader=yaml.SafeLoader)

    assembly_build = configuration["Assembly-Build"]
    annotation_dict = configuration["Annotation-Resources"]

    if args.threads is not None:
        num_threads = int(args.threads)
    else:
        num_threads = 1

    basic_annotation = args.basic
    expanded_annotation = args.expanded
    inhouse_sample = args.inhouse_sample

    try:
        # create intermediate temp files
        tmp_sorted = tempfile.NamedTemporaryFile(mode="w+b", suffix="_sorted.vcf", delete=False)
        tmp_vep_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_vepAnnot.vcf", delete=False)
        tmp_vcf_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_vcfAnnot.vcf", delete=False)
        tmp_bed_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_bedAnnot.vcf", delete=False)
        tmp_bigwig_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_bigwigAnnot.vcf", delete=False)

        # perform annotations
        left_normalize_and_sort_vcf(input_vcf_file, tmp_sorted.name, annotation_dict, args.reference, inhouse_sample)
        call_vep_and_annotate_vcf(tmp_sorted.name, tmp_vep_annot.name, annotation_dict, assembly_build, basic_annotation, expanded_annotation, num_threads)
        annotate_from_vcf(tmp_vep_annot.name, tmp_vcf_annot.name, annotation_dict, expanded_annotation, basic_annotation, num_threads)

        if basic_annotation:
            annotate_from_bed(tmp_vcf_annot.name, tmp_bed_annot.name, annotation_dict, num_threads)
            filter_regions(tmp_bed_annot.name, output_vcf_file, annotation_dict)

        elif expanded_annotation:
            annotate_from_bigwig(tmp_vcf_annot.name, output_vcf_file, annotation_dict, num_threads)

        else:
            annotate_from_bed(tmp_vcf_annot.name, tmp_bed_annot.name, annotation_dict, num_threads)
            annotate_from_bigwig(tmp_bed_annot.name, tmp_bigwig_annot.name, annotation_dict, num_threads)
            filter_regions(tmp_bigwig_annot.name, output_vcf_file, annotation_dict)

    finally:
        # clean up
        os.remove(tmp_sorted.name)
        os.remove(tmp_vep_annot.name)
        os.remove(tmp_vcf_annot.name)
        os.remove(tmp_bed_annot.name)
        os.remove(tmp_bigwig_annot.name)

