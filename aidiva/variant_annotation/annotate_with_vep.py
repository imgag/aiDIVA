import subprocess
import argparse
import tempfile
import os
import logging
import yaml


logger = logging.getLogger(__name__)


# TODO: restrict annotation to SIFT and PolyPhen (get rid of VEP in the long run)
def call_vep_and_annotate_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict, build="GRCh37", basic=False, expanded=False, num_cores=1):
    # the path to the executable
    vep_command = vep_annotation_dict["vep"] + "/" + "vep "

    # set the correct paths to the needed perl modules
    if "PERL5LIB" in os.environ:
        #os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/:" + os.environ["PERL5LIB"]
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:/mnt/storage2/share/opt/perl_cpan_ubuntu20.04/lib/perl5/:" + os.environ["PERL5LIB"]
    else:
        #os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/"
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:/mnt/storage2/share/opt/perl_cpan_ubuntu20.04/lib/perl5/"

    cache_path = vep_annotation_dict['vep-cache'] + "/"

    # add essential parameters
    vep_command = f"{vep_command} --species homo_sapiens --assembly {build} "
    vep_command = f"{vep_command} --offline"
    vep_command = f"{vep_command} --cache"
    vep_command = f"{vep_command} --dir_cache {cache_path}"
    vep_command = f"{vep_command} --gencode_basic"
    vep_command = f"{vep_command} --symbol"
    vep_command = f"{vep_command} --biotype"
    vep_command = f"{vep_command} --variant_class"

    if not expanded:
        pass # deprecated
        # allele frequencies to include
        #vep_command = f"{vep_command} --max_af"

        # the following AF annotations could be dropped since we only need the max AF
        #vep_command = f"{vep_command} --af"
        #vep_command = f"{vep_command} --af_1kg"
        #vep_command = f"{vep_command} --af_esp"
        if build == "GRCh37":
            vep_command = f"{vep_command} --af_gnomad"
        elif build == "GRCh38":
            vep_command = f"{vep_command} --af_gnomadg"

    # vep plugins to use
    if not basic:
        vep_command = f"{vep_command} --sift s"
        vep_command = f"{vep_command} --polyphen s"

    vep_command = f"{vep_command} -i " + input_vcf_file + " "
    vep_command = f"{vep_command} -o " + output_vcf_file + " "
    vep_command = f"{vep_command} --fork " + str(num_cores) + " "
    vep_command = f"{vep_command} --format vcf" + " " # we need this to prevent vep from not working if the VCF file has no variant entries
    vep_command = f"{vep_command} --vcf" + " "
    vep_command = f"{vep_command} --no_stats" + " "
    vep_command = f"{vep_command} --force_overwrite"

    subprocess.run(vep_command, shell=True, check=True)
    logger.debug("The VEP annotated VCF is saved as %s" % (output_vcf_file))


# TODO rewrite method to use consequence annotation from ngs-bits instead of VEP
def annotate_consequence_information(input_vcf_file, output_vcf_file, vep_annotation_dict, build="GRCh37", num_cores=1):
    # the path to the executable
    vep_command = vep_annotation_dict["vep"] + "/" + "vep "

    # set the correct paths to the needed perl modules
    if "PERL5LIB" in os.environ:
        #os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/:" + os.environ["PERL5LIB"]
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:/mnt/storage2/share/opt/perl_cpan_ubuntu20.04/lib/perl5/:" + os.environ["PERL5LIB"]
    else:
        #os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:" + vep_annotation_dict["vep"] + "/" + "cpan/lib/perl5/"
        os.environ["PERL5LIB"] = vep_annotation_dict["vep"] + "/" + "Bio/:/mnt/storage2/share/opt/perl_cpan_ubuntu20.04/lib/perl5/"

    cache_path = vep_annotation_dict['vep-cache'] + "/"

    # add essential parameters
    vep_command = f"{vep_command} --species homo_sapiens --assembly {build}"
    vep_command = f"{vep_command} --offline"
    vep_command = f"{vep_command} --cache"
    vep_command = f"{vep_command} --dir_cache {cache_path}"
    vep_command = f"{vep_command} --gencode_basic"

    vep_command = f"{vep_command} -i {input_vcf_file}"
    vep_command = f"{vep_command} -o {output_vcf_file}"
    vep_command = f"{vep_command} --fields Consequence"
    vep_command = f"{vep_command} --vcf_info_field CONS"
    vep_command = f"{vep_command} --fork {num_cores}"
    vep_command = f"{vep_command} --vcf"
    vep_command = f"{vep_command} --no_stats"
    vep_command = f"{vep_command} --force_overwrite"

    subprocess.run(vep_command, shell=True, check=True)
    logger.debug(f"The consequence annotated VCF is saved as {output_vcf_file}")


def annotate_from_vcf(input_vcf_file, output_vcf_file, annotation_dict, expanded=False, basic=False, num_cores=1):
    tmp = tempfile.NamedTemporaryFile(mode="w+b", suffix=".config", delete=False)

    vcf_annotation = annotation_dict['vcf-files']
    command = f"{annotation_dict['ngs-bits']}/VcfAnnotateFromVcf"

    try:
        if not basic:
            tmp.write(f"{vcf_annotation['CONDEL']}\t\tCONDEL\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['EIGEN_PHRED']}\t\tEIGEN_PHRED\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['FATHMM_XF']}\t\tFATHMM_XF\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['MutationAssessor']}\t\tMutationAssessor\t\ttrue\n".encode())
            # add allele frequencies to gnomAD annotation (instead of annotation from VEP)
            tmp.write(f"{vcf_annotation['gnomAD']}\tgnomAD\tAN,Hom\t\ttrue\n".encode())
            #tmp.write(f"{vcf_annotation['gnomAD']}\tgnomAD\tAFR_AF,AMR_AF,EAS_AF,NFE_AF,SAS_AF\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['CAPICE']}\t\tCAPICE\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['dbscSNV']}\t\tADA_SCORE,RF_SCORE\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['CADD']}\t\tCADD\t\ttrue\n".encode())
            tmp.write(f"{vcf_annotation['REVEL']}\t\tREVEL\t\ttrue\n".encode())

        if not expanded:
            tmp.write(f"{vcf_annotation['clinvar']}\tCLINVAR\tDETAILS\t\ttrue\n".encode())

            # HGMD needs a valid license, therefor we check if the file exists otherwise this annotation is skipped
            if os.path.isfile(f"{vcf_annotation['hgmd']}"):
                tmp.write(f"{vcf_annotation['hgmd']}\tHGMD\tCLASS,RANKSCORE\t\ttrue\n".encode())
            else:
                logger.warn("HGMD file is not found! Skip HGMD annotation!")

        # close temporary file to make it accessible
        tmp.close()

        command = f"{command} -config_file {tmp.name} -in {input_vcf_file} -out {output_vcf_file} -threads {num_cores}"
        subprocess.run(command, shell=True, check=True)

    finally:
        # clean up
        os.remove(tmp.name)


def annotate_from_bed(input_vcf_file, output_vcf_file, annotation_dict, num_cores=1):
    bed_annotation = annotation_dict['bed-files']
    command = f"{annotation_dict['ngs-bits']}/VcfAnnotateFromBed"

    try:
        tmp_segDup = tempfile.NamedTemporaryFile(mode="w+b", suffix="_segDup.vcf", delete=False)
        tmp_simpleRepeat = tempfile.NamedTemporaryFile(mode="w+b", suffix="_simpleRepeat.vcf", delete=False)
        tmp_oe_lof = tempfile.NamedTemporaryFile(mode="w+b", suffix="_oe_lof.vcf", delete=False)
        tmp_repeatmasker = tempfile.NamedTemporaryFile(mode="w+b", suffix="_repeatmasker.vcf", delete=False)

        # close temporary files to make them accessible
        tmp_segDup.close()
        tmp_simpleRepeat.close()
        tmp_oe_lof.close()
        tmp_repeatmasker.close()

        subprocess.run(f"{command} -bed {bed_annotation['segmentDuplication']} -name SegDup -sep '&' -in {input_vcf_file} -out {tmp_segDup.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bed {bed_annotation['simpleRepeat']} -name SimpleRepeats -sep '&' -in {tmp_segDup.name} -out {tmp_simpleRepeat.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bed {bed_annotation['oe_lof']} -name oe_lof -sep '&' -in {tmp_simpleRepeat.name} -out {tmp_oe_lof.name} -threads {num_cores}", shell=True, check=True)

        if os.path.isfile(f"{bed_annotation['omim']}"):
            subprocess.run(f"{command} -bed {bed_annotation['repeatMasker']} -name REPEATMASKER -sep '&' -in {tmp_oe_lof.name} -out {tmp_repeatmasker.name} -threads {num_cores}", shell=True, check=True)
            subprocess.run(f"{command} -bed {bed_annotation['omim']} -name OMIM -sep '&' -in {tmp_repeatmasker.name} -out {output_vcf_file} -threads {num_cores}", shell=True, check=True)
        else:
            subprocess.run(f"{command} -bed {bed_annotation['repeatMasker']} -name REPEATMASKER -sep '&' -in {tmp_oe_lof.name} -out {output_vcf_file} -threads {num_cores}", shell=True, check=True)
            logger.warn("OMIM file is not found! Skip OMIM annotation!")

    finally:
        # clean up
        os.remove(tmp_segDup.name)
        os.remove(tmp_simpleRepeat.name)
        os.remove(tmp_oe_lof.name)
        os.remove(tmp_repeatmasker.name)


def annotate_from_bigwig(input_vcf_file, output_vcf_file, annotation_dict, num_cores=1):
    bigwig_annotation = annotation_dict['bigwig-files']
    command = f"{annotation_dict['ngs-bits']}/VcfAnnotateFromBigWig"

    try:
        tmp_phyloP_primate = tempfile.NamedTemporaryFile(mode="w+b", suffix="_phyloP_primate.vcf", delete=False)
        tmp_phyloP_mammal = tempfile.NamedTemporaryFile(mode="w+b", suffix="_phyloP_mammal.vcf", delete=False)
        tmp_phyloP_vertebrate = tempfile.NamedTemporaryFile(mode="w+b", suffix="_phyloP_vertebrate.vcf", delete=False)

        tmp_phastCons_primate = tempfile.NamedTemporaryFile(mode="w+b", suffix="_phastCons_primate.vcf", delete=False)
        tmp_phastCons_mammal = tempfile.NamedTemporaryFile(mode="w+b", suffix="_phastCons_mammal.vcf", delete=False)
        tmp_phastCons_vertebrate = tempfile.NamedTemporaryFile(mode="w+b", suffix="_phastCons_vertebrate.vcf", delete=False)

        # close temporary files to make them accessible
        tmp_phyloP_primate.close()
        tmp_phyloP_mammal.close()
        tmp_phyloP_vertebrate.close()
        tmp_phastCons_primate.close()
        tmp_phastCons_mammal.close()
        tmp_phastCons_vertebrate.close()

        subprocess.run(f"{command} -bw {bigwig_annotation['phyloP_primate']} -name phyloP_primate -desc 'phyloP primate dataset' -mode avg -in {input_vcf_file} -out {tmp_phyloP_primate.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bw {bigwig_annotation['phyloP_mammal']} -name phyloP_mammal -desc 'phyloP mammalian dataset' -mode avg -in {tmp_phyloP_primate.name} -out {tmp_phyloP_mammal.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bw {bigwig_annotation['phyloP_vertebrate']} -name phyloP_vertebrate -desc 'phyloP vertebrate dataset' -mode avg -in {tmp_phyloP_mammal.name} -out {tmp_phyloP_vertebrate.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bw {bigwig_annotation['phastCons_primate']} -name phastCons_primate -desc 'phastCons primate dataset' -mode avg -in {tmp_phyloP_vertebrate.name} -out {tmp_phastCons_primate.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bw {bigwig_annotation['phastCons_mammal']} -name phastCons_mammal -desc 'phastCons mammalian dataset' -mode avg -in {tmp_phastCons_primate.name} -out {tmp_phastCons_mammal.name} -threads {num_cores}", shell=True, check=True)
        subprocess.run(f"{command} -bw {bigwig_annotation['phastCons_vertebrate']} -name phastCons_vertebrate -desc 'phastCons vertebrate dataset' -mode avg -in {tmp_phastCons_mammal.name} -out {output_vcf_file} -threads {num_cores}", shell=True, check=True)

    finally:
        # clean up
        os.remove(tmp_phyloP_primate.name)
        os.remove(tmp_phyloP_mammal.name)
        os.remove(tmp_phyloP_vertebrate.name)
        os.remove(tmp_phastCons_primate.name)
        os.remove(tmp_phastCons_mammal.name)
        os.remove(tmp_phastCons_vertebrate.name)


def filter_regions(input_vcf_file, output_vcf_file, annotation_dict):
    command = f"{annotation_dict['ngs-bits']}/VariantFilterRegions"

    low_confidence_filter = annotation_dict['low-confidence']
    command = f"{command} -in {input_vcf_file} -mark low_conf_region -inv -reg {low_confidence_filter} -out {output_vcf_file}"
    subprocess.run(command, shell=True, check=True)


def sort_vcf(input_vcf_file, output_vcf_file, vep_annotation_dict):
    command = f"{vep_annotation_dict['ngs-bits']}/VcfSort"

    command = f"{command} -in {input_vcf_file} -out {output_vcf_file}"
    subprocess.run(command, shell=True, check=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "AIdiva -- Annotation with VEP")
    parser.add_argument("--in_data", type=str, dest="in_data", metavar="data.vcf", required=True, help="VCF file containing the data, you want to annotate with VEP\n")
    parser.add_argument("--out_data", type=str, dest="out_data", metavar="out.vcf", required=True, help="Specifies the extended output file\n")
    parser.add_argument("--config", type=str, dest="config", metavar="config.yaml", required=True, help="Config file specifying the annotation parameters")
    parser.add_argument("--basic", dest="basic", action="store_true", required=False, help="Flag to perform basic annotation on InDels")
    parser.add_argument("--expanded", dest="expanded", action="store_true", required=False, help="Flag to perform annotation on expanded InDels")
    parser.add_argument("--threads", dest="threads", metavar=1, required=False, help="Number of threads to use during annotation")
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

    try:
        # create intermediate temp files
        tmp_sorted = tempfile.NamedTemporaryFile(mode="w+b", suffix="_sorted.vcf", delete=False)
        tmp_vep_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_vepAnnot.vcf", delete=False)
        tmp_vcf_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_vcfAnnot.vcf", delete=False)
        tmp_bed_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_bedAnnot.vcf", delete=False)
        tmp_bigwig_annot = tempfile.NamedTemporaryFile(mode="w+b", suffix="_bigwigAnnot.vcf", delete=False)

        # perform annotations
        sort_vcf(input_vcf_file, tmp_sorted.name, annotation_dict)
        call_vep_and_annotate_vcf(tmp_sorted.name, tmp_vep_annot.name, annotation_dict, assembly_build, basic_annotation, expanded_annotation, num_threads)
        annotate_from_vcf(tmp_vep_annot.name, tmp_vcf_annot.name, annotation_dict, expanded_annotation, basic_annotation, num_threads)

        if basic_annotation:
            annotate_from_bed(tmp_vcf_annot.name, tmp_bed_annot.name, annotation_dict, num_threads)
            filter_regions(tmp_bed_annot.name, output_vcf_file, annotation_dict)

        elif expanded_annotation:
            annotate_from_bigwig(tmp_vcf_annot.name, output_vcf_file, annotation_dict, num_threads)

        else:
            annotate_from_bed(tmp_vcf_annot.name, tmp_bed_annot.name, annotation_dict, num_threads)
            annotate_from_bigwig(tmp_vcf_annot.name, tmp_bigwig_annot.name, annotation_dict, num_threads)
            filter_regions(tmp_bigwig_annot.name, output_vcf_file, annotation_dict)

    finally:
        # clean up
        os.remove(tmp_sorted.name)
        os.remove(tmp_vep_annot.name)
        os.remove(tmp_vcf_annot.name)
        os.remove(tmp_bed_annot.name)
        os.remove(tmp_bigwig_annot.name)

