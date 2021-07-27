import argparse
import gzip
import tempfile
import logging


coding_variants = ["splice_acceptor_variant",
                   "splice_donor_variant",
                   "stop_gained",
                   "frameshift_variant",
                   "stop_lost",
                   "start_lost",
                   "inframe_insertion",
                   "inframe_deletion",
                   "missense_variant",
                   "protein_altering_variant",
                   "splice_region_variant",
                   "incomplete_terminal_codon_variant",
                   "start_retained_variant",
                   "stop_retained_variant",
                   "synonymous_variant",
                   "coding_sequence_variant",
                   "5_prime_UTR_variant",
                   "3_prime_UTR_variant"]

logger = logging.getLogger(__name__)


def filter_coding_variants(filepath, filepath_out, annotation_field_name):
    if filepath.endswith(".gz"):
        vcf_file_to_reformat = gzip.open(filepath, "rt")
    else:
        vcf_file_to_reformat = open(filepath, "r")
    outfile = open(filepath_out, "w")
    consequence_index = 0

    # make sure that there are no unwanted linebreaks in the variant entries
    tmp = tempfile.NamedTemporaryFile(mode="w+")
    tmp.write(vcf_file_to_reformat.read().replace(r"(\n(?!((((((chr)?[0-9]{1,2}|(chr)?[xXyY]{1}|(chr)?(MT|mt){1})\t)(.+\t){6,}(.+(\n|\Z))))|(#{1,2}.*(\n|\Z))|(\Z))))", ""))
    tmp.seek(0)

    # extract header from vcf file
    indel_ID = 0
    for line in tmp:
        splitted_line = line.split("\t")
        if line.strip().startswith("##"):
            if line.strip().startswith("##fileformat="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##filedate="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##source="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##reference="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##contig="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##FORMAT="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##FILTER="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##ANALYSISTYPE="):
                outfile.write(line)
                continue
            elif line.strip().startswith("##SAMPLE="):
                outfile.write(line)
                continue

            if line.strip().startswith("##INFO=<ID=" + annotation_field_name):
                annotation_header = line.strip().replace("\">", "").split(": ")[1].split("|")
                for i in range(len(annotation_header)):
                    if annotation_header[i] == "Consequence":
                        consequence_index = i
                        break
                continue
            continue

        if line.strip().startswith("#CHROM"):
            outfile.write(line)
            continue

        # skip empty lines if the VCF file is not correctly formatted (eg. if there are multiple blank lines in the end of the file)
        if line == "\n":
            continue

        # check if variant is coding and write to outfile
        annotation_field = ""
        for field in splitted_line[7].split(";"):
            if field.startswith(annotation_field_name + "="):
                annotation_field = field.replace(annotation_field_name + "=", "")

        if annotation_field:
            # check all annotated transcripts
            transcript_annotations = [annotation for annotation in annotation_field.split(",")]
            consequence_list = [transcript.split("|")[consequence_index] for transcript in transcript_annotations]
            consequences = []

            for consequence in consequence_list:
                consequences += consequence.split("&")

            if any(term for term in coding_variants if term in consequences):
                splitted_line[7] = "."
                outfile.write("\t".join(splitted_line))
        else:
            logger.error("Annotation field missing!")
            logger.warn("Variant filtering will be skipped!")

    vcf_file_to_reformat.close()
    outfile.close()
    tmp.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_file", type=str, dest="in_file", metavar="input.vcf", required=True, help="VCF file to filter\n")
    parser.add_argument("--out_file", type=str, dest="out_file", metavar="output.vcf", required=True, help="VCF file containing only the filtered coding variants\n")
    parser.add_argument("--annotation_field", type=str, dest="annotation_field", metavar="CSQ", required=True, help="ID of the annotation field with the Consequence information\n")
    args = parser.parse_args()

    filter_coding_variants(args.in_file, args.out_file, args.annotation_field)
