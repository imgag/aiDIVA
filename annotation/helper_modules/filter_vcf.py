import argparse
import gzip
import logging
import tempfile


logger = logging.getLogger(__name__)


def filter_coding_variants(filepath, filepath_out, annotation_field_name, CONSTANT_DICTIONARY):
    # get CONSTANTS
    CODING_VARIANTS = CONSTANT_DICTIONARY["CODING_VARIANTS"]
    SPLICE_VARIANTS = CONSTANT_DICTIONARY["SPLICE_VARIANTS"]
    SYNONYMOUS_VARIANTS = CONSTANT_DICTIONARY["SYNONYMOUS_VARIANTS"]
    
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
    for line in tmp:
        splitted_line = line.replace("\n", "").split("\t")
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

        # extract annotation header
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

            # check if variant is coding and write to outfile
            ## TODO: check if the for loop can be dropped maybe use any() in the if condition instead
            for term in [*CODING_VARIANTS, *SPLICE_VARIANTS, *SYNONYMOUS_VARIANTS]:
                if term in consequences:
                    splitted_line[7] = "CODING_FILTERED=TRUE"
                    outfile.write(str("\t".join(splitted_line) + "\n"))
                    break

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
