import argparse
import numpy as np
import pandas as pd
import warnings

from operator import itemgetter


SYNONYMOUS_VARIANTS = ["synonymous_variant",
                       "start_retained_variant",
                       "stop_retained_variant"]

SPLICE_VARIANTS = ["splice_acceptor_variant",
                   "splice_donor_variant",
                   "splice_donor_5th_base_variant",
                   "splice_region_variant",
                   "splice_donor_region_variant",
                   "splice_polypyrimidine_tract_variant"]

VARIANT_CONSEQUENCES = {"transcript_ablation": 1,
                        "splice_acceptor_variant": 2,
                        "splice_donor_variant": 3,
                        "stop_gained": 4,
                        "frameshift_variant": 5,
                        "stop_lost": 6,
                        "start_lost": 7,
                        "transcript_amplification": 8,
                        "inframe_insertion": 9,
                        "inframe_deletion": 10,
                        "missense_variant": 11,
                        "protein_altering_variant": 12,
                        "splice_region_variant": 13,
                        "splice_donor_5th_base_variant": 14,
                        "splice_donor_region_variant": 15,
                        "splice_polypyrimidine_tract_variant": 16,
                        "incomplete_terminal_codon_variant": 17,
                        "start_retained_variant": 18,
                        "stop_retained_variant": 19,
                        "synonymous_variant": 20,
                        "coding_sequence_variant": 21,
                        "mature_miRNA_variant": 22,
                        "5_prime_UTR_variant": 23,
                        "3_prime_UTR_variant": 24,
                        "non_coding_transcript_exon_variant": 25,
                        "intron_variant": 26,
                        "NMD_transcript_variant": 27,
                        "non_coding_transcript_variant": 28,
                        "upstream_gene_variant": 29,
                        "downstream_gene_variant": 30,
                        "TFBS_ablation": 31,
                        "TFBS_amplification": 32,
                        "TF_binding_site_variant": 33,
                        "regulatory_region_ablation": 34,
                        "regulatory_region_amplification": 35,
                        "feature_elongation": 36,
                        "regulatory_region_variant": 37,
                        "feature_truncation": 38,
                        "intergenic_variant": 39,
                        # use "unknown" consequence as default if new consequence terms are added to the database that are not yet implemented (this prevents the program from exiting with an error)
                        "unknown": 40}


def annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, feature):
    with warnings.catch_warnings():
        # we expect to see RuntimeWarnings if the values for a certain variant are missing (the following filter will prevent them from bloating the log file)
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')

        if grouped_expanded_vcf[feature].get_group(row["INDEL_ID"]).empty:
            logger.error(f"Could not combine expanded InDels, INDEL_ID {row['INDEL_ID']} missing in data!")
            return np.nan

        else:
            current_group = grouped_expanded_vcf.get_group(row["INDEL_ID"])

            # we only use non-splicing und non-synonymous variants with impact High or Moderate
            current_group = current_group[(~(current_group["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SYNONYMOUS_VARIANTS))) & ~(current_group["MOST_SEVERE_CONSEQUENCE"].str.contains("|".join(SPLICE_VARIANTS)))) & ((current_group["IMPACT"] == "HIGH") | (current_group["IMPACT"] == "MODERATE"))]

            # to prevent underscoring of High impact variants
            current_group.loc[(current_group["IMPACT"] == "High") & (current_group[feature].isna()), feature] = 1.0

            return current_group[feature].mean()


def combine_vcf_dataframes(grouped_expanded_vcf, vcf_as_dataframe):
    vcf_as_dataframe["phyloP_vertebrate"] = vcf_as_dataframe.apply(lambda row : pd.Series(annotate_indels_with_combined_snps_information(row, grouped_expanded_vcf, "phyloP_vertebrate")), axis=1)

    return vcf_as_dataframe


def extract_columns(cell, process_indel):
    info_fields = str(cell).split(";")

    indel_ID = np.nan
    phyloP_vertebrate = np.nan
    gnomAD_afr_af = np.nan
    gnomAD_amr_af = np.nan
    gnomAD_eas_af = np.nan
    gnomAD_nfe_af = np.nan
    gnomAD_sas_af = np.nan
    max_af = np.nan
    spliceAI = np.nan
    spliceAI_raw = []
    oe_lof = np.nan
    oe_mis = np.nan
    oe_syn = np.nan
    clinvar_details = ""
    hgmd_class = ""
    omim_details = ""
    annotation = ""

    for field in info_fields:
        field_splitted = field.split("=")

        # skip empty INFO field annotations to prevent problems
        if field_splitted[1] != "nan" and field_splitted[1] != "."  and field_splitted[1] != "":

                if field_splitted[0] == "INDEL_ID":
                    indel_ID = field_splitted[1]

                elif field_splitted[0] == "CSQ":
                    annotation = field_splitted[1]

                elif field_splitted[0] == "phyloP_vertebrate":
                    if "&" in field_splitted[1]:
                        phyloP_vertebrate = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)

                    else:
                        phyloP_vertebrate = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_AFR_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_afr_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)

                    else:
                        gnomAD_afr_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_AMR_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_amr_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)

                    else:
                        gnomAD_amr_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_EAS_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_eas_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)

                    else:
                        gnomAD_eas_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_NFE_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_nfe_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)

                    else:
                        gnomAD_nfe_af = float(field_splitted[1])

                elif field_splitted[0] == "gnomAD_SAS_AF":
                    if "&" in field_splitted[1]:
                        gnomAD_sas_af = max([float(value) for value in field_splitted[1].split("&") if (value != "." and value != "nan" and value !="")], default=np.nan)

                    else:
                        gnomAD_sas_af = float(field_splitted[1])

                elif field_splitted[0] == "SpliceAI":
                    spliceAI_raw = field_splitted[1].split("|")
                    spliceAI = max([float(spliceAI_raw[2]), float(spliceAI_raw[3]), float(spliceAI_raw[4]), float(spliceAI_raw[5])], default=np.nan)

                elif field_splitted[0] == "oe_lof":
                    oe_lof = min([float(value) for value in field_splitted[1].split("&")  if (value != "." and value != "nan" and value !="")], default=np.nan)

                ## currently not used
                elif field_splitted[0] == "oe_mis":
                    oe_mis = min([float(value) for value in field_splitted[1].split("&")  if (value != "." and value != "nan" and value !="")], default=np.nan)

                ## currently not used
                elif field_splitted[0] == "oe_syn":
                    oe_syn = min([float(value) for value in field_splitted[1].split("&")  if (value != "." and value != "nan" and value !="")], default=np.nan)

                elif field_splitted[0] == "CLINVAR_DETAILS":
                    clinvar_details = field_splitted[1]

                elif field_splitted[0] == "HGMD_CLASS":
                    hgmd_class = field_splitted[1]

                elif field_splitted[0] == "OMIM":
                    omim_details = field_splitted[1]

                else:
                    print(f"ERROR: Could not recognize INFO field {field}")

        else:
            print(f"DEBUG: Skip empty (NaNs are handled as empty) INFO field {field}")

    max_af = max([gnomAD_afr_af, gnomAD_amr_af, gnomAD_eas_af, gnomAD_nfe_af, gnomAD_sas_af], default=np.nan)

    if process_indel:
        extracted_columns = [indel_ID, annotation, phyloP_vertebrate, gnomAD_afr_af, gnomAD_amr_af, gnomAD_eas_af, gnomAD_nfe_af, gnomAD_sas_af, max_af, oe_lof, oe_mis, oe_syn, spliceAI, clinvar_details, hgmd_class, omim_details]

    else:
        extracted_columns = [annotation, phyloP_vertebrate, gnomAD_afr_af, gnomAD_amr_af, gnomAD_eas_af, gnomAD_nfe_af, gnomAD_sas_af, max_af, oe_lof, oe_mis, oe_syn, spliceAI, clinvar_details, hgmd_class, omim_details]

    return extracted_columns


def extract_vep_annotation(cell, annotation_header, canonical_transcripts=[]):
    annotation_fields = str(cell["CSQ"]).split(",")
    new_cols = []

    if (len(annotation_fields) >= 1) and (annotation_fields[0] != ""):

        if canonical_transcripts:
            for annotation in annotation_fields:
                transcript_index = annotation_header.index("Feature")

                if annotation.split("|")[transcript_index] in canonical_transcripts:
                    new_cols = annotation.split("|")
                    break

        if new_cols == []:
            # choose the most severe annotation variant
            # if new consequence terms were added to the database that are not yet handled from aiDIVA use default consequence "unknown" with lowest severity value
            consequences = [min([VARIANT_CONSEQUENCES.get(consequence) if consequence in VARIANT_CONSEQUENCES.keys() else VARIANT_CONSEQUENCES.get("unknown") for consequence in field.split("|")[annotation_header.index("Consequence")].split("&")]) for field in annotation_fields]
            target_index = min(enumerate(consequences), key=itemgetter(1))[0]
            new_cols = annotation_fields[target_index].strip().split("|")

    else:
        # can happen with the expanded InDels
        print("WARNING: Empty VEP annotation!")

    return new_cols


# Find and return common prefix of two strings
def find_common_prefix(string_a, string_b):
    prefix_end_position = 0

    length_string_a = len(string_a)
    length_string_b = len(string_b)
    max_length_prefix = min(length_string_a, length_string_b)

    while prefix_end_position < max_length_prefix and string_a[prefix_end_position] == string_b[prefix_end_position]:
        prefix_end_position += 1

    shared_prefix = string_a[:prefix_end_position]

    return shared_prefix


# Find and return common suffix of two strings
def find_common_suffix(string_a, string_b):
    suffix_start_position = 0

    length_string_a = len(string_a)
    length_string_b = len(string_b)
    max_length_suffix = min(length_string_a, length_string_b)

    while suffix_start_position < max_length_suffix and string_a[length_string_a-suffix_start_position-1] == string_b[length_string_b-suffix_start_position-1]:
        suffix_start_position += 1

    if suffix_start_position == 0:
        return ""

    common_suffix = string_a[suffix_start_position:]

    return common_suffix


def convert_variant_representation(row):
    chrom = row["#CHROM"]
    start_position = row["POS"]
    ref = row["REF"].upper()
    alt = row["ALT"].upper()

    # remove common first base
    if ref != "" and alt != "" and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        start_position +=1;

    # remove common suffix
    suffix_length = len(find_common_suffix(ref, alt));
    if suffix_length > 0:
        ref = ref[:-suffix_length]
        alt = alt[:-suffix_length]

    # remove common prefix
    prefix_length = len(find_common_prefix(ref, alt))
    if prefix_length > 0:
        ref = ref[prefix_length:]
        alt = alt[prefix_length:]
        start_position += prefix_length

    # determine start and end
    end_position = start_position
    ref_length = len(ref)
    alt_length = len(alt)

    if alt_length == 1 and ref_length == 1: # SNV
        # nothing to do for SNVs
        pass

    elif ref_length == 0: # insertion
        ref = "-"
        # change insertions from before the coordinate to after the coordinate!!!!
        start_position -= 1
        end_position -= 1

    elif alt_length == 0: # deletion
        end_position = start_position + ref_length - 1
        alt="-"

    elif alt_length >= 1 and ref_length > 1: #complex indel
        end_position = start_position + ref_length - 1

    normalized_variant = f"{chrom}:{start_position}-{end_position}_{ref}>{alt}"

    return normalized_variant


def extract_sample_information(row, sample, sample_header=None):
    if str(row["FORMAT"]) == "MULTI":
        if sample_header is not None:
            sample_dict = {}
            sample_information = []
            for sample_entry in row[sample].split(","):
                splitted_entry = sample_entry.split("=")
                sample_dict[splitted_entry[0]] = splitted_entry[1]

            for sample_id in sample_header:
                splitted_sample_information = sample_dict[sample_id].split("|")
                if splitted_sample_information[0] == "wt":
                    sample_information.append("0/0")

                elif splitted_sample_information[0] == "hom":
                    sample_information.append("1/1")

                elif splitted_sample_information[0] == "het":
                    sample_information.append("0/1")

                else:
                    print("WARNING: Genotype not recognized! (%s)" % (splitted_sample_information[0]))

                sample_information.append(splitted_sample_information[1])
                sample_information.append(splitted_sample_information[2])
                sample_information.append(sample_id + "=" + sample_dict[sample_id])

        else:
            print("ERROR: Format is MULTI but no sample_header is given!")

    else:
        format_entries = str(row["FORMAT"]).strip().split(":")
        sample_fields = str(row[sample + "_full"]).strip().split(":")

        if len(format_entries) != len(sample_fields):
            num_missing_entries = abs(len(format_entries) - len(sample_fields))
            for i in range(num_missing_entries):
                sample_fields.append(".")

        if "GT" in format_entries:
            sample_gt_information = sample_fields[format_entries.index("GT")]
        
        else:
            sample_gt_information = "./."

    return sample_gt_information


def add_INFO_fields_to_dataframe(process_indel, expanded_indel, vcf_as_dataframe):
    indel_annotation_columns = ["INDEL_ID", "CSQ", "phyloP_vertebrate", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "MAX_AF", "oe_lof", "oe_mis", "oe_syn", "SpliceAI", "CLINVAR_DETAILS", "HGMD_CLASS", "OMIM"]
    snp_annotation_columns = ["CSQ", "phyloP_vertebrate", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_EAS_AF", "gnomAD_NFE_AF", "gnomAD_SAS_AF", "MAX_AF", "oe_lof", "oe_mis", "oe_syn", "SpliceAI", "CLINVAR_DETAILS", "HGMD_CLASS", "OMIM"]

    if process_indel:
        if not expanded_indel:
            vcf_as_dataframe["GSVAR_VARIANT"] = vcf_as_dataframe.apply(lambda row: convert_variant_representation(row), axis=1)
            vcf_as_dataframe[indel_annotation_columns] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))
            vcf_as_dataframe["IS_INDEL"] = 1

        else:
            vcf_as_dataframe[indel_annotation_columns] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))

    else:
        vcf_as_dataframe["GSVAR_VARIANT"] = vcf_as_dataframe.apply(lambda row: convert_variant_representation(row), axis=1)
        vcf_as_dataframe[snp_annotation_columns] = vcf_as_dataframe["INFO"].apply(lambda x: pd.Series(extract_columns(x, process_indel)))
        vcf_as_dataframe["IS_INDEL"] = 0

    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["INFO"])

    return vcf_as_dataframe


def add_VEP_annotation_to_dataframe(annotation_header, canonical_transcripts, vcf_as_dataframe):
    vcf_as_dataframe[annotation_header] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_vep_annotation(x, annotation_header, canonical_transcripts)), axis=1)
    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["CSQ"])

    vcf_as_dataframe["MOST_SEVERE_CONSEQUENCE"] = vcf_as_dataframe.apply(lambda row: get_most_severe_consequence(row), axis=1)

    return vcf_as_dataframe


def add_sample_information_to_dataframe(sample_ids, sample_header, vcf_as_dataframe):
    for sample in sample_ids:
        if (sample == "trio") or (sample == "multi"):
            sample_header_multi = []
            for sample_id in sample_header:
                sample_header_multi.append("GT_" + sample_id)

            vcf_as_dataframe[sample_header_multi] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_sample_information(x, sample, sample_header)), axis=1)

        else:
            vcf_as_dataframe.rename(columns={sample: sample + "_full"}, inplace=True)
            sample_header_default = ["GT_" + sample]
            vcf_as_dataframe[sample_header_default] = vcf_as_dataframe.apply(lambda x: pd.Series(extract_sample_information(x, sample)), axis=1)
            vcf_as_dataframe = vcf_as_dataframe.drop(columns=[sample + "_full"])

    vcf_as_dataframe = vcf_as_dataframe.drop(columns=["FORMAT"])

    return vcf_as_dataframe


def convert_vcf_to_pandas_dataframe(input_file, process_indel, expanded_indel):
    header, vcf_as_dataframe = reformat_vcf_file_and_read_into_pandas_and_extract_header(input_file)

    #print("DEBUG: Convert VCF file")

    canonical_transcripts = []

    sample_ids = []
    # FORMAT column has index 8 (counted from 0) and sample columns follow afterwards (sample names are unique)
    # Check if FORMAT column exists
    if len(vcf_as_dataframe.columns) > 8:
        for i in range(9, len(vcf_as_dataframe.columns)):
            sample_ids.append(vcf_as_dataframe.columns[i])

        sample_header = extract_sample_header(header)

    annotation_header = extract_annotation_header(header)

    if not vcf_as_dataframe.empty:
        vcf_as_dataframe = add_INFO_fields_to_dataframe(process_indel, expanded_indel, vcf_as_dataframe)
        vcf_as_dataframe = add_VEP_annotation_to_dataframe(annotation_header, canonical_transcripts, vcf_as_dataframe)

        if len(vcf_as_dataframe.columns) > 8:
            if "FORMAT" in vcf_as_dataframe.columns:
                vcf_as_dataframe = add_sample_information_to_dataframe(sample_ids, sample_header, vcf_as_dataframe)

            else:
                # This warning is always triggered when the expanded indel vcf file is processed. If it is only triggered once in this case it can be ignored.
                print("WARNING: It seems that your VCF file does contain sample information but no FORMAT description!")

        else:
            print("ERROR: MISSING SAMPLE INFORMATION!")

        # replace empty strings or only spaces with NaN
        vcf_as_dataframe = vcf_as_dataframe.replace(r"^\s*$", np.nan, regex=True)

    else:
        print("ERROR: The given VCF file is empty!")

    return vcf_as_dataframe


def reformat_vcf_file_and_read_into_pandas_and_extract_header(filepath):
    header_line = ""
    comment_lines = []

    if filepath.endswith(".gz"):
        vcf_file_to_reformat = gzip.open(filepath, "rt")

    else:
        vcf_file_to_reformat = open(filepath, "r")

    # extract header from vcf file
    for line in vcf_file_to_reformat:
        if line.strip().startswith("##"):
            comment_lines.append(line.strip())

        elif line.strip().startswith("#CHROM"):
            header_line = line.strip()
            comment_lines.append(header_line)
            break # now the variant entries are coming

        else:
            continue

    if header_line == "":
        print("WARNING: The VCF file seems to be corrupted")

    vcf_file_to_reformat.close()

    vcf_header = header_line.strip().split("\t")
    vcf_as_dataframe = pd.read_csv(filepath, names=vcf_header, sep="\t", comment="#", low_memory=False)

    return comment_lines, vcf_as_dataframe


def extract_annotation_header(header):
    annotation_header = [entry.strip().replace("\">", "").split(": ")[1].split("|") for entry in header if entry.startswith("##INFO=<ID=CSQ,")][0]

    return annotation_header


def extract_sample_header(header):
    sample_header = []

    for line in header:
        if line.startswith("##SAMPLE=<"):
            sample_entry = line.strip().replace("##SAMPLE=<", "").replace(">", "")
            sample_id = sample_entry.split(",")[0]
            sample_header.append(sample_id.split("=")[1])

    return sample_header


def get_sample_information(row):
    if str(row) == "1/0":
        sample_information = "het"

    elif str(row) == "0/1":
        sample_information = "het"

    elif str(row) == "0/0":
        sample_information = "wt"

    elif str(row) == "1/1":
        sample_information = "hom"

    else:
        sample_information = "n/a"

    return sample_information


def get_gnomad_allele_frequency(row):
    gnomAD = row["MAX_AF"]

    gnomad_afr = str(row[f"gnomAD_AFR_AF"])
    gnomad_amr = str(row[f"gnomAD_AMR_AF"])
    gnomad_eas = str(row[f"gnomAD_EAS_AF"])
    gnomad_nfe = str(row[f"gnomAD_NFE_AF"])
    gnomad_sas = str(row[f"gnomAD_SAS_AF"])

    if gnomad_afr == "nan" and gnomad_amr == "nan" and gnomad_eas == "nan" and gnomad_nfe == "nan" and gnomad_sas == "nan":
        gnomAD_sub = np.nan

    else:
        gnomAD_sub = ",".join([gnomad_afr, gnomad_amr, gnomad_eas, gnomad_nfe, gnomad_sas])

    return [gnomAD, gnomAD_sub]


def extract_gene_information(row):
    gene = row["SYMBOL"]
    inheritance = "n/a"

    if str(row["oe_lof"]) == "nan":
        oe_lof = "n/a"

    else:
        oe_lof = row["oe_lof"]

    #if str(row["oe_mis"]) == "nan":
    #    oe_mis = "n/a"

    #else:
    #    oe_mis = row["oe_mis"]

    #if str(row["oe_syn"]) == "nan":
    #    oe_syn = "n/a"

    #else:
    #    oe_syn = row["oe_syn"]

    oe_mis = "n/a"
    oe_syn = "n/a"

    gene_information = f"{gene} (inh={inheritance} oe_syn={oe_syn} oe_mis={oe_mis} oe_lof={oe_lof})"

    return gene_information


def extract_variant_information(row):
    gene = str(row["SYMBOL"])

    if str(row["Feature"]) != "nan":
        transcript_id = str(row["Feature"])

    else:
        exon_intron_number = ""

    if str(row["Consequence"]) != "nan":
        consequence = str(row["Consequence"])

    else:
        exon_intron_number = ""

    if str(row["IMPACT"]) != "nan":
        impact = str(row["IMPACT"])

    else:
        exon_intron_number = ""

    if str(row["EXON"]) != "nan":
        exon_intron_number = str(row["EXON"])

    elif str(row["INTRON"]) != "nan":
        exon_intron_number = str(row["INTRON"])

    else:
        exon_intron_number = ""

    if str(row["HGVSc"]) != "nan":
        hgvsc = str(row["HGVSc"])

    else:
        hgvsc = ""

    if str(row["HGVSp"]) != "nan":
        hgvsp = str(row["HGVSp"])

    else:
        hgvsp = ""

    pfam_domain = ""

    variant_information = ":".join([gene, transcript_id, consequence, impact, exon_intron_number, hgvsc, hgvsp, pfam_domain])

    return variant_information


def get_clinvar_information(row):
    clinvar_information = row["CLINVAR_DETAILS"]
    if clinvar_information == "Conflicting_classifications_of_pathogenicity":
        clinvar_details = np.nan

    elif "benign" in clinvar_information and "pathogenic" in clinvar_information:
        clinvar_details = np.nan

    elif "benign" in clinvar_information and "uncertain" in clinvar_information:
        clinvar_details = np.nan

    elif "pathogenic" in clinvar_information and "uncertain" in clinvar_information:
        clinvar_details = np.nan

    else:
        clinvar_details = clinvar_information

    return clinvar_details


def extract_variant_position(row):
    gsvar_variant = row["GSVAR_VARIANT"]
    chrom = gsvar_variant.split(":")[0]
    position = gsvar_variant.split(":")[1].split("_")[0]
    variation = gsvar_variant.split("_")[1]
    start = position.split("-")[0]
    end = position.split("-")[1]
    ref = variation.split(">")[0]
    obs = variation.split(">")[1]

    return [chrom, start, end, ref, obs]


# extract the most severe consequence if overlapping consequences were found
def get_most_severe_consequence(row):
    consequences = str(row["Consequence"])
    found_consequences = consequences.split("&")

    # use most severe consequence for filtering if overlapping consequences are present
    if len(found_consequences) > 1:
        most_severe_consequence = found_consequences[0]
        
        for consequence in found_consequences:
            if VARIANT_CONSEQUENCES[consequence] < VARIANT_CONSEQUENCES[most_severe_consequence]:
                most_severe_consequence = consequence
    else:
        most_severe_consequence = found_consequences[0]
    
    return most_severe_consequence


def main(snp_vcf, indel_vcf, indel_expanded_vcf, out_file, sample_id):
    snp_data = convert_vcf_to_pandas_dataframe(snp_vcf, False, False)
    indel_data = convert_vcf_to_pandas_dataframe(indel_vcf, True, False)
    indel_expanded_data = convert_vcf_to_pandas_dataframe(indel_expanded_vcf, True, True)

    indel_expanded_data["phyloP_vertebrate"] = indel_expanded_data["phyloP_vertebrate"].apply(lambda row: max([float(value) for value in str(row).split("&") if ((value != ".") & (value != "nan"))], default=np.nan))
    grouped_expanded_data = indel_expanded_data.groupby("INDEL_ID")
    combined_indel_data = combine_vcf_dataframes(grouped_expanded_data, indel_data)

    input_table = pd.concat([snp_data, combined_indel_data])
    input_table.sort_values(["#CHROM", "POS"], ascending=[True, True], inplace=True)
    input_table.reset_index(inplace=True, drop=True)
    input_table = input_table[snp_data.columns]

    with open(out_file, "w") as outfile:
        #input_table = pd.read_csv(in_file, sep="\t", low_memory=False)
        input_table[["#chr", "start", "end", "ref", "obs"]] = input_table.apply(lambda x: pd.Series(extract_variant_position(x)), axis=1)

        if f"GT_{sample_id}" in input_table.columns:
            input_table[sample_id] = input_table[f"GT_{sample_id}"].apply(lambda x: get_sample_information(x))

        else:
            sample_id = "sample_info"
            input_table["sample_info"] = "n/a"

        outfile.write(f"""##GENOME_BUILD=GRCh38
##CREATION_DATE=YEAR-MONTH-DAY
##CALLING_DATE=YEAR-MONTH-DAY
##SAMPLE=<ID={sample_id},DiseaseStatus=affected>
##FILTER=off-target=Variant marked as 'off-target'.
##FILTER=low_conf_region=Variant marked as 'low_conf_region'.
##DESCRIPTION={sample_id}=Genotype of variant in sample.
##DESCRIPTION=filter=Annotations for filtering and ranking variants.
##DESCRIPTION=quality=Quality parameters - variant quality (QUAL)
##DESCRIPTION=coding_and_splicing=Coding and splicing details (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain).
##DESCRIPTION=OMIM=OMIM database annotation.
##DESCRIPTION=ClinVar=ClinVar database annotation.
##DESCRIPTION=HGMD=HGMD database annotation.
##DESCRIPTION=gnomAD=Maximum allele frequency from the main sub populations (AFR,AMR,EAS,NFE,SAS).
##DESCRIPTION=gnomAD_sub=Sub-population allele frequenciens (AFR,AMR,EAS,NFE,SAS) in gnomAD project.
##DESCRIPTION=phyloP=phyloP (100way vertebrate) annotation..
##DESCRIPTION=MaxEntScan=MaxEntScan reference score and alternate score for (1) native splice site, (2) acceptor gain and (3) donor gain. Comma-separated list if there are different predictions for several transcripts.
##DESCRIPTION=SpliceAI=SpliceAI prediction of splice-site variations. Probability of the variant being splice-altering (range from 0-1). The score is the maximum value of acceptor/donor gain/loss of all effected genes.
##DESCRIPTION=NGSD_hom=Homozygous variant count in NGSD.
##DESCRIPTION=NGSD_het=Heterozygous variant count in NGSD.
##DESCRIPTION=classification=Classification from the NGSD.
##DESCRIPTION=gene_info=Gene information from NGSD (inheritance mode, gnomAD o/e scores).
""")

        input_table["quality"] = input_table["QUAL"].apply(lambda x: "QUAL=" + str(x))
        input_table["filter"] = input_table["FILTER"].apply(lambda x: np.nan if x == "." else x)

        input_table["coding_and_splicing"] = input_table.apply(lambda x: extract_variant_information(x), axis=1)

        input_table[["gnomAD", "gnomAD_sub"]] = input_table.apply(lambda x: pd.Series(get_gnomad_allele_frequency(x)), axis=1)

        if "CLINVAR_DETAILS" in input_table.columns:
            input_table = input_table.rename(columns={"CLINVAR_DETAILS": "ClinVar"})

        else:
            input_table["ClinVar"] = np.nan

        if "HGMD" in input_table.columns:
            input_table["HGMD"] = input_table["HGMD_CLASS"]

        else:
            input_table = input_table.rename(columns={"HGMD_CLASS": "HGMD"})

        input_table = input_table.rename(columns={"phyloP_vertebrate": "phyloP"})

        input_table["oe_lof"] = input_table.apply(lambda row: min([float(value) for value in str(row["oe_lof"]).split("&") if ((value != ".") and (value != "nan") and (value != "") and (":" not in value) and ("-" not in value))], default=np.nan), axis=1)
        input_table["gene_info"] = input_table.apply(lambda x: extract_gene_information(x), axis=1)
        input_table["NGSD_hom"] = np.nan
        input_table["NGSD_het"] = np.nan
        input_table["classification"] = np.nan
        input_table["MaxEntScan"] = np.nan

        input_table[["#chr", "start", "end", "ref", "obs", "filter", "coding_and_splicing", "OMIM", "ClinVar", "HGMD", "gnomAD", "gnomAD_sub", "phyloP", "SpliceAI", "gene_info", "NGSD_hom", "NGSD_het", "classification","MaxEntScan", "quality", sample_id]].to_csv(outfile, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp_vcf", type=str, dest="snp_vcf", metavar="snp_annot.vcf", required=True, help="Annotated SNP VCF.")
    parser.add_argument("--indel_vcf", type=str, dest="indel_vcf", metavar="indel_annot.vcf", required=True, help="Annotated InDel VCF.")
    parser.add_argument("--indel_expanded_vcf", type=str, dest="indel_expanded_vcf", metavar="indel_expanded_annot.vcf", required=True, help="Annotated expanded InDel VCF.")
    parser.add_argument("--out_file", type=str, dest="out_file", metavar="output.csv", required=True, help="Converted basic GSvar file.")
    parser.add_argument("--sample_id", type=str, dest="sample_id", metavar="NA12878", required=False, help="Sample ID to get the sample information from the table.")
    args = parser.parse_args()

    if args.sample_id is not None:
        sample_id = args.sample_id

    else:
        sample_id = "UNKNOWN"

    main(args.snp_vcf, args.indel_vcf, args.indel_expanded_vcf, args.out_file, sample_id)
