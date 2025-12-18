import argparse
import numpy as np
import pandas as pd


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


def get_gnomad_allele_frequency(row, legacy_mode):
    gnomAD = row["MAX_AF"]

    if legacy_mode:
        genome_identifier = ""

    else:
        genome_identifier = "g"

    gnomad_afr = str(row[f"gnomAD{genome_identifier}_AFR_AF"])
    gnomad_amr = str(row[f"gnomAD{genome_identifier}_AMR_AF"])
    gnomad_eas = str(row[f"gnomAD{genome_identifier}_EAS_AF"])
    gnomad_nfe = str(row[f"gnomAD{genome_identifier}_NFE_AF"])
    gnomad_sas = str(row[f"gnomAD{genome_identifier}_SAS_AF"])

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


def main(in_file, out_file, sample_id, legacy_mode=False):
    with open(out_file, "w") as outfile:
        input_table = pd.read_csv(in_file, sep="\t", low_memory=False)
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

        input_table[["gnomAD", "gnomAD_sub"]] = input_table.apply(lambda x: pd.Series(get_gnomad_allele_frequency(x, legacy_mode)), axis=1)

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
    parser.add_argument("--in_file", type=str, dest="in_file", metavar="input.tsv", required=True, help="Annotated table that should be converted to a basic GSvar file.")
    parser.add_argument("--out_file", type=str, dest="out_file", metavar="output.GSvar", required=True, help="GSvar file containing the basic information needed to run the VariantRanking tool.")
    parser.add_argument("--sample_id", type=str, dest="sample_id", metavar="NA12878", required=False, help="Sample ID to get the sample information from the table.")
    parser.add_argument("--legacy", dest="legacy", action="store_true", required=False, help="Flag to use legacy mode. Needed to support tables from versions 1.1.0 to 1.2.0.")
    args = parser.parse_args()

    if args.sample_id is not None:
        sample_id = args.sample_id

    else:
        sample_id = "UNKNOWN"

    legacy_mode = args.legacy

    main(args.in_file, args.out_file, sample_id, legacy_mode)
