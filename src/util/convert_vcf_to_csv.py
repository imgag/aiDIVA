import pandas as pd
import sys

import vcfpy

variant_consequences = {'transcript_ablation': 1,
                        'splice_acceptor_variant': 2,
                        'splice_donor_variant': 3,
                        'stop_gained': 4,
                        'frameshift_variant': 5,
                        'stop_lost': 6,
                        'start_lost': 7,
                        'transcript_amplification': 8,
                        'inframe_insertion': 9,
                        'inframe_deletion': 10,
                        'missense_variant': 11,
                        'protein_altering_variant': 12,
                        'splice_region_variant': 13,
                        'incomplete_terminal_codon_variant': 14,
                        'start_retained_variant': 15,
                        'stop_retained_variant': 16,
                        'synonymous_variant': 17,
                        'coding_sequence_variant': 18,
                        'mature_miRNA_variant': 19,
                        '5_prime_UTR_variant': 20,
                        '3_prime_UTR_variant': 21,
                        'non_coding_transcript_exon_variant': 22,
                        'intron_variant': 23,
                        'NMD_transcript_variant': 24,
                        'non_coding_transcript_variant': 25,
                        'upstream_gene_variant': 26,
                        'downstream_gene_variant': 27,
                        'TFBS_ablation': 28,
                        'TFBS_amplification': 29,
                        'TF_binding_site_variant': 30,
                        'regulatory_region_ablation': 31,
                        'regulatory_region_amplification': 32,
                        'feature_elongation': 33,
                        'regulatory_region_variant': 34,
                        'feature_truncation': 35,
                        'intergenic_variant': 36}


reader = vcfpy.Reader.from_path(sys.argv[1])
# for entry in reader.header.lines:
#     print(entry)

# print(reader.header.get_info_field_info("CSQ").description)

header = reader.header.get_info_field_info("CSQ").description
header = header.split(":")[1]
splitted_header = header.split('|')

# print(header)

# print(reader.parser.parse_next_record().INFO.get("CSQ"))

unique_consequence = set()
unique_impact = set()
unique_featuretype = set()
unique_biotype = set()

header_entries = ['Chr', 'Pos', 'Ref', 'Alt', 'Rank']
header_entries.extend(splitted_header)

test_variant_pandas = pd.DataFrame(columns=header_entries)

print("Extract data from VCF")
for record in reader:
    if not record.is_snv():
        continue

    # print(record.INFO.get('CSQ'))
    # print(record.INFO.get('RANK'))
    # print(len(record.INFO.get('CSQ')))

    record_entries = record.INFO.get('CSQ')

    record_entry = [record.CHROM, record.POS, record.REF, ','.join([alt.serialize() for alt in record.ALT]), record.INFO.get('RANK')[0]]
    target_entry = []
    current_consequence = 100

    for entry in record_entries:
        splitted_entry = entry.split('|')

        consequences = [variant_consequences[consequence] for consequence in splitted_entry[1].split('&')]

        if min(consequences) < current_consequence:
            current_consequence = min(consequences)
            target_entry = splitted_entry

    record_entry.extend(target_entry)
    # print(len(record_entry))
    record_dict = {header_entries[index]: record_entry[index] for index in range(len(header_entries))}

    test_variant_pandas = test_variant_pandas.append(record_dict, ignore_index=True)
    # print(test_variant_pandas)

print("Write data to CSV")
test_variant_pandas.to_csv(sys.argv[2], sep='\t', encoding='utf-8', index=False)
