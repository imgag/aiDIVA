from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
import pandas as pd
import numpy as np
import pickle
import vcfpy
import sys
import time


variant_consequences = {"transcript_ablation" : 1,
                        "splice_acceptor_variant" : 2,
                        "splice_donor_variant" : 3,
                        "stop_gained" : 4,
                        "frameshift_variant" : 5,
                        "stop_lost" : 6,
                        "start_lost" : 7,
                        "transcript_amplification" : 8,
                        "inframe_insertion" : 9,
                        "inframe_deletion" : 10,
                        "missense_variant" : 11,
                        "protein_altering_variant" : 12,
                        "splice_region_variant" : 13,
                        "incomplete_terminal_codon_variant" : 14,
                        "start_retained_variant" : 15,
                        "stop_retained_variant" : 16,
                        "synonymous_variant" : 17,
                        "coding_sequence_variant" : 18,
                        "mature_miRNA_variant" : 19,
                        "5_prime_UTR_variant" : 20,
                        "3_prime_UTR_variant" : 21,
                        "non_coding_transcript_exon_variant" : 22,
                        "intron_variant" : 23,
                        "NMD_transcript_variant" : 24,
                        "non_coding_transcript_variant" : 25,
                        "upstream_gene_variant" : 26,
                        "downstream_gene_variant" : 27,
                        "TFBS_ablation" : 28,
                        "TFBS_amplification" : 29,
                        "TF_binding_site_variant" : 30,
                        "regulatory_region_ablation" : 31,
                        "regulatory_region_amplification" : 32,
                        "feature_elongation" : 33,
                        "regulatory_region_variant" : 34,
                        "feature_truncation" : 35,
                        "intergenic_variant" : 36}


print(sys.argv[1])

reader = vcfpy.Reader.from_path(sys.argv[1])

#for entry in reader.header.lines:
    #print(entry)

#print(reader.header.get_info_field_info("CSQ").description)

header = reader.header.get_info_field_info("CSQ").description
header = header.split(":")[1]
splitted_header = header.split("|")
#print(header)

#print(reader.parser.parse_next_record().INFO.get("CSQ"))

unique_consequence = set()
unique_impact = set()
unique_featuretype = set()
unique_biotype = set()

test_list = ['Chr','Pos','Ref','Alteration']
#print(test_list)
#print(type(splitted_header))
test_list.extend(splitted_header)
#print(test_list)
variant_data = pd.DataFrame(columns=test_list)

counter = 0
for record in reader:
    #t = time.process_time()
    if not record.is_snv():
        continue

    #print("\nRecord:")
    
    #print(header)
    #print(record.CHROM)
    #print(record.POS)
    #print(record.REF)
    #print(record.ALT)
    
    record_entries = record.INFO.get("CSQ")
    
    record_entry = [record.CHROM, record.POS, record.REF, record.ALT]
    target_entry = []
    currentConsequence = 100
    
    for entry in record_entries:
        entry_line = ""

        splitted_entry = entry.split("|")
        
        #print(splitted_entry)
        
        #print(splitted_entry[1])
        
        consequences = [variant_consequences[consequence] for consequence in splitted_entry[1].split('&')]
        #print(consequences)
        if min(consequences) < currentConsequence:
            currentConsequence = min(consequences)
            target_entry = splitted_entry
        #print(min(consequences))
        #for consequence in splitted_entry[1].split("&"):
        #    print(variant_consequences[consequence])
        
        unique_consequence.add(splitted_entry[1])
        
        #if splitted_entry[2] == "HIGH":
        #    print(splitted_entry[1].split("&"))
            
        unique_impact.add(splitted_entry[2])
        unique_featuretype.add(splitted_entry[6])
        unique_biotype.add(splitted_entry[24])
        
        for i in range(len(splitted_header)):
            entry_line += splitted_header[i] + ": " + splitted_entry[i]

            if i < len(splitted_header):
                entry_line += ";  "

        break
        #print(entry)
        #print(entry_line)
        
    record_entry.extend(target_entry)
    #print(record_entry)
    test_dict = {test_list[x] : record_entry[x] for x in range(len(test_list))}
    #print(test_dict)
    variant_data = variant_data.append(test_dict, ignore_index=True)
    
    #print("\n")
    counter += 1
    print(counter)
    #elapsed_time = time.process_time() - t
    #print(elapsed_time)

variant_data.replace('', np.nan, inplace=True)
variant_data.to_csv('variant_data_test_full.csv', sep='\t', encoding='utf-8', index=False)
print(variant_data)
print(unique_consequence)
#print(unique_impact)
#print(unique_featuretype)
#print(unique_biotype)

# vcf_file = VariantFile("/home/dominic/PycharmProjects/Masterarbeit/res/trio_variant_call_simulated_variant_dominant-denovo_vep.vcf")

count = 0
# for rec in vcf_file.fetch():
#     print(rec.info)
#     # print(rec.header)
#     print(list(rec.header.records))
#
#     for x in vcf_file.header.records:
#         print(x)
#         print(x.type, x.key)
#
#     print(rec.info.keys())
#     print(rec.info["CSQ"])
#     count += 1
#     if count == 1:
#         break

# Allele|Consequence|IMPACT|SYMBOL|HGNC_ID|Feature|Feature_type|EXON|INTRON|HGVSc|HGVSp|DOMAINS|SIFT|PolyPhen|Existing_variation|AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_EAS_AF|gnomAD_NFE_AF|gnomAD_SAS_AF|EA_AF|AA_AF|BIOTYPE|CADD_PHRED|REVEL|FATHMM_MKL_C|FATHMM_MKL_NC|MaxEntScan_ref|MaxEntScan_alt|GeneSplicer|ada_score|rf_score|gnomADg_AF|gnomADg_Hom|gnomADg_Hemi|REPEATMASKER|CLINVAR|CLINVAR_DETAILS|PHYLOP|OMIM|HGMD|HGMD_CLASS|HGMD_MUT|HGMD_GENE|HGMD_PHEN

train_data = pd.read_csv("/home/dominic/PycharmProjects/Masterarbeit/res/train_set.csv", sep="\t")

train_data['Cadd2'].fillna(0)
train_data['Condel'].fillna(0)
train_data['SegMentDup'].fillna(0)
train_data['PrimatesPhyloP'].fillna(0)
train_data['PlacentalMammalPhyloP'].fillna(0)
train_data['PrimatesPhastCons'].fillna(0)
train_data['PlacentalMammalPhastCons'].fillna(0)
train_data['Eigen_Phred'].fillna(0)
train_data['MaxAF'].fillna(0)
train_data['MutAss'].fillna(-5)
train_data['ABB_score'].fillna(1)

labels = np.asarray(train_data['rank'])
training_features = np.asarray(train_data[['Cadd2','Condel','SegMentDup','PrimatesPhyloP','PlacentalMammalPhyloP',
                                           'PrimatesPhastCons','PrimatesPhastCons','PlacentalMammalPhastCons',
                                           'Eigen_Phred','MaxAF','MutAss','ABB_score']])

# print(train_data)
# print(labels)
# print(training_features)

# print(training_features[0:5])

# clf = RandomForestClassifier(n_estimators=1000, random_state=SyncRNG(14038).seed)
# clf.fit(training_features, labels)

# print(clf.feature_importances_)

#export_file = "/home/dominic/PycharmProjects/Masterarbeit/res/test-model.pkl"

#if not export_file.endswith(".pkl"):
#        export_file = export_file + ".pkl"

# pickle.dump(clf, open(export_file, 'wb'))
