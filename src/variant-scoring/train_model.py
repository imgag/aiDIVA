from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
import pandas as pd
import numpy as np
import pickle
import vcfpy



# exit(0)
#
#
# for record in reader:
#     if not record.is_snv():
#         continue

    # print("\nRecord:")
    # record_entries = record.INFO.get("CSQ")
    #
    # for entry in record_entries:
    #     entry_line = ""
    #
    #     splitted_entry = entry.split("|")
    #
    #     unique_consequence.add(splitted_entry[1])
    #     unique_impact.add(splitted_entry[2])
    #     unique_featuretype.add(splitted_entry[6])
    #     unique_biotype.add(splitted_entry[24])
    #
    #     for i in range(len(splitted_header)):
    #         entry_line += splitted_header[i] + ": " + splitted_entry[i]
    #
    #         if i < len(splitted_header):
    #             entry_line += ";  "

        # print(entry)
        # print(entry_line)
    # print('\n')

# print(unique_consequence)
# print(unique_impact)
# print(unique_featuretype)
# print(unique_biotype)

# vcf_file = VariantFile('/home/dominic/PycharmProjects/Masterarbeit/res/trio_variant_call_simulated_variant_dominant-denovo_vep.vcf')

# count = 0
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
#     print(rec.info['CSQ'])
#     count += 1
#     if count == 1:
#         break

train_data = pd.read_csv('/home/dominic/PycharmProjects/Masterarbeit/res/train_set.csv', sep='\t')

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

export_file = '/home/dominic/PycharmProjects/Masterarbeit/res/test-model.pkl'

if not export_file.endswith('.pkl'):
    export_file = export_file + '.pkl'

# pickle.dump(clf, open(export_file, 'wb'))
