from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
from pprint import pprint
import pickle
#import vcfpy

random_seed = 14038

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

train_data = pd.read_csv('/home/dominic/Masterarbeit/eDiVA_new/data/train_set.csv', sep='\t')

train_data['Cadd2'].fillna(0, inplace=True)
train_data['Condel'].fillna(0, inplace=True)
train_data['SegMentDup'].fillna(0, inplace=True)
train_data['PrimatesPhyloP'].fillna(0, inplace=True)
train_data['PlacentalMammalPhyloP'].fillna(0, inplace=True)
train_data['PrimatesPhastCons'].fillna(0, inplace=True)
train_data['PlacentalMammalPhastCons'].fillna(0, inplace=True)
train_data['Eigen_Phred'].fillna(0, inplace=True)
train_data['MaxAF'].fillna(0, inplace=True)
train_data['MutAss'].fillna(-5, inplace=True)
train_data['ABB_score'].fillna(1, inplace=True)

training_labels = np.asarray(train_data['rank'])
training_features = np.asarray(train_data[['Cadd2','Condel','SegMentDup','PrimatesPhyloP','PlacentalMammalPhyloP', 'PrimatesPhastCons', 'PlacentalMammalPhastCons', 'Eigen_Phred','MaxAF','MutAss','ABB_score']])

train_features, test_features, train_labels, test_labels = train_test_split(training_features, training_labels, test_size = 0.2, random_state = random_seed)

test_data_ind = pd.read_csv('/home/dominic/Masterarbeit/eDiVA_new/data/test_set.csv', sep='\t')
test_data_ind['Cadd2'].fillna(0, inplace=True)
test_data_ind['Condel'].fillna(0, inplace=True)
test_data_ind['SegMentDup'].fillna(0, inplace=True)
test_data_ind['PrimatesPhyloP'].fillna(0, inplace=True)
test_data_ind['PlacentalMammalPhyloP'].fillna(0, inplace=True)
test_data_ind['PrimatesPhastCons'].fillna(0, inplace=True)
test_data_ind['PlacentalMammalPhastCons'].fillna(0, inplace=True)
test_data_ind['Eigen_Phred'].fillna(0, inplace=True)
test_data_ind['MaxAF'].fillna(0, inplace=True)
test_data_ind['MutAss'].fillna(-5, inplace=True)
test_data_ind['ABB_score'].fillna(1, inplace=True)

test_ind_labels = np.asarray(test_data_ind['rank'])
test_ind_features = np.asarray(test_data_ind[['Cadd2','Condel','SegMentDup','PrimatesPhyloP','PlacentalMammalPhyloP', 'PrimatesPhastCons', 'PlacentalMammalPhastCons', 'Eigen_Phred','MaxAF','MutAss','ABB_score']])

# print(train_data)
# print(labels)
# print(training_features)

# print(training_features[0:5])

#clf = RandomForestClassifier(n_estimators=1000, random_state=random_seed)
#clf.fit(training_features, labels)

#base_model = RandomForestClassifier(random_state = random_seed, verbose = 1)
#base_model.fit(train_features, train_labels)
#base_accuracy = evaluate(base_model, test_features, test_labels)

parameter_grid = {"n_estimators": [1000]}

rf_clf = RandomForestClassifier(random_state = random_seed)
grid_search = GridSearchCV(estimator = rf_clf, param_grid = parameter_grid, cv = 5, n_jobs = -1, verbose = 2)
grid_search.fit(training_features, training_labels)

best_grid = grid_search.best_estimator_
#grid_accuracy =evaluate(best_grid, test_features, test_labels)

#print("Improvement of {:0.2f}%.".format(100 * (grid_accuracy - base_accuracy) / base_accuracy))

#print(clf.feature_importances_)

export_file = '/home/dominic/Masterarbeit/eDiVA_new/data/rf-model_origTrainSet.pkl'

if not export_file.endswith('.pkl'):
    export_file = export_file + '.pkl'

pickle.dump(best_grid, open(export_file, 'wb'))
