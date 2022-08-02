import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from tqdm import trange

def make_dataset(**kwargs):
	dataset_drug_list = []
	col_list = []
	for key, value in kwargs.items():
		if(key == 'target'):
			target = pd.read_csv(value, sep = '\t')
			target_val = pd.read_csv('../data/target_value.tsv', sep = '\t')
			drug_list = []
			print('filtering molecule')
			for i in trange(len(target)):
				dt = target.iloc[i,1].split('|')
				if(len(dt) * 0.9 <= len(target_val.loc[target_val['gene_id'].isin(dt) == True, :])):
					drug_list.append(target.iloc[i,0])
			target_filter = target.loc[target.iloc[:,0].isin(drug_list) == True, :]
			if(len(dataset_drug_list) == 0):
				dataset_drug_list.extend(drug_list)
			else:
				dataset_drug_list = list(set(dataset_drug_list) & set(drug_list))
			col_list.extend(['OGE','CGE','Network','Expression'])
		elif(key == 'chemical'):
			chemical = pd.read_csv(value, sep = '\t')
			chemical = chemical.dropna()
			if(len(dataset_drug_list) == 0):
				dataset_drug_list.extend(chemical.iloc[:,0].tolist())
			else:
				dataset_drug_list = list(set(dataset_drug_list) & set(chemical.iloc[:,0].tolist()))
			col_list.extend(['MolecularWeight','XLogP','HydrogenBondDonorCount','HydrogenBondAcceptorCount','PolarSurfaceArea','FormalCharge','NumRings','RotatableBondCount','Refractivity','Ro5','Ghose','Veber','wQED'])

	df_dataset = pd.DataFrame(data = {'Molecule' : dataset_drug_list}, columns = ['Molecule'] + col_list)
	print('making dataset')
	for i in trange(len(df_dataset)):
		if(len(set(['OGE','CGE','Network','Expression']) & set(df_dataset.columns))):
			dt_temp = target.loc[target.iloc[:,0] == df_dataset.iloc[i,0], :].iloc[0,1].split('|')
			for col in ['OGE','CGE','Network','Expression']:
				df_dataset.iloc[i,list(df_dataset.columns).index(col)] = target_val.loc[target_val['gene_id'].isin(dt_temp) == True, col].mean()
		if(len(set(['MolecularWeight','XLogP','HydrogenBondDonorCount','HydrogenBondAcceptorCount','PolarSurfaceArea','FormalCharge','NumRings','RotatableBondCount','Refractivity','Ro5','Ghose','Veber','wQED']) & set(df_dataset.columns))):
			for col in ['MolecularWeight','XLogP','HydrogenBondDonorCount','HydrogenBondAcceptorCount','PolarSurfaceArea','FormalCharge','NumRings','RotatableBondCount','Refractivity','Ro5','Ghose','Veber','wQED']:
				df_dataset.iloc[i,list(df_dataset.columns).index(col)] = chemical.loc[chemical.iloc[:,0] == df_dataset.iloc[i,0], [col]].iloc[0,0]
	df_dataset.to_csv('../../result/new_molecule_prediction/New_molecule_dataset.tsv', sep = '\t', header = True, index = False)


def prediction(dataset_path):
	df = pd.read_csv(dataset_path, sep = '\t')
	feature_list = list(set(df.columns) & set(['OGE','CGE','Network','Expression'] + ['MolecularWeight','XLogP','HydrogenBondDonorCount','HydrogenBondAcceptorCount','PolarSurfaceArea','FormalCharge','NumRings','RotatableBondCount','Refractivity','Ro5','Ghose','Veber','wQED']))

	RF = RandomForestClassifier(n_estimators = 1000)
	
	train = pd.read_csv('../data/ML_drug_set.tsv', sep = '\t')
	train.loc[train['Drug status'] == 'unapproved', 'Drug status'] = 1
	train.loc[train['Drug status'] == 'approved', 'Drug status'] = 0
	train['Drug status'] = train['Drug status'].astype(int)
	X_train = train.loc[:,feature_list]
	y_train = train.loc[:,'Drug status']

	RF.fit(X_train, y_train)
	y_proba = RF.predict_proba(df.loc[:,feature_list])
	
	df_proba = pd.DataFrame(data = {'Molecule' : df['Molecule'].tolist(), 'Approval probability' : list(y_proba[:,0])}, columns = ['Molecule','Approval probability', 'Predicted drug status'])
	df_proba.loc[df_proba['Approval probability'] >= 0.5, 'Predicted drug status'] = 'Approved'
	df_proba.loc[df_proba['Approval probability'] < 0.5, 'Predicted drug status'] = 'Unapproved'
	df_proba.to_csv('../../result/new_molecule_prediction/New_molecule_approval_probability.tsv', sep = '\t', header = True, index = False)
