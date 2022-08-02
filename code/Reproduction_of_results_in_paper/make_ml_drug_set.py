import pandas as pd
from tqdm import trange

drug = pd.read_csv('../../data/Drug_info.tsv', delimiter = '\t', header = 0)

target = pd.read_csv('../../data/target_value.tsv', sep = '\t', header = 0)

link = pd.read_csv('../../data/stitch_link_700_final.txt', sep = '\t', header = 0)
link = link.loc[link['chemical'].isin(drug['STITCH ID'].tolist()) == True, :]

chem = pd.read_csv('../../data/Chemical_properties.tsv', sep = '\t', header = 0)
chem = chem.dropna()

drug_list = []
print('filtering drug with target information coverage over than 90% & chemical information')
for i in trange(len(drug)):
	target_list = link.loc[link['chemical'] == drug.iloc[i,list(drug.columns).index('STITCH ID')], 'ensgid'].tolist()
	if(len(target_list) * 0.9 <= len(target.loc[target['gene_id'].isin(target_list) == True, :]) and drug.iloc[i,list(drug.columns).index('STITCH ID')] in chem['STITCH ID'].tolist()):
		drug_list.append(drug.iloc[i,list(drug.columns).index('STITCH ID')])
drug_filter = drug.loc[drug['STITCH ID'].isin(drug_list) == True, :]

print('making DrugApprovalPrediction dataset')
df = pd.DataFrame(data = {'STITCH ID' : drug_filter['STITCH ID'].tolist(), 'Drug status' : drug_filter['Drug status'].tolist(), 'Phase' : drug_filter['Phase'].tolist()}, columns = ['STITCH ID','Drug status','Phase','CGE','OGE','Network','Expression'] + list(chem.columns[2:]))
for i in trange(len(df)):
	target_temp = target.loc[target['gene_id'].isin(link.loc[link['chemical'] == df.iloc[i,list(df.columns).index('STITCH ID')], 'ensgid'].tolist()) == True, :]
	for j in ['CGE','OGE','Network','Expression']:
		df.iloc[i,list(df.columns).index(j)] = target_temp.loc[:,j].mean()
	for j in list(chem.columns)[2:]:
		df.iloc[i,list(df.columns).index(j)] = chem.loc[chem['STITCH ID'] == df.iloc[i,list(df.columns).index('STITCH ID')], [j]].iloc[0,0]
df.to_csv('../../data/ML_drug_set.tsv', sep = '\t', header = True, index = False)
