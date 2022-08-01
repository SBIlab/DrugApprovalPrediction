import pandas as pd
import numpy as np
import copy as cp
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import *
from sklearn.metrics import *
from tqdm import trange
from itertools import combinations
from multiprocessing import *
import os

def rf(feature):
	drug_list = []; status_list = []; pred_list = []; prob_list = []
	roc_auc_list = []; pr_auc_list = []; accuracy_list = []; sens_list = []; spec_list = []; prec_list = []; recall_list = []; f1_list = []; #f_dic = dict(zip(feature, [[] for x in range(len(feature))]))
	df_proba = pd.DataFrame(columns = ['STITCH ID','Drug status','Predicted status','Approval probability'])
	df_perf = pd.DataFrame(columns = ['ROC_AUC','AUPRC','Accuracy','Sensitivity','Specificity','Precision','Recall','F1_score'])
	
	for i in trange(1000):
		data_un = data.loc[data['Drug status'] == 1, :]; data_ap = data.loc[data['Drug status'] == 0, :]
		if(len(data_un) > len(data_ap)):
			data_un = data_un.sample(n = len(data_ap), random_state = i)
		else:
			data_ap = data_ap.sample(n = len(data_un), random_state = i)
		data_sum = pd.concat([data_un, data_ap])

		X_train, X_test, y_train, y_test = train_test_split(data_sum, data_sum['Drug status'], test_size = 0.1, random_state = i+1, stratify = data_sum['Drug status'])

		RF = RandomForestClassifier(n_estimators = 1000, n_jobs = 1)
		RF.fit(X_train.loc[:,feature], y_train)
		y_pred = RF.predict(X_test.loc[:,feature])
		y_prob = RF.predict_proba(X_test.loc[:,feature])

		drug_list.extend(list(X_test.iloc[:,list(X_test.columns).index('STITCH ID')])); status_list.extend(list(X_test.iloc[:,list(X_test.columns).index('Drug status')])); pred_list.extend(list(y_pred)); prob_list.extend(list(y_prob[:,1]))

		roc_auc = roc_auc_score(y_test, y_prob[:,1]); roc_auc_list.append(roc_auc)

		precision, recall, thresholds = precision_recall_curve(y_test, y_prob[:,1])
		pr_auc = auc(recall, precision); pr_auc_list.append(pr_auc);

		cm = confusion_matrix(y_test, y_pred)
		accuracy = accuracy_score(y_test, y_pred); sensitivity = cm[1,1] / (cm[1,0] + cm[1,1]); specificity = cm[0,0] / (cm[0,0] + cm[0,1]); precision = cm[1,1] / (cm[0,1] + cm[1,1]); recall = cm[1,1] / (cm[1,0] + cm[1,1]); f1 = 2 * precision * recall / (precision + recall)
		
		accuracy_list.append(accuracy); sens_list.append(sensitivity); spec_list.append(specificity); prec_list.append(precision); recall_list.append(recall); f1_list.append(f1)
	
	df_proba['STITCH ID'] = drug_list; df_proba['Drug status'] = status_list; df_proba['Predicted status'] = pred_list; df_proba['Approval probability'] = list(1-np.array(prob_list))
	df_proba.loc[df_proba['Drug status'] == 1, 'Drug status'] = 'unapproved'; df_proba.loc[df_proba['Predicted status'] == 1, 'Predicted status'] = 'unapproved'; df_proba.loc[df_proba['Drug status'] == 0, 'Drug status'] = 'approved'; df_proba.loc[df_proba['Predicted status'] == 0, 'Predicted status'] = 'approved'
	df_perf['ROC_AUC'] = roc_auc_list; df_perf['AUPRC'] = pr_auc_list; df_perf['Accuracy'] = accuracy_list; df_perf['Sensitivity'] = sens_list; df_perf['Specificity'] = spec_list; df_perf['Precision'] = prec_list; df_perf['Recall'] = recall_list; df_perf['F1_score'] = f1_list
	
	f_name = ''
	n = 0; m = 0
	for f in feature:
		if(f in list(data.columns)[7:]):
			m += 1
		else:
			if(n == 0):
				f_name = f
				n += 1
			else:
				f_name  = f_name + '_' + f
	if(m > 0):
		if(n == 0):
			f_name = 'Chemical'
		else:
			f_name = f_name + '_' + 'Chemical'
	
	df_proba['Feature'] = f_name
	df_proba_pivot = pd.pivot_table(df_proba, index = 'STITCH ID', columns = 'Feature', aggfunc = 'mean')
	df_proba_final = pd.DataFrame(data = {'STITCH ID' : list(df_proba_pivot.index), f_name : df_proba_pivot[('Approval probability',f_name)].tolist()}, columns = ['STITCH ID',f_name])
	df_proba_final.to_csv('../../result/df_proba_' + f_name + '.tsv', sep = '\t', header = True, index = False)

	df_perf['Feature'] = f_name
	df_perf.to_csv('../../result/df_performance_' + f_name + '.tsv', sep = '\t', header = True, index = False)


data = pd.read_csv('../../data/ML_drug_set.tsv', delimiter = '\t', header = 0)
feature_comb_list_temp = [['OGE','CGE'], ['OGE'], ['CGE'], ['OGE','CGE','Network','Expression'], ['Chemical'], ['OGE','CGE','Network','Expression','Chemical']]
feature_comb_list = []
for i in feature_comb_list_temp:
	if('Chemical' in i):
		temp_list = list(i); temp_list.remove('Chemical'); temp_list.extend(list(data.columns)[7:])
		feature_comb_list.append(tuple(temp_list))
	else:
		feature_comb_list.append(i)

data.loc[data.iloc[:,1] == 'unapproved', 'Drug status'] = 1
data.loc[data.iloc[:,1] == 'approved', 'Drug status'] = 0
data['Drug status'] = data['Drug status'].astype(int)

print('Machine learning start')
procs = []
for i in range(len(feature_comb_list)):
	procs.append(Process(target = rf, args = (feature_comb_list[i], )))
for proc in procs:
	proc.start()
for proc in procs:
	proc.join()

proba_file_list = [x for x in os.listdir('../../result/') if 'proba' in x]
perf_file_list = [x for x in os.listdir('../../result/') if 'performance' in x]

n = 0
for file_ in proba_file_list:
	if(n == 0):
		df_sum = pd.read_csv('../../result/' + file_, sep = '\t')
		df_sum = pd.merge(df_sum, data.loc[:,['STITCH ID','Drug status','Phase']], on = 'STITCH ID', how = 'left')
		df_sum.loc[df_sum['Drug status'] == 1, 'Drug status'] = 'unapproved'
		df_sum.loc[df_sum['Drug status'] == 0, 'Drug status'] = 'approved'
		df_sum = df_sum.iloc[:,[0,-2,-1,1]]
		n += 1
	else:
		df_temp = pd.read_csv('../../result/' + file_, sep = '\t')
		df_sum = pd.merge(df_sum, df_temp, on = 'STITCH ID', how = 'left')
	os.system('rm -rf ../../result/' + file_)
df_sum.to_csv('../../result/df_approval_probability.tsv', sep = '\t', header = True, index = False)

n = 0
for file_ in perf_file_list:
	if(n == 0):
		df_sum = pd.read_csv('../../result/' + file_, sep = '\t')
		n += 1
	else:
		df_temp = pd.read_csv('../../result/' + file_, sep = '\t')
		df_sum = pd.concat([df_sum, df_temp])
	os.system('rm -rf ../../result/' + file_)
df_sum.to_csv('../../result/df_prediction_performance.tsv', sep = '\t', header = True, index = False)
