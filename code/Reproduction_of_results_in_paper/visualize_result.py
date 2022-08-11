import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


perf = pd.read_csv('../../result/Results_in_paper/df_prediction_performance.tsv', sep = '\t')
proba = pd.read_csv('../../result/Results_in_paper/df_approval_probability.tsv', sep = '\t')


### Box plot of prediction performance (AUPRC) for Monte-Carlo cross-validation ###
fig = sns.boxplot(x = 'Feature', y = 'AUPRC', order = ['OGE_CGE','OGE','CGE','Network','Expression','OGE_CGE_Network_Expression','Chemical','OGE_CGE_Network_Expression_Chemical'], data = perf)
plt.xticks(rotation = 45)
plt.xlabel('Feature')
plt.ylabel('AUPRC')
plt.tight_layout()
plt.savefig('../../result/Results_in_paper/prediction_performance_box_plot.eps', dpi = 300, format = 'eps')
plt.savefig('../../result/Results_in_paper/prediction_performance_box_plot.jpg', format = 'jpg')
plt.close()


### Box plots of approval probabilities for drugs by clinical phase ###
for feature in ['OGE_CGE','OGE_CGE_Network_Expression_Chemical']:
	fig = sns.boxplot(x = 'Phase', y = feature, order = ['Phase 1','Phase 2','Phase 3','Phase 4'], data = proba)
	plt.xlabel('Clinical phase')
	plt.ylabel('Approval probability')
	plt.tight_layout()
	plt.savefig('../../result/Results_in_paper/approval_probability_' + feature + '.eps', dpi = 300, format = 'eps')
	plt.savefig('../../result/Results_in_paper/approval_probability_' + feature + '.jpg', format = 'jpg')
	plt.close()

