import pandas as pd
from tqdm import trange
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from gprofiler import GProfiler

ess = pd.read_csv('../../data/Gene_info.tsv', sep = '\t')
ess = ess.loc[(ess['CGE'].isnull() != True) & (ess['OGE'].isnull() != True), :]
gp = GProfiler(return_dataframe = True)
conv = gp.convert(organism = 'hsapiens', query = ess['gene_id'].tolist(), target_namespace = 'ENSG')
conv.loc[conv['converted'] == 'None', 'converted'] = float('nan')
ess['hg38_ensgid'] = conv['converted'].tolist()
ess = ess.loc[ess['hg38_ensgid'].isnull() != True, :]
net = pd.read_csv('../../data/string_physical_centrality.tsv', sep = '\t', header = 0)
exp = pd.read_csv('../../data/consensus_tissue_epxression.tsv', sep = '\t', header = 0)
exp = exp.loc[exp['Exp_level_mean'] != 0, :]

scaler = MinMaxScaler()

net['Network_rank_score'] = (scaler.fit_transform(pd.DataFrame(net['Degree Centrality'].rank(method = 'dense'))) + scaler.fit_transform(pd.DataFrame(net['Betweenness Centrality'].rank(method = 'dense')))) / 2

exp['Expression_rank_score'] = (scaler.fit_transform(pd.DataFrame(exp['Exp_level_mean'].rank(method = 'dense'))) + scaler.fit_transform(pd.DataFrame(exp['Exp_broadness'].rank(method = 'dense')))) / 2

co_gene = list(set(ess['hg38_ensgid'].tolist()) & set(net['hg38_ensgid'].tolist()) & set(exp['hg38_ensgid'].tolist()))
df_score = pd.DataFrame(data = {'hg38_ensgid' : co_gene})

df_score = pd.merge(df_score, ess.loc[:,['hg38_ensgid','gene_id','CGE','OGE']], on = 'hg38_ensgid', how = 'left')
df_score = pd.merge(df_score, net.loc[:,['hg38_ensgid','Network_rank_score']], on = 'hg38_ensgid', how = 'left')
df_score = pd.merge(df_score, exp.loc[:,['hg38_ensgid','Expression_rank_score']], on = 'hg38_ensgid', how = 'left')

df_score.to_csv('../../data/target_value.tsv', sep = '\t', header = True, index = False)

