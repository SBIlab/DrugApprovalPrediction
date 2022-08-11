import pandas as pd
from tqdm import trange
from sklearn.preprocessing import MinMaxScaler
from gprofiler import GProfiler

### Load genes with CGE and OGE, filtering out genes without CGE and OGE ###  
ess = pd.read_csv('../../data/gene_essentiality.tsv', sep = '\t')
ess = ess.loc[(ess['CGE'].isnull() != True) & (ess['OGE'].isnull() != True), :]

### For genes with CGE and OGE, hg19 ensg ids were converted to hg38 ensg ids ###
gp = GProfiler(return_dataframe = True)
conv = gp.convert(organism = 'hsapiens', query = ess['gene_id'].tolist(), target_namespace = 'ENSG')
conv.loc[conv['converted'] == 'None', 'converted'] = float('nan')
ess['hg38_ensgid'] = conv['converted'].tolist()
ess = ess.loc[ess['hg38_ensgid'].isnull() != True, :]

### Load genes with centrality of protein interaction network ###
net = pd.read_csv('../../data/string_physical_centrality.tsv', sep = '\t')

### Load genes with tissue expression profiles, filtering out genes, which were not expressed in any tissues ###
exp = pd.read_csv('../../data/consensus_tissue_epxression.tsv', sep = '\t')
exp = exp.loc[exp['Exp_level_mean'] != 0, :]

### Two network centrality values and two expression profile values were integrated by averaging values and scaled from 0 to 1. ###
scaler = MinMaxScaler()

net['Network'] = (scaler.fit_transform(pd.DataFrame(net['Degree Centrality'].rank(method = 'dense'))) + scaler.fit_transform(pd.DataFrame(net['Betweenness Centrality'].rank(method = 'dense')))) / 2

exp['Expression'] = (scaler.fit_transform(pd.DataFrame(exp['Exp_level_mean'].rank(method = 'dense'))) + scaler.fit_transform(pd.DataFrame(exp['Exp_broadness'].rank(method = 'dense')))) / 2

### Generate a data table organizing CGE, OGE, network, and expression information ###
co_gene = list(set(ess['hg38_ensgid'].tolist()) & set(net['hg38_ensgid'].tolist()) & set(exp['hg38_ensgid'].tolist()))
df_score = pd.DataFrame(data = {'hg38_ensgid' : co_gene})

df_score = pd.merge(df_score, ess.loc[:,['hg38_ensgid','gene_id','CGE','OGE']], on = 'hg38_ensgid', how = 'left')
df_score = pd.merge(df_score, net.loc[:,['hg38_ensgid','Network']], on = 'hg38_ensgid', how = 'left')
df_score = pd.merge(df_score, exp.loc[:,['hg38_ensgid','Expression']], on = 'hg38_ensgid', how = 'left')

df_score.to_csv('../../data/Gene_info.tsv', sep = '\t', header = True, index = False)

