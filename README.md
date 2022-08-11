# DrugApprovalPrediction
## Description
Source codes for generating results of "Drug approval prediction based on the discrepancy of gene perturbation effects between cells and the human population".


## Requirements
### python
- python (3.7.9)
- pandas (1.3.4)
- numpy (1.21.6)
- scipy (1.6.0)
- scikit-learn (0.24.2)
- networkx (2.4)
- gprofiler-official (1.0.0)
- matplotlib (3.3.1)
- seaborn (0.11.0)
- tqdm (4.59.0)

### R
- R (4.0.5)
- ChemmineR (3.42.2)
- ChemmineOB (1.28.4)
- rcdk (3.6.0)
- rcdklibs (2.3)
- rJava (1.0-6)
- Rcpi (1.26.0)
- randomForest (4.6-14)
- ggplot2 (3.3.6)
- gridExtra (2.3)
- data.table (1.14.2)

## Installation
All python packages can be installed via pip (https://pypi.org/project/pip/).
pip install [package name]
e.g. pip install pandas

All R packages can be installed via BiocManager (https://github.com/Bioconductor/BiocManager).
Details in PrOCTOR (https://github.com/kgayvert/PrOCTOR).

Gnerally, a couple of minutes is needed for installing each package.

## DrugApprovalPrediction
### Codes for results reproduction
All codes for results reproduction is provided under the'./code/Reproduction_of_results_in_paper' folder.
- "make_gene_info.py" to generate data table organizing drug target information values for genes
- "make_ml_dataset.py" to generate datasets of machine learning for DrugApprovalPrediction 
- "run_ml.py" to make prediction for drug approval probability using Monte-Carlo cross-validation
- "visualize_result.py" to visualize the results from "run_ml.py"
- "run_all.py" to run the above code at once.

- Expected results (prediction performance & approval probability) are provided under './result/Results_in_paper' folder.
- The expected run time is under 90 minutes for running 1,000 times Monte-Carlo cross-validation.

```
To make Monte-Carlo cross-validation predictions, run 'run_all.py' under the './code/Reproduction_of_results_in_paper' folder.

Use the following commnad line to run 'run_all.py' in linux.

   $ python3 run_all.py
   
Prediction results are generated under the './result/Results_in_paper' folder.
```

### New molecule approval prediction
- Code (run_prediction_new_molecule.py) for predicting new molecule approval probability is provided under the './code/New_molecule_prediction' folder.
- Result of predicted approval probability is provided under './result/New_molecule_prediction' folder.
- Examples of input files (target.tsv & chemical.tsv) were provided under the './code/New_molecule_prediction' folder. 
- Either target information file or chemical information file should be necessary.
- Target genes of new molecule should be annotated with Ensembl gene id (see './code/New_molecule_prediction/target.tsv').

```
from run_prediction_new_molecule import *

make_new_molecule_dataset(target = './target.tsv', chemical = './chemical.tsv')     # 'New_molecule_ML_dataset.tsv' is saved under './data' folder.

prediction('./data/New_molecule_ML_dataset.tsv')     # 'New_molecule_approval_probability.tsv' is saved under './result/New_molecule_prediction' folder.
```

- R scripts for extracting chemical information, which are provided under the './code/New_molecule_prediction/PrOCTOR', are downloaded from Gayvert et al. (https://github.com/kgayvert/PrOCTOR).
```
source('PrOCTOR.R')
getStructureFeatures(SMILE = 'CC(=O)Nc1ccc(O)cc1')
```
