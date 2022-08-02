# DrugApprovalPrediction
## Description
Source codes for generating results of "Drug approval prediction based on the discrepancy of gene perturbation effects between cells and the human population".


## Requirements
### python
- python (3.7.9)
- pandas (1.3.4)
- numpy (1.18.5)
- scikit-learn (0.24.2)
- gprofiler-official (1.0.0)
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

All R packages can be installed via BiocManager (https://github.com/Bioconductor/BiocManager).

Gnerally, a couple of minutes is needed for installing each package.

## DrugApprovalPrediction
### Result reproduction
- Code (run_all.py) for reproducing drug target information, machine learning dataset, Monte-Carlo cross-validation is provided under the './code/Reproduction_of_results_in_paper' folder.
- Expected results (prediction performance & approval probability) are provided under './result/result_in_paper' folder.
- The expected run time is under 90 minutes for running 1,000 times Monte-Carlo cross-validation.

```
To make Monte-Carlo cross-validation predictions, run 'run_all.py' under the './code/Reproduction_of_results_in_paper' folder.

Use the following commnad line to run 'run_all.py' in linux.

   $ python3 run_all.py
   
Prediction results are generated under the './result/result_in_paper' folder.
```

### New molecule approval prediction
- Code (run_prediction_new_molecule.py) for predicting new molecules' approval probabilities is provided under the './code/New_prediction' folder.
- predicted approval probabilities are provided under './result/new_molecule_prediction' folder.
- Examples of input files (target.tsv & chemical.tsv) were provided under the './code/New_prediction' folder. 
- Target genes of new molecule should be annotated with Ensembl gene id (see './code/New_prediction/target.tsv').
- R scripts for extracting chemical information, which are provided under the './code/New_prediction/PrOCTOR', were downloaded from Gayvert et al. (https://github.com/kgayvert/PrOCTOR).
```
source('PrOCTOR.R')
getStructureFeatures(SMILE = 'CC(=O)Nc1ccc(O)cc1')
```
