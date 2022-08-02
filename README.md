# DrugApprovalPrediction
## Description
Source codes for generating results of "Drug approval prediction based on the discrepancy of gene perturbation effects between cells and the human population".


## Requirements
- python (3.7.9)
- pandas (1.3.4)
- numpy (1.18.5)
- scikit-learn (0.24.2)
- gprofiler-official (1.0.0)
- tqdm (4.59.0)

## Installation
All packages can be installed via pip (https://pypi.org/project/pip/). Gnerally, a couple of minutes is needed for installing each package.

**e.g.** pip install pandas


## DrugApprovalPrediction
- Code (run_all.py) for reproducing drug target information, machine learning dataset, Monte-Carlo cross-validation is provided under the './code/Reproduction_of_results_in_paper' folder.
- Expected results (prediction performance & approval probability) are provided under './result/result_in_paper' folder.
- The expected run time is under 90 minutes for running 1,000 times Monte-Carlo cross-validation.



```
To make Monte-Carlo cross-validation predictions, run 'run_all.py' under the './code/Reproduction_of_results_in_paper' folder.

Use the following commnad line to run 'run_all.py' in linux.

   $ python3 run_all.py
   
Prediction results are generated under the './result/result_in_paper' folder.
```
