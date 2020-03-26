# 'bagb'
### R code to accompany "Bias analysis gone bad"

## File descriptions

- **confounding.R** implements a simple bias analysis for an unmeasured confounder based on Mørch _et. al._

- **misclassification.R** implements [Greenland's approach](https://pubmed.ncbi.nlm.nih.gov/3043623/?from_single_result=3043623) to adjust for exposure misclassification based on Chien _et. al._ This analysis can also be performed in a Microsoft Excel spreadsheet available [here](https://sites.google.com/site/biasanalysis/exposure-misclassification-variance-correction).

- **selection.R** implements probabilistic bias analysis to adjust for selection bias based on DiForti _et. al._, using OR<sub>sel</sub> as the bias factor (based on selection proportions).

## Required R packages
- none (base R only)

## Contact:
<02tahern@med.uvm.edu>