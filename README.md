# 'bagb'
### R code to accompany "Bias analysis gone bad"

## File descriptions
- **misclassification.R** implements probabilistic bias analysis to adjust for exposure misclassification based on Chien _et. al._

- **selection_orsel.R** implements probabilistic bias analysis to adjust for selection bias based on DiForti _et. al._, using exposure prevalences in non-participating cases and controls as the bias parameters.

- **selection_prevs.R** implements probabilistic bias analysis to adjust for selection bias based on DiForti _et. al._, using OR<sub>sel</sub> as the bias factor (based on selection proportions).

## Required R packages
- tidyverse, reshape2, MASS, RColorBrewer

## Contact:
<02tahern@med.uvm.edu>