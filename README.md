# Do your exposures need supervision?

Jenna R Krall, Howard H. Chang, W. Michael Caudle, Katherine M. Gass, Matthew J. Strickland

These files contain all functions and analysis code used to submit to the NIEHS epi-stats workshop in July 2015

# Code

All code was written using R version 3.1

- `niehs-epistats-results.Rmd` contains all the analysis, including results and figures, generated for the poster and oral presentation.
- `niehs-epistats-fn.R` contains all the functions needed to run the analysis.  

## Our results

To obtain our results, first install the `rmarkdown` package from CRAN.

`library(rmarkdown)`

`render("niehs-epistats-results.Rmd")`

## Using our functions

- use `pca1` to perform the PCA approach
- use `run_crt` to run the C&RT model
- use `niehs_outer` to run both approaches and obtain relevant output

 
