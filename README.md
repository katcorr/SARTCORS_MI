# SARTCORS_MI
This repository contains files relevant for implementing multiple imputation by chained equations (mice) in R on the [Society for Assisted Reproductive Technology Clinical Outcome Reporting System](https://www.sart.org/professionals-and-providers/research/) (SART CORS) database for a project investigating the impact of state insurance mandates on racial and ethnic disparities in assisted reproductive technology utilization and outcomes.

Fields imputed:
- patient race and ethnicity (Hispanic/Latinx, non-Hispanic white, non-Hispanic Black/African American, non-Hispanic Asian, other race, multiracial)
- BMI category (<18.5, 18.5-24.9, 25-25.9, 30+)
- FSH category (<=10, >10)
- AMH category (<1, 1-4, >4)
- number of embryos transferred
- parous
- any prior spontaneous abortion

Workflow:
- run `create_imputed_datasets.R` x20 in parallel 
- run `assess_imputed_convergence.R` and review plots and summaries to assess convergence
- run `fit_modelsOR_MI.R` to fit the GEE models on each imputed dataset
    - note that the gee-1 function used to fit the models is available in [this GitHub repository](https://github.com/katcorr/GEE1-Rfunction)
- run `compute_pooledOR_MI.R` to compute pooled estimates (ORs, 95% CIs)

Stef van Buuren's [Flexible Imputation of Missing Data](https://stefvanbuuren.name/fimd/) is an excellent resource for multiple imputation.
Documentation on R's [mice package](https://cran.r-project.org/web/packages/mice/mice.pdf) and the [corresponding vignettes](https://www.gerkovink.com/miceVignettes/) are useful resources for implementing multiple imputation by chained equations in R.
