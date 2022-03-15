# SARTCORS_MI
Implementing multiple imputation by chained equations (mice) in R on the [Society for Assisted Reproductive Technology Clinical Outcome Reporting System](https://www.sart.org/professionals-and-providers/research/) (SART CORS) database for a project investigating the impact of state insurance mandates on racial and ethnic disparities in assisted reproductive technology utilization and outcomes.

Workflow:
- run `create_imputed_datasets.R` x20 in parallel 
- run `assess_imputed_convergence.R` and review plots and summaries to assess convergence
- run `fit_modelsOR_MI.R` to fit the GEE models on each imputed dataset
- run `create_pooledOR_MI.R` to compute pooled estimates (ORs, 95% CIs)

Stef van Buuren's [Flexible Imputation of Missing Data](https://stefvanbuuren.name/fimd/) is an excellent resource for multiple imputation.
Documentation on R's [mice package](https://cran.r-project.org/web/packages/mice/mice.pdf) and the [corresponding vignettes](https://www.gerkovink.com/miceVignettes/Ad_hoc_and_mice/Ad_hoc_methods.html) are excellent resources for implementing multiple imputation by chained equations in R.
