## How to use the package

### Installation

The package runs C++ behind the scenes, and therefore the computer on which this package will be used must have C++ compilers installed beforehand. Otherwise, the installation might fail.

* **For Windows users** - Download and install [Rtools from the website](https://cran.r-project.org/bin/windows/Rtools/).
* **For Mac users** - Please install the latest Xcode. It comes with all related compilers.
* **For Linux** - Install `r-base-dev` through the package manager of your system. For example, on Ubuntu 18.04, run `apt install r-base-dev`. 

After the compilers are ready, run the following to install the dependencies:
```r
install.packages("Rcpp", "RcppArmadillo", "RcppProgress", "BH")
```
When the packages above are successfully installed, run the following to install the `BayesMeta` packages:
```r
install.packages("BayesMeta_1.0.tar.gz", type="source", repo=NULL)
```

Instead of building the source file yourself, you can use the `devtools` package to install directly from the version control services like GitHub or GitLab. When using `devtools::install_github()`, sometimes `R` complains of the following:
```
tar: This does not look like a tar archive
gzip: stdin: unexpected end of file
```
In that case, check if `getOption("download.file.method")` is `curl`. Then, run `options("download.file.method" = "libcurl")` and then try installing again through `devtools::install_github("daeyounglim/BayesMeta")`.

### Example Code
```r
##############
# BNMR-HTMRe #
##############
groupinfo <- list(c(0,1), c(2,3), c(4)) # define the variance structure
fit <- bayesnmr(df$y, df$sd, x, df$ids, df$iarm, df$npt, groupinfo, prior = list(c01=1.0e05, c02=4, nu=3), mcmc=list(ndiscard=2500,nskip=1,nkeep=10000), init = list(beta = c(beta_true, gamma_true), sig2 = sig2_true))

## Get the estimates
##  (will be made into 'fitted' function soon)
beta.est <- rowMeans(fit$mcmc.draws$beta)
phi.est <- rowMeans(fit$mcmc.draws$phi)

## Get 'goodness of fit (gof)'
dic <- gof(fit, type = "dic") # calculate DIC
lpml <- gof(fit, type = "lpml") # calculate LPML
```

### Inputs explained for `bayesnmr`

* `y` - aggregate mean of the responses for each arm of each study
* `sd` - standard deviation of the responses for each arm of each study
* `x` - aggregate covariates for the mean component
* `ids` - study number in integers
* `iarm` - arm number of each study
* `groupinfo` - list of grouping information; the control(baseline) group must start from 0; the aggregate covariates `z` explaining the variance of the random effect of the t-th treatment will be construct based on this grouping information
* `npt` - number of observations per trial
* `nT` - number of treatments
* `prior` - list of hyperparameters; when not given, algorithm will run in default setting
* `mcmc` - list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning (nskip), and posterior sample size (nkeep)
* `add.z` - additional covariates other than the grouping vectors that should be column-concatenated to `z`. This should have the same number of rows as `y`, and `x`
* `scale.x` - logical variable for scaling x. Defaulting to `TRUE`. If not scaled, the `gamma[1]` cannot be interpreted as placebo
* `verbose` - logical variable for printing progress bar. Default to `FALSE`.
* `init` - initial values for beta (`ns + nT` dimensional) and phi. Dimensions must be conformant.

## `R` Package Development Reminders

* When creating an R package, don't forget to register the native symbols by running the following:
```r
tools::package_native_routine_registration_skeleton(".", "src/init.c", character_only=FALSE)
```

* When put datasets in the package, use `usethis::use_data_raw()` to create `data-raw/` folder and after preparing the dataset, say `df`, run the following:
```r
devtools::use_data(df)
```
Then, the convention is to document the datasets in `R/data.R` in the format of
```r
#' Data from Experiment 1
#'
#' This is data from the first experiment ever to try XYZ using Mechanical
#' Turk workers.
#'
#' @format A data frame with NNNN rows and NN variables:
#' \describe{
#'   \item{subject}{Anonymized Mechanical Turk Worker ID}
#'   \item{trial}{Trial number, from 1..NNN}
#'   ...
#' }
"df"
```

* Create `R/global.R` and include all global variables that will be used in the documentation. This includes all 

## To do

- [x] MCMC for ***Bayesian Network Meta-Regression Hierarchical Models Using Heavy-Tailed Multivariate Random Effects with Covariate-Dependent Variances (BNMR-HTMRe)***
- [x] DIC and LPML for BNMR-HTMRe
- [ ] R code for network metadata diagram
- [ ] R functions such as `fitted`, `print`, `summary`, `plot` for each class
- [ ] MCMC for ***Li, H., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., Lin, J., & Tershakovec, A. M. (2019). Bayesian inference for network meta-regression using multivariate random effects with applications to cholesterol lowering drugs. Biostatistics, 20(3), 499-516.***
- [ ] MCMC for ***Bayesian Flexible Hierarchical Skew Heavy-Tailed Multivariate Meta Regression Models for Individual Patient Data with Applications***



