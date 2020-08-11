# 📝 Table of Contents
+ [Installation](#installation)
+ [Example Code](#example_code)
+ [Inputs explained for `bayesnmr`](#inputs)
+ [`R` Package Development Reminders](#dev_reminders)
+ [To-do list](#to_do)

## 🔨 Installation <a name = "installation"></a>
### Simple
Run
```r
devtools::install_github("daeyounglim/BayesMeta")
```

### Complex
Click the "Clone or download" button and download ZIP. Run the following commands:
```
unzip("BayesMeta-master.zip")
file.rename("BayesMeta-master", "BayesMeta")
system("R CMD build BayesMeta")
install.packages(list.files(pattern="(BayesMeta)(.*)(.tar.gz)"), repos = NULL, type = "source")
```

### 🔧 Troubleshooting for potential errors
The package runs C++ behind the scenes, and therefore the computer on which this package will be used must have C++ compilers installed beforehand. Otherwise, the installation might fail.

* **For Windows users** - Download and install [Rtools from the website](https://cran.r-project.org/bin/windows/Rtools/).
* **For Mac users** - Please install the latest Xcode. It comes with all related compilers.
* **For Linux** - Install `r-base-dev` through the package manager of your system. For example, on Ubuntu 18.04, run `apt install r-base-dev`. While installing the `devtools` package, you might encounter error messages that flash by in fleeting moments that say the following:

```
------------------------- ANTICONF ERROR ---------------------------
Configuration failed because libcurl was not found.
Try installing:
 * deb: libcurl4-openssl-dev (Debian, Ubuntu, etc)
 * rpm: libcurl-devel (Fedora, CentOS, RHEL)
 * csw: libcurl_dev (Solaris)
If libcurl is already installed, check that 'pkg-config' is in your PATH and PKG_CONFIG_PATH contains a libcurl.pc file. If pkg-config is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:
R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'                                                                                       --------------------------------------------------------------------  
```
and
```
---------------------------------------------                                                                                    Unable to find the LibSSH2 (ver>=v1.8) library on this system. Building git2r without support for SSH transport.
To build with SSH support, please install:
    libssh2-1-dev (package on e.g. Debian and Ubuntu)
	libssh2-devel (package on e.g. Fedora, CentOS and RHEL)
    libssh2 (Homebrew package on OS X)
and try again.

If the LibSSH2 library is installed on your system but the git2r configuration is unable to find it, you can specify the include and lib path to LibSSH2 with:
 R CMD INSTALL git2r --configure-vars='LIBS=-L/path/to/libs CPPFLAGS=-I/path/to/headers'
---------------------------------------------
```
and
```
------------------------- ANTICONF ERROR ---------------------------
Configuration failed because libxml-2.0 was not found. Try installing:                            
 * deb: libxml2-dev (Debian, Ubuntu, etc)
 * rpm: libxml2-devel (Fedora, CentOS, RHEL)
 * csw: libxml2_dev (Solaris)
If libxml-2.0 is already installed, check that 'pkg-config' is in your PATH and PKG_CONFIG_PATH contains a libxml-2.0.pc file.
If pkg-config is unavailable you can set INCLUDE_DIR and LIB_DIR manually via:
R CMD INSTALL --configure-vars='INCLUDE_DIR=... LIB_DIR=...'                                                                                       --------------------------------------------------------------------    
```
Install the Linux libraries above according to your Linux distributions.

After the compilers are ready, try installing the package again. The dependency installation is supposed to be automated but if `R` complains, then run the following:
```r
install.packages("Rcpp", "RcppArmadillo", "RcppProgress", "BH")
```
When the packages above are successfully installed, run `devtools::install_github('daeyounglim/BayesMeta')` or following the complex version.

#### Linux distros
When using the `devtools::install_github`, some Linux distros fire the following error message:
```
tar: This does not look like a tar archive
gzip: stdin: unexpected end of file
```
In that case, check if `getOption("download.file.method")` is `curl`. Then, run `options("download.file.method" = "libcurl")` and then try installing again through `devtools::install_github("daeyounglim/BayesMeta")`.

## 💻 Example Code <a name="example_code"></a>
```r
##############
# BNMR-HTMRe #
##############
groupinfo <- list(c(0,1), c(2,3), c(4)) # define the variance structure
fit <- bayesnmr(df$y, df$sd, x, df$ids, df$iarm, df$npt, groupinfo, prior = list(c01=1.0e05, c02=4, nu=3), mcmc=list(ndiscard=2500,nskip=1,nkeep=10000), verbose=TRUE)

## Get the estimates
##  (will be made into 'fitted' function soon)
beta.est <- rowMeans(fit$mcmc.draws$beta)
phi.est <- rowMeans(fit$mcmc.draws$phi)

## Get 'goodness of fit (gof)'
dic <- gof(fit, type = "dic") # calculate DIC
lpml <- gof(fit, type = "lpml") # calculate LPML
```

## 🔣 Inputs explained for `bayesnmr` <a name="inputs"></a>

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

## 👨🏻‍💻 `R` Package Development Reminders <a name="dev_reminders"></a>

* When creating an R package, don't forget to register the native symbols by running the following:
```r
tools::package_native_routine_registration_skeleton(".", "src/init.c", character_only=FALSE)
```
However, if you're using `Rcpp`, do **NOT** use the above command. Simply place the `[[Rcpp::export]]` directive above the `C++` function you want exported. It will streamline the native symbol registration for you. Running the above command will only end up conflicting with the configuration generated by `Rcpp`.

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

## ✔️ To do <a name="to_do"></a>

- [x] MCMC for ***Bayesian Network Meta-Regression Hierarchical Models Using Heavy-Tailed Multivariate Random Effects with Covariate-Dependent Variances (BNMR-HTMRe)***
- [x] DIC and LPML for BNMR-HTMRe
- [ ] R code for network metadata diagram
- [ ] R code for SUCRA plot
- [x] R functions such as `fitted`, `print` for each class
- [ ] R functions such as `summary`, `plot` for each class
- [ ] MCMC for ***Li, H., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., Lin, J., & Tershakovec, A. M. (2019). Bayesian inference for network meta-regression using multivariate random effects with applications to cholesterol lowering drugs. Biostatistics, 20(3), 499-516.***
- [ ] MCMC for ***Bayesian Flexible Hierarchical Skew Heavy-Tailed Multivariate Meta Regression Models for Individual Patient Data with Applications***



