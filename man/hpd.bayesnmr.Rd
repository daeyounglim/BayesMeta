% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hpd.bayesnmr.R
\name{hpd.bayesnmr}
\alias{hpd.bayesnmr}
\title{get the highest posterior density (HPD) interval}
\usage{
\method{hpd}{bayesnmr}(object, prob = 0.95)
}
\arguments{
\item{object}{the output model from fitting a meta analysis/regression model}

\item{prob}{the probability which the HPD interval will cover}
}
\value{
dataframe containing HPD intervals for the parameters
}
\description{
get the highest posterior density (HPD) interval
}
