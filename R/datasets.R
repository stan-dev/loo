#' Datasets for loo examples and vignettes
#'
#' Small datasets for use in **loo** examples and vignettes. The `Kline`
#' and `milk` datasets are also included in the **rethinking** package
#' (McElreath, 2016a), but we include them here as **rethinking** is not
#' on CRAN.
#'
#' @name loo-datasets
#' @aliases Kline milk
#'
#' @details
#' Currently the data sets included are:
#' * `Kline`:
#'   Small dataset from Kline and Boyd (2010) on tool complexity and demography
#'   in Oceanic islands societies. This data is discussed in detail in
#'   McElreath (2016a,2016b). [(Link to variable descriptions)](https://www.rdocumentation.org/packages/rethinking/versions/1.59/topics/Kline)
#' * `milk`:
#'   Small dataset from Hinde and Milligan (2011) on primate milk
#'   composition.This data is discussed in detail in McElreath (2016a,2016b).
#'   [(Link to variable descriptions)](https://www.rdocumentation.org/packages/rethinking/versions/1.59/topics/milk)
#'
#' @references
#' Hinde and Milligan. 2011. *Evolutionary Anthropology* 20:9-23.
#'
#' Kline, M.A. and R. Boyd. 2010. *Proc R Soc B* 277:2559-2564.
#'
#' McElreath, R. (2016a). rethinking: Statistical Rethinking book package.
#' R package version 1.59.
#'
#' McElreath, R. (2016b). *Statistical rethinking: A Bayesian course with
#' examples in R and Stan*. Chapman & Hall/CRC.
#'
#' @examples
#' str(Kline)
#' str(milk)
#'
NULL
