# Datasets for loo examples and vignettes

Small datasets for use in **loo** examples and vignettes. The `Kline`
and `milk` datasets are also included in the **rethinking** package
(McElreath, 2016a), but we include them here as **rethinking** is not on
CRAN.

## Details

Currently the data sets included are:

- `Kline`: Small dataset from Kline and Boyd (2010) on tool complexity
  and demography in Oceanic islands societies. This data is discussed in
  detail in McElreath (2016a,2016b). [(Link to variable
  descriptions)](https://www.rdocumentation.org/packages/rethinking/versions/1.59/topics/Kline)

- `milk`: Small dataset from Hinde and Milligan (2011) on primate milk
  composition.This data is discussed in detail in McElreath
  (2016a,2016b). [(Link to variable
  descriptions)](https://www.rdocumentation.org/packages/rethinking/versions/1.59/topics/milk)

- `voice`: Voice rehabilitation data from Tsanas et al. (2014).

## References

Hinde and Milligan. 2011. *Evolutionary Anthropology* 20:9-23.

Kline, M.A. and R. Boyd. 2010. *Proc R Soc B* 277:2559-2564.

McElreath, R. (2016a). rethinking: Statistical Rethinking book package.
R package version 1.59.

McElreath, R. (2016b). *Statistical rethinking: A Bayesian course with
examples in R and Stan*. Chapman & Hall/CRC.

A. Tsanas, M.A. Little, C. Fox, L.O. Ramig: Objective automatic
assessment of rehabilitative speech treatment in Parkinson's disease,
IEEE Transactions on Neural Systems and Rehabilitation Engineering, Vol.
22, pp. 181-190, January 2014

## Examples

``` r
str(Kline)
#> 'data.frame':    10 obs. of  5 variables:
#>  $ culture    : Factor w/ 10 levels "Chuuk","Hawaii",..: 4 7 6 10 3 9 1 5 8 2
#>  $ population : int  1100 1500 3600 4791 7400 8000 9200 13000 17500 275000
#>  $ contact    : Factor w/ 2 levels "high","low": 2 2 2 1 1 1 1 2 1 2
#>  $ total_tools: int  13 22 24 43 33 19 40 28 55 71
#>  $ mean_TU    : num  3.2 4.7 4 5 5 4 3.8 6.6 5.4 6.6
str(milk)
#> 'data.frame':    29 obs. of  8 variables:
#>  $ clade         : Factor w/ 4 levels "Ape","New World Monkey",..: 4 4 4 4 4 2 2 2 2 2 ...
#>  $ species       : Factor w/ 29 levels "A palliata","Alouatta seniculus",..: 11 8 9 10 16 2 1 6 28 27 ...
#>  $ kcal.per.g    : num  0.49 0.51 0.46 0.48 0.6 0.47 0.56 0.89 0.91 0.92 ...
#>  $ perc.fat      : num  16.6 19.3 14.1 14.9 27.3 ...
#>  $ perc.protein  : num  15.4 16.9 16.9 13.2 19.5 ...
#>  $ perc.lactose  : num  68 63.8 69 71.9 53.2 ...
#>  $ mass          : num  1.95 2.09 2.51 1.62 2.19 5.25 5.37 2.51 0.71 0.68 ...
#>  $ neocortex.perc: num  55.2 NA NA NA NA ...
```
