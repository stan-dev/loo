# Helper functions for K-fold cross-validation

These functions can be used to generate indexes for use with K-fold
cross-validation. See the **Details** section for explanations.

## Usage

``` r
kfold_split_random(K = 10, N = NULL)

kfold_split_stratified(K = 10, x = NULL)

kfold_split_grouped(K = 10, x = NULL)
```

## Arguments

- K:

  The number of folds to use.

- N:

  The number of observations in the data.

- x:

  A discrete variable of length `N` with at least `K` levels (unique
  values). Will be coerced to a
  [factor](https://rdrr.io/r/base/factor.html).

## Value

An integer vector of length `N` where each element is an index in `1:K`.

## Details

`kfold_split_random()` splits the data into `K` groups of equal size (or
roughly equal size).

For a categorical variable `x` `kfold_split_stratified()` splits the
observations into `K` groups ensuring that relative category frequencies
are approximately preserved.

For a grouping variable `x`, `kfold_split_grouped()` places all
observations in `x` from the same group/level together in the same fold.
The selection of which groups/levels go into which fold (relevant when
when there are more groups than folds) is randomized.

## Examples

``` r
ids <- kfold_split_random(K = 5, N = 20)
print(ids)
#>  [1] 3 5 2 4 3 3 5 1 4 1 2 4 2 2 5 1 1 4 3 5
table(ids)
#> ids
#> 1 2 3 4 5 
#> 4 4 4 4 4 


x <- sample(c(0, 1), size = 200, replace = TRUE, prob = c(0.05, 0.95))
table(x)
#> x
#>   0   1 
#>  12 188 
ids <- kfold_split_stratified(K = 5, x = x)
print(ids)
#>   [1] 3 4 1 4 5 5 2 2 2 4 2 1 5 5 1 3 5 1 4 5 1 3 1 3 1 2 2 3 1 2 2 3 4 4 2 1 3
#>  [38] 2 2 5 2 3 3 2 5 3 4 4 2 5 4 2 1 3 3 3 1 4 2 3 2 5 3 5 5 5 1 1 5 1 5 1 1 3
#>  [75] 5 2 5 5 5 1 4 3 4 1 4 3 2 4 2 5 3 3 4 3 2 5 2 2 5 5 4 4 2 3 4 1 2 2 5 1 3
#> [112] 3 1 1 5 1 1 3 2 2 3 3 5 4 4 3 2 5 4 4 5 3 5 4 5 4 3 1 5 2 2 2 3 1 4 1 2 3
#> [149] 5 4 4 3 1 3 5 1 1 5 5 4 1 3 4 5 1 1 2 4 2 4 4 1 2 1 3 1 4 4 1 5 1 2 4 5 1
#> [186] 4 4 4 4 2 3 4 1 3 5 2 2 5 3 3
table(ids, x)
#>    x
#> ids  0  1
#>   1  3 37
#>   2  3 37
#>   3  2 38
#>   4  2 38
#>   5  2 38

grp <- gl(n = 50, k = 15, labels = state.name)
length(grp)
#> [1] 750
head(table(grp))
#> grp
#>    Alabama     Alaska    Arizona   Arkansas California   Colorado 
#>         15         15         15         15         15         15 

ids_10 <- kfold_split_grouped(K = 10, x = grp)
(tab_10 <- table(grp, ids_10))
#>                 ids_10
#> grp               1  2  3  4  5  6  7  8  9 10
#>   Alabama         0  0 15  0  0  0  0  0  0  0
#>   Alaska          0  0  0  0  0  0  0 15  0  0
#>   Arizona         0  0  0 15  0  0  0  0  0  0
#>   Arkansas        0  0  0  0  0 15  0  0  0  0
#>   California      0  0  0  0  0 15  0  0  0  0
#>   Colorado        0  0  0  0  0  0 15  0  0  0
#>   Connecticut     0 15  0  0  0  0  0  0  0  0
#>   Delaware        0  0  0  0 15  0  0  0  0  0
#>   Florida         0  0  0  0  0  0  0 15  0  0
#>   Georgia         0  0  0  0  0  0  0  0 15  0
#>   Hawaii         15  0  0  0  0  0  0  0  0  0
#>   Idaho           0  0  0 15  0  0  0  0  0  0
#>   Illinois        0  0  0  0  0  0  0 15  0  0
#>   Indiana         0 15  0  0  0  0  0  0  0  0
#>   Iowa            0  0  0  0  0  0  0  0  0 15
#>   Kansas         15  0  0  0  0  0  0  0  0  0
#>   Kentucky        0  0  0  0  0  0  0  0 15  0
#>   Louisiana       0  0  0  0  0  0  0  0  0 15
#>   Maine           0  0  0  0 15  0  0  0  0  0
#>   Maryland       15  0  0  0  0  0  0  0  0  0
#>   Massachusetts   0 15  0  0  0  0  0  0  0  0
#>   Michigan        0  0  0 15  0  0  0  0  0  0
#>   Minnesota       0  0  0 15  0  0  0  0  0  0
#>   Mississippi     0 15  0  0  0  0  0  0  0  0
#>   Missouri        0  0  0  0 15  0  0  0  0  0
#>   Montana         0  0  0  0  0  0 15  0  0  0
#>   Nebraska        0  0 15  0  0  0  0  0  0  0
#>   Nevada          0  0  0  0  0  0  0  0  0 15
#>   New Hampshire   0  0  0  0  0  0  0  0  0 15
#>   New Jersey      0  0  0  0  0  0 15  0  0  0
#>   New Mexico      0  0  0  0  0  0  0  0 15  0
#>   New York        0  0 15  0  0  0  0  0  0  0
#>   North Carolina  0  0  0  0  0 15  0  0  0  0
#>   North Dakota    0  0  0  0  0  0  0 15  0  0
#>   Ohio            0 15  0  0  0  0  0  0  0  0
#>   Oklahoma        0  0 15  0  0  0  0  0  0  0
#>   Oregon          0  0  0  0  0  0 15  0  0  0
#>   Pennsylvania    0  0  0  0  0  0  0 15  0  0
#>   Rhode Island    0  0  0  0  0  0  0  0 15  0
#>   South Carolina  0  0  0  0  0  0  0  0  0 15
#>   South Dakota    0  0  0  0  0  0  0  0 15  0
#>   Tennessee       0  0  0  0 15  0  0  0  0  0
#>   Texas           0  0  0  0  0  0 15  0  0  0
#>   Utah           15  0  0  0  0  0  0  0  0  0
#>   Vermont         0  0  0  0 15  0  0  0  0  0
#>   Virginia        0  0  0  0  0 15  0  0  0  0
#>   Washington      0  0  0 15  0  0  0  0  0  0
#>   West Virginia   0  0  0  0  0 15  0  0  0  0
#>   Wisconsin       0  0 15  0  0  0  0  0  0  0
#>   Wyoming        15  0  0  0  0  0  0  0  0  0
colSums(tab_10)
#>  1  2  3  4  5  6  7  8  9 10 
#> 75 75 75 75 75 75 75 75 75 75 

ids_9 <- kfold_split_grouped(K = 9, x = grp)
(tab_9 <- table(grp, ids_9))
#>                 ids_9
#> grp               1  2  3  4  5  6  7  8  9
#>   Alabama         0  0  0  0  0  0 15  0  0
#>   Alaska          0 15  0  0  0  0  0  0  0
#>   Arizona         0  0  0  0  0 15  0  0  0
#>   Arkansas        0  0  0  0 15  0  0  0  0
#>   California      0 15  0  0  0  0  0  0  0
#>   Colorado        0  0  0  0  0  0  0 15  0
#>   Connecticut     0  0  0  0 15  0  0  0  0
#>   Delaware        0  0  0  0  0  0  0  0 15
#>   Florida         0  0  0  0  0  0 15  0  0
#>   Georgia         0  0  0  0  0  0  0  0 15
#>   Hawaii          0 15  0  0  0  0  0  0  0
#>   Idaho           0  0  0  0 15  0  0  0  0
#>   Illinois       15  0  0  0  0  0  0  0  0
#>   Indiana        15  0  0  0  0  0  0  0  0
#>   Iowa            0  0  0  0  0  0  0 15  0
#>   Kansas          0  0  0 15  0  0  0  0  0
#>   Kentucky        0  0 15  0  0  0  0  0  0
#>   Louisiana      15  0  0  0  0  0  0  0  0
#>   Maine           0  0  0 15  0  0  0  0  0
#>   Maryland        0  0  0  0  0  0 15  0  0
#>   Massachusetts   0  0  0  0  0  0  0  0 15
#>   Michigan        0  0  0  0  0  0  0 15  0
#>   Minnesota       0  0  0  0  0 15  0  0  0
#>   Mississippi     0  0 15  0  0  0  0  0  0
#>   Missouri       15  0  0  0  0  0  0  0  0
#>   Montana         0  0  0  0  0  0  0  0 15
#>   Nebraska        0 15  0  0  0  0  0  0  0
#>   Nevada          0  0  0 15  0  0  0  0  0
#>   New Hampshire   0  0  0 15  0  0  0  0  0
#>   New Jersey      0  0  0  0  0  0 15  0  0
#>   New Mexico      0 15  0  0  0  0  0  0  0
#>   New York        0  0 15  0  0  0  0  0  0
#>   North Carolina  0  0  0 15  0  0  0  0  0
#>   North Dakota    0  0  0  0  0  0  0 15  0
#>   Ohio            0  0  0  0  0  0 15  0  0
#>   Oklahoma        0  0  0  0  0 15  0  0  0
#>   Oregon          0  0  0  0  0 15  0  0  0
#>   Pennsylvania    0  0 15  0  0  0  0  0  0
#>   Rhode Island    0  0  0  0  0  0  0 15  0
#>   South Carolina  0  0  0  0  0 15  0  0  0
#>   South Dakota    0  0  0  0 15  0  0  0  0
#>   Tennessee       0  0  0  0  0  0  0  0 15
#>   Texas           0  0  0 15  0  0  0  0  0
#>   Utah            0  0 15  0  0  0  0  0  0
#>   Vermont        15  0  0  0  0  0  0  0  0
#>   Virginia        0  0  0  0 15  0  0  0  0
#>   Washington     15  0  0  0  0  0  0  0  0
#>   West Virginia   0 15  0  0  0  0  0  0  0
#>   Wisconsin       0  0  0  0 15  0  0  0  0
#>   Wyoming         0  0 15  0  0  0  0  0  0
colSums(tab_9)
#>  1  2  3  4  5  6  7  8  9 
#> 90 90 90 90 90 75 75 75 75 
```
