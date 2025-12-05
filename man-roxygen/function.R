#' @describeIn <%= fn %>
#'  A function `f()` that takes arguments `data_i` and `draws` and returns a
#'  vector containing the log-likelihood for a single observation `i` evaluated
#'  at each posterior draw. The function should be written such that, for each
#'  observation `i` in `1:N`, evaluating
#'
#'     f(data_i = data[i,, drop=FALSE], draws = draws)
#'
#'  results in a vector of length `S` (size of posterior sample). The
#'  log-likelihood function can also have additional arguments but `data_i` and
#'  `draws` are required.
#'
#'  If using the function method then the arguments `data` and `draws` must also
#'  be specified in the call to `loo()`:
#'  * `data`: A data frame or matrix containing the data (e.g.
#'    observed outcome and predictors) needed to compute the pointwise
#'    log-likelihood. For each observation `i`, the `i`th row of
#'    `data` will be passed to the `data_i` argument of the
#'    log-likelihood function.
#'  * `draws`: An object containing the posterior draws for any
#'    parameters needed to compute the pointwise log-likelihood. Unlike
#'    `data`, which is indexed by observation, for each observation the
#'    entire object `draws` will be passed to the `draws` argument of
#'    the log-likelihood function.
#'  * The `...` can be used if your log-likelihood function takes additional
#'    arguments. These arguments are used like the `draws` argument in that they
#'    are recycled for each observation.
#'
