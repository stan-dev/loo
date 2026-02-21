# psis throws correct errors and warnings

    Code
      psis(-LLarr, r_eff = r_eff_arr)
    Message
      Replacing NAs in `r_eff` with 1s
    Output
      Computed from 1000 by 32 log-weights matrix.
      MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.0]).
      
      All Pareto k estimates are good (k < 0.67).
      See help('pareto-k-diagnostic') for details.

---

    Code
      psis(-LLarr[1:5, , ])
    Condition
      Warning:
      Not enough tail samples to fit the generalized Pareto distribution in some or all columns of matrix of log importance ratios. Skipping the following columns: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ... [22 more not printed].
      Warning:
      Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
    Output
      Computed from 10 by 32 log-weights matrix.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      Pareto k diagnostic values:
                             Count Pct.    Min. ESS
      (-Inf, 0]   (good)      0      0.0%  <NA>    
         (0, 1]   (bad)       0      0.0%  <NA>    
       (1, Inf)   (very bad) 32    100.0%  <NA>    
      See help('pareto-k-diagnostic') for details.
