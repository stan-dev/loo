# tis throws correct errors and warnings

    Code
      psis(-LLarr, r_eff = r_eff_arr)
    Message
      Replacing NAs in `r_eff` with 1s
    Output
      Computed from 1000 by 32 log-weights matrix.
      MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.0]).
      
      All Pareto k estimates are good (k < 0.67).
      See help('pareto-k-diagnostic') for details.

# tis_loo and sis_loo are returned

    Code
      print(loo_tis)
    Output
      
      Computed from 1000 by 32 log-likelihood matrix using tis_loo .
      
               Estimate  SE
      elpd_loo    -83.6 4.3
      p_loo         3.3 1.2
      looic       167.2 8.6
      ------

---

    Code
      print(loo_sis)
    Output
      
      Computed from 1000 by 32 log-likelihood matrix using sis_loo .
      
               Estimate  SE
      elpd_loo    -83.6 4.3
      p_loo         3.3 1.2
      looic       167.2 8.6
      ------

