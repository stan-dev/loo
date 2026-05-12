# print.waic output is ok

    Code
      print(waic1)
    Output
      
      Computed from 1000 by 32 log-likelihood matrix.
      
                Estimate  SE
      elpd_waic    -83.5 4.3
      p_waic         3.3 1.1
      waic         167.1 8.5
      
      3 (9.4%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

# print.psis_loo and print.psis output ok

    Code
      print(psis1)
    Output
      Computed from 1000 by 32 log-weights matrix.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.67).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(loo1)
    Output
      
      Computed from 1000 by 32 log-likelihood matrix.
      
               Estimate  SE
      elpd_loo    -83.6 4.3
      p_loo         3.3 1.2
      looic       167.2 8.6
      ------
      MCSE of elpd_loo is 0.1.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.67).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(loo1_r_eff)
    Output
      
      Computed from 1000 by 32 log-likelihood matrix.
      
               Estimate  SE
      elpd_loo    -83.6 4.3
      p_loo         3.3 1.2
      looic       167.2 8.6
      ------
      MCSE of elpd_loo is 0.1.
      MCSE and ESS estimates assume MCMC draws (r_eff in [0.6, 1.0]).
      
      All Pareto k estimates are good (k < 0.67).
      See help('pareto-k-diagnostic') for details.

# mcse_loo extractor gives correct value

    WAoAAAACAAQFAAACAwAAAAAOAAAAAT+2J8YDcP5s

# print.loo supports kfold with pareto-k diagnostics - calibrated

    Code
      print(kfold1)
    Output
      
      Based on 10-fold cross-validation.
      
                 Estimate   SE
      elpd_kfold   -285.0  9.2
      p_kfold         2.5  0.6
      kfoldic       570.0 18.4
      ------
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

# print.loo supports kfold with pareto-k diagnostics - miscalibrated

    Code
      print(kfold1)
    Output
      
      Based on 10-fold cross-validation.
      
                 Estimate     SE
      elpd_kfold  -5556.6  701.0
      p_kfold       358.2  108.5
      kfoldic     11113.1 1401.9
      ------
      
      Pareto k diagnostic values:
                               Count Pct.    Min. ESS
      (-Inf, 0.7]   (good)     245   93.5%   24      
         (0.7, 1]   (bad)        8    3.1%   <NA>    
         (1, Inf)   (very bad)   9    3.4%   <NA>    
      See help('pareto-k-diagnostic') for details.

# print.loo supports kfold without pareto-k diagnostics

    Code
      print(kfold1)
    Output
      
      Based on 10-fold cross-validation.
      
                 Estimate     SE
      elpd_kfold  -5556.6  701.0
      p_kfold       358.2  108.5
      kfoldic     11113.1 1401.9

