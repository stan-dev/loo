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

