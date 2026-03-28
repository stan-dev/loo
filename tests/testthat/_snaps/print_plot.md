# mcse_loo extractor gives correct value

    WAoAAAACAAQFAAACAwAAAAAOAAAAAT+2J8YDcP5s

# print.loo supports kfold with pareto-k diagnostics - calibrated

    Code
      print(kfold1)
    Output
      
      Based on 10-fold cross-validation.
      
                 Estimate   SE
      elpd_kfold   -284.7 10.0
      p_kfold         2.4  0.6
      kfoldic       569.3 20.1
      ------
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

# print.loo supports kfold with pareto-k diagnostics - miscalibrated

    Code
      print(kfold1)
    Output
      
      Based on 10-fold cross-validation.
      
                 Estimate     SE
      elpd_kfold  -5521.0  713.1
      p_kfold       318.5   97.9
      kfoldic     11042.0 1426.3
      ------
      
      Pareto k diagnostic values:
                               Count Pct.    Min. ESS
      (-Inf, 0.7]   (good)     249   95.0%   <NA>    
         (0.7, 1]   (bad)        6    2.3%   <NA>    
         (1, Inf)   (very bad)   7    2.7%   <NA>    
      See help('pareto-k-diagnostic') for details.

# print.loo supports kfold without pareto-k diagnostics

    Code
      print(kfold1)
    Output
      
      Based on 10-fold cross-validation.
      
                 Estimate     SE
      elpd_kfold  -5521.0  713.1
      p_kfold       318.5   97.9
      kfoldic     11042.0 1426.3

