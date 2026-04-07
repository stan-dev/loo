# Test the vignette

    Code
      print(looss_1)
    Output
      
      Computed from 4000 by 100 subsampled log-likelihood
      values from 3020 total observations.
      
               Estimate   SE subsampling SE
      elpd_loo  -1968.5 15.6            0.3
      p_loo         3.1  0.1            0.4
      looic      3936.9 31.2            0.6
      ------
      MCSE of elpd_loo is 0.0.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(looss_1b)
    Output
      
      Computed from 4000 by 200 subsampled log-likelihood
      values from 3020 total observations.
      
               Estimate   SE subsampling SE
      elpd_loo  -1968.3 15.6            0.2
      p_loo         3.2  0.1            0.4
      looic      3936.7 31.2            0.5
      ------
      MCSE of elpd_loo is 0.0.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(aploo_1)
    Output
      
      Computed from 2000 by 3020 log-likelihood matrix.
      
               Estimate   SE
      elpd_loo  -1968.4 15.6
      p_loo         3.2  0.2
      looic      3936.8 31.2
      ------
      Posterior approximation correction used.
      MCSE of elpd_loo is 0.0.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(looapss_1)
    Output
      
      Computed from 2000 by 100 subsampled log-likelihood
      values from 3020 total observations.
      
               Estimate   SE subsampling SE
      elpd_loo  -1968.2 15.6            0.4
      p_loo         2.9  0.1            0.5
      looic      3936.4 31.1            0.8
      ------
      Posterior approximation correction used.
      MCSE of elpd_loo is 0.0.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(looss_2)
    Output
      
      Computed from 4000 by 100 subsampled log-likelihood
      values from 3020 total observations.
      
               Estimate   SE subsampling SE
      elpd_loo  -1952.0 16.2            0.2
      p_loo         2.6  0.1            0.3
      looic      3903.9 32.4            0.4
      ------
      MCSE of elpd_loo is 0.0.
      MCSE and ESS estimates assume independent draws (r_eff=1).
      
      All Pareto k estimates are good (k < 0.7).
      See help('pareto-k-diagnostic') for details.

---

    Code
      print(comp)
    Output
             elpd_diff se_diff subsampling_se_diff
      model2  0.0       0.0     0.0               
      model1 16.5      22.5     0.4               

---

    Code
      print(comp)
    Output
             elpd_diff se_diff subsampling_se_diff
      model2  0.0       0.0     0.0               
      model1 16.1       4.4     0.1               

---

    Code
      print(comp2)
    Output
             elpd_diff se_diff subsampling_se_diff
      model2  0.0       0.0     0.0               
      model1 16.3       4.4     0.1               

---

    Code
      print(comp3)
    Output
             elpd_diff se_diff subsampling_se_diff
      model2  0.0       0.0     0.0               
      model1 16.5       4.4     0.3               

