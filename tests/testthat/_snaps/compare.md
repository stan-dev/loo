# model_compare works with three loo_pred_measure models

    Code
      print(comp)
    Output
      Each measure compared against its own best model (elpd: B, r2: B, mae: C).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
           B     2
           C     1
           A     1
      
       model elpd_diff se_diff p_worse diag_diff
           B       0.0     0.0      NA          
           C     -25.5   129.1    0.58          
           A    -850.3   372.3    0.99          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

# print marks and explains flipped measure signs

    Code
      print(comp, measures = "all")
    Output
      Each measure compared against its own best model (r2: m2, mse: m2).
      PSIS-LOO unreliable for both models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m1     1
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2   0.000      0.000
          m1  -0.105      0.217
      
      -- mse (vs m2, sign flipped) --
       model mse_diff mse_se_diff
          m2      0.0         0.0
          m1   -212.5       448.0
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mse.

# print.compare.loo works for loo_pred_measure comparisons

    Code
      print(comp)
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
       model elpd_diff se_diff p_worse diag_diff
          m2       0.0     0.0      NA          
          m3     -25.5   129.1    0.58          
          m1    -850.3   372.3    0.99          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

---

    Code
      print(comp, measures = "all", digits = 2)
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
      -- elpd (vs m2) --
       model elpd_diff se_diff p_worse diag_diff
          m2      0.00    0.00      NA          
          m3    -25.47  129.10    0.58          
          m1   -850.29  372.31    0.99          
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2    0.00       0.00
          m3   -0.09       0.18
          m1   -0.10       0.22
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff
          m3     0.00        0.00
          m2    -0.07        1.24
          m1    -6.34        3.08
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mae.
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp, measures = "all", digits = c(r2 = 1))
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
      -- elpd (vs m2) --
       model elpd_diff se_diff p_worse diag_diff
          m2       0.0     0.0      NA          
          m3     -25.5   129.1    0.58          
          m1    -850.3   372.3    0.99          
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2     0.0        0.0
          m3    -0.1        0.2
          m1    -0.1        0.2
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff
          m3      0.0         0.0
          m2     -0.1         1.2
          m1     -6.3         3.1
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mae.
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp, measures = c("r2", "mae"))
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2   0.000      0.000
          m3  -0.093      0.178
          m1  -0.105      0.217
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff
          m3      0.0         0.0
          m2     -0.1         1.2
          m1     -6.3         3.1
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mae.

---

    Code
      print(comp, simplify = FALSE)
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
       model elpd_diff se_diff p_worse diag_diff    elpd se_elpd    p se_p    r2
          m2       0.0     0.0      NA           -2071.4   468.9 61.6 20.7 0.153
          m3     -25.5   129.1    0.58           -2096.8   438.7 96.4 46.2 0.060
          m1    -850.3   372.3    0.99           -2921.7   449.8 75.8 21.1 0.048
       se_r2  mae se_mae
       0.214 22.0    3.4
       0.300 21.9    3.6
       0.040 28.2    3.2
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

---

    Code
      print(comp, measures = "all", simplify = FALSE)
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
      -- elpd (vs m2) --
       model elpd_diff se_diff p_worse diag_diff    elpd se_elpd    p se_p
          m2       0.0     0.0      NA           -2071.4   468.9 61.6 20.7
          m3     -25.5   129.1    0.58           -2096.8   438.7 96.4 46.2
          m1    -850.3   372.3    0.99           -2921.7   449.8 75.8 21.1
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff    r2 se_r2
          m2   0.000      0.000 0.153 0.214
          m3  -0.093      0.178 0.060 0.300
          m1  -0.105      0.217 0.048 0.040
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff  mae se_mae
          m3      0.0         0.0 21.9    3.6
          m2     -0.1         1.2 22.0    3.4
          m1     -6.3         3.1 28.2    3.2
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mae.
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp, measures = "r2", simplify = FALSE)
    Output
      Each measure compared against its own best model (elpd: m2, r2: m2, mae: m3).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m2     2
          m3     1
          m1     1
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff    r2 se_r2
          m2   0.000      0.000 0.153 0.214
          m3  -0.093      0.178 0.060 0.300
          m1  -0.105      0.217 0.048 0.040

# model_compare returns expected results (2 models)

    WAoAAAACAAQGAQACAwAAAAMTAAAADAAAABAAAAACAAQACQAAAAZtb2RlbDEABAAJAAAABm1v
    ZGVsMgAAAA4AAAACAAAAAAAAAAAAAAAAAAAAAAAAAA4AAAACAAAAAAAAAAAAAAAAAAAAAAAA
    AA4AAAACf/AAAAAAB6J/8AAAAAAHogAAABAAAAACAAQACQAAAAAABAAJAAAAAAAAABAAAAAC
    AAQACQAAAAAABAAJAAAAAAAAAA4AAAACwFTh8N3JQljAVOHw3clCWAAAAA4AAAACQBEIPbMR
    cF9AEQg9sxFwXwAAAA4AAAACQAoowGHVuVNACijAYdW5UwAAAA4AAAACP/H9Zexy814/8f1l
    7HLzXgAAAA4AAAACQGTh8N3JQlhAZOHw3clCWAAAAA4AAAACQCEIPbMRcF9AIQg9sxFwXwAA
    BAIAAAABAAQACQAAAAVuYW1lcwAAABAAAAAMAAQACQAAAAVtb2RlbAAEAAkAAAAJZWxwZF9k
    aWZmAAQACQAAAAdzZV9kaWZmAAQACQAAAAdwX3dvcnNlAAQACQAAAAlkaWFnX2RpZmYABAAJ
    AAAACWRpYWdfZWxwZAAEAAkAAAAJZWxwZF93YWljAAQACQAAAAxzZV9lbHBkX3dhaWMABAAJ
    AAAABnBfd2FpYwAEAAkAAAAJc2VfcF93YWljAAQACQAAAAR3YWljAAQACQAAAAdzZV93YWlj
    AAAEAgAAAAEABAAJAAAABWNsYXNzAAAAEAAAAAIABAAJAAAAC2NvbXBhcmUubG9vAAQACQAA
    AApkYXRhLmZyYW1lAAAEAgAAAAEABAAJAAAACXJvdy5uYW1lcwAAAA0AAAACgAAAAP////4A
    AAQCAAAAAQAEAAkAAAARY29tcGFyZV9yZWZlcmVuY2UAAAIQAAAAAQAEAAkAAAAGbW9kZWwx
    AAAEAgAAAf8AAAAQAAAAAQAEAAkAAAAEZWxwZAAAAP4AAAD+

---

    Code
      print(comp1)
    Output
        model elpd_diff se_diff p_worse diag_diff diag_elpd
       model1       0.0     0.0      NA                    
       model2       0.0     0.0      NA                    

---

    WAoAAAACAAQGAQACAwAAAAMTAAAADAAAABAAAAACAAQACQAAAAZtb2RlbDEABAAJAAAABm1v
    ZGVsMgAAAA4AAAACAAAAAAAAAADAEDpTX5xF7gAAAA4AAAACAAAAAAAAAAA/tmpHtC8TAQAA
    AA4AAAACf/AAAAAAB6I/8AAAAAAAAAAAABAAAAACAAQACQAAAAAABAAJAAAAB04gPCAxMDAA
    AAAQAAAAAgAEAAkAAAAAAAQACQAAAAAAAAAOAAAAAsBU4fDdyUJYwFXllhPDBrkAAAAOAAAA
    AkARCD2zEXBfQBEalRIN2T8AAAAOAAAAAkAKKMBh1blTQCZnlesA0IoAAAAOAAAAAj/x/WXs
    cvNeP/GbYJxtZ8cAAAAOAAAAAkBk4fDdyUJYQGXllhPDBrkAAAAOAAAAAkAhCD2zEXBfQCEa
    lRIN2T8AAAQCAAAAAQAEAAkAAAAFbmFtZXMAAAAQAAAADAAEAAkAAAAFbW9kZWwABAAJAAAA
    CWVscGRfZGlmZgAEAAkAAAAHc2VfZGlmZgAEAAkAAAAHcF93b3JzZQAEAAkAAAAJZGlhZ19k
    aWZmAAQACQAAAAlkaWFnX2VscGQABAAJAAAACWVscGRfd2FpYwAEAAkAAAAMc2VfZWxwZF93
    YWljAAQACQAAAAZwX3dhaWMABAAJAAAACXNlX3Bfd2FpYwAEAAkAAAAEd2FpYwAEAAkAAAAH
    c2Vfd2FpYwAABAIAAAABAAQACQAAAAVjbGFzcwAAABAAAAACAAQACQAAAAtjb21wYXJlLmxv
    bwAEAAkAAAAKZGF0YS5mcmFtZQAABAIAAAABAAQACQAAAAlyb3cubmFtZXMAAAANAAAAAoAA
    AAD////+AAAEAgAAAAEABAAJAAAAEWNvbXBhcmVfcmVmZXJlbmNlAAACEAAAAAEABAAJAAAA
    Bm1vZGVsMQAABAIAAAH/AAAAEAAAAAEABAAJAAAABGVscGQAAAD+AAAA/g==

---

    Code
      print(comp2)
    Output
        model elpd_diff se_diff p_worse diag_diff diag_elpd
       model1       0.0     0.0      NA                    
       model2      -4.1     0.1    1.00   N < 100          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp2, p_worse = FALSE)
    Output
        model elpd_diff se_diff
       model1       0.0     0.0
       model2      -4.1     0.1

---

    Code
      print(comp2, simplify = FALSE)
    Output
        model elpd_diff se_diff p_worse diag_diff diag_elpd elpd_waic se_elpd_waic
       model1       0.0     0.0      NA                         -83.5          4.3
       model2      -4.1     0.1    1.00   N < 100               -87.6          4.3
       p_waic se_p_waic  waic se_waic
          3.3       1.1 167.1     8.5
         11.2       1.1 175.2     8.6
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp2, simplify = FALSE, p_worse = FALSE)
    Output
        model elpd_diff se_diff elpd_waic se_elpd_waic p_waic se_p_waic  waic se_waic
       model1       0.0     0.0     -83.5          4.3    3.3       1.1 167.1     8.5
       model2      -4.1     0.1     -87.6          4.3   11.2       1.1 175.2     8.6

# model_compare returns expected result (3 models)

    WAoAAAACAAQGAQACAwAAAAMTAAAADAAAABAAAAADAAQACQAAAAZtb2RlbDEABAAJAAAABm1v
    ZGVsMgAEAAkAAAAGbW9kZWwzAAAADgAAAAMAAAAAAAAAAMAQOlNfnEXuwDANypG2BBgAAAAO
    AAAAAwAAAAAAAAAAP7ZqR7QvEwE/y6/t4TTtYAAAAA4AAAADf/AAAAAAB6I/8AAAAAAAAD/w
    AAAAAAAAAAAAEAAAAAMABAAJAAAAAAAEAAkAAAAHTiA8IDEwMAAEAAkAAAAHTiA8IDEwMAAA
    ABAAAAADAAQACQAAAAAABAAJAAAAAAAEAAkAAAAAAAAADgAAAAPAVOHw3clCWMBV5ZYTwwa5
    wFjlY4I2w2IAAAAOAAAAA0ARCD2zEXBfQBEalRIN2T9AEPIF3GigEwAAAA4AAAADQAoowGHV
    uVNAJmeV6wDQikBByNhSGt0KAAAADgAAAAM/8f1l7HLzXj/xm2CcbWfHP/GA0JJnyV8AAAAO
    AAAAA0Bk4fDdyUJYQGXllhPDBrlAaOVjgjbDYgAAAA4AAAADQCEIPbMRcF9AIRqVEg3ZP0Ag
    8gXcaKATAAAEAgAAAAEABAAJAAAABW5hbWVzAAAAEAAAAAwABAAJAAAABW1vZGVsAAQACQAA
    AAllbHBkX2RpZmYABAAJAAAAB3NlX2RpZmYABAAJAAAAB3Bfd29yc2UABAAJAAAACWRpYWdf
    ZGlmZgAEAAkAAAAJZGlhZ19lbHBkAAQACQAAAAllbHBkX3dhaWMABAAJAAAADHNlX2VscGRf
    d2FpYwAEAAkAAAAGcF93YWljAAQACQAAAAlzZV9wX3dhaWMABAAJAAAABHdhaWMABAAJAAAA
    B3NlX3dhaWMAAAQCAAAAAQAEAAkAAAAFY2xhc3MAAAAQAAAAAgAEAAkAAAALY29tcGFyZS5s
    b28ABAAJAAAACmRhdGEuZnJhbWUAAAQCAAAAAQAEAAkAAAAJcm93Lm5hbWVzAAAADQAAAAKA
    AAAA/////QAABAIAAAABAAQACQAAABFjb21wYXJlX3JlZmVyZW5jZQAAAhAAAAABAAQACQAA
    AAZtb2RlbDEAAAQCAAAB/wAAABAAAAABAAQACQAAAARlbHBkAAAA/gAAAP4=

---

    Code
      print(comp1)
    Output
        model elpd_diff se_diff p_worse diag_diff diag_elpd
       model1       0.0     0.0      NA                    
       model2      -4.1     0.1    1.00   N < 100          
       model3     -16.1     0.2    1.00   N < 100          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

# model_compare with simplify=FALSE returns expected result

    Code
      print(comp, simplify = FALSE)
    Output
        model elpd_diff se_diff p_worse diag_diff diag_elpd elpd_loo se_elpd_loo
       model3       0.0     0.0      NA                        -19.6         4.3
       model2     -32.0     0.0    1.00   N < 100              -51.6         4.3
       model1     -64.0     0.0    1.00   N < 100              -83.6         4.3
       p_loo se_p_loo looic se_looic
         3.3      1.2  39.2      8.6
         3.3      1.2 103.2      8.6
         3.3      1.2 167.2      8.6
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

# compare returns expected result (3 models)

    WAoAAAACAAQFAAACAwAAAAMOAAAAGAAAAAAAAAAAwBA6U1+cRe7AMA3KkbYEGAAAAAAAAAAA
    P7ZqR7QvEwE/y6/t4TTtXsBU4fDdyUJYwFXllhPDBrnAWOVjgjbDYkARCD2zEXBfQBEalRIN
    2T9AEPIF3GigE0AKKMBh1blTQCZnlesA0IpAQcjYUhrdCj/x/WXscvNeP/GbYJxtZ8c/8YDQ
    kmfJX0Bk4fDdyUJYQGXllhPDBrlAaOVjgjbDYkAhCD2zEXBfQCEalRIN2T9AIPIF3GigEwAA
    BAIAAAABAAQACQAAAANkaW0AAAANAAAAAgAAAAMAAAAIAAAEAgAAAAEABAAJAAAACGRpbW5h
    bWVzAAAAEwAAAAIAAAAQAAAAAwAEAAkAAAACdzEABAAJAAAAAncyAAQACQAAAAJ3MwAAABAA
    AAAIAAQACQAAAAllbHBkX2RpZmYABAAJAAAAB3NlX2RpZmYABAAJAAAACWVscGRfd2FpYwAE
    AAkAAAAMc2VfZWxwZF93YWljAAQACQAAAAZwX3dhaWMABAAJAAAACXNlX3Bfd2FpYwAEAAkA
    AAAEd2FpYwAEAAkAAAAHc2Vfd2FpYwAABAIAAAABAAQACQAAAAVjbGFzcwAAABAAAAAEAAQA
    CQAAAAtjb21wYXJlLmxvbwAEAAkAAAAGbWF0cml4AAQACQAAAAVhcnJheQAEAAkAAAAPb2xk
    X2NvbXBhcmUubG9vAAAA/g==

