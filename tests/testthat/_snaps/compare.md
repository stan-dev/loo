# model_compare works with three loo_pred_measure models

    Code
      print(comp)
    Output
      Each measure compared against its own best model (elpd: B, mae: C, r2: B).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
           B     6
           C     4
           A     4
      
       model elpd_diff se_diff p_worse diag_diff
           B       0.0     0.0      NA          
           C     -22.4   129.6    0.57          
           A    -841.5   373.2    0.99          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

# print marks and explains flipped measure signs

    Code
      print(comp, measures = "all")
    Output
      Each measure compared against its own best model (mse: m2, r2: m2).
      PSIS-LOO unreliable for both models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m1     4
      
      -- mse (vs m2, sign flipped) --
       model mse_diff mse_se_diff
          m2      0.0         0.0
          m1   -199.6       460.2
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2   0.000      0.000
          m1  -0.098      0.223
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mse.

# print.compare.loo works for loo_pred_measure comparisons

    Code
      print(comp)
    Output
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
       model elpd_diff se_diff p_worse diag_diff
          m2       0.0     0.0      NA          
          m3     -22.4   129.6    0.57          
          m1    -841.5   373.2    0.99          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

---

    Code
      print(comp, measures = "all", digits = 2)
    Output
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
      -- elpd (vs m2) --
       model elpd_diff se_diff p_worse diag_diff
          m2      0.00    0.00      NA          
          m3    -22.44  129.62    0.57          
          m1   -841.47  373.24    0.99          
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff
          m3     0.00        0.00
          m2    -0.15        1.17
          m1    -6.35        3.07
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2    0.00       0.00
          m3   -0.07       0.16
          m1   -0.10       0.22
      
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
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
      -- elpd (vs m2) --
       model elpd_diff se_diff p_worse diag_diff
          m2       0.0     0.0      NA          
          m3     -22.4   129.6    0.57          
          m1    -841.5   373.2    0.99          
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff
          m3      0.0         0.0
          m2     -0.2         1.2
          m1     -6.3         3.1
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2     0.0        0.0
          m3    -0.1        0.2
          m1    -0.1        0.2
      
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
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2   0.000      0.000
          m3  -0.073      0.160
          m1  -0.098      0.223
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff
          m3      0.0         0.0
          m2     -0.2         1.2
          m1     -6.3         3.1
      
      All differences: 0 = best model, negative = worse.
      Signs flipped for loss measures: mae.

---

    Code
      print(comp, simplify = FALSE)
    Output
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
       model elpd_diff se_diff p_worse diag_diff    elpd se_elpd    p se_p  mae
          m2       0.0     0.0      NA           -2074.2   469.5 67.9 22.5 22.0
          m3     -22.4   129.6    0.57           -2096.7   438.5 89.9 41.0 21.9
          m1    -841.5   373.2    0.99           -2915.7   448.1 68.8 19.5 28.2
       se_mae    r2 se_r2
          3.4 0.144 0.221
          3.6 0.071 0.290
          3.2 0.046 0.040
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

---

    Code
      print(comp, measures = "all", simplify = FALSE)
    Output
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
      -- elpd (vs m2) --
       model elpd_diff se_diff p_worse diag_diff    elpd se_elpd    p se_p
          m2       0.0     0.0      NA           -2074.2   469.5 67.9 22.5
          m3     -22.4   129.6    0.57           -2096.7   438.5 89.9 41.0
          m1    -841.5   373.2    0.99           -2915.7   448.1 68.8 19.5
      
      -- mae (vs m3, sign flipped) --
       model mae_diff mae_se_diff  mae se_mae
          m3      0.0         0.0 21.9    3.6
          m2     -0.2         1.2 22.0    3.4
          m1     -6.3         3.1 28.2    3.2
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff    r2 se_r2
          m2   0.000      0.000 0.144 0.221
          m3  -0.073      0.160 0.071 0.290
          m1  -0.098      0.223 0.046 0.040
      
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
      Each measure compared against its own best model (elpd: m2, mae: m3, r2: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff    r2 se_r2
          m2   0.000      0.000 0.144 0.221
          m3  -0.073      0.160 0.071 0.290
          m1  -0.098      0.223 0.046 0.040

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

# print names only the ranking reference with more than four measures

    Code
      print(comp)
    Output
      Each measure compared against its own best model (elpd: m2, ...).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.5); measures may be biased.
       model bad_k
          m2     6
          m3     4
          m1     4
      
       model elpd_diff se_diff p_worse diag_diff
          m2       0.0     0.0      NA          
          m3     -22.4   129.6    0.57          
          m1    -841.5   373.2    0.99          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Use print(x, measures = "all") to see all measures.

