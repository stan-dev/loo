# loo_compare works with three loo_pred_measure models

    Code
      print(comp)
    Output
      Models ranked by mae (reference: B).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
           A    23
           C    19
           B     9
      
       model mae_diff mae_se_diff
           B      0.0         0.0
           C      0.0         0.7
           A     -7.6         1.5
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Other measures compared: elpd, r2. Use print(x, measures = "all").

# loo_compare informs when measure signs are converted

    Code
      comp <- loo_compare(pm1, pm2)
    Message
      For model comparison, differences for mse is
      reported on a utility scale (higher is better).

# print.compare.loo works for loo_pred_measure comparisons

    Code
      print(comp)
    Output
      Each measure compared against its own best model (elpd: m3, r2: m2, mae: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m1    23
          m3    19
          m2     9
      
       model elpd_diff se_diff p_worse diag_diff
          m3       0.0     0.0      NA          
          m2    -205.5   221.8    0.82          
          m1   -2739.1   630.4    1.00          
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Other measures compared: r2, mae. Use print(x, measures = "all").

---

    Code
      print(comp, measures = "all", digits = 2)
    Output
      Each measure compared against its own best model (elpd: m3, r2: m2, mae: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m1    23
          m3    19
          m2     9
      
      -- elpd (vs m3) --
       model elpd_diff se_diff p_worse diag_diff
          m3      0.00    0.00      NA          
          m2   -205.45  221.82    0.82          
          m1  -2739.08  630.38    1.00          
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2    0.00       0.00
          m3   -0.03       0.08
          m1   -0.22       0.09
      
      -- mae (vs m2) --
       model mae_diff mae_se_diff
          m2     0.00        0.00
          m3    -0.02        0.72
          m1    -7.60        1.54
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp, measures = c("r2", "mae"))
    Output
      Each measure compared against its own best model (elpd: m3, r2: m2, mae: m2).
      PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
       model bad_k
          m1    23
          m3    19
          m2     9
      
      -- r2 (vs m2) --
       model r2_diff r2_se_diff
          m2     0.0        0.0
          m3     0.0        0.1
          m1    -0.2        0.1
      
      -- mae (vs m2) --
       model mae_diff mae_se_diff
          m2      0.0         0.0
          m3      0.0         0.7
          m1     -7.6         1.5
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.

---

    Code
      print(comp_mae)
    Output
      Models ranked by mae (reference: m2).
      PSIS-LOO unreliable for both models (k_psis > 0.62); measures may be biased.
       model bad_k
          m1    23
          m2     9
      
       model mae_diff mae_se_diff
          m2      0.0         0.0
          m1     -7.6         1.5
    Message
      
      Diagnostic flags present.
      See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
      or https://mc-stan.org/loo/reference/loo-glossary.html.
      
      Other measures compared: elpd, r2. Use print(x, measures = "all").

# loo_compare returns expected results (2 models)

    WAoAAAACAAQEAgACAwAAAAMTAAAADAAAABAAAAACAAQACQAAAAZtb2RlbDEABAAJAAAABm1v
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
    AAD+

---

    Code
      print(comp1)
    Output
        model elpd_diff se_diff p_worse diag_diff diag_elpd
       model1       0.0     0.0      NA                    
       model2       0.0     0.0      NA                    

---

    WAoAAAACAAQEAgACAwAAAAMTAAAADAAAABAAAAACAAQACQAAAAZtb2RlbDEABAAJAAAABm1v
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
    AAD////+AAAA/g==

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

# loo_compare returns expected result (3 models)

    WAoAAAACAAQEAgACAwAAAAMTAAAADAAAABAAAAADAAQACQAAAAZtb2RlbDEABAAJAAAABm1v
    ZGVsMgAEAAkAAAAGbW9kZWwzAAAADgAAAAMAAAAAAAAAAMAQOlNfnEXuwDANypG2BBgAAAAO
    AAAAAwAAAAAAAAAAP7ZqR7QvEwE/y6/t4TTtXgAAAA4AAAADf/AAAAAAB6I/8AAAAAAAAD/w
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
    AAAA/////QAAAP4=

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

