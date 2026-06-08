# loo_compare_subsample

    Code
      lcss <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2, lss3))
    Condition
      Warning:
      Different subsamples in 'model3' and 'model2'. Naive diff SE is used.
      Warning:
      Different subsamples in 'model3' and 'model1'. Naive diff SE is used.

---

    Code
      lcssapi <- loo_compare(lss1, lss2, lss3)
    Condition
      Warning:
      Different subsamples in 'model3' and 'model2'. Naive diff SE is used.
      Warning:
      Different subsamples in 'model3' and 'model1'. Naive diff SE is used.

