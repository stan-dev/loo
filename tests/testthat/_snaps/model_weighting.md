# loo_model_weights (stacking and pseudo-BMA) gives expected result

    WAoAAAACAAQFAgACAwAAAAAOAAAAAz/KEXnf1DM9P+l7oXTYyi4+YzIpAwAAAA==

---

    Code
      print(w1)
    Output
      Method: stacking
      ------
             weight
      model1 0.204 
      model2 0.796 
      model3 0.000 

---

    WAoAAAACAAQFAgACAwAAAAAOAAAAAz+xA6UGtqDDP+3eFS5zKzY/J2MLYAsc4A==

---

    Code
      print(w2)
    Output
      Method: pseudo-BMA+ with Bayesian bootstrap
      ------
             weight
      model1 0.066 
      model2 0.933 
      model3 0.000 

---

    Code
      print(w3)
    Output
      Method: pseudo-BMA
      ------
             weight
      model1 0.000 
      model2 1.000 
      model3 0.000 

