loo_compare accepts as input the outputs of:
- loo()
- kfold()
- waic()

They all produce outputs that inherit the "loo" class.

Deprectation
============
- remove "compare()" which is already deprecated for a while.
- deprecate "loo_compare()"
- introduce "model_compare()"

New implementation
==================
model_compare(type = "kfold"/"loo"/"test"/"insample") -> consistent with <type>_pred_measure



Trivia
======

