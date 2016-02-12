## Submission summary

This update adds additional unit tests, removes model weights from the output
given by the compare function, and makes minor improvements to 
the package documentation and messages issued by the waic function.

## Test environments
* local OS X install, R 3.2.3
* win-builder (release and devel)
* ubuntu 12.04 (on travis-ci)

## R CMD check results
No ERRORs, WARNINGs, or NOTEs

## Downstream dependencies
One very minor test fails in the brms package. I have discussed this with the 
author of brms and he will make the necessary fix, but the functionality of brms
is not affected by this loo update.
