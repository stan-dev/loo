## Test environments
* local OS X install, R 3.4.4
* ubuntu 14.04 (on travis-ci)
* win-builder (devel and release)

## R CMD check results
No ERRORs, WARNINGs, or NOTEs

## Downstream dependencies
This is a major release and all package maintainers with packages depending on
'loo' were notified two months ago of a few changes that affect their packages.
We also sent another reminder more recently and will also notify the maintainers
again after this update is up on CRAN so that they know to go ahead and submit
their own updates if necessary. At a minimum, new versions of 'rstanarm' and
'brms' will be submitted immediately following the acceptance of this 'loo'
update.
