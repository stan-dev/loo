## Updating DESCRIPTION, NEWS
+ changed in DESCRIPTION: version and date
+ update NEWS.md

## Running reverse dependency checks
+ flocker 0.1.0 still fails but on GitHub is already an updated release 
(it is just not yet in Cran)


## Running `devtools::check_win_devel()`
+ Jonah, did you got a message with the tar file?

```
> withr::with_envvar(c("NOT_CRAN" = "true"), devtools::check_win_devel())
Building windows version of loo (2.10.0)
ℹ Using R-devel with win-builder.r-project.org.
Email results to jgabry@gmail.com?

1: I forget
2: Nope
3: Yeah

Selection: 
3
── R CMD build ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/u/21/wa.bocktif1/unix/GitHub2/loo/DESCRIPTION’ ...
─  preparing ‘loo’: (3.3s)
✔  checking DESCRIPTION meta-information ...
─  installing the package to build vignettes
✔  creating vignettes (49m 41.6s)
─  checking for LF line-endings in source and make files and shell scripts (358ms)
─  checking for empty or unneeded directories
─  building ‘loo_2.10.0.tar.gz’
   Warning: invalid uid value replaced by that for user 'nobody'
   Warning: invalid gid value replaced by that for user 'nobody'
   
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
ℹ Check <jgabry@gmail.com> for the results in 15-30 mins (~06:59 ).
```

## Running `devtools::check_mac_release()`
+ seems to be all fine
+ tar file is in this repo: `results.tar.bz2`

```
> withr::with_envvar(c("NOT_CRAN" = "true"), devtools::check_mac_release())
macOS builder results
Build system: r-release-macosx-arm64|4.6.0|macosx|macOS 26.2 (25C56)|Mac mini|Apple M1||en_US.UTF-8|macOS 14.4|clang-1700.6.3.2|GNU Fortran (GCC) 14.2.0

You can download a full tar ball results.tar.bz2 which contains the check results, tests, installed packages and binaries. A subset of the files are shown below.

* using log directory ‘/Volumes/PkgBuild/work/1781534234-1930b57b62ecdbc5/packages/sonoma-arm64/results/4.6/loo.Rcheck’
* using R version 4.6.0 Patched (2026-04-24 r89963)
* using platform: aarch64-apple-darwin23
* R was compiled by
    Apple clang version 17.0.0 (clang-1700.3.19.1)
    GNU Fortran (GCC) 14.2.0
* running under: macOS Tahoe 26.2
* using session charset: UTF-8
* current time: 2026-06-15 14:37:37 UTC
* checking for file ‘loo/DESCRIPTION’ ... OK
* checking extension type ... Package
* this is package ‘loo’ version ‘2.10.0’
* package encoding: UTF-8
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... NOTE
Found the following hidden files and directories:
  .git-blame-ignore-revs
  .agents
  .agents/skills/r-cli-app/.evals
These were most likely included in error. See section ‘Package
structure’ in the ‘Writing R Extensions’ manual.
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking whether package ‘loo’ can be installed ... [2s/2s] OK
* checking installed package size ... OK
* checking package directory ... OK
* checking ‘build’ directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking code files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... [0s/0s] OK
* checking whether the package can be loaded with stated dependencies ... [0s/0s] OK
* checking whether the package can be unloaded cleanly ... [0s/0s] OK
* checking whether the namespace can be loaded with stated dependencies ... [0s/0s] OK
* checking whether the namespace can be unloaded cleanly ... [0s/0s] OK
* checking loading without being on the library search path ... [0s/0s] OK
* checking whether startup messages can be suppressed ... [0s/0s] OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... [2s/2s] OK
* checking Rd files ... [0s/0s] OK
* checking Rd metadata ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of ‘data’ directory ... OK
* checking data for non-ASCII characters ... [0s/0s] OK
* checking LazyData ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking R/sysdata.rda ... OK
* checking installed files from ‘inst/doc’ ... OK
* checking files in ‘vignettes’ ... OK
* checking examples ... [1s/1s] OK
* checking for unstated dependencies in ‘tests’ ... OK
* checking tests ... [15s/8s] OK
  Running ‘testthat.R’ [15s/8s]
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... [10s/16s] OK
* checking PDF version of manual ... [3s/3s] OK
* DONE
Status: 1 NOTE
* using check arguments '--no-clean-on-error '
* elapsed time (check, wall clock): 0:51
* result classification: OK, errors:no, warnings:no, notes:yes
```