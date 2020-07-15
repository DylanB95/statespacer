## Test environments
* local R installation on Windows, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)
* win-builder (old release)
* win-builder (release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on rhub)
* Fedora Linux, R-devel, clang, gfortran (on rhub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (on rhub)

## R CMD check results

0 errors | 0 warnings | 1 note

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Dylan Beijers <dylanbeijers@gmail.com>'
  Days since last update: 2

  Added a patch for macOS. On macOS, assignment of nested lists is handled a bit differently than on other platforms. For instance, if `first <- list()`, then assigning `first$second[[["third"]] <- 1` returns on Windows (and other platforms) a list named `first` containing a list named `second` containing a named element `third` equal to an unnamed length 1 numerical vector. On macOS though, it returns a list named `first` containing a named element `second` equal to a named (`third`) length 1 numerical vector. So `second` is a list on the other platforms, while being a named numerical vector on macOS. This caused a bug on macOS while computing standard errors. My apologies for not noticing this before submitting v0.2.0.