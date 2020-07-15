## Test environments
* local R installation on Windows, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* Fedora Linux, R-devel, clang, gfortran (fedora-clang-devel)
* Debian Linux, R-devel, GCC ASAN/UBSAN (linux-x86_64-rocker-gcc-san)
* Oracle Solaris 10, x86, 32 bit, R-release (solaris-x86-patched)
* Ubuntu Linux 16.04 LTS, R-release, GCC (ubuntu-gcc-release)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (windows-x86_64-devel)
* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit (windows-x86_64-oldrel)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit (windows-x86_64-release)

## R CMD check results

0 errors | 0 warnings | 1 note

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Dylan Beijers <dylanbeijers@gmail.com>'
  Days since last update: 2

  Added a patch for macOS. On macOS, assignment of nested lists is handled a bit differently than on other platforms. For instance, if `first <- list()`, then assigning `first$second[[["third"]] <- 1` returns on Windows (and other platforms) a list named `first`, that contains another list named `second`, that contains a named element `third` equal to an unnamed length 1 numerical vector. On macOS though, it returns a list named `first`, that contains a named element `second` equal to a named (`third`) length 1 numerical vector. So `second` is a list on the other platforms, while being a named numerical vector on macOS. This caused a bug on macOS while computing standard errors. My apologies for not noticing this before submitting v0.2.0.