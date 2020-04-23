.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("Welcome to the", pkgname, "package!"))
}