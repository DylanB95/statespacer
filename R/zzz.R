.onUnload <- function (libpath) {
  library.dynam.unload("statespacer", libpath)
}
