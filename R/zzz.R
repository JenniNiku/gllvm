.onLoad <- function(libname, pkgname) {
  # by default tape sequentially to avoid memory spikes
  # TMB::config(tape.parallel = FALSE, DLL = "gllvm")
}

.onUnLoad <- function(libpath) {
  # reset TMB config settings
  # TMB::config(tape.parallel = TRUE, DLL = "gllvm")
}