.onAttach <- function(libname, pkgname) {
    packageStartupMessage("\n*** Genoset API Changes ***\nThe genoset class has transitioned to the
RangedSummarizedExperiment API from the eSet API (e.g. use colnames instead of sampleNames). ***")
}
