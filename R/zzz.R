.onAttach <- function(libname, pkgname) {
    packageStartupMessage("\n*** Genoset API Changes ***\nThe genoset class is transitioning to extending 
RangedSummarizedExperiment rather than eSet. For this release, 
please use the RSE API as the eSet API has been deprecated\n (e.g. colnames instead of sampleNames). ***")
}
