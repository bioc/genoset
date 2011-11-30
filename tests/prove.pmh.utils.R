#!/usr/bin/env Rscript
#####
#####  A script to run tests on the genoset package

suppressMessages(library(genoset))
suppressMessages(library(RUnit))
suppressMessages(library(bigmemory))

options(warn=1)
RUnit = options("RUnit")$RUnit
RUnit$silent = TRUE
options("RUnit"=RUnit)

testsuite.genoset <- defineTestSuite("genoset.test.suite", 
                                       dirs = "unit",
                                       testFileRegexp = "^runit.+\\.R", 
                                       testFuncRegexp = "^test.+", 
                                       rngKind = "default",
                                       rngNormalKind = "default")

results <- capture.output(runTestSuite(testsuite.genoset))

## Finish up
quit(runLast=F)
