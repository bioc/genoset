#!/usr/bin/env Rscript
#####
#####  A script to run tests on the genoset package

suppressMessages(library(genoset))
suppressMessages(library(RUnit))
suppressMessages(library(BSgenome))

options(warn=1)

testsuite.genoset <- defineTestSuite("genoset.test.suite", 
                                       dirs = "unit",
                                       testFileRegexp = "^runit.+\\.R", 
                                       testFuncRegexp = "^test.+", 
                                       rngKind = "default",
                                       rngNormalKind = "default")

results <- capture.output(runTestSuite(testsuite.genoset))

## Finish up
quit(runLast=F)
