#!/usr/bin/env Rscript

library(GenomicRanges)
library(SummarizedExperiment)
library(genoset)

sample.names = LETTERS[11:13]
probe.names = paste("p", 1:1000, sep = "")
num.samples = length(sample.names)
num.probes = length(probe.names)

locs = GRanges(ranges = IRanges(start = c(seq(from = 1.25e+08, by = 30000, length = 400), 
    seq(from = 1, length = 400, by = 32500), seq(from = 3e+07, length = 200, by = 30000)), 
    width = 1, names = probe.names), seqnames = factor(c(rep("chr8", 400), rep("chr12", 
    400), rep("chr17", 200)), levels = c("chr8", "chr12", "chr17")))

genome(locs) = "hg19"
fake.baf = sample(c(0, 0.5, 1), length(probe.names), replace = TRUE) + rnorm(length(probe.names), 
    0, 0.01)
fake.baf[fake.baf > 1] = fake.baf[fake.baf > 1] - 1
fake.baf[fake.baf < 0] = fake.baf[fake.baf < 0] + 1
hets = fake.baf < 0.75 & fake.baf > 0.25
fake.baf[101:200][hets[101:200]] = fake.baf[101:200][hets[101:200]] + rep(c(-0.2, 
    0.2), 50)[hets[101:200]]
fake.cn = matrix(c(c(rnorm(200, 0.4, 0.05), rnorm(200, 0.2, 0.05), rnorm(200, 0.23, 
    0.05), rnorm(200, 0.15, 0.05), rnorm(200, 2, 0.05)), c(rnorm(200, 0, 0.05), rnorm(200, 
    3, 0.05), rnorm(200, 14, 0.05), rnorm(200, 0.1, 0.05), rnorm(200, -0.05, 0.05)), 
    c(rnorm(200, 0.1, 0.05), rnorm(200, 1, 0.05), rnorm(200, -0.5, 0.05), rnorm(200, 
        3, 0.05), rnorm(200, 3, 0.05))), nrow = num.probes, ncol = num.samples, dimnames = list(probe.names, 
    sample.names))
fake.baf = matrix(fake.baf, nrow = num.probes, ncol = num.samples, dimnames = list(probe.names, 
    sample.names))
fake.pData = data.frame(matrix(LETTERS[1:15], nrow = 3, ncol = 5, dimnames = list(sample.names, 
    letters[1:5])))

assaylist = list(lrr = fake.cn, baf = fake.baf)

genoset.ds = GenoSet(locs, assaylist, fake.pData)

save(genoset.ds, file = "genoset.RData", compression_level = 9, compress = "xz")
