#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
fes <- read.table(args[1])
outputname <- args[2]
require(plot3D)
x <- unique(fes$V1)
y <- unique(fes$V2)
z <- fes$V3
binx <- length(x)
biny <- length(y)
matz <- matrix(z, nrow=binx, byrow=TRUE)
png(filename = outputname, res = 300, height = 2000, width = 2000, bg = "transparent")
image2D(matz, x=x, y=y, resfac=5)
dev.off()
