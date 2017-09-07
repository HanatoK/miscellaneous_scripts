#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
fes <- read.table(args[1])
outputname <- args[2]
require(akima)
require(plot3D)
x <- fes$V1*180/pi
y <- fes$V2*180/pi
z <- fes$V3/4.184
s <- interp(x, y, z, nx=720, ny=720)
png(filename = outputname, res = 300, height = 2000, width = 2000, bg = "transparent")
image2D(s)
dev.off()
