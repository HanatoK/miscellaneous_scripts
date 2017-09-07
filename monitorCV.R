#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) 
colvars <- read.table(args[1])
ncol <- args[2]
outputname <- "CVevolution"
i <- 1
interval <- 100000
start <- 0
x <- colvars[,as.numeric(ncol)]
end <- length(x) - 1
repeat {
  fn = paste(outputname, as.character(i), ".png", sep = "")
  png(filename = fn, res = 300, width = 10000, height = 1000, bg = "transparent")
  ep = start + interval
  plot(x[start:ep], pch = ".", type = "l", ylim = c(-1.0, 1.0))
  dev.off()
  start = start + interval
  i = i + 1
  if(start + interval >= end) {
    break
  }
}