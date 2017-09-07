#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
rmsd <- function(a, b) {
  print(sqrt(mean((a - b)^2)))
}

ref <- read.table(args[1])
i <- 2
resultvector <- character(length(args) - 1)
filevector <- character(length(args) - 1)
for(i in 2:length(args)) {
  print(args[i])
  x <- read.table(args[i])
  filevector[i - 1] <- args[i]
  resultvector[i - 1] <- rmsd(ref$V3, x$V3)
}

print(resultvector)
out <- data.frame(filevector, resultvector)
write.table(out, file = "fesrmsd.dat", quote = FALSE, col.names = FALSE, row.names = FALSE)