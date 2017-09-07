#!/usr/bin/env Rscript 
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1])
fn = args[2]
x = data$V1 * 10
y = data$V2
shift = min(y)
y = (y - shift)/4.184
png(filename=fn,res=300,width=2000,height=2000,bg="transparent")
par(mar=c(5, 5, 4, 2)+0.1)
plot(x,y,xlab="CV",ylab="△G(kcal/mol)",type="l",lwd=2,cex.lab=1.5,cex.axis=1.2,tcl=(0.3))
box(lwd=2)
dev.off()
comparepmf <- data.frame(x,y)
write.table(comparepmf, file = "comp.pmf", col.names = FALSE, row.names = FALSE, sep = '\t')
