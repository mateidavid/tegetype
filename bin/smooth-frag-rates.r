#!/usr/bin/Rscript
filename <- commandArgs(TRUE)[1]
t <- read.table(filename)
rg <- as.character(unique(t$V1))
t$V5 <- t$V3 / (t$V4 + 1)
m <- matrix(nrow=dim(t[t$V1==rg[1],])[1], ncol=0)
for (i in 1:length(rg)) {
    a <- array(unlist(lowess(t[t$V1==rg[i],][,c('V2','V5')], f=0.3)[2]))
    for (i in 1:dim(a)) {
        a[i] <- max(a[i], 1e-10)
    }
    m <- cbind(m, a)
}
t$V6 <- c(m)
write.table(t[,c("V1", "V2", "V6")], quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
