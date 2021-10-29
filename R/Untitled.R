p <- seq(from=0.001, to =0.05, by=0.001)
library(lattice)
df <- data.frame(p=p, BFupper=-1/(exp(1)*p*log(p)))
labp <- c(0.05, 0.01, 0.005, 0.001)
xyplot(BFupper~log(p), type='l', scales=list(x=list(at=log(labp),
                                                    labels=paste(labp))), data=df)


