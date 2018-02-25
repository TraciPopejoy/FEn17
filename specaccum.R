library(xlsx)
SACdata<-read.xlsx("./FEn17_data/w19FirstSamples.xlsx",sheetIndex = 1)
SACdata$SampleT<-paste(SACdata$Enc, SACdata$SampleSplit, sep=".")

CMT<-t(table(SACdata[,c(7,10)]))

library(vegan)
spE9 <- specaccum(CMT[grep("E9", rownames(CMT)),], "collector")
spD9 <- specaccum(CMT[grep("D9", rownames(CMT)),], "collector")
spA5 <- specaccum(CMT[grep("A5", rownames(CMT)),], "collector")
spB6 <- specaccum(CMT[grep("B6", rownames(CMT)),], "collector")
spA3 <- specaccum(CMT[grep("A3", rownames(CMT)),], "collector")

plot(spD9, col="pink", lwd=2)
plot(spE9, col="red", lwd=2, add=T)
plot(spB6,  col="orange", lwd=2, add=T)
plot(spA5, col="yellow", lwd=2, add=T)
plot(spA3, col="green", lwd=2, add=T)

spE92 <- specaccum(CMT[grep("E9", rownames(CMT)),])
spD92 <- specaccum(CMT[grep("D9", rownames(CMT)),])
spA52 <- specaccum(CMT[grep("A5", rownames(CMT)),])
spB62 <- specaccum(CMT[grep("B6", rownames(CMT)),])
spA32 <- specaccum(CMT[grep("A3", rownames(CMT)),])

plot(spD92, col="pink", lwd=2, ci.type="bar")
plot(spE92, col="red", lwd=2, add=T,ci.type="bar")
plot(spB62,  col="orange", lwd=2, add=T,ci.type="bar")
plot(spA52, col="yellow", lwd=2, add=T,ci.type="bar")
plot(spA32, col="green", lwd=2, add=T,ci.type="bar")

S <- specnumber(CMT) # observed number of species
(raremax <- min(rowSums(CMT)))
Srare <- rarefy(CMT, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(CMT, step = 20, sample = raremax, col = "blue", cex = 0.6)

