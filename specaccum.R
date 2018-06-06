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
spB5 <- specaccum(CMT[grep("B5", rownames(CMT)),], "collector")
spC10 <- specaccum(CMT[grep("C10", rownames(CMT)),], "collector")

plot(spD9, col="pink", lwd=2, ylab="Species")
plot(spE9, col="red", lwd=2, add=T)
plot(spB6,  col="orange", lwd=2, add=T)
plot(spA5, col="yellow", lwd=2, add=T)
plot(spA3, col="green", lwd=2, add=T)
plot(spC10, col="blue", lwd=2, add=T)
plot(spB5, col="purple", lwd=2, add=T)

spE92 <- specaccum(CMT[grep("E9", rownames(CMT)),])
spD92 <- specaccum(CMT[grep("D9", rownames(CMT)),])
spA52 <- specaccum(CMT[grep("A5", rownames(CMT)),])
spB62 <- specaccum(CMT[grep("B6", rownames(CMT)),])
spA32 <- specaccum(CMT[grep("A3", rownames(CMT)),])
spC102 <- specaccum(CMT[grep("C10", rownames(CMT)),])
spB52 <- specaccum(CMT[grep("B5", rownames(CMT)),])


plot(spD92, col="pink", lwd=2, ci.type="bar")
plot(spE92, col="red", lwd=2, add=T,ci.type="bar")
plot(spB62,  col="orange", lwd=2, add=T,ci.type="bar")
plot(spA52, col="yellow", lwd=2, add=T,ci.type="bar")
plot(spA32, col="green", lwd=2, add=T,ci.type="bar")
plot(spC102, col="blue", lwd=2, add=T,ci.type="bar")
plot(spB52, col="purple", lwd=2, add=T,ci.type="bar")


#pulling in week 12 data for those samples
InvSPA1<-rbind(Inv[Inv$Enc=="D09",], Inv[Inv$Enc=="E09",],Inv[Inv$Enc=="B06",],
               Inv[Inv$Enc=="A05",],Inv[Inv$Enc=="A03",],Inv[Inv$Enc=="C10",],
               Inv[Inv$Enc=="B05",])
CMTInv<-t(table(InvSPA1[,c(9,6)]))
CMTall<-t(table(Inv[,c(9,6)]))

S <- specnumber(CMT) # observed number of species
(raremax <- min(rowSums(CMT)))
Srare <- rarefy(CMT, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(CMT, step = 20, sample = raremax, col = "blue", cex = 0.6)
(raremax2 <- min(rowSums(CMTInv)))
rarecurve(CMTInv, step = 20, sample = raremax2)
(raremax3 <- min(rowSums(CMTall)))
rarecurve(CMTall, step = 20)

SPall<-specaccum(CMTall, "exact")
plot(SPall)

colnames(SACdata)[7]<-"Taxa"
SACdata$Family<-as.character(TaxaList$Family[match(SACdata$Taxa,TaxaList$Taxa)])
SACdata$Order<-as.character(TaxaList$Order[match(SACdata$Taxa,TaxaList$Taxa)])
SACdata$Length<-SACdata$Length*10

#removing taxa that give me trouble 
SACdataA<-SACdata[!SACdata$Order=="misc",]
SACdataB<-na.omit(SACdataA)


SACdataBM<-ddply(.data=SACdataB, .var=c("Taxa"), .fun=function(x) {
  idx<-x[1,c("Family","Order")]
  if(idx$Family=="misc"){
    plcoeffs<-BiomassReg[BiomassReg$Order == idx$Order &
                           !is.na(BiomassReg$Order == idx$Order),]
  }else{
    plcoeffs<-BiomassReg[BiomassReg$Family==idx$Family&
                           BiomassReg$Order == idx$Order&
                           !is.na(BiomassReg$Family==idx$Family), ]  
  }
  ldply(lapply(1:dim(x)[[1]],FUN=function(i) {
    idx2<- c(x$Length[i]>=plcoeffs[,19] & x$Length[i]<=plcoeffs[,20]) 
    idx2[is.na(idx2)]<-T
    d1<-plcoeffs[idx2,]
    indmassest<-d1$a*(x$Length[i]^d1$b) 
    data.frame(
      Enc=x$Enc[i],
      Week=9,
      SampleT=x$SampleT[i],
      length=x$Length[i],
      neq=length(idx2),
      ninR=sum(idx2),
      meanBM.mg=mean(indmassest),
      median=median(indmassest),
      stDev=sd(indmassest))}), 
    data.frame)
})
library(plyr)
SACdataTotalBM<-ddply(SACdataBM, .var=c("SampleT", "Taxa"),
                  .fun=function(x){data.frame(mean.mg=mean(x$meanBM.mg, na.rm=T), 
                                              median.mg=median(x$meanBM.mg, na.rm=T), 
                                              SD.mg=sd(x$meanBM.mg, na.rm=T),
                                              Sum.mg=sum(x$meanBM.mg,na.rm=T),
                                              SampleT=x$SampleT[1],
                                              Enc=x$Enc[1],
                                              mean.length=mean(x$length, na.rm=T)
                  )}) 
library(ggplot2)
ggplot(na.omit(SACdataTotalBM), aes(x=Enc, y=log10(Sum.mg), color=Enc))+
  geom_point()
SlurryTemp<-SlurryData[SlurryData$Week=="w09",]
SACdataTotalBM$Density<-SACdataTotalBM$Sum.mg/(SlurryTemp[match(SACdataTotalBM$Enc, SlurryTemp$Enclosure),5]*.03315)
SACdataSUM<-ddply(SACdataTotalBM,.variables=c("SampleT"),
                  .fun=function(x) data.frame(Enc=x$Enc[1],
                                              TotalBM.mg=sum(x$Sum.mg, na.rm=T),
                                              BMDensity.mgpm2=sum(x$Density, na.rm=T),
                                              BMcorrected=(sum(x$Density, na.rm=T)*4)))

ggplot(SACdataSUM, aes(x=Enc, y=BMcorrected))+geom_point(cex=2)+theme_krementz()+ylab("BM mgpm2 corrected")
ggplot(SACdataSUM, aes(x=Enc, y=TotalBM.mg))+geom_point(cex=2)+theme_krementz()

SACdataTotalBM<-SACdataTotalBM[order(SACdataTotalBM$SampleT),]

SACsumEn<-ddply(SACdataSUM,.variables=c("Enc"),.fun=function(x) data.frame(REALBMdensity.mg=sum(x$BMDensity.mgpm2)))
SACsumEn$variable<-"RealBMDensity.mgpm2"
SACsumEn$SampleT<-NA
colnames(SACsumEn)[2]<-"value"
SACsumEn<-SACsumEn[,c(4,1,3,2)]
library(reshape2)
SacGraph<-melt(SACdataSUM[,c(1,2,5)])
SacGraph<-rbind(SacGraph, SACsumEn)

ggplot(SacGraph, aes(x=Enc, y=value, color=variable))+geom_point(cex=3)+theme_krementz()
