#libraries
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(xlsx)
library(vegan)

#############     Invert Data in Field Surber Samples    ##############
Treat<-read.csv("FEn17OKTreatments.csv", sep=",")#bring in enclosure and treatment data
TaxaList<-read.xlsx("NSFMacroInvertTaxaList.xlsx", sheetIndex = 1)
BiomassReg<-read.xlsx("Macroinv Power Law Coeffs TBP.xlsx", sheetIndex = 1, stringsAsFactors=F)

FEn17Inv<-read.csv("FEn17InvMeas.csv")
Inv<-FEn17Inv

#clean the data frame
Inv<-Inv[-c(1:3),-c(3:5)]
Inv$TEid<-substring(Inv$Label, 6,11)
Inv$Enc<-substring(Inv$Label,9,11)
Inv$Week<-substring(Inv$Label,6,8)
Inv$Treatment<-Treat[match(Inv$Enc, Treat$Enclosure2), 3]
colnames(Inv)[4]<-"Taxa"

#summarize data
Counts<-ddply(Inv, .variables = c("TEid","Taxa"), .fun=function(x) {
  data.frame(TEid=x[1,5],
            Enc=x[1,6],
            Week=x[1,7],
            Treatment=x[1,8],
            n=count(x,x[1,4])[,2])
  })

Counts$Family<-as.character(TaxaList$Family[match(Counts$Taxa,TaxaList$Taxa)])
Counts$Order<-as.character(TaxaList$Order[match(Counts$Taxa,TaxaList$Taxa)])

ggplot(data=Counts, aes(x=Treatment, y=log10(n), color=Order)) + geom_point(cex=4)

#############     Field Invert Biomass Calculation     #############
Inv$Family<-as.character(TaxaList$Family[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Order<-as.character(TaxaList$Order[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Length<-Inv$Length.cm*10

#removing taxa that give me trouble 
InvA<-Inv[!Inv$Order=="misc",]
InvB<-na.omit(InvA)


InvBM<-ddply(.data=InvB, .var=c("Taxa"), .fun=function(x) {
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
      TEid=x$TEid[i],
      Enc=x$Enc[i],
      Week=x$Week[i],
      Treatment=x$Treatment[i],
      length=x$Length[i],
      neq=length(idx2),
      ninR=sum(idx2),
      meanBM.mg=mean(indmassest),
      median=median(indmassest),
      stDev=sd(indmassest))}), 
    data.frame)
})

InvTotalBM<-ddply(InvBM, .var=c("TEid", "Taxa"),
                  .fun=function(x){data.frame(mean.mg=mean(x$meanBM.mg, na.rm=T), 
                                              median.mg=median(x$meanBM.mg, na.rm=T), 
                                              SD.mg=sd(x$meanBM.mg, na.rm=T),
                                              Sum.mg=sum(x$meanBM.mg,na.rm=T),
                                              Treatment=x$Treatment[1],
                                              Enc=x$Enc[1],
                                              Week=x$Week[1]
                                              )}) 

ggplot(na.omit(InvTotalBM), aes(x=Taxa, y=log10(Sum.mg), color=Treatment))+
  geom_point()+coord_flip()


#converting mean biomass to biomass/meter
SlurryData<-read.csv("SubSamFEn17OK.csv", stringsAsFactors = F)
SlurryData$Enc2<-Treat[match(SlurryData$Enclosure, Treat$ï..Enclosure),2]
SlurryData$TEid<-paste("w",SlurryData$Week,SlurryData$Enc2, sep="")
InvTotalBM$Density<-InvTotalBM$Sum.mg/(SlurryData[match(InvTotalBM$TEid, SlurryData$TEid),5]*.03315)

###would use density because sampling was not constant (not always full basket recovery)


###Total Summary of Data####
InvSumA<-ddply(Counts,.variables=c("TEid"),.fun=function(x) {count(x,x[1,1])[,2]})
commat<-dcast(Counts[,-c(3,4,5,7,8)], TEid~...)
commat[is.na(commat)]<-0
rownames(commat)<-commat[,1]
commat<-commat[,-1]
InvSumA$Shannon<-diversity(commat)
InvSumB<-ddply(InvTotalBM,.variables=c("TEid"),.fun=function(x) data.frame(TotalBM.mg=sum(x$Sum.mg),
                                                                           BMDensity.mgpm2=sum(x$Density)))
InvSum<-merge(InvSumA,InvSumB, by="TEid")
colnames(InvSum)[2]<-"richness"
InvSum$Enc<-Inv[match(InvSum$TEid,Inv$TEid),6]
InvSum$Week<-Inv[match(InvSum$TEid,Inv$TEid),7]
InvSum$Treatment<-Treat[match(InvSum$Enc, Treat$Enclosure2), 3]
InvSum$basketn<-SlurryData[match(InvSum$TEid, SlurryData$TEid),5]

ggplot(InvSum, aes(x=Treatment, y=log10(BMDensity.mgpm2)))+
  geom_point()


############## writing the data for Caryn ######################

totalBMx<-field.totalBM[,-c(1,5,8)]
fTuncast<-dcast(totalBMx, Site+SamplingSeason+Surber+Treatment~Taxa, mean)
write.table(fTuncast, file="FieldTotalBM.csv",sep=",")

avgBMx<-totalBMx[,c(1:3,5:6,4)]
fTuncast<-dcast(temp, Site+SamplingSeason+Surber+Treatment~Taxa, mean)
write.table(fTuncast, file="Field 15 Average Biomass.csv",sep=",")



