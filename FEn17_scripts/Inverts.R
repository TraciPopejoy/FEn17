#libraries
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(xlsx)
library(vegan)

#############     Invert Data in Field Surber Samples    ##############
#bring in enclosure and treatment data, trait/taxa list, and length weight regressions
Treat<-read.csv("./FEn17_data/FEn17OKTreatments.csv", sep=",", stringsAsFactors = F) 
TaxaList<-read.csv("./FEn17_data/TaxaTable.csv", sep=",", stringsAsFactors = F)
BiomassReg<-read.xlsx("./FEn17_data/Macroinv Power Law Coeffs TBP.xlsx", sheetIndex = 1, stringsAsFactors=F)
#THIS IS THE INSECT DATA
FEn17Inv<-read.csv("./FEn17_data/FEn17InvMeas.csv", stringsAsFactors = F)
Inv<-FEn17Inv

#clean the data frame
Inv<-Inv[-c(1:3),-c(3:5)]
Inv$TEid<-substring(Inv$Label, 6,11)
Inv$Enc<-substring(Inv$Label,9,11)
Inv$Week<-substring(Inv$Label,6,8)
Inv$Treatment<-Treat[match(Inv$Enc, Treat$Enclosure2), 3]
colnames(Inv)[4]<-"Taxa"
sort(unique(Inv$Taxa)) #check to make sure no misspellings
Inv[Inv$Taxa=="Col.Elm.A",4]<-"Col.ElmA"

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

#converting to density to compensate for different sampling effort
head(SlurryData) #found in Slurry Analysis sheet
Counts$Density.npm<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5]*.03315)


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
                                              Week=x$Week[1],
                                              mean.length=mean(x$length, na.rm=T)
                                              )}) 

ggplot(na.omit(InvTotalBM), aes(x=Taxa, y=log10(Sum.mg), color=Treatment))+
  geom_point()+coord_flip()

#converting mean biomass to biomass/meter
InvTotalBM$Density<-InvTotalBM$Sum.mg/(SlurryData[match(InvTotalBM$TEid, SlurryData$TEid),5]*.03315)

###would use density because sampling was not constant (not always full basket recovery)
InvGraph<-merge(InvTotalBM[,-c(4,5)],Counts, by=c("TEid","Taxa","Treatment","Enc","Week"))
InvGraph<-merge(InvGraph, TaxaList[,c(1,5:9)], by="Taxa")
InvGraph[InvGraph$Treatment=="AMBL","Type"]<-"Live"
InvGraph[InvGraph$Treatment=="ACTL","Type"]<-"Live"
InvGraph[InvGraph$Treatment=="AMBS","Type"]<-"Sham"
InvGraph[InvGraph$Treatment=="ACTS","Type"]<-"Sham"
InvGraph[InvGraph$Treatment=="CTRL","Type"]<-"Ctrl"
InvGraph[InvGraph$Treatment=="AMBL","Spp"]<-"AMB"
InvGraph[InvGraph$Treatment=="ACTL","Spp"]<-"ACT"
InvGraph[InvGraph$Treatment=="AMBS","Spp"]<-"AMB"
InvGraph[InvGraph$Treatment=="ACTS","Spp"]<-"ACT"
InvGraph[InvGraph$Treatment=="CTRL","Spp"]<-"Ctrl"
InvGraph$Type<-factor(InvGraph$Type, levels=c("Live","Sham","Ctrl", ordered=T))

InvGraph<-InvGraph[InvGraph$Taxa!= "Tri.Leptoceridae",]
InvGraph<-InvGraph[InvGraph$Taxa!= "Dip.Other",]

ggplot(InvGraph, 
       aes(x=T.Trop, y=T.Habit, color=Type))+
  #scale_y_log10() +
  geom_point(aes(size=Density.npm), position="jitter")+
  theme_bw()#+facet_grid(.~Type, space="free", scales = "free")

ggplot(InvGraph, 
       aes(x=Order, y=Density.npm, fill=T.Trop))+
  scale_y_log10() +coord_flip()+
  geom_boxplot()+
  facet_wrap(~Type)+theme_classic()

ggplot(na.omit(InvGraph), 
       aes(x=Order, y=Density.npm, fill=Order))+
  scale_y_log10() + coord_flip() +
  geom_boxplot()+
  facet_wrap(~Type)+theme_classic()


###Total Summary of Data####
InvSumA<-ddply(Counts,.variables=c("TEid"),.fun=function(x) {count(x,x[1,1])[,2]})
commat<-dcast(Counts[,-c(3,4,5,7,8,9)], TEid~...)
commat[is.na(commat)]<-0
rownames(commat)<-commat[,1]
commat<-commat[,-1]
InvSumA$Shannon<-diversity(commat)
InvSumA$TotalN<-rowSums(commat)
InvSumB<-ddply(InvTotalBM,.variables=c("TEid"),.fun=function(x) data.frame(TotalBM.mg=sum(x$Sum.mg),
                                                                           BMDensity.mgpm2=sum(x$Density)))
InvSum<-merge(InvSumA,InvSumB, by="TEid")
colnames(InvSum)[2]<-"richness"
InvSum$Enc<-Inv[match(InvSum$TEid,Inv$TEid),6]
InvSum$Week<-Inv[match(InvSum$TEid,Inv$TEid),7]
InvSum$Treatment<-Treat[match(InvSum$Enc, Treat$Enclosure2), 3]
InvSum$basketn<-SlurryData[match(InvSum$TEid, SlurryData$TEid),5]
InvSum$Density.npm<-InvSum$TotalN/(InvSum$basketn*.03315)

InvSum[InvSum$Treatment=="AMBL","Type"]<-"Live"
InvSum[InvSum$Treatment=="ACTL","Type"]<-"Live"
InvSum[InvSum$Treatment=="AMBS","Type"]<-"Sham"
InvSum[InvSum$Treatment=="ACTS","Type"]<-"Sham"
InvSum[InvSum$Treatment=="CTRL","Type"]<-"Ctrl"
InvSum[InvSum$Treatment=="AMBL","Spp"]<-"AMB"
InvSum[InvSum$Treatment=="ACTL","Spp"]<-"ACT"
InvSum[InvSum$Treatment=="AMBS","Spp"]<-"AMB"
InvSum[InvSum$Treatment=="ACTS","Spp"]<-"ACT"
InvSum[InvSum$Treatment=="CTRL","Spp"]<-"Ctrl"
InvSum$Type<-factor(InvSum$Type, levels=c("Live","Sham","Ctrl", ordered=T))

RAcommat<-commat/rowSums(commat)
RAcommat$sites<-rownames(RAcommat)
RAcommat$Type<-InvSum[match(RAcommat$sites, InvSum$TEid),12]
ComGraph<-melt(RAcommat)
ComGraph$Order<-InvGraph[match(ComGraph$variable, InvGraph$Taxa),12]
ComGraph$Trop<-TaxaList[match(ComGraph$variable, TaxaList$Taxa),"T.Trop"]

commat1<-commat
commat1$sites<-rownames(commat1)
commat1$Type<-InvSum[match(commat1$sites, InvSum$TEid),12]
commat2<-melt(commat1)
commat2$Trop<-TaxaList[match(commat2$variable, TaxaList$Taxa),"T.Trop"]

ggplot(data = ComGraph, aes(x = sites, y = value, fill = Trop)) + 
  geom_bar(stat="identity") + coord_flip()+
  labs(x="Sites", y="Relative Abundance") +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))+
  facet_grid(Type~., space="free", scales="free")

library(wesanderson)
colors<-c(wes_palette("Zissou")[c(1,3,5)])
colors<-c("#3B9AB2","#EBCC2A","#F21A00")
ggplot(InvSum, aes(x=Treatment, y=BMDensity.mgpm2, color=Type))+
  geom_point(cex=5)+ylim(0,max(InvSum$BMDensity.mgpm2))+
  ylab("Biomass Density (mg/m sq.)")+
  scale_color_manual(values=colors)+stat_summary(color="black")+
  theme_light()

plot(nmds, display = c("sites", "species"), choices = c(1, 2),
     type = "n", shrink = FALSE)
points(nmds, display = c("sites", "species"),
       choices = c(1,2), shrink = FALSE)
## S3 method for class 'metaMDS'
text(nmds, display = c( "species"), labels=colnames(commat), 
     choices = c(1,2), cex=.5)
text(nmds, display = c( "sites"), labels=rownames(commat), 
     choices = c(1,2), cex=1, col="red")

##### Functional Diversity Analysis #####
trait<-TaxaList[,-c(1:4,10,11)]
rownames(trait)<-TaxaList[,1]
ordtrait<-trait[order(rownames(trait),decreasing=F),]
trait1<-as.matrix(ordtrait[,c(3,4)])

AbMatrixT<-dcast(Counts, TEid + Treatment ~ Taxa, value.var = "Density.npm")
AbMatrixT[is.na(AbMatrixT)]<-0
Abundances<-AbMatrixT[,-c(1,2)]
rownames(Abundances)<-AbMatrixT[,1]

rownames(ordtrait)==colnames(Abundances[,-25]) #need it to be true to run the function
Abundances1<-Abundances[,-25]

ex <- dbFD(trait1,Abundances1)

FuncGraph<-data.frame(FDis=ex$FDis)
FuncGraph$TEid<-rownames(FuncGraph)
FuncGraph$Treatment<-InvGraph[match(FuncGraph$TEid, InvGraph$TEid),"Treatment"]
FuncGraph<-FuncGraph[,c(2,3,1)]

ggplot(FuncGraph, aes(x=Treatment, y=FDis))+geom_boxplot()
