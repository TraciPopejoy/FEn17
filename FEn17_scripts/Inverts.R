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
Counts<-Counts[-743,]

#converting to density to compensate for different sampling effort
head(SlurryData) #found in Slurry Analysis sheet
Counts$Density.npm<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5]*.03315)
Counts$Density.npb<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5])


ggplot(data=Counts, aes(x=Treatment, y=Density.npm, color=Order)) + 
  geom_point(position="jitter") +
  scale_y_log10()

#############     Field Invert Biomass Calculation     #############
Inv$Family<-as.character(TaxaList$Family[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Order<-as.character(TaxaList$Order[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Length<-Inv$Length.cm*10

Inv$FFG<-InvGraph[match(Inv$Taxa, InvGraph$Taxa), "FFG"]
Inv$Type<-InvGraph[match(Inv$Treatment, InvGraph$Treatment), "Type"]

ggplot(Inv[!is.na(Inv$Treatment),], 
       aes(x=FFG, y=Length))+
  geom_violin()+scale_y_sqrt()+facet_wrap(~Type)+fungraph

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

ggplot(na.omit(InvTotalBM), aes(x=Taxa, y=Sum.mg, color=Treatment))+
  geom_point()+coord_flip()+scale_y_log10()

#converting mean biomass to biomass/meter
InvTotalBM$Density<-InvTotalBM$Sum.mg/(SlurryData[match(InvTotalBM$TEid, SlurryData$TEid),5]*.03315)

###would use density because sampling was not constant (not always full basket recovery)
InvGraph<-merge(InvTotalBM[,-c(4,5)],Counts, by=c("TEid","Taxa","Treatment","Enc","Week"))
InvGraph<-merge(InvGraph, TaxaList[,c(1,2,5:9)], by="Taxa")
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
InvGraph$Treatment<-factor(InvGraph$Treatment, levels=c("ACTL","ACTS","AMBL","AMBS","CTRL"))
TropTable<-data.frame(TropN=seq(1:6),
                      FFG=c("C-Gatherer","C-Filterer",
                            "Herbivore","Predator","Shredder","Parasite"))
InvGraph$FFG<-TropTable[match(InvGraph$T.Trop, TropTable$TropN),2]

fungraph<-theme(axis.text.x=element_text(angle = 90,size=12,color="black"),
      axis.text.y = element_text(size=12,color="black"),
      axis.title.y=element_text(size=20),
      plot.background = element_blank(),
      panel.border=element_blank(),
      panel.grid.major= element_line(colour=NA), 
      panel.grid.minor=element_line(colour=NA),
      title=element_text(size=20),
      panel.background = element_rect(fill = "white"),
      axis.line.x=element_line(colour="black"),
      axis.line.y=element_line(colour="black"),
      strip.background=element_rect(fill="white", color="black"),
      strip.text=element_text(size=15))

ggplot(InvGraph, 
       aes(x=FFG, y=mean.length, color=Treatment))+
  scale_y_sqrt() +
  geom_point(aes(size=Density.npm), alpha=.5, position="jitter")+
  ylab("Mean Length (mm)")+xlab("Functional Feeding Group")+
  facet_grid(.~Type, space="free", scales = "free")+fungraph


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

sizegraph<-merge(InvB, TaxaList)
sizegraph$Type<-InvGraph[match(sizegraph$Treatment, InvGraph$Treatment), "Type"]
ggplot(sizegraph[sizegraph$T.Trop==2,], aes(x=Treatment, y=Length))+
  geom_violin()+facet_wrap(~Type, scales="free")

ggplot(InvGraph, 
       aes(x=Type, y=mean.length))+
  scale_y_sqrt() + geom_violin(aes(fill=Treatment))+
  geom_boxplot(alpha=.3, aes(group=Treatment))+
  #geom_point(aes(size=Density.npb), alpha=.3)+
  ylab("Mean Length (mm)")+xlab("Functional Feeding Group")+
  facet_grid(.~FFG, space="free", scales = "free")+fungraph+
  scale_fill_brewer(palette="Paired")

ggplot(InvGraph, 
       aes(x=Type, y=Density.npm, fill=Treatment))+
  scale_y_log10() +
  geom_boxplot()+
  ylab("Density (n/meter.sq)")+xlab("Treatments")+
  facet_grid(.~FFG, space="free", scales = "free")+
  scale_fill_brewer(palette="Paired")+fungraph

(p1<-ggplot(InvGraph[InvGraph$T.Trop==2,], 
       aes(x=Type, y=Density.npm, fill=Treatment))+
  geom_boxplot()+
  ylab("Density (n/meter.sq)")+xlab("Treatments")+ylim(c(0,160))+
  facet_grid(.~FFG, space="free", scales = "free")+
  scale_fill_brewer(palette="Paired")+fungraph)

(p2<-ggplot(InvGraph[InvGraph$T.Trop==2,], 
       aes(x=Type, y=mean.length, fill=Treatment))+
  geom_boxplot()+
  ylab("Mean Length (mm)")+xlab("Treatments")+
  facet_grid(.~FFG, space="free", scales = "free")+
  scale_fill_brewer(palette="Paired")+fungraph)
grid.arrange(p1,p2,ncol=2)

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
InvSum$Enclosure<-Treat[match(InvSum$Enc, Treat$Enclosure2),1]

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

test<-InvGraph %>% group_by(TEid)%>%filter(T.Trop==2) %>%
  mutate(CFiltDen=sum(Density.npm))
InvSum$CFiltDen<-test[match(InvSum$TEid, test$TEid), "CFiltDen"]
InvSum[is.na(InvSum$CFiltDen),"CFiltDen"]<-0
InvSum$CFiltDen<-unlist(InvSum$CFiltDen)
mathss<-as.data.frame(InvSum)


testing<-aov(CFiltDen~Type, data=mathss)
summary(testing)
plot(testing)
TukeyHSD(testing)
library(lsmeans)
leastm<-lsmeans(testing, "Type",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

test<-InvGraph %>% group_by(TEid)%>%filter(T.Trop==2) %>%
  mutate(CFiltBM=mean(mean.length))
InvSum$CFiltBM<-test[match(InvSum$TEid, test$TEid), "CFiltBM"]
InvSum[is.na(InvSum$CFiltBM),"CFiltBM"]<-0
InvSum$CFiltBM<-unlist(InvSum$CFiltBM)
mathss2<-as.data.frame(InvSum)

testing<-aov(CFiltBM~Type, data=mathss2)
summary(testing)
plot(testing)
TukeyHSD(testing)
library(lsmeans)
leastm<-lsmeans(testing, "Type",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)


write.csv(InvSum, "FEn17week12insects.csv")


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
  fungraph


library(vegan)
nmds<-metaMDS(commat)
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
trait1<-ordtrait[,c(3,4)]

library(reshape2)
AbMatrixT<-dcast(Counts, TEid + Treatment ~ Taxa, value.var = "Density.npm")
AbMatrixT[is.na(AbMatrixT)]<-0
Abundances<-AbMatrixT[,-c(1,2)]
rownames(Abundances)<-AbMatrixT[,1]

rownames(ordtrait)==colnames(Abundances) #need it to be true to run the function

library(FD)
ex <- dbFD(trait1,Abundances)

FunciGraph<-data.frame(FDis=ex$FDis,
                       FRich=ex$FRic,
                       FEven=ex$FEve,
                       RaoQ=ex$RaoQ)
FunciGraph$TEid<-rownames(FunciGraph)
FunciGraph$Treatment<-InvGraph[match(FunciGraph$TEid, InvGraph$TEid),"Treatment"]
FunciGraph[FunciGraph$Treatment=="AMBL","Type"]<-"Live"
FunciGraph[FunciGraph$Treatment=="ACTL","Type"]<-"Live"
FunciGraph[FunciGraph$Treatment=="AMBS","Type"]<-"Sham"
FunciGraph[FunciGraph$Treatment=="ACTS","Type"]<-"Sham"
FunciGraph[FunciGraph$Treatment=="CTRL","Type"]<-"Ctrl"
FunciGraph[FunciGraph$Treatment=="AMBL","Spp"]<-"AMB"
FunciGraph[FunciGraph$Treatment=="ACTL","Spp"]<-"ACT"
FunciGraph[FunciGraph$Treatment=="AMBS","Spp"]<-"AMB"
FunciGraph[FunciGraph$Treatment=="ACTS","Spp"]<-"ACT"
FunciGraph[FunciGraph$Treatment=="C?ATRL","Spp"]<-"Ctrl"
FunciGraph$Treatment<-factor(FunciGraph$Treatment, c("ACTL","ACTS","AMBL","AMBS","CTRL"))
FunciGraph$Type<-factor(FunciGraph$Type, c("Live","Sham","Ctrl"))

fit <- aov(FDis ~ Type, data=FunciGraph)
summary(fit)
plot(fit)
TukeyHSD(fit)
leastm<-lsmeans(fit, "Type",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)



mFGraph<-melt(FunciGraph)

ggplot(mFGraph, aes(x=Type, y=value, fill=Treatment))+geom_boxplot()+
  facet_wrap(~variable, scales = "free")+
  scale_fill_brewer(palette = "Paired")+fungraph

ModelData<-merge(InvSum, FunciGraph, by=c("TEid","Treatment","Type","Spp"))
ModelData$depth<-EncDV[match(ModelData$Enclosure, EncDV$ï..Enclosure), "Depth.m"]
ModelData$velocity<-EncDV[match(ModelData$Enclosure, EncDV$ï..Enclosure), "V.mps"]
ModelData$Treatment<-factor(ModelData$Treatment, 
                            levels=c("ACTL","ACTS","AMBL","AMBS","CTRL"))
#ModelData$Treatment<-relevel(ModelData$Treatment, ref="CTRL")

library(car)
library(compute.es)
library(effects)
library(multcomp) 
#Explore the data & check ANOVA assumptions

#Check for homogeneity
leveneTest(ModelData$FDis, ModelData$Treatment, center=median)

#Check for independence between treatment and covariate 
testco1<-aov(depth~Treatment, data=ModelData)
summary(testco1)
testco2<-aov(depth~Type, data=ModelData)
summary(testco2)
testco3<-aov(velocity~Type, data=ModelData)
summary(testco3)
#Run the ANCOVA 
mod8 <- lm(FDis~depth+velocity+Type, data=ModelData)
Anova(mod8, type="III")
adjustedMeans<-effect("Type", mod8)
postHocs<-glht(mod8, linfct=mcp(Type="Tukey"))
summary(postHocs)
#Run contrasts/comparisons 
#Check for homogeneity of regression slopes
hist(residuals(mod8), col="darkgray")
plot(fitted(mod8), residuals(mod8))

####MacroAbund####
leveneTest(ModelData$Density.npm, ModelData$Type, center=median)
#Check for independence between treatment and covariate 

#Run the ANCOVA 
mod9 <- lm(Density.npm~depth+Type, data=ModelData)
Anova(mod9, type="III")
adjustedMeans<-effect("Type", mod9)
postHocs<-glht(mod9, linfct=mcp(Type="Tukey"))
summary(postHocs)


##NMDS
DenC<-dcast(Counts[,-c(2,4,5,6,7,8,9)], Enc~...)
DenC[is.na(DenC)]<-0
rownames(DenC)<-DenC[,1]
DenC<-DenC[,-1]
comNMDS<-metaMDS(DenC)
plot(comNMDS, type="t")
library(tidyverse)
library(plyr)
Fcom<-ddply(InvGraph, .variables = c("Enc","FFG"), .fun = function(x){sum(x$Density.npb)})
FunC<-dcast(Fcom, Enc~FFG)
FunC[is.na(FunC)]<-0
rownames(FunC)<-FunC[,1]
FunC<-FunC[,-1]
funNMDS<-metaMDS(FunC)
plot(funNMDS, type="t")
funNMDS1<-funNMDS$points

EnclosureRaster$funORD<-funNMDS1[match(EnclosureRaster$enc, rownames(funNMDS1)),1]
plot(EnclosureRaster["funORD"])
text(cc[,1],cc[,2],zc)
