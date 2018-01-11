##### Basket Slurry Analysis #####
SlurryData<-read.csv("SubSamFEn17OK.csv", stringsAsFactors = F) #bring in volume data; bucket volume, subsample, filter, enclosure
SlurryData$BucketVol<-as.numeric(paste(SlurryData$BucketVol))
Treat<-read.csv("FEn17OKTreatments.csv", sep=",")#bring in enclosure and treatment data

##### Ash Free Dry Mass Data #####
AFDMraw<-read.csv("AFDMFEn17OK.csv", stringsAsFactors = F) #datasheet with dried and burned results
FilterPre<-read.csv("FiltPreWFen17OK.csv") #filter preweights
colnames(FilterPre)[1]<-"Filter."

AFDM2<-merge.data.frame(AFDMraw, FilterPre, by="Filter.")
AFDM<-merge.data.frame(AFDM2, SlurryData, by="Filter.")
AFDM$Tin.filter.combusted<-as.numeric(paste(AFDM$Tin.filter.combusted))

AFDM$SubstrateVol<-AFDM$Basket.*.03315
AFDM$Tare<-AFDM$Tin.Weight+AFDM$PreWeight
AFDM$DryMass<-AFDM$Tin.filter.dry-AFDM$Tare
AFDM$OrganicM<-AFDM$Tin.filter.dry-AFDM$Tin.filter.combusted
AFDM$InOrganicM<-AFDM$DryMass-AFDM$OrganicM
AFDM$AFDMdensity<-(AFDM$OrganicM/AFDM$VolFilrwe)*(AFDM$BucketVol/AFDM$SubstrateVol)
AFDM$Treat<-Treat[match(AFDM$Enclosure, Treat$ï..Enclosure),2]

library(ggplot2)
ggplot(AFDM, aes(x=Treat, y=AFDMdensity))+geom_boxplot()+facet_wrap(~AFDM$Week)
  
#### Chlorophyll ####
#this is from the algal tiles generally
ChlAraw<-read.csv("ChlDataFEn17OK.csv")#bring in raw chlorophyll data
ChlAtile<-ChlAraw[ChlAraw$Type=="Tiles",]
ChlAtile$Treatment<-Treat[match(ChlAtile$Enclosure.,Treat$ï..Enclosure),2]
ChlAtile$Area=0.0005905

ChlAtile$ChlAdensity<-26.7*((ChlAtile$X664nm-ChlAtile$fir750nm)-(ChlAtile$X665nm-ChlAtile$sec750nm))*ChlAtile$Vacetone/(ChlAtile$Area*1)/1000

ggplot(ChlAtile, aes(x=Treatment, y=ChlAdensity))+geom_boxplot()+facet_wrap(~Week)

ChlAfil<-ChlAraw[ChlAraw$Type=="Filter",]
ChlAfil$BucketVol<-SlurryData[match(ChlAfil$Order, SlurryData$FilterN), ]
ChlAfil$FilterVol<-SlurryData[match(ChlAfil$Order, SlurryData$FilterN), ]
ChlAfil$Area<-SlurryData[match(ChlAfil$Order, SlurryData$FilterN), ]*


ChlAfil$ChlAdensity<-26.7*((ChlAtile$664nm-ChlAtile$fir750nm)-(ChlAtile$665nm-ChlAtile$sec750nm))*ChlAtile$VAcetone/(ChlAtile$Area*1)/1000

#### Combining all the data ####
ChlAtile$EnW<-paste(ChlAtile$Enclosure., ChlAtile$Week, sep="")
AFDM$EnW<-paste(AFDM$Enclosure, AFDM$Week, sep="")

basalres<-merge.data.frame(ChlAtile, AFDM,by="EnW")
basalres$depth<-EnclosureRaster[match(basalres$Enclosure, EnclosureRaster$enc),8]
View(basalres)

test<-lm(AFDMdensity~Treat*Week.x, data=basalres)
anova(test)
summary(test)
leastm<-lsmeans(test, "Treat",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)
