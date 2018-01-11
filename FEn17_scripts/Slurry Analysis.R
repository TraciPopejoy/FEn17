##### Basket Slurry Analysis #####
SlurryData<-read.csv("./FEn17_data/SubSamFEn17OK.csv", stringsAsFactors = F) #bring in volume data; bucket volume, subsample, filter, enclosure
SlurryData$BucketVol<-as.numeric(paste(SlurryData$BucketVol)) #tell volume to stop being a factor/text
colnames(SlurryData)[6]<-"WeekBad"
SlurryData[SlurryData$WeekBad=="8",8]<-"w09"
SlurryData[SlurryData$WeekBad=="4",8]<-"w04"
SlurryData[SlurryData$WeekBad=="12",8]<-"w12"
colnames(SlurryData)[8]<-"Week"

Treat<-read.csv("./FEn17_data/FEn17OKTreatments.csv", sep=",")#bring in enclosure and treatment data

##### Ash Free Dry Mass Data #####
AFDMraw<-read.csv("./FEn17_data/AFDMFEn17OK.csv", stringsAsFactors = F) #datasheet with dried and burned results
FilterPre<-read.csv("./FEn17_data/FiltPreWFen17OK.csv") #filter preweights
colnames(FilterPre)[1]<-"Filter."

AFDM2<-merge.data.frame(AFDMraw, FilterPre, by="Filter.")
AFDM1<-merge.data.frame(AFDM2, SlurryData, by="Filter.")
AFDM1$Tin.filter.combusted<-as.numeric(paste(AFDM$Tin.filter.combusted))
AFDM1$BucketVol<-as.numeric(paste(AFDM$BucketVol))

AFDM1$SubstrateVol<-AFDM$Basket.*.03315
AFDM1$Tare<-AFDM$Tin.Weight+AFDM$PreWeight
AFDM1$DryMass<-AFDM$Tin.filter.dry-AFDM$Tare
AFDM1$OrganicM<-AFDM$Tin.filter.dry-AFDM$Tin.filter.combusted
AFDM1$InOrganicM<-AFDM$DryMass-AFDM$OrganicM
AFDM1$AFDMdensity<-(AFDM$OrganicM/AFDM$VolFilrwe)*(AFDM$BucketVol/AFDM$SubstrateVol)

##NOTE: need to deal with the NAs here - from lost buckets, lost filters, water column

library(plyr)
AFDM<-ddply(AFDM1, .variables=c("TEid"),.fun=function(x) mean(x[,20]))
AFDM$Treat<-Treat[match(AFDM$Enclosure, Treat$ï..Enclosure),3]

library(ggplot2)
ggplot(AFDM, aes(x=Treat, y=AFDMdensity))+geom_boxplot()+facet_wrap(~AFDM$Week)
  
#### Chlorophyll ####
#this is from the algal tiles generally
ChlAraw<-read.csv("./FEn17_data/ChlDataFEn17OK.csv")#bring in raw chlorophyll data
colnames(ChlAraw)[2]<-"WeekBad"
ChlAraw[ChlAraw$WeekBad=="9",12]<-"w09"
ChlAraw[ChlAraw$WeekBad=="4",12]<-"w04"
ChlAraw[ChlAraw$WeekBad=="12",12]<-"w12"
colnames(ChlAraw)[12]<-"Week"

ChlAtile<-ChlAraw[ChlAraw$Type=="Tiles",]
ChlAtile$Treatment<-Treat[match(ChlAtile$Enclosure.,Treat$ï..Enclosure),3]
ChlAtile$Area=0.0005905

ChlAtile$ChlAdensity<-26.7*((ChlAtile$X664nm-ChlAtile$fir750nm)-(ChlAtile$X665nm-ChlAtile$sec750nm))*ChlAtile$Vacetone/(ChlAtile$Area*1)/1000

ggplot(ChlAtile, aes(x=Treatment, y=ChlAdensity))+geom_boxplot()+facet_wrap(~Week)

ChlAfil<-ChlAraw[ChlAraw$Type=="Filter",]
ChlAfil$BucketVol<-SlurryData[match(ChlAfil$Order, SlurryData$Filter.), 2]
ChlAfil$FilterVol<-SlurryData[match(ChlAfil$Order, SlurryData$Filter.), 3]
ChlAfil$Area<-SlurryData[match(ChlAfil$Order, SlurryData$Filter.), 5] * .03315
ChlAfil$Treatment<-Treat[match(ChlAfil$Enclosure.,Treat$ï..Enclosure),3]

ChlAfil$ChlAdensity<-26.7*((ChlAfil$X664nm-ChlAfil$fir750nm)-(ChlAfil$X665nm-ChlAfil$sec750nm))*ChlAfil$Vacetone/(ChlAfil$Area*1)/1000

ggplot(ChlAfil, aes(x=Treatment, y=ChlAdensity))+geom_boxplot()+facet_wrap(~Week)

####NOTE: we have chlA for the water column, noted by #N/A enclosure, filter volumes are somewhere

#### Combining all the data ####
ChlAtile$Enc2<-Treat[match(ChlAtile$Enclosure., Treat$ï..Enclosure),2]
ChlAtile$TEid<-paste(ChlAtile$Week, ChlAtile$Enc2,sep="")
AFDM$Enc2<-Treat[match(AFDM$Enclosure, Treat$ï..Enclosure),2]
AFDM$TEid<-paste(AFDM$Week, AFDM$Enc2, sep="")

basalres<-merge.data.frame(ChlAtile, AFDM,by="TEid")
##NOTE: reduce number of column by removing superflous columns

View(basalres)