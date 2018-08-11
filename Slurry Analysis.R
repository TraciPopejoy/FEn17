##### Basket Slurry Analysis #####
#bring in volume data; bucket volume, subsample, filter, enclosure
SlurryData<-read.csv("./FEn17_data/SubSamFEn17OK.csv", stringsAsFactors = F) 
SlurryData$BucketVol<-as.numeric(paste(SlurryData$BucketVol)) #tell volume to stop being a factor/text
colnames(SlurryData)[6]<-"WeekBad"
SlurryData[SlurryData$WeekBad=="8",8]<-"w09"
SlurryData[SlurryData$WeekBad=="4",8]<-"w04"
SlurryData[SlurryData$WeekBad=="12",8]<-"w12"
colnames(SlurryData)[8]<-"Week"
#bring in enclosure and treatment data
Treat<-read_excel("./FEn17_data/FEn17OKTreatments.xlsx")

SlurryData$Enc2<-Treat[match(SlurryData$Enclosure, Treat$ï..Enclosure),"Enclosure2"]
SlurryData$TEid<-paste(SlurryData$Week,SlurryData$Enc2, sep="")
#Some filters didn't have bucket volumes associated with them
#5L is written next to them and is the mean bucket volume (considering sig dig)
SlurryData[is.na(SlurryData$BucketVol),"BucketVol"]<-5000
#Some filters didn't have baskets associated with them
#assumed 3 baskets
SlurryData[is.na(SlurryData$Basket.),"Basket."]<-3


##### Ash Free Dry Mass Data #####
AFDMraw<-read.csv("./FEn17_data/AFDMFEn17OK.csv", stringsAsFactors = F) #datasheet with dried and burned results
FilterPre<-read.csv("./FEn17_data/FiltPreWFen17OK.csv") #filter preweights
colnames(FilterPre)[1]<-"Filter."

AFDM2<-merge.data.frame(AFDMraw, FilterPre, by="Filter.")
AFDM<-merge.data.frame(AFDM2, SlurryData, by="Filter.")
AFDM$Tin.filter.combusted<-as.numeric(paste(AFDM$Tin.filter.combusted))
AFDM$BucketVol<-as.numeric(paste(AFDM$BucketVol))

AFDM$SubstrateVol<-AFDM$Basket.*331.5 #centimeter squared
AFDM$Tare<-AFDM$Tin.Weight+AFDM$PreWeight
AFDM$DryMass<-AFDM$Tin.filter.dry-AFDM$Tare
AFDM$OrganicM<-AFDM$Tin.filter.dry-AFDM$Tin.filter.combusted
AFDM$InOrganicM<-AFDM$DryMass-AFDM$OrganicM

#AFDMdensity units g/cm3
AFDM$AFDMdensity<-(AFDM$OrganicM/AFDM$VolFilrwe)*(AFDM$BucketVol/AFDM$SubstrateVol)

#following filters were water column samples, AFDM density is g / ml
for(i in 1:length(WC<-c("FA062","FA063","FA132","FA133","FA284","FA285"))){
  AFDM[AFDM$Filter.==WC[i],"Enclosure"]<-"WaterColumn"
  AFDM[AFDM$Filter.==WC[i],"AFDMdensity"]<-AFDM[AFDM$Filter.==WC[i],"OrganicM"]/AFDM[AFDM$Filter.==WC[i],"VolFilrwe"]
}

#adding graphing variables on the table
AFDM$Enc2<-Treat[match(AFDM$Enclosure, Treat$ï..Enclosure),2]
AFDM$TEid<-paste(AFDM$Week, AFDM$Enc2, sep="")
AFDM$Treatment<-Treat[match(AFDM$Enclosure, Treat$ï..Enclosure),3]

library(ggplot2)
ggplot(AFDM, aes(x=Treatment, y=AFDMdensity))+geom_boxplot()+
  labs(y="AFDM g/m3")+facet_wrap(~AFDM$Week)
  
#### Chlorophyll ####
#this is from the algal tiles generally
ChlAraw<-read.csv("./FEn17_data/ChlDataFEn17OK.csv", stringsAsFactors = F)#bring in raw chlorophyll data
colnames(ChlAraw)[2]<-"WeekBad"
ChlAraw[ChlAraw$WeekBad=="9",12]<-"w09"
ChlAraw[ChlAraw$WeekBad=="4",12]<-"w04"
ChlAraw[ChlAraw$WeekBad=="12",12]<-"w12"
colnames(ChlAraw)[12]<-"Week"

ChlAtile<-ChlAraw[ChlAraw$Type=="Tiles",]
ChlAtile$Treatment<-Treat[match(ChlAtile$Enclosure.,Treat$ï..Enclosure),3]
((2.75/2)^2)*pi
ChlAtile$Area=5.939574 # centimeter squared

#chl mg/cm2
ChlAtile$ChlAdensity<-26.7*((ChlAtile$X664nm-ChlAtile$fir750nm)-(ChlAtile$X665nm-ChlAtile$sec750nm))*(ChlAtile$Vacetone/ChlAtile$Area)*1

ggplot(ChlAtile, aes(x=Treatment, y=ChlAdensity))+geom_boxplot()+facet_wrap(~Week)

#when tiles were unavailable, we instead used filters to capture ChlA
#these have NAs that need fixing (bleh)
#myster filters FC010, FC017 and FC018 (17 & 18 might be E7)
ChlAfil<-ChlAraw[ChlAraw$Type=="Filter",]
ChlAfil[ChlAfil$Order=="D76","Order"]<-ChlAfil[ChlAfil$Order=="D76","Enclosure."]
ChlAfil[ChlAfil$Order=="D75","Order"]<-ChlAfil[ChlAfil$Order=="D75","Enclosure."]
ChlAfil$BucketVol<-SlurryData[match(ChlAfil$Order, SlurryData$Filter.), 2]
ChlAfil$FilterVol<-SlurryData[match(ChlAfil$Order, SlurryData$Filter.), 3]
ChlAfil$Area<-SlurryData[match(ChlAfil$Order, SlurryData$Filter.), "Basket."] * 331.5 #cm2
ChlAfil$Enc3<-SlurryData[match(ChlAfil$Order,SlurryData$Filter.),1]
ChlAfil$Treatment<-Treat[match(ChlAfil$Enc3, Treat$ï..Enclosure),3]
#5L is written next to them and is the mean bucket volume (considering sig dig)
ChlAfil[is.na(ChlAfil$BucketVol),"BucketVol"]<-5000

ChlAfil$ChlAdensity<-26.7*((ChlAfil$X664nm-ChlAfil$fir750nm)-(ChlAfil$X665nm-ChlAfil$sec750nm))*(ChlAfil$Vacetone/ChlAfil$FilterVol)*(ChlAfil$BucketVol/ChlAfil$Area)*1

for(k in which(ChlAfil$Enclosure.=="WC")){
  ChlAfil[k,"ChlAdensity"]<-26.7*((ChlAfil$X664nm[k]-ChlAfil$fir750nm[k])-(ChlAfil$X665nm[k]-ChlAfil$sec750nm[k]))*(ChlAfil$Vacetone[k]/ChlAfil$FilterVol[k])*1
}

ggplot(ChlAfil, aes(x=Treatment, y=ChlAdensity))+geom_boxplot()+facet_wrap(~Week)

#### Combining all the data ####
ChlAtile$Enc2<-Treat[match(ChlAtile$Enclosure., Treat$ï..Enclosure),2]
ChlAtile$TEid<-paste(ChlAtile$Week, ChlAtile$Enc2,sep="")
ChlAfil$Enc2<-Treat[match(ChlAfil$Enclosure., Treat$ï..Enclosure),2]
ChlAfil$TEid<-paste(ChlAfil$Week, ChlAfil$Enc2,sep="")
#USE ME TO FILL IN DATA GAPS
ChlAfil2<-ddply(ChlAfil, c("TEid","Week"), summarize, mean.Chl=mean(ChlAdensity, na.rm=T), Chlsd=sd(ChlAdensity, na.rm=T))
library(plyr)
AFDMT<-ddply(AFDM, .variables = c("TEid"), .fun=function(x) {data.frame(TEid=x$TEid[1],
                                                                        Treatment=x$Treatment[1],
                                                                        Week=x$Week[1],
                                                                        Enc2=x$Enc2[1],
                                                                        AFDMg.cm2=mean(x$AFDMdensity, na.rm=T),
                                                                        sdAFDM=sd(x$AFDMdensity,na.rm=T))})

basalres<-merge.data.frame(ChlAtile, AFDMT,by=c("TEid","Enc2","Week","Treatment"), all=T)
basalres<-basalres[,-c(5:16)]
colnames(basalres)[5]<-"ChlA.ug.cm2"
basalres$ChlA.g.cm2<-basalres$ChlA.ug.cm2/1e6
basalres$AutoIndex<-(basalres$AFDMg.cm2*1000)/(basalres$ChlA.ug.cm2*.001)
View(basalres)

##### Chlorophyll on Shells week 12 #####
ChlAShell<-ChlAraw[ChlAraw$Type=="shellsampling",]
ChlAShell$BucketID<-SlurryData[match(ChlAShell$Enclosure., SlurryData$Filter.), "Type"]
ChlAShell$Enc<-substring(ChlAShell$BucketID, 1,2)
ChlAShell$Spp<-substring(ChlAShell$BucketID, 3,5)
ChlAShell[ChlAShell$Order=="D69",15]<-'AMB'
ChlAShell[ChlAShell$Order=="D65",15]<-'ACT'
ChlAShell[ChlAShell$Order=="D84",15]<-'AMB'
ChlAShell[ChlAShell$Order=="D86",15]<-'ACT'

ChlAShell$FilterVol<-SlurryData[match(ChlAShell$Enclosure., SlurryData$Filter.), 3]
ChlAShell$BucketVol<-SlurryData[match(ChlAShell$Enclosure., SlurryData$Filter.), 2]

MusselData<-read.csv("./FEn17_data/MusselBMExpFEn17OK.csv")
MusselData$ShellArea<-(MusselData$L*MusselData$H*2)*0.01
shellareaag<-ddply(MusselData, .variables = c('Genus','Unit'), .fun=function(x) sum(x$ShellArea,na.rm=T))
shellareaag$ty<-paste(shellareaag$Unit, shellareaag$Genus, sep="")

ChlAShell$ShellArea<-shellareaag[match(ChlAShell$BucketID, shellareaag$ty), "V1"] #cm2
ChlAShell$Enc2<-SlurryData[match(ChlAShell$Enclosure.,SlurryData$Filter.),9]
ChlAShell$Treatment<-Treat[match(ChlAShell$Enc2, Treat$Enclosure2),3]

ChlAShell$ChlAdensity<-26.7*((ChlAShell$X664nm-ChlAShell$fir750nm)-(ChlAShell$X665nm-ChlAShell$sec750nm))*(ChlAShell$Vacetone/ChlAShell$FilterVol)*(ChlAShell$BucketVol/ChlAShell$ShellArea)*1
library(ggplot2)
ggplot(ChlAShell, aes(x=Spp, y=ChlAdensity, group=Enc2))+
  geom_line()+facet_wrap(~Treatment)+theme_bw()



install.packages('profileR')
library('profileR')


