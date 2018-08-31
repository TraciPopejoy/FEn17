#libraries
library(sp)
library(xlsx)
library(plyr)
library(RColorBrewer)

EnclosureRaster = data.frame(enc=
                      c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10",
                      "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10",
                      "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10",
                      "D1","D2","D3","D4","D5","D6","D7","D8","D9","D10",
                      "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10"),
                z = seq(1:50),
                xc = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10)),
                yc = c(rep(seq(1:10),5)))
coordinates(EnclosureRaster) = ~xc+yc #promotes to SpatialPointsDataFrame
gridded(EnclosureRaster) = TRUE #promotes to SpatialPixelsDataFrame
df = as(EnclosureRaster, "SpatialGridDataFrame") # to raster data

EnclosureRaster$Treatment<-Treat[match(EnclosureRaster$enc, Treat$ï..Enclosure), "TreatA"]

# draw labels to verify:
cc = coordinates(df)
a=df[["enc"]]
zc=as.character(a)
zc[is.na(zc)]="NA"
text(cc[,1],cc[,2],zc)

#### Spatial data ####
EncDVw09<-read.csv("./FEn17_data/EncPhysDisFEn17OK.csv") #table has depth and velocity for each enclosure
EncDVw09$Enc2<-Treat[match(EncDVw09$ï..Enclosure, Treat$ï..Enclosure), "Enclosure2"]
EnclosureRaster$depthw9<-EncDVw09[match(EnclosureRaster$enc, EncDVw09$ï..Enclosure),2]
EnclosureRaster$velocityw9<-EncDVw09[match(EnclosureRaster$enc, EncDVw09$ï..Enclosure),3]

plot(EnclosureRaster["depthw9"], col=brewer.pal(5,"Blues"))
text(cc[,1],cc[,2],zc)

image(EnclosureRaster["velocityw9"])
text(cc[,1],cc[,2],zc)

Treat$TreatN<-as.numeric(as.factor(Treat$TreatA))
EnclosureRaster$TreatN<-Treat[match(EnclosureRaster$enc, Treat$ï..Enclosure),5]
image(EnclosureRaster["TreatN"])
text(cc[,1],cc[,2],zc)

CHltemp<-ChlAtile[ChlAtile$Week=="w12",]
EnclosureRaster$Chl12<-CHltemp[match(EnclosureRaster$enc, CHltemp$Enclosure.),"ChlAdensity"]
EnclosureRaster@data[13,"Chl12"]<-NA
image(EnclosureRaster["Chl12"])
text(cc[,1],cc[,2],zc)

### Velocity Transects ###

Rvel<-read_excel("./FEn17_data/FieldEncDataSum2017V2.xlsx", sheet = 1)
Rvel %>% group_by(Transect) %>% summarize(meandepth=mean(depth.m, na.rm=T),
                                          depthsd=sd(depth.m, na.rm=T),
                                          velmean=mean(`velocity.m/s`, na.rm=T),
                                          sdvel=sd(`velocity.m/s`, na.rm=T))


####Water Depth Data####
library(tidyverse)
PD1<-read.csv("./FEN17_data/10699194.csv", header=F, stringsAsFactors = F)
PD2<-read.csv("./FEN17_data/beg10699194_5.csv", header=F, stringsAsFactors = F)

PDa<-PD1[-c(1:2),]
names(PDa)<-c("number","DateTimeGMT","AbsP.psi","TempF")
PDb<-PD2[-c(1:2),]
names(PDb)<-c("number","DateTimeGMT","AbsP.psi","TempF")

pressuredata<-rbind(PDb[,1:4],PDa[,1:4])
pressuredata<-pressuredata[pressuredata$TempF!="",]
pressuredata<-pressuredata[-c(10554:10557),]

library(lubridate)
pressuredata$GoodDate<-mdy_hms(pressuredata[,2], tz="GMT")
pressuredata$Date<-date(pressuredata$GoodDate)
pressuredata$AbsP.psi<-as.numeric(pressuredata$AbsP.psi)
pressuredata$AbsP.Nm2<-pressuredata$AbsP.psi*6894.744825
pressuredata$TempF<-as.numeric(pressuredata$TempF)

head(pressuredata)
ggplot(pressuredata, aes(x=GoodDate, y=AbsP.psi))+geom_point()

#pressure = density * gravity * height  
#pressure / (9.81m/s*density)=height
#1 psi = 6894.744825 N/m2
#on 7/13/2017, sediment to midpoint=18cm, midpoint to top=68cm total height is 86cm
#7/20/2017	13.7	42.8	56.5
#10/8/2017 water depth was 48.5cm
Jul13<-mean(pressuredata[pressuredata$Date=="2017-07-13","AbsP.Nm2"], na.rm=T)
(WdensityJ13<-Jul13/(9.81*.86))
(Jul20<-mean(pressuredata[pressuredata$Date=="2017-07-20","AbsP.Nm2"], na.rm=T))
(WdensityJ20<-Jul20/(9.81*.565))
Oct8<-mean(pressuredata[pressuredata$Date=="2017-10-07","AbsP.Nm2"], na.rm=T)
(WdensityO<-Oct8/(9.81*.485))
(WDensity<-mean(WdensityJ13, WdensityJ20, WdensityO))

pressuredata$Height<-pressuredata$AbsP.Nm2/(9.81*WDensity)

library(scales)
fmt_dcimals <- function(decimals=0){
  function(x) as.character(round(x,decimals))
}

(SiteDepth<-ggplot(pressuredata[pressuredata$Height>=.58,], aes(x=GoodDate, y=Height))+
  geom_line()+
  geom_vline(xintercept = ymd_hms("2017-09-18 12:00:00"), color="red",alpha=.3,cex=5)+
  geom_vline(xintercept = ymd_hms("2017-10-07 12:00:00"), color="red",alpha=.3,cex=5)+
  ylab("Water Depth (m)")+xlab("Date")+
  scale_x_datetime(date_breaks="2 weeks", date_labels="%b %d")+
  scale_y_continuous(labels=fmt_dcimals(2)))+theme_bw()
ggsave("SiteDepth.tiff",SiteDepth,width=7, height=4, dpi=300)

##### MUSSEL BIOMASS ESTIMATE ####
head(MusselData)
mreg<-read.xlsx("./FEn17_data/LENGTH-MASS_CLA_CCV_20161208-reg coefficients.xlsx", sheetIndex = 2)
head(mreg)
mreg$Spp<-c("boot","all","ACT","AMB","OREF","POCC","QVER","FFLA",NA,NA,NA)
MusselData$BiomassE<-mreg[match(MusselData$Genus, mreg$Spp),"a"]*(MusselData$L^mreg[match(MusselData$Genus, mreg$Spp),"b"]) 
ACTbm<-mean(MusselData[MusselData$Genus=="AMB","BiomassE"], na.rm=T)
#replacing not available mussel data with average
MusselData[MusselData$Genus=="ACT" & is.na(MusselData$BiomassE),"BiomassE"]<-mean(MusselData[MusselData$Genus=="ACT","BiomassE"], na.rm=T)
MusselData[MusselData$Genus=="AMB" & is.na(MusselData$BiomassE),"BiomassE"]<-mean(MusselData[MusselData$Genus=="AMB","BiomassE"], na.rm=T)

MusBiomass<-MusselData %>% group_by(Enc2,Genus, Treatment) %>% 
  summarize(sumBM=sum(BiomassE), meanBM=mean(BiomassE), sdBM=sd(BiomassE))
MBM<-MusselData %>% group_by(Enc2, Treatment) %>% 
  summarize(sumBM=sum(BiomassE), meanBM=mean(BiomassE), sdBM=sd(BiomassE),
            ACTbm=mean(BiomassE[Genus=="ACT"]), 
            AMBbm=mean(BiomassE[Genus=="AMB"]))
ggplot(MusBiomass, aes(x=Treatment, y=sumBM))+geom_boxplot()
ggplot(MBM, aes(x=Treatment, y=sumBM))+geom_boxplot()
names(MBM)[1]<-"Enc"
ggsave("./Figures/biomass.png")

#### Structure data ####
MusselData<-read.csv("./FEn17_data/MusselBMExpFEn17OK.csv")
colnames(MusselData)[1]<-"Enclosure"
MusselData$Enc2<-Treat[match(MusselData$Enclosure,Treat$ï..Enclosure),2]
expos<-ddply(MusselData, .variables =c("Enclosure"), .fun=function(x) mean(na.omit(x$Exposure)))
EncDVw09$Avg.Exposure<-expos[match(EncDVw09$ï..Enclosure, expos$Enclosure),2]
EnclosureRaster$exposure<-EncDVw09[match(EnclosureRaster$enc, EncDVw09$ï..Enclosure),"Avg.Exposure"]
plot(EnclosureRaster["exposure"])
text(cc[,1],cc[,2],zc)