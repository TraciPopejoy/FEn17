#libraries
library(sp)
library(xlsx)
library(plyr)

EnclosureRaster = data.frame(enc=
                      c("A01","A02","A03","A04","A05","A06","A07","A08","A09","A10",
                      "B01","B02","B03","B04","B05","B06","B07","B08","B09","B10",
                      "C01","C02","C03","C04","C05","C06","C07","C08","C09","C10",
                      "D01","D02","D03","D04","D05","D06","D07","D08","D09","D10",
                      "E01","E02","E03","E04","E05","E06","E07","E08","E09","E10"),
                z = seq(1:50),
                xc = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10)),
                yc = c(rep(seq(1:10),5)))
coordinates(EnclosureRaster) = ~xc+yc #promotes to SpatialPointsDataFrame
gridded(EnclosureRaster) = TRUE #promotes to SpatialPixelsDataFrame
df = as(EnclosureRaster, "SpatialGridDataFrame") # to raster data


EnclosureRaster$Treatment<-Treat[match(EnclosureRaster$enc, Treat$Enclosure2), 3]

# draw labels to verify:
cc = coordinates(df)
a=df[["enc"]]
zc=as.character(a)
zc[is.na(zc)]="NA"
text(cc[,1],cc[,2],zc)

#### Spatial data ####
EncDV<-read.csv("./FEn17_data/EncPhysDisFEn17OK.csv") #table has depth and velocity for each enclosure
EncDV$Enc2<-Treat[match(EncDV$ï..Enclosure, Treat$ï..Enclosure), "Enclosure2"]
EnclosureRaster$depth<-EncDV[match(EnclosureRaster$enc, EncDV$ï..Enclosure),2]
EnclosureRaster$velocity<-EncDV[match(EnclosureRaster$enc, EncDV$ï..Enclosure),3]

Vidtemp<-Viddata[1:50,]
EnclosureRaster$fishpres<-Vidtemp[match(EnclosureRaster$enc, Vidtemp$Row.Labels),9]

EnclosureRaster$fishN<-Vidtemp[match(EnclosureRaster$enc, Vidtemp$Row.Labels),7]

image(EnclosureRaster["fishN"])
text(cc[,1],cc[,2],zc)

image(EnclosureRaster["depth"])
text(cc[,1],cc[,2],zc)

image(EnclosureRaster["velocity"])
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

#### Structure data ####
MusselData<-read.csv("./FEn17_data/MusselBMExpFEn17OK.csv")
colnames(MusselData)[1]<-"Enclosure"
MusselData$Enc2<-Treat[match(MusselData$Enclosure,Treat$ï..Enclosure),2]
expos<-ddply(MusselData, .variables =c("Enclosure"), .fun=function(x) mean(na.omit(x$Exposure)))
EncDV$Avg.Exposure<-expos[match(EncDV$ï..Enclosure, expos$Enclosure),2]
EnclosureRaster$exposure<-EncDV[match(EnclosureRaster$enc, EncDV$ï..Enclosure),4]
image(EnclosureRaster["exposure"])
text(cc[,1],cc[,2],zc)

###Water Depth Data###
PD1<-read.csv("./FEN17_data/10699194.csv", header=F, stringsAsFactors = F)
PD2<-read.csv("./FEN17_data/beg10699194_5.csv", header=F, stringsAsFactors = F)

PDa<-PD1[-c(1:2),]
names(PDa)<-c("number","DateTimeGMT","AbsP.psi","TempF")
PDb<-PD2[-c(1:2),]
names(PDb)<-c("number","DateTimeGMT","AbsP.psi","TempF")

pressuredata<-rbind(PDb[,1:4],PDa[,1:4])
pressuredata<-pressuredata[pressuredata$TempF!="",]
pressuredata<-pressuredata[-c(10554:10557),]

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
  geom_vline(xintercept = ymd_hms("2017-08-07 12:00:00"), color="red",alpha=.3,cex=5)+
  geom_vline(xintercept = ymd_hms("2017-09-18 12:00:00"), color="red",alpha=.3,cex=5)+
  geom_vline(xintercept = ymd_hms("2017-10-07 12:00:00"), color="red",alpha=.3,cex=5)+
  fungraph+ylab("Water Depth (m)")+xlab("Date")+
  scale_x_datetime(date_breaks="2 weeks", date_labels="%b %d")+
  scale_y_continuous(labels=fmt_dcimals(2)))
ggsave("SiteDepth.tiff",SiteDepth,width=7, height=4, dpi=300)
