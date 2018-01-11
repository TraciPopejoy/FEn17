library(sp)
library(xlsx)

EnclosureRaster = data.frame(enc=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","A10",
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


EnclosureRaster$Treatment<-Treat[match(EnclosureRaster$enc, Treat$ï..Enclosure), 2]

#### Spatial data ####
EncDV<-read.csv("EncPhysDisFEn17OK.csv") #table has depth and velocity for each enclosure
EnclosureRaster$depth<-EncDV[match(EnclosureRaster$enc, EncDV$ï..Enclosure),2]
EnclosureRaster$velocity<-EncDV[match(EnclosureRaster$enc, EncDV$ï..Enclosure),3]
image(EnclosureRaster["depth"])
image(EnclosureRaster["velocity"])

# draw labels to verify:
cc = coordinates(df)
a=df[["enc"]]
zc=as.character(a)
zc[is.na(zc)]="NA"
text(cc[,1],cc[,2],zc)
