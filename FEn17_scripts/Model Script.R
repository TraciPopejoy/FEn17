#script to bring models in, run last

library(xlsx)
fishdata<-read.xlsx("./FEn17_data/videowithabiotics.xlsx", sheetIndex =1 )
colnames(fishdata)[c(3,10)]<-c("GTreat","Log(probFish+1)")
#adding columns that match other data
fishdata$Enc2<-Treat[match(fishdata$Unit, Treat$ï..Enclosure),2] #Treat found in Slurry Analysis script
fishdata$Treatment<-Treat[match(fishdata$Unit, Treat$ï..Enclosure),3]
fishdata$Week<-rep(NA, nrow(fishdata))
for(i in 1:nrow(fishdata)){
  if(fishdata$Month[i]=="Oct"){fishdata$Week[i]<-"w12"}else{fishdata$Week[i]<-"w09"}}
fishdata$TEid<-paste(fishdata$Week,fishdata$Enc2, sep="")

####losing treatments when merge basal and fish - bc basal not complete
#adding basal resources (chlA and AFDM)
head(basalres) #result of Slurry Analysis script
fishbas<-merge(fishdata, basalres, by=c("TEid","Treatment", "Enc2","Week"), all.x=T)

#adding structure stuff ()
head(MusselData) #result of spatialanalysis script
####NEED to ask Thomas for code/length-biomass regressions for mussels

head(EncDV) #result of spatialanalysis script
EncDV$Enc2<-Treat[match(EncDV$ï..Enclosure, Treat$ï..Enclosure),2]
EncDV$Treatment<-Treat[match(EncDV$ï..Enclosure, Treat$ï..Enclosure),3]
EncDV[EncDV$Treatment=="CTRL","Avg.Exposure"]<-0
#following forloop adds velocity collected at week 9
for(k in 1:nrow(fishbas)){
if(fishbas$Week[k]=="w09"){fishbas$Velocity[k]<-EncDV[fishbas[k,"Enc2"],3]}
if(fishbas$Week[k]=="w09"){fishbas$Depth[k]<-EncDV[fishbas[k,"Enc2"],2]}}

fishbas$Structure<-EncDV[match(fishbas$Enc2, EncDV$Enc2),4]
               
#adding invertebrate biomass
head(InvSum) #from Inverts script
AllmodelD<-merge(fishbas,InvSum, by=c("TEid","Treatment","Week"))


####Note, since merge removes rows that lack data in one of the data sets
#this is a conservative data set MEANING if data didn't exist for any of these variables
#the row/observation was removed
#to include all rows from one data set, add all=TRUE in the merge commands

#exporting to csv for Garrett
write.csv(AllmodelD, "preliminarydata.csv")
write.csv(fishbas, "fish&basal2.csv")
