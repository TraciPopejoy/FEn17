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
AllmodelD<-merge(fishbas,ModelData, by=c("TEid","Treatment","Week"), all=T)
AllmodelD<-AllmodelD[AllmodelD$Week=="w12",]

#Explore the data & check ANOVA assumptions
#Check for homogeneity
leveneTest(AllmodelD$Prob...Fish, AllmodelD$Treatment, center=median)
#Check for independence between treatment and covariate 
testco1<-aov(depth~Treatment, data=AllmodelD)
summary(testco1)
testco2<-aov(depth~Type, data=ModelData)
summary(testco2)
testco3<-aov(velocity~Type, data=ModelData)
summary(testco3)
#Run the ANCOVA 
mod1 <- lm(Total.Fish~depth+AFDMg.cm2+BMDensity.gpm2+Type, data=AllmodelD)
Anova(mod1, type="III")
adjustedMeans<-effect("Type", mod1)
postHocs<-glht(mod1, linfct=mcp(Type="Tukey"))
summary(postHocs)
#Run contrasts/comparisons 
#Check for homogeneity of regression slopes
hist(residuals(mod1), col="darkgray")
plot(fitted(mod1), residuals(mod1))

AllmodelD$BMDensity.gpm2<-AllmodelD$BMDensity.mgpm2/1000
trophic<-melt(AllmodelD[,c(1:4,9,11,16,39)])

afdm<-ggplot(trophic[trophic$variable=="AFDMg.cm2",], 
       aes(x=Treatment, y=value, color=Treatment))+
  geom_point()+facet_wrap(~variable, shrink=T,scales="free")+
  ylab("AFDM (g/cm2)")+xlab("")+fungraph+theme(legend.position="none")

depthG<-ggplot(trophic[trophic$variable=="Depth",], 
             aes(x=Treatment, y=value, color=Treatment))+
  geom_point()+facet_wrap(~variable, shrink=T,scales="free")+
  ylab("Enclosure depth (m)")+xlab("")+fungraph+theme(legend.position="none")

fish<-ggplot(trophic[trophic$variable=="Total.Fish",], 
             aes(x=Treatment, y=value, color=Treatment))+
  geom_point()+facet_wrap(~variable, shrink=T,scales="free")+
  ylab("N fish")+xlab("")+fungraph+theme(legend.position="none")

bug<-ggplot(trophic[trophic$variable=="BMDensity.gpm2",], 
             aes(x=Treatment, y=value, color=Treatment))+
  geom_point()+facet_wrap(~variable, shrink=T,scales="free")+
  ylab("Macroinv. Biomass g/m2")+xlab("")+fungraph+theme(legend.position="none")

grid.arrange(depthG,afdm,bug,fish, nrow=1)

####Note, since merge removes rows that lack data in one of the data sets
#this is a conservative data set MEANING if data didn't exist for any of these variables
#the row/observation was removed
#to include all rows from one data set, add all=TRUE in the merge commands

#exporting to csv for Garrett
write.csv(AllmodelD, "preliminarydata.csv")
write.csv(fishbas, "fish&basal.csv")

ggplot(test[test$ChlA.mg.m2>0,], aes(x=Treatment, y=ChlA.mg.m2))+geom_boxplot()+facet_wrap(~Week)
