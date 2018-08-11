#script to bring models in, run last
library(readxl)

fishdata<-read_excel("./FEn17_data/videowithabiotics.xlsx")
colnames(fishdata)[c(3,10)]<-c("GTreat","Log(probFish+1)")
#adding columns that match other data
fishdata[,11:12]<-Treat[match(fishdata$Unit, Treat$Enclosure),c("Enc2","TreatA") ] #Treat found in Slurry Analysis script
fishdata$Week<-rep(NA, nrow(fishdata))
for(i in 1:nrow(fishdata)){
  if(fishdata$Month[i]=="Oct"){fishdata$Week[i]<-"w12"}else{fishdata$Week[i]<-"w09"}}
fishdata$TEid<-paste(fishdata$Week,fishdata$Enc2, sep="")

#adding basal resources (chlA and AFDM)
head(basalres) #result of Slurry Analysis script
fishbas<-merge(fishdata, basalres, by=c("TEid", "Week"), all.x=T)
fishbas<-fishbas[,-15]
names(fishbas)[13]<-"Enc2"
###### add in filter derived ChlA?

#adding mussel structure
head(EncDVw09) #result of spatialanalysis script
fishbas$Structure<-EncDVw09[match(fishbas$Enc2, EncDVw09$Enc2),"Avg.Exposure"]
fishbas$Depth[1:50]<-EncDVw09[match(fishbas$Enc2, EncDVw09$Enc2), "Depth.m"][51:100]
fishbas$Velocity[1:50]<-EncDVw09[match(fishbas$Enc2, EncDVw09$Enc2), "V.mps"][51:100]
               
#adding invertebrate biomass
head(InvSum) #from Inverts script
AllmodelD<-merge(fishbas,InvSum, by=c("TEid","Week"), all=T)
names(MBM)[1]<-"Enc2"
AllDATA<-merge(AllmodelD, MBM, by=c("Enc2"), all=T)
AllDATA<-AllDATA[,-c(14,27:30,33,34)]
write_csv(AllDATA, "./Results/modeldata.csv")
ChlAfil2
#### cleaned up in excel - added filter estimates, removed superflous columns

EncDF<-read_excel("./Results/EncModelTableBEST.xlsx", 
                  col_types=c("text","text","text", "skip","text","text","text","numeric",
                              "numeric","numeric","numeric","numeric","numeric","text","numeric",
                              "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                              "numeric","numeric","numeric","numeric","numeric","numeric","numeric",
                              "numeric","numeric","numeric"))
str(EncDF)
EncDF$BMDensity.gpm2<-EncDF$BMDensity.mgpm2/1000
TreatCOV<-EncDF %>% group_by(Treatment.x,Week) %>% summarize(Depth=round(mean(Depth, na.rm=T),2),
                                                           Velocity=round(mean(Velocity, na.rm=T),3),
                                                           Chl=round(mean(as.numeric(ChlA.ug.cm2), na.rm=T),4), 
                                                           ChlSD=round(sd(ChlA.ug.cm2, na.rm=T),4), 
                                                           AFDM=round(mean(AFDMg.cm2, na.rm=T),4),
                                                           AFDMSD=round(sd(AFDMg.cm2, na.rm=T),4),
                                                           MStructure=round(mean(Structure, na.rm=T),2),
                                                           MStructureSD=round(sd(Structure, na.rm=T),2),
                                                           MSum=round(mean(sumBM, na.rm=T),2),
                                                           MSumSD=round(sd(sumBM, na.rm=T),2))
write.csv(TreatCOV, "./Results/covariate.csv")
MusselT<-EncDF %>% group_by(Treatment.x) %>% summarize(MStructure=round(mean(Structure, na.rm=T),2),
                                                       MStructureSD=round(sd(Structure, na.rm=T),2),
                                                       MSum=round(mean(sumBM, na.rm=T),2),
                                                       MSumSD=round(sd(sumBM, na.rm=T),2),
                                                       ACT=round(mean(ACTbm, na.rm=T),2),
                                                       ACTsd=round(sd(ACTbm, na.rm=T),2),
                                                       AMB=round(mean(AMBbm, na.rm=T),2),
                                                       AMBsd=round(sd(AMBbm, na.rm=T),2),
                                                       meanBM=round(mean(meanBM, na.rm=T),2),
                                                       meanBNsd=round(sd(meanBM, na.rm=T),2))
write.csv(MusselT, "./Results/musselBM.csv")

EncGRAPH<-melt(EncDF[,c(1:6,15,17,32,12)])
ggplot(EncGRAPH,aes(x=variable, y=value, fill=GTreat))+geom_boxplot()+
  fungraph




###### Build GLMM models #####

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
