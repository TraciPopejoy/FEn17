library(xlsx)
TR203<-read.xlsx("./Documents/41TR203.xlsx", sheetIndex=1)
head(TR203)
TR203$Length.mm<-round(TR203$Length*25.4,1)
TR203$PSL.mm<-round(TR203$PSL*25.4,1)

TXtax<-read.xlsx("./Documents/41TR203.xlsx",sheetIndex=3)
TR203$Tribe<-TXtax[match(TR203$Taxa,TXtax$Taxa),"Tribe"]

library(tidyverse)
AnU3 <- TR203 %>% group_by(AU,Lot,Taxa) %>% 
  summarize(Tribe=unique(Tribe), n=n(), meanLength=mean(Length.mm, na.rm=T)) 
head(AnU)

library(colorspace)

ggplot(AnU, aes(x = as.character(Level), y = n, fill=Tribe)) + 
  geom_bar(stat="identity") + 
  labs(x="Analytical Unit", y="Non-Repetitive Elements") +theme_bw()+
  facet_grid(~AU, switch="x",space="free", scales="free") 

ggplot(AnU[AnU$Taxa!="fossilized marine shell",], aes(x=AU, y=n, fill=Taxa))+
  geom_bar(stat="identity")+xlim(breaks=paste(seq(1:8)))


### Quality Check ###
reid<- read.xlsx("./Documents/41TR203.xlsx",sheetIndex=5)
reidcom<-reid %>% 
  group_by(Lot,Taxa) %>% summarize(n=n(), name="ReID")%>% spread(Taxa, n)

comtab<-TR203 %>% group_by(Lot,Taxa) %>% 
  summarize(n=n(), name="first") %>% spread(Taxa,n)

CHECK<-rbind(reidcom,comtab[match(reidcom$Lot, comtab$Lot),])

SumTOT <- TR203 %>% group_by(Taxa) %>% 
  summarize(n=n(), RA=n/564*100)
Sum$SiteNRE<-rowSums(Sum, na.rm=T)
Sum[9,]<-colSums(Sum, na.rm=T)
rownames(Sum)[9]<-"SpeciesNRE"
Sum[9,1]<-NA
ggplot(SumTOT, aes(x=Taxa, y=n))+geom_bar(stat="identity")+coord_flip()+theme_bw()

###Habitat Quality ###
UnioOUT<- read.xlsx("./Documents/41TR203.xlsx",sheetIndex=6)
colnames(UnioOUT)<-paste(as.character(UnioOUT[1,]))
library(reshape2)
UnioHAB<-melt(UnioOUT, id="AU")
UnioHAB$Habitat<-"Water-body Type"
UnioHAB$Habitat[49:176]<-"Water Depth (dm)"
UnioHAB$Habitat[177:208]<-"Current Velocity"
UnioHAB$Habitat[209:256]<-"Substrate"
UnioHAB$AU<-as.character(UnioHAB$AU)
ggplot(UnioHAB, aes(x=variable, y=value, color=AU,group=AU))+
  geom_path(size=2)+geom_point(color="black")+
  ylab("Percentage NRE")+
  facet_grid(~Habitat, space="free",scales="free")+theme_bw()
ggplot(UnioHAB, aes(x=variable, y=value))+
  geom_boxplot()+
  ylab("Percentage NRE")+
  facet_grid(~Habitat, space="free",scales="free")+theme_bw()+
  theme(axis.text.x=element_text(angle = 35,color="black", hjust=1))

#### Taphonomy ####
taph<-read.xlsx("./Documents/41TR203.xlsx",sheetIndex = 7)
ggplot(taph, aes(x=Density, y=Sphericity))+
  geom_point(aes(size=RA))+
  geom_text(aes(label=Taxa, hjust="inward", vjust="inward"))+
  scale_size_continuous(name="Relative Abundance")+theme_bw()+
  theme(legend.position="top")

ggplot(taph, aes(x=Estimated.Rank, y=RA))+
  geom_point()+
  geom_text(aes(label=Taxa, hjust="inward", vjust="inward"))+
  ylab("Relative Abundance (%)")+xlab("Estimated Robusticity Rank")+
  theme_bw()

TXtax<-TXtax[order(TXtax$NRE, decreasing=T),]
TXtax$Fact<-factor(TXtax$TableName, levels=TXtax$TableName[c(2,21:3,1)]) 
ggplot(TXtax, aes(x=Fact, y=NRE, reorder_size(NRE)))+
  xlab("Taxa")+
  geom_bar(stat="identity")+coord_flip()+theme_bw()
