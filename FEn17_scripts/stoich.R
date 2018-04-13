library(xlsx)
stoich<-read.csv("./FEn17_data/171204_PratherPopejoy_processed_data.csv", stringsAsFactors = F)
stoich$perC<-stoich[,"C_mg"]/stoich[,"Amount_mg"]*100
stoich$perN<-stoich[,"N_mg"]/stoich[,"Amount_mg"]*100
stoich$CtN<-stoich[,"perC"]/stoich[,"perN"]

stoichgraph<-stoich[stoich$Enc!="",]
library(reshape2)
stoichgraph2<-melt(stoichgraph[,-c(1:9)])
library(ggplot2)
ggplot(stoichgraph2, aes(x=Treatment, y=value))+geom_boxplot()+facet_wrap(~variable, scales="free")

phos<-read.csv("./FEn17_data/PeriPhosphorus.csv", stringsAsFactors = F)
PStand<-phos[1:24,]
PstandC<-lm(Abs885~mgPpL, PStand)
summary(PstandC)
ggplot(PStand, aes(x=mgPpL, y=Abs885))+geom_point()+
  geom_smooth(method='lm',formula=y~x)

Pb=-0.17671
Pa=0.62089
phos$mgPpLCAL<-((phos$Abs885-Pb)/Pa)*phos$Dilution
phos$DM<-(phos$mgPpLCAL/phos$DryWeight)*.012
phos$perP<-phos$DM*100
colnames(phos)[5]<-"DryWeight.mg"
phos$Treatment<-Treat[match(phos$Enc,Treat$Enclosure2),3]

ggplot()+geom_point(data=PStand, aes(x=mgPpL, y=Abs885), size=2.5)+
  geom_smooth(data=PStand, aes(x=mgPpL, y=Abs885),method='lm',formula=y~x)+
  geom_point(data=phos[!is.na(phos$Treatment),], size=4, alpha=.3,
             aes(x=mgPpLCAL, y=Abs885, color=Treatment))+
  theme_bw()

ggplot()+geom_boxplot(data=phos, aes(x=Treatment, y=Abs885))+
  geom_boxplot(data=PStand, aes(x=Type, y=Abs885))

OKStoich<-merge(stoich[stoich$Enc!="",-c(1,2)],phos[,-c(6:8)], by="Enc")
OKStoich$CtP<-OKStoich[,"perC"]/OKStoich[,"perP"]
OKStoich$NtP<-OKStoich[,"perN"]/OKStoich[,"perP"]

ggplot(OKStoich, aes(x=Treatment, y=NtP))+geom_boxplot()

OKStoichPer<-OKStoich[,c("Enc","Treatment","perC","perN","perP")]
mperNUT<-melt(OKStoichPer)

ggplot(data = mperNUT, aes(x = Enc, y = value, fill = variable)) + 
  geom_bar(stat="identity") + coord_flip()+
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))+
  facet_grid(Treatment~., space="free", scales="free")
ggplot(mperNUT, aes(x=variable, y=value))+geom_boxplot()+
  facet_wrap(~Treatment, scales="free")



testing<-aov(NtP~Treatment, data=OKStoich)
summary(testing)
plot(testing)
TukeyHSD(testing)
library(lsmeans)
leastm<-lsmeans(testing, "Treatment",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

