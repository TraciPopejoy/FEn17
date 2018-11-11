#libraries
library(reshape2);library(plyr);library(dplyr);library(ggplot2);library(xlsx)
library(vegan)

#############     Invert Data in Field Surber Samples    ##############
#bring in enclosure and treatment data, trait/taxa list, and length weight regressions
Treat<-read.xlsx("./FEn17_data/FEn17OKTreatments.xlsx", sheetIndex = 1) 
TaxaList<-read.xlsx("./FEn17_data/TaxaTable.xlsx", sheetIndex = 1)
BiomassReg<-read.xlsx("./FEn17_data/Macroinv Power Law Coeffs TBP.xlsx", sheetIndex = 1)
#THIS IS THE INSECT DATA
FEn17Inv12<-read.csv("./FEn17_data/FEn17InvMeas.csv", stringsAsFactors = F)
FEn17Inv12<-FEn17Inv12[FEn17Inv12$Label!="rulerxocc.tif",-c(3:5)]
FEn17Inv12$TEid<-substring(FEn17Inv12$Label, 6,11)
FEn17Inv12$Enc<-substring(FEn17Inv12$Label,9,11)
FEn17Inv12$Week<-substring(FEn17Inv12$Label,6,8)
FEn17Inv09<-read.csv("./FEn17_data/FEn17w09.csv", stringsAsFactors = F)
FEn17Inv09<-FEn17Inv09[,-c(3:5,9)]
FEn17Inv09$Week<-"w09"
FEn17Inv09$TEid<-paste(FEn17Inv09$Week, FEn17Inv09$Enc, sep="")
FEn17Inv09<-FEn17Inv09[,c(1:4,7,5,6)]
Inv<-rbind(FEn17Inv12, FEn17Inv09) #contains every insect identified from baskets

#clean the data frame
Inv[,8:10]<-Treat[match(Inv$Enc, Treat$Enc2), c("TreatA","Type","Spp")]
sort(unique(Inv$Taxa)) #check to make sure no misspellings
Inv$Family<-as.character(TaxaList$Family[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Order<-as.character(TaxaList$Order[match(Inv$Taxa,TaxaList$Taxa)])
Inv$Length<-Inv$Length.cm*10

ggplot(Inv[Inv$Order=="Odonata",], aes(x=Length.cm))+geom_histogram()
#removing all odonates >1.4cm
remOD<-Inv[(Inv$Order=="Odonata" & Inv$Length.cm>1),]
FEn17ODcount<-remOD %>% group_by(Enc, Week) %>% 
  mutate(biomass=case_when(Taxa=="Od.Dragonflies"~0.0082*((Length.cm*10)^2.813),
                           Taxa=="Od.Damselflies"~0.0086*((Length.cm*10)^2.666))) %>% #Smock 1980
  summarize(totalODbiomass.mg=sum(biomass),
            avg.length.cm=mean(Length.cm),
            nLargeOd=n()) %>% 
  left_join(Treat, by=c("Enc"="Enc2")) %>% select(-TreatF)
ggplot(remOD, aes(x=TreatA, y=Length.cm, color=Week))+
  geom_point(size=2, position=position_dodge(width=.6))+
  stat_summary(color="black")+facet_wrap(~Type, scales="free")+
  theme_bw()
remOD %>% group_by(TreatA, Week) %>% tally()
#Inv<-Inv[!(Inv$Order=="Odonata" & Inv$Length.cm>1.4),]

#how many individuals of each taxa in each enclosure/time; long format
Counts<-ddply(Inv, .variables = c("TEid","Taxa"), .fun=function(x) {
  data.frame(TEid=x[1,5],
            Enc=x[1,6],
            Week=x[1,7],
            Treatment=x[1,8],
            n=count(x,x[1,4])[,2])
  })
Counts$Family<-as.character(TaxaList$Family[match(Counts$Taxa,TaxaList$Taxa)])
Counts$Order<-as.character(TaxaList$Order[match(Counts$Taxa,TaxaList$Taxa)])

#converting to density to compensate for different sampling effort
head(SlurryData) #found in Slurry Analysis sheet
Counts$Density.npm<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5]*.03315)
Counts$Density.npb<-Counts$n/(SlurryData[match(Counts$TEid, SlurryData$TEid),5])

#############     Field Invert Biomass Calculation     #############
#removing taxa that give me trouble 
InvA<-Inv[!Inv$Order=="misc",]
InvB<-na.omit(InvA)
#apply appropriate biomass regressions to each length
InvBMsize<-ddply(.data=InvB, .var=c("Taxa"), .fun=function(x) {
  idx<-x[1,c("Family","Order")] #what family/order are we on
  if(idx$Family=="misc"){ #if not ID'd to family, use order level regressions
    plcoeffs<-BiomassReg[BiomassReg$Order == idx$Order &
                         !is.na(BiomassReg$Order == idx$Order),]
  }else{ #pull all the regressions for that family
    plcoeffs<-BiomassReg[BiomassReg$Family==idx$Family&
                         BiomassReg$Order == idx$Order&
                         !is.na(BiomassReg$Family==idx$Family), ]  
  }
  #which regressions were actually built for insects this size
  ldply(lapply(1:dim(x)[[1]],FUN=function(i) { 
    #idx2<-c(x$Length[i]>=plcoeffs[,19] & x$Length[i]<=plcoeffs[,20]) 
    idx2<-T #if no size range listed, use the regression anyways
    d1<-plcoeffs[idx2,]
    indmassest<-d1$a*(x$Length[i]^d1$b) #power law to determine biomass
    data.frame(
      TEid=x$TEid[i],
      Enc=x$Enc[i],
      Week=x$Week[i],
      Treatment=x$TreatA[i],
      length=x$Length[i],
      neq=length(idx2), #number of possible equations used
      ninR=sum(idx2), #number of equations used
      meanBM.mg=mean(indmassest),
      median=median(indmassest),
      stDev=sd(indmassest))}), 
    data.frame)
})
#get mean biomass and sum of each taxa for each enclosure/time
InvTotalBM<-ddply(InvBM, .var=c("TEid", "Taxa"),
                  .fun=function(x){data.frame(mean.mg=mean(x$meanBM.mg, na.rm=T), 
                                              median.mg=median(x$meanBM.mg, na.rm=T), 
                                              SD.mg=sd(x$meanBM.mg, na.rm=T),
                                              Sum.mg=sum(x$meanBM.mg,na.rm=T),
                                              Treatment=x$Treatment[1],
                                              Enc=x$Enc[1],
                                              Week=x$Week[1],
                                              mean.length=mean(x$length, na.rm=T)
                                              )}) 
#check it worked
#ggplot(na.omit(InvTotalBM), aes(x=Taxa, y=Sum.mg, color=Treatment))+
#  geom_point()+coord_flip()+scale_y_log10()

#converting mean biomass to biomass/meter
InvTotalBM$BMDensity<-InvTotalBM$Sum.mg/(SlurryData[match(InvTotalBM$TEid, SlurryData$TEid),"Basket."]*.03315)
###would use density because sampling was not constant (not always full basket recovery)

#get each taxa and teid with counts and biomass, and density of both
InvGraph<-merge(InvTotalBM[,-c(4,5)],Counts, by=c("TEid","Taxa","Treatment","Enc","Week"))
#get traits into the graph
InvGraph<-merge(InvGraph, TaxaList[,c(1,2,5:25)], by=c("Taxa", "Family"))
InvGraph[,36:37]<-Treat[match(InvGraph$Enc, Treat$Enc2), c("Type","Spp")]
TraitDef<-read.xlsx("./FEn17_data/TaxaTable.xlsx",sheetIndex = 2)
trophic<-TraitDef[TraitDef$Trait=="T.TropP",]
InvGraph$FFGp<-trophic[match(InvGraph$T.TropP, trophic$Num),"T.state"]

#### graph city ####
fungraph<-theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1),
                axis.text.y = element_text(size=12,color="black"),
                axis.title.y=element_text(size=20),
                plot.background = element_blank(),
                panel.border=element_blank(),
                panel.grid.major= element_line(colour=NA), 
                panel.grid.minor=element_line(colour=NA),
                title=element_text(size=20),
                panel.background = element_rect(fill = "white"),
                legend.key=element_rect(colour=NA), 
                axis.line.x=element_line(colour="black"),
                axis.line.y=element_line(colour="black"),
                strip.background=element_rect(fill="white", color="black"),
                strip.text=element_text(size=15))
library(colorspace)
CP<-diverge_hcl(5, h=c(180,70), c = 100, l = c(50, 90), power = 1)
CP[3]<-"black"
CP2<-data.frame(colorss=CP[c(1,2,5,4)], Treat=unique(Treat$TreatA))

Traitplot<-InvGraph %>% group_by(Treatment,Family,FFGp) %>% summarize(meanDenMeter=mean(Density.npm),
                                                                      meansize=mean(mean.length),
                                                                      meanDenBask=mean(Density.npb))
Traitplot$Type<-Treat[match(Traitplot$Treatment, Treat$Treatment), "Type"]

file<-c("1.tiff","2.tiff","3.tiff","4.tiff","5.tiff")

for(i in 1:length(unique(Traitplot$Treatment))){
  ggplot(Traitplot, 
         aes(x=FFGp, y=meansize, color=Treatment))+
    geom_point(data=subset(Traitplot, 
                           Treatment==c(paste(unique(Traitplot$Treatment)[i]))),
      aes(size=sqrt(meanDenMeter)), position=position_dodge(width=.1))+
    ylim(c(0.5,5.25))+
    labs(y="Mean Size (mm)", x="Functional Feeding Group")+
    scale_color_manual(values=paste(CP2[CP2$Treat==paste(unique(Traitplot$Treatment)[i]),"colorss"]),
                       name="Treatment")+
    scale_size_area(name=expression(sqrt(Individuals/m^{2}))) +
    theme_bw()+theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1),
                     axis.text.y = element_text(size=12,color="black"),
                     axis.title.y=element_text(size=20),
                     title=element_text(size=20),
                     panel.background = element_rect(fill = "white"),
                     legend.key=element_rect(colour=NA), 
                     axis.line.x=element_line(colour="black"),
                     axis.line.y=element_line(colour="black"),
                     strip.background=element_rect(fill="white", color="black"),
                     strip.text=element_text(size=15))
  ggsave(file[i])
}

ggplot(Traitplot, 
       aes(x=FFG, y=meansize, color=Treatment))+
  geom_point(aes(size=sqrt(meanDenMeter)), alpha=.6)+
  ylim(c(.5,5.25))+
  labs(y="Mean Size (mm)", x="Functional Feeding Group")+
  scale_color_manual(breaks=c("CTRL","ACTL","ACTS","AMBL","AMBS"),
                     labels=c("Control","ACT live","ACT sham","AMB live","AMB sham"),
                     values=CP[c(3,1,2,5,4)],
                     name="Treatment")+
  scale_size_area(name=expression(sqrt(Individuals/m^{2}))) +
  theme_bw()+theme(axis.text.x=element_text(angle = 35,size=12,color="black", hjust=1),
                   axis.text.y = element_text(size=12,color="black"),
                   axis.title.y=element_text(size=20),
                   title=element_text(size=20),
                   panel.background = element_rect(fill = "white"),
                   legend.key=element_rect(colour=NA), 
                   axis.line.x=element_line(colour="black"),
                   axis.line.y=element_line(colour="black"),
                   strip.background=element_rect(fill="white", color="black"),
                   strip.text=element_text(size=15))
ggsave("sfsplot.png")

ggplot(Traitplot[Traitplot$FFG=="Herbivore",], 
       aes(x=Family, y=meansize, color=Treatment, label=FFG))+
  geom_point(aes(size=sqrt(meanDenBask)), shape=21)+
  labs(y="Mean Size (mm)", x="Functional Feeding Group")+
  scale_color_manual(values=CP[c(3,1,2,5,4)],
                     name="Treatment")+
  scale_size_area(name=expression(individuals/basket))+
  fungraph + facet_grid(~FFG, switch="x", scales="free", space="free")

ggplot(InvGraph[InvGraph$Family=="Heptageniidae" | 
                   InvGraph$Family=="Polycentropidae"|
                   InvGraph$Family=="Chironomidae",], 
       aes(x=Family, y=Density.npm, fill=Treatment))+
  geom_boxplot()+
  labs(y=expression(Individuals/m^{2}), x="Functional Feeding Group")+
  scale_fill_manual(values=CP[c(3,1,2,5,4)],
                     name="Treatment")+
  theme_bw() +  theme(axis.text.x=element_text(angle = 0,size=12,color="black", 
                                               hjust=.5),
                      axis.text.y = element_text(size=12,color="black"),
                      axis.title.y=element_text(size=20),
                      title=element_text(size=20),
                      panel.background = element_rect(fill = "white"),
                      legend.key=element_rect(colour=NA), 
                      axis.line.x=element_line(colour="black"),
                      axis.line.y=element_line(colour="black"),
                      strip.background=element_rect(fill="white", color="black"),
                      strip.text=element_text(size=15))+
  facet_grid(~FFG+Week, switch="x", scales="free", space="free")
ggsave("aund.png")
ggplot(remOD, aes(x=Family, y=Length.cm, color=Treatment))+
  geom_point(size=3, alpha=.5, position=position_dodge(width=.2))+
  labs(y="Mean Size (cm)", x="Functional Feeding Group")+
  scale_color_manual(values=CP[c(3,1,2,5,4)],
                     name="Treatment")+fungraph
 

ggplot(Traitplot[Traitplot$FFG=="Herbivore" | 
                   Traitplot$FFG=="Predator",], 
       aes(x=FFG, y=meansize, color=Treatment))+
  geom_point(aes(size=meanDenMeter))+
  labs(y="Mean Size (mm)", x="Functional Feeding Group")+
  scale_color_manual(values=CP[c(3,1,2,5,4)],
                     name="Treatment")+
  scale_size_area(name=expression(individuals/m^{2})) +
  fungraph+facet_wrap(~Treatment, ncol=5)
ggsave("ComboTrait.tiff")

ggplot(Traitplot[Traitplot$FFG=="Herbivore" | 
                   Traitplot$FFG=="Predator",|
                   Traitplot$FFG=="C-Filterer"], 
       aes(x=Family, y=meansize, color=Treatment))+
  geom_point(aes(size=meanDenMeter))+
  labs(y="Mean Size (mm)", x="Functional Feeding Group")+
  scale_color_manual(values=CP[c(3,1,2,5,4)],
                     name="Treatment")+
  scale_size_area(name=expression(individuals/m^{2})) +
  fungraph+facet_wrap(~Treatment, ncol=5)
ggsave("ComboTrait.tiff")

ggplot(InvGraph, 
       aes(x=FFG, y=mean.length, color=Treatment))+
  scale_y_sqrt() +scale_color_manual(values=CP[c(3,1,2,5,4)])+
  scale_size_area(name=expression(individuals/m^{2}))+
  geom_point(aes(size=Density.npm), alpha=.5, position="jitter")+
  ylab("Mean Length (mm)")+xlab("Functional Feeding Group")+
  facet_grid(.~Type, space="free", scales = "free")+fungraph
ggsave("./Figures/EnclosureTraitSpace.tiff")

ggplot(InvGraph, 
       aes(x=Order, y=Density.npm, fill=FFG))+
  scale_y_log10() +coord_flip()+
  geom_boxplot()+
  facet_wrap(~Type)+theme_classic()

ggplot(na.omit(InvGraph), 
       aes(x=Order, y=Density.npm, fill=Order))+
  scale_y_log10() + coord_flip() +
  geom_boxplot()+
  facet_wrap(~Type)+theme_classic()

sizegraph<-merge(InvB, TaxaList)
sizegraph$Type<-InvGraph[match(sizegraph$Treatment, InvGraph$Treatment), "Type"]
ggplot(sizegraph[sizegraph$T.Trop==2,], aes(x=Treatment, y=Length))+
  geom_violin()+facet_wrap(~Type, scales="free")

ggplot(InvGraph, 
       aes(x=Type, y=mean.length))+
  scale_y_sqrt() + geom_violin(aes(fill=Treatment))+
  geom_boxplot(alpha=.3, aes(group=Treatment))+
  #geom_point(aes(size=Density.npb), alpha=.3)+
  ylab("Mean Length (mm)")+xlab("Functional Feeding Group")+
  facet_grid(.~FFG, space="free", scales = "free")+fungraph+
  scale_fill_brewer(palette="Paired")

ggplot(InvGraph, aes(x=FFG))

(p1<-ggplot(InvGraph[InvGraph$T.Trop==2,], 
       aes(x=Type, y=Density.npm, fill=Treatment))+
  geom_boxplot()+
  ylab("Density (n/meter.sq)")+xlab("Treatments")+ylim(c(0,160))+
  facet_grid(.~FFG, space="free", scales = "free")+
  scale_fill_brewer(palette="Paired")+fungraph)

(p2<-ggplot(InvGraph[InvGraph$T.Trop==2,], 
       aes(x=Type, y=mean.length, fill=Treatment))+
  geom_boxplot()+
  ylab("Mean Length (mm)")+xlab("Treatments")+
  facet_grid(.~FFG, space="free", scales = "free")+
  scale_fill_brewer(palette="Paired")+fungraph)
library(gridExtra)
grid.arrange(p1,p2,ncol=2)
ggsave("./Figures/Filterer.tiff")

###Total Summary of Data####

InvSumA<-ddply(Counts,.variables=c("TEid"),.fun=function(x) {count(x,x[1,1])[,2]})
#make a typical community matrix (col:species, row: obs)
commat<-dcast(Counts[,-c(3,4,5,6,7,8,9)], TEid~...)
commat[is.na(commat)]<-0
rownames(commat)<-commat[,1]
commat<-commat[,-1]
InvSumA$Shannon<-diversity(commat) #calculate Shannon diversity index
InvSumA$TotalN<-rowSums(commat) #calculate total, density standardized insects
#sum all the biomass measures for each TEid
InvSumB<-ddply(InvTotalBM,.variables=c("TEid"),.fun=function(x) data.frame(TotalBM.mg=sum(x$Sum.mg),
                                                                           BMDensity.mgpm2=sum(x$BMDensity)))
InvSum<-merge(InvSumA,InvSumB, by="TEid")
colnames(InvSum)[2]<-"richness"
InvSum[,7:11]<-Inv[match(InvSum$TEid, Inv$TEid), c("Enc","Week","TreatA","Type","Spp")]
InvSum$basketn<-SlurryData[match(InvSum$TEid, SlurryData$TEid),5]
InvSum$Density.npm<-InvSum$TotalN/(InvSum$basketn*.03315)
InvSum$Enclosure<-Treat[match(InvSum$Enc, Treat$Enc2),1]

####testing CLEANME ####
test<-InvGraph %>% group_by(TEid)%>%filter(T.TropP==2) %>%
  mutate(CFiltDen=sum(Density.npm))
InvSum$CFiltDen<-test[match(InvSum$TEid, test$TEid), "CFiltDen"]
InvSum[is.na(InvSum$CFiltDen),"CFiltDen"]<-0
InvSum$CFiltDen<-unlist(InvSum$CFiltDen)
mathss<-as.data.frame(InvSum)

testing<-aov(CFiltDen~Type, data=mathss)
summary(testing)
plot(testing)
TukeyHSD(testing)

library(lsmeans)
leastm<-lsmeans(testing, "Type",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

test<-InvGraph %>% group_by(TEid)%>%filter(T.TropP==2) %>%
  mutate(CFiltBM=mean(mean.length))
InvSum$CFiltBM<-test[match(InvSum$TEid, test$TEid), "CFiltBM"]
InvSum[is.na(InvSum$CFiltBM),"CFiltBM"]<-0
InvSum$CFiltBM<-unlist(InvSum$CFiltBM)
mathss2<-as.data.frame(InvSum)

testing<-aov(CFiltBM~Type, data=mathss2)
summary(testing)
plot(testing)
TukeyHSD(testing)
leastm<-lsmeans(testing, "Type",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)


write.csv(InvSum, "FEn17week12insects.csv")

RAcommat<-commat/rowSums(commat)
RAcommat$sites<-rownames(RAcommat)
RAcommat$Type<-InvSum[match(RAcommat$sites, InvSum$TEid),12]
ComGraph<-melt(RAcommat)
ComGraph$Order<-InvGraph[match(ComGraph$variable, InvGraph$Taxa),12]
ComGraph$Trop<-TaxaList[match(ComGraph$variable, TaxaList$Taxa),"T.Trop"]

commat1<-commat
commat1$sites<-rownames(commat1)
commat1$Type<-InvSum[match(commat1$sites, InvSum$TEid),12]
commat2<-melt(commat1)
commat2$Trop<-TaxaList[match(commat2$variable, TaxaList$Taxa),"T.TropP"]

ggplot(InvSum, aes(x=TreatA, y=BMDensity.mgpm2, color=TreatA))+
  geom_point(cex=5)+ylim(0,max(InvSum$BMDensity.mgpm2))+
  ylab(expression(Biomass~Density~~mg/m^{2}))+
  scale_color_manual(values=CP[c(1,2,5,4,3)])+
  stat_summary(color="black")+
  fungraph
ggsave("./Figures/BMDensity.tiff")

#####NMDS analysis####
#plotting advice from Christopher Chizinski GitHub
library(vegan)
nmds<-metaMDS(commat[-51,-1])
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- Inv[match(data.scores$site, Inv$TEid), "TreatA"]  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

grp.a <- data.scores[data.scores$grp == "CTRL", ][chull(data.scores[data.scores$grp == 
                                                                   "CTRL", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "ACTL", ][chull(data.scores[data.scores$grp == 
                                                                   "ACTL", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- data.scores[data.scores$grp == "ACTS", ][chull(data.scores[data.scores$grp == 
                                                                    "ACTS", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.d <- data.scores[data.scores$grp == "AMBL", ][chull(data.scores[data.scores$grp == 
                                                                    "AMBL", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.e <- data.scores[data.scores$grp == "AMBS", ][chull(data.scores[data.scores$grp == 
                                                                    "AMBS", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
hull.data <- rbind(grp.a, grp.b, grp.c, grp.d, grp.e)  #combine grp.a and grp.b

ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=2) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=4) + # add the point markers
  scale_colour_manual(values=CP[c(3,1,2,5,4)]) +
  scale_fill_manual(values=CP[c(3,1,2,5,4)])+
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
ggsave("./Figures/NMDSwhole.png")

colSums(commat[-51,-1])
ComThin<-commat[-51,colSums(commat[-51,])>5]

nmds<-metaMDS(ComThin)
data.scores <- as.data.frame(scores(nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- Inv[match(data.scores$site, Inv$TEid), "TreatA"]  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

grp.a <- data.scores[data.scores$grp == "CTRL", ][chull(data.scores[data.scores$grp == 
                                                                      "CTRL", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "ACTL", ][chull(data.scores[data.scores$grp == 
                                                                      "ACTL", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.c <- data.scores[data.scores$grp == "ACTS", ][chull(data.scores[data.scores$grp == 
                                                                      "ACTS", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.d <- data.scores[data.scores$grp == "AMBL", ][chull(data.scores[data.scores$grp == 
                                                                      "AMBL", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.e <- data.scores[data.scores$grp == "AMBS", ][chull(data.scores[data.scores$grp == 
                                                                      "AMBS", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
hull.data <- rbind(grp.a, grp.b, grp.c, grp.d, grp.e)  #combine grp.a and grp.b
hull.data

ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=2) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=4) + # add the point markers
  scale_colour_manual(values=CP[c(3,1,2,5,4)]) +
  scale_fill_manual(values=CP[c(3,1,2,5,4)])+
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
ggsave("./Figures/NmdsThin.png")

head(EncDF) #result of Model Script, has environmental varialbes in there
pts<-envfit(nmds,EncDF[EncDF$TEid!="w09E02",c(7,8,14,15,27)], na.rm=F) #not working 
pts.df<-as.data.frame(pts$vectors$arrows*sqrt(pts$vectors$r))
pts.df$species<-rownames(pts.df)

ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=2) +  # add the species labels
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=2) + # add the point markers
  geom_segment(data=pts.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey") + 
  geom_text(data=pts.df,aes(x=NMDS1,y=NMDS2,label=species),size=5) +
  scale_colour_manual(values=CP[c(3,1,2,5,4)]) +
  scale_fill_manual(values=CP[c(3,1,2,5,4)])+
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
ggsave("./Figures/nmdschar.png")

##### Functional Diversity Analysis #####
invtrait<-read.csv("./FEn17_data/InvertTraitsTable_v1.txt", sep="\t")

Invmisc<-Inv[Inv$Family!="misc",]
oldtraits<-NULL
for(j in 1:length(unique(Invmisc$Family))){
  k<-unique(Invmisc$Family)[!is.na(unique(Invmisc$Family))]
  idx<-k[j] #what family/order are we on
  traits<-invtrait[invtrait$Family==idx & 
                   !is.na(invtrait$Family==idx), ] 
  oldtraits<-rbind(oldtraits,traits)
}
which(duplicated(oldtraits))#checking for duplicate rows within the trait matrix

FeedTraits<-oldtraits[oldtraits$Feed_mode_prim!="" &
                      !is.na(oldtraits$Feed_mode_prim==""),
                      c("Family","Feed_mode_prim","Feed_mode_sec","Feed_mode_comments",
                        "TraitRecord_ID")]
library(tidyverse)
library(dplyr)
getmode<- function(v) {
  uniqv<-unique(v)
  uniqv[which.max(tabulate(match(v,uniqv)))]
}

FFGtrait<-FeedTraits %>% group_by(Family) %>% summarize(n=n(),
                                              primFFG=getmode(Feed_mode_prim),
                                              secFFG=getmode(Feed_mode_sec))
### replace Others with consensus
### replace factors with Leroy Poff numbers
### look at Schnider for new method
### explicitly looking at differences between scale. 


for(j in 1:nrow(invtrait)){
  if(!is.na(match(invtrait$Family[1:5], Inv$Family))==T){print(invtrait$Family[j])}}

test<-invtrait[match(invtrait$Family, unique(Inv$Family)),]
test2<-test[!is.na(test$Family),]
nrow(test2)

for(j in 1:nrow(invtrait)){if(match(invtrait$Family[j], Inv$Family)==T){print(invtrait$Family[j])}}

trait<-TaxaList[,-c(1:4,26,27,28)]
rownames(trait)<-TaxaList[,1]
ordtrait<-trait[order(rownames(trait),decreasing=F),]
trait1<-ordtrait[-33,-21]

AbMatrixT<-dcast(Counts, TEid + Treatment ~ Taxa, value.var = "Density.npm")
AbMatrixT[is.na(AbMatrixT)]<-0
Abundances<-AbMatrixT[,-c(1:2)]
rownames(Abundances)<-AbMatrixT[,1]
AbundR<-Abundances[-c(19,49)]

rownames(trait1)==colnames(AbundR) #need it to be true to run the function

library(FD)
ex <- dbFD(trait1,AbundR)

FunciGraph<-data.frame(FDis=ex$FDis,
                       FRich=ex$FRic,
                       FEven=ex$FEve,
                       RaoQ=ex$RaoQ)
FunciGraph$TEid<-rownames(FunciGraph)
FunciGraph$Treatment<-InvGraph[match(FunciGraph$TEid, InvGraph$TEid),"Treatment"]
FunciGraph[,c("Type","Spp")]<-treattype[match(FunciGraph$Treatment, treattype$Treatment),c("Type","Spp")]


fit <- aov(FDis ~ Type, data=FunciGraph)
summary(fit)
plot(fit)
TukeyHSD(fit)
leastm<-lsmeans(fit, "Type",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

mFGraph<-melt(FunciGraph)

ggplot(mFGraph, aes(x=Type, y=value, fill=Treatment))+geom_boxplot()+
  facet_wrap(~variable, scales = "free")+
  scale_fill_manual(values)+fungraph
ggsave("./Figures/AllFunctionInd.tiff")

ggplot(mFGraph[mFGraph$variable=="FDis",], 
       aes(x=Type, y=value, fill=Treatment))+geom_boxplot()+
  scale_fill_manual(values=CP[c(3,1,2,5,4)])+
  labs(y="Functional Distance")+theme_bw() +  
  theme(axis.text.x=element_text(angle = 0,size=12,color="black"), 
              axis.text.y = element_text(size=12,color="black"),
              axis.title.y=element_text(size=20),
              title=element_text(size=20),
              panel.background = element_rect(fill = "white"),
              legend.key=element_rect(colour=NA), 
              axis.line.x=element_line(colour="black"),
              axis.line.y=element_line(colour="black"),
              strip.background=element_rect(fill="white", color="black"),
              strip.text=element_text(size=15))
ggsave("./Figures/FunctionalDis.tiff")


ModelData<-merge(InvSum, FunciGraph, by=c("TEid","Treatment","Type","Spp"))
ModelData$depth<-EncDV[match(ModelData$Enclosure, EncDV$ï..Enclosure), "Depth.m"]
ModelData$velocity<-EncDV[match(ModelData$Enclosure, EncDV$ï..Enclosure), "V.mps"]
ModelData$Treatment<-factor(ModelData$Treatment, 
                            levels=c("ACTL","ACTS","AMBL","AMBS","CTRL"))
#ModelData$Treatment<-relevel(ModelData$Treatment, ref="CTRL")

library(car)
library(compute.es)
library(effects)
library(multcomp) 
#Explore the data & check ANOVA assumptions

#Check for homogeneity
leveneTest(ModelData$FDis, ModelData$Treatment, center=median)

#Check for independence between treatment and covariate 
testco1<-aov(depth~Treatment, data=ModelData)
summary(testco1)
testco2<-aov(depth~Type, data=ModelData)
summary(testco2)
testco3<-aov(velocity~Type, data=ModelData)
summary(testco3)
#Run the ANCOVA 
mod8 <- lm(FDis~depth+velocity+Type, data=ModelData)
Anova(mod8, type="III")
adjustedMeans<-effect("Type", mod8)
postHocs<-glht(mod8, linfct=mcp(Type="Tukey"))
summary(postHocs)
#Run contrasts/comparisons 
#Check for homogeneity of regression slopes
hist(residuals(mod8), col="darkgray")
plot(fitted(mod8), residuals(mod8))

####MacroAbund####
leveneTest(ModelData$Density.npm, ModelData$Type, center=median)
#Check for independence between treatment and covariate 

#Run the ANCOVA 
mod9 <- lm(Density.npm~depth+Type, data=ModelData)
Anova(mod9, type="III")
adjustedMeans<-effect("Type", mod9)
postHocs<-glht(mod9, linfct=mcp(Type="Tukey"))
summary(postHocs)




install.packages("DecomposingFD")
library("DecomposingFD")

