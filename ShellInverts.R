#### Invertebrates on Shells ####

library(readxl)
library(tidyverse)

Treat<-read_excel("./FEn17_data/FEn17OKTreatments.xlsx") 

MusselData<-read_csv("./FEn17_data/MusselBMExpFEn17OK.csv")
names(MusselData)[1]<-"Enclosure"
MDat<-MusselData %>% full_join(Treat) %>% 
  select(Enc2,Genus, TreatA, Type,Spp,L, H,W) %>% 
  mutate(ShellSpecies=recode(Genus, "AMBL"="AMB", "ACT"="ACT"),
         ShellSurArea.mm2=2*(L*H)+2*(L*W)+2*(H*W),
         SamID=paste(Enc2, ShellSpecies, sep=".")) %>%
  group_by(SamID) %>% summarize(TShellSurArea.cm2=sum(ShellSurArea.mm2, na.rm=T)*.01)

shellInv<-read_excel("./FEn17_data/ShellInv.xlsx")
unique(shellInv$Taxa)[order(unique(shellInv$Taxa))] #check taxa for correct spelling
SInv<-shellInv %>% filter(Taxa!="Spider") %>% 
  select(Enc2, ShellSpecies, SamType, Taxa, Length.cm)

SCounts<-SInv %>% group_by(Enc2, ShellSpecies) %>% 
  summarize(Nsam=n(),
            richness=length(unique(Taxa)))%>%
  mutate(SamID=paste(Enc2,ShellSpecies, sep="."))%>% full_join(MDat) %>%
  mutate(InvertDensity.npcm2=Nsam/TShellSurArea.cm2) %>% full_join(Treat)

ggplot(na.omit(SCounts[,c("ShellSpecies","InvertDensity.npcm2","Enc2","TreatA")]),
       aes(x=ShellSpecies, y=InvertDensity.npcm2, group=Enc2, color=TreatA))+
  geom_line(size=2)+geom_point()+facet_wrap(~TreatA)+
  ylab("Number of Invertebrates per shell surface area (cm2)")+
  xlab("Species of Shell sampled")+
  scale_color_manual(name="Cage Treatment", values=CP[c(1,2,5,4)])
ggplot(SCounts, aes(x=Type, y=InvertDensity.npcm2, color=TreatA))+geom_boxplot()

SCountsT<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  mutate(SamID=paste(Enc2,ShellSpecies, sep=".")) %>% summarize(N=n()) %>%
  full_join(SCounts) %>% mutate(RA=N/Nsam*100)

ggplot(na.omit(SCountsT[,c("SamID","Taxa","RA","ShellSpecies")]), aes(x=SamID, y=RA, fill=Taxa))+
  geom_bar(stat="identity")+
  facet_grid(~ShellSpecies, scales="free")+
  fungraph

ShellInvD<- SCountsT %>%full_join(MDat)

TaxaTable<-read_excel("./FEn17_data/TaxaTable.xlsx")
match(unique(shellInv$Taxa), TaxaTable$Taxa)