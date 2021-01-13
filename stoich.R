library(readr); library(tidyverse)
#bring in N & C data, processed through Allen's lab late 2017
stoich<-read_excel("data/Enclosure_Stoichiometry.xlsx", sheet="AlgalCN") %>% 
  mutate(perC=C_mg/Amount_mg*100,
         perN=N_mg/Amount_mg*100,
         C_mg_mgDW=C_mg/Amount_mg,
         N_mg_mgDW=N_mg/Amount_mg,
         moleC=C_mg_mgDW/12,
         moleN=N_mg_mgDW/14,
         moleCtN=moleC/moleN)
#bring in phosphorus data, collected by TPD early 2018
phos<-read_excel("data/Enclosure_Stoichiometry.xlsx", sheet="AlgalP")
PStand<-phos[1:24,]
PstandC<-lm(Abs885~mgPpL, PStand[PStand$mgPpL!=.2,]) #.2 are outliers
summary(PstandC) #need R2 to be >0.99
#plot of standard curve
ggplot(PStand[PStand$mgPpL!=.2,], aes(x=mgPpL, y=Abs885))+geom_point()+
  geom_smooth(method='lm',formula=y~x)
#save the coefficients for calculation of P in samples
Pb=coefficients(PstandC)[1] 
Pa=coefficients(PstandC)[2]
phos <- phos %>% mutate(mgPpLCAL=((phos$Abs885-Pb)/Pa)*Dilution, #mg P per Liter
                        P_mg_mgDW=(mgPpLCAL/DryWeight)*.012, #mg P per mg DW
                        moleP=P_mg_mgDW/31, #moles of P per DW
                        perP=P_mg_mgDW*100) #percentage of P per DW
colnames(phos)[5]<-"DryWeight.mg"
ggplot(phos[phos$Type.1!="invert",], aes(x=mgPpLCAL, y=Abs885))+geom_point(aes(color=Type.1))
#what do my non-periphyton samples look like (Reach scale)
phos %>% filter(Type.1 =="invert")
#join N, C, and P data and calculate molar ratios per DW
OKStoich<-stoich %>% filter(Identifier2=="sample") %>% left_join(phos) %>%
  select(Identifier1, SampleID, Enc, C_mg_mgDW, N_mg_mgDW,
         perC, perN, moleC, moleN, moleCtN, P_mg_mgDW, moleP, perP) %>%
  filter(Identifier1!="mystery 1" & Identifier1!="mystery 2") %>% 
  mutate(CtP=moleC/moleP, NtP=moleN/moleP) %>% left_join(Treat, by=c("Enc"="Enc2"))

biostoich<-OKStoich %>% left_join(MusBiomass, by=c("Enc"="Enc2")) %>%
  mutate(ACTC=replace_na(ACT,0),
         AMBC=replace_na(AMB,0),
         TreatF2=factor(TreatA, levels=c("ACTL","AMBL","ACTS","AMBS","CTRL")))
names(biostoich)

#### stats ####
ggplot(biostoich, aes(x=TreatF2,y=NtP))+geom_boxplot(aes(fill=Type))+
  ylab("molar N : P ratio")+
  theme_bw()
ggplot(biostoich, aes(x=TreatF2,y=moleCtN))+geom_boxplot(aes(fill=Type))+
  ylab("molar C : N ratio")+
  theme_bw()
ggplot(biostoich, aes(x=TreatF2,y=CtP))+geom_boxplot(aes(fill=Type))+
  ylab("molar C : P ratio")+
  theme_bw()
ggplot(biostoich, aes(x=ACTC,y=NtP))+geom_point(aes(shape=Type, color=TreatA), size=3)+
  xlab("Actinonaias biomass (mg)") + ylab("molar N : P ratio")+
  geom_hline(yintercept=mean(biostoich[biostoich$TreatA=="CTRL","NtP"]))+
  theme_bw()
ggplot(biostoich, aes(x=AMBC,y=NtP))+geom_point(aes(shape=Type, color=TreatA), size=3)+
  xlab("Amblema biomass (mg)") + ylab("molar N : P ratio")+
  geom_hline(yintercept=mean(biostoich[biostoich$TreatA=="CTRL","NtP"]))+
  theme_bw()

# Nitrogen to Phosphorus
library(ggplot2)
ggplot(biostoich, aes(y=NtP))+geom_point(aes(x=ACTC), color="blue")+
  geom_point(aes(x=AMBC), color="red")
NtpAnova<-lm(NtP~TreatA, data=biostoich)
summary(NtpAnova)
plot(NtpAnova)
TukeyHSD(NtpAnova)
library(lsmeans)
leastmNTP<-lsmeans(NtpAnova, "TreatA",adjust="tukey")
cld(leastmNTP, alpha=.05, Letters=letters)
# Carbon to Nitrogen
ggplot(biostoich, aes(y=moleCtN))+geom_point(aes(x=ACTC), color="blue")+
  geom_point(aes(x=AMBC), color="red")
CtnAnova<-aov(moleCtN~TreatA, data=biostoich)
summary(CtnAnova) # marginal
plot(CtnAnova)
TukeyHSD(CtnAnova)
leastmCTN<-lsmeans(CtnAnova, "Treatment.x",adjust="tukey")
cld(leastmCTN, alpha=.05, Letters=letters)
# Carbon to Phosphorus
ggplot(biostoich, aes(y=CtP))+geom_point(aes(x=ACTC), color="blue")+
  geom_point(aes(x=AMBC), color="red")
CtpAnova<-aov(CtP~TreatA, data=biostoich)
summary(CtpAnova)
plot(CtpAnova)
TukeyHSD(CtpAnova)
leastmCTP<-lsmeans(CtpAnova, "Treatment.x",adjust="tukey")
cld(leastmCTP, alpha=.05, Letters=letters)