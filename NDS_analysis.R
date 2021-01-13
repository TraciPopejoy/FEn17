library(plyr); library(readr)

nds1<-read_excel("data/Enclosure_Stoichiometry.xlsx", sheet="NDSrawdata")
  read.xlsx("./FEn17_data/Field Encl Stoich.xlsx", sheetIndex = 3, stringAsFactors=F)
nds<-nds1[-c(167,168),c(2,3,5,6,14)] #remove extra rows and columns
nds<-nds[,c(2:5,1)]#longform bar, treatment, chlorophyll
nd<-nds[nds$Bar.!="047",] #had to remove this because only had two C discs
  
#define some summary functions
SE<-function(x){sem<-sd(x)/sqrt(length(x))}

trtlist<-unique(nd[,c(1:4)]) #site,mussel reach, treatment
ls2<-trtlist[trtlist$Treat!="C",]

a<-list() #start with an empty list
#creates a list that isolates each type of disc (bar,treat) and then gives outer table
for(i in unique(na.exclude(nd$Bar.))){
    for(j in c("N", "P", "NP")){
n<-length(a)+1
a[[n]]<-list(i,j,(outer(nd$Chl.A..mg.m2.[nd$Bar.==i & nd$Treat==j],
                        nd$Chl.A..mg.m2.[nd$Bar.==i & nd$Treat=="C"],FUN="/")))
    }
  }

#decompose list of lists into a list of N, P, NP
df_list <- tapply(a, 
                  sapply(a,`[[`,2), FUN=function(x) do.call(rbind,lapply(a,`[[`,3)))

sapply(df_list, function(x){
  list("mean"=mean(unlist(x), na.rm = T),
       "median"=median(unlist(x), na.rm = T),
       "se"=SE(x[is.finite(x)]),
       "sd"=sd(unlist(x), na.rm = T),
       "n"=length(x[is.finite(x)]))})

q<-ldply(a, .fun=function(x){data.frame(barID=x[[1]],
                                        trt=x[[2]],
                                        as.data.frame(sapply(x[[3]],unlist)))})
dnova1<-ddply(q, .variables=c("barID", "trt"),
             .fun=function(x){data.frame(mean=mean(x[,3], na.rm = T),
                                         median=median(x[,3], na.rm = T),
                                         se=SE(x[is.finite(x[,3]),3]),
                                         sd=sd(x[,3], na.rm = T),
                                         n=length(x[is.finite(x[,3]),3]))})
dnova1$Type<-trtlist[match(dnova1$barID, trtlist$Bar.),2]
dnova1$Week<-trtlist[match(dnova1$barID, trtlist$Bar.),1]

mod1<-lm(log(mean)~Type*trt,data=dnova1)

mod4<-lm(mean~Type*trt,data=dnova1[!dnova1$Week==4,])
summary(mod4)
anova(mod4)


#####################################################################################
###Analysis with mean subtracted Chl-a samples
#####################################################################################
#'normalize' each bar by subtracting out the average chl-a control on that bar

dnova_norm<-ddply(nds, .variables=c("Week","Type","Bar."), 
                  .fun=function(x){data.frame(RR=x$Chl.A..mg.m2.-mean(x$Chl.A..mg.m2.[x$Treat=="C"]), Treat=x$Treat)})

###Conduct anova at the site level after subtracting out the control values for each bar
mod2<-lm(RR~Bar.*Treat,data=dnova_norm[dnova_norm$Week!=4 & dnova_norm$Type=="AL",])
mod<-lm(Chl.A..mg.m2.~Type*Treat, data=nds[nds$Week!=4,])
mod3<-aov(RR~Type*Treat,data=dnova_norm[dnova_norm$Week!=4,])
mod4<-lm(RR~Type*Treat*Week,data=dnova_norm)
tuklet<-TukeyHSD(mod3)

HSm1<-lm(RR~Treat, data=dnova_norm[dnova_norm$Week!=4 & dnova_norm$Type=="HS",])
ALm1<-lm(RR~Treat, data=dnova_norm[dnova_norm$Week!=4 & dnova_norm$Type=="AL",])

library(car)
Anova(mod2,type=2)
summary(mod2)

hist(residuals(mod), col="darkgray") #normal distribution assumption
plot(fitted(mod), residuals(mod2)) #homoscedastic assumption

library(lsmeans)
leastm<-lsmeans(mod2, "Treat",adjust="tukey")
cld(leastm, alpha=.05, Letters=letters)

#### plots ####
library(ggplot2)
ggplot(nds[nds$Week!=4,], aes(x=Treat, y=Chl.A..mg.m2.)) + geom_boxplot() +facet_wrap(~Type)
