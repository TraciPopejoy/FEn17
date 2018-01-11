#dnova<-read.table(file="clipboard", header=TRUE, sep="\t")
mod1<-lm(log(ChlA_mgM2)~Reach*Treat,data=dnova)

mod4<-lm(ChlA_mgM2~Reach*Treat,data=dnova[!dnova$Site%in%c("LY"),])
summary(mod4)
anova(mod4)

Anova(mod4, type=2, white.adjust='hc3')
Anova(mod4, type=3, white.adjust='hc3')

#assumes equal variance
lsmeans(mod4, list(pairwise ~ Reach | Treat, pairwise ~ Treat | Reach))
#allows unequal variance
library(multcomp)
typing.lsm = lsmeans(typing.lm, pairwise ~ type, glhargs=list())
typing.lsm[[2]]




windows()
boxplot(ChlA_mgM2~Site*Treat*Reach, data=dnova, las=2)



##Test is chl-a controls are different among sites
Mod_con<-lm(ChlA_mgM2~Reach, data=dnova[dnova$Treat%in%"C"&dnova$Site!="LY",])
windows()
boxplot(ChlA_mgM2~Reach*Site,data=dnova[dnova$Treat%in%"C"&dnova$Site!="LY",])
library(car)
leveneTest(log(ChlA_mgM2)~Reach,data=dnova[dnova$Treat%in%"C",])

paired_res<-pairwise.t.test(x=log(dnova$ChlA_mgM2), g=paste0(dnova$Site,":",dnova$Reach), p.adjust.method="none")
paired_res2<-as.data.frame(paired_res$p.value)
paired_res2[seq(1,13,2),seq(1,13,2)]

#####################################################################################
###Analysis and plotting with mean subtracted Chl-a samples
#####################################################################################
#'normalize' each bar by subtracting out the average chl-a control on that bar

dnova_norm<-ddply(dnova, .variables=c("Site","Reach","Bar"), .fun=function(x){data.frame(RR=x$ChlA_mgM2-mean(x$ChlA_mgM2[x$Treat%in%c("C")]), Treat=x$Treat)})

###Conduct anova at the site level after subtracting out the control values for each bar
mod2<-lm(RR~Reach*Treat,data=dnova_norm)
mod3<-aov(RR~Reach*Treat,data=dnova_norm)
mod4<-lm(RR~Reach*Treat*Site,data=dnova_norm[dnova_norm$Site!="LY",])
tuklet<-TukeyHSD(mod3)
library(multcompView)
multcompLetters(tuklet$Reach[,4])
multcompLetters(tuklet$Treat[,4])
a<-multcompLetters(tuklet$`Reach:Treat`[,4])
a2<-a$Letters[c(8,1,2,3,4,5,6,7)]
labs<-data.frame(labs=a2,
                 y=ddply(dnova[dnova$Site!="LY",], .variables=c("Treat", "Reach"), .fun=function(x){max(x$ChlA_mgM2, na.rm=TRUE)+3.5})$V1,
                 x=c(1,2, 4,5, 7,8, 10,11))

windows()
#pdf(file=file.choose(), width=6, height=6 )
boxplot(ChlA_mgM2~Reach*Treat, 
        data=dnova[dnova$Site!="LY",],
        ylab = expression(Chlorophyll-a ~"("*response ~-~ control ~ "("*mu*"g" ~ Chlorophyll-a ~~~ m^{-2} ~"))"),
        at = c(1,2, 4,5, 7,8, 10,11),
        names = rep(c("MR", "NM"),4),
        main="Fall 2016 NDS experiment",
        cex.main=2)
axis(c("C","N","NP","P"), 
     side = 1,
     at = c(1.5, 4.5, 7.5, 10.5),
     tck= 0,
     col=NA,
     line=1.5,
     cex.axis=1.5)
abline(h=0, lty=3, col="gray")
text(x=labs$x, y=labs$y, labels=labs$labs, cex=1.5)
dev.off()







