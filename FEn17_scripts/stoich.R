library(xlsx)
stoich<-read.csv("./FEn17_data/171204_PratherPopejoy_processed_data.csv", stringsAsFactors = F)
stoich$perC<-stoich[,"C_mg"]/stoich[,"Amount_mg"]*100
stoich$perN<-stoich[,"N_mg"]/stoich[,"Amount_mg"]*100
stoich$CtN<-stoich[,"perC"]/stoich[,"perN"]

stoichgraph<-stoich[stoich$X.1!="",]
library(reshape2)
stoichgraph2<-melt(stoichgraph[,-c(1:9)])
library(ggplot2)
ggplot(stoichgraph2, aes(x=X.2, y=value))+geom_boxplot()+facet_wrap(~variable, scales="free")
