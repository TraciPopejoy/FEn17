require(plyr)
#d<-read.table(file="clipboard", header=FALSE, as.is=TRUE, sep="\t")
#dfull<-read.table(file="clipboard", header=FALSE, as.is=TRUE, sep="\t")
d2<-d
d2$V3<-41
dtest<-rbind(d,d2)
#define some summary functions
SE<-function(x){sem<-sd(x)/sqrt(length(x))}

ndssummary<-function(y){
a<-list()
for(i in unique(y$V3)){
  for(j in c("N", "P", "NP")){
    n<-length(a)+1
    a[[n]]<-list(i,j,(outer(y$V5[y$V3==i&y$V4==j],y$V5[y$V3==i&y$V4=="C"],FUN="/")))
  }
}

#decompose list of lists into a list of N, P, NP
df_list <- tapply(a, 
                  sapply(a, `[[`,2), 
                  FUN=function(x) do.call(rbind,lapply(x, `[[`,3)))

data.frame("trt"=c("N","NP","P"), as.data.frame(t(sapply(df_list, function(x){
  list("mean"=mean(unlist(x), na.rm = T),
       "median"=median(unlist(x), na.rm = T),
       "se"=SE(x[is.finite(x)]),
       "sd"=sd(unlist(x), na.rm = T),
       "n"=length(x[is.finite(x)]))}))))
}



ddply(dfull, c("V1","V2","V6"), .fun=ndssummary)



dply(dfull, c("V1","V2"), .fun=function(x){mean(x$V5, na.rm = TRUE)})
