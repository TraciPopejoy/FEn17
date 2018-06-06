library(lme4)
#setwd("~/Ph.D. stuff/VideoR")
#load video dataset
Viddata<-read.csv(file="./FEn17_data/VideoData.csv", header = TRUE, sep = ",", row.names = c(1))
#Create data frame
Viddata<-(data.frame(Viddata))
#check that headers are right
head(Viddata)
attach(Viddata)
#make treatment and month factors
as.factor(Treatment)
as.factor(Month)
#use lattice to visualize data 
library(lattice)
xyplot(Nsampledetect~Treatment|Month)
xyplot(TotalFishDetected~Treatment|Month)
xyplot(prop.sampl.detect~Treatment|Month)
xyplot(Avg.Fish.Detect~Treatment|Month)
xyplot(richness~Treatment|Month)
head(Viddata)

#bivariate plots of different response metrics.
plot(Depth,Nsampledetect)
plot(Depth,Velocity)
plot(Nsampledetect, TotalFishDetected)
plot(Nsampledetect, Avg.Fish.Detect)
plot(Nsampledetect, prop.sampl.detect)
plot(Nsampledetect, richness)

##glm with poisson distribuiton, use AIC to check support for model with and without treatment
##The response variable is the number of samples (30 sec windows) that a fish was detected.There
##are 24, 30 sec windows possible for each enclosure. Random effect month is just the sampling date. 
glm1<-glmer(Nsampledetect~ (1|Month), family=poisson, data=Viddata);summary(glm1)
glm2<-glmer(Nsampledetect~Treatment+(1|Month), family=poisson, data=Viddata);summary(glm2)
AIC(glm1, glm2)

#load lsmeans to get pairwise comparison for treatment effect while controlling for random effect of sampling period (month)
library(lsmeans)
library(multcompView)
least1<-lsmeans(glm2, list(pairwise~Treatment), adjust ="Tukey")
cld(least1, alpha=.05, Letters=letters)

##Same model instead the response variable is the total fish detected in all 24, 30 second
##time windows for each enclosure. 
glm3<-glmer(TotalFishDetected~ (1|Month), family=poisson, data=Viddata);summary(glm3)
glm4<-glmer(TotalFishDetected~Treatment+(1|Month), family=poisson, data=Viddata);summary(glm4)
AIC(glm3, glm4)
##pairwise comparison of treatment means while controlling for sampling date. 
least2<-lsmeans(glm4, list(pairwise~Treatment), adjust ="Tukey")
cld(least2, alpha=.05, Letters=letters)

##### assumptions
par(mfrow=c(1,1))
plot(fitted(glm4),resid(glm4))
abline(h=0,lty=2,col="red")
qqnorm(resid(glm4))
qqline(resid(glm4))
qqnorm(ranef(glm4))
qqline(ranef(glm4))
qqnorm(ranef(glm4))
qqline(ranef(glm4))
scatter.smooth(fitted(glm4),sqrt(abs(resid(glm4))))

