
---
  title: "Organic Contaminants in Great Lakes Tributaries: Watersheds and Chemicals of Greatest Concern"
author: "<b>USGS</b>: Steve Corsi, Laura DeCicco, Austin Baldwin, David Alvarez, Peter Lenaker <br/><b>North Carolina State University</b>: David Reif"
date: "Wednesday, May 27, 2015"
output:
  ioslides_presentation:
  smaller: yes
css: styles.css    
---
  
## Passive Sampler Regression Tree Model: 
```{r echo=FALSE}
library(rpart)
library(party)
library(partykit)

#################
dat <- mapDataPassive
mainTitle <- "Passive Samplers"
response <- "nChem"

set.seed(seed = 10)
form <- formula(paste(response,"~","Urban.pct + Ag.pct"))
m1 <- rpart(form,control=rpart.control(xval=10, cp=0.01,minBucket=3, minsplit = 10),data=dat)
#plot(as.party(m1), tp_args=list(id=FALSE),main=mainTitle,font=2)

print(m1$cptable)

m2 <- prune(m1,cp=m1$cptable[which.min(m1$cptable[,"xerror"])])
plot(as.party(m2), tp_args=list(id=FALSE))


#summary(m2)


#m <- lm(nChem ~ Urban.pct + Ag.pct,data=dat)
#summary(m)

#plot(dat$nChem,predict(m,newdata=dat))

##################
```


```{r echo=FALSE}
#################
dat <- mapData
mainTitle <- "Water Samples"

set.seed(seed = 10)
m1 <- rpart(nChem ~ Urban.pct + Ag.pct,control=rpart.control(xval=10,cp = 0.01, minBucket=3),data=dat)
plot(as.party(m1), tp_args=list(id=FALSE),main=mainTitle,font=2)

print(m1$cptable)

m2 <- prune(m1,cp=m1$cptable[which.min(m1$cptable[,"xerror"])])
plot(as.party(m2), tp_args=list(id=FALSE))


summary(m2)


m <- lm(nChem ~ Urban.pct + Ag.pct,data=dat)
summary(m)

plot(dat$nChem,predict(m,newdata=dat))

##################
```

```{r echo=FALSE}
mapDataBM <- left_join(stationINFO, summ2_BM, by=c("shortName"="site"))
mapDataBM <- mapDataBM[!is.na(mapDataBM$dec.lat.va),]
mapDataBM <- mapDataBM[!is.na(mapDataBM$maxEAR),]
dat <- mapDataBM
mainTitle <- "Water Quality Guidelines"

response <- "nChem"


set.seed(seed = 13)
form <- formula(paste(response,"~","Urban.pct + Ag.pct"))
m1 <- rpart(form,control=rpart.control(xval=5, cp=0.01,minBucket=10,minsplit = 15,),data=dat)
plot(as.party(m1), tp_args=list(id=FALSE),main=mainTitle,font=2)

print(m1$cptable)

m2 <- prune(m1,cp=m1$cptable[which.min(m1$cptable[,"xerror"])])
plot(as.party(m2), tp_args=list(id=FALSE),main=mainTitle,font=2)


summary(m2)


m <- lm(nChem ~ Urban.pct + Ag.pct,data=dat)
summary(m)

plot(dat$nChem,predict(m,newdata=dat))
abline(0,1,col="blue",lty=2)
##################
```
                       
