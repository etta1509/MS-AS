### ### ### ### ### ### ### ### ### 
### Si vuole determinare le variabili che influenzano 
### il salario
### ### ### ### ### ### ### ### ### 

## 
rm(list =ls())
## 


## librerie
library(contrast)
library(deldir)
library(tidyverse)
library(aod)
##


# carichiamo i dati
Data = read.csv("guadagno.csv",sep=",")
summary(Data)
head(Data)


## puliamo e sistemiamo il dataset 
#Data = Data[,-1]
recast_Data <- Data %>% select(-x)%>%select(-workclass) %>%
  mutate(education = factor(ifelse(education == "Preschool" | education == "10th" | education == "11th" | education == "12th" | education == "1st-4th" | education == "5th-6th" | education == "7th-8th" | education == "9th", "dropout", ifelse(education == "HS-grad", "HighGrad", ifelse(education == "Some-college" | education == "Assoc-acdm" | education == "Assoc-voc", "Community", ifelse(education == "Bachelors", "Bachelors", ifelse(education == "Masters" | education == "Prof-school", "Master", "PhD")))))))

recast_data_2 <-  recast_Data %>% mutate(marital.status = factor(ifelse(marital.status == "Never-married" | marital.status == "Married-spouse-absent", "Not_married", ifelse(marital.status == "Married-AF-spouse" | marital.status == "Married-civ-spouse", "Married", ifelse(marital.status == "Separated" | marital.status == "Divorced", "Separated", "Widow")))))

Data = recast_data_2


## ## ## ## ## ## ## ##  
## Descrittive
## ## ## ## ## ## ## ## 


###  Numeriche
Data$income = as.factor(Data$income)
plot(Data$age ~ Data$income, col=2:3)
plot(Data$hours.per.week ~ Data$income, col=2:3)
plot(Data$educational.num ~ Data$income, col=2:3)
plot(Data$educational.num~Data$education)
# riscaliamo le variabili
Data_rescale <- Data 
for(i in 1:ncol(Data_rescale))
{
  if(is.numeric(Data_rescale[,i]))
  {
    Data_rescale[,i] = scale(Data_rescale[,i])
  }
}

###  Fattori

attributes(as.factor(Data$gender))
table(Data$gender)
table(Data$gender, Data$income)/nrow(Data)
table(Data$gender, Data$income)/matrix(table(Data$income), ncol=2, nrow=2, byrow=T)
table(Data$gender, Data$income)/matrix(table(Data$gender), ncol=2, nrow=2, byrow=F)

attributes(Data$ marital.status)
table(Data$marital.status, Data$income)/matrix(table(Data$income), ncol=2, nrow=4, byrow=T)
table(Data$marital.status, Data$income)/matrix(table(Data$marital.status), ncol=2, nrow=4, byrow=F)


attributes(as.factor(Data$race))
table(Data$race, Data$income)/matrix(table(Data$income), ncol=2, nrow=5, byrow=T)
table(Data$race, Data$income)/matrix(table(Data$race), ncol=2, nrow=5, byrow=F)

attributes(Data$education )
table(Data$education, Data$income)/matrix(table(Data$income), ncol=2, nrow=6, byrow=T)
table(Data$education, Data$income)/matrix(table(Data$education), ncol=2, nrow=6, byrow=F)



## ## ## ## ## ## ## ##  
## Modello logistico
## ## ## ## ## ## ## ## 

?glm
?family
?make.link
utils::str(make.link("logit"))

# se si voglino cambiare i reference dei fattori
# usare relevel
Data$education = relevel(Data$education,ref="PhD")
Mod1 = glm(income ~ age+education+educational.num+marital.status+race+gender+hours.per.week, data=Data, family=binomial(link = "logit"))

summary(Mod1)


XX = model.matrix(Mod1)
colnames(XX)
W = which(rowSums(XX[,-c(1,2,8,17)])==0 )
# Possiamo fare dei contrasti
Cont1 = contrast(Mod1 , 
              list(education="PhD" ,marital.status = "Not_married", race="White", gender="Female",age=26, educational.num= 7, hours.per.week=40),
              list(education="Bachelors" ,marital.status = "Not_married", race="White", gender="Female",age=26, educational.num= 7, hours.per.week=40)) # type="individual"
print(Cont1, X=T)
print(Cont1)



# e possiamo calcolarli manualmente 
BetaPar = matrix(coef(Mod1),ncol=1)
JpowerLess1= summary(Mod1)$cov.unscaled #

VecConf = Cont1$X

Diff = VecConf%*%BetaPar # contrasto (differenza*varianza)
Var  = VecConf%*%JpowerLess1%*%t(VecConf)

zvalue = Diff/(Var)^0.5
zvalue
pnorm(zvalue,0,1,lower.tail=F)*2





## null deviance

ydic = as.numeric(Data$income)-1
pred = rep(mean(ydic), length(ydic))

W1 = which(ydic==1)
W0 = which(ydic==0)

### il modello è binomiale per la y (mi chiedo se è > o < di 50)
# devianza del modello nullo (dal modello saturo)
D0 =  2*(
  sum(log(ydic/pred)[W1])+sum(log((1-ydic  )/(1-pred))[W0])
)
D0

gdl0 = length(ydic)-1 ### togliamo l'intercetta
gdl0

## Model DEviance (residuals)
?predict
pred = predict(Mod1,type="response")

W1 = which(ydic==1)
W0 = which(ydic==0)

### devianza del mio modello (dal modello saturo)
D1 = 2*(
  sum(log(ydic/pred)[W1])+sum(log((1-ydic  )/(1-pred))[W0])
)
gdl1 = length(ydic)-length(BetaPar)

D1
gdl1

## testiamo se il modello è meglio del nullo
deltaD = D0-D1
# probabilità che la chiquadro con gdl0-gdl1 gradi di libertà
# sia maggiore del valore deltaD
pchisq(deltaD, gdl0-gdl1, lower.tail=F)

## vediamo se il modello è buono come il saturo
pchisq(D1, gdl1, lower.tail=F)
# rappresentazione grafica
xseq = seq(35000,50000, length.out=100)
plot(xseq,dchisq(xseq,gdl1), type="l")
abline(v=D1,col=2)

xseq = seq(0,20000, length.out=100)
plot(xseq,dchisq(xseq,gdl0-gdl1), type="l")
abline(v=D0-D1,col=2)

## analisi dei residui
?residuals.glm

p.res   = residuals.glm(Mod1, type="pearson")
dev.res =  residuals.glm(Mod1, type="deviance")
p.res.stand = rstandard(Mod1, type="pearson")
# Possiamo vedere come si relazionano alle covariate
plot(Data$age, p.res)
# o ai valori predetti -  attenzione che le indipendenze non sono valide qui
plot(pred, p.res)

# standard plots
plot(Mod1)

# oppur possiamo usare il predittore lineare
plot(predict(Mod1,type="link"), p.res)

## Test quadrato residui pearson
Z2 = sum(p.res^2)
pchisq(Z2, gdl1, lower.tail=T)


## possiamo vedere anche le misure di influenza
# attenzione che sono time-consuming
#InfMeas1 = influence.measures(Mod1)




## ## ## ## ## ## ## ## ## ## ## ## 
## Scelta del modello
## ## ## ## ## ## ## ## ## ## ## ## 

## con il test ANOVA possiamo verificare
## se le componenti del modello
# sono significative - 
# ATTENZIONE che le variabili sono
# testate nell'ordine dato dal modello



anova(Mod1,test="Chisq")

# se vogliamo testare modelli diversi facciamo l'update 
Mod2 = update(Mod1,update(Mod1$formula,    ~ . - age))
anova(Mod2,Mod1 ,test="Chisq")
Dconf = deviance(Mod2)-deviance(Mod1)

pchisq(Dconf, 1, lower.tail=F)
xseq = seq(0,600, length.out=100)
plot(xseq,dchisq(xseq,1), type="l")
abline(v=Dconf,col=2)

