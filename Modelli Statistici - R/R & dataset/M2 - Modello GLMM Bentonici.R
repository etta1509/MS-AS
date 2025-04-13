##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####   ##### #####
#In Questo dataset ci sono dati bentonici marini (organismi  che vivono in stretto contatto
#con il fondo o fissati ad un substrato solido), di 9 spiagge tedesche, osservati vicino alla
#costa.
#Per ogni spiaggia sono stati raccolti 5 campioni, e sono stati misurati Richness
#(Numero di specie diverse) il NAP (Differenza tra punto di campionamento e livello del mare) e
#l'Esposure (un indice composto dall'azione delle onde, inclinazione, tipo di sabbia etc).
#L'interesse � capire come la richness dipende dalle altre variabili
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####   ##### #####

rm(list =ls())


# LIBRARIES  -------------------------------------------------------------

library(lme4) 
library(nlme)
library(lattice)
# Functions
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

# Settiamo la directory
# DIR = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/2024 - 01VJWNG - Modelli statistici_Apprendimento statistico (Model - 273975/script per esame/Mastrantonio/datasets/"
# setwd(DIR)

# Carichiamo i dati -  Ricordatevi di controllare separatori e decimali
Data = read.csv("Dati Bentonici.txt", sep="", dec=".", header=T)
summary(Data)

# qualche plot
pairs(Data[,2:4], upper.panel=panel.cor,diag.panel=panel.hist)

par(mfrow=c(1,2))
scatter.smooth(x=Data[,3], y=Data[,2], pch=16, col="red", xlab=colnames(Data)[3], ylab=colnames(Data)[2])
scatter.smooth(x=Data[,4], y=Data[,2], pch=16, col="red", xlab=colnames(Data)[4], ylab=colnames(Data)[2])

##### ##### ##### ##### ##### #####
# Modello GLM semplice
##### ##### ##### ##### ##### #####

GLM1 = glm(Richness  ~ Exposure+NAP, data= Data, family=poisson)
summary(GLM1)

# Modello GLM semplice con interazione: l'effetto di exposure varia con nap
GLM2 = glm(Richness  ~ Exposure+NAP+Exposure:NAP, data= Data, family=poisson)
summary(GLM2)

# Possiamo ipotizzare che gli effetti siano diversi per spiaggia
# l'effetto di exposure varia con nap e quello di nap con beach
GLM3wrong = glm(Richness  ~ Exposure:Beach+NAP:Beach, data= Data, family=poisson)
summary(GLM3wrong)
Data$Beach = as.factor(Data$Beach)

# Modello Diviso per spiaggie
GLM3wrong2 = glm(Richness  ~ Exposure+NAP+Beach+Exposure:Beach+NAP:Beach, data= Data, family=poisson)
summary(GLM3wrong2)
table(Data$Exposure, Data$Beach)

GLM3 = glm(Richness  ~ Exposure+NAP+Beach+NAP:Beach, data= Data, family=poisson)
summary(GLM3)


# Per verificarlo calcoliamo la matrice del disegno
Modmat = model.matrix( ~ Exposure+NAP+Beach+NAP:Beach, data= Data)
summary(Modmat)
# e modellizziamo, tramite lm, la settima spiaggia in funzione delle precedenti
TestDep = lm(Modmat[,"Beach7"] ~ -1+Modmat[,1:8])
summary(TestDep)



# Un modello alternativo: togliamo l'intercetta e mettiamo un'interazione l'effetto di nap varia con beach
GLM3_alt1 = glm(Richness  ~-1+ Exposure+NAP+Beach+NAP:Beach, data= Data, family=poisson)
summary(GLM3_alt1)

# scelta modello: AIC minore
AIC(GLM1,GLM2,GLM3)

# calcoliamo residui e fitted values del modello scelto
# scegliere il tipo di residui
# ?predict.glm - type
resGLM3     = residuals(GLM3)
fitGLM3     = fitted(GLM3)

# plot residui colorati per spiaggia
par(mfrow=c(2,2))
plot(resGLM3, col=c(1:9)[Data$Beach], pch=20)
plot(resGLM3, fitGLM3, col=c(1:9)[Data$Beach], pch=20)
plot(Data$NAP,resGLM3, col=c(1:9)[Data$Beach], pch=20)

# diagnostiche generali
par(mfrow=c(2,2))
plot(GLM3)





##### ##### ##### ##### ##### #####
# Modello a effetti misti
##### ##### ##### ##### ##### #####
# introdurre spiaggia "costa" 8 gradi di libert�
# se non siamo interessati a sapere il valore di preciso per ogni spiaggia
# possiamo trattarla come un effetto random

GLM1_app    = glm(Richness ~ NAP+Exposure+Beach, data=Data, family=poisson)
summary(GLM1_app )


# modello a effetti misti con effetto random su spiaggia: la media di y varia al variare di beach
GLMM1    = glmer(Richness ~ NAP+Exposure+(1|Beach), data=Data, family=poisson)
summary(GLMM1)

GLMM1_2    = glmer(Richness ~ 1- NAP+Exposure+(1|Beach), data=Data, family=poisson)
summary(GLMM1_2)

AIC(GLMM1,GLMM1_2)

# estraiamo i coefficienti random dell'intercetta
coef(GLMM1)
# estraiamo i valori specifici della spiaggia
ranef(GLMM1)
# differenza tra i due effetti
coef(GLMM1)$Beach[,1]-ranef(GLMM1)$Beach
# la differenza da l'intercetta
summary(GLMM1)


# Plottiamo la matrice Z di y= X*beta+Z*b+epsilon, b = effetto random
Zmat1=getME(GLMM1,"Z")
image(Zmat1)

# calcoliamo la matrice di covarianza degli effetti random
# la varianza
Cov = as.numeric((summary(GLMM1))$varcor$Beach)
# MatCov1 = Cov(b)
MatCov1 = diag(Cov,9)
#MatCov = Cov(Z*b) = ZCov(b)Z'
MatCov = Zmat1%*%MatCov1%*%t(Zmat1)



# testiamo una random slope su NAP: la variabilità tra i gruppi si vede anche al di fuori della media
GLMM2    = glmer(Richness ~ NAP+(-1+NAP|Beach), data=Data, family=poisson)
summary(GLMM2)

# estraiamo i coefficienti random dell'intercetta
coef(GLMM2)
# estraiamo i valori specifici della spiaggia
ranef(GLMM2)
# differenza tra i due effetti
coef(GLMM2)$Beach[,2]-ranef(GLMM2)$Beach
# la differenza da l'intercetta
summary(GLMM2)


# Plottiamo la matrice Z di y= X*beta+Z*b+epsilon, b = effetto random
Zmat2=getME(GLMM2,"Z")
image(Zmat2)

# calcoliamo la matrice di covarianza degli effetti random
# la varianza
Cov = as.numeric((summary(GLMM2))$varcor$Beach)
# MatCov1 = Cov(b)
MatCov1 = diag(Cov,9)
#MatCov = Cov(Z*b) = ZCov(b)Z'
MatCov = Zmat2%*%MatCov1%*%t(Zmat2)



# modello con random slope e intercetta
GLMM3    = glmer(Richness ~ NAP+(-1+NAP|Beach) + (1|Beach), data=Data, family=poisson)
summary(GLMM3)
coef(GLMM3)

Zmat3=getME(GLMM3,"Z")
image(Zmat3)


# calcoliamo la matrice di covarianza degli effetti random
# la varianza
CovNAP      = as.numeric((summary(GLMM3))$varcor[1])
CovBeach    = as.numeric((summary(GLMM3))$varcor[2])
# MatCov1 = Cov(b)
MatCovNAP1       = diag(CovNAP,9)
MatCovBeach1     = diag(CovBeach,9)
#MatCov = Cov(Z*b) = ZCov(b)Z'
MatCovNAP       = Zmat3[,1:9]%*%MatCovNAP1%*%t(Zmat3[,1:9])
MatCovBeach     = Zmat3[,9+1:9]%*%MatCovBeach1%*%t(Zmat3[,9+1:9])

D = (MatCovNAP +MatCovBeach)




# modello con random slope e intercetta
GLMM4    = glmer(Richness ~ NAP+(NAP|Beach), data=Data, family=poisson)
summary(GLMM4)
coef(GLMM4)

Zmat4=getME(GLMM4,"Z")
image(Zmat4)


# calcoliamo la matrice di covarianza degli effetti random
# la varianza
Cov = as.numeric((summary(GLMM4))$varcor$Beach)
Cov = matrix(Cov, ncol=2, nrow=2)
# MatCov1 = Cov(b)
# La funzione kronecker crea la matrice a blocchi diagonali
MatCov1 = kronecker(diag(1,9),Cov)
#MatCov = Cov(Z*b) = ZCov(b)Z'
MatCov = Zmat4%*%MatCov1%*%t(Zmat4)




# Scelta del modello: modello a effetti misti con effetto random su spiaggia
AIC(GLMM1,GLMM2,GLMM3,GLMM4)