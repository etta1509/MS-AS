#### #### #### #### #### #### ####
#### Eyes
#### #### #### #### #### #### ####

rm(list=ls())

#### librerie
library(R2jags)
library(sjmisc)
library(rjags)

#### #### #### #### #### ####
#### The data
#### #### #### #### #### ####
#### #### DESCRIZIONE
#### Bowmaker et al (1985) analyse data on the peak sensitivity wavelengths for individual
#### microspectrophotometric records on a small set of monkey's eyes.
#### Data for one monkey (S14 in the paper) are given below
#### (500 has been subtracted from each of the 48 measurements).
#### #### #### #### #### ####

#setwd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/didattica/Modelli Statistici/Codici R/WinBugs/")

### ### ### ### ### ### ###
### I dati
### ### ### ### ### ### ###

# valori y
Y = c(529.0, 530.0, 532.0, 533.1, 533.4, 533.6, 533.7, 534.1, 534.8, 535.3,
         535.4, 535.9, 536.1, 536.3, 536.4, 536.6, 537.0, 537.4, 537.5, 538.3,
         538.5, 538.6, 539.4, 539.6, 540.4, 540.8, 542.0, 542.8, 543.0, 543.5,
         543.8, 543.9, 545.3, 546.2, 548.8, 548.7, 548.9, 549.0, 549.4, 549.9,
         550.6, 551.2, 551.4, 551.5, 551.6, 552.8, 552.9,553.2)

# vediamo la densità
plot(density(Y))

#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####

# lista con tutti gli elementi che servono al modello
dataList = list(
	Y   = Y,
	N   = length(Y)
)

#### Modello Bayesiano
mod1_string <- '
model
{
	for(i in 1:N)
	{
		Y[i] ~ dnorm(mu[i], tau)
		mu[i] <- lambda[T[i]]
		T[i] ~ dcat(P[]) # equivalente a Bern(P[1])
	}
	alpha[1] = 1
	alpha[2] = 1
	P[1:2] ~ ddirich(alpha[]) # Equivalente a B(alpha1],alpha[2])
	theta ~ dunif(0.0, 1000)
	lambda[2] <- lambda[1] + theta
	lambda[1] ~ dnorm(0.0, 1.0E-6)
	tau ~ dgamma(0.001, 0.001)
	sigma <- 1 / sqrt(tau)
}
'


# parametri da salvare
SavePar = c("P", "lambda","sigma","T")

# inits
inits =  list(
  list(lambda = c(100, NA), theta = 50, tau = 1)
)
# fittiamo il modello
set.seed(1)
model = jags(
			data 	= dataList,
			parameters.to.save 	= SavePar,
			inits					= inits,
			model.file 		= textConnection(mod1_string),
			n.chains 						= 1,
			n.iter 						  = 25000,
			n.burnin 						= 15000,
			n.thin 							= 5,
			DIC 								= T
)

#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####

# alpha
Wpar = which(substr(names(model $BUGSoutput$sims.array[1,1,]),1,1)%in%c("P","l","s"))
##
parPOST = as.mcmc(model$BUGSoutput$sims.array[,1,Wpar])
summary(parPOST)
plot(parPOST)
acf(parPOST)

## classificazione

Wclass = which(substr(names(model $BUGSoutput$sims.array[1,1,]),1,1)%in%c("T"))
##
classPOST = as.mcmc(model$BUGSoutput$sims.array[,1,Wclass])
summary(classPOST)



#### #### #### #### #### #### ####
#### stime con missing
#### #### #### #### #### #### ####
indexNA = c(27,30,35)
Y[indexNA] = NA

# lista con tutti gli elementi che servono al modello
dataList = list(
	Y   = Y,
	N   = length(Y)
)

#### Modello Bayesiano
mod2_string <- '
model
{
	for(i in 1:N)
	{
		Y[i] ~ dnorm(mu[i], tau)
		mu[i] <- lambda[T[i]]
		T[i] ~ dcat(P[])
	}
	alpha[1] = 1
	alpha[2] = 1
	P[1:2] ~ ddirich(alpha[])
	theta ~ dunif(0.0, 1000)
	lambda[2] <- lambda[1] + theta
	lambda[1] ~ dnorm(0.0, 1.0E-6)
	tau ~ dgamma(0.001, 0.001)
	sigma <- 1 / sqrt(tau)

}
'

# parametri da salvare
SavePar = c("P", "lambda","sigma","T","Y")

# inits
inits =  list(
  list(lambda = c(100, NA), theta = 50, tau = 1)
)
# fittiamo il modello
set.seed(1)
model_NA = jags(
			data 	= dataList,
			parameters.to.save 	= SavePar,
			inits					= inits,
			model.file 		= textConnection(mod2_string),
			n.chains 						= 1,
			n.iter 						  = 25000,
			n.burnin 						= 15000,
			n.thin 							= 5,
			DIC 								= T
)

#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####

# alpha
Wpar = which(substr(names(model_NA$BUGSoutput$sims.array[1,1,]),1,1)%in%c("P","l","s"))
##
parPOST = as.mcmc(model_NA$BUGSoutput$sims.array[,1,Wpar])
summary(parPOST)
plot(parPOST)
acf(parPOST)

## classificazione

Wclass = which(substr(names(model_NA$BUGSoutput$sims.array[1,1,]),1,1)%in%c("T"))
##
classPOST = as.mcmc(model_NA$BUGSoutput$sims.array[,1,Wclass])
summary(classPOST)

## previsione
Wy = which(substr(names(model_NA$BUGSoutput$sims.array[1,1,]),1,1)%in%c("Y"))
##
yPOST = as.mcmc(model_NA$BUGSoutput$sims.array[,1,Wy])

## plot missing
plot(yPOST[,indexNA])

## esempio plot non missingj
plot(yPOST[,c(1,15,40)])

