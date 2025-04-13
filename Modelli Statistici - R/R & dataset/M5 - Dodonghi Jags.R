#### #### #### #### #### #### ####
#### Dugongs
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
#### The data are length (Y) and age (x) measurements for 27 captured dugongs (seacows- mucche di mare).
#### Carlin and Gelfand (1991) model this data using a nonlinear growth curve with no inflection
#### point and an asymptote as X i tends to infinity:
#### #### #### #### #### ####

#setwd("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/didattica/Modelli Statistici/Codici R/WinBugs/")

### ### ### ### ### ### ###
### I dati
### ### ### ### ### ### ###
# valori x
x = c( 1.0,  1.5,  1.5,  1.5, 2.5,   4.0,  5.0,  5.0,  7.0,
	            8.0,  8.5,  9.0,  9.5, 9.5,  10.0, 12.0, 12.0, 13.0,
	           13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5)
# valori y
Y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
	           2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
	           2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57)


plot(x,Y)
#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####

# lista con tutti gli elementi che servono al modello
dataList = list(
	Y   = Y,
	x   = x,
	N   = length(Y)
)

#### Modello Bayesiano
mod1_string <- 'model {
	for(i in 1:N)
	{
		Y[i] ~ dnorm(mu[i], tau)
		mu[i] <- alpha - beta*pow(gamma,x[i])
  }
	alpha ~ dunif(0,200)
	beta ~ dunif(0,200)
	gamma ~ dunif(0.5, 1.0)
	tau ~ dgamma(0.001, 0.001)
	sigma <- 1/sqrt(tau)
	U3 <- logit(gamma)

	for(i in 1:N)
	{
		res[i] = Y[i]-mu[i]
	}

}'


# parametri da salvare
SavePar = c("alpha", "beta", "gamma", "res")

#inits
inits =  list(
  list("alpha" = 1, "beta" = 1)
)
# fittiamo il modello
set.seed(1)
model = jags(
			data 	= dataList,
			parameters.to.save 	= SavePar,
			inits					= inits,
			model.file 		= textConnection(mod1_string),
			n.chains 						= 1,
			n.iter 						= 15000,
			n.burnin 						= 5000,
			n.thin 							= 5,
			DIC 								= T
)

#### #### #### #### #### #### ####
#### stime
#### #### #### #### #### #### ####

# alpha
Wpar = which(substr(names(model $BUGSoutput$sims.array[1,1,]),1,3)%in%c("alp","bet","gam"))
##
parPOST = as.mcmc(model$BUGSoutput$sims.array[,1,Wpar])
summary(parPOST)
plot(parPOST)
acf(parPOST)

### funzione media - alpha - beta*pow(gamma,x[i])

# valori di x da usare
xseq = seq(min(x), max(x), by=0.5)

PosterioFunction = matrix(NA, nrow=nrow(parPOST), ncol=length(xseq))

for(imcmc in 1:nrow(parPOST))
{
	PosterioFunction[imcmc,] = parPOST[imcmc,"alpha"]-parPOST[imcmc,"beta"]*parPOST[imcmc,"gamma"]^xseq
}

## esempi delle realizzazioni della funzione a posteriori

# tutte
par(mfrow=c(1,1))
plot(xseq,PosterioFunction[1,], type="l", col=1,ylim=c(1.6,2.8))
for(i in 1:nrow(PosterioFunction))
{
	lines(xseq,PosterioFunction[i,],  col=i)
}

# solo alcune
plot(xseq,PosterioFunction[1,], type="l", col=1,ylim=c(1.6,2.8))
lines(xseq,PosterioFunction[150,],  col=3)
lines(xseq,PosterioFunction[250,],  col=4)

### calcoliamo la media a posteriori e intervalli di credibilità
PostMeam = colMeans(PosterioFunction)
PostQ1   = apply(PosterioFunction,2,quantile, prob=0.025)
PostQ2   = apply(PosterioFunction,2,quantile, prob=1-0.025)


plot(xseq,PostMeam, type="l", col=1, ylim=c(1.6,2.8))
lines(xseq,PostQ1, col=2, lty=3)
lines(xseq,PostQ2, col=2, lty=3)
points(x,Y)
### Vediamo qualche residuo
Wres = which(substr(names(model $BUGSoutput$sims.array[1,1,]),1,3)=="res")
##
resPOST = model$BUGSoutput$sims.array[,1,Wres]
resMeans = colMeans(resPOST)

par(mfrow=c(1,2))
plot(resMeans)
plot(x,resMeans)
par(mfrow=c(1,1))
