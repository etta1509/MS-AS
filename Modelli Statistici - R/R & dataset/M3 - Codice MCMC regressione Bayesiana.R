#### #### #### #### #### #### #### ####
#### Stima modelli Lineare - MCMC
#### #### #### #### #### #### #### ####

rm(list =ls())
   
# definiamo la funzione che ci permette di simulare da una
# normale multivariata e calcolarne la densità
rmnorm=function(n = 1, mean = rep(0, d), varcov)
{
  d <- if (is.matrix(varcov))
    ncol(varcov)
  else 1
  z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
  y <- t(mean + t(z))
  return(y)
}
dmnorm=function (x, mean = rep(0, d), varcov, log = FALSE)
{
  d <- if (is.matrix(varcov))
    ncol(varcov)
  else 1
  if (d > 1 & is.vector(x))
    x <- matrix(x, 1, d)
  n <- if (d == 1)
    length(x)
  else nrow(x)
  X <- t(matrix(x, nrow = n, ncol = d)) - mean
  Q <- apply((solve(varcov) %*% X) * X, 2, sum)
  logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
  logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
  if (log)
    logPDF
  else exp(logPDF)
}

library(coda)

## ## ## ## ## ##
## Modello Lineare - MCMC
## ## ## ## ## ##

# beta ~ N()
#sigma2 ~ IG()

set.seed(10)
n       = 100
beta    = c(1,2,1,-1)
x1      = runif(n, 0,1)
x2      = runif(n, -1,1)
x3      = rnorm(n, 0,1)
sigma2  = 1
y       = beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x3+rnorm(n, 0,sigma2^0.5)

Data = data.frame(y =y, x1=x1,x2=x2,x3=x3)


# si assume che beta ~ N(beta.mean, beta.variance) e sigma^2 ~ IG(sigma2.a, sigma2.b)
ModLin_Conj = function(formula,  beta.mean, beta.variance,sigma2.a, sigma2.b, start.beta, start.sigma2, iter, burnin, thin, Data)
{

 


  ### Definiamo l'osservazioni e la matrice X
  Cov     = model.matrix(formula,Data)
  NameObs = as.character(terms(formula)[[2]])
  Obs     = Data[,NameObs,drop=F] # drop serve a fare in modo che il risultato sia ancora una matrice

  nbeta = ncol(Cov)
  nobs  = nrow(Cov)

  ### definisco quanti campioni a posteriori devo salvare
  ### e gi oggetti che ritorna la funzione
  NsampleSave     = floor((iter-burnin)/thin)
  BetaSave        = matrix(NA,ncol=nbeta, nrow=NsampleSave )
  sigma2Save      = matrix(NA,ncol=1, nrow=NsampleSave )

  ### definisco i valori correnti di beta e sigma2
  betaMCMC   = start.beta
  sigma2MCMC = start.sigma2



  ## oggetti utili per l'MCMC
  appSamp = burnin
  i = 1
  XtX = t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
  Xty = t(Cov[i,,drop=F])*Obs[i,1]
  for(i in 2:nobs)
  {
    XtX = XtX+t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
    Xty = Xty+t(Cov[i,,drop=F])*Obs[i,1]
  }
  for(iMCMC in 1:NsampleSave)
  {
    for(jMCMC in 1:appSamp)
    {
      ## ## ## ## ## ## ## ##
      ## campioniamo Beta
      ## ## ## ## ## ## ## ##

      # prior parameters
      Vpost = solve(beta.variance)
      Mpost = solve(beta.variance)%*%matrix(beta.mean,ncol=1)

      # contributi delle osservazioni
      Vpost = Vpost+XtX/sigma2MCMC
      Mpost = Mpost+Xty/sigma2MCMC

      Vpost = solve(Vpost)
      Mpost = Vpost%*%Mpost

      #simlazione Gibbs
      betaMCMC = rmnorm(1,Mpost,Vpost)

      ## ## ## ## ## ## ## ##
      ## campioniamo sigma2
      ## ## ## ## ## ## ## ##

      a_post = sigma2.a+nobs/2
      
      ObsMinisXbeta =  Obs -Cov%*%matrix(betaMCMC,ncol=1)

      b_post = sigma2.b+sum(ObsMinisXbeta^2)/2

      sigma2MCMC = 1/rgamma(1,shape =  a_post, rate=b_post)


    }
    appSamp = thin

    BetaSave[iMCMC,]    = betaMCMC
    sigma2Save[iMCMC,]  = sigma2MCMC
  }
  return(list(Beta =BetaSave, sigma2 = sigma2Save ))
}




## ## ## ## ## ##
## Modello Lineare - MCMC
## ## ## ## ## ##

# si assume che beta ~ N(beta.mean, beta.variance) e sigma^2 ~ G(sigma2.a, sigma2.b)
ModLin = function(formula,  beta.mean, beta.variance,sigma2.a, sigma2.b, start.beta, start.sigma2, iter, burnin, thin, Data, sd.prop = 1)
{

  ### Definiamo l'osservazioni e la matrice X
  Cov     = model.matrix(formula,Data)
  NameObs = as.character(terms(formula)[[2]])
  Obs     = Data[,NameObs,drop=F] # drop serve a fare in modo che il risultato sia ancora una matrice

  nbeta = ncol(Cov)
  nobs  = nrow(Cov)

  ### definisco quanti campioni a posteriori devo salvare
  ### e gi oggetti che ritorna la funzione
  NsampleSave     = floor((iter-burnin)/thin)
  BetaSave        = matrix(NA,ncol=nbeta, nrow=NsampleSave )
  sigma2Save      = matrix(NA,ncol=1, nrow=NsampleSave )

  ### definisco i valori correnti di beta e sigma2
  betaMCMC   = start.beta
  sigma2MCMC = start.sigma2



  ## oggetti utili per l'MCMC
  appSamp = burnin
  i = 1
  XtX = t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
  Xty = t(Cov[i,,drop=F])*Obs[i,1]
  for(i in 2:nobs)
  {
    XtX = XtX+t(Cov[i,,drop=F])%*%Cov[i,,drop=F]
    Xty = Xty+t(Cov[i,,drop=F])*Obs[i,1]
  }
  for(iMCMC in 1:NsampleSave)
  {
    for(jMCMC in 1:appSamp)
    {
      ## ## ## ## ## ## ## ##
      ## campioniamo Beta
      ## ## ## ## ## ## ## ##

      # prior parameters
      Vpost = solve(beta.variance)
      Mpost = solve(beta.variance)%*%matrix(beta.mean,ncol=1)

      # contributi delle osservazioni
      Vpost = Vpost+XtX/sigma2MCMC
      Mpost = Mpost+Xty/sigma2MCMC

      Vpost = solve(Vpost)
      Mpost = Vpost%*%Mpost

      #simlazione Gibbs
      betaMCMC = rmnorm(1,Mpost,Vpost)

      ## ## ## ## ## ## ## ##
      ## campioniamo sigma2 _ Metropolis
      ## ## ## ## ## ## ## ##

      # propongo  un nuovo valore
      

      # calcolo l'acceptance rate
      # potete scrivelo in due modi
      # 1) la distribuzione proposal è un log normale (esponenziale di una normale) e la prior gamma
      # 2) invece di sigma2 lavoriamo con tau=log(sigma2), che ha proposal simmetrica, e la prior la dobbiamo trovare
      # con la regola di trasformazioni di variabili:
      # f(tau) = f(sigma2(tau)) |d sigma2(tau) / d tau|
      # il second metodo,  ha il vantaggio
      # che essendo la proposa simmetrica, si semplifica nel
      # calcolo di alpha
      # solo uno dei due deve essere non commentato


      
      # #### PROPOSTA SU tau #####
      # proposta normale e prior f(tau) = f(sigma2(tau)) |d sigma2(tau) / d tau|
      
      tauMCMC = log(sigma2MCMC)
      tau2prop = rnorm(1,tauMCMC,sd.prop)
      sigma2prop = exp(tau2prop)

      logalpha = 0
      logalpha = logalpha+(-(sigma2.a+1)*log(sigma2prop)-sigma2.b/sigma2prop+log(sigma2prop))
    
      logalpha = logalpha-(-(sigma2.a+1)*log(sigma2MCMC)-sigma2.b/sigma2MCMC+log(sigma2MCMC))

      ## DENSITà PROPOSTA - in realtà non serve che si semplifica
 
      

      # #### PROPOSTA SU sigma #####
      # # proposta log-normale e prior Gamma
      # tauMCMC = log(sigma2MCMC)
      # tau2prop = rnorm(1,tauMCMC,sd.prop)
      # sigma2prop = exp(tau2prop)

      # logalpha = 0
      # logalpha = logalpha+(-(sigma2.a+1)*log(sigma2prop)-sigma2.b/sigma2prop
    
      # logalpha = logalpha-(-(sigma2.a+1)*log(sigma2MCMC)-sigma2.b/sigma2MCMC

      # ## DENSITà PROPOSTA
      # logalpha = logalpha+dlnorm(tauMCMC,tau2prop, sd.prop)
      # logalpha = logalpha-dlnorm(tau2prop,tauMCMC, sd.prop) 


      # likelihood contribution
      for(i in 1:nobs)
      {
        logalpha = logalpha+dnorm(Obs[i,],Cov[i,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2prop^0.5, log=T)
        logalpha = logalpha-dnorm(Obs[i,],Cov[i,,drop=F]%*%matrix(betaMCMC,ncol=1), sigma2MCMC^0.5, log=T)
      }

      alpha = min(1,exp(logalpha))

      u = runif(1,0,1)
      
      if(u<alpha)
      {
        sigma2MCMC = sigma2prop
      }
    }
    appSamp = thin

    BetaSave[iMCMC,]    = betaMCMC
    sigma2Save[iMCMC,]  = sigma2MCMC
  }
  return(list(Beta =BetaSave, sigma2 = sigma2Save ))
}


#### #### #### #### #### #### ####
#### TESTIAMO IL CODICE CON DATI SIMULATI
#### #### #### #### #### #### ####
set.seed(10)
n       = 100
beta    = c(1,2,1,-1)
x1      = runif(n, 0,1)
x2      = runif(n, -1,1)
x3      = rnorm(n, 0,1)
sigma2  = 1
y       = beta[1]+beta[2]*x1+beta[3]*x2+beta[4]*x3+rnorm(n, 0,sigma2^0.5)

Data = data.frame(y =y, x1=x1,x2=x2,x3=x3)


Results = ModLin(
  y ~ x1+x2+x3,
  beta.mean = rep(0,4),
  beta.variance=diag(100,4),
  sigma2.a=1,
  sigma2.b=1,
  start.beta=rep(0,4),
  start.sigma2=1,
  iter=3000,
  burnin=1,
  thin=1,
  Data=Data,
  sd.prop = 1
)


ResultsGIBBS = ModLin_Conj(
  y ~ x1+x2+x3,
  beta.mean = rep(0,4),
  beta.variance=diag(100,4),
  sigma2.a=1,
  sigma2.b=1,
  start.beta=rep(0,4),
  start.sigma2=1,
  iter=3000,
  burnin=1,
  thin=1,
  Data=Data
)

par(mfrow=c(2,2))
plot(Results$Beta[,1], type="l")
lines(ResultsGIBBS$Beta[,1], col=2)
acf(Results$Beta[,1])

plot(Results$Beta[,2], type="l")
lines(ResultsGIBBS$Beta[,2], col=2)
acf(Results$Beta[,2])

plot(Results$Beta[,3], type="l")
lines(ResultsGIBBS$Beta[,3], col=2)
acf(Results$Beta[,3])

plot(Results$Beta[,4], type="l")
lines(ResultsGIBBS$Beta[,4], col=2)
acf(Results$Beta[,4])

plot(Results$sigma2[,1], type="l")
lines(ResultsGIBBS$sigma2[,1], type="l", col=2)
acf(Results$sigma2[,1])
par(mfrow=c(1,1))

ResultsMCMC = as.mcmc(cbind(Results$Beta,Results$sigma2))

colnames(ResultsMCMC) = c("beta0","beta1","beta2","beta3", "sigma2")
summary(ResultsMCMC)


acf(ResultsMCMC)
plot(Results$Beta[,1],Results$Beta[,2])
