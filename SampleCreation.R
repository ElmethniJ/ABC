####################################################################################
# Definition and initialization of variables and parameters
####################################################################################

FonctCreerEch<-function(ChoixLoi)
{
  Sample=matrix(0,nr=N,nc=n)	# empty matrix of N samples of size n
  EchU=matrix(runif(n*N,0,1),nr=N,nc=n)   # matrix of N samples of size n drawn according to a uniform distribution (0,1)
  dfech=1/EchU  # Matrix of the inverses of the uniforms to be given to FonctQ

  if(ChoixLoi == 1) {
    
  ######## RPD ######## 
  # Quantile function
  FonctQ<-function(x,gammaSecond,RhoSecond,betaSecond)
  {
    (x^(gammaSecond))*(exp(betaSecond*gammaSecond*((x^(RhoSecond)-1)/RhoSecond)))
  }
  
  # M samples of size n generate with FonctQ
  Sample=t(sapply(seq_len(nrow(dfech)), function(i) FonctQ(dfech[i,], gammaSecond, RhoSecond, betaSecond)))
  
  # Theoretical quantile to compare with
  QuantSecond=FonctQ(1/pn,gammaSecond,RhoSecond,betaSecond)	
  Quant = QuantSecond

  } else if(ChoixLoi == 2){
    
    ######## Burr ########
    # Quantile function
    FonctQB<-function(y,gammaB,RhoB)
    {
      ((y^(-RhoB))-1)^(-gammaB/RhoB)
    }
    
    # M samples of size n generate with FonctQB
    Sample=FonctQB(dfech,gammaB,RhoB)
   
    # Theoretical quantile to compare with
    QuantBurr=FonctQB(1/pn,gammaB,RhoB)
    Quant = QuantBurr
    
  } else if(ChoixLoi == 3){
    
    ######## Frechet ######## 
    # Quantile function
    FonctQF<-function(y,gammaF)
    {
      (log(1/(1-y)))^(-gammaF)
    }
    
    # M samples of size n generate with FonctQF
    Sample=FonctQF(EchU,gammaF)

    # Theoretical quantile to compare with
    QuantFrechet=FonctQF(pn,gammaF)
    Quant = QuantFrechet
    
  } else if(ChoixLoi == 4){
    
    ######## Student ########
    # M samples of size n generate with Student quantile function
    Sample=matrix(abs(rt(n*N,NuS)),nr=N,nc=n)
    
    # Theoretical quantile to compare with
    QuantS=qt((1-(pn/2)),NuS) 
    Quant = QuantS
    
  } else if(ChoixLoi == 5){
    
    ######## Fisher ########
    # M samples of size n generate with Fisher quantile function
    Sample=matrix(rf(n*N,NuF1,NuF2),nr=N,nc=n) 
    
    # Theoretical quantile to compare with
    QuantFi=qf(1-pn,NuF1,NuF2)
    Quant = QuantFi
    
  } else if(ChoixLoi == 6){
    
    ######## Inverse Gamma ########
    # M samples of size n generate with Inverse Gamma quantile function
    Sample=matrix(rinvgamma(n*N,AlphaInvg,rate = BetaInvg, scale = 1/BetaInvg),nr=N,nc=n)
    
    # Theoretical quantile to compare with
    QuantInvg=qinvgamma(1-pn,AlphaInvg, rate=BetaInvg,scale=1/rate)
    Quant = QuantInvg
    
  } else {
    
    ######## GPD ########
    # M samples of size n generate with GPD quantile function
    Sample=matrix(revd(n*N,loc = 0, scale = 1, shape = gammaG, threshold = 0,type =  "GP"),nr=N,nc=n)
    
    # Theoretical quantile to compare with
    QuantGPD=qevd(1-pn,loc=0,scale=1,shape=gammaG,type="GP")
    Quant = QuantGPD
    
}
  ListeRetour <- list(Sample, Quant)  
  # Returns the matrix of samples created and the theoretical quantile.
  return(ListeRetour) 
}

