
############################# ABC function #################################

# Reference : reference sample for comparison
# alphan : level of the quantile 
# n : fixed size of samples to be drawn
# M:  number of samples to be drawn for a fixed size
# tau : percentage of samples to be retained
# sigmaR : standard deviation of the a priori gamma distribution on rho
# sigmaB : standard deviation of gamma distribution a priori on beta

FonctABC<-function(Reference,taille,M,tau,sigmaR,sigmaB,alphan,gammagomes,rhogomes,betagomes)
{
  nk=length(Reference)                  # size of the reference sample
  SampleT=matrix(0,nr=M,nc=nk)	        # empty matrix of NbEch samples sorted
  Delta=matrix(0,nr=M,nc=1)	            # empty error matrix in the sense of optimal transport
  QuantileABC=matrix(0,nr=(tau*M),nc=1)	# empty matrix of simulated quantiles
  CTEABC=matrix(0,nr=(tau*M),nc=1)      # empty matrix of simulated CTE
  
  # Level of quantile creation according to a uniform distribution
  EchU=matrix(runif(nk*M,0,1),nr=M,nc=nk)
  
  # A priori distributions
  
  # Gamma sampled according to a uniform distribution
  gammaC=matrix(runif(M,0,max(gammagomes[1:(taille/2)])),nr=M,nc=1)

  # Rho sampled according to atruncated normal distribution
  rhoC=matrix(rtruncnorm(M, a=-Inf, b=-0.001, mean = rhogomes, sd = sigmaR),nr=M,nc=1)
  
  # Beta sampled according to atruncated normal distribution
  betaC=matrix(rtruncnorm(M, a=-1, b=Inf, mean = betagomes, sd = sigmaB),nr=M,nc=1)
  
  # Matrix of the inverses of the uniforms to be given to FonctQ
  dfech=1/EchU 
  # M samples of size n generate with FonctQ
  SampleT=rowSort(t(sapply(seq_len(nrow(dfech)), function(i) FonctQ(dfech[i,], gamma = gammaC[i], rho = rhoC[i], beta=betaC[i]))))
  
  # Computation of the errors between Reference and the samples
  ReferenceT=rowSort(matrix(Reference, nrow=M, ncol=length(Reference), byrow=TRUE))
  Delta=rowSums((abs( log(abs(SampleT)) - log(abs(ReferenceT)) )))
  
  # Dataframe storage of samples and corresponding parameters
  df=data.frame(Gamma=gammaC,Rho=rhoC,Beta=betaC,Delta)
  
  # Sorts all dataframe columns by errors
  df_T=df[order(df$Delta),] 
  
  # We only keep the tau*M best
  df_final=df_T[1:(tau*M),]

  # Generation of estimated quantiles for tau*M best samples
  QuantileABC=sapply(seq_len(nrow(df_final)), function(i) log(FonctQ(1/alphan, gamma = df_final$Gamma[i], rho = df_final$Rho[i], beta=df_final$Beta[i])))
  
  # Generation of estimated CTE for tau*M best samples
  CTEABC=sapply(seq_len(nrow(df_final)), function(i) FonctCTE(1/alphan, gamma = df_final$Gamma[i], rho = df_final$Rho[i], beta=df_final$Beta[i]))
  
  # Final results in data frame
  results=data.frame(df_final,QuantileABC,CTEABC)
}

# ________________________________________________
# Function giving exces Y
Exces<-function(Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(RefT[1:k]/RefT[k+1])
}

# ________________________________________________
# Function giving threshold X n-k,n
Seuil<-function(Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(RefT[k+1])
}

# ________________________________________________
# Hill function
GammaHill<-function(Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(mean(log(RefT[1:k]/RefT[k+1])))
}

# ________________________________________________
# PWM function
FonctPWM<-function(a,Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(  ((a+1)^2)*mean( ( ( (1:k)/k )^a )*log(RefT[1:k]/RefT[k+1]) )  )
}

# ________________________________________________
# Quantile function denoted U
FonctQ<-function(x,gamma,rho,beta)
{
  (x^(gamma))*(exp(beta*gamma*((x^(rho)-1)/rho)))
}

#________________________________________________
## CTE function denoted FonctCTE
FonctCTE<-function(x,gamma,rho,beta)
{
  integrand2<-function(y) {((y)^(-gamma))*exp((gamma*beta/rho)*((y)^(-rho)))}
  (-gamma*beta/rho)+log(x)+log(integrate(integrand2,lower=0,upper=1/x)$value)
}

