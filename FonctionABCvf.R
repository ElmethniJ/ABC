
############################# Fonction ABC #################################

# Reference : echantillon de reference auquel se comparer
# alphan : ordre du quantile 
# n : Taille fixee des echantillons a tirer
# M:  Nombre d'echantillons a tirer pour une taille fixee
# tau : Pourcentage d'echantillons a garder
# sigmaR : ecart-type de la loi gamma a priori sur rho
# sigmaB : ecart-type de la loi gamma a priori sur beta

FonctABC<-function(Reference,taille,M,tau,sigmaR,sigmaB,alphan,gammagomes,rhogomes,betagomes)
{
  # Matrice vide des echantillons et des delta
  nk=length(Reference)
  SampleT=matrix(0,nr=M,nc=nk)	# matrice des NbEch echantillons tries
  Delta=matrix(0,nr=M,nc=1)	# matrice des erreurs au sens du transport optimal
  QuantileABC=matrix(0,nr=(tau*M),nc=1)	# matrice des quantiles simules
  CTEABC=matrix(0,nr=(tau*M),nc=1)	# matrice des CTE simules
  
  # Creation des echantillons
  EchU=matrix(runif(nk*M,0,1),nr=M,nc=nk)
  
  # Lois a priori
  
  # Gamma tire selon une loi uniforme
  gammaC=matrix(runif(M,0,max(gammagomes[1:(taille/2)])),nr=M,nc=1)

  # Rho tire selon une loi normale tronquee
  rhoC=matrix(rtruncnorm(M, a=-Inf, b=-0.001, mean = rhogomes, sd = sigmaR),nr=M,nc=1)
  
  # Beta tire selon une loi normale tronquee
  betaC=matrix(rtruncnorm(M, a=-1, b=Inf, mean = betagomes, sd = sigmaB),nr=M,nc=1)
  
  # Matrice des Uniforme a donner a FonctQ
  dfech=1/EchU 
  # Les M echantillons de taille n selon FonctQ
  SampleT=rowSort(t(sapply(seq_len(nrow(dfech)), function(i) FonctQ(dfech[i,], gamma = gammaC[i], rho = rhoC[i], beta=betaC[i]))))
  
  # Calcul de lecart au sens du transport optimal
  ReferenceT=rowSort(matrix(Reference, nrow=M, ncol=length(Reference), byrow=TRUE))
  Delta=rowSums((abs( log(abs(SampleT)) - log(abs(ReferenceT)) )))
  
  # Stockage dans un dataframe des echantillons et parametres correspondants
  df=data.frame(Gamma=gammaC,Rho=rhoC,Beta=betaC,Delta)
  
  # Trie de toutes les colonnes du dataframe selon Delta
  df_T=df[order(df$Delta),] 
  
  # On ne retient que les meilleurs 
  df_final=df_T[1:(tau*M),]

  # Generation des quantiles estimes pour les tau*M meilleurs echantillons
  QuantileABC=sapply(seq_len(nrow(df_final)), function(i) log(FonctQ(1/alphan, gamma = df_final$Gamma[i], rho = df_final$Rho[i], beta=df_final$Beta[i])))
  
  # Generation des quantiles estimes pour les tau*M meilleurs echantillons
  CTEABC=sapply(seq_len(nrow(df_final)), function(i) FonctCTE(1/alphan, gamma = df_final$Gamma[i], rho = df_final$Rho[i], beta=df_final$Beta[i]))
  
  # Data frame des resultats
  results=data.frame(df_final,QuantileABC,CTEABC)
}

# ________________________________________________
# Fonction donnant les exces Y
Exces<-function(Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(RefT[1:k]/RefT[k+1])
}

# ________________________________________________
# Fonction donnant le seuil X n-k,n
Seuil<-function(Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(RefT[k+1])
}

# ________________________________________________
# Fonction donnant Hill
GammaHill<-function(Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(mean(log(RefT[1:k]/RefT[k+1])))
}

# ________________________________________________
# Fonction Moments notee PWM
FonctPWM<-function(a,Ref,k)
{
  RefT=sort(Ref,decreasing = T)
  return(  ((a+1)^2)*mean( ( ( (1:k)/k )^a )*log(RefT[1:k]/RefT[k+1]) )  )
}

# ________________________________________________
# Fonction Quantile notee U
FonctQ<-function(x,gamma,rho,beta)
{
  (x^(gamma))*(exp(beta*gamma*((x^(rho)-1)/rho)))
}

#________________________________________________
#Fonction CTE notee FonctCTE
FonctCTE<-function(x,gamma,rho,beta)
{
  integrand2<-function(y) {((y)^(-gamma))*exp((gamma*beta/rho)*((y)^(-rho)))}
  (-gamma*beta/rho)+log(x)+log(integrate(integrand2,lower=0,upper=1/x)$value)
}

