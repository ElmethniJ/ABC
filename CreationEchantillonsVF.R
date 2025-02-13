# _________________________________________________________________________
# Definition et initialisation des variables et des parametres

##################### Creation #####################

FonctCreerEch<-function(ChoixLoi)
{
  Sample=matrix(0,nr=N,nc=n)	# matrice des N echantillons de taille n
  # Creation des N echantillons de taille n tires selon une loi uniforme 0, 1
  EchU=matrix(runif(n*N,0,1),nr=N,nc=n)   # matrice des N echantillons de taille n tires selon une loi uniforme 0, 1
  
  # Matrice des Uniformes a donner a FonctQ
  dfech=1/EchU       # matrice des inverses des elemnts de la matrice EchU !

  if(ChoixLoi == 1) {
  #### Loi du second ordre   
  # Parametres
  # gammaSecond=0.25  
  # RhoSecond=-0.5   
  # betaSecond=1
  # ________________________________________________
  # Fonction quantile du second ordre
  FonctQ<-function(x,gammaSecond,RhoSecond,betaSecond)
  {
    (x^(gammaSecond))*(exp(betaSecond*gammaSecond*((x^(RhoSecond)-1)/RhoSecond)))
  }
  
  #for(j in 1:N)
  #{
  # Tirage d'un echantillon de taille n de la loi choisie
  #  Sample[j,]=FonctQ(1/EchU[j,],gammaSecond,RhoSecond,betaSecond)
  #}
  
  # Les N echantillons (quantiles) de taille n selon FonctQ
  Sample=t(sapply(seq_len(nrow(dfech)), function(i) FonctQ(dfech[i,], gammaSecond, RhoSecond, betaSecond)))
  
  
  # ________________________________________________
  # Quantile theorique auquel se comparer
  QuantSecond=FonctQ(1/pn,gammaSecond,RhoSecond,betaSecond)	# Quantile theorique
  Quant = QuantSecond

  # apply(Sample,1,max)
  } else if(ChoixLoi == 2){
    ######## Loi de Burr ########
    
    # gammaB=0.5    # Loi de Burr
    # RhoB=-0.25    # rho, Loi de Burr
    # ________________________________________________
    # Fonction quantile de Burr : 
    FonctQB<-function(y,gammaB,RhoB)
    {
      ((y^(-RhoB))-1)^(-gammaB/RhoB)
    }
    
    # curve(FonctQB(x,gammaB,RhoB),from=1,to=100)
    
    # ________________________________________________
    # Creation des echantillons
    Sample=FonctQB(dfech,gammaB,RhoB)
    # ________________________________________________
    # Quantile theorique auquel se comparer
    QuantBurr=FonctQB(1/pn,gammaB,RhoB)     # Quantile theorique : Loi de Burr
    Quant = QuantBurr
    
    # apply(Sample,1,max)
  } else if(ChoixLoi == 3){
    ###### Loi de Frechet #######
    # gammaF=0.25       # Loi de Frechet implique Rho
    # RhoF=-1        # rho, Loi de Frechet
    # ________________________________________________
    # Fonction quantile de Frechet : Inverse de la fonction de survie
    FonctQF<-function(y,gammaF)
    {
      (log(1/(1-y)))^(-gammaF)
    }
    # ________________________________________________
    # Creation des echantillons
    Sample=FonctQF(EchU,gammaF)
    # ________________________________________________
    # Quantile theorique auquel se comparer
    QuantFrechet=FonctQF(pn,gammaF)          # Quantile theorique : Loi de Frechet
    Quant = QuantFrechet
    
    # apply(Sample,1,max)
    
    # apply(Sample,1,max)
  } else if(ChoixLoi == 4){
    ###### Loi de Student #######
    # NuS=2            # Ddl de la Loi de Student
    # gammaS=1/NuS     # Loi de Student
    # RhoS=-2/NuS      # rho, Loi de Student
    Sample=matrix(abs(rt(n*N,NuS)),nr=N,nc=n)  # Tirage d'un echantillon de taille n de la loi de Student
    QuantS=qt((1-(pn/2)),NuS)           # Quantile theorique : Loi de Student
    Quant = QuantS
    
    # apply(Sample,1,max)
  } else if(ChoixLoi == 5){
    ####### Loi de Fisher #######
    # NuF1=2            # Ddl1 de la Loi de Fisher
    # NuF2=3            # Ddl2 de la Loi de Fisher
    # gammaFi=2/NuF2        # Loi de Fisher
    # RhoFi=-2/NuF2       # rho, Loi de Fisher
    Sample=matrix(rf(n*N,NuF1,NuF2),nr=N,nc=n) # Tirage d'un echantillon de taille n de la loi de Fisher
    QuantFi=qf(1-pn,NuF1,NuF2)       # Quantile theorique : Loi de Fisher
    Quant = QuantFi
    
    # apply(Sample,1,max)
  } else if(ChoixLoi == 6){
    ###### Loi de InvGamma ######
    #############################
    # BetaInvg=1      # ou 2, arbitraire ! determine le rate
    # AlphaInvg=8     # AlphaInvg=-1/RhoInvg
    # RhoInvg=-1/AlphaInvg
    # gammaInvg=1/AlphaInvg
    # Tirage d'un echantillon de taille n de la loi InvGamma
    Sample=matrix(rinvgamma(n*N,AlphaInvg,rate = BetaInvg, scale = 1/BetaInvg),nr=N,nc=n)
    QuantInvg=qinvgamma(1-pn,AlphaInvg, rate=BetaInvg,scale=1/rate)  # Quantile theorique : Loi InvGamma
    Quant = QuantInvg
    
    # apply(Sample,1,max)
  } else {
    ########## Loi GPD ##########
    # RhoG=-0.5        # rho GPD
    # gammaG=-RhoG       # Loi GPD
    Sample=matrix(revd(n*N,loc = 0, scale = 1, shape = gammaG, threshold = 0,type =  "GP"),nr=N,nc=n)
    QuantGPD=qevd(1-pn,loc=0,scale=1,shape=gammaG,type="GP")         # Quantile theorique : GPD
    Quant = QuantGPD
    
    # apply(Sample,1,max)
}
  ListeRetour <- list(Sample, Quant)  
  return(ListeRetour)        # Retourne la matrice des echantillons crees ainsi que le quantile theorique.
}
# ________________________________________________
