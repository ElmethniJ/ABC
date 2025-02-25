# Remove all
rm(list = ls())
gc()

# Useful packages
require(evir) 
require(evt0)    
require(invgamma) 
require(extRemes) 
require(Rfast)
require(ReIns)
require(truncnorm)
require(ggplot2)
require(pracma)
require(rmutil)

# Load reference sample creation program 
source("CreationEchantillonsVF.R")

# Load ABC function
source("FonctionABCvf.R")


# Step 1: choose sample size 
# n=500   # to reproduce simulations results
n=371     # for the real data set : besecura
# n= 555  # for the real data set : flood

# Step 2: choose number of replications
N=1    # for real data set
#N=500  # to reproduce simulations results

# Step 3: level of the quantile 
# pn=1/(2*n)  # to reproduce simulations results
pn=1/n        # for real data set
# pn=0.005    # for real data set

# Extrapolation factor
dn=(2:(n-1))/(n*pn)

# Set the seed
set.seed(12)

# Step 4: choose one of the 7 distribution by initializing the variable ChoixLoi
# 1=RPD, 2=Burr, 3=Frechet, 4=Student, 5=Fisher, 6 InvGamma, 7=GPD
ChoixLoi=4	

# Step 5: choose the parameters of the distribution

# RPD parameters   
gammaSecond=1/2 # gamma
RhoSecond=-1/2  # rho
betaSecond=1    # beta

# Burr parameters
gammaB=1/2  # gamma
RhoB=-1/2   # rho

# Frechet parameters
gammaF=1/2  # gamma
RhoF=-1     # rho fixed

# Student parameters
NuS=2         # Degrees of freedom
gammaS=1/NuS  # gamma
RhoS=-2/NuS   # rho

# Fisher parameters
NuF1=3           # First degrees of freedom
NuF2=4           # Second degrees of freedom
gammaFi=2/NuF2   # gamma
RhoFi=-2/NuF2    # rho

# Inverse Gamma parameters
BetaInvg=1      # Rate
AlphaInvg=2     # Shape
gammaInvg=1/AlphaInvg # gamma
RhoInvg=-1/AlphaInvg  # rho 


# GPD parameters
RhoG=-1/2       # rho 
gammaG=-RhoG    # gamma 


#### not to be modified but must be executed even for real data ####

LoiChoisie=function(choixLoi)
{
  switch(choixLoi,"Second Ordre ","Burr ", "Frechet ", "Student","Fisher", "InvGamma","GPD")
}

loi=LoiChoisie(ChoixLoi)	# choosen distribution

ChoixRho=function(choixLoi)
{
  switch(choixLoi,RhoSecond,RhoB, RhoF, RhoS,RhoFi, RhoInvg,RhoG)
}

rho=ChoixRho(ChoixLoi)	# choosen rho

ChoixGamma=function(choixLoi)
{
  switch(choixLoi,gammaSecond,gammaB, gammaF, gammaS,gammaFi, gammaInvg,gammaG)
}

gamma=ChoixGamma(ChoixLoi)	# choosen gamma

# Creation of the matrix of N samples of size n
ListeRetour=FonctCreerEch(ChoixLoi)
Reference=ListeRetour[[1]] # Matrix of N rows and n columns
QuantTheo=ListeRetour[[2]]
QuantTheo      # Theoretical quantile to compare with

#### end ####


### Start : for simulations do not execute this section 

#### For real Data #### 
# Step 6 : choose the real data set

# Besecura article Weissman
require(CASdatasets)
data("besecura")
n=length(besecura$Loss)
Reference=matrix(0,nr=N,nc=n)
Reference[N,]=besecura$Loss

# Uncomment lines 130 to 134 if you want to execute on flood

# Flood article Expected Shortfall
# library(Expectrem)
# data("flood_data")
# n=length(flood_data$Area_A_2012)
# Reference=matrix(0,nr=N,nc=n)
# Reference[N,]=flood_data$Area_A_2012

dn=(2:(n-1))/(n*pn)

# Step 6 bis: be vigilant on the level of the quantile to compare with
QuantTheo=max(Reference) # if pn=1/n
# QuantTheo=quantile(Reference[1,],1-pn) # if pn=0.005
CTETheo=mean(Reference[Reference>QuantTheo])
CTETheo

### End

############################### ABC procedure ###############################

# Caution for these two lists which go from 1 to n-2 corresponds to k which goes from 2 to n-1

# Arguments of the ABC function 

# RefCour : reference sample for comparison
# alphan : level of the quantile 
# n : fixed size of samples to be drawn
# M:  number of samples to be drawn for a fixed size
# tau : percentage of samples to be retained
# sigmaR : standard deviation of the a priori gamma distribution on rho
# sigmaB : standard deviation of gamma distribution a priori on beta

############### Loop over the number N of replications ###############################

MatQGomes=matrix(0,nr=N,nc=n-2) # Matrix initialization of N vectors of k QGomez estimates by Gomes
MatQWeiss=matrix(0,nr=N,nc=n-2) # Matrix initialization of N vectors of k QWeiss estimates by Weiss 
MatQPWM=matrix(0,nr=N,nc=n-2)   # Matrix initialization of N vectors of k QPWM estimates by PWM
MatQABC=matrix(0,nr=N,nc=n-2)   # Matrix initialization of N vectors of k QABC estimates by ABC
MatCTEABC=matrix(0,nr=N,nc=n-2) # Matrix initialization of N vectors of k CTEABC estimates by ABC

for(Nrep in 1:N)# Start loop on number N of replications
{ 
  LseuilBis=NULL     # Reset the vector of thresholds
  LexcesBis=list()   # Rest the list of excesses 
  LabcBis=list()     # Reset the dataframe which contains Gamma, Rho, Beta, Delta, QuantileABC
  EstimHill=NULL     # Reset the vector of Hill estimations
  EstimGomes=NULL    # Reset the vector of Gomes estimations
  EstimPWM0=NULL     # Reset the vector of PWM0 estimations 
  EstimPWM1=NULL     # Reset the vector of PWM1 estimations 
  EstimPWM2=NULL     # Reset the vector of PWM2 estimations 
  EstimR=NULL        # Reset the vector of moments ratio estimations 
  GammaPWM=NULL      # Reset the vector of GammaPWM estimations
  BetaPWM=NULL       # Reset the vector of BetaPWM estimations
  RhoPWM=NULL        # Reset the vector of RhoPWM estimations
  QGomes=NULL        # Reset the vector of QGomes estimations
  QWeiss=NULL        # Reset the vector of QWeiss estimations
  QPWM=NULL          # Reset the vector of QPWM estimations
  MeanQbis=NULL      # Reinitialiser le vecteur des n-2 estimations des moyennes des quantiles sans les seuils
  MeanLogQbis=NULL   # Reinitialiser le vecteur des n-2 estimations des moyennes des logquantiles avec les seuils
  
  RefCour=Reference[Nrep,] # Reference Courant de numero : Nrep de 1 a N


# Estimation of Gamma for Hill and CH 
RhoChap=mop(RefCour,p=0,k=1,"RBMOP")$rho   # Estimation of rho 
BetaChap=mop(RefCour,p=0,k=1,"RBMOP")$beta # Estimation of beta 
EstimHill=sapply(2:(n-1), function(i) GammaHill(RefCour,i))          # Hill estimation
EstimGomes=EstimHill*(1-((BetaChap*(dn*pn)^(-RhoChap))/(1-RhoChap))) # CH estimation

  
# Calcul les k (de 2 a n-1) seuils, les exces correspondants et ABC sur les exces pour RefCour = Reference Courant (1 echantillon)
LseuilBis=sapply(2:(n-1), function(i) Seuil(RefCour,i))
# LseuilBis : vecteur des n-2 seuils
LexcesBis=lapply(2:(n-1), function(i) Exces(RefCour,i))
# LexcesBis : liste de n-2  vecteurs  de taille croissante de 2 a n-1
LabcBis=lapply(2:(n-1), function(i) FonctABC(LexcesBis[[i-1]],taille=n,M=2000,tau=0.05,sigmaR=1,sigmaB=1,alpha=(n*pn/i),gammagomes=EstimGomes,rhogomes=RhoChap,betagomes=BetaChap))
# LabcBis : liste de n-2  dataframe Gamma, Rho, Beta, Delta, QuantileABC



# PWM
EstimPWM0=sapply(2:(n-1), function(i) FonctPWM(a=0,RefCour,i))   # Calcul des n-2 estimations des moments PWM0 (meme que Hill)
EstimPWM1=sapply(2:(n-1), function(i) FonctPWM(a=1,RefCour,i))   # Calcul des n-2 estimations des moments PWM1
EstimPWM2=sapply(2:(n-1), function(i) FonctPWM(a=2,RefCour,i))   # Calcul des n-2 estimations des moments PWM2
EstimR=pmin(pmax((EstimPWM2-EstimPWM0)/(EstimPWM1-EstimPWM0),42/31),1.99)   # Calcul du ratio des moments

#GammaPWM2=EstimPWM0*(1-RhoPWM)/(1-RhoPWM+BetaPWM)
GammaPWM=EstimPWM0+((EstimPWM2-EstimPWM0)/(4-3*EstimR))   # Calcul des n-2 estimations de gamma par PWM
# MatGammaPWM[Nrep,]=GammaPWM  # Stoker le vecteur des estimations de PWM pour l echantillon courant
RhoPWM=3+(2/(EstimR-2))   # Calcul des n-2 estimations de rho par RhoPWM
# MatRhoPWM[Nrep,]=RhoPWM  # Stoker le vecteur des estimations de Rho par PWM pour l echantillon courant
BetaPWM=(2*(1-EstimR)/(EstimR-2))*((EstimPWM0-EstimPWM2)/(EstimPWM2+3*EstimPWM0*(1-EstimR)))   
# Calcul des n-2 estimations de BetaPWM
# MatBetaPWM[Nrep,]=BetaPWM  # Stoker le vecteur des estimations de Beta par PWM pour l echantillon courant

# Quantile Weissman, Gomes et PWM 
QGomes=log(LseuilBis)+EstimGomes*(log(dn))+EstimGomes*(BetaChap*((dn*pn)^(-RhoChap)))*((dn^(RhoChap)-1)/RhoChap)   
# Calcul des n-2 estimations des log quantiles par Gomes
MatQGomes[Nrep,]=QGomes  # Stoker le vecteur des estimations de Q par Gomes pour l echantillon courant
QWeiss=log(LseuilBis)+EstimHill*(log(dn))   # Calcul des n-2 estimations des log quantiles par QWeiss

# IC
QWeissIC1=QWeiss-EstimHill*(log(dn))*(qnorm(0.95)/sqrt(2:(n-1)))
QWeissIC2=QWeiss+EstimHill*(log(dn))*(qnorm(0.95)/sqrt(2:(n-1)))

b1=((qnorm(0.95)/sqrt(2:(n-1)))-BetaChap*(n/2:(n-1))^(RhoChap)/(1-RhoChap))
b2=((qnorm(0.95)/sqrt(2:(n-1)))+BetaChap*(n/2:(n-1))^(RhoChap)/(1-RhoChap))

LCLG=EstimGomes/(1+b2)
UCLG=EstimGomes/(1-b1)

QGomesIC1=QGomes-LCLG*(log(dn))*((qnorm(0.95)/sqrt(2:(n-1))))   
QGomesIC2=QGomes+UCLG*(log(dn))*((qnorm(0.95)/sqrt(2:(n-1))))   
# Fin IC

MatQWeiss[Nrep,]=QWeiss # Stoker le vecteur des estimations de Q par Weiss pour l echantillon courant
QPWM=log(LseuilBis)+GammaPWM*(log(dn))+GammaPWM*BetaPWM*((dn)^(RhoPWM)-1)/RhoPWM   # Calcul des log quantiles par PWM
MatQPWM[Nrep,]=QPWM  # Stoker le vecteur des estimations de Q par PWM pour l echantillon courant
  
# Boucle les moyennes des Quantiles, Gamma, Rho et Beta pour les resultats ABC sur les exces
MeanQbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$QuantileABC))))  
MeanCTE=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$CTEABC))))

# Calcul les n-2 moyennes des log quantiles par ABC
IC5=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$QuantileABC,0.05)))))   # Calcul les n-2 bornes inf
IC95=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$QuantileABC,0.95)))))   # Calcul les n-2 bornes sup

Cte5=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$CTEABC,0.05)))))   # Calcul les n-2 bornes inf
Cte95=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$CTEABC,0.95)))))   # Calcul les n-2 bornes sup

MeanGbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$Gamma))))   # Calcul les n-2 moyennes de Gamma
MeanRbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$Rho))))     # Calcul les n-2 moyennes de Rho
MeanBbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$Beta))))    # Calcul les n-2 moyennes de Beta
DataFrameBis=data.frame(MoyenneGamma=MeanGbis,MoyenneRho=MeanRbis,MoyenneBeta=MeanBbis) # Constitue le DataFrame des 3 vecteurs des moyennes precedentes

# Ajout du log du seuil au log de la moyenne des quantiles
MeanLogQbis=log(LseuilBis)+MeanQbis   # Calcul du log quantile final en ajoutant le seuil
MeanLogCTE=log(LseuilBis)+MeanCTE   # Calcul du log CTE final en ajoutant le seuil

MatQABC[Nrep,]=MeanLogQbis            # Stoker le vecteur des estimations de Q par ABC pour l echantillon courant
MatCTEABC[Nrep,]=MeanLogCTE            # Stoker le vecteur des estimations de CTE par ABC pour l echantillon courant

IC5seuil=log(LseuilBis)+IC5           # Calcul IC 5% du log quantile final en ajoutant le seuil
IC95seuil=log(LseuilBis)+IC95         # Calcul IC 95% du log quantile final en ajoutant le seuil

Cte5seuil=log(LseuilBis)+Cte5           # Calcul IC 5% du log quantile final en ajoutant le seuil
Cte95seuil=log(LseuilBis)+Cte95         # Calcul IC 95% du log quantile final en ajoutant le seuil

DataFrameFinal=data.frame(DataFrameBis,MoyenneLogQuantile=MeanLogQbis,IC005=IC5seuil,IC095=IC95seuil,MoyenneLogCTE=MeanLogCTE,CTE005=Cte5seuil,CTE095=Cte95seuil) # Complete le dataframe precedent par ces 3 vecteurs

} ############### FIN boucle sur le nombre N de replications (d echantillons) ###############################

################################################
################ FOR SIMULATIONS ###############
################################################

# ___________________________ Sauvegarde __________________________________
# ATTENTION : Changer de repertoire ou changer de nom

# Enregistrement de MatQGomes : tableau des estimations des quantiles par Gomes
write.table(MatQGomes,"MatQGomes", row.names = FALSE, col.names = FALSE)

# Enregistrement de MatQWeiss : tableau des estimations des quantiles par Weiss
write.table(MatQWeiss,"MatQWeiss", row.names = FALSE, col.names = FALSE)

# Enregistrement de MatQPWM : tableau des estimations des quantiles par PWM
write.table(MatQPWM,"MatQPWM", row.names = FALSE, col.names = FALSE)

# Enregistrement de MatQABC : tableau des estimations des quantiles par ABC
write.table(MatQABC,"MatQABC", row.names = FALSE, col.names = FALSE)


###########################################################################
# ______________________________   MSE   __________________________________

MatMse=matrix(0,nr=(n-2),nc=4)  # matrice des 4 vecteurs des estimations des 4 MSE 

##########  Lecture des tableaux
# Choisir le repertoire du travail courant
# setwd(dirname(file.choose()))
# Pour verification
# getwd()


# Lecture de MatQWeiss : tableau des estimations des quantiles par Weiss
MatQWeiss=read.table("MatQWeiss", header = FALSE)
# Lecture de MatQABC : tableau des estimations des quantiles par ABC
MatQABC=read.table("MatQABC", header = FALSE)
# Lecture de MatQGomes : tableau des estimations des quantiles par Gomes
MatQGomes=read.table("MatQGomes", header = FALSE)
# Lecture de MatQPWM : tableau des estimations des quantiles par PWM
MatQPWM=read.table("MatQPWM", header = FALSE)


# Weissman
BiaisWeiss=colMeans(MatQWeiss-log(QuantTheo))
Biais2Weiss=(colMeans(MatQWeiss-log(QuantTheo)))^2
VarWeiss=colVars(MatQWeiss-log(QuantTheo))
MseWeiss=Biais2Weiss+VarWeiss
MatMse[,1]=MseWeiss
# ABC 
BiaisABC=colMeans(MatQABC-log(QuantTheo))
Biais2ABC=(colMeans(MatQABC-log(QuantTheo)))^2
VarABC=colVars(MatQABC-log(QuantTheo))
MseABC=Biais2ABC+VarABC
MatMse[,2]=MseABC
# Gomes 
BiaisGomes=colMeans(MatQGomes-log(QuantTheo))
Biais2Gomes=(BiaisGomes)^2
VarGomes=colVars(MatQGomes-log(QuantTheo))
MseGomes=Biais2Gomes+VarGomes
MatMse[,3]=MseGomes
# PWM 
BiaisPWM=colMeans(MatQPWM-log(QuantTheo))
Biais2PWM=(colMeans(MatQPWM-log(QuantTheo)))^2
VarPWM=colVars(MatQPWM-log(QuantTheo))
MsePWM=Biais2PWM+VarPWM
MatMse[,4]=MsePWM


################ Graphiques Simulations ################

MSE=data.frame(MSEABC=MseABC,MSEWeiss=MseWeiss,MSEGomes=MseGomes,MSEPWM=MsePWM,BiasABC=BiaisABC,BiasWeiss=BiaisWeiss,BiasGomes=BiaisGomes,BiasPWM=BiaisPWM)
library("ggplot2")
xValue=1:(n-2)
DataFrame=data.frame(MSE,xValue)


GraphMSE=ggplot(DataFrame,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(0,1)
MseFinal=GraphMSE+geom_line(aes(y=MSEABC),col="black",lwd=1.1)+geom_line(aes(y=MSEWeiss),col="blue",lwd=1.1) +geom_line(aes(y=MSEGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=MSEPWM),col="Purple",lwd=1.1)+theme(axis.text=element_text(size=30))                
MseFinal

ggsave("MseStudentG12R1.pdf",MseFinal)

GraphBiais=ggplot(DataFrame,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(-0.5,1.5)
BiasFinal=GraphBiais+geom_line(aes(y=BiasABC),col="black",lwd=1.1)+geom_line(aes(y=BiasWeiss),col="blue",lwd=1.1) +geom_line(aes(y=BiasGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=BiasPWM),col="Purple",lwd=1.1)+ geom_hline(yintercept=0,col="red")+theme(axis.text=element_text(size=30))               
BiasFinal

ggsave("BiasStudentG12R1.pdf",BiasFinal)

IC=data.frame(IC05=DataFrameFinal$IC005,IC95=DataFrameFinal$IC095,LogQabc=DataFrameFinal$MoyenneLogQuantile)
xValue=1:(n-2)
DataFrameIC=data.frame(IC,xValue)

GraphIC=ggplot(DataFrameIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(0,13)#min(DataFrameIC$IC05),max(DataFrameIC$IC95))
ICFinal=GraphIC+geom_line(aes(y=LogQabc),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICStudentG12R1.pdf",ICFinal)


# Enregistrement de MatMse : tableau des estimations des 4 MSE
 write.table(MatMse,"MatMse", row.names = FALSE, col.names = FALSE)


# Calcul du min MSE apres calcul des matrice dans le cas de la loi du second ordre
# Calcul des moyennes autour du min MSE
# Ecart=5    # Amplitude autour du min MSE pour calculer la moyenne
# Calcul du min de la MSE pour ABC 
MinMSEABC=min(MseABC)
# Calcul de la moyenne de la MSE autour du min pour ABC 
IndMinMseABC=which.min(MseABC)     # Indice du minimum
# MoyMSEABC=mean(MseABC[(IndMinMseABC-Ecart):(IndMinMseABC+Ecart)])

# Calcul du min de la MSE pour Weissman
MinMSEWeiss=min(MseWeiss)
# Calcul de la moyenne de la MSE autour du min pour Weissman
IndMinMseWeiss=which.min(MseWeiss)     # Indice du minimum
# MoyMSEWeiss=mean(MseWeiss[(IndMinMseWeiss-Ecart):(IndMinMseWeiss+Ecart)])
# Calcul du min de la MSE pour Gomes 
MinMSEGomes=min(MseGomes)
# Calcul de la moyenne de la MSE autour du min pour Gomes
IndMinMseGomes=which.min(MseGomes)     # Indice du minimum
# MoyMSEGomes=mean(MseGomes[(IndMinMseGomes-Ecart):(IndMinMseGomes+Ecart)])
# Calcul du min de la MSE pour PWM 
MinMSEPWM=min(MsePWM)
# Calcul de la moyenne de la MSE autour du min pour PWM
IndMinMsePWM=which.min(MsePWM)     # Indice du minimum
# MoyMSEPWM=mean(MsePWM[(IndMinMsePWM-Ecart):(IndMinMsePWM+Ecart)])

# Affichage des min des ’MSE
cat( "Weiss : ", round(MinMSEWeiss,4), "ABC : ", round(MinMSEABC,4), "Gomes : ", round(MinMSEGomes,4), "PWM : ", round(MinMSEPWM,4) )
# Affichage des indices des min des MSE
cat( "Weiss : ", IndMinMseWeiss,"ABC : ", IndMinMseABC, "Gomes : ", IndMinMseGomes, "PWM : ", IndMinMsePWM )
# Affichage des min des MSE pour verification
cat( "Weiss : ", MseWeiss[IndMinMseWeiss],"ABC : ", MseABC[IndMinMseABC], "Gomes : ", MseGomes[IndMinMseGomes], "PWM : ", MsePWM[IndMinMsePWM] )


###########################################################################

################################################
################ FOR REAL DATA  ################
################################################

######### Besecura  ylim(15.5,17.2) #########

DataHist=data.frame(ref=as.numeric(Reference))
Histo=ggplot(DataHist, aes(x=ref)) +theme_light()+labs(x=" ")+labs(y=" ") + geom_histogram(fill = "blue",color = "blue")+theme(axis.text=element_text(size=30))
Histo

ggsave("HistoBesecura.pdf",Histo)


xValue=1:(n-2)
DataFrameFinalIC=data.frame(xValue,LogQuantile=DataFrameFinal$MoyenneLogQuantile,IC05=DataFrameFinal$IC005,IC95=DataFrameFinal$IC095)

# IC Log Quantile ABC
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICbesecura.pdf",ICFinal)


DataFrameFinalQ=data.frame(DataFrameFinalIC,QGomes,QWeiss,QPWM)
# Log Quantile
QuantileS=ggplot(DataFrameFinalQ,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.5)#+ylim(min(DataFrameFinalQ$IC05),max(DataFrameFinalQ$IC95))
QFinal=QuantileS+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_line(aes(y=QGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=QWeiss),col="blue",lwd=1.1)+geom_line(aes(y=QPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
QFinal

ggsave("LogQbesecura.pdf",QFinal)


# linetype="dashed"
DataFrameFinalG=data.frame(DataFrameFinalQ,EstimGomes,EstimHill,GammaPWM,GammaABC=DataFrameFinal$MoyenneGamma)
# Gamma
GammaS=ggplot(DataFrameFinalG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(0,0.5)
GFinal=GammaS+geom_line(aes(y=GammaABC),col="black",lwd=1.1)+geom_line(aes(y=EstimGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=EstimHill),col="blue",lwd=1.1)+geom_line(aes(y=GammaPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
GFinal

ggsave("Gbesecura.pdf",GFinal)


# IC Log Quantile W
DataFrameFinalICW=data.frame(xValue,QWeiss,IC05Q=QWeissIC1,IC95Q=QWeissIC2)
ICLogQuantileW=ggplot(DataFrameFinalICW,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.5)#ylim(min(DataFrameFinalICW$IC05Q),max(DataFrameFinalICW$IC95Q))
ICFinalW=ICLogQuantileW+geom_line(aes(y=QWeiss),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05Q, ymax=IC95Q), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05Q),col="blue")+geom_line(aes(y=IC95Q),col="blue")+theme(axis.text=element_text(size=30))
ICFinalW

ggsave("ICbesecuraW.pdf",ICFinalW)


# IC Log Quantile G
DataFrameFinalICG=data.frame(xValue,QGomes,IC05G=QGomesIC1,IC95G=QGomesIC2)
ICLogQuantileG=ggplot(DataFrameFinalICG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.2)#ylim(min(DataFrameFinalICG$IC05G),max(DataFrameFinalICG$IC95G))
ICFinalG=ICLogQuantileG+geom_line(aes(y=QGomes),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05G, ymax=IC95G), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05G),col="blue")+geom_line(aes(y=IC95G),col="blue")+theme(axis.text=element_text(size=30))
ICFinalG

ggsave("ICbesecuraG.pdf",ICFinalG)

# IC Log Quantile ABC
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.2)#ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICbesecuraC.pdf",ICFinal)

xValue=1:(n-2)
DataFrameFinalCTE=data.frame(xValue,LogCTE=DataFrameFinal$MoyenneLogCTE,ICcte05=DataFrameFinal$CTE005,ICcte95=DataFrameFinal$CTE095)

# IC Log CTE ABC
ICLogCTE=ggplot(DataFrameFinalCTE,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.2)
ICcte=ICLogCTE+geom_line(aes(y=LogCTE),col="black",lwd=1.1)+theme(axis.text=element_text(size=30))+geom_ribbon(aes(ymin=ICcte05, ymax=ICcte95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=ICcte05),col="blue")+geom_line(aes(y=ICcte95),col="blue")+geom_hline(yintercept=log(CTETheo),col="red")
ICcte

ggsave("CTEbesecura.pdf",ICcte)




######### Floodylim(11,24) #########


DataHist=data.frame(ref=as.numeric(Reference))
Histo=ggplot(DataHist, aes(x=ref)) +theme_light()+labs(x=" ")+labs(y=" ") + geom_histogram(fill = "blue",color = "blue")+theme(axis.text=element_text(size=30))
Histo

ggsave("HistoFlood.pdf",Histo)

xValue=1:(n-2)
xValue=1:((n-1)/2)
DataFrameFinalIC=data.frame(xValue,LogQuantile=DataFrameFinal$MoyenneLogQuantile,IC05=DataFrameFinal$IC005,IC95=DataFrameFinal$IC095)

# IC Log Quantile ABC
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICflood.pdf",ICFinal)

 
DataFrameFinalQ=data.frame(DataFrameFinalIC,QGomes,QWeiss,QPWM)
# Log Quantile
QuantileS=ggplot(DataFrameFinalQ,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#+ylim(min(DataFrameFinalQ$IC05),max(DataFrameFinalQ$IC95))
QFinal=QuantileS+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_line(aes(y=QGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=QWeiss),col="blue",lwd=1.1)+geom_line(aes(y=QPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
QFinal

ggsave("LogQflood.pdf",QFinal)


# linetype="dashed"
DataFrameFinalG=data.frame(DataFrameFinalQ,EstimGomes,EstimHill,GammaPWM,GammaABC=DataFrameFinal$MoyenneGamma)
# Gamma
GammaS=ggplot(DataFrameFinalG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(0,1)
GFinal=GammaS+geom_line(aes(y=GammaABC),col="black",lwd=1.1)+geom_line(aes(y=EstimGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=EstimHill),col="blue",lwd=1.1)+geom_line(aes(y=GammaPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
GFinal

ggsave("Gflood.pdf",GFinal)


# IC Log Quantile W
DataFrameFinalICW=data.frame(xValue,QWeiss,IC05Q=QWeissIC1,IC95Q=QWeissIC2)
ICLogQuantileW=ggplot(DataFrameFinalICW,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#ylim(min(DataFrameFinalICW$IC05Q),max(DataFrameFinalICW$IC95Q))
ICFinalW=ICLogQuantileW+geom_line(aes(y=QWeiss),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05Q, ymax=IC95Q), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05Q),col="blue")+geom_line(aes(y=IC95Q),col="blue")+theme(axis.text=element_text(size=30))
ICFinalW

ggsave("ICfloodW.pdf",ICFinalW)


# IC Log Quantile G
DataFrameFinalICG=data.frame(xValue,QGomes,IC05G=QGomesIC1,IC95G=QGomesIC2)
ICLogQuantileG=ggplot(DataFrameFinalICG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#ylim(min(DataFrameFinalICG$IC05G),max(DataFrameFinalICG$IC95G))
ICFinalG=ICLogQuantileG+geom_line(aes(y=QGomes),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05G, ymax=IC95G), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05G),col="blue")+geom_line(aes(y=IC95G),col="blue")+theme(axis.text=element_text(size=30))
ICFinalG

ggsave("ICfloodG.pdf",ICFinalG)

# IC Log Quantile ABC
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICfloodC.pdf",ICFinal)

xValue=1:(n-2)
DataFrameFinalCTE=data.frame(xValue,LogCTE=DataFrameFinal$MoyenneLogCTE,ICcte05=DataFrameFinal$CTE005,ICcte95=DataFrameFinal$CTE095)

# IC Log CTE ABC
ICLogCTE=ggplot(DataFrameFinalCTE,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)
ICcte=ICLogCTE+geom_line(aes(y=LogCTE),col="black",lwd=1.1)+theme(axis.text=element_text(size=30))+geom_ribbon(aes(ymin=ICcte05, ymax=ICcte95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=ICcte05),col="blue")+geom_line(aes(y=ICcte95),col="blue")+geom_hline(yintercept=log(CTETheo),col="red")
ICcte

ggsave("CTEflood.pdf",ICcte)

