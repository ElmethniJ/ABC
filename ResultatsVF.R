# Effacer tout
rm(list = ls())
gc()

# Packages utiles
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

# Module de creation des echantillons de reference 
source("CreationEchantillonsVF.R")

# Fonction ABC
source("FonctionABCvf.R")  # Nouvelle version VF

# Taille de l echantillon
n=500 # simulations
# n=371 # besecura
# n= 555 # flood

# Nombre de replications
N=1
  
# Ordre du quantile 
pn=1/(2*n)  # Ordre du quantile
# pn=1/n  # Ordre du quantile
# pn=0.005    # Ordre du quantile

# Facteur dextrapolation
dn=(2:(n-1))/(n*pn)

# Fixe la seed
set.seed(12)

# Choisir une loi parmi les 6 en initialisant la variable : ChoixLoi  par un nombre entier de 1 a 7
ChoixLoi=4	# Switch choix loi:  1=Loi du second ordre, 2=Burr,  3=Frechet, 4=Student, 5=Fisher, 6 InvGamma, 7=GPD

# Parametres Loi du second ordre  
gammaSecond=1/2  
RhoSecond=-1/2  
betaSecond=1

# Parametres Loi de Burr
gammaB=1/2    # Loi de Burr
RhoB=-1/2  # rho, Loi de Burr

# Parametres Loi de Frechet
gammaF=1/2       # Loi de Frechet implique Rho
RhoF=-1           # rho, Loi de Frechet

# Parametres Loi de Student
NuS=2            # Ddl de la Loi de Student
gammaS=1/NuS     # Loi de Student
RhoS=-2/NuS      # rho, Loi de Student

# Parametres Loi de Fisher
NuF1=3            # Ddl1 de la Loi de Fisher
NuF2=4            # Ddl2 de la Loi de Fisher
gammaFi=2/NuF2        # Loi de Fisher
RhoFi=-2/NuF2       # rho, Loi de Fisher

# Parametres Loi InvGamma
BetaInvg=1      # ou 2, arbitraire ! determine le rate
AlphaInvg=2     # AlphaInvg=-1/RhoInvg
RhoInvg=-1/AlphaInvg
gammaInvg=1/AlphaInvg

# Parametres Loi GPD
RhoG=-1/2        # rho GPD
gammaG=-RhoG       # Loi GPD

LoiChoisie=function(choixLoi)
{
  switch(choixLoi,"Second Ordre ","Burr ", "Frechet ", "Student","Fisher", "InvGamma","GPD")
}

loi=LoiChoisie(ChoixLoi)	# loi choisie

ChoixRho=function(choixLoi)
{
  switch(choixLoi,RhoSecond,RhoB, RhoF, RhoS,RhoFi, RhoInvg,RhoG)
}

rho=ChoixRho(ChoixLoi)	# rho choisi

ChoixGamma=function(choixLoi)
{
  switch(choixLoi,gammaSecond,gammaB, gammaF, gammaS,gammaFi, gammaInvg,gammaG)
}

gamma=ChoixGamma(ChoixLoi)	# gamma choisi

# Matrice des N echantillons de taille n
ListeRetour=FonctCreerEch(ChoixLoi)
Reference=ListeRetour[[1]] 
# Reference      # Matrice de N lignes (Nb echantillons) et n colonnes (taille echantillons)
QuantTheo=ListeRetour[[2]]
QuantTheo      # Quantile theorique auquel se comparer


# Real Data pour tester 
# Ne pas oublier de changer le n avant

# Besecura article Weissman
# require(CASdatasets)
# data("besecura")
# n=length(besecura$Loss)
# Reference=matrix(0,nr=N,nc=n)	# matrice des N(=1) echantillons de taille n
# Reference[N,]=besecura$Loss

# Flood article Expected Shortfall
# library(Expectrem)
# data("flood_data")
# n=length(flood_data$Area_A_2012)
# Reference=matrix(0,nr=N,nc=n)	# matrice des N(=1) echantillons de taille n
# Reference[N,]=flood_data$Area_A_2012

# dn=(2:(n-1))/(n*pn)
# QuantTheo=max(Reference) 
# QuantTheo=quantile(Reference[1,],1-pn)
# CTETheo=mean(Reference[Reference>QuantTheo])
# CTETheo

############################### Procedure ###############################
 # _______________________________________________________________________
# Procedure pour choisir le quantile theorique selon une loi choisie

# Attention pour ces deux listes qui vont de 1 a n-2 correspond a k qui va de  2 a n-1

# Arguments Fonction ABC
# RefCour : echantillon de reference auquel se comparer
# alphan : ordre du quantile 
# M:  Nombre d'echantillons a tirer pour une taille fixee
# tau : Pourcentage d'echantillons a garder
# sigmaR : ecart-type de la loi gamma a priori sur rho
# sigmaB : ecart-type de la loi gamma a priori sur beta

############### Boucle sur le nombre N de replications (d echantillons) ###############################

# MatGamHill=matrix(0,nr=N,nc=n-2)     # matrice des N vecteurs des estimations des k Gamma par Hill
# MatGamGomes=matrix(0,nr=N,nc=n-2)    # matrice des N vecteurs des estimations des k Gamma par Gomes
# MatGamPWM0=matrix(0,nr=N,nc=n-2)     # matrice des N vecteurs des estimations des k Gamma par PWM0
# MatGamPWM1=matrix(0,nr=N,nc=n-2)     # matrice des N vecteurs des estimations des k Gamma par PWM1
# MatGamPWM2=matrix(0,nr=N,nc=n-2)     # matrice des N vecteurs des estimations des k Gamma par PWM2
# MatGamR=matrix(0,nr=N,nc=n-2)        # matrice des N vecteurs des estimations des k Gamma corrigees ??
# MatGammaPWM=matrix(0,nr=N,nc=n-2)    # matrice des N vecteurs des estimations des k Gamma par PWM
# MatBetaPWM=matrix(0,nr=N,nc=n-2)     # matrice des N vecteurs des estimations des k Beta par PWM
# MatRhoPWM=matrix(0,nr=N,nc=n-2)      # matrice des N vecteurs des estimations des k Rho par PWM
MatQGomes=matrix(0,nr=N,nc=n-2)        # matrice des N vecteurs des estimations des k QGomez par Gomes
MatQWeiss=matrix(0,nr=N,nc=n-2)        # matrice des N vecteurs des estimations des k QWeiss par Weiss
MatQPWM=matrix(0,nr=N,nc=n-2)          # matrice des N vecteurs des estimations des k QGomez par PWM
MatQABC=matrix(0,nr=N,nc=n-2)          # matrice des N vecteurs des estimations des k QABC par ABC
MatCTEABC=matrix(0,nr=N,nc=n-2)          # matrice des N vecteurs des estimations des k QABC par ABC

for(Nrep in 1:N)	# Debut boucle sur le nombre N de replications
{ 
  LseuilBis=NULL     # Reinitialiser le vecteur des n-2 seuils 
  LexcesBis=list()   # Reinitialiser la liste de n-2  vecteurs  de taille croissante de 2 a n-1 
  LabcBis=list()     # Reinitialiser la liste de n-2  dataframe Gamma, Rho, Beta, Delta, QuantileABC
  EstimHill=NULL     # Reinitialiser le vecteur des n-2 estimations de Hill
  EstimGomes=NULL    # Reinitialiser le vecteur des n-2 estimations de Gomes
  EstimPWM0=NULL     # Reinitialiser le vecteur des n-2 estimations de PWM0  
  EstimPWM1=NULL     # Reinitialiser le vecteur des n-2 estimations de PWM1 
  EstimPWM2=NULL     # Reinitialiser le vecteur des n-2 estimations de PWM2 
  EstimR=NULL        # Reinitialiser le vecteur des n-2 estimations du ratio des moments
  GammaPWM=NULL      # Reinitialiser le vecteur des n-2 estimations de GammaPWM
  BetaPWM=NULL       # Reinitialiser le vecteur des n-2 estimations de BetaPWM 
  RhoPWM=NULL        # Reinitialiser le vecteur des n-2 estimations de RhoPWM
  QGomes=NULL        # Reinitialiser le vecteur des n-2 estimations de QGomes
  QWeiss=NULL        # Reinitialiser le vecteur des n-2 estimations de QWeiss
  QPWM=NULL          # Reinitialiser le vecteur des n-2 estimations de QPWM
  MeanQbis=NULL      # Reinitialiser le vecteur des n-2 estimations des moyennes des quantiles sans les seuils
  MeanLogQbis=NULL   # Reinitialiser le vecteur des n-2 estimations des moyennes des logquantiles avec les seuils
  
  RefCour=Reference[Nrep,]    # Reference Courant de numero : Nrep de 1 a N

  
# Gamma Hill et Gomes 
# EstimHill=as.vector(unlist(lapply(2:(n-1), function(i) mean(log(LexcesBis[[i-1]])))))
RhoChap=mop(RefCour,p=0,k=1,"RBMOP")$rho    # Calcul RhoChap pour Gomes
BetaChap=mop(RefCour,p=0,k=1,"RBMOP")$beta   # Calcul BetaChap  pour Gomes
EstimHill=sapply(2:(n-1), function(i) GammaHill(RefCour,i))   # Calcul des n-2 estimations de gamma par Hill
# MatGamHill[Nrep,]=EstimHill   # Stoker le vecteur des estimations de gamma par Hill pour l echantillon courant
EstimGomes=EstimHill*(1-((BetaChap*(dn*pn)^(-RhoChap))/(1-RhoChap)))   # Calcul des n-2 estimations de gamma par Gomes
# MatGamGomes[Nrep,]=EstimGomes  # Stoker le vecteur des estimations de gamma par Gomes pour l echantillon courant
  
  
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

