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
require(ggplot2)

# Load reference sample creation program 
source("SampleCreation.R")

# Load ABC function et others useful functions
source("FunctionABC.R")

# Step 1: choose sample size n
# n=500   # to reproduce simulations results
n=371     # for the real data set : besecura
# n= 555  # for the real data set : flood

# Step 2: choose number of replications N
N=1    # for real data set
# N=500  # to reproduce simulations results

# Step 3: choose level of the quantile pn
# pn=1/(2*n)  # to reproduce simulations results
pn=1/n        # for real data set
# pn=0.005    # for real data set

# Extrapolation factor dn
dn=(2:(n-1))/(n*pn)

# Set the seed
set.seed(12)

# Step 4: choose one of the 7 distributions by initializing the variable ChoixLoi
# 1=RPD, 2=Burr, 3=Frechet, 4=Student, 5=Fisher, 6 InvGamma, 7=GPD
ChoixLoi=1	

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


#### Not to be modified and must be executed even for real data ####

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

# Theoretical quantile to compare with
QuantTheo=ListeRetour[[2]]
QuantTheo     

### Start : for simulations do not execute this section 

#### For Real Data #### 
# Step 6: choose the real data set

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

# Step 6 bis: be vigilant on the level of the quantile to compare with
QuantTheo=max(Reference) # if pn=1/n
# QuantTheo=quantile(Reference[1,],1-pn) # if pn=0.005
CTETheo=mean(Reference[Reference>QuantTheo])
CTETheo

### End

############################### ABC procedure ###############################

# Caution for these two lists which goes from 1 to n-2 corresponds to k which goes from 2 to n-1

# Arguments of the ABC function 

# RefCour : reference sample for comparison
# alphan : level of the quantile 
# n : fixed size of samples to be drawn
# M:  number of samples to be drawn for a fixed size
# tau : percentage of samples to be retained
# sigmaR : standard deviation of the a priori gamma distribution on rho
# sigmaB : standard deviation of gamma distribution a priori on beta

############### START : Loop over the number N of replications ###############################

MatQGomes=matrix(0,nr=N,nc=n-2) # Matrix initialization of N vectors of k QGomez estimates by Gomes
MatQWeiss=matrix(0,nr=N,nc=n-2) # Matrix initialization of N vectors of k QWeiss estimates by Weissman
MatQPWM=matrix(0,nr=N,nc=n-2)   # Matrix initialization of N vectors of k QPWM estimates by PWM
MatQABC=matrix(0,nr=N,nc=n-2)   # Matrix initialization of N vectors of k QABC estimates by ABC
MatCTEABC=matrix(0,nr=N,nc=n-2) # Matrix initialization of N vectors of k CTEABC estimates by ABC

for(Nrep in 1:N)     # Start loop on number N of replications
{ 
  LseuilBis=NULL     # Reset the vector of thresholds
  LexcesBis=list()   # Rest the list of excesses 
  LabcBis=list()     # Reset the data frame which contains Gamma, Rho, Beta, Delta, QuantileABC and CTEABC
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
  MeanQbis=NULL      # Reset the vector of the Quantile ABC estimations without threshold
  MeanLogQbis=NULL   # Reset the vector of the log Quantile ABC estimations without threshold
  
  RefCour=Reference[Nrep,] # Current replication : Nrep goes from 1 to N

  # Estimation of Gamma by Hill and CH 
  RhoChap=mop(RefCour,p=0,k=1,"RBMOP")$rho   # Estimation of rho 
  BetaChap=mop(RefCour,p=0,k=1,"RBMOP")$beta # Estimation of beta 
  EstimHill=sapply(2:(n-1), function(i) GammaHill(RefCour,i))          # Hill estimation
  EstimGomes=EstimHill*(1-((BetaChap*(dn*pn)^(-RhoChap))/(1-RhoChap))) # CH estimation
  
  # Compute the k (from 2 to n-1) thresholds
  LseuilBis=sapply(2:(n-1), function(i) Seuil(RefCour,i))
  
  # Compute the k (from 2 to n-1) corresponding excesses
  LexcesBis=lapply(2:(n-1), function(i) Exces(RefCour,i))
  
  # Compute the k (from 2 to n-1) ABC on the excesses
  # LabcBis : list of n-2  data frame Gamma, Rho, Beta, Delta, QuantileABC and CTEABC
  # Be careful is you run on Real Data you have to use EstimGomes instead of EstimHill. The latter is dedicated to simulations results. 
  LabcBis=lapply(2:(n-1), function(i) FonctABC(LexcesBis[[i-1]],taille=n,M=2000,tau=0.05,sigmaR=1,sigmaB=1,alpha=(n*pn/i),gammagomes=EstimGomes,rhogomes=RhoChap,betagomes=BetaChap))
  
  # Estimation of PWM of orders 0, 1 and 2
  EstimPWM0=sapply(2:(n-1), function(i) FonctPWM(a=0,RefCour,i))   # Compute PWM0 (Hill)
  EstimPWM1=sapply(2:(n-1), function(i) FonctPWM(a=1,RefCour,i))   # Compute PWM1
  EstimPWM2=sapply(2:(n-1), function(i) FonctPWM(a=2,RefCour,i))   # Compute PWM2
  EstimR=pmin(pmax((EstimPWM2-EstimPWM0)/(EstimPWM1-EstimPWM0),42/31),1.99)   # Compute moments ratio
  
  # Estimation of Gamma by PWM
  GammaPWM=EstimPWM0+((EstimPWM2-EstimPWM0)/(4-3*EstimR))
  # Estimation of Rho by PWM
  RhoPWM=3+(2/(EstimR-2))
  # Estimation of Beta by PWM
  BetaPWM=(2*(1-EstimR)/(EstimR-2))*((EstimPWM0-EstimPWM2)/(EstimPWM2+3*EstimPWM0*(1-EstimR))) 

  # Estimation of the log Quantile by Gomes 
  QGomes=log(LseuilBis)+EstimGomes*(log(dn))+EstimGomes*(BetaChap*((dn*pn)^(-RhoChap)))*((dn^(RhoChap)-1)/RhoChap)  
  MatQGomes[Nrep,]=QGomes  # Store the vector of estimates by Gomes for the current sample
  
  # Estimation of the log Quantile by Weissman
  QWeiss=log(LseuilBis)+EstimHill*(log(dn)) 
  MatQWeiss[Nrep,]=QWeiss # Store the vector of estimates by Weissman for the current sample
  
  # Start 90% confidence intervals CI
  QWeissIC1=QWeiss-EstimHill*(log(dn))*(qnorm(0.95)/sqrt(2:(n-1)))
  QWeissIC2=QWeiss+EstimHill*(log(dn))*(qnorm(0.95)/sqrt(2:(n-1)))

  b1=((qnorm(0.95)/sqrt(2:(n-1)))-BetaChap*(n/2:(n-1))^(RhoChap)/(1-RhoChap))
  b2=((qnorm(0.95)/sqrt(2:(n-1)))+BetaChap*(n/2:(n-1))^(RhoChap)/(1-RhoChap))

  LCLG=EstimGomes/(1+b2)
  UCLG=EstimGomes/(1-b1)

  QGomesIC1=QGomes-LCLG*(log(dn))*((qnorm(0.95)/sqrt(2:(n-1))))   
  QGomesIC2=QGomes+UCLG*(log(dn))*((qnorm(0.95)/sqrt(2:(n-1))))   
  # End 90% confidence intervals CI

  # Estimation of the log Quantile by PWM 
  QPWM=log(LseuilBis)+GammaPWM*(log(dn))+GammaPWM*BetaPWM*((dn)^(RhoPWM)-1)/RhoPWM
  MatQPWM[Nrep,]=QPWM  # Store the vector of estimates by PWM for the current sample
  
  # Store average log Quantiles for ABC results on exces
  MeanQbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$QuantileABC)))) 
  # Store average log CTE for ABC results on exces
  MeanCTE=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$CTEABC))))

  # Compute 90% credible intervals for log Quantiles obtained thanks to ABC
  IC5=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$QuantileABC,0.05))))) 
  IC95=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$QuantileABC,0.95)))))
  
  # Compute 90% credible intervals for log CTE obtained thanks to ABC
  Cte5=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$CTEABC,0.05))))) 
  Cte95=as.vector(unlist(lapply(2:(n-1), function(i) unname(quantile(LabcBis[[i-1]]$CTEABC,0.95)))))

  # Compute average estimations of Gamma for ABC
  MeanGbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$Gamma)))) 
  # Compute average estimations of Rho for ABC
  MeanRbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$Rho)))) 
  # Compute average estimations of Beta for ABC
  MeanBbis=as.vector(unlist(lapply(2:(n-1), function(i) mean(LabcBis[[i-1]]$Beta))))
  
  # Creation of the DataFrame of the 3 previous average vectors
  DataFrameBis=data.frame(MoyenneGamma=MeanGbis,MoyenneRho=MeanRbis,MoyenneBeta=MeanBbis) 

  # Adding the log of the threshold to the log of the quantile and the CTE estimations
  MeanLogQbis=log(LseuilBis)+MeanQbis
  MeanLogCTE=log(LseuilBis)+MeanCTE   

  # Store the vector of log Quantiles and CTE estimates by ABC for the current sample
  MatQABC[Nrep,]=MeanLogQbis            
  MatCTEABC[Nrep,]=MeanLogCTE 

  # Adding the log of the threshold to 90% credible intervals for log Quantiles obtained thanks to ABC
  IC5seuil=log(LseuilBis)+IC5 
  IC95seuil=log(LseuilBis)+IC95 

  # Adding the log of the threshold to 90% credible intervals for log CTE obtained thanks to ABC
  Cte5seuil=log(LseuilBis)+Cte5   
  Cte95seuil=log(LseuilBis)+Cte95 

  # Creation of the Final Data Frame which complete the previous data frame
  DataFrameFinal=data.frame(DataFrameBis,MoyenneLogQuantile=MeanLogQbis,IC005=IC5seuil,IC095=IC95seuil,MoyenneLogCTE=MeanLogCTE,CTE005=Cte5seuil,CTE095=Cte95seuil) 

} 

############### END : Loop over the number N of replications ###############################





################################################
############## FOR SIMULATIONS ONLY ############
################################################

# Save of Load : it depends on whether you have done the calculations or are loading the matrices.

# SAVE
# MatQWeiss save : table of quantile estimates by Weissman
write.table(MatQWeiss,"MatQWeiss", row.names = FALSE, col.names = FALSE)
# MatQGomes save : table of quantile estimates by Gomes
write.table(MatQGomes,"MatQGomes", row.names = FALSE, col.names = FALSE)
# MatQPWM save : table of quantile estimates by PWM 
write.table(MatQPWM,"MatQPWM", row.names = FALSE, col.names = FALSE)
# MatQABC save : table of quantile estimates by ABC
write.table(MatQABC,"MatQABC", row.names = FALSE, col.names = FALSE)


# LOAD
# ATTENTION : Change directory or change name
# Select current directory
# setwd(dirname(file.choose()))

# Loading MatQWeiss : table of quantile estimates by Weissman
MatQWeiss=read.table("MatQWeiss", header = FALSE)
# Loading MatQGomes : table of quantile estimates by Gomes
MatQGomes=read.table("MatQGomes", header = FALSE)
# Loading MatQPWM : table of quantile estimates by PWM
MatQPWM=read.table("MatQPWM", header = FALSE)
# Loading MatQABC : table of quantile estimates by ABC
MatQABC=read.table("MatQABC", header = FALSE)


############ MSE ############
MatMse=matrix(0,nr=(n-2),nc=4) # Empty matrix of the 4 vectors of the 4 MSE estimates

# Bias and MSE for Weissman
BiaisWeiss=colMeans(MatQWeiss-log(QuantTheo))
Biais2Weiss=(colMeans(MatQWeiss-log(QuantTheo)))^2
VarWeiss=colVars(MatQWeiss-log(QuantTheo))
MseWeiss=Biais2Weiss+VarWeiss
MatMse[,1]=MseWeiss

# Bias and MSE for Gomes 
BiaisGomes=colMeans(MatQGomes-log(QuantTheo))
Biais2Gomes=(BiaisGomes)^2
VarGomes=colVars(MatQGomes-log(QuantTheo))
MseGomes=Biais2Gomes+VarGomes
MatMse[,3]=MseGomes

# Bias and MSE for PWM 
BiaisPWM=colMeans(MatQPWM-log(QuantTheo))
Biais2PWM=(colMeans(MatQPWM-log(QuantTheo)))^2
VarPWM=colVars(MatQPWM-log(QuantTheo))
MsePWM=Biais2PWM+VarPWM
MatMse[,4]=MsePWM

# Bias and MSE for ABC 
BiaisABC=colMeans(MatQABC-log(QuantTheo))
Biais2ABC=(colMeans(MatQABC-log(QuantTheo)))^2
VarABC=colVars(MatQABC-log(QuantTheo))
MseABC=Biais2ABC+VarABC
MatMse[,2]=MseABC

# Save MatMse : table of MSE estimates by Weissman, Gomes, PWM and ABC
write.table(MatMse,"MatMse", row.names = FALSE, col.names = FALSE)


###############################################
#### To reproduce the figures : Simulations ###
###############################################

# MSE as a function of k

MSE=data.frame(MSEABC=MseABC,MSEWeiss=MseWeiss,MSEGomes=MseGomes,MSEPWM=MsePWM,BiasABC=BiaisABC,BiasWeiss=BiaisWeiss,BiasGomes=BiaisGomes,BiasPWM=BiaisPWM)
xValue=1:(n-2)
DataFrame=data.frame(MSE,xValue)

GraphMSE=ggplot(DataFrame,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(0,1)
MseFinal=GraphMSE+geom_line(aes(y=MSEABC),col="black",lwd=1.1)+geom_line(aes(y=MSEWeiss),col="blue",lwd=1.1) +geom_line(aes(y=MSEGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=MSEPWM),col="Purple",lwd=1.1)+theme(axis.text=element_text(size=30))                
MseFinal

ggsave("MseRPDG12R12.pdf",MseFinal) # save figure in pdf format

# Biases as a function of k

GraphBiais=ggplot(DataFrame,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(-0.5,1.5)
BiasFinal=GraphBiais+geom_line(aes(y=BiasABC),col="black",lwd=1.1)+geom_line(aes(y=BiasWeiss),col="blue",lwd=1.1) +geom_line(aes(y=BiasGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=BiasPWM),col="Purple",lwd=1.1)+ geom_hline(yintercept=0,col="red")+theme(axis.text=element_text(size=30))               
BiasFinal

ggsave("BiasRPDG12R12.pdf",BiasFinal) # save figure in pdf format

# Credible Intervals associated with ABC as a function of k

IC=data.frame(IC05=DataFrameFinal$IC005,IC95=DataFrameFinal$IC095,LogQabc=DataFrameFinal$MoyenneLogQuantile)
xValue=1:(n-2)
DataFrameIC=data.frame(IC,xValue)

GraphIC=ggplot(DataFrameIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(0,13)#min(DataFrameIC$IC05),max(DataFrameIC$IC95))
ICFinal=GraphIC+geom_line(aes(y=LogQabc),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICRPDG12R12.pdf",ICFinal) # save figure in pdf format


###############################################
##### To reproduce the tables : Simulations ###
###############################################


# Tables corresponding to MRMSE 1

min(MseWeiss)
min(MseABC)
min(MseGomes)
min(MsePWM)

# Tables corresponding to MRMSE 1

Diff2WeissMRMSE=as.matrix((MatQWeiss-log(QuantTheo))^2)
Diff2AbcMRMSE=as.matrix((MatQABC-log(QuantTheo))^2)
Diff2GomesMRMSE=as.matrix((MatQGomes-log(QuantTheo))^2)
Diff2PwmMRMSE=as.matrix((MatQPWM-log(QuantTheo))^2)

Ecart=0.05*n # neighborhood of the optimal sample fraction 

MinEcatWeiss=matrix(0,nrow=1,ncol=N)
IndEcatWeiss=matrix(0,nrow=1,ncol=N)
for (j in 1:N) {
  IndMinMseWeiss=which.min(Diff2WeissMRMSE[j,1:(0.9*n)])
  IndMin1=max(1,IndMinMseWeiss-Ecart)
  IndMin2=min(n-2,IndMinMseWeiss+Ecart)
  IndEcatWeiss[,j]=IndMinMseWeiss
  MinEcatWeiss[,j]=mean(Diff2WeissMRMSE[j,IndMin1:IndMin2])
}

MinEcatABC=matrix(0,nrow=1,ncol=N)
IndEcatABC=matrix(0,nrow=1,ncol=N)
for (j in 1:N) {
  IndMinMseABC=which.min(Diff2AbcMRMSE[j,1:(0.9*n)])
  IndMin1=max(1,IndMinMseABC-Ecart)
  IndMin2=min(n-2,IndMinMseABC+Ecart)
  IndEcatABC[,j]=IndMinMseABC
  MinEcatABC[,j]=mean(Diff2AbcMRMSE[j,IndMin1:IndMin2])
}

MinEcatGomes=matrix(0,nrow=1,ncol=N)
IndEcatGomes=matrix(0,nrow=1,ncol=N)
for (j in 1:N) {
  IndMinMseGomes=which.min(Diff2GomesMRMSE[j,1:(0.9*n)])
  IndMin1=max(1,IndMinMseGomes-Ecart)
  IndMin2=min(n-2,IndMinMseGomes+Ecart)
  IndEcatGomes[,j]=IndMinMseGomes
  MinEcatGomes[,j]=mean(Diff2GomesMRMSE[j,IndMin1:IndMin2])
}

MinEcatPWM=matrix(0,nrow=1,ncol=N)
IndEcatPWM=matrix(0,nrow=1,ncol=N)
for (j in 1:N) {
  IndMinMsePWM=which.min(Diff2PwmMRMSE[j,1:(0.9*n)])
  IndMin1=max(1,IndMinMsePWM-Ecart)
  IndMin2=min(n-2,IndMinMsePWM+Ecart)
  IndEcatPWM[,j]=IndMinMsePWM
  MinEcatPWM[,j]=mean(Diff2PwmMRMSE[j,IndMin1:IndMin2])
}

# Tables corresponding to MRMSE 2

mean(MinEcatWeiss)
mean(MinEcatABC)
mean(MinEcatGomes)
mean(MinEcatPWM)

############# END of Simualtions part #############




################################################
############## FOR REAL DATA ONLY ##############
################################################

################################################
########### Besecura  ylim(15.5,17.5) ##########
################################################

# Histogram of the dataset
DataHist=data.frame(ref=as.numeric(Reference))
Histo=ggplot(DataHist, aes(x=ref)) +theme_light()+labs(x=" ")+labs(y=" ") + geom_histogram(fill = "blue",color = "blue")+theme(axis.text=element_text(size=30))
Histo

ggsave("HistoBesecura.pdf",Histo) # save figure in pdf format

# Useful data frame for the following
xValue=1:(n-2)
DataFrameFinalIC=data.frame(xValue,LogQuantile=DataFrameFinal$MoyenneLogQuantile,IC05=DataFrameFinal$IC005,IC95=DataFrameFinal$IC095)

# Additional figure not in the paper : plot of the 90% credible intervals of Log Quantile by ABC but not with the same ylim
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICbesecura.pdf",ICFinal) # save figure in pdf format

# Log Quantile estimations of Weissman, Gomes, PWM and ABC 
DataFrameFinalQ=data.frame(DataFrameFinalIC,QGomes,QWeiss,QPWM)
QuantileS=ggplot(DataFrameFinalQ,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.5)#+ylim(min(DataFrameFinalQ$IC05),max(DataFrameFinalQ$IC95))
QFinal=QuantileS+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_line(aes(y=QGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=QWeiss),col="blue",lwd=1.1)+geom_line(aes(y=QPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
QFinal

ggsave("LogQbesecura.pdf",QFinal) # save figure in pdf format

# Plot of the estimations of Gamma
DataFrameFinalG=data.frame(DataFrameFinalQ,EstimGomes,EstimHill,GammaPWM,GammaABC=DataFrameFinal$MoyenneGamma)
GammaS=ggplot(DataFrameFinalG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(0,0.5)
GFinal=GammaS+geom_line(aes(y=GammaABC),col="black",lwd=1.1)+geom_line(aes(y=EstimGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=EstimHill),col="blue",lwd=1.1)+geom_line(aes(y=GammaPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
GFinal

ggsave("Gbesecura.pdf",GFinal) # save figure in pdf format

# Plot of the 90% CI of Log Quantile by Weissman
DataFrameFinalICW=data.frame(xValue,QWeiss,IC05Q=QWeissIC1,IC95Q=QWeissIC2)
ICLogQuantileW=ggplot(DataFrameFinalICW,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.5)#ylim(min(DataFrameFinalICW$IC05Q),max(DataFrameFinalICW$IC95Q))
ICFinalW=ICLogQuantileW+geom_line(aes(y=QWeiss),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05Q, ymax=IC95Q), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05Q),col="blue")+geom_line(aes(y=IC95Q),col="blue")+theme(axis.text=element_text(size=30))
ICFinalW

ggsave("ICbesecuraW.pdf",ICFinalW) # save figure in pdf format

# Plot of the 90% confidence intervals of Log Quantile by Gomes
DataFrameFinalICG=data.frame(xValue,QGomes,IC05G=QGomesIC1,IC95G=QGomesIC2)
ICLogQuantileG=ggplot(DataFrameFinalICG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.5)#ylim(min(DataFrameFinalICG$IC05G),max(DataFrameFinalICG$IC95G))
ICFinalG=ICLogQuantileG+geom_line(aes(y=QGomes),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05G, ymax=IC95G), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05G),col="blue")+geom_line(aes(y=IC95G),col="blue")+theme(axis.text=element_text(size=30))
ICFinalG

ggsave("ICbesecuraG.pdf",ICFinalG) # save figure in pdf format

# Plot of the 90% credible intervals of Log Quantile by ABC
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.5)#ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICbesecuraC.pdf",ICFinal) # save figure in pdf format

# Plot of the 90% credible intervals of Log CTE by ABC
xValue=1:(n-2)
DataFrameFinalCTE=data.frame(xValue,LogCTE=DataFrameFinal$MoyenneLogCTE,ICcte05=DataFrameFinal$CTE005,ICcte95=DataFrameFinal$CTE095)

ICLogCTE=ggplot(DataFrameFinalCTE,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,n-2)+ylim(15.5,17.2)
ICcte=ICLogCTE+geom_line(aes(y=LogCTE),col="black",lwd=1.1)+theme(axis.text=element_text(size=30))+geom_ribbon(aes(ymin=ICcte05, ymax=ICcte95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=ICcte05),col="blue")+geom_line(aes(y=ICcte95),col="blue")+geom_hline(yintercept=log(CTETheo),col="red")
ICcte

ggsave("CTEbesecura.pdf",ICcte) # save figure in pdf format

################################################
##############@ Flood ylim(11,24) ##############
################################################

# Histogram of the dataset
DataHist=data.frame(ref=as.numeric(Reference))
Histo=ggplot(DataHist, aes(x=ref)) +theme_light()+labs(x=" ")+labs(y=" ") + geom_histogram(fill = "blue",color = "blue")+theme(axis.text=element_text(size=30))
Histo

ggsave("HistoFlood.pdf",Histo) # save figure in pdf format

# Useful data frame for the following
xValue=1:((n-1)/2)
DataFrameFinalIC=data.frame(xValue,LogQuantile=DataFrameFinal$MoyenneLogQuantile,IC05=DataFrameFinal$IC005,IC95=DataFrameFinal$IC095)

# Additional figure not in the paper : plot of the 90% credible intervals of Log Quantile by ABC but not with the same ylim
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICflood.pdf",ICFinal) # save figure in pdf format

# Log Quantile estimations of Weissman, Gomes, PWM and ABC 
DataFrameFinalQ=data.frame(DataFrameFinalIC,QGomes,QWeiss,QPWM)
QuantileS=ggplot(DataFrameFinalQ,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#+ylim(min(DataFrameFinalQ$IC05),max(DataFrameFinalQ$IC95))
QFinal=QuantileS+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_line(aes(y=QGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=QWeiss),col="blue",lwd=1.1)+geom_line(aes(y=QPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
QFinal

ggsave("LogQflood.pdf",QFinal) # save figure in pdf format

# Plot of the estimations of Gamma
DataFrameFinalG=data.frame(DataFrameFinalQ,EstimGomes,EstimHill,GammaPWM,GammaABC=DataFrameFinal$MoyenneGamma)
GammaS=ggplot(DataFrameFinalG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(0,1)
GFinal=GammaS+geom_line(aes(y=GammaABC),col="black",lwd=1.1)+geom_line(aes(y=EstimGomes),col="darkgreen",lwd=1.1)+geom_line(aes(y=EstimHill),col="blue",lwd=1.1)+geom_line(aes(y=GammaPWM),col="purple",lwd=1.1)+theme(axis.text=element_text(size=30))
GFinal

ggsave("Gflood.pdf",GFinal) # save figure in pdf format

# Plot of the 90% confidence intervals of Log Quantile by Weissman
DataFrameFinalICW=data.frame(xValue,QWeiss,IC05Q=QWeissIC1,IC95Q=QWeissIC2)
ICLogQuantileW=ggplot(DataFrameFinalICW,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#ylim(min(DataFrameFinalICW$IC05Q),max(DataFrameFinalICW$IC95Q))
ICFinalW=ICLogQuantileW+geom_line(aes(y=QWeiss),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05Q, ymax=IC95Q), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05Q),col="blue")+geom_line(aes(y=IC95Q),col="blue")+theme(axis.text=element_text(size=30))
ICFinalW

ggsave("ICfloodW.pdf",ICFinalW) # save figure in pdf format

# Plot of the 90% confidence intervals of Log Quantile by Gomes
DataFrameFinalICG=data.frame(xValue,QGomes,IC05G=QGomesIC1,IC95G=QGomesIC2)
ICLogQuantileG=ggplot(DataFrameFinalICG,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#ylim(min(DataFrameFinalICG$IC05G),max(DataFrameFinalICG$IC95G))
ICFinalG=ICLogQuantileG+geom_line(aes(y=QGomes),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05G, ymax=IC95G), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05G),col="blue")+geom_line(aes(y=IC95G),col="blue")+theme(axis.text=element_text(size=30))
ICFinalG

ggsave("ICfloodG.pdf",ICFinalG) # save figure in pdf format

# Plot of the 90% credible intervals of Log Quantile by ABC
ICLogQuantile=ggplot(DataFrameFinalIC,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)#ylim(min(DataFrameFinalIC$IC05),max(DataFrameFinalIC$IC95))
ICFinal=ICLogQuantile+geom_line(aes(y=LogQuantile),col="black",lwd=1.1)+geom_hline(yintercept=log(QuantTheo),col="red")+geom_ribbon(aes(ymin=IC05, ymax=IC95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=IC05),col="blue")+geom_line(aes(y=IC95),col="blue")+theme(axis.text=element_text(size=30))
ICFinal

ggsave("ICfloodC.pdf",ICFinal) # save figure in pdf format

# Plot of the 90% credible intervals of Log CTE by ABC
xValue=1:((n-1)/2)
DataFrameFinalCTE=data.frame(xValue,LogCTE=DataFrameFinal$MoyenneLogCTE,ICcte05=DataFrameFinal$CTE005,ICcte95=DataFrameFinal$CTE095)

ICLogCTE=ggplot(DataFrameFinalCTE,aes(x=xValue))+theme_light()+labs(x=" ")+labs(y=" ")+xlim(1,(n-1)/2)+ylim(12,16)
ICcte=ICLogCTE+geom_line(aes(y=LogCTE),col="black",lwd=1.1)+theme(axis.text=element_text(size=30))+geom_ribbon(aes(ymin=ICcte05, ymax=ICcte95), alpha=0.1, fill = "blue",color = "blue")+geom_line(aes(y=ICcte05),col="blue")+geom_line(aes(y=ICcte95),col="blue")+geom_hline(yintercept=log(CTETheo),col="red")
ICcte

ggsave("CTEflood.pdf",ICcte) # save figure in pdf format

