#Import necessary libraries

## @knitr library
library(R2WinBUGS)
library(R2jags)
library(ggmcmc)
library(matrixStats)
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(beepr)
library(knitr)

## @knitr cores
cores<-detectCores()

## @knitr wd
setwd("~/")

## @knitr source
sourceCpp("derivedvalues20.cpp")
sourceCpp("probest.cpp")
sourceCpp("probest_multinomial.cpp")

## @knitr dataload
h2d<-read.csv("htod.final.csv")  #susceptible to dead data
colnames(h2d)[1] <- "sheep"
h2i<-read.csv("htoi.final.csv")  #susceptible to infected data
i2dd<-read.csv("itodd.final.csv") #infected to disease death data
i2d<-read.csv("itod.final.csv")  #infected to non-disease death data
i2r<-read.csv("itor.final.csv")  #infected to recovered data
r2i<-read.csv("rtori.final.csv")  #recovered to infected data
r2d<-read.csv("rtod.final.csv")  #recovered to non-disease death data
ritod<-read.csv("ritod.final.csv") #reinfected to non-disease death data
ritodd<-read.csv("ritodd.final.csv") #reinfected to disease death data 
hcens<-read.csv("hcens.final.csv") #censored healthy (no clinical signs of disease)
icens<-read.csv("icens.final.csv") #censored infected
rcens<-read.csv("rcens.final.csv") #censored recovered

colnames(ritod) <- colnames(i2d)  #put reinfected with infecteds
i2d <- rbind(i2d,ritod)

colnames(ritodd) <- colnames(i2dd)
i2dd <- rbind(i2dd,ritodd)

## @knitr standardize
##Standardize covariates - for each state's data
covar.elisa <- read.csv("covariates.rawelisa.final.csv")

#Initial Elisa values
covar.elisa$std.initial <- (covar.elisa$initialELISA - 
                              mean(covar.elisa$initialELISA,na.rm=T))/
  sd(covar.elisa$initialELISA,na.rm=T)

covar.elisa$std.initial[is.na(covar.elisa$std.initial)] <- 0 #set missing = mean

#Elisa values during infection
covar.elisa$infectELISA <- as.numeric(as.character(covar.elisa$infectELISA))

covar.elisa$std.infectELISA <- (covar.elisa$infectELISA -
                                  mean(covar.elisa$infectELISA,na.rm=T))/
  sd(covar.elisa$infectELISA,na.rm=T)
covar.elisa$std.infectELISA [is.na(covar.elisa$std.infectELISA )] <- 0 #set missing = mean

#Elisa values during recovery
covar.elisa$std.recovELISA <- (covar.elisa$recovELISA - 
                                 mean(covar.elisa$recovELISA,na.rm=T))/
  sd(covar.elisa$recovELISA,na.rm=T)

covar.elisa$std.recovELISA [is.na(covar.elisa$std.recovELISA )] <- 0 #set missing = mean

colnames(covar.elisa)[colnames(covar.elisa) %in% c("std.initial","std.infectELISA")] <-
  c("elisaint","elisaa")
levels(covar.elisa$sheep)[2] <- levels(covar$sheep)[2]

###Function to add covariate values to each dataset
covmerge<-function(x,y){
  out<-merge(x,y, by="sheep")
  return(out)
}

h2d<-covmerge(h2d,covar.elisa) #susceptible to dead
h2i<-covmerge(h2i, covar.elisa) #susceptible to infected
i2dd<-covmerge(i2dd,covar.elisa) #infected to disease-dead
i2d<-covmerge(i2d,covar.elisa)   #infected to dead
i2r<-covmerge(i2r,covar.elisa)   #infected to recovered
r2i<-covmerge(r2i,covar.elisa)   #recovered to reinfected
r2d<-covmerge(r2d,covar.elisa)   #recovered to dead
hcens<-covmerge(hcens,covar.elisa) #susceptible - censored
icens<-covmerge(icens,covar.elisa) #infected - censored
rcens<-covmerge(rcens,covar.elisa) #recovered - censored


## @knitr dataorg
###Function organize data for creation of data file
extrdata<-function(x,namelist){
  nameout<-vector(length=11,mode="list")
  names(nameout)<-namelist
  nameout[[1]]<-nrow(x)
  nameout[[2]]<-length(unique(x$sheep))
  nameout[[3]]<- x[,grep("left",colnames(x))]
  nameout[[4]]<- x[,grep("right",colnames(x))]
  nameout[[5]]<- x[,grep("censor",colnames(x))]
  nameout[[6]]<- x[,grep("elisaint",colnames(x))]
  nameout[[7]]<- x[,grep("elisaa" ,colnames(x))]
  nameout[[8]]<- x[,grep("dist",colnames(x))]
  nameout[[9]]<- x[,grep("strain",colnames(x))]
  nameout[[10]]<- x[,grep("sex",colnames(x))]-1
  nameout[[11]]<- x[,grep("igs400",colnames(x))] 
  
  return(nameout)
}


#Susceptible to dead #1
h2da<-extrdata(h2d,c("records1","n1","left1","right1","censored1",
                     "elisaint","elisaa","dist","intstr","sex","igs400"))

#Sucsceptible to infected #2
h2ia<-extrdata(h2i,c("records2","n2","left2","right2","censored2",
                     "elisaint","elisaa","dist","intstr","sex","igs400"))

#Infected to infected disease dead #3
i2dda<-extrdata(i2dd,c("records3","n3","left3","right3","censored3",
                       "elisaint","elisaa","dist","intstr","sex","igs400"))

#Infected to infected recovered #4
i2ra<-extrdata(i2r,c("records4","n4","left4","right4","censored4",
                     "elisaint","elisaa","dist","intstr","sex","igs400"))

#Infected to infected non-disease dead #5
i2da<-extrdata(i2d,c("records5","n5","left5","right5","censored5",
                     "elisaint","elisaa","dist","intstr","sex","igs400"))

#Recovered to infected #6
r2ia<-extrdata(r2i,c("records6","n6","left6","right6","censored6",
                     "elisaint","elisaa","dist","intstr","sex","igs400"))

#Recovered to non-disease dead #7 
r2da<-extrdata(r2d,c("records7","n7","left7","right7","censored7",
                     "elisaint","elisaa","dist","intstr","sex","igs400"))

## @knitr datain1
T1<-800  #study length
T2<-666  #max time after infection for analysis
T3<-533  #max time after recovery for analysis

lefth<-hcens$lefth  #Censored individuals for susceptible/healthy state
lefti<-icens$lefti  #Censored individuals for infected state
leftr<-rcens$leftr  #Censored individuals for recovered state
righth<-hcens$righth
righti<-icens$righti
rightr<-rcens$rightr
nh1<-length(unique(hcens$sheep)) 
ns1<-length(unique(icens$sheep))
nr1<-length(unique(rcens$sheep))

## @knitr distmat
##Distance matrix for time-kernel convolution
intvl<-1 # day interval for knots
knots<-seq(1,800,intvl) #daily knots
distmat<-matrix(0,T1,length(knots))

for(i in 1:nrow(distmat)){
  for(j in 1:length(knots)){  
    distmat[i,j]<-abs(i-knots[j]) #absolute value
  }
}

## @knitr indices1
col.index<-matrix(0,nrow(distmat),2)  #create indices of knots to include in kernel effect  
##for each time point -- truncated at 120 days

col.index<-apply(distmat,1,function(x){  #for infection
  x[x>120]<--1
  out<-c(min(which(x>=0)),
         max(which(x>=0)))
  return(out)
})

col.index<-t(col.index)

col.index2<-matrix(0,nrow(distmat),2) #for disease death

col.index2<-apply(distmat[1:T2,1:length(seq(1,T2,intvl))],1,function(x){
  x[x>120]<--1
  out<-c(min(which(x>=0)),
         max(which(x>=0)))
  return(out)
})

col.index2<-t(col.index2)

## @knitr zeros
z<-0  #Zeros for trick in BUGS/JAGS

## @knitr dataorg2
###create infected matrix for dcat
ninfected2<-cbind(rep(1,i2dda$n3+i2ra$n4+i2da$n5),c(rep(1,i2dda$n3),rep(2,i2ra$n4),
                                                    rep(3,i2da$n5)))

###Create covariate matrix for ninfected
covarcomb<-cbind(c(i2dda$elisaint[i2dda$left3==1],i2ra$elisaint[i2ra$left4==1],
                   i2da$elisaint[i2da$left5==1]),
                 c(i2dda$elisaa[i2dda$left3==1],i2ra$elisaa[i2ra$left4==1],
                   i2da$elisaa[i2da$left5==1]),
                 c(i2dda$dist[i2dda$left3==1],i2ra$dist[i2ra$left4==1],
                   i2da$dist[i2da$left5==1]),
                 c(i2dda$intstr[i2dda$left3==1],i2ra$intstr[i2ra$left4==1],
                   i2da$intstr[i2da$left5==1]),
                 c(i2dda$sex[i2dda$left3==1],i2ra$sex[i2ra$left4==1],
                   i2da$sex[i2da$left5==1]),
                 c(i2dda$igs400[i2dda$left3==1],i2ra$igs400[i2ra$left4==1],
                   i2da$igs400[i2da$left5==1]))

#covariates for i2dda, i2ra, and i2da
colnames(covarcomb)<-c("elisaint","elisaa", "dist","intstr","sex","igs400") 

#change column names
colnames(rcens)[colnames(rcens)=="initial.strain"] <- "strain"

colnames(icens)[colnames(icens)=="initial.strain"] <- "strain"

colnames(hcens)[colnames(hcens)=="initial.strain"] <- "strain"

## @knitr Bugsdata
#R list for creating datafile for use with BUGS and used directly in JAGS
datain<-list(records1=h2da$records1,n1 = h2da$n1, left1 = h2da$left1, right1 =h2da$right1, 
             censored1 = h2da$censored1, elisaint1 = h2da$elisaint, intstr1 = h2da$intstr,
             
             records2=h2ia$records2,n2 = h2ia$n2, left2 = h2ia$left2, right2 =h2ia$right2, 
             censored2 = h2ia$censored2, elisaint2 = h2ia$elisaint,  
             dist2 = h2ia$dist, intstr2 = h2ia$intstr,
            
             records3=i2dda$records3, left3 = i2dda$left3, right3 = i2dda$right3, 
             censored3 = i2dda$censored3, elisaint3 = i2dda$elisaint,elisaa3=i2dda$elisaa, 
             
             records4=i2ra$records4,left4 = i2ra$left4, right4 = i2ra$right4, 
             censored4 = i2ra$censored4, elisaint4 = i2ra$elisaint, elisaa4 = i2ra$elisaa, 
            intstr4 = i2ra$intstr, 
             
             records5=i2da$records5,left5 = i2da$left5, right5 = i2da$right5, 
             censored5 = i2da$censored5, 
             
             records6=r2ia$records6,n6 = r2ia$n6, left6 = r2ia$left6, 
             right6 = r2ia$right6, censored6 = r2ia$censored6, 
             intstr6 = r2ia$intstr, 
             
             records7=r2da$records7,n7 = r2da$n7, left7 = r2da$left7,right7 = r2da$right7, 
             censored7 = r2da$censored7, intstr7 = r2da$intstr,
             
             elisab = hcens$elisaint, distb = hcens$dist, 
             intstrb = hcens$strain,
             
             elisac = icens$elisaint, elisaca=icens$elisaa, 
             intstrc = icens$strain, igs400c = icens$igs400,
             
             intstrr = rcens$strain, 
             
             covarcomb = covarcomb,T1=T1, T2=T2, lefth=lefth, lefti=lefti, 
             righti=righti, leftr=leftr, rightr=rightr, z=z, nh1=nh1, ns1=ns1, nr1=nr1, 
             distmat=distmat,
             
             NI =i2dda$n3+i2ra$n4+i2da$n5, ninfected=ninfected2, nconst=1/sqrt(2*pi),
             nknots=length(knots), nknots2=length(seq(1,666,intvl)),
             colindex=col.index, colindex2=col.index2)

## @knitr temp2

####################################################BEGIN OF MODEL STATEMENT
####Note:  Same as model statement for BUGS - but create an R function instead by 
####pasting BUGS code into R.  Must remove "model{" and the closing "}" for "model"
####from BUGS code.

## @knitr model1


model <- function(){
  
  ####Susceptible, infected, recovered model 
  
  ## 1 - Susceptible to dead transition
  ## 2 - Susceptible to infected transition (i.e., clinical disease)
  ## 3 - Infected to disease related dead transition
  ## 4 - Infected to recovered transition (i.e., no more/background clinical signs)
  ## 5 - Infected to non-disease related death transition
  ## 6 - Recovered to infected transition
  ## 7 - Recovered to dead transition

## @knitr itime   
  ##Kernel convolution- time to infection
  # lnsdk ~ dunif(0,4.1)
  # sdk <- exp(lnsdk)
  tauk ~ dgamma(1,1)#<-1/(sdk*sdk)
  sdk <- pow(1/tauk,0.5)
  stauk<-sqrt(tauk)
  
  taua ~ dgamma(1,1) #precision of time random effect
  sda <- 1/sqrt(taua)
  
  
  for(i in 1:nknots){
    alpha[i] ~ dnorm(0,1)
    alphau[i] <- sda*alpha[i]
  }
  
  ratioinf<-sdk/sda #ratio of variability

## @knitr itimeeffects  
  for(i in 1:T1){
    for(j in 1:nknots){ 
      temp1[i,j] <-(stauk*nconst*exp(-0.5*distmat[i,j]^2*tauk))
    }
  }

  for(i in 1:T1){
    for(j in colindex[i,1]:colindex[i,2]){
      temp[i,j] <- temp1[i,j]*alphau[j]
      
    }
    KA[i] <- sum(temp[i,colindex[i,1]:colindex[i,2]])
  }
  
  
## @knitr ddtime
  
  ##Kernel convolution- time to infection to disease-death
  taukdd ~ dgamma(1,1)
  sdkdd <- pow(1/taukdd,0.5)

  tauadd ~ dgamma(1,1) #precision of time random effect
  sdadd <-  1/sqrt(tauadd)

  staukdd<-sqrt(taukdd)
  ratiodd<-sdkdd/sdadd #ratio of variability
  
  for(i in 1:nknots2){
    alphadd[i] ~ dnorm(0,1)
    alphaddu[i]<-sdadd*alphadd[i]
  }
  
## @knitr ddtimeeffects  
  for(i in 1:T2){
    for(j in 1:nknots2){ 
      tempdd1[i,j] <-(staukdd*nconst*exp(-0.5*distmat[i,j]^2*taukdd))
    }
  }
 
  
  for(i in 1:T2){
    for(j in colindex2[i,1]:colindex2[i,2]){
      tempdd[i,j] <- tempdd1[i,j]*alphaddu[j]# (tempdd1[i,j]/tempdd2[j])*alphaddu[j]
      
    }
    KAdd[i] <- sum(tempdd[i,colindex2[i,1]:colindex2[i,2]])
  }
  
  
## @knitr deadhazard
  ##1
    #Priors
    gamma1 ~ dunif(-100,100) #baseline-log hazard
    
    #Event-time Log-Likelihood
    for(j in 1:records1){
      for(k in left1[j]:(right1[j]-1)){
        UCH1[j,k] <- exp(gamma1)  # constant hazard model for non-disease related deaths
      }
      surv1[j] <- (-sum(UCH1[j,left1[j]:(right1[j]-1)]))*censored1[j] + 
        (1-censored1[j])*log(1-exp(-sum(UCH1[j,left1[j]:(right1[j]-1)])) )
    }

## @knitr probdead    
    #Priors
    for(i in 1:2){
      beta01[i] ~ dnorm(0,0.01)
    }
    
    beta01str[1] <- 0 #set to zero for base-line strain
    for(i in 2:4){
      beta01str[i] ~ dnorm(0,0.01) # strain parm- categorical	
    }
    
    
    #Risk-type Log-Likelihood
    for(i in 1:n1){
      logit(p1[i])<-beta01[1] + beta01[2]*elisaint1[i] + beta01str[intstr1[i]]
                                                             
                                                  
      risk1[i] <- log(1 - p1[i]) #probability not transitioning to infected
    }
    
    
    logl1<-sum(surv1[]) + sum(risk1[])  #log likelihood for all individuals in ##1
    
## @knitr infechazard    
    ##2
    #Priors
    gamma2 ~ dunif(-100,100) #baseline-log hazard
    
    for(i in 1:2){
      betah1[i] ~ dnorm(0,0.01)
    }
    
    
    #Event-time Log-Likelihood
    for(j in 1:records2){
      for(k in left2[j]:(right2[j]-1)){
        UCH2[j,k] <- exp(gamma2 + betah1[1] *elisaint2[j]+betah1[2]*dist2[j] + KA[k])
      }
        surv2[j] <- (-sum(UCH2[j,left2[j]:(right2[j]-1)]))*censored2[j] + 
        (1-censored2[j])*log(1-exp(-sum(UCH2[j,left2[j]:(right2[j]-1)])) )
    }

## @knitr probinfec    
    #Risk-type Log-Likelihood
    for(i in 1:n2){
      logit(p2[i])<-beta01[1]+beta01[2]*elisaint2[i]+beta01str[intstr2[i]]
      risk2[i]<-log(p2[i])  #probability transitioning to infected
    }
    
    logl2<-sum(risk2[])+sum(surv2[])  #log likelihood for all individuals in ##2

## @knitr scens        
    ##Stay Susceptible - Censored Individuals
    for(i in 1:nh1){
      for(k in lefth[i]:(T1-1)){
        UCHa[i,k] <- exp(gamma1)   # constant hazard model for non-disease related deaths
      }
      for(k in lefth[i]:(T1-1)){
        UCHb[i,k] <- exp(gamma2 + betah1[1]*elisab[i] + betah1[2]*distb[i] + KA[k]) 
        # constant hazard model for infection
      }
      
      # transitioning to infected
      logit(p1a[i]) <-beta01[1] + beta01[2]*elisab[i] + beta01str[intstrb[i]] 
     
      #log likelihood for censored 
      loglh[i] <- log(p1a[i] *exp(-sum(UCHb[i,lefth[i]:(T1-1)]))+(1-p1a[i])*
                        exp(-sum(UCHa[i,lefth[i]:(T1-1)])))  
    }

## @knitr slike    
    logl1all<-sum(loglh[])+logl1+logl2  #Susceptible contribution to log-likelihood	
    
## @knitr infhazpriors
   
    ##Priors
    gamma3~ dunif(-100,100) #baseline-log hazard
    gamma4 ~ dunif(-100,100) #baseline-log hazard
    
    for(i in 1:2){
      betah3[i] ~ dnorm(0,0.01)
      betah4[i] ~ dnorm(0,0.01)
    }
    
    
    betah4str[1] <-0
    betah4str[2] <-0
    betah4str[3] ~ dnorm(0,0.01) #only 2 and 3 strains recovered
    betah4str[4] <-0

## @knitr ddhaz       
    ##3
    
    #Event-time Log-Likelihood
    for(j in 1:records3){
      for(k in left3[j]:(right3[j]-1)){
        UCH3[j,k] <- exp(gamma3 + betah3[1]*elisaint3[j] + betah3[2]*elisaa3[j] +
                           KAdd[k]) # hazard model for disease-death
      }
      surv3[j] <- (-sum(UCH3[j,left3[j]:(right3[j]-1)]))*censored3[j] + 
        (1-censored3[j])*log(1-exp(-sum(UCH3[j,left3[j]:(right3[j]-1)])) )
    }
    
    logl3<-sum(surv3[])  #log likelihood for all individuals in ##3

## @knitr rechaz    
    
    ##4
    #Event-time Log-Likelihood
    for(j in 1:records4){
      for(k in left4[j]:(right4[j]-1)){
        UCH4[j,k] <- exp(gamma4 + betah4[1]*elisaint4[j] + betah4[2]*elisaa4[j] + 
                           betah4str[intstr4[j]]) # hazard model for recovery
      }
      surv4[j] <- (-sum(UCH4[j,left4[j]:(right4[j]-1)]))*censored4[j] + 
        (1-censored4[j])*log(1-exp(-sum(UCH4[j,left4[j]:(right4[j]-1)])) )
    }
    
    
    logl4<-sum(surv4[])  #log likelihood for all individuals in ##4
    
## @knitr idhaz     
    ##5 - constrained to be equal to ##1 - healthy deaths
    
    #Event-time Log-Likelihood
    for(j in 1:records5){
      for(k in left5[j]:(right5[j]-1)){
        UCH5[j,k] <- exp(gamma1)# constant hazard model for death - same as healthy death
      }
      surv5[j] <- (-sum(UCH5[j,left5[j]:(right5[j]-1)]))*censored5[j] + 
        (1-censored5[j])*log(1-exp(-sum(UCH5[j,left5[j]:(right5[j]-1)])) )
    }
    
    
    logl5<-sum(surv5[])  #log likelihood for all individuals in ##5
    
    
      
## @knitr infprobpriors 
    #Risk-type Likelihood
    ##Priors
    for(i in 1:4){
      beta0dd[i] ~ dnorm(0,0.01)
      beta0rec[i] ~ dnorm(0,0.01)
    }
    
    
    ## @knitr ddhaz 
    
    for(i in 2:4){
      beta0ddstr[i] ~ dnorm(0,0.01)
      beta0recstr[i] ~ dnorm(0,0.01)
    }	
    
    beta0ddstr[1]<-0  #set to zero - 400 set as base-line strain - probability
    beta0recstr[1]<-0
    
   
## @knitr infprobs        
    for(i in 1:NI){
      #probability of disease death
      P[i,1] <- exp(beta0dd[1] + beta0dd[2]*covarcomb[i,1] + beta0dd[3]*covarcomb[i,2] + 
                    beta0dd[4]*covarcomb[i,6] + beta0ddstr[covarcomb[i,4]])/
        (1 + exp(beta0dd[1] + beta0dd[2]*covarcomb[i,1] + beta0dd[3]*covarcomb[i,2] + 
                   beta0dd[4]*covarcomb[i,6] + beta0ddstr[covarcomb[i,4]]) + 
             exp(beta0rec[1] + beta0rec[2]*covarcomb[i,1] + beta0rec[3]*covarcomb[i,2]+
                  beta0rec[4]*covarcomb[i,6] + beta0recstr[covarcomb[i,4]])) 
      
      #probability of recovery
      P[i,2] <- exp(beta0rec[1] + beta0rec[2]*covarcomb[i,1] + 
                    beta0rec[3]*covarcomb[i,2] + beta0rec[4]*covarcomb[i,6] + 
                    beta0recstr[covarcomb[i,4]])/
        (1+exp(beta0dd[1] + beta0dd[2]*covarcomb[i,1] + beta0dd[3]*covarcomb[i,2] + 
              beta0dd[4]*covarcomb[i,6] + beta0ddstr[covarcomb[i,4]]) + 
           exp(beta0rec[1] + beta0rec[2]*covarcomb[i,1] + beta0rec[3]*covarcomb[i,2] +
              beta0rec[4]*covarcomb[i,6] + beta0recstr[covarcomb[i,4]]))  
      
      P[i,3]<-1-sum(P[i,1:2])   #probability of non-disease death
      
      ninfected[i,2] ~ dcat(P[i,])  #likelihood contribution
    }
    
## @knitr icens     
    ##Stay Infected - Censored Individuals
    for(i in 1:ns1){
      for(k in lefti[i]:(righti[i]-1)){
        UCHc[i,k] <- exp(gamma3 + betah3[1]*elisac[i] + betah3[2]*elisaca[i]+
                          KAdd[k])
            # hazard model for disease related deaths
      }
      for(k in lefti[i]:(righti[i]-1)){
        UCHd[i,k] <- exp(gamma4 + betah4[1]*elisac[i] + betah4[2]*elisaca[i] +
                           betah4str[intstrc[i]])  # hazard model for recovery
      }
      for(k in lefti[i]:(righti[i]-1)){
        UCHe[i,k] <- exp(gamma1)   # constant hazard model for non-disease related deaths
      }
      
      # probability of disease death
      pc[i,1] <-exp(beta0dd[1] + beta0dd[2]*elisac[i] + beta0dd[3]*elisaca[i] +
                   beta0dd[4]*igs400c[i]+ beta0ddstr[intstrc[i]])/
        (1 + exp(beta0dd[1] + beta0dd[2]*elisac[i] + beta0dd[3]*elisaca[i] +
                   beta0dd[4]*igs400c[i] + beta0ddstr[intstrc[i]]) + 
              exp(beta0rec[1]+beta0rec[2]*elisac[i] + beta0rec[3]*elisaca[i] + 
                    beta0rec[4]*igs400c[i] + beta0recstr[intstrc[i]]))     
                
      #probability of infected to recovered censored
      pc[i,2] <- exp(beta0rec[1] + beta0rec[2]*elisac[i] + beta0rec[3]*elisaca[i] +
                      beta0rec[4]*igs400c[i] + beta0recstr[intstrc[i]])/
        (1+exp(beta0dd[1] + beta0dd[2]*elisac[i] + beta0dd[3]*elisaca[i] + 
                 beta0dd[4]*igs400c[i] + beta0ddstr[intstrc[i]]) + 
           exp(beta0rec[1] + beta0rec[2]*elisac[i] + beta0rec[3]*elisaca[i] +
                 beta0rec[4]*igs400c[i] + beta0recstr[intstrc[i]]))  
               
 
      pc[i,3]<-1-sum(pc[i,1:2])
      
      #log likelihood for censored 
      logls[i] <- log(pc[i,1] *exp(-sum(UCHd[i,lefti[i]:(righti[i]-1)]))+pc[i,2]*
                        exp(-sum(UCHc[i,lefti[i]:(righti[i]-1)]))
                      +pc[i,3]*exp(-sum(UCHe[i,lefti[i]:(righti[i]-1)]))) 
    }

## @knitr ilike
    #Infected contribution to log-likelihood	
    logl2all<-sum(logls[]) + logl3 + logl4 + logl5  
    

## @knitr recpriors    
   ##6 
    #Priors
    gamma6 ~ dunif(-100,100) #baseline-log hazard
    betah6 ~ dnorm(0,0.01)
    beta6igs[1] <- 0
    beta6igs[2] ~ dnorm(0,0.01)

## @knitr rec2ihaz  
    #Event-time Log-Likelihood
    for(j in 1:records6){
      for(k in left6[j]:(right6[j]-1)){
        UCH6[j,k] <- exp(gamma6) # hazard model for recovered to infected
      }
      surv6[j] <- (-sum(UCH6[j,left6[j]:(right6[j]-1)]))*censored6[j] + 
        (1-censored6[j])*log(1-exp(-sum(UCH6[j,left6[j]:(right6[j]-1)])) )
    }
  
## @knitr recprobpriors     
    beta06[1] ~ dnorm(0,0.01)
    beta06str[1] <- 0 #set to zero for base-line strain
    beta06str[2] <- 0 #only 2 and 3 strains recover to infected
    beta06str[3] ~ dnorm(0,0.01) # strain parm- categorical	
    beta06str[4] <- 0

## @knitr reinfprob
   for(i in 1:n6){
      logit(p6[i])<-beta06[1] + beta06str[intstr6[i]]
      risk6[i] <- log(p6[i]) #probability recovered transitioning to infected
    }
    
 
    logl6<-sum(surv6[])+sum(risk6[])  #log likelihood for all individuals in ##6

## @knitr r2d    
    ##7
    #Event-time Log-Likelihood
    for(j in 1:records7){
      for(k in left7[j]:(right7[j]-1)){
        UCH7[j,k] <- exp(gamma1) # hazard model for recovered to infected
      }
      surv7[j] <- (-sum(UCH7[j,left7[j]:(right7[j]-1)]))*censored7[j] + 
        (1-censored7[j])*log(1-exp(-sum(UCH7[j,left7[j]:(right7[j]-1)])) )
    }
    
## @knitr rdprob    
    for(i in 1:n7){
      logit(p7[i])<-beta06[1] +  beta06str[intstr7[i]]
      risk7[i] <- log(1-p7[i]) #probability recovered to non-disease death
    }

    logl7<-sum(risk7[]) + sum(surv7[])  #log likelihood for all individuals in ##7
    
## @knitr rcens    
    ##Stay Recovered - Censored Individuals
    for(i in 1:nr1){
      for(k in leftr[i]:(rightr[i]-1)){
        UCHf[i,k] <- exp(gamma1)   # constant hazard model for non-disease related deaths
      }
      for(k in leftr[i]:(rightr[i]-1)){
        UCHg[i,k] <- exp(gamma6) # constant hazard model for infection
      }  
      
      #probability transitioning to infected
      logit(p7a[i]) <-beta06[1] + beta06str[intstrr[i]] 
      
      #log likelihood for censored 
      loglh7[i] <- log(p7a[i] *exp(-sum(UCHg[i,leftr[i]:(rightr[i]-1)])) + (1-p7a[i])*
                         exp(-sum(UCHf[i,leftr[i]:(rightr[i]-1)])))  
    }

## @knitr rlike    
    logl3all<-sum(loglh7[]) + logl6 + logl7  #Recovered contribution to log-likelihood	
    
## @knitr zerostrick
    #"Zeros" trick (Lunn et al. 2013, p.204)
    const<-10000 #ensures phi[i] >0
    phi <--(logl1all + logl2all +logl3all ) + const
      z ~ dpois(phi)
}

## @knitr temp2
####################################################################################END OF MODEL STATEMENT

## @knitr inits
###Create initial values - use random number generator so don't have to 
####specify starting values for each chain

init1 <- function(){list("gamma1"=runif(1,-3,0), "gamma2"=runif(1,-3,0), 
              "gamma3"=runif(1,-3,0), "gamma4"=runif(1,-3,0), "gamma6"=runif(1,-3,0), 
              "beta01"=runif(2,-2,0.5), "betah1"=runif(2,-2,0.5),"betah3"=runif(2,-2,0.5),
              "betah4"=runif(2,-2,0.5), "beta0dd"=c(-1,0,0,0),
              "beta0rec"=c(-1,0,0,0), "beta01str" = c(NA,runif(3,-4,0.5)), 
              "betah4str" = c(NA,NA,runif(1,-4,0),NA), "beta0ddstr"=c(NA,runif(3,-4,0.5)),
              "beta0recstr"=c(NA,runif(3,-4,0.5)),
              "beta06str"= c(NA,NA,runif(1,-4,0.5),NA),
              "beta06"=runif(1,-2,0.5),"alpha"=rep(0, length(knots)),"alphadd"=rep(0,length(seq(1,666,intvl))))}

## @knitr parms
#identify parms to monitor
parms <- c("beta0rec","beta0dd","beta0ddstr","beta0recstr","gamma1","gamma2","gamma3",
         "gamma4","gamma6","beta01","betah1","betah3","betah4","beta01str","betah4str",
         "beta06str","sdk", "sda","sdkdd", "sdadd","KA", 
         "KAdd","beta06","ratioinf","ratiodd") 
  
  ## @knitr fit
  jagsfit <- jags(datain, init1, parameters.to.save=parms, 
                jags.module = c("dic"), n.chains=3, n.iter=10, 
                model.file=model, n.burnin=1, n.thin=1, 
                DIC=FALSE, working.directory=NULL)  #fit the model in JAGS
  
  ## @knitr temp3
  for(i in 1:100){
    beep(2)
  }

## @knitr mcmcobj
####Summarize Posterior Distributions
jags.mcmc <- as.mcmc(jagsfit)  #convert to mcmc object for plotting below

## @knitr temp4
saveRDS(jags.mcmc,file="results_final19.rds")
jags.mcmc <- readRDS("results_final19.rds")


plot(jags.mcmc[,"sdk"],main="sdk")
plot(jags.mcmc[,"sdkdd"],main="sdkdd")
plot(jags.mcmc[,"beta06str[3]"],main="beta06str[3]")
plot(jags.mcmc[,"beta06"],main="beta06")
plot(jags.mcmc[,"KA[1]"],main="KA[1]")
plot(jags.mcmc[,"gamma3"],main="gamma3")
plot(jags.mcmc[,"beta0dd[1]"],main="beta0dd[1]")


## @knitr timeorder
##Obtain colnames for proper sorting of random time effects
orderfunc <- function(name){
  temp<-strsplit(colnames(jags.mcmc[[1]])[grep(name,colnames(jags.mcmc[[1]]))],"\\[")
  temp2<-vector(length=length(temp),mode="numeric")
  for(i in 1:length(temp)){
    temp2[i]<-as.numeric(strsplit(temp[[i]][2],"\\]")[[1]]  ) #colnames as numeric
  }
  return(order(temp2))
}


orderkadd <- orderfunc("KAdd")
orderka <- orderfunc("KA\\[")
test <- as.matrix(jags.mcmc)
test2 <- test[,grep("KA\\[",colnames(test))]
test2 <- test2[,orderka]


test3 <- test[,grep("KAdd",colnames(test))]
test3 <- test3[,orderkadd]


strain.design <- matrix(0,nrow(covar),4)
for(i in 1:nrow(covar)){
  strain.design[i,covar[i,"strain"]]<-1
}

## @knitr sumresults1
results <- deriveval(test[,"gamma1"],test[,"gamma2"],test[,"gamma3"],test[,"gamma4"],
          test[,"gamma6"], test[,grep("beta01",colnames(test))],
          test[,grep("beta06",colnames(test))],test[,grep("beta0dd",colnames(test))],
          test[,grep("beta0rec",colnames(test))],test[,grep("betah1",colnames(test))],
          test[,c(grep("betah3",colnames(test)),grep("betaigs3",colnames(test)))],
          test[,grep("betah4",colnames(test))],
          test2,
          test3,150,seq(0,720,by=30),
          cbind(1,covar.elisa[,"elisaint"], strain.design, covar.elisa[,"elisaa"], 
                covar.elisa[,"dist"], 0, covar.elisa[,"igs400"]), 
          cores)


colnames(results$Nmu) <- seq(0,720,by=30)
rownames(results$Nmu) <- c("Susceptible","Dead","Infected","Recovered","Disease-dead")
rownames(results$quantile025) <- c("Susceptible","Dead",
                                   "Infected","Recovered","Disease-dead")
rownames(results$quantile975) <- c("Susceptible","Dead",
                                   "Infected","Recovered","Disease-dead")

## @knitr temp6
saveRDS(results,file="popsimresults20.rds")

## @knitr summary
#remove random time effects for summary
jags.mcmc2 <- jags.mcmc[,-grep("KA\\[",colnames(jags.mcmc[[1]]))] 
jags.mcmc2 <- jags.mcmc2[,-grep("KAdd",colnames(jags.mcmc2[[1]]))]

for(i in 1:length(jags.mcmc2)){
  colnames(jags.mcmc2[[i]])[c(1:48)]<-c(
            "probinfec-int","probinfec-elisai","probinfec-400","probinfec-404",
            "probinfec-393", "probinfec-398", 
            
            "probr2i-int","probr2i-400", "probr2i-404", "probr2i-393", 
            "probr2i-398", 
            
            "probdd-int", "probdd-elisai", "probdd-elisaa", "probdd-igs400",
            "probdd-400", "probdd-404", "probdd-393", "probdd-398",
            
            "probrec-int","probrec-elisai", "probrec-elisaa", "probrec-igs400",
            "probrec-400", "probrec-404", "probrec-393", "probrec-398",
            
            "timeinfec-elisai", "timeinfec-dist", 
            
            "timei2dd-elisai", "timei2ddelisaa", 
            
            "timerecov-elisai", "timerecov-elisaa", 
            "timerecov-400", "timerecov-404", "timerecov-393", "timerecov-398", 
            
            "time2d-int", "time2i-int", "timei2dd-int", "timei2r-int", 
            "timer2i-int", "ratiokernelsd-dd","ratiokernelsd-inf","precision-inf",
            "precision-dd", "precision-time-inf","precision-time-dd")
  
}

sumresults <- summary(jags.mcmc2,quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))
write.csv(cbind(sumresults[[1]],sumresults[[2]]),"results_final19.csv")

## @knitr createplots
##Create plots of hazard/survival functions

#hazard of infected to disease death
data.p1 <- as.data.frame(results$hazdd)
colnames(data.p1) <- c("median","LCB","UCB")
hazdd.plot <- ggplot(data=data.p1, aes(x=1:nrow(data.p1),y=median,group=1)) + 
  geom_line() +
  geom_ribbon(aes(ymin=LCB,ymax=UCB),alpha=0.3) + 
  xlab("days") + 
  ylab("death due to disease hazard")
hazdd.plot

#hazard of susceptible to infected
data.p2 <- as.data.frame(results$hazinf)
colnames(data.p2) <- c("median","LCB","UCB")
hazinf.plot <- ggplot(data=data.p2,aes(x=1:nrow(data.p2),y=median,group=1)) +
  geom_line() + 
  geom_ribbon(aes(ymin=LCB,ymax=UCB),alpha=0.3) + 
  xlab("days") + 
  ylab("infection hazard")
hazinf.plot


#survival curve of susceptible
data.p3 <- as.data.frame(results$Ssucept)
colnames(data.p3) <- c("median","LCB","UCB")
Ss.plot <- ggplot(data=data.p3,aes(x=1:nrow(data.p3),y=median)) +
  geom_line() +  
  geom_ribbon(aes(ymin=LCB,ymax=UCB),alpha=0.3) +
  xlab("days") +
  ylab("Transition curve - susceptible")
Ss.plot


#survival curve of infected
data.p4 <- as.data.frame(results$Sinfect)
colnames(data.p4) <- c("median","LCB","UCB")
Si.plot <- ggplot(data=data.p4, aes(x=1:nrow(data.p4),y=median)) +
  geom_line() +  
  geom_ribbon(aes(ymin=LCB,ymax=UCB),alpha=0.3) + 
  xlab("days") +
  ylab("Transition curve - infected")
Si.plot

#survival curve of recovered
data.p5 <- as.data.frame(results$Srec)
colnames(data.p5) <- c("median","LCB","UCB")
Sr.plot <- ggplot(data=data.p5,aes(x=1:nrow(data.p5),y=median)) +
  geom_line() + 
  geom_ribbon(aes(ymin=LCB,ymax=UCB),alpha=0.3) + 
  xlab("days") +
  ylab("Transition curve - recovered")
Sr.plot

#plot of elisa effect for prob h2i

beta.in <- test[,grep("beta01\\[",colnames(test))]
covar.range <- round(range(c(datain$elisaint1,datain$elisaint2)),2)
X.in <- cbind(1,seq(covar.range[1],covar.range[2], by = 0.05))

prob.effect1 <- probest(X.in, beta.in, cores)

data.p6 <- data.frame(x = X.in[,2], mu = t(prob.effect1$mu), 
                      lcl = prob.effect1$quantile025,
                      ucl = prob.effect1$quantile975)

peffect1.plot <- ggplot(data=data.p6,aes(x=x,y=mu)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl,ymax=ucl),alpha=0.3) + 
  xlab("Initial standardized ELISA value") +
  ylab("Probability of infection")
peffect1.plot


#plot of average elisa effect for prob i2dd
beta.in1 <- test[,grep("beta0dd\\[",colnames(test))]
beta.in2 <- test[,grep("beta0rec\\[",colnames(test))]

beta.in1 <- beta.in1[,1:3] #remove strain covariate
beta.in2 <- beta.in2[,1:3] #remove strain covariate

covar.range1 <- round(range(c(datain$covarcomb[,2],datain$covarcomb[,2])),2)

X.in1 <- cbind(1, mean(datain$covarcomb[,1]), 
               seq(covar.range1[1],covar.range1[2], by = 0.05))

#Infected to disease dead
prob.effect2 <- probestmulti(X.in1, beta.in1, beta.in2, cores)

data.p7 <- data.frame(x = X.in1[,3], mu = t(prob.effect2$mu), 
                      lcl = prob.effect2$quantile025,
                      ucl = prob.effect2$quantile975)

peffect2.plot <- ggplot(data=data.p7,aes(x=x,y=mu)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl,ymax=ucl),alpha=0.3) + 
  xlab("Average standardized ELISA value") +
  ylab("Probability of dying of disease")
peffect2.plot


#Infected to recovered
prob.effect3 <- probestmulti(X.in1, beta.in2, beta.in1, cores)


data.p8 <- data.frame(x = X.in1[,3], mu = t(prob.effect3$mu), 
                      lcl = prob.effect3$quantile025,
                      ucl = prob.effect3$quantile975)

peffect3.plot <- ggplot(data=data.p8,aes(x=x,y=mu)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl,ymax=ucl),alpha=0.3) + 
  xlab("Average standardized ELISA value") +
  ylab("Probability of recovering from disease")
peffect3.plot



