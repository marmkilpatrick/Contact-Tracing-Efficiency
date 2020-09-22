setwd("") #setwd to location where Fig 4_5 ModSimReps.r lives

#common parameters--------------
shape=function(mu,sig) mu^2/sig^2;scale = function(mu,sig) sig^2/mu;mn=3.1;sig2=2.15
presympInf=pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2))/2.3
postsympInf=(pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/5
ssInf=(pgamma(13.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/6
d_msTest=5;f_asymp=0.2; f_hosp=0.036; d_EPS=3.2; d_EA=1/( f_asymp*(1/d_EPS)/(1-f_asymp) ); d_PSMS=2.3; d_MSSS=8;
beta   = .37; cTraceMx = 1; msTest  = 1/d_msTest; ssTest  = 1.; 
qEA  = 1/d_EA; qEPS  = 1/d_EPS; qPSMS  = 1/d_PSMS; qMSSS  = 1/d_MSSS; gammaA  = 1/8;
gammaMS=1/5; #gammaMS  = 1/d_MSSS*(1/f_hosp - 1)-1/d_EA;
gammaSS  = 1/10.7; gammaQ  = 1/14; sigmaA  = 1; sigmaPS  = presympInf/postsympInf; sigmaMS  = 1 ;sigmaSS  = ssInf/postsympInf
alpha  = 0.0066 * (gammaMS + gammaSS + 1/d_EA)/(1-0.0066) #based on IFR verity

Ninf=1 # # of infections for initial R0 w/ contact tracing
Ncont_p_inf=15 # contacts per infection 15
NCCPD=12 # # of contacts a contact tracer can contact per day to isolate 12
NCT=15 # 15# of contact tracers for population - based on data for SCC was 42; Roadmap 15/100,000: https://www.natlawreview.com/article/california-s-resilience-roadmap-and-guidance-to-employers-stage-two-reopening
#-----------------


#--- Fig 4---------------------
#Fig 4A,B
f_mstr=1 #optimistic fraction traced
kappa = 0.65 # moderate social distancing to match R0 w/ CT; kappa = .791 if f_mstr=0.5 to match w/ CT
FrCaseTr=f_mstr*msTest/(msTest+qMSSS+gammaMS)/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)

cTrace=0 #No contact tracing
R01noCT =((kappa * beta) / (cTrace  + qEA + qEPS)) * #initial R0 w/out CT
  (((qEA * sigmaA) / (cTrace + gammaA)) + 
     ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
     ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
     ((qEPS * qPSMS * qMSSS * sigmaSS) /
        ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
           (cTrace + ssTest + alpha + gammaSS))));R01noCT

#Initial conditions (Et was 1)
nt=225;TS = 1 #Time step = 1 / fraction of a day #time steps
St=100000;Et=5;IAt=0;IPSt=0;IMSt=0;ISSt=0;N_I_tr=0;Rt=0;Qt=0;In=Dt=0;Nt=St+Et+IAt+IPSt+IMSt+ISSt+Rt+Qt

R0t=R01noCT;#Fig 4A,B
NSIM=1;Et=1;CTR=F;R0t=R01noCT;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) 
source("Fig 4_5 ModSimReps.r") 
outF=data.frame(day=seq(0,nt,by=1/TS),St,Et,In,Rt)

options(scipen=5);titloc=175
par(mfcol=c(2,3),mar=c(0, 3, 0, 0),mgp=c(1.5,.5,0),oma = c(.5, .5, 0.5, 0.5),cex=1.25)
matplot(x=outF[,1],y=outF[,2:ncol(outF)],col=1:4,type="l",xaxt="n",lty=1,
        xlab="Day",ylab="Number of people")
legend(-10,.8e5,c("S","E","I","R"),col=1:4,lty=1,bty="n",seg.len=.75)
text(0,95000,"A")
text(titloc,95000,"Social distancing")
text(titloc,90000,"No contact tracing")
par(mar=c(3, 3, 0, 0) )#b,l,t,r
plot(R0t[-1]~outF[-1,1],xlab="Day",ylab=expression(R[t]),type="l",ylim=c(0,2.7)) #Rt over time
text(0,2.6,"B")

#Fig 4C,D - Unlimited CT, no social distancing
kappa=1 #no social distancing
NCT=1500 # Number of contact tracers: effectively Unlimited
FrCaseTr=f_mstr*msTest/(msTest+qMSSS+gammaMS)/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)
cTraceF1 = function(FrCaseTr,Delay,Ninf,Ncont_p_inf,NCT,NCCPD) { #function estimating cTrace parameter
  dTracing=Delay+0.5*Ninf*Ncont_p_inf/(NCT*NCCPD) #duration until traced (2=1 for test results 1 for symptoms)
  return(FrCaseTr*qEPS/(qEPS+qEA)*1/dTracing)} #if tracing all symptomatic cases
cTrace=cTraceF1(FrCaseTr,d_msTest,Ninf,Ncont_p_inf,NCT,NCCPD) #tracing removal rate w/ Ninf=1
R01 =((kappa * beta) / (cTrace  + qEA + qEPS)) * #initial R0
  (((qEA * sigmaA) / (cTrace + gammaA)) + 
     ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
     ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
     ((qEPS * qPSMS * qMSSS * sigmaSS) /
        ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
           (cTrace + ssTest + alpha + gammaSS))));R01

NSIM=1;Et=1;CTR=T;R0t=R01;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM)
source("Fig 4_5 ModSimReps.r") 
outF=data.frame(day=seq(0,nt,by=1/TS),St,Et,In,Rt)

par(mar=c(0, 0, 0, 0),mgp=c(1.5,.5,0),cex=1.25)
#par(mfrow=c(2,1),mar=c(0, 3, 0, 0), oma = c(3.5, .5, 0.5, 0.5),mgp=c(1.5,.5,0),cex=1.25)#b,l,t,r
matplot(x=outF[,1],y=outF[,2:ncol(outF)],col=1:4,type="l",xaxt="n",yaxt="n",lty=1,
        xlab="Day",ylab="Number of people")
#legend(50,.8e5,c("S","E","I","R"),col=1:4,lty=1,bty="n")

text(0,95000,"C")
text(titloc,95000,"No social distancing")
text(titloc,90000,"Unlimited contact tracing")
par(mar=c(3, 0, 0, 0) )#b,l,t,r
plot(R0t[-1]~outF[-1,1],xlab="Day",type="l",yaxt="n",ylim=c(0,2.7)) #Rt over time
#abline(h=R0t[2],lty=2)
text(0,2.6,"D")

#Fig 4E,F
kappa=1 #no social distancing
NCT=5 # Number of contact tracers - limited pre-scale up
FrCaseTr=f_mstr*msTest/(msTest+qMSSS+gammaMS)/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)
cTraceF1 = function(FrCaseTr,Delay,Ninf,Ncont_p_inf,NCT,NCCPD) { #function estimating cTrace parameter
  dTracing=Delay+0.5*Ninf*Ncont_p_inf/(NCT*NCCPD) #duration until traced (2=1 for test results 1 for symptoms)
  #  return(parsR0$qMSSS/(parsR0$qMSSS+parsR0$gammaMS+parsR0$qEA)*1/dTracing)} #if only severe cases traced
  return(FrCaseTr*qEPS/(qEPS+qEA)*1/dTracing)} #if tracing all symptomatic cases
cTrace=cTraceF1(FrCaseTr,d_msTest,Ninf,Ncont_p_inf,NCT,NCCPD) #tracing removal rate w/ I=1
R01 =((kappa * beta) / (cTrace  + qEA + qEPS)) * #initial R0
  (((qEA * sigmaA) / (cTrace + gammaA)) + 
     ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
     ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
     ((qEPS * qPSMS * qMSSS * sigmaSS) /
        ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
           (cTrace + ssTest + alpha + gammaSS))));R01

NSIM=1;Et=1;CTR=T;R0t=R01;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM)
source("Fig 4_5 ModSimReps.r") 
outF=data.frame(day=seq(0,nt,by=1/TS),St,Et,In,Rt)

par(mar=c(0, 0, 0, 0),mgp=c(1.5,.5,0),cex=1.25)
matplot(x=outF[,1],y=outF[,2:ncol(outF)],col=1:4,type="l",xaxt="n",yaxt="n",lty=1,
        xlab="Day",ylab="Number of people")
text(0,95000,"E")
text(titloc,95000,"No social distancing")
text(titloc,90000,"Limited contact tracing")
par(mar=c(3, 0, 0, 0) )#b,l,t,r
plot(R0t[-1]~outF[-1,1],xlab="Day",type="l",yaxt="n",ylim=c(0,2.7)) #Rt over time
#abline(h=R0t[2],lty=2)
text(0,2.6,"F")


#Fig 5------------------------------
kappa=0.7 #Limited social distancing
NCT=15 # of contact tracers for population - based on data for SCC was 42; Roadmap 15/100,000: https://www.natlawreview.com/article/california-s-resilience-roadmap-and-guidance-to-employers-stage-two-reopening
f_mstr=.5 #Fraction of mildly symptomatic cases traced 1 for Fig 4; 0.5 for fig 5
FrCaseTr=f_mstr*msTest/(msTest+qMSSS+gammaMS)/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)
cTraceF1 = function(FrCaseTr,Delay,Ninf,Ncont_p_inf,NCT,NCCPD) { #function estimating cTrace parameter
  dTracing=Delay+0.5*Ninf*Ncont_p_inf/(NCT*NCCPD) #duration until traced (2=1 for test results 1 for symptoms)
  #  return(parsR0$qMSSS/(parsR0$qMSSS+parsR0$gammaMS+parsR0$qEA)*1/dTracing)} #if only severe cases traced
  return(FrCaseTr*qEPS/(qEPS+qEA)*1/dTracing)} #if tracing all symptomatic cases
cTrace=cTraceF1(FrCaseTr,d_msTest,Ninf,Ncont_p_inf,NCT,NCCPD) #tracing removal rate w/ I=1
R01 =((kappa * beta) / (cTrace  + qEA + qEPS)) *
         (((qEA * sigmaA) / (cTrace + gammaA)) + 
            ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
            ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
            ((qEPS * qPSMS * qMSSS * sigmaSS) /
               ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
                  (cTrace + ssTest + alpha + gammaSS))));R01

cTrace=0 #zero tracing
R01noCT =((kappa * beta) / (cTrace  + qEA + qEPS)) *
  (((qEA * sigmaA) / (cTrace + gammaA)) + 
     ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
     ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
     ((qEPS * qPSMS * qMSSS * sigmaSS) /
        ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
           (cTrace + ssTest + alpha + gammaSS))));R01noCT


#Fig 5 Initial conditions
set.seed(1)
St=100000;Et=5;IAt=0;IPSt=0;IMSt=0;ISSt=0;N_I_tr=0;Rt=0;Qt=0;R0t=0;In=Dt=0;Nt=St+Et+IAt+IPSt+IMSt+ISSt+Rt+Qt
nt=365 #simulation length in days
TS = 1 #Inverse of time step in days; TS=1 => 24 hr time step; TS=6 => 4 hr time step

# Run simulations via Fig 4_5 ModSimReps.r and plot
par(mfrow=c(2,2),mar=c(0, 0, 0, 0), oma = c(3.5, 3.5, 0.5, 0.5),mgp=c(1.5,.5,0),cex=1.25)#b,l,t,r
ht=5000 #ylim max for all panels
outD=data.frame(day=seq(0,nt,by=1/TS)) #x-axis time label

#Fig 5A: Run NSIM=50 simulations w/ E=5, CT, stoch
NSIM=50;Et=5;CTR=T;R0t=R01;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",xaxt="n",ylim=c(0,ht))
text(nt*.9,ht*.95,expression('CT, E'[0]*'=5'));text(nt*.05,ht*.95,"A")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(nt*.9,ht*.85,paste("Frac. est. ",sum(mxI>3*Et[1])/NSIM))
#Fig 5A: black (deterministic) line
Et=5;CTR=T;R0t=R01;stoch=F;NSIM=1;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r") #Run w/ E=5, CT, NO stoch
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)

#Fig 5B: Run w/ E=50, CT, stoch
NSIM=50;Et=50;CTR=T;R0t=R01;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",xaxt="n",yaxt="n",ylab="",ylim=c(0,ht))
text(nt*.9,ht*.95,expression('CT, E'[0]*'=50'));text(nt*.05,ht*.95,"B")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(nt*.9,ht*.85,paste("Frac. est. ",sum(mxI>Et[1])/NSIM))
#Fig 5B: black (deterministic) line
NSIM=1;Et=50;CTR=T;R0t=R01;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r")
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)

#Fig 5c: Run w/ E=5, NO CT,  stoch
NSIM=50;Et=5;CTR=F;R0t=R01noCT;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r") 
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",ylab="",ylim=c(0,ht))
text(nt*.9,ht*.95,expression('No CT, E'[0]*'=5'));text(nt*.05,ht*.95,"C")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(nt*.9,ht*.85,paste("Frac. est. ",sum(mxI>Et[1])/NSIM))
#Fig 5C: black (deterministic) line
NSIM=1;Et=5;CTR=F;R0t=R01noCT;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r") #Run w/ E=5, NO CT, NO stoch
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)

#Fig 5d: Run w/ E=50, NO CT,  stoch
NSIM=50;Et=50;CTR=F;R0t=R01noCT;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4_5 ModSimReps.r")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",ylab="",yaxt="n",ylim=c(0,ht))
text(nt*.9,ht*.95,expression('No CT, E'[0]*'=50'));text(nt*.05,ht*.95,"D")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(nt*.9,ht*.85,paste("Frac. est. ",sum(mxI>Et[1])/NSIM))
#Fig 5D: black (deterministic) line
NSIM=1;Et=50;CTR=F;R0t=R01noCT;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
stoch=F;source("Fig 4_5 ModSimReps.r") #Run w/ E=5, NO CT, NO stoch
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)
mtext("Day", side = 1, outer = TRUE, line = 1.7,cex=1.7)
mtext("Number Latently Infected", side = 2, outer = TRUE, line = 1.7,cex=1.7)


