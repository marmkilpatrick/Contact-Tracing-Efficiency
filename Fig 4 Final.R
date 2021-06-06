setwd("c:/marm/research/covid-19/modeling T-CT-I/") #setwd so source("Fig 4 ModSimRepsFinal.r") works
#params
#Infectiousness vs day from He et al 1c gamma distribution parameterized w/ mean, variance
shape=function(mu,sig) mu^2/sig^2;scale = function(mu,sig) sig^2/mu;mn=3.1;sig2=2.15
presympInf=pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2))/2.3 #relative daily infectiousness pre-symptomatic
postsympInf=(pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/5 #relative daily infectiousness mildly-symptomatic
ssInf=(pgamma(13.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/6 #relative daily infectiousness severely-symptomatic
pars = data.frame(N=100001,
                  alpha=0.2*(0.0066/0.8)/(1-(0.0066/0.8)),
                  beta = 0.35,
                  kappa=.6,
                  sigma_Ips=presympInf/postsympInf,
                  sigma_Is=1,
                  gamma_Is=1/5,
                  q_E_Ips=1/3.2,
                  q_Ips_Is=1/2.3,
                  q_Q_R=1/14,
                  tau_Is=1/10,
                  tracetime=0.5,
                  CTR=T,
                  Ncont_p_inf= 10, # contacts per infection
                  NCCTD=12, # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
                  NCT=15 # # of contact tracers for population
                  
)


#Fig 4---------------------------
set.seed(1);NSIM=50;nt=365;stoch=T #time steps
TS = 2 #Time step = 1 / fraction of a day
S=100000;E=5;Ips=0;Is=0;Iss=0;R=0;Q=0; In=0; N=S+E+Ips+Is+R+Q;R0t=NA;ft=NA
Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM)
outD=data.frame(day=seq(0,nt,by=1/TS)) #x-axis time label

tiff("Fig 4.tif",width=12,height=8,units="in",res=300,compression = "lzw")
#pdf("Fig 4.pdf",width=12,height=8)
par(mfrow=c(2,2),mar=c(0, 0, 0, 0), oma = c(3.5, 3.5, 0.5, 0.5),mgp=c(1.5,.5,0),cex=1.25)#b,l,t,r
ht=4700 #ylim max for all panels
xl=0.75*nt
#Fig 4A: Run NSIM=50 simulations w/ E=5, CT, stoch
S=100000;E=5;Ips=0;Is=0;Iss=0;R=0;Q=0; In=0; N=S+E+Ips+Is+R+Q;R0t=NA;ft=NA
NSIM=50;E=5;CTR=T;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",xaxt="n",ylim=c(0,ht))
text(xl,ht*.95,expression('Contact tracing, E'[0]*'=5'));text(nt*.05,ht*.95,"A")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(xl,ht*.85,paste("Fraction established ",sum(mxI>3*E[1])/NSIM))
#Fig 4A: black (deterministic) line
NSIM=1;Ea=1;Es=4;CTR=T;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R") #Run w/ E=5, CT, NO stoch
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)
R0t[1]
#Fig 4B: Run w/ E=50, CT, stoch
S=100000;E=5;Ips=0;Is=0;Iss=0;R=0;Q=0; In=0; N=S+E+Ips+Is+R+Q;R0t=NA;ft=NA
NSIM=50;E=50;CTR=T;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",xaxt="n",yaxt="n",ylab="",ylim=c(0,ht))
text(xl,ht*.95,expression('Contact tracing, E'[0]*'=50'));text(nt*.05,ht*.95,"B")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(xl,ht*.85,paste("Fraction established ",sum(mxI>E[1])/NSIM))
#Fig 4B: black (deterministic) line
NSIM=1;Ea=10;Es=40;CTR=T;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)

#Fig 4c: Run w/ E=5, NO CT,  stoch
S=100000;E=5;Ips=0;Is=0;Iss=0;R=0;Q=0; In=0; N=S+E+Ips+Is+R+Q;R0t=NA;ft=NA
NSIM=50;E=5;CTR=F;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",ylab="",ylim=c(0,ht))
text(xl,ht*.95,expression('No contact tracing, E'[0]*'=5'));text(nt*.05,ht*.95,"C")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(xl,ht*.85,paste("Fraction established ",sum(mxI>E[1])/NSIM))
#Fig 4C: black (deterministic) line
NSIM=1;Ea=1;Es=4;CTR=F;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)

#Fig 4d: Run w/ E=50, NO CT,  stoch
S=100000;E=5;Ips=0;Is=0;Iss=0;R=0;Q=0; In=0; N=S+E+Ips+Is+R+Q;R0t=NA;ft=NA
NSIM=50;E=50;CTR=F;stoch=T;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
matplot(x=outD[,1],y=Ens,type="l",col="gray",xlab="",ylab="",yaxt="n",ylim=c(0,ht))
text(xl,ht*.95,expression('No contact tracing, E'[0]*'=50'));text(nt*.05,ht*.95,"D")
mxI=apply(Ens, 2, function(x) max(x, na.rm = TRUE))
text(xl,ht*.85,paste("Fraction established ",sum(mxI>E[1])/NSIM))
#Fig 4D: black (deterministic) line
NSIM=1;Ea=10;Es=40;CTR=F;stoch=F;Ens=matrix(NA,nrow=(nt*TS+1),ncol=NSIM) #Output matrix - Total infected over time for each sim
source("Fig 4 ModSimsRepFinal.R")
lines(x=outD[,1],y=Ens[,1],type="l",col="black",xlab="",ylab="",yaxt="n",ylim=c(0,ht),lwd=1.5)
mtext("Day", side = 1, outer = TRUE, line = 1.7,cex=1.7)
mtext("Number of latently infected people, E", side = 2, outer = TRUE, line = 1.7,cex=1.7)
R0t[1]
dev.off()

