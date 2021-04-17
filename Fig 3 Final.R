setwd("c:/marm/research/covid-19/modeling T-CT-I/")
require(deSolve);require(tidyverse)
R0a = function(parsR) { #function for R0 expression
  with(parsR,
       kappa*beta*(S/N)*(sigma_Is*q_Ips_Is/(alpha_Is+gamma_Is+tau_Is+eps_Is_Ips*ftraced+eps_Is_Is*ftraced)*(q_E_Ips/(q_Ips_Is+eps_Ips_Ips*ftraced+eps_Ips_Is*ftraced))
                         +(sigma_Ips*q_E_Ips/(q_Ips_Is+eps_Ips_Ips*ftraced+eps_Ips_Is*ftraced)))*(1/(q_E_Ips+eps_E_Ips*ftraced+eps_E_Is*ftraced))) }

SIRb=function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {    
    capcontratio=(NCT*NCCTD)/(Is*tau_Is*Ncont_p_inf);capcontratio[capcontratio>1]=1
    ftraced=ifelse(CTR,capcontratio,0)
    
    N=S+E+Ips+Is+R+Q
    dS=-kappa*beta*S*(sigma_Ips*Ips+sigma_Is*Is)/N
    dE=kappa*beta*S*(sigma_Ips*Ips+sigma_Is*Is)/N-q_E_Ips*E-eps_E_Ips*ftraced*E-eps_E_Is*ftraced*E
    dIps=q_E_Ips*E-q_Ips_Is*Ips-eps_Ips_Ips*ftraced*Ips-eps_Ips_Is*ftraced*Ips
    dIs=q_Ips_Is*Ips-gamma_Is*Is-tau_Is*Is-eps_Is_Ips*ftraced*Is-eps_Is_Is*ftraced*Is-alpha_Is*Is
    dQ=(eps_E_Ips*ftraced+eps_E_Is*ftraced)*E+(eps_Ips_Ips*ftraced+eps_Ips_Is*ftraced)*Ips+
      (eps_Is_Ips*ftraced+eps_Is_Is*ftraced+tau_Is)*Is-q_Q_R*Q-alpha_Q*Q
    dR=gamma_Is*Is+q_Q_R*Q
    
    return(list(c(dS,dE,dIps,dIs,dQ,dR)))
  })
}

#params
#Infectiousness vs day from He et al gamma distribution parameterized w/ mean, variance
shape=function(mu,sig) mu^2/sig^2;scale = function(mu,sig) sig^2/mu;mn=3.1;sig2=2.15
presympInf=pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2))/2.3 #relative daily infectiousness pre-symptomatic
postsympInf=(pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/5 #relative daily infectiousness mildly-symptomatic
ssInf=(pgamma(13.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/6 #relative daily infectiousness severely-symptomatic

pars3 = data.frame(N=100001,
                   beta = 0.35,
                   kappa=1,
                   sigma_Ips=presympInf/postsympInf,
                   sigma_Is=1,
                   gamma_Is=1/5,
                   q_E_Ips=1/3.2,
                   q_Ips_Is=1/2.3,
                   q_Q_R=1/14,
                   tau_Is=1/5,
                   tracetime=0.5,
                   CTR=T,
                   Ncont_p_inf= 10, # contacts per infection
                   NCCTD=12, # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
                   NCT=15 # # of contact tracers for population
)

pars3 = pars3 %>% mutate (
  alpha_Is=gamma_Is*(0.0066/0.8)/(1-(0.0066/0.8)),
  alpha_Q=q_Q_R*(0.0066/0.8)/(1-(0.0066/0.8)),
  f_test=tau_Is/(tau_Is+gamma_Is),
  f_Ips=(sigma_Ips/q_Ips_Is)/
    (sigma_Ips/q_Ips_Is+sigma_Is/(gamma_Is+tau_Is)),
  f_Is=(sigma_Is/(gamma_Is+tau_Is))/
    (sigma_Ips/q_Ips_Is+sigma_Is/(gamma_Is+tau_Is)),
  ttt_Ips=1/tau_Is+1/q_Ips_Is+tracetime, #time to trace infect caused by Ips
  ttt_Is=1/tau_Is+tracetime, #time to trace infect caused by Ips
  f_E_Ips=exp(-q_E_Ips*ttt_Ips), #fraction of traced infections caused while Ips still in Ea
  f_E_Is=exp(-q_E_Ips*ttt_Is), #fraction of traced infections caused while Is still in Ea
  f_Ips_Ips=(1-f_E_Ips)*exp(-(1/(1/q_E_Ips+1/q_Ips_Is)*ttt_Ips)),
  f_Ips_Is=(1-f_E_Is)*exp(-(1/(1/q_E_Ips+1/q_Ips_Is)*ttt_Is)),
  f_Is_Ips=(1-f_E_Ips-f_Ips_Ips)*
    exp(-(1/(1/q_E_Ips+1/q_Ips_Is+1/gamma_Is)*ttt_Ips)),
  f_Is_Is=(1-f_E_Is-f_Ips_Is)*
    exp(-(1/(1/q_E_Ips+1/q_Ips_Is+1/gamma_Is)*ttt_Is)),
  eps_E_Ips=f_test*f_Ips*f_E_Ips, #traced infections while Ips to be removed
  eps_E_Is=f_test*f_Is*f_E_Is, #traced infections while Is to be removed
  eps_Ips_Ips=f_test*f_Ips*f_Ips_Ips, #traced infections while Ips to be removed
  eps_Ips_Is=f_test*f_Is*f_Ips_Is,#traced infections while Ips to be removed
  eps_Is_Ips=f_test*f_Ips*f_Is_Ips,#traced infections while Ips to be removed
  eps_Is_Is=f_test*f_Is*f_Is_Is)#traced infections while Ips to be removed

times=seq(0, 350, by = .5);yini=c(S=100000,E=0,Ips=5,Is=5,Q=0,R=0)
options(scipen=5);titlocH=225;titlocV1=99200;titlocV2=93200

#Social distancing panels A,B
parsSD= pars3 %>% mutate (kappa=0.5799608, CTR=F)
#parsSD= pars3 %>% mutate (kappa=1,tau_Is=0, CTR=F) #run to get R0 w/ no testing or tracing 3.2

dynSD=as.data.frame(ode(yini, times, SIRb, parsSD));#dynSD$I=dynSD$Ips+dynSD$Is
tiff("Fig 3.tif",width=14,height=8,units="in",res=300,compression = "lzw")
#pdf("Fig 3.pdf",width=14,height=8)
cl1=c(3,1,2,7,4,5)
par(mfcol=c(2,3),mar=c(0, 3, 0, 0),mgp=c(1.5,.5,0),oma = c(.5, .5, 0.5, 0.5),cex=1.25)
matplot(x=dynSD[,1],y=dynSD[,-1],type="l",xaxt="n",lty=1,
        xlab="Day",ylab="Number of people",col=cl1)
legend(20,.8e5,c("S","E",expression(I[ps]),expression(I[s]),"Q","R"),col=cl1,lty=1,bty="n",seg.len=1,lwd=2,x.intersp=.5)
text(0,95000,"A")
text(titlocH+5,titlocV1,"Social distancing")
text(titlocH+5,titlocV2,"No contact tracing")

parsSD2=data.frame(dynSD, parsSD)
parsSD2$ftraced=0
R0sd=R0a(parsSD2);
par(mar=c(3, 3, 0, 0) )#b,l,t,r
plot(R0sd~dynSD[,1],xlab="",ylab=expression(R[t]),type="l",ylim=c(0,2.0),col=8,lwd=2) #Rt over time
text(0,2.0,"B")

parsUCT= pars3 %>% mutate (kappa=1,CTR=T,NCT=1500)
dynUCT=ode(yini, times, SIRb, parsUCT);
par(mar=c(0, 0, 0, 0),mgp=c(1.5,.5,0),cex=1.25)
matplot(x=dynUCT[,1],y=dynUCT[,2:ncol(dynUCT)],col=cl1,type="l",xaxt="n",yaxt="n",lty=1,
        xlab="Day",ylab="Number of people")
parsUCT=data.frame(dynUCT, parsUCT)
parsUCT$capcontratio=(parsUCT$NCT*parsUCT$NCCTD)/(parsUCT$Is*parsUCT$tau_Is*parsUCT$Ncont_p_inf)
parsUCT$ftraced=parsUCT$capcontratio;parsUCT$ftraced[parsUCT$ftraced>1]=1
R0uct=R0a(parsUCT)
text(0,95000,"C")
text(titlocH,titlocV1,"No social distancing")
text(titlocH,titlocV2,"Eff. unlimited contact tracing")
par(mar=c(3, 0, 0, 0) )#b,l,t,r
plot(R0uct~dynUCT[,1],xlab="Day",cex.lab=1.2,type="l",yaxt="n",ylim=c(0,2.0),col=8,lwd=2) #Rt over time
text(0,2.0,"D")

parsLCT=pars3 %>% mutate(kappa=1,NCT=15) 
dynLCT=ode(yini, times, SIRb, parsLCT);
par(mar=c(0, 0, 0, 0),mgp=c(1.5,.5,0),cex=1.25)
matplot(x=dynLCT[,1],y=dynLCT[,2:ncol(dynLCT)],col=cl1,type="l",xaxt="n",yaxt="n",lty=1,
        xlab="Day",ylab="Number of people")
text(0,95000,"E")
text(titlocH,titlocV1,"No social distancing")
text(titlocH,titlocV2,"Limited contact tracing")

parsLCT2=data.frame(dynLCT, parsLCT)
parsLCT2$capcontratio=(parsLCT2$NCT*parsLCT2$NCCTD)/(parsLCT2$Is*parsLCT2$tau_Is*parsLCT2$Ncont_p_inf)
parsLCT2$ftraced=replace(parsLCT2$capcontratio, parsLCT2$capcontratio>1,1)
R0lct=R0a(parsLCT2)
par(mar=c(3, 0, 0, 0) )#b,l,t,r
plot(R0lct~dynLCT[,1],xlab="",type="l",yaxt="n",ylim=c(0,2.0),col=8,lwd=2) #Rt over time
text(0,2.0,"F")
dev.off()
#----------------------
R0uct[1]/R0sd[1] #ratio of starting Rt values
R0lct[1]/R0sd[1] #ratio of starting Rt values
