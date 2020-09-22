lapply(c("ggplot2","ggthemes","deSolve","dplyr","gridExtra"),require,character.only=T) #load multiple packages

R0 = function(parsR01,cTrace) { #function for R0 expression
  with(parsR01,  ((kappa * beta) / (cTrace + qEA + qEPS)) *
         (((qEA * sigmaA) / (cTrace + gammaA)) + 
            ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
            ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
            ((qEPS * qPSMS * qMSSS * sigmaSS) /
               ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
                  (cTrace + ssTest + alpha + gammaSS))  )    )             )
}
cTraceF = function(FrCaseTr,Delay,Ninf,Ncont_p_inf,NCT,NCCPD,parsR0) { #function estimating cTrace parameter
  dTracing=Delay+0.5*Ninf*Ncont_p_inf/(NCT*NCCPD) #duration until traced
  return(FrCaseTr*parsR0$qEPS/(parsR0$qEPS+parsR0$qEA)*1/dTracing)}

#Parameters
f_asymp=0.20 #frac asymptomatic from Buitrago-Garcia
d_EPS=3.2; #Duration latent period for symptomatic cases He et al 2020 Nat Med in days
d_EA=1/( f_asymp*(1/d_EPS)/(1-f_asymp) ) #set duration asympt to give f_asymp fraction
d_PSMS=2.3 #duration infectious period before symptom onset for 5.5 d incubation He et al 2020 Nat Med
d_MSSS=8 #duration mild to severe symptoms diff breathing; Zhou et al 2020 Lancet; also Lewnard et al 2020 BMJ
#Infectiousness vs day from He et al 1c gamma distribution parameterized w/ mean, variance
shape=function(mu,sig) mu^2/sig^2;scale = function(mu,sig) sig^2/mu;mn=3.1;sig2=2.15
presympInf=pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2))/2.3 #relative daily infectiousness pre-symptomatic
postsympInf=(pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/5 #relative daily infectiousness mildly-symptomatic
ssInf=(pgamma(13.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/6 #relative daily infectiousness severely-symptomatic

#Contact tracing efficiency w/ infections, contact tracers, contacts; only ratio important
Ninf=c(10^seq(0,4,by=.01)) # # of cases/infections traced
Ncont_p_inf=c(5,10,20,30) # contacts per infection - value before lockdowns 30
NCCPD=12 # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
NCT=15 # # of contact tracers for population - CA target

#-----------Fig 1
dat2=data.frame(Delay=rep(1:5,each=length(Ninf)),Ncont_p_inf=10,Ninf=Ninf,
                cases_per_tracer_calls=c(Ninf/(NCT*NCCPD)))

d_msTest=rep(1:5,each=length(Ninf)) #Variable delays
parsR0_F1 = data.frame(kappa = 1,       #Social distancing factor: set to 1/2.6 to get R0=1
                       beta   = .37, # transmission rate - set to give R0 = ~2.4 
                       msTest  = 1/d_msTest ,   # Removal by testing of mildly symptomatic infected toggle w/ active surv
                       ssTest  = 1. ,   # Removal by testing of severely symptomatic infected
                       alpha  = 0.0066 ,   # IFR - Verity et al 2020 Lancet
                       qEA  = 1/d_EA ,   # Transition rate from exposed to asymptomatic infected
                       qEPS  = 1/d_EPS ,   # Transition rate from exposed to presymptomatic infected
                       qPSMS  = 1/d_PSMS ,   # Transition rate from presymptomatic to mild symptoms
                       qMSSS  = 1/d_MSSS,   # Transition rate from mildly to severely infected
                       gammaA  = 1/8 ,   # Recovery rate of asymptomatic infected - approx 8d infect period He et al
                       gammaMS  = 1/5,   # Recovery rate of mildly symptomatic infected 
                       gammaSS  = 1/10.7 ,   # Recovery rate of severely symptomatic infected - Lewnard et al 2020
                       gammaQ  = 0.14 ,   # Recovery rate of quarantined not important but 1/7d
                       sigmaA  = 1 ,   # Infectiousness of asymptomatic infected #unknown
                       sigmaPS  = presympInf/postsympInf ,   # Infectiousness of pre-symptomatic infected
                       sigmaMS  = 1 ,   # Infectiousness of mildly symptomatic infected
                       sigmaSS  = ssInf/postsympInf)     # Infectiousness of severely symptomatic infected

f_mstr=1 #fraction of mildly symptomatic cases traced
FrCaseTr=f_mstr*parsR0_F1$msTest/(parsR0_F1$msTest+parsR0_F1$qMSSS+parsR0_F1$gammaMS)/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)
dat2$R0v=R0(parsR0_F1,cTraceF(FrCaseTr,Delay=dat2$Delay,dat2$Ninf,dat2$Ncont_p_inf,NCT,NCCPD,parsR0_F1))

p1=ggplot(dat2, aes(y=R0v, x=cases_per_tracer_calls,col=R0v)) + 
  theme_bw()+geom_point(size=2)+scale_x_continuous(trans='log10')+
  coord_cartesian(xlim=c(.005, 50))+
  scale_color_gradientn(colours=c("green","yellow","orange","red"),limits=c(0.2,2.5))+
  theme(axis.title=element_text(size=25),legend.position = c(.8, .50),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+labs(colour = expression(R[t]))+
  geom_hline(yintercept=1)+
  xlab(expression(Cases~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(R[t]))+
  annotate(geom="text",x=0.01, y=aggregate(R0v~Delay,data=dat2,FUN=min)$R0v,
           label=unique(dat2$Delay),size=5)+
  annotate(geom="text",x=0.016, y=1.1*aggregate(R0v~Delay,data=dat2,FUN=min)$R0v[5],
         label="Delay: symptoms->test/tracing",size=6);p1

#Fig 2--------------------------
d_msTest=5 #delay from symptom onset to test/trace ~5d from recent data in CA
parsR0 = data.frame(kappa = 1,       #Social distancing factor
         beta   = .37, # transmission rate - set to give R0 = ~2.4
         msTest  = 1/d_msTest ,   # Removal rate by testing of mildly symptomatic infected
         ssTest  = 1. ,   # Removal by testing of severely symptomatic infected
         alpha  = 0.0066 * (1/5 + 1/10.7 + 1/d_EA)/(1-0.0066), # set to give IFR=0.66% - Verity et al 2020 Lancet
         qEA  = 1/d_EA ,   # Transition rate from exposed to asymptomatic infected
         qEPS  = 1/d_EPS ,   # Transition rate from exposed to presymptomatic infected
         qPSMS  = 1/d_PSMS ,   # Transition rate from presymptomatic to mild symptoms
         qMSSS  = 1/d_MSSS,   # Transition rate from mildly to severely infected
         gammaA  = 1/8 ,   # Recovery rate of asymptomatic infected - approx 8d infect period He et al
         gammaMS  = 1/5,   # Recovery rate of mildly symptomatic infected 
         gammaSS  = 1/10.7 ,   # Recovery rate of severely symptomatic infected - Lewnard et al 2020
         gammaQ  = 1/14 ,   # Recovery rate of quarantined not important but 1/7d
         sigmaA  = 1 ,   # Infectiousness of asymptomatic infected #unknown
         sigmaPS  = presympInf/postsympInf ,   # Infectiousness of pre-symptomatic infected
         sigmaMS  = 1 ,   # Infectiousness of mildly symptomatic infected
         sigmaSS  = ssInf/postsympInf)     # Infectiousness of severely symptomatic infected

R0(parsR0,cTrace=0) #R0 w/ no tracing

f_mstr=1 #fraction of mildly symptomatic cases traced
FrCaseTr=f_mstr*(parsR0$msTest/(parsR0$msTest+parsR0$qMSSS+parsR0$gammaMS))/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)

dat=data.frame(Delay=d_msTest,Ncont_p_inf=rep(Ncont_p_inf,each=length(Ninf)),Ninf=Ninf,
               cases_per_tracer_calls=c(Ninf/(NCT*NCCPD)),contacts_per_tracer_calls_per_day=c(Ninf*rep(Ncont_p_inf,each=length(Ninf))/(NCT*NCCPD)))
dat$CTR_eps=cTraceF(FrCaseTr,Delay=dat$Delay,dat$Ninf,dat$Ncont_p_inf,NCT,NCCPD,parsR0) #contact tracing removal rate
dat$R0v=R0(parsR0,cTraceF(FrCaseTr,Delay=dat$Delay,dat$Ninf,dat$Ncont_p_inf,NCT,NCCPD,parsR0)) #R0 values
ylh=R0(parsR0,cTraceF(FrCaseTr,Delay=dat$Delay[1],10^2.25,Ncont_p_inf,NCT,NCCPD,parsR0)) #heights for labels of # contacts/case

p2=ggplot(dat, aes(y=R0v, x=cases_per_tracer_calls,col=R0v)) + 
  theme_bw()+geom_point(size=2)+scale_x_continuous(trans='log10')+
  coord_cartesian(xlim=c(.005, 50), ylim = c(0.9,2.6))+
  scale_color_gradientn(colours=c("green","yellow","orange","red"),limits=c(0.2,2.5))+
  theme(axis.title=element_text(size=25),legend.position = c(.8, .30),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+labs(colour = expression(R[t]))+
  geom_hline(yintercept=1) +
  xlab(expression(Cases~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(R[t]))+
  annotate(geom="text",x=1., y=ylh,
           label=Ncont_p_inf,size=5)+
  annotate(geom="text",x=.9, y=1.07*ylh[4],
           label="Contacts/case",size=5);p2

#Fig 3
f_mstr=0.5 #fraction of mildly symptomatic cases traced
FrCaseTr=f_mstr*(parsR0$msTest/(parsR0$msTest+parsR0$qMSSS+parsR0$gammaMS))/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)

dat=data.frame(Delay=d_msTest,Ncont_p_inf=rep(Ncont_p_inf,each=length(Ninf)),Ninf=Ninf,
               cases_per_tracer_calls=c(Ninf/(NCT*NCCPD)),contacts_per_tracer_calls_per_day=c(Ninf*rep(Ncont_p_inf,each=length(Ninf))/(NCT*NCCPD)))
dat$CTR_eps=cTraceF(FrCaseTr,Delay=dat$Delay,dat$Ninf,dat$Ncont_p_inf,NCT,NCCPD,parsR0) #contact tracing removal rate
dat$R0v=R0(parsR0,cTraceF(FrCaseTr,Delay=dat$Delay,dat$Ninf,dat$Ncont_p_inf,NCT,NCCPD,parsR0)) #R0 values
ylh=R0(parsR0,cTraceF(FrCaseTr,Delay=dat$Delay[1],10^2.25,Ncont_p_inf,NCT,NCCPD,parsR0)) #heights for labels of # contacts/case

p3=ggplot(dat, aes(y=R0v, x=cases_per_tracer_calls,col=R0v)) + 
  theme_bw()+geom_point(size=2)+scale_x_continuous(trans='log10')+
  coord_cartesian(xlim=c(.005, 50), ylim = c(0.9,2.6))+
  scale_color_gradientn(colours=c("green","yellow","orange","red"),limits=c(0.2,2.5))+
  theme(axis.title=element_text(size=25),legend.position = c(.8, .30),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+labs(colour = expression(R[t]))+
  geom_hline(yintercept=1) +
  xlab(expression(Cases~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(R[t]))+
  annotate(geom="text",x=1., y=ylh,
           label=Ncont_p_inf,size=5)+
  annotate(geom="text",x=1., y=1.05*ylh[4],
           label="Contacts/case",size=5);p3

#grid.arrange(p1,p2, p3, nrow = 2,as.table=F)
#--------------
#Not in paper: Plot w/ contacts per tracer calls as X-axis - all data on one line/curve
ggplot(dat, aes(y=R0v, x=contacts_per_tracer_calls_per_day,col=R0v)) + 
  theme_bw()+geom_point(size=2)+scale_x_continuous(trans='log10')+
  scale_color_gradientn(colours=c("green","yellow","orange","red"))+
  theme(axis.title=element_text(size=25),legend.position = c(.8, .30),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+labs(colour = expression(R[0]))+
  geom_hline(yintercept=1) +
  xlab(expression(Case-contacts~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(R[0]))

#Not in paper: Plot epsilon (contact tracing removal rate) vs cases/CT calls/day
ggplot(dat, aes(y=CTR_eps, x=cases_per_tracer_calls,col=CTR_eps)) + 
  theme_bw()+geom_point(size=2)+scale_x_continuous(trans='log10')+
  scale_color_gradientn(colours=c("green","yellow","orange","red"))+
  theme(axis.title=element_text(size=25),legend.position = c(.8, .70),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+labs(colour = expression(epsilon))+
  xlab(expression(Cases~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(epsilon))+
  annotate(geom="text",x=1, y=c(.019, .025, .038, .051),
           label=Ncont_p_inf,size=5)# # contacts per case

#Not in paper: Plot epsilon (contact tracing removal rate)w/ contacts per tracer calls as X-axis - all data on one line/curve
ggplot(dat, aes(y=CTR_eps, x=contacts_per_tracer_calls_per_day,col=CTR_eps)) + 
  theme_bw()+geom_point(size=2)+scale_x_continuous(trans='log10')+
  scale_color_gradientn(colours=c("green","yellow","orange","red"))+
  theme(axis.title=element_text(size=25),legend.position = c(.8, .70),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+labs(colour = expression(epsilon))+
  xlab(expression(Case-contacts~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(epsilon))

#-----------------------Sensitivity analysis Fig S2
Ninf= 18# # of cases/infections traced
Ncont_p_inf=10 # contacts per infection - value before lockdowns 30
NCCPD=12 # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
NCT=15 # # of contact tracers for population - CA target
f_mstr=1
FrCaseTr=f_mstr*(parsR0$msTest/(parsR0$msTest+parsR0$qMSSS+parsR0$gammaMS))/(1-f_asymp) #fraction of symptomatic cases traced (numerator is frac of total cases)

Rt_sens=data.frame(matrix(NA,nrow=length(parsR0),ncol=2));colnames(Rt_sens)=c("dec","inc");rownames(Rt_sens)=names(parsR0)
for (k in 1:length(parsR0)) {
  parsR0_sens=parsR0
  for (l in 1:2) {
    parsR0_sens[k]=parsR0[k]*(.7+0.2*l)
      Rt_sens[k,l]=R0(parsR0_sens,
      cTrace=cTraceF(FrCaseTr,Delay=1/parsR0_sens$msTest,Ninf=Ninf,Ncont_p_inf,NCT,NCCPD,parsR0_sens)) 
  }
}

parlabs=c(expression(kappa),
          expression(beta),
          expression(tau[ms]),
          expression(tau[ss]),
          expression(alpha),
          expression(q[E-a]),
          expression(q[E-ps]),
          expression(q[ps-ms]),
          expression(q[ms-ss]),
          expression(gamma[a]),
          expression(gamma[ms]),
          expression(gamma[ss]),
          expression(gamma[Q]),
          expression(sigma[a]),
          expression(sigma[ps]),
          expression(sigma[ms]),
          expression(sigma[ss]) )

Rt_sens=Rt_sens/R0(parsR0,cTrace=cTraceF(FrCaseTr,Delay=1/parsR0$msTest,Ninf=Ninf,Ncont_p_inf,NCT,NCCPD,parsR0))-1
par(mar=c(4, 3, .5, 0.5),mgp=c(2.5,.5,0),cex=2)
barplot(Rt_sens[,1],horiz=T,names.arg=parlabs,las=2,col="blue",xlim=c(-.11,.11),#xlab="Change in R0")
        xlab=expression(Change~"in"~R[0]))
barplot(Rt_sens[,2],horiz=T,names.arg=parlabs,las=2,col="red",add=T)
box();grid(ny=nrow(Rt_sens)+1)

