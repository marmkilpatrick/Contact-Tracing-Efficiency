setwd("c:/marm/research/covid-19/modeling T-CT-I/")
lapply(c("ggplot2","tidyverse","deSolve"),require,character.only=T) #load packages

R0a = function(parsR) { #function for R0 expression
  with(parsR,
       kappa*beta*(S/N)*(sigma_Is*q_Ips_Is/(alpha+gamma_Is+tau_Is+eps_Is_Ips*ftraced+eps_Is_Is*ftraced)*(q_E_Ips/(q_Ips_Is+eps_Ips_Ips*ftraced+eps_Ips_Is*ftraced))
                         +(sigma_Ips*q_E_Ips/(q_Ips_Is+eps_Ips_Ips*ftraced+eps_Ips_Is*ftraced)))*(1/(q_E_Ips+eps_E_Ips*ftraced+eps_E_Is*ftraced))) }
#params
#Infectiousness vs day from He et al gamma distribution parameterized w/ mean, variance
shape=function(mu,sig) mu^2/sig^2;scale = function(mu,sig) sig^2/mu;mn=3.1;sig2=2.15
presympInf=pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2))/2.3 #relative daily infectiousness pre-symptomatic
postsympInf=(pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(2.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/5 #relative daily infectiousness mildly-symptomatic
ssInf=(pgamma(13.3,shape=shape(mn,sig2),scale=scale(mn,sig2))-pgamma(7.3,shape=shape(mn,sig2),scale=scale(mn,sig2)))/6 #relative daily infectiousness severely-symptomatic

Ninf=c(10^seq(0.2,3.4,by=.01)) # # of cases/infections traced
Ncont_p_inf= 10 #c(5,10,20,30) # contacts per infection - value before lockdowns 30
NCCTD=12 # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
NCT=15 # # of contact tracers for population - CA target
capcontratio=(NCT*NCCTD)/(Ninf*Ncont_p_inf);capcontratio[capcontratio>1]=1
tdelay=c(1:5,10);testdelay=rep(tdelay,each=length(Ninf))

pars1 = data.frame(N=100001,S=100001,
                   beta = 0.35,
                   kappa=1,
                   sigma_Ips=presympInf/postsympInf,
                   sigma_Is=1,
                   gamma_Is=1/5,
                   q_E_Ips=1/3.2,
                   q_Ips_Is=1/2.3,
                   q_Q_R=1/14,
                   tau_Is=1/testdelay,
                   ftraced=capcontratio,
                   tracetime=0.5,
                   CTR=T,
                   Ncont_p_inf= 10, # contacts per infection
                   NCCTD=12, # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
                   NCT=15 # # of contact tracers for population
)

pars1 = pars1 %>% mutate (
  alpha=gamma_Is*(0.0066/0.8)/(1-(0.0066/0.8)),
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

dat3=data.frame(Delay=testdelay,Ncont_p_inf=Ncont_p_inf,Ninf=rep(Ninf,length(tdelay)*length(Ncont_p_inf)),
                contacts_per_tracer_calls=c(Ninf*Ncont_p_inf/(NCT*NCCTD)),R0v=R0a(pars1))
cbp1 <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p1=ggplot(dat3, aes(y=R0v, x=contacts_per_tracer_calls,colour=factor(Delay))) + 
  theme_bw()+geom_point(size=2)+coord_cartesian(ylim=c(0.4,3.05))+
  scale_x_continuous(trans = 'log10')+
  theme(axis.title=element_text(size=25),legend.position = "none",
        axis.text=element_text(size=25))+
  xlab(expression(Ratio~of~Contacts~needing~tracing~to~Tracing~Capacity))+
  ylab(expression(R[t]))+
  annotate(geom="text",x=0.25, y=aggregate(R0v~Delay,data=dat3,FUN=min)$R0v+.05,
           label=unique(dat3$Delay),size=7)+
  annotate(geom="text",x=0.5, y=.88*aggregate(R0v~Delay,data=dat3,FUN=min)$R0v[1],
           label="Delay: symptoms->testing/tracing",size=7)+
  theme(panel.grid.minor = element_blank())+scale_colour_manual(values=cbp1);p1
#ggsave("Fig 1.pdf",plot=p1,width=12,height=8)
ggsave("Fig 1.tiff",plot=p1,width=12,height=8,units="in",compression = "lzw")

#Fig 2-------------------------------------------
Ncont_p_inf= c(5,10,20,30) # contacts per infection - value before lockdowns 30
capcontratio=(NCT*NCCTD)/(rep(Ninf,length(Ncont_p_inf))*rep(Ncont_p_inf,each=length(Ninf)))
capcontratio[capcontratio>1]=1;tracetime=0.5#fixing avg time to trace to 1/2 day
tdelay=5;testdelay=rep(tdelay,each=length(Ninf)*length(Ncont_p_inf))

pars2 = data.frame(N=100001,S=100001,
                   beta = 0.4,
                   kappa=1,
                   sigma_Ips=presympInf/postsympInf,
                   sigma_Is=1,
                   gamma_Is=1/5,
                   q_E_Ips=1/3.2,
                   q_Ips_Is=1/2.3,
                   q_Q_R=1/14,
                   tau_Is=1/testdelay,
                   ftraced=capcontratio,
                   tracetime=0.5,
                   CTR=T,
                   Ncont_p_inf= 10, # contacts per infection
                   NCCTD=12, # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
                   NCT=15 # # of contact tracers for population
)

pars2 = pars2 %>% mutate (
  alpha=gamma_Is*(0.0066/0.8)/(1-(0.0066/0.8)),
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

dat4=data.frame(Delay=testdelay,Ncont_p_inf=rep(Ncont_p_inf,each=length(Ninf)),Ninf=rep(Ninf,length(Ncont_p_inf)*length(tdelay)),
                cases_per_tracer_calls=c(rep(Ninf,length(Ncont_p_inf)*length(tdelay))/(NCT*NCCTD)))
dat4$R0v=R0a(pars2)
cbp2 <- c("#F0E442","#009E73","#0072B2", "#D55E00", "#CC79A7")
p2=ggplot(dat4, aes(y=R0v, x=cases_per_tracer_calls,color=factor(Ncont_p_inf))) + 
  theme_bw()+geom_point(size=2)+coord_cartesian(ylim=c(0.4,3.05))+
  scale_x_continuous(trans = 'log10')+
  theme(axis.title=element_text(size=25),legend.position = "none",
        axis.text=element_text(size=25))+
  xlab(expression(Cases~per~Contact~Tracer~Calls~Per~Day))+
  ylab(expression(R[t]))+
  annotate(geom="text",x=c(10^(c(-.4,-.7,-1,-1.2)-.35)), y=1.6,
           label=Ncont_p_inf,size=7)+
  annotate(geom="text",x=.09, y=1.4,
           label="Contacts/case",size=7)+
  theme(panel.grid.minor = element_blank())+scale_colour_manual(values=cbp2);p2
#ggsave("Fig 2.pdf",plot=p2,width=12,height=8)
ggsave("Fig 2.tiff",plot=p2,width=12,height=8,units="in",compression = "lzw")

#-----------------------Sensitivity analysis Fig S3
Ninf= 18# # of cases/infections traced
Ncont_p_inf=10 # contacts per infection - value before lockdowns 30
NCCTD=12 # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
NCT=15 # # of contact tracers for population - CA target
capcontratio=(NCT*NCCTD)/(Ninf*Ncont_p_inf);capcontratio[capcontratio>1]=1
pars5 = data.frame(
                   kappa=1,
                   beta = 0.35,
                   sigma_Ips=presympInf/postsympInf,
                   sigma_Is=1,
                   gamma_Is=1/5,
                   q_E_Ips=1/3.2,
                   q_Ips_Is=1/2.3,
                   tau_Is=1/5,
                   alpha=0.0025)

Rt_sens=data.frame(matrix(NA,nrow=length(pars5),ncol=2));colnames(Rt_sens)=c("dec","inc");rownames(Rt_sens)=names(pars5)
for (k in 1:length(pars5)) {
  parsR0_sens=pars5
  for (l in 1:2) {
    parsR0_sens[k]=pars5[k]*(.7+0.2*l)
    pars6 = parsR0_sens %>% mutate (N=100001,S=100001,
                              q_Q_R=1/14,
                              ftraced=capcontratio,
                              tracetime=0.5,
                              CTR=T,
                              Ncont_p_inf= 10, # contacts per infection
                              NCCTD=12, # # of contacts a contact tracer can contact per day to isolate; 1.5/hr
                              NCT=15, # # of contact tracers for population
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
    
    Rt_sens[k,l]=R0a(pars6)   
}}

parlabs=c(expression(kappa),
          expression(beta),
          expression(sigma[Ips]),
          expression(sigma[Is]),
          expression(gamma[Is]),
          expression(q[E-Ips]),
          expression(q[Ips-Is]),
          expression(tau[Is]),
          expression(alpha[Is])
)
Rtbase=sum(Rt_sens[1,])/2
Rt_sensR=(Rt_sens-Rtbase)/Rtbase
#pdf("Fig S3.pdf",width=10,height=8)
tiff("Fig S3.tif",width=10,height=8,units="in",res=300,compression = "lzw")
par(mar=c(4, 3, .5, 0.5),mgp=c(2.5,.5,0),cex=2)
barplot(Rt_sensR[,1],horiz=T,names.arg=parlabs,las=2,col="blue",xlim=c(-.11,.11),#xlab="Change in R0")
        xlab=expression(Change~"in"~R[t]))
barplot(Rt_sensR[,2],horiz=T,names.arg=parlabs,las=2,col="red",add=T)
abline(h=c(.75,1.85,3.05,4.25,5.45,6.65,7.85,9.05,10.25),col="grey",lty=3)
box();#grid(ny=nrow(Rt_sens)+1)
dev.off()
