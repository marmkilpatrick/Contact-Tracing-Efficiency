for(j in 1:NSIM){
  for (t in 1:(nt*TS)) { #run loop to move through time
    capcontratio=(pars$NCT*pars$NCCTD)/(Is[t]*pars$tau_Is*pars$Ncont_p_inf);capcontratio[capcontratio>1]=1
    ftraced=ifelse(CTR,capcontratio,0)
    tracetime=0.5
    ft[t]=ftraced
    
    N=S+E+Ips+Is+R+Q;#N=100006 for all timesteps if alpha=0
    
    f_test=pars$tau_Is/(pars$tau_Is+pars$gamma_Is)
    f_Ips=(pars$sigma_Ips/pars$q_Ips_Is)/
      (pars$sigma_Ips/pars$q_Ips_Is+pars$sigma_Is/(pars$gamma_Is+pars$tau_Is))
    f_Is=(pars$sigma_Is/(pars$gamma_Is+pars$tau_Is))/
      (pars$sigma_Ips/pars$q_Ips_Is+pars$sigma_Is/(pars$gamma_Is+pars$tau_Is))
    ttt_Ips=1/pars$tau_Is+1/pars$q_Ips_Is+tracetime #time to trace infect caused by Ips
    ttt_Is=1/pars$tau_Is+tracetime #time to trace infect caused by Ips
    f_E_Ips=exp(-pars$q_E_Ips*ttt_Ips) #fraction of traced infections caused while Ips still in Ea
    f_E_Is=exp(-pars$q_E_Ips*ttt_Is) #fraction of traced infections caused while Is still in Ea
    f_Ips_Ips=(1-f_E_Ips)*exp(-(1/(1/pars$q_E_Ips+1/pars$q_Ips_Is)*ttt_Ips))
    f_Ips_Is=(1-f_E_Is)*exp(-(1/(1/pars$q_E_Ips+1/pars$q_Ips_Is)*ttt_Is))
    f_Is_Ips=(1-f_E_Ips-f_Ips_Ips)*
      exp(-(1/(1/pars$q_E_Ips+1/pars$q_Ips_Is+1/pars$gamma_Is)*ttt_Ips))
    f_Is_Is=(1-f_E_Is-f_Ips_Is)*
      exp(-(1/(1/pars$q_E_Ips+1/pars$q_Ips_Is+1/pars$gamma_Is)*ttt_Is))
    eps_E_Ips=f_test*f_Ips*f_E_Ips #traced infections while Ips to be removed
    eps_E_Is=f_test*f_Is*f_E_Is #traced infections while Is to be removed
    eps_Ips_Ips=f_test*f_Ips*f_Ips_Ips #traced infections while Ips to be removed
    eps_Ips_Is=f_test*f_Is*f_Ips_Is #traced infections while Ips to be removed
    eps_Is_Ips=f_test*f_Ips*f_Is_Ips #traced infections while Ips to be removed
    eps_Is_Is=f_test*f_Is*f_Is_Is #traced infections while Ips to be removed
    
    
    R0t[t]=pars$kappa*pars$beta*(S[t]/N[t])*
      (pars$q_E_Ips/(pars$q_E_Ips+eps_E_Ips*ftraced + eps_E_Is*ftraced))*
      ((pars$sigma_Ips/pars$q_Ips_Is+eps_Ips_Ips*ftraced + eps_Ips_Is*ftraced)+
         (pars$sigma_Is/(pars$alpha+pars$tau_Is+pars$gamma_Is+eps_Is_Ips*ftraced + eps_Is_Is*ftraced))*
         (pars$q_Ips_Is/(pars$q_Ips_Is+eps_Ips_Ips*ftraced + eps_Ips_Is*ftraced)) )
    
    
    mCPC=(pars$kappa*pars$beta/R0t[t]*S[t]/N[t]*
            (pars$sigma_Ips*Ips[t] + pars$sigma_Is*Is[t]))#mean # of new infections/Rt
    actualTT=(mCPC-floor(mCPC))*rnbinom(1,mu=R0t[t],size=0.16)+ #stoch # of new infections
      sum(rnbinom(floor(mCPC),mu=R0t[t],size=0.16))
    if(stoch==F) actualTT=mCPC*R0t[t] #For deterministic simulation use mean/Rt * Rt
    
    S[t+1]=S[t]-actualTT/TS
    E[t+1]=E[t]+(actualTT - pars$q_E_Ips*E[t] -eps_E_Is*ftraced*E[t]-eps_E_Ips*ftraced*E[t])/TS
    Ips[t+1]=Ips[t]+(pars$q_E_Ips*E[t]-pars$q_Ips_Is*Ips[t]- eps_Ips_Is*ftraced*Ips[t] - eps_Ips_Ips*ftraced*Ips[t])/TS
    Is[t+1]=Is[t]+(pars$q_Ips_Is*Ips[t] - (pars$alpha+pars$tau_Is + pars$gamma_Is)* Is[t] -
                       eps_Is_Is*ftraced*Is[t] - eps_Is_Ips*ftraced*Is[t])/TS
    In[t+1]=Ips[t+1]+Is[t+1]
    R[t+1]=R[t]+(pars$gamma_Is*Is[t]+pars$q_Q_R*Q[t])/TS
    Q[t+1]=Q[t]+(E[t]*ftraced*(eps_E_Is+eps_E_Ips)+Ips[t]*ftraced*(eps_Ips_Is+eps_Ips_Ips)+Is[t]*ftraced*(eps_Is_Is+eps_Is_Ips)+pars$tau_Is*Is[t]-pars$q_Q_R*Q[t])/TS
    N[t+1]=S[t+1]+E[t+1]+Ips[t+1]+Is[t+1]+R[t+1]+Q[t+1]
  }
  
  outD=data.frame(day=seq(0,nt,by=1/TS),S,E,Ips,Is,R,Q)
  Ens[,j]=E }
