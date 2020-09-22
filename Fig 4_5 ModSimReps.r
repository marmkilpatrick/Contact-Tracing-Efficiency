for(j in 1:NSIM){
  for (i in 1:(nt*TS)) { #run loop to move through time - nt = # days; TS = #Inverse of time step in days; TS=1 => 24 hr time step; TS=6 => 4 hr time step
    N_I_tr[i]=FrCaseTr*IPSt[i]*qPSMS/TS #number cases being traced = new symptomatic cases
    cTrace=ifelse(CTR==F,0,cTraceF1(FrCaseTr,d_msTest,N_I_tr[i],Ncont_p_inf,NCT,NCCPD/TS)) #trace rate using function
    mCPC=(kappa * beta/R0t[i] * St[i] / Nt[i] *  #mean # of new infections/Rt
            (sigmaA * IAt[i] + sigmaPS * IPSt[i] + sigmaMS * IMSt[i] + sigmaSS * ISSt[i]))
    actualTT=(mCPC-floor(mCPC))*rnbinom(1,mu=R0t[i],size=0.16)+ #stoch # of new infections
      sum(rnbinom(floor(mCPC),mu=R0t[i],size=0.16))
    if(stoch==F) actualTT=mCPC*R0t[i] #For deterministic simulation use mean/Rt * Rt
    St[i+1]=St[i] - actualTT/TS
    Et[i+1]= Et[i] + (actualTT - Et[i] * (cTrace + qEA + qEPS)) /TS
    IAt[i+1]=IAt[i] + ((Et[i] * qEA) - (IAt[i] * (cTrace + gammaA)))/TS
    IPSt[i+1]=IPSt[i] + ((Et[i] * qEPS) - (IPSt[i] * (cTrace + qPSMS)))/TS
    IMSt[i+1]= IMSt[i] + ((IPSt[i] * qPSMS) - (IMSt[i] * (cTrace + msTest + qMSSS + gammaMS)))/TS
    ISSt[i+1]= ISSt[i] + ((IMSt[i] * qMSSS) - (ISSt[i] * (cTrace + ssTest + alpha + gammaSS)))/TS
#    In[i+1]= Et[i+1]
    In[i+1]= IAt[i+1]+IPSt[i+1]+IMSt[i+1]+ISSt[i+1]
    Rt[i+1]= Rt[i] + ((IAt[i] * gammaA) + (IMSt[i] * gammaMS) + (ISSt[i] * gammaSS) + (Qt[i] * gammaQ))/TS
    Qt[i+1]= Qt[i] + (IAt[i] * cTrace + Et[i] * cTrace + IPSt[i] * cTrace + IMSt[i] * (cTrace + msTest) + ISSt[i] * (cTrace + ssTest) - Qt[i] * gammaQ)/TS
    Dt[i+1]= Dt[i] + (ISSt[i] * alpha)/TS
    Nt[i+1]= St[i+1]+Et[i]+IAt[i+1]+IPSt[i+1]+IMSt[i+1]+ISSt[i+1]+Rt[i+1]+Qt[i]
    R0t[i+1] =((kappa * beta)*St[i]/Nt[i] / (cTrace + qEA + qEPS)) *
      (((qEA * sigmaA) / (cTrace + gammaA)) + 
         ((qEPS * sigmaPS) / (cTrace + qPSMS)) +
         ((qEPS * qPSMS * sigmaMS) / ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS))) +
         ((qEPS * qPSMS * qMSSS * sigmaSS) /
            ((cTrace + qPSMS) * (cTrace + msTest + qMSSS + gammaMS) * 
               (cTrace + ssTest + alpha + gammaSS))))
  } 
  Ens[,j]=Et }

