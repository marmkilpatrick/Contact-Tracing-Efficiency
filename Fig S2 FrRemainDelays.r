t_ms=1/c(1:15) #testing delays
q_E_Ips=1/3.2
q_Ips_Is=1/2.3
gamma_Is=1/5
f_Eps=exp(-q_E_Ips*1/t_ms);#f_Eps #frac remaining in E that hasn't left; decays w/ time
#     frac left E*    frac hasn't left Ips
f_Ips=(1-f_Eps)*(exp(-(1/t_ms)/(1/q_E_Ips+1/q_Ips_Is)));#f_Ips
f_ms1=(1-f_Ips-f_Eps)*(exp(-(1/t_ms)/(1/q_E_Ips+1/q_Ips_Is+1/gamma_Is)));#f_ms1
f_R=(1-f_Ips-f_Eps-f_ms1);#f_R
tab1=data.frame(delay=1/t_ms,E=f_Eps,I_ps=f_Ips,I_s=f_ms1,R=f_R);#tab1
#pdf("Fig S2.pdf",width=10,height=8)
tiff("Fig S2.tif",width=10,height=8,units="in",res=300,compression = "lzw")
par(mfcol=c(1,1),mar=c(2.5, 2.5, 0, 0),mgp=c(1.5,.5,0),oma = c(rep(0.1,4)),cex=2)
matplot(x=1/t_ms,y=cbind(f_Eps,f_Ips,f_ms1,f_R),type="b",ylim=c(0,1),pch=19,col=c(1,2,7,4),
        xlab="Delay from infection to tracing",ylab="Fraction of infected contacts in each class")
legend(5,1,c("E",expression(I[ps]),expression(I[s]),"R"),col=c(1,2,7,4),lty=1:4,bty="n",pch=19,x.intersp=.25,lwd=2)
dev.off()