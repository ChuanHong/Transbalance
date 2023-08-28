ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
  {
    out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
    if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
    mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
    for(k in 1:pp)
      {
       yy = yy0; 
       if(!is.null(fpr0)){
         tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
         fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
         TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
         TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
          yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
         FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
       }else{
         TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
         FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
       }
       out.yy = cbind(out.yy, yy)
       out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
       out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
       PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
       out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
       #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
       #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
       AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
       out.AUC <- c(out.AUC, AUC)
     }
    out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
    out
  }