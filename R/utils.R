#' @export
g.logit = function(xx){exp(xx)/(exp(xx)+1)}
#' @export
logit = function(xx){log(xx/(1-xx))}
#' @export
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}
#' @export
Intercept.GLM = function(yi,bxi){glm(yi~bxi,family="binomial")}
#' @export
ROC.FUN.Sup.Par = function(St,Yt,fpr=jump.u){
  mhati = g.logit(St); ss = unique(sort(St)); mu1 = mean(Yt); mu0 = 1-mu1; nt = length(Yt)
  S1.ss = sum.I(ss, "<=", St, mhati)/sum(mhati)
  S0.ss = sum.I(ss, "<=", St, 1-mhati)/sum(1-mhati); auc = sum(S1.ss[-1]*(S0.ss[-nt]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=g.logit(ss),"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr)$y}); colnames(out)=tmpnm
  list(auc,out)
}

#' @export
ROC.FUN.SSL.NP = function(Sv,St,Yt,fpr=jump.u,yes.CV=F,Xt=NULL,Xv=NULL,rep=10,regularize=yes.regularize){
  nv = length(Sv); nt = length(Yt); Pt = sum.I(St, ">=", Sv)/nv
  Pv = sum.I(Sv,">=",Sv)/nv; bw = sd(Pt)/nt^0.3; ss = sort(Pv)
  mhat.v = predict(locfit(Yt~lp(Pt, deg=0, h=bw),ev=Pv))
  mu1 = mean(Yt); mu0 = 1-mu1; 
  if(!yes.CV){
    S1.ss = sum.I(ss, "<=", Pv, mhat.v)/sum(mhat.v)
    S0.ss = sum.I(ss, "<=", Pv, 1-mhat.v)/sum(1-mhat.v); 
  }else{
    K.fold = 2; nt.t = floor(nt/K.fold); S1.ss = S0.ss = matrix(NA,nrow=length(ss),ncol=K.fold*rep)
    for(rr in 1:rep){
      dat.t = dat.t[sample(1:nt),]
      for(kk in 1:K.fold){
        indt.v = 1:nt.t + (kk-1)*nt.t; indt.t = setdiff(1:nt, indt.v)
        bhat.t = Est.ALASSO.GLM(dat.t[indt.t,],rtn="EST",regularize=regularize)[1+0:ncol(Xv)];
        bhat.v = Est.ALASSO.GLM(dat.t[indt.v,],rtn="EST",regularize=regularize)[1+0:ncol(Xv)];
        mvi = g.logit(cbind(1,Xv)%*%bhat.t); Pvi = sum.I(cbind(1,Xv)%*%bhat.v,">=",Sv)/nv
        S1.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Pvi,mvi)/sum(mvi)
        S0.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Pvi,1-mvi)/sum(1-mvi)
      }}
    S1.ss = apply(S1.ss,1,mean); S0.ss = apply(S0.ss,1,mean)
  }  
  auc = sum(S1.ss[-1]*(S0.ss[-nv]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=ss,"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr,rule=2)$y}); colnames(out)=tmpnm
  list(auc,out)
}
#' @export
ROC.FUN.SSL.Par = function(Xv,dat.t,fpr=jump.u,yes.CV=F,rep=10,wgt=NULL,rtn="list",bptb=NULL,bptb2=NULL,regularize=yes.regularize){
  if(is.null(wgt)){wgt=rep(1,nrow(dat.t))}
  if(is.null(bptb)){bptb = Est.ALASSO.GLM(dat.t,rtn="EST",regularize=regularize,Wi=wgt)[1+0:ncol(Xv)]}
  Yt = dat.t[,1]; if(is.null(bptb2)){bptb2 = bptb}; 
  Sv = c(cbind(1,Xv)%*%bptb); mhat.v = g.logit(c(cbind(1,Xv)%*%bptb2)); 
  ss = unique(sort(Sv)); mu1 = mean(mhat.v); mu0 = 1-mu1; nv = length(Sv); nt = length(Yt)
  if(!yes.CV){
    S1.ss = sum.I(ss, "<=", Sv, mhat.v)/sum(mhat.v)
    S0.ss = sum.I(ss, "<=", Sv, 1-mhat.v)/sum(1-mhat.v); 
  }else{
    K.fold = 2; nt.t = floor(nt/K.fold); S1.ss = S0.ss = matrix(NA,nrow=length(ss),ncol=K.fold*rep)
    for(rr in 1:rep){
      dat.t = dat.t[sample(1:nt),]
      for(kk in 1:K.fold){
        indt.v = 1:nt.t + (kk-1)*nt.t; indt.t = setdiff(1:nt, indt.v)
        bhat.t = Est.ALASSO.GLM(dat.t[indt.t,],rtn="EST",regularize=regularize,BIC.power=0.1)[1+0:ncol(Xv)];
        bhat.v = Est.ALASSO.GLM(dat.t[indt.v,],rtn="EST",regularize=regularize,BIC.power=0.1)[1+0:ncol(Xv)];
        mvi = g.logit(cbind(1,Xv)%*%bhat.t); Svi = cbind(1,Xv)%*%bhat.v
        S1.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Svi,mvi)/sum(mvi)
        S0.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Svi,1-mvi)/sum(1-mvi)
    }}
    S1.ss = apply(S1.ss,1,mean); S0.ss = apply(S0.ss,1,mean)
  }  
  auc = sum(S1.ss[-1]*(S0.ss[-nv]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=g.logit(ss),"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr,rule=2)$y}); colnames(out)=tmpnm
  if(rtn=="vec"){return(c(auc,out))}else{return(list(auc,out))}
}


## ================================================================================================================== ##
## estimate beta w/ adaptive LASSO regularization (if regularize=T) or standard logistic regression (if regularize=F) ##
## data: 1st column y; remaining x; nopen.ind indexes which subset of x should not be penalized                       ##
## Wi: weights for resampling if interested in obtaining standard errors                                              ## 
## ================================================================================================================== ##
#' @export
Est.ALASSO.GLM = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
    bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family="binomial"))); ##print(bini)
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    ## glmpath provides solution path for a range of penalty parameters ##
    tmpfit = glmpath(x.t,y,nopenalty.subset=nopen.ind,family=binomial,weight=Wi,standardize=F,min.lambda=0,lambda2=lam.ridge)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda")
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(predict(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+min(sum(y)^0.1,log(sum(y)))*df.all 
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=binomial,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(out)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}
#' @export
SIM.FUN = function(nn,rtn="data.t")
  {
	xx = mvrnorm(nn,mu=rep(0,p.x),Sigma=Sig0.X); 
	icd.B = rbinom(nn,size=1,prob=g.logit(xx[,1]*3+2))
	xx[,1] = (xx[,1]*(xx[,1]>0)+(rexp(nn,rate=0.1)+5)*rbinom(nn,size=1,prob=0.1))*icd.B
	prob.x = g.logit(-alp0+c(xx%*%beta0)-3*(1-icd.B)+0.2*(xx[,1]>15))
	yy = rbinom(nn,prob=prob.x,size=1); dat = cbind(yy,xx)
	if(rtn=="data.t"){return(dat)}else{
		zz = rbinom(nn, size=2, prob=g.logit(log(maf)+gam.z*yy))
		return(list(dat,cbind("D"=yy,"P.x"=prob.x,"G"=zz)))}
  }


#' @export
ROC.FUN.ALASSO.632boot  <- function(data, ind=NULL, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,
                            FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T){
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}    
  n.set <- 1; pp = ncol(data);  yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  ## === Apparent Accuracy === ##
  betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
  if(length(betahat)>1){
    yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp])    
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV&q.nonzero>0){
      ## === 0.632 bootstrap === ##
      roc.cv = NULL; 
      for(i in 1:rep){
        tmpind=sample(1:nn,replace=T); ind.v = setdiff(1:nn,unique(tmpind))
        dat.t = data[tmpind,]; wgti.t = wgti[tmpind]; dat.v = data[ind.v,];  wgti.v = wgti[ind.v]
        beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti.t,rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
        if(length(beta.t)>1){
          beta.t = beta.t[1:pp]; yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t); yyi.vmat[ind.v,i]=yyi.v
          roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti.v),silent=T)
          if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}  } }
       roc.cv = apply(roc.cv,1,mean)*0.632+rochat*0.368  
     }else{roc.cv = NULL}
   }else{roc.cv=beta.hat=NULL}
   if(rtn=="EST"){
     if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
     return(out)
   }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, "beta"=betahat, "Si"=yyi, "yyi.v"=yyi.vmat))}
}
#' @export
ROC.FUN.ALASSO.cv  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,
                            FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T,lam.ridge=0)
  {
    ## ========================================================================= ##
    ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
    ## ==== 1st column is binary disease status;                            ==== ##
    ## ==== 2nd column and on are individual markes                         ==== ##
    ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
    ## ========================================================================= ##

    nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)}; 
    if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}    
    n.set <- 1; pp = ncol(data)

    ## ========================= ##
    ## === Apparent Accuracy === ##
    ## ========================= ##

  	betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind,lam.ridge=lam.ridge),silent=T)
    if(length(betahat)>1)
      {
        yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp]); 
        q.nonzero = sum(abs(betahat[2:(pp)])!=0)
        if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
        if(yes.seed){set.seed(315)}
        if(yes.CV & q.nonzero >0)
        {
            ## =============================== ##
            ## === K-fold cross validation === ##
            ## =============================== ##
            K = 2; nv = floor(nn/K); bias.cv = NULL
            for(i in 1:rep)
              {
                tmpind=sample(1:nn); data = data[tmpind,]; wgti = wgti[tmpind]
                for(k in 1:K)
                 {
                    ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
                    dat.t = data[ind.t,]; dat.v = data[ind.v,]            
                    ## ============================================================================== ##
                    ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
                    ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
                    ## ============================================================================== ##
                    beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti[ind.t],rtn="EST",regularize=regularize,nopen.ind=nopen.ind,lam.ridge=lam.ridge),silent=T)
                    if(length(beta.t)>1)
                      {
                        beta.t = beta.t[1:pp]
                        yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t)
                        yyi.t = g.logit(cbind(1,as.matrix(dat.t[,-1]))%*%beta.t)
                        bias.k = try(ROC.Est.FUN(dat.t[,1],yyi.t,yy0,FPR0,wgti[ind.t])-ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti[ind.v]),silent=T)
                        if(length(bias.k)>1){bias.cv = cbind(bias.cv, bias.k)}
                      }
                }
              }
            print(ncol(bias.cv)); bias.cv = apply(bias.cv,1,mean,trim=0.05,na.rm=T)/2; roc.cv = rochat - bias.cv  
        }else{
            roc.cv = NULL
            }
     }else{
      roc.cv=beta.hat=NULL
     }
   if(rtn=="EST"){
     if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
     return(out)
   }else{return(list(rochat, roc.cv, betahat, yyi))}
  }
#' @export
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
#' @export
AUC.FUN = function(data)
  {
	dd = data[,1]; xx = data[,2]; n0 = sum(1-dd); n1 = sum(dd) 
	x0 = xx[dd==0]; x1 = xx[dd==1]
  	sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
  }
#' @export
chat.FUN.632boot <- function(data,fpr0=0.05,betahat,regularize=T,nopen.ind=NULL)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data); data = as.matrix(data); nn = length(yy)
    set.seed(1202); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:20)
      { 
      	ind.t = sample(1:nn,rep=T); ind.v = setdiff(1:nn,ind.t)
        datat = data[ind.t,]; datav = data[ind.v,]
        beta = try(Est.ALASSO.GLM(datat,nopen.ind=nopen.ind,regularize=regularize,rtn="EST")[1:p0],silent=T)
        if(length(beta)>1){
            phatv = (cbind(1,as.matrix(datav[,-1]))%*%beta)
            c.fpr   = as.numeric(try(quantile(phatv[datav[,1]==0],1-fpr0),silent=T))
            c.all = cbind(c.all,c.fpr); ##c.t.all = c(c.t.all, c.fpr.t);  
         }
      }
    phat0 = c(cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat0[data[,1]==0],1-fpr0)
    c.fpr.bc  = 0.368*c.fpr.0 + 0.632*apply(as.matrix(c.all),1,mean,na.rm=T)     
    c(g.logit(c.fpr.0),g.logit(c.fpr.bc))
  }
#' @export
chat.FUN <- function(data,fpr0=0.05,betahat,regularize=T,nopen.ind=NULL)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data)
    phat.0 = (cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat.0[data[,1]==0],1-fpr0)
    set.seed(1202); nt = floor(nrow(data)/2); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:10)
      { 
        data = data[sample(1:nrow(data)),]
        for(k in 1:2)
          {
            tmpind = 1:nt+(k-1)*nt
            datat = data[-tmpind,]; datav = data[tmpind,]
        	beta = try(Est.ALASSO.GLM(datat,nopen.ind=nopen.ind,regularize=regularize,rtn="EST")[1:p0],silent=T)
            if(length(beta)>1){
                phatv = (cbind(1,as.matrix(datav[,-1]))%*%beta)
                phatt = (cbind(1,as.matrix(datat[,-1]))%*%beta)
                c.fpr   = quantile(phatv[datav[,1]==0],1-fpr0)
                c.fpr.t = quantile(phatt[datat[,1]==0],1-fpr0) 
                c.all = rbind(c.all,c.fpr); c.t.all = rbind(c.t.all, c.fpr.t)
                tmpb.all = cbind(tmpb.all,beta)}
          }
      }
    c.fpr.bc  = c.fpr.0 + apply(c.all-c.t.all,2,mean,trim=0.05)
    g.logit(c(c.fpr.0,c.fpr.bc))
  }
#' @export
Predict.FUN <- function(data,newdata,fpr0=0.05,betahat)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data)
    set.seed(1202); nt = floor(nrow(data)/3); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:10)
      { 
        data = data[sample(1:nrow(data)),]
        for(k in 1:3)
          {
            tmpind = 1:nt+(k-1)*nt
            datat = data[-tmpind,]; datav = data[tmpind,]
            beta = try(Est.ALASSO.GLM(datat,nopen.ind=c(1,2))[[1]][1:p0])
            if(length(beta)>1){
                phatv = g.logit(cbind(1,as.matrix(datav[,-1]))%*%beta)
                phatt = g.logit(cbind(1,as.matrix(datat[,-1]))%*%beta)
                c.fpr   = quantile(phatv[datav[,1]==0],1-fpr0)
                c.fpr.t = quantile(phatt[datat[,1]==0],1-fpr0) 
                c.all = c(c.all,c.fpr); c.t.all = c(c.t.all, c.fpr.t)
                tmpb.all = cbind(tmpb.all,beta)}
          }
      }
    #pnew.all = g.logit(cbind(1,as.matrix(newdata[,-1]))%*%tmpb.all)
    pnew.0 = g.logit(cbind(1,as.matrix(newdata[,-1]))%*%betahat[1:p0])
    phat.0 = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat.0[data[,1]==0],1-fpr0)
    c.fpr.bc  = c.fpr.0 + mean(c.all-c.t.all)
    data.frame("patient_num"=newdata$patient,"Prob.RA"=round(pnew.0,5),
    "Cut-off.0"=round(c.fpr.0,5), "Cut-off.bc"=round(c.fpr.bc,5))
  }
#' @export
ProbTestPos <- function(tpr,fpr,prev)
  {
    tpr*prev + fpr*(1-prev)
  }
      
#' @export    
L2Norm <- function(data,coef,intercept=F)
  {
    yy = data[,1]; xx = data[,-1]; if(intercept){xx=cbind(1,xx)}
    apply((yy - xx%*%t(coef))^2,2,mean)
  }

#' @export                        
PPV.FUN <- function(fpr,SE, mu0){ 1/(1+fpr/SE*(1-mu0)/mu0)}
#' @export
NPV.Project <- function(npv.e0y0,p.e1,p.y0.e0,npv.e1=1)
  {
  	## P(D=0 | E=1 or E=0&Y=0) = {P(D=0|E=1)P(E=1)+P(D=0|E=0&Y=0)P(E=0&Y=0)}/{P(E=1)+P(E=0&Y=0)
  	npv.all = (npv.e1 * p.e1 + npv.e0y0*p.y0.e0*(1-p.e1))/(p.e1+p.y0.e0*(1-p.e1))
  	npv.all
  }
#' @export
S.FUN <- function(yy,Yi,Di,yes.smooth=F)
  {
  	if(yes.smooth){
		Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
		c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  	}else{
		return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  	}
    ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
  }
#' @export
Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
  {
    yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
    return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
  }

#' @export
CV.FUN <- function(data)
  {
    y = data[,1]; x = data[,-1]; nn=length(y)
    w = 1/abs(lm(y~x)$coef[-1]); x.w = x/VTM(w,nrow(x))
    s0.all = (1:100)/100
    fit.all = l1ce(y~x.w, standardize=F, bound=s0.all)
    K = 5; L2cv = NULL
    for(k in 1:K)
      {
        indv = 1:floor(nn/K) + (k-1)*floor(nn/K)
        indt = setdiff(1:nn,indv)
        fitk = l1ce(y~x.w,subset=indt, standardize=F, bound=s0.all)
        L2cv = rbind(L2cv, L2Norm(cbind(y,x.w)[indv,],coef(fitk),intercept=T))
      }
    L2cv = apply(L2cv,2,mean)
    s0 = min(s0.all[L2cv==min(L2cv)])
    bhat  = l1ce(y~x.w, standardize=F, bound=s0)
    list("b"=coef(bhat)/w, "s0"=s0, "lam0" =bhat$Lagrangian,"b0"=bhat$bound)    
  }
#' @export
VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}
#' @export
sum.I <- function(yy,FUN,Yi,Vi=NULL)
## sum_i I(yy FUN Yi)Vi
# Vi weight
  {
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    # for each distinct ordered failure time t[j], number of Xi < t[j]
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
    if (!is.null(Vi)) {
	   ## if FUN contains '=', tmpind is the order of decending
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
}
#' @export
convert <- function(fit) {
	rochat.auc = fit$rochat[1]
	rochat.values = matrix(fit$rochat[-1],ncol=6)
	colnames(rochat.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
	roc.cv.auc = fit$roc.cv[1]
	roc.cv.values = matrix(fit$roc.cv[-1],ncol=6)
	colnames(roc.cv.values) = c("cutoff","pos.rate","FPR","TPR","PPV","NPV")
	betahat = fit$beta[1:(grep("\\<bini",names(fit$beta))[1]-1)]
	names(betahat) = gsub("\\<b\\.","",names(betahat))
	return(list(ROC.hat.auc=rochat.auc, ROC.hat.values=rochat.values, ROC.cv.auc=roc.cv.auc,ROC.cv.values=roc.cv.values, beta=betahat, Y.hat=fit$Si))
}
#' @export
AUC = function(D,score,wgt=NULL){
	if(is.null(wgt)) {wgt=rep(1,length(D))}
	auc = sum(S.FUN(score,Yi=score,D*wgt,yes.smooth=F)*(1-D)*wgt)/sum((1-D)*wgt)
	return(auc)
}
#' @export
ROC = function(D,score){
	roc = ROC.Est.FUN(D,score,0.5,seq(.01,.99,by=.01))
	roc = matrix(roc[-1],ncol=6)
	colnames(roc) = c("cutoff","est.pos.rate","FPR","TPR","PPV","NPV")
	return(roc)
}


#' @export
intercept.fun=function(n, p, beta0, mean.X, Sig.X, prev){
  x <- rmvnorm(n, rep(mean.X, p), Sig.X)
  beta0(betaX=beta0, X=cbind(1,x), N=rep(1,dim(x)[1]), rhoY=prev, expandX="all")
  
  y=unlist(lapply(1:n, function(ll) rbinom(1,1,g.logit(c(1,x[ll,]) %*% beta0))))
}
#' @export
Data.gen.sim=function(M, n.s, n.t,n.t.other, n.v, p, beta0.s, beta0.t, mean.X.s, mean.X.t, Sig.X.s, Sig.X.t, inter){
  dat.s.list=dat.t.list=vector("list", length = M)
  dat.v=Data.gen2(n.v, p, beta0.t, mean.X.t, Sig.X.t,inter)
  
  for(m in 1:M){
    if(m==1){
      dat.t.list[[m]]=Data.gen2(n.t, p, beta0.t, mean.X.t, Sig.X.t,inter)}else{
        dat.t.list[[m]]=Data.gen2(n.t.other, p, beta0.t, mean.X.t, Sig.X.t,inter) ### this is to make sure other sites have larger sample size, which is consistent with the real practice
      }
    dat.s.list[[m]]=Data.gen2(n.s, p, beta0.s, mean.X.s, Sig.X.s,inter)
    if(sum(dat.t.list[[m]][,1])<2){dat.t.list[[m]][1:2,1]=1}
  }
  names(dat.t.list)=paste0("site", 1:M)
  names(dat.s.list)=paste0("site", 1:M)
  
  return(list(dat.v=dat.v, dat.t.list=dat.t.list, dat.s.list=dat.s.list))
}
#' @export
Data.gen1=function(n, p, prev, mean.X1, mean.X0, Sig.X){
  y=rbinom(n,1,prev)
  x=matrix(NA, ncol=p, nrow=n)
  x[y==1,] = rmvnorm(sum(y==1), rep(mean.X1, p), Sig.X)
  x[y==0,] = rmvnorm(sum(y==0), rep(mean.X0, p), Sig.X)
  data.frame(y, x)
}

#' @export
Coef.gen<- function(s, h, sig.beta, sig.delta, p){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- beta0
  samp0<- sample(1:s, h, replace=F)
  W[samp0] <-W[samp0] + rep(-sig.delta, h)
  return(list(W=W, beta0=beta0))
}
#' @export
Coef.gen.new<- function(s, h, sig.beta, sig.delta, p){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- beta0
  samp0<- sample(1:s, h, replace=F)
  W[samp0] <-W[samp0] + rep(sig.delta, h)
  return(list(W=W, beta0=beta0))
}
#' @export
table.sim.fun=function(mytable, ave.method=mean, M.list=4, n.s.list=5000, n.t.list, prev.t.list, setting.list, 
                       sampling.list, model.list, fed.list, by.cols=c("Model")){
  cols.rm=unique(c(setdiff(c("M","n.s", "n.t", "prev.t","p","s", "sig.beta", "setting", "Sampling", "Model", "Fed", "FRM"),by.cols),
                   "fpr", "cut", "p.pos", "npv"))
  mytable.new=mytable[mytable$M%in%M.list & 
                        mytable$n.s%in%n.s.list & 
                        mytable$n.t%in%n.t.list & 
                        mytable$prev.t%in%prev.t.list &
                        mytable$setting%in%setting.list & 
                        mytable$Sampling%in%sampling.list & 
                        mytable$Model%in%model.list &
                        mytable$Fed%in%fed.list,]          
  mytable.new=data.frame(setDT(mytable.new[,setdiff(colnames(mytable.new), cols.rm)])[, lapply(.SD, ave.method, na.rm=T), by= by.cols])
  tryCatch({mytable.new$Sampling=factor(mytable.new$Sampling,levels=sampling.list)},error=function(e){})
  tryCatch({mytable.new$Model=factor(mytable.new$Model,levels=model.list)},error=function(e){})
  tryCatch({mytable.new$Fed=factor(mytable.new$Fed,levels=fed.list)},error=function(e){})
  mytable.new=mytable.new[ do.call(order, as.list(as.data.frame(mytable.new[,by.cols])) ), ]
  mytable.new
}

#' @export
table.sim.FindBest.fun=function(mytable, ave.method=mean, M.list=4, n.s.list=5000, n.t.list, prev.t.list, setting.list, 
                                sampling.list, model.list, fed.list, by.cols=c("Model")){
  cols.rm=unique(c(setdiff(c("M","n.s", "n.t", "prev.t","p","s", "sig.beta", "setting", "Sampling", "Model", "Fed", "FRM"),by.cols),
                   "fpr", "cut", "p.pos", "npv"))
  mytable.new=mytable[mytable$M%in%M.list & 
                        mytable$n.s%in%n.s.list & 
                        mytable$n.t%in%n.t.list & 
                        mytable$prev.t%in%prev.t.list &
                        mytable$setting%in%setting.list & 
                        mytable$Sampling%in%sampling.list & 
                        mytable$Model%in%model.list &
                        mytable$Fed%in%fed.list,]          
  mytable.new=data.frame(setDT(mytable.new[,setdiff(colnames(mytable.new), cols.rm)])[, lapply(.SD, ave.method, na.rm=T), by= by.cols])
  tryCatch({mytable.new$Sampling=factor(mytable.new$Sampling,levels=sampling.list)},error=function(e){})
  tryCatch({mytable.new$Model=factor(mytable.new$Model,levels=model.list)},error=function(e){})
  tryCatch({mytable.new$Fed=factor(mytable.new$Fed,levels=fed.list)},error=function(e){})
  mytable.new=mytable.new[ do.call(order, as.list(as.data.frame(mytable.new[,by.cols])) ), ]
  mytable.new[which.max(mytable.new$auc.est),]
}

#' @export
table.sim.print.fun=function(mytable, ave.method=mean, M.list=4, n.s.list=5000, n.t.list, prev.t.list, setting.list, 
                             sampling.list, model.list, fed.list, by.cols=c("Model")){
  cols.rm=unique(c(setdiff(c("M","n.s", "n.t", "prev.t","p","s", "sig.beta", "setting", "Sampling", "Model", "Fed", "FRM"),by.cols),
                   "fpr", "cut", "p.pos", "npv"))
  mytable.new=mytable[mytable$M%in%M.list & 
                        mytable$n.s%in%n.s.list & 
                        mytable$n.t%in%n.t.list & 
                        mytable$prev.t%in%prev.t.list &
                        mytable$setting%in%setting.list & 
                        mytable$Sampling%in%sampling.list & 
                        mytable$Model%in%model.list &
                        mytable$Fed%in%fed.list,]          
  mytable.new=data.frame(setDT(mytable.new[,setdiff(colnames(mytable.new), cols.rm)])[, lapply(.SD, ave.method, na.rm=T), by= by.cols])
  tryCatch({mytable.new$Sampling=factor(mytable.new$Sampling,levels=sampling.list)},error=function(e){})
  tryCatch({mytable.new$Model=factor(mytable.new$Model,levels=model.list)},error=function(e){})
  tryCatch({mytable.new$Fed=factor(mytable.new$Fed,levels=fed.list)},error=function(e){})
  mytable.new=mytable.new[ do.call(order, as.list(as.data.frame(mytable.new[,by.cols])) ), ]
  mytable.new
}
#' @export
table.data.print.fun=function(mytable, ave.method=mean, nm.study.list, year.cut.list, sampling.list, model.list, fed.list, frm.list, by.cols=c("Cohort","Model")){
  cols.rm=unique(c(setdiff(c("Cohort", "Year", "Sampling", "Model", "Fed", "FRM"),by.cols),
                   "fpr", "cut", "p.pos", "npv"))
  mytable.new=mytable[mytable$Year%in%year.cut.list & 
                        grepl(paste0(nm.study.list,collapse="|"),mytable$Cohort) & 
                        mytable$Sampling%in%sampling.list & 
                        mytable$Model%in%model.list &
                        mytable$Fed%in%fed.list & 
                        mytable$FRM%in%frm.list,]          
  mytable.new=data.frame(setDT(mytable.new[,setdiff(colnames(mytable.new), cols.rm)])[, lapply(.SD, ave.method, na.rm=T), by= by.cols])
  tryCatch({mytable.new$Sampling=factor(mytable.new$Sampling,levels=sampling.list)},error=function(e){})
  tryCatch({mytable.new$Model=factor(mytable.new$Model,levels=model.list)},error=function(e){})
  tryCatch({mytable.new$Fed=factor(mytable.new$Fed,levels=fed.list)},error=function(e){})
  mytable.new=mytable.new[ do.call(order, as.list(as.data.frame(mytable.new[,by.cols])) ), ]
  mytable.new
}
#' @export
table.fun=function(mytable, ave.method=mean, nm.study.list, year.cut.list, sampling.list, model.list, fed.list, frm.list, by.cols=c("Cohort","Model")){
  cols.rm=unique(c(setdiff(c("Cohort", "Year", "Sampling", "Model", "Fed", "FRM"),by.cols),
                   "fpr", "cut", "p.pos", "npv"))
  mytable.new=mytable[mytable$Year%in%year.cut.list & 
                        grepl(paste0(nm.study.list,collapse="|"),mytable$Cohort) & 
                        mytable$Sampling%in%sampling.list & 
                        mytable$Model%in%model.list &
                        mytable$Fed%in%fed.list & 
                        mytable$FRM%in%frm.list,]          
  mytable.new=data.frame(setDT(mytable.new[,setdiff(colnames(mytable.new), cols.rm)])[, lapply(.SD, ave.method, na.rm=T), by= by.cols])
  tryCatch({mytable.new$Sampling=factor(mytable.new$Sampling,levels=sampling.list)},error=function(e){})
  tryCatch({mytable.new$Model=factor(mytable.new$Model,levels=model.list)},error=function(e){})
  tryCatch({mytable.new$Fed=factor(mytable.new$Fed,levels=fed.list)},error=function(e){})
  mytable.new=mytable.new[ do.call(order, as.list(as.data.frame(mytable.new[,by.cols])) ), ]
  mytable.new
}

#' @export
myrose.fun=function(dat,N.rose=5000, myseed){
  ROSE(Y~., data=dat, seed = myseed, N=max(N.rose, max(table(dat[,"Y"]))*2))$data
}
#' @export
myrose.fed.fun=function(dat,N.rose=5000, myseed){
  ROSE(Y~., data=dat, seed = myseed, N=max(N.rose, max(table(dat[,"Y"]))*2))$data
}
#' @export
myrose.list.fun=function(dat.train, dat.list, nm.study, N.rose=5000, myseed){
  dat.list.new=NULL
  for(ll in setdiff(ls(dat.list), nm.study)){
    dat.tmp=dat.list[[ll]]
    dat.list.new[[ll]]=data.matrix(ROSE(Y~., data=dat.tmp, seed = myseed, N=max(N.rose, max(table(dat.tmp[,"Y"]))*2))$data)
  }
  dat.list.new[[nm.study]]=data.matrix(dat.train)
  dat.list.new
}
#' @export
val.fun.new=function(mybeta, dat.val,nm.out){
  tmp.x=dat.val[,names(mybeta)]
  mys=data.matrix(tmp.x)%*%as.numeric(mybeta)
  list(auc=ROC.Est.FUN(dat.val[,nm.out], mys, yy0=0.5)[1],S=mys)
}
#' @export
val.fun=function(mybeta, dat.val){
  tmp.x=dat.val[,names(mybeta)]
  mys=data.matrix(tmp.x)%*%as.numeric(mybeta)
  list(auc=ROC.Est.FUN(dat.val$stroke, mys, yy0=0.5)[1],S=mys)
}
#' @export
myroc.fun=function(junk, cut.nm="fpr", cut.val=0.1){
  mt=junk[-1]
  mt=data.frame(matrix(mt, ncol=6, byrow=F))
  colnames(mt)=c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
  mt$f1=2*mt$tpr*mt$ppv/(mt$tpr+mt$ppv)
  indx=which.min(abs(mt[,cut.nm]-cut.val))
  mt[indx,]
}
#' @export
Trans.update.fun=function(beta.ini, x.t, y.t, family="binomial"){
  p=dim(x.t)[2]
  n.t=length(y.t)
  lam.const <- 0.0001/sqrt(2*log(p)/n.t)
  myoffset=x.t%*%beta.ini
  delta.fit <- as.numeric(glmnet(x=x.t,y=y.t, offset=myoffset, lambda=lam.const*sqrt(2*log(p)/n.t), family=family)$beta)
  delta.fit<-delta.fit*(abs(delta.fit)>=lam.const*sqrt(2*log(p)/n.t))
  beta.fit <- beta.ini + delta.fit
  beta.fit
}
#' @export
Trans.fun=function(x.s, y.s, x.t, y.t, family="binomial"){
  n.s=length(y.s)
  n.t=length(y.t)
  p=dim(x.t)[2]
  cv.init<-cv.glmnet(x=x.s, y=y.s, nfolds=8, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/n.s), family=family)
  lam.const <- cv.init$lambda.min/sqrt(2*log(p)/n.s)
  w.fit <- as.numeric(glmnet(x.s, y.s, lambda=lam.const*sqrt(2*log(p)/n.s), family=family)$beta)
  w.fit<-w.fit*(abs(w.fit)>=lam.const*sqrt(2*log(p)/n.s))
  myoffset=x.t%*%w.fit
  delta.fit <- as.numeric(glmnet(x=x.t,y=y.t, offset=myoffset, lambda=lam.const*sqrt(2*log(p)/n.t), family=family)$beta)
  delta.fit<-delta.fit*(abs(delta.fit)>=lam.const*sqrt(2*log(p)/n.t))
  beta.fit <- w.fit + delta.fit
  beta.fit
}
#' @export
iteration.fun2=function (dat.list, bini, kk.list) {
  K = length(dat.list)
  Uini.list = sapply(kk.list, function(kk) {
    U.fun(bini, dat = dat.list[[kk]])
  })
  Aini.list = lapply(1:K, function(kk) {
    A.fun(bini, dat = dat.list[[kk]])
  })
  Ahat.ini = Reduce("+", Aini.list)/K
  bhat.list = -ginv(Ahat.ini) %*% Uini.list + bini
  list(b.k = bhat.list, Ahat = Ahat.ini)
}
#' @export
DCOS.FUN=function(dat.list, niter, ridge=F, lambda=NULL){
  K=length(dat.list)
  lambda.grid=10^seq(-100,3,0.1)
  N=sum(unlist(lapply(dat.list, function(xx) dim(xx)[1])))
  ### Round 1
  dat.list.00=dat.list[[1]]
  y1=dat.list.00[,1]
  x1=dat.list.00[,-1]
  
  ##initial
  if(ridge==F){
    bini = glm(y1~x1,family="binomial")$coef}
  if(ridge==T){
    bini = as.vector(coef(glmnet(x1,y1,alpha=0,standardize=F,lambda=lambda,family="binomial")))
  }
  ##update
  update=iteration.fun2(dat.list=dat.list,bini=bini,kk.list=1:K)
  bnew.update=apply(update$b.k,1,mean)
  
  for(ii in 1:(niter-1)){
    bnew=bnew.update
    update.DCOS = iteration.fun2(dat.list=dat.list,bini=bnew,kk.list=1:K)
    bnew.update = apply(update.DCOS$b.k,1,mean) 
  }
  
  Ahalf.DCOS= svd(-update.DCOS$Ahat); Ahalf.DCOS = Ahalf.DCOS$u%*%diag(sqrt(Ahalf.DCOS$d))%*%t(Ahalf.DCOS$v)
  betahat = Est.ALASSO.Approx.GLMNET(ynew=Ahalf.DCOS%*%bnew.update,xnew=Ahalf.DCOS,bini=bnew.update,N.adj=N,lambda.grid)$bhat.BIC
  betahat
}
#' @export
DCOS.Trans.FUN=function(bini, dat.list){
  K=length(dat.list)
  lambda.grid=10^seq(-100,3,0.1)
  N=sum(unlist(lapply(dat.list, function(xx) dim(xx)[1])))
  
  ##update
  update.DCOS=iteration.fun2(dat.list=dat.list,bini=bini,kk.list=1:K)
  bnew.update=apply(update.DCOS$b.k,1,mean)
  
  Ahalf.DCOS= svd(-update.DCOS$Ahat); Ahalf.DCOS = Ahalf.DCOS$u%*%diag(sqrt(Ahalf.DCOS$d))%*%t(Ahalf.DCOS$v)
  betahat = Est.ALASSO.Approx.GLMNET(ynew=Ahalf.DCOS%*%bnew.update,xnew=Ahalf.DCOS,bini=bnew.update,N.adj=N,lambda.grid)$bhat.BIC
  betahat
}

#' @export
Trans_DCOS.FUN=function(dat.list, niter, ridge=F, lambda=NULL){
  K=length(dat.list)
  lambda.grid=10^seq(-100,3,0.1)
  N=sum(unlist(lapply(dat.list, function(xx) dim(xx)[1])))
  ### Round 1
  dat.list.00=dat.list[[1]]
  y1=dat.list.00[,1]
  x1=dat.list.00[,-1]
  
  ##initial
  if(ridge==F){
    bini = glm(y1~x1,family="binomial")$coef}
  if(ridge==T){
    bini = as.vector(coef(glmnet(x1,y1,alpha=0,standardize=F,lambda=lambda,family="binomial")))
  }
  ##update
  update=iteration.fun2(dat.list=dat.list,bini=bini,kk.list=1:K)
  bnew.update=apply(update$b.k,1,mean)
  
  for(ii in 1:(niter-1)){
    bnew=bnew.update
    update.DCOS = iteration.fun2(dat.list=dat.list,bini=bnew,kk.list=1:K)
    bnew.update = apply(update.DCOS$b.k,1,mean) 
  }
  
  Ahalf.DCOS= svd(-update.DCOS$Ahat); Ahalf.DCOS = Ahalf.DCOS$u%*%diag(sqrt(Ahalf.DCOS$d))%*%t(Ahalf.DCOS$v)
  betahat = Est.ALASSO.Approx.GLMNET(ynew=Ahalf.DCOS%*%bnew.update,xnew=Ahalf.DCOS,bini=bnew.update,N.adj=N,lambda.grid)$bhat.BIC
  betahat
}

#' @export
density_fit_fun=function(y.s, x.s, y.t, x.t, is.pca, fit.type="all"){
  y.dr.fit=c(rep(0,length(y.s)), rep(1,length(y.t)))
  x.dr.fit=rbind(x.s, x.t)
  if(fit.type=="source"){y.fit=y.s; x.fit=x.s}
  if(fit.type=="target"){y.fit=y.t; x.fit=x.t}
  if(fit.type=="all"){y.fit=c(y.s, y.t); x.fit=rbind(x.s,x.t)}
  
  if(is.pca){
    fit.pca=prcomp(x.dr.fit)
    junk=glm(y.dr.fit~fit.pca$x[,1:5], family="binomial")
    pp1=as.numeric(g.logit(cbind(1,predict(fit.pca, x.fit)[,1:5])%*%coef(junk)))}else{
      junk=lasso.fun2(x.dr.fit, y.dr.fit,family="binomial")
      pp1=as.numeric(g.logit(cbind(1,x.fit)%*%junk))}
  fit.ds=Est.ALASSO.GLMNET(data=cbind(y.fit,x.fit),fam0="binomial", Wi=pp1/(1-pp1))$bhat.BIC
  fit.ds
}
#' @export
density_fed_fit_fun=function(y.s.list, x.s.list, y.t.list, x.t.list, y.s, x.s, y.t, x.t, is.pca){
  y.dr.fit.list=lapply(1:length(y.s.list),function(ll) c(rep(0,length(y.s.list[[ll]])), rep(1,length(y.t.list[[ll]]))))
  x.dr.fit.list=lapply(1:length(x.s.list),function(ll) rbind(x.s.list[[ll]], x.t.list[[ll]]))
  x.fit=rbind(x.s,x.t)
  y.fit=c(y.s,y.t)
  
  pp1.list=list()
  for(ll in 1:length(y.dr.fit.list)){
    if(is.pca){
      fit.pca=prcomp(x.dr.fit.list[[ll]])
      junk=glm(y.dr.fit.list[[ll]]~fit.pca$x[,1:5], family="binomial")
      pp1.list[[ll]]=g.logit(cbind(1,predict(fit.pca, x.fit)[,1:5])%*%coef(junk))
    }else{
      junk=lasso.fun2(x.dr.fit.list[[ll]], y.dr.fit.list[[ll]],family="binomial")
      pp1.list[[ll]]=g.logit(cbind(1,x.fit)%*%junk)}
  }
  pp1=rowMeans(do.call(cbind,pp1.list))
  fit.ds=Est.ALASSO.GLMNET(data=cbind(y.fit,x.fit),fam0="binomial", Wi=pp1/(1-pp1))$bhat.BIC
  fit.ds
}
#' @export
cohort.sel.fun=function(dat.split, nm.cohort.t, nm.cov,rr=2){ 
  x.t=dat.split[dat.split$cohort==nm.cohort.t ,nm.cov]
  y.t=dat.split[dat.split$cohort==nm.cohort.t,"stroke"]
  
  cov.t=cov(x.t)
  cov.diff=r.case=NULL
  for(nm.cohort in unique(dat.split$cohort)){
    x.s=dat.split[dat.split$cohort==nm.cohort,nm.cov]
    y.s=dat.split[dat.split$cohort==nm.cohort,"stroke"]
    cov.s=cov(x.s)
    cov.diff=c(cov.diff,sqrt(sum((cov.t-cov.s)^2)))
    r.case=c(r.case, sum(y.s)/sum(y.t))
  }
  names(cov.diff)=names(r.case)=unique(dat.split$cohort)
  ## if the effect size of target sample is small, sample size will be more important
  z=rr*length(nm.cov)/sum(y.t)*r.case/((max(r.case)-min(r.case)))-cov.diff/((max(cov.diff)-min(cov.diff)))
  z=sort(z, decreasing=T)
  beta.sel=names(z)
  beta.sel=beta.sel[grepl("ARIC|MESA|REGARDS|OFFSPRING",beta.sel)==1]
  beta.sel=beta.sel[1:5]
  beta.sel
}
#' @export
evaluate.fun=function(y,s, cut.metric="fpr", metric.cut=0.5){
  junk=ROC.Est.FUN(y,s, yy0=0.5,fpr0=seq(0,1,0.01))
  junk.mtx=data.frame(matrix(junk[-1], ncol=6))
  colnames(junk.mtx)=c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
  junk.mtx=junk.mtx[-1,]
  junk.mtx$f1=2*junk.mtx$tpr*junk.mtx$ppv/(junk.mtx$tpr+junk.mtx$ppv)
  if(cut.metric!="f1"){
    myroc=junk.mtx[which.min(abs(junk.mtx[,cut.metric]-metric.cut)),]
  }else{
    myroc=junk.mtx[which.max(junk.mtx[,cut.metric]),]  
  }
  myroc=as.numeric(c(junk[1], myroc))
  pr <- pr.curve(scores.class0 = s, weights.class0 = y)
  myroc=data.frame(cut.metric, metric.cut,pr$auc.integral,t(myroc), max(junk.mtx[,"f1"],na.rm=T))
  colnames(myroc)=c("cut.metric", "metric.cut","auprc", "aucroc", colnames(junk.mtx), "f1.max")
  myroc
}
#' @export
lasso.fun=function(x,y, family, offset=NULL, weights=NULL, penalty.factor=NULL, lambda=NULL,alpha=1){
  if(is.null(weights)){weights=rep(1,dim(x)[1])}
  if(is.null(penalty.factor)){penalty.factor=rep(1,dim(x)[1])}
  if(is.null(lambda)){
    cvfit = cv.glmnet(x, y, family=family, offset=offset, weights=weights, alpha=alpha,penalty.factor=penalty.factor)
    coef.1se = coef(cvfit, s = "lambda.min")
  }else{
    fit = glmnet(x, y, family=family, offset=offset, weights=weights,lambda=lambda, alpha=alpha,penalty.factor=penalty.factor)
    coef.1se = coef(fit)
    
  }
  as.numeric(coef.1se)
}
#' @export
weight.fun=function(x.source, y.source, x.target, y.target){
  #y.tmp=c(rep(0,length(y.target)), rep(1, length(y.target)))
  y.tmp=c(rep(1,length(y.source)), rep(0, length(y.target)))
  x.tmp=rbind(cbind(x.source), cbind(x.target))
  #x.tmp=rbind(cbind(y.source,x.source)[1:length(y.target),],cbind(y.target,x.target))
  junk=lasso.fun(x=x.tmp, y=y.tmp, family="binomial")
  p.all=g.logit(cbind(1,x.tmp)%*%junk)
  p.source=p.all[1:length(y.source)]
  p.target=p.all[-c(1:length(y.source))]
  auc=ROC.Est.FUN(y.tmp,p.all,yy0=0.5)[1]
  list(beta=junk,auc=auc,p.all=p.all, p.source=p.source, p.target=p.target)
}
#' @export
data.gen.fun=function(n.source, n.target, p, setting,
                      setting1.delta.mean=1,
                      setting1.delta.sd=1,
                      setting2.h=4, 
                      setting2.delta=0.3){
  if(setting==1){
    ### three race related variables: x1, x2, x3
    x.source=matrix(rnorm(n.source*p, 0, 1),ncol=p)
    x.target=matrix(rnorm(n.target*p, 0, 1),ncol=p)
    x.source[,c(1:10)]=matrix(rnorm(n.source*10, 1, 1),ncol=10)
    x.target[,c(1:10)]=matrix(rnorm(n.target*10, 1+setting1.delta.mean, 1+setting1.delta.sd),ncol=10)
    
    b.target=b.source=c(0.4,rep(0.3,16),rep(0,p-16))
  }
  
  if(setting==2){
    x.source=matrix(rnorm(n.source*p, 0, 1),ncol=p)
    x.target=matrix(rnorm(n.target*p, 0, 1),ncol=p)
    b.target=c(0.4,rep(0.3,16),rep(0,p-16))
    #set.seed(1234)
    #index=sample(1:p, setting2.h, replace=F)
    index=11:(setting2.h+10)
    
    ind=rep(0,p)
    ind[index]=1
    b.source=b.target[-1]+setting2.delta*ind
    b.source=c(0.4,b.source)
  }
  
  if(setting==3){
    x.source=matrix(rnorm(n.source*p, 0, 1),ncol=p)
    x.target=matrix(rnorm(n.target*p, 0, 1),ncol=p)
    x.source[,c(1:3)]=matrix(rnorm(n.source*3, 1, 1),ncol=3)
    x.target[,c(1:3)]=matrix(rnorm(n.target*3, 1+setting1.delta.mean, 1+setting1.delta.sd),ncol=3)
    
    b.target=c(0.4,rep(0.3,16),rep(0,p-16))
    #set.seed(1234)
    #index=sample(1:p, setting2.h, replace=F)
    index=11:(setting2.h+10)
    
    ind=rep(0,p)
    ind[index]=1
    b.source=b.target[-1]+setting2.delta*ind
    b.source=c(0.4,b.source)
  }
  y.source=unlist(lapply(g.logit(cbind(1,x.source)%*%b.source), function(pp) rbinom(1,1,pp)))
  y.target=unlist(lapply(g.logit(cbind(1,x.target)%*%b.target), function(pp) rbinom(1,1,pp)))
  return(list(b.source=b.source, b.target=b.target,y.source=y.source, y.target=y.target, x.source=x.source, x.target=x.target))
}

#' @export
data.gen.fun2=function(n.s, n.t, p, a, sigma2.ee){
  x.s=matrix(rnorm(n.s*p, 0, 1),ncol=p)
  x.t=matrix(rnorm(n.t*p, 0, 1),ncol=p)
  theta.s=c(1,rnorm(p,0,1))
  theta.w=c(1,rnorm(p,0,1))
  theta.t=a*theta.s+(1-a)*theta.w
  #ee.t=rnorm(n.t, 0, sigma2.ee)
  y.s=unlist(lapply(g.logit(cbind(1,x.s)%*%theta.s), function(pp) rbinom(1,1,pp)))
  y.t=unlist(lapply(g.logit(cbind(1,x.t)%*%theta.t), function(pp) rbinom(1,1,pp)))
  return(list(y.s=y.s, x.s=x.s, y.t=y.t, x.t=x.t, theta.s=theta.s, theta.w=theta.w, theta.t=theta.t))
}
#' @export
lasso.fun2=function(x,y, family, offset=NULL, weights=NULL, lambda=NULL,alpha=1, method="lasso"){
  if(is.null(weights)){weights=rep(1,dim(x)[1])}
  if(method=="lasso"){
    if(is.null(lambda)){
      cvfit = cv.glmnet(x, y, family=family, offset=offset, weights=weights, alpha=alpha)
      coef.1se = coef(cvfit, s = "lambda.min")
    }else{
      fit = glmnet(x, y, family=family, offset=offset, weights=weights,lambda=lambda, alpha=alpha)
      coef.1se = coef(fit)
    }
    res=as.numeric(coef.1se)}
  
  if(method=="glm"){
    fit = glm(as.numeric(y)~x, family=family)
    res=coef(fit)
  }
  res
}
#' @export
Est.ALASSO.GLMNET=function (data, BIC.factor = 0.1, fam0 = "binomial", w.b = NULL, Wi=NULL, lambda.grid = 10^seq(-4,0,0.01)) 
{  if(is.null(Wi)==1){Wi=rep(1,dim(data)[1])}
  if (fam0 != "Cox") {
    data = as.matrix(data)
    y = data[, 1]
    x = data[, -1, drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
    if (is.null(w.b)) {
      bini = as.vector(coef(glmnet(x,y, family=fam0, weights=Wi,alpha=0,standardize=F,lambda=0.1)))
      w.b = 1/abs(bini[-1])
      gc()
    }
    tmpfit = glmnet(x = x, y = y, family = fam0, penalty.factor = w.b, weights=Wi,
                    alpha = 1, lambda = lambda.grid, intercept = T)
    gc()
    N.adj = dim(x)[1]
  }
  else {
    data = as.matrix(data)
    y = cbind(time = data[, 1], status = data[, 2])
    x = data[, -c(1, 2), drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
    if (is.null(w.b)) {
      bini = as.vector(coef(glmnet(x, y, family = "cox",  weights=Wi,
                                   alpha = 0, lambda = 0.1)))
      w.b = 1/abs(bini)
      gc()
    }
    tmpfit = glmnet(x, y, family = "cox", penalty.factor = w.b, weights=Wi,
                    alpha = 1 , lambda = lambda.grid)
    gc()
    N.adj = sum(y[, 2])
  }
  dev = deviance(tmpfit)
  BIC.lam = dev + min(N.adj^0.1, log(N.adj)) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.modBIC = c(tmpfit$a0[m.opt], tmpfit$beta[, m.opt])
  lamhat.modBIC = tmpfit$lambda[m.opt]
  BIC.lam = dev + log(N.adj) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.BIC = c(tmpfit$a0[m.opt], tmpfit$beta[, m.opt])
  lamhat.BIC = tmpfit$lambda[m.opt]
  
  
  tLL <- tmpfit$nulldev - deviance(tmpfit)
  k <- tmpfit$df
  n <- tmpfit$nobs
  AIC.lam <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  
  m.opt=which.min(AIC.lam)
  bhat.AIC=c(tmpfit$a0[m.opt], tmpfit$beta[,m.opt])
  lamhat.AIC = tmpfit$lambda[m.opt]
  
  return(list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC, bhat.AIC=bhat.AIC,
              lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC, lambda.AIC=lamhat.AIC))
}
#' @export
Est.ALASSO.GLMNET.offset=function (data, BIC.factor = 0.1, fam0 = "binomial", w.b = NULL, Wi=NULL, lambda.grid = 10^seq(-4,0,0.01), myoffset) 
{  if(is.null(Wi)==1){Wi=rep(1,dim(data)[1])}
  if (fam0 != "Cox") {
    data = as.matrix(data)
    y = data[, 1]
    x = data[, -1, drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
    if (is.null(w.b)) {
      bini = as.vector(coef(glmnet(x,y, family=fam0, offset=myoffset,weights=Wi,alpha=0,standardize=F,lambda=0.1)))
      w.b = 1/abs(bini[-1])
      gc()
    }
    tmpfit = glmnet(x = x, y = y, family = fam0, offset=myoffset,penalty.factor = w.b, weights=Wi,
                    alpha = 1, lambda = lambda.grid, intercept = T)
    gc()
    N.adj = dim(x)[1]
  }
  else {
    data = as.matrix(data)
    y = cbind(time = data[, 1], status = data[, 2])
    x = data[, -c(1, 2), drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
    if (is.null(w.b)) {
      bini = as.vector(coef(glmnet(x, y, family = "cox", offset=myoffset, weights=Wi,
                                   alpha = 0, lambda = 0.1)))
      w.b = 1/abs(bini)
      gc()
    }
    tmpfit = glmnet(x, y, family = "cox", offset=myoffset,penalty.factor = w.b, weights=Wi,
                    alpha = 1 , lambda = lambda.grid)
    gc()
    N.adj = sum(y[, 2])
  }
  dev = deviance(tmpfit)
  BIC.lam = dev + min(N.adj^0.1, log(N.adj)) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.modBIC = c(tmpfit$a0[m.opt], tmpfit$beta[, m.opt])
  lamhat.modBIC = tmpfit$lambda[m.opt]
  BIC.lam = dev + log(N.adj) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.BIC = c(tmpfit$a0[m.opt], tmpfit$beta[, m.opt])
  lamhat.BIC = tmpfit$lambda[m.opt]
  
  
  tLL <- tmpfit$nulldev - deviance(tmpfit)
  k <- tmpfit$df
  n <- tmpfit$nobs
  AIC.lam <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  
  m.opt=which.min(AIC.lam)
  bhat.AIC=c(tmpfit$a0[m.opt], tmpfit$beta[,m.opt])
  lamhat.AIC = tmpfit$lambda[m.opt]
  
  return(list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC, bhat.AIC=bhat.AIC,
              lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC, lambda.AIC=lamhat.AIC))
}
