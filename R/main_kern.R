
main_kern=function(dat.train.s, dat.train.t, 
                   dat.s.list, dat.t.list,
                   nm.study,
                   nm.cov,
                   nm.out,
                   sampling,
                   model,
                   fed,
                   nm.study.rm.fed=NULL,
                   N.rose=5000, myseed=1){

dat.train.s.new=dat.train.s[,c(nm.out, nm.cov)]
dat.train.t.new=dat.train.t[,c(nm.out, nm.cov)]
colnames(dat.train.s.new)[colnames(dat.train.s.new)==nm.out]="Y"
colnames(dat.train.t.new)[colnames(dat.train.t.new)==nm.out]="Y"

dat.s.list.new=lapply(dat.s.list, function(xx){xx.new=xx[,c(nm.out,nm.cov)];colnames(xx.new)[colnames(xx.new)==nm.out]="Y"; return(xx.new)})
dat.t.list.new=lapply(dat.t.list, function(xx){xx.new=xx[,c(nm.out,nm.cov)];colnames(xx.new)[colnames(xx.new)==nm.out]="Y"; return(xx.new)})
nm.study.rm.s=names(dat.s.list.new)[which(unlist(lapply(dat.s.list.new, function(xx) dim(xx)[1]))==0)]
nm.study.rm.t=names(dat.t.list.new)[which(unlist(lapply(dat.t.list.new, function(xx) dim(xx)[1]))==0)]

dat.s.list.new=dat.s.list.new[setdiff(names(dat.s.list.new), nm.study.rm.s)]
dat.t.list.new=dat.t.list.new[setdiff(names(dat.t.list.new), nm.study.rm.t)]
#print("tag: sampling")
### sampling
if(sampling=="Rose"){
dat.train.s.new <- myrose.fun(dat.train.s.new, myseed=myseed)
dat.train.t.new <- myrose.fun(dat.train.t.new, myseed=myseed)

dat.s.list.new <- myrose.list.fun(dat.train=dat.train.s.new, dat.list=dat.s.list.new, nm.study=nm.study, myseed=myseed)
dat.t.list.new <- myrose.list.fun(dat.train=dat.train.t.new, dat.list=dat.t.list.new, nm.study=nm.study, myseed=myseed)
}

### change data.frame to data.matrix 
dat.train.s.new=data.matrix(dat.train.s.new)
dat.train.t.new=data.matrix(dat.train.t.new)

dat.list.new=list()
for(ll in unique(c(ls(dat.s.list.new), ls(dat.t.list.new)))){
dat.list.new[[ll]]=data.matrix(rbind(dat.s.list.new[[ll]], dat.t.list.new[[ll]]))
if(is.null(dat.s.list.new[[ll]])!=1){dat.s.list.new[[ll]]=data.matrix(dat.s.list.new[[ll]][,setdiff(colnames(dat.s.list.new[[ll]]), "black")])}
if(is.null(dat.t.list.new[[ll]])!=1){dat.t.list.new[[ll]]=data.matrix(dat.t.list.new[[ll]][,setdiff(colnames(dat.t.list.new[[ll]]), "black")])}
}

### remove pre-specified study for federated learning
dat.s.list.new=dat.s.list.new[setdiff(names(dat.s.list.new), nm.study.rm.fed)]
dat.t.list.new=dat.t.list.new[setdiff(names(dat.t.list.new), nm.study.rm.fed)]
dat.list.new=dat.list.new[setdiff(names(dat.list.new), nm.study.rm.fed)]

#print("tag: data")

### prepare data for modeling
if(model=="Source"){
  dat.train.new=dat.train.s.new[,setdiff(colnames(dat.train.s.new), "black")]
  dat.train.list.new=dat.s.list.new}
if(model=="Target"){
  dat.train.new=dat.train.t.new[,setdiff(colnames(dat.train.t.new), "black")]
  dat.train.list.new=dat.t.list.new}
if(model=="Pool"){
  dat.train.new=rbind(dat.train.s.new,dat.train.t.new)
  dat.train.list.new=dat.list.new}

#print("tag: modeling")
### non-transfer modeling
if(model!="Translasso"){
  if(fed=="nFed"){
    fit=Est.ALASSO.GLMNET(dat.train.new, fam0="binomial")$bhat.BIC
    names(fit)=c("intercept",colnames(dat.train.new)[-1])}
  if(fed=="Dac"){
    fit=DCOS.FUN(dat.train.list.new, niter=2, ridge=T, lambda=0.0001)
    names(fit)=c("intercept",colnames(dat.train.list.new[[1]])[-1])}
  if(fed=="Ave"){
    fit=lapply(dat.train.list.new, function(xx) Est.ALASSO.GLMNET(xx, fam0="binomial")$bhat.BIC)
    fit=do.call(cbind, fit)
    fit.rm=colnames(fit)[which(colSums(fit[-1,])==0)]
    fit=fit[,setdiff(colnames(fit), fit.rm)]
    fit=apply(fit,2,function(x) x/max(abs(x)))
    fit=as.numeric(rowMeans(fit))
    names(fit)=c("intercept",colnames(dat.train.list.new[[1]])[-1])}
  if(fed=="Meta"){
    fit=lapply(dat.train.list.new, function(xx) Est.ALASSO.GLMNET(xx, fam0="binomial")$bhat.BIC)
    fit=do.call(cbind,fit)
    fit.rm=colnames(fit)[which(colSums(fit[-1,])==0)]
    fit=fit[,setdiff(colnames(fit), fit.rm)]
    fit=apply(fit,2,function(x) x/max(abs(x)))
    meta.wt=unlist(lapply(dat.train.list.new,function(xx) sum(xx[,"Y"])))
    meta.wt=meta.wt[setdiff(names(meta.wt), fit.rm)]
    fit=as.numeric(fit%*%meta.wt/sum(meta.wt))
    names(fit)=c("intercept",colnames(dat.train.list.new[[1]])[-1])}
  if(fed=="Refit"){
    bini=lapply(dat.train.list.new, function(xx) Est.ALASSO.GLMNET(xx, fam0="binomial")$bhat.BIC)
    bini.rm=do.call(cbind,bini)
    bini.rm=colnames(bini.rm)[which(colSums(bini.rm[-1,])==0)]
    bini=bini[setdiff(names(bini), bini.rm)]
    dat.refit=dat.train.t.new[,c("Y",names(bini[[1]])[-1])]
    fit.score=do.call(cbind,lapply(bini, function(bb) dat.refit[,-1]%*%bb[-1]))
    refit=coef(glm(dat.refit[,"Y"]~fit.score, family="binomial"))
    refit[is.na(refit)]=0
    fit=c(refit[1],rowSums(do.call(cbind,lapply(1:length(bini), function(ll) refit[ll+1]*bini[[ll]][-1]))))
    names(fit)=c("intercept", colnames(dat.train.list.new[[1]])[-1])}
}

### transfer modeling
if(model=="Translasso"){
  if(fed=="nFed"){
    fit=Trans.fun(x.s=data.matrix(dat.train.s.new[,setdiff(colnames(dat.train.s.new),c("Y", "black"))]), 
                  y.s=dat.train.s.new[,"Y"], 
                  x.t=data.matrix(dat.train.t.new[,setdiff(colnames(dat.train.t.new),c("Y", "black"))]), 
                  y.t=dat.train.t.new[,"Y"], family="binomial")
    names(fit)=setdiff(colnames(dat.train.s.new),c("Y", "black"))}
  if(fed=="Dac"){
    nm.study.intersect=intersect(ls(dat.s.list.new), ls(dat.t.list.new))
    bini=lapply(nm.study.intersect, function(xx){
      Trans.fun(x.s=data.matrix(dat.s.list.new[[xx]][,setdiff(colnames(dat.s.list.new[[xx]]),c("Y", "black"))]), 
                y.s=dat.s.list.new[[xx]][,"Y"], 
                x.t=data.matrix(dat.t.list.new[[xx]][,setdiff(colnames(dat.t.list.new[[xx]]),c("Y", "black"))]), 
                y.t=dat.t.list.new[[xx]][,"Y"], family="binomial")
    })
    bini=do.call(cbind,bini)
    bini=c(0,as.numeric(rowMeans(bini)))
    fit=DCOS.Trans.FUN(bini, dat.t.list.new)
    names(fit)=c("intercept",colnames(dat.t.list.new[[1]])[-1])}
  if(fed=="Ave"){
    nm.study.intersect=intersect(ls(dat.s.list.new), ls(dat.t.list.new))
    fit=lapply(nm.study.intersect, function(xx){
      Trans.fun(x.s=data.matrix(dat.s.list.new[[xx]][,setdiff(colnames(dat.s.list.new[[xx]]),c("Y", "black"))]), 
                  y.s=dat.s.list.new[[xx]][,"Y"], 
                  x.t=data.matrix(dat.t.list.new[[xx]][,setdiff(colnames(dat.t.list.new[[xx]]),c("Y", "black"))]), 
                  y.t=dat.t.list.new[[xx]][,"Y"], family="binomial")
    })
    fit=do.call(cbind, fit)
    colnames(fit)=nm.study.intersect
    fit.rm=colnames(fit)[which(colSums(fit[-1,])==0)]
    fit=fit[,setdiff(colnames(fit), fit.rm)]
    if(is.null(dim(fit))!=1){
    fit=apply(fit,2,function(x) x/max(abs(x)))
    fit=as.numeric(rowMeans(fit))}
    names(fit)=setdiff(colnames(dat.s.list.new[[1]]),c("Y", "black"))}
  if(fed=="Meta"){
    nm.study.intersect=intersect(ls(dat.s.list.new), ls(dat.t.list.new))
    fit=lapply(nm.study.intersect, function(xx){
      Trans.fun(x.s=data.matrix(dat.s.list.new[[xx]][,setdiff(colnames(dat.s.list.new[[xx]]),c("Y", "black"))]), 
                y.s=dat.s.list.new[[xx]][,"Y"], 
                x.t=data.matrix(dat.t.list.new[[xx]][,setdiff(colnames(dat.t.list.new[[xx]]),c("Y", "black"))]), 
                y.t=dat.t.list.new[[xx]][,"Y"], family="binomial")
    })
    names(fit)=nm.study.intersect
    fit=do.call(cbind,fit)
    fit.rm=colnames(fit)[which(colSums(fit[-1,])==0)]
    fit=fit[,setdiff(colnames(fit), fit.rm)]
    if(is.null(dim(fit))!=1){
    fit=apply(fit,2,function(x) x/max(abs(x)))
    meta.wt=unlist(lapply(setdiff(nm.study.intersect,fit.rm),function(ll) sum(dat.t.list[[ll]][,nm.out])))
    fit=as.numeric(fit%*%meta.wt/sum(meta.wt))
    }
    names(fit)=setdiff(colnames(dat.s.list.new[[1]]),c("Y", "black"))}
  if(fed=="Refit"){
    nm.study.intersect=intersect(ls(dat.s.list.new), ls(dat.t.list.new))
    bini=lapply(nm.study.intersect, function(xx){
      Trans.fun(x.s=data.matrix(dat.s.list.new[[xx]][,setdiff(colnames(dat.s.list.new[[xx]]),c("Y", "black"))]), 
                y.s=dat.s.list.new[[xx]][,"Y"], 
                x.t=data.matrix(dat.t.list.new[[xx]][,setdiff(colnames(dat.t.list.new[[xx]]),c("Y", "black"))]), 
                y.t=dat.t.list.new[[xx]][,"Y"], family="binomial")
    })
    names(bini)=nm.study.intersect
    bini.rm=do.call(cbind,bini)
    bini.rm=colnames(bini.rm)[which(colSums(bini.rm[-1,])==0)]
    bini=bini[setdiff(names(bini), bini.rm)]
    dat.refit=dat.train.t.new[,setdiff(colnames(dat.train.t.new), "black")]
    fit.score=do.call(cbind,lapply(bini, function(bb) dat.refit[,-1]%*%bb))
    refit=coef(glm(dat.refit[,"Y"]~fit.score, family="binomial"))
    refit[is.na(refit)]=0
    fit=c(refit[1],rowSums(do.call(cbind,lapply(1:length(bini), function(ll) refit[ll+1]*bini[[ll]]))))
    names(fit)=c("intercept",setdiff(colnames(dat.s.list.new[[1]]),c("Y", "black")))}
}

betahat=rep(0, length(nm.cov)+1)
names(betahat)=c("intercept",nm.cov)
betahat[names(fit)]=fit
betahat
}
