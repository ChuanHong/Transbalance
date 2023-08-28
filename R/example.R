
#source("code/Est.Survfit.COXNET.R")
source("code/Est.ALASSO.GLMNET.R")
source("code/library_v2.R")
source("code/function.R")
source("code/main_kern.R")

### parameter settings
Xmean.s=Xmean.t+sig.Xmean.delta
Xsd.s=Xsd.t+sig.Xsd.delta

Sig.X.t <- diag(Xsd.t, p)
Sig.X.s <- diag(Xsd.s, p)


set.seed(1234)
coef.all <-Coef.gen(s, h = h, sig.beta = sig.beta, sig.delta = sig.beta.delta, p = p) 
beta0.t=c(int.beta, coef.all$beta0)
beta0.s=c(int.beta, coef.all$W)

nsim=500
res=NULL
for(isim in 1:nsim){
  print(isim)
  myseed=1234+isim
  set.seed(myseed)
  tryCatch({
### other sites need to increase the sample size of t
dat.junk=Data.gen.sim(M=M, n.s=n.s, n.t=n.t,n.t.other, n.v=n.v, p=p, beta0.s=beta0.s, beta0.t=beta0.t, mean.X.s=Xmean.s, mean.X.t=Xmean.t, Sig.X.s=Sig.X.s, Sig.X.t=Sig.X.t, inter)
dat.val=dat.junk$dat.v
dat.s.list=dat.junk$dat.s.list
dat.t.list=dat.junk$dat.t.list
dat.train.s=dat.s.list[["site1"]]
dat.train.t=dat.t.list[["site1"]]
auc.00=ROC.Est.FUN(dat.val$y, data.matrix(dat.val[,-1])%*%beta0.t[-1], yy0=0.5)[1]

beta.tmp=main_kern(dat.train.s, dat.train.t, 
                   dat.s.list, dat.t.list,
                   nm.study="site1",
                   nm.cov=colnames(dat.train.s)[-1],
                   nm.out="y",
                   sampling=sampling,
                   model=model,
                   fed=fed,
                   nm.study.rm.fed="None", myseed=myseed)

junk=val.fun.new(beta.tmp[-1], dat.val,nm.out="y")
if(length(unique(junk$S))==1){junk$S[1]=1}
tmp.table=data.frame(
  M=M,
  n.s=n.s,
  n.v=n.v, 
  n.t=n.t,
  n.t.other=n.t.other,
  prev.s=mean(dat.train.s[,1]),
  prev.t=prev.t, 
  neffect.t=n.t*prev.t,
  n.t=n.t,
  p=p, 
  s=s, 
  inter=inter,
  int.beta=int.beta,
  sig.beta=sig.beta,
  auc=auc,
  auc.00=auc.00,
  Xmean.t=Xmean.t,
  Xsd.t=Xsd.t,
  setting=setting, 
  h=h, 
  sig.beta.delta=sig.beta.delta,
  sig.Xmean.delta=sig.Xmean.delta,
  sig.Xsd.delta=sig.Xsd.delta,
  sampling=sampling, 
  model=model, 
  fed=fed, isim=isim,
  (rbind(
  evaluate.fun(y=dat.val$y, s=junk$S, cut.metric="fpr", metric.cut=0.1),
  evaluate.fun(y=dat.val$y, s=junk$S, cut.metric="fpr", metric.cut=0.15),
  evaluate.fun(y=dat.val$y, s=junk$S, cut.metric="ppv", metric.cut=0.9),
  evaluate.fun(y=dat.val$y, s=junk$S, cut.metric="f1", metric.cut=0.9)
                     )))
res=rbind(res,tmp.table)
}, error=function(e) NA)
if(nsim>=100){
if(isim%in%seq(100,nsim,100)){save(res, file=paste0("simulation/result/sim.", 
                        "M", M,
                        ".n.s", n.s,
                        ".n.v", n.v,
                        ".n.t", n.t,
                        ".n.t.other", n.t.other,
                        ".prev.t", prev.t,
                        ".p",p,
                        ".s", s,
                        ".inter", inter,
                        ".int.beta", int.beta,
                        ".sig.beta.", sig.beta,
                        ".Xmean.t",Xmean.t,
                        ".Xsd.t",Xsd.t,
                        ".setting.", setting,
                        ".h",h,
                        ".sig.beta.delta", sig.beta.delta,
                        ".sig.Xmean.delta", sig.Xmean.delta,
                        ".sig.Xsd.delta", sig.Xsd.delta,
                        ".", sampling, 
                        ".", model, 
                        ".", fed,
                        ".Rdata"))}
}
}

save(res, file=paste0("simulation/result/sim.", 
                      "M", M,
                      ".n.s", n.s,
                      ".n.v", n.v,
                      ".n.t", n.t,
                      ".n.t.other", n.t.other,
                      ".prev.t", prev.t,
                      ".p",p,
                      ".s", s,
                      ".inter", inter,
                      ".int.beta", int.beta,
                      ".sig.beta.", sig.beta,
                      ".Xmean.t",Xmean.t,
                      ".Xsd.t",Xsd.t,
                      ".setting.", setting,
                      ".h",h,
                      ".sig.beta.delta", sig.beta.delta,
                      ".sig.Xmean.delta", sig.Xmean.delta,
                      ".sig.Xsd.delta", sig.Xsd.delta,
                      ".", sampling, 
                      ".", model, 
                      ".", fed,
                      ".Rdata"))


