
M=4
n.s=5000
n.v=5000
n.t=200
n.t.other=2000
prev.t=0.05
p=100
s=20
inter=1
int.beta=-5.025
sig.beta=0.296
auc=0.835
Xmean.t=0.218
Xsd.t=1.003
setting='shift.cov'
h=0
sig.beta.delta=0
sig.Xmean.delta=0.2
sig.Xsd.delta=0.2
sampling='Rose'
model='Source'
fed='Meta'

### parameter settings
Xmean.s=Xmean.t+sig.Xmean.delta
Xsd.s=Xsd.t+sig.Xsd.delta

Sig.X.t <- diag(Xsd.t, p)
Sig.X.s <- diag(Xsd.s, p)


set.seed(1234)
coef.all <-Coef.gen(s, h = h, sig.beta = sig.beta, sig.delta = sig.beta.delta, p = p) 
beta0.t=c(int.beta, coef.all$beta0)
beta0.s=c(int.beta, coef.all$W)

dat.junk=Data.gen.sim(M=M, n.s=n.s, n.t=n.t,n.t.other, n.v=n.v, p=p, beta0.s=beta0.s, beta0.t=beta0.t, mean.X.s=Xmean.s, mean.X.t=Xmean.t, Sig.X.s=Sig.X.s, Sig.X.t=Sig.X.t, inter)
dat.val=dat.junk$dat.v
dat.s.list=dat.junk$dat.s.list
dat.t.list=dat.junk$dat.t.list
dat.train.s=dat.s.list[["site1"]]
dat.train.t=dat.t.list[["site1"]]
auc.00=ROC.Est.FUN(dat.val$y, data.matrix(dat.val[,-1])%*%beta0.t[-1], yy0=0.5)[1]

betahat=main_kern(dat.train.s, dat.train.t, 
                   dat.s.list, dat.t.list,
                   nm.study="site1",
                   nm.cov=colnames(dat.train.s)[-1],
                   nm.out="y",
                   sampling=sampling,
                   model=model,
                   fed=fed,
                   nm.study.rm.fed="None", myseed=myseed)

shat=val.fun.new(betahat[-1], dat.val,nm.out="y")



