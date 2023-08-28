# Trans-Balance: Reducing Demographic Disparity for Prediction Models in the Presence of Class Imbalance

`Transbalance` is an R package for estimating coefficients for risk prediction models in the presence of class imbalance and demographi disparity. It provides a comprehensive of functions for generating synthetic data, fitting Trans-Balance models, and evaluating their performance. The package aims to help researchers and practitioners in the field of risk prediction to deal with class imbalance and demographic disparity to improve the accuracy of their models, especially when working with limited labeled data.

## Features

- Synthetic data generation for imbalance class and demographic disparity
- Implements various models, including:
  - Proposed Trans-Balance Model with ROSE oversampling and federated learning 
  - Naive local
  - Naive global
  - Local model with ROSE oversampling
  - Global model with ROSE oversampling
  - Density ratio estimator
- Functions for calculating the predicted probabilities
- Performance evaluation metrics
- Extensive documentation and examples

## Installation

To install the latest version of `Transbalance` from GitHub, run the following commands in your R console:

```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("michaelyanwang/divideconquer")
devtools::install_github("ChuanHong/Transbalance")
```

## Usage

To get started with `Transbalance`, load the package in your R script or console:

```R
library(Transbalance)
```

Next, you can use the provided functions to generate synthetic data, fit survival models, and evaluate their performance. Here's a brief example:

```R
### Data generation
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
```

This example demonstrates how to generate synthetic data, fit different prediction models, and obtain final predicted probabilities estimates using the `Transbalance` package.
