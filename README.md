# SCANER: Semi-Supervised Calibration of Noisy Event Risk with Electronic Health Records

`SCANER` is an R package for estimating survival curves using both labeled and unlabeled data in a semi-supervised learning framework. It provides a comprehensive of functions for generating synthetic data, fitting survival models, and evaluating their performance. The package aims to help researchers and practitioners in the field of survival analysis to leverage the power of semi-supervised learning to improve the accuracy of their models, especially when working with limited labeled data.

## Features

- Synthetic data generation for survival analysis
- Implements various survival models, including:
  - Proposed semi-supervised estimator
  - Kaplan-Meier
  - Proposed SCANER estimator which ensembles the semi-supervised estimator and Kaplan-Meier
  - Deep Learning-based models (using the `deepsurv` package)
  - Non-parametric estimator
  - Density ratio estimator
- Functions for estimating survival curves
- Performance evaluation metrics
- Extensive documentation and examples

## Installation

To install the latest version of `SCANER` from GitHub, run the following commands in your R console:

```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("ChuanHong/SCANER")
```

## Usage

To get started with `SCANER`, load the package in your R script or console:

```R
library(SCANER)
```

Next, you can use the provided functions to generate synthetic data, fit survival models, and evaluate their performance. Here's a brief example:

```R
### Data generation
t0.all = seq(0.02,2,0.02); n.t0 = length(t0.all)
n     = 200
N     = 2000
lam   = 5

set.seed(1234)
data0 = data_generation(n, N, lam)
Xi = data0$Xi
Di = data0$Di
Ci = data0$Ci
Zi = data0$Zi
Ci.UL = data0$Ci.UL
Zi.UL = data0$Zi.UL

### Fit various survival models
KM = get.KM(Xi, Di, t0.all)
Semi = get.semi(Xi, Ci, Di, Zi, Ci.UL, Zi.UL, t0.all)
SCANER = get.SCANER(n.t0, n, Semi, KM)
DR = get.DR(Zi, Zi.UL, Xi, Ci, Di, t0.all)
DL = get.DL(Zi, Zi.UL, Xi, Ci, Di, t0.all)

### Calculate AUC of the proposed semi-supervised model
myAUC = get.semi.auc(Xi,Ci,Di,Zi,Ci.UL,Zi.UL, t0.all)
```

This example demonstrates how to generate synthetic data, fit different survival models, and obtain survival curve estimates using the `SCANER` package.
