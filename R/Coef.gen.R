Coef.gen<- function(s, h, sig.beta, sig.delta, p){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- beta0
  samp0<- sample(1:s, h, replace=F)
  W[samp0] <-W[samp0] + rep(-sig.delta, h)
  return(list(W=W, beta0=beta0))
}