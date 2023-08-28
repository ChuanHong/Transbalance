Data.gen2=function(n, p, beta.x, mean.X, Sig.X,inter){
  x <- rmvnorm(n, rep(mean.X, p), Sig.X)
  if(inter==0){
    y=unlist(lapply(1:n, function(ll) rbinom(1,1,g.logit(c(1,x[ll,]) %*% beta.x))))}
  if(inter==1){
    y=unlist(lapply(1:n, function(ll) rbinom(1,1,g.logit(c(1,x[ll,]) %*% beta.x+x[ll,1]*x[ll,2]))))}
  data.frame(y, x)
}