val.fun.new=function(mybeta, dat.val,nm.out){
  tmp.x=dat.val[,names(mybeta)]
  mys=data.matrix(tmp.x)%*%as.numeric(mybeta)
  list(auc=ROC.Est.FUN(dat.val[,nm.out], mys, yy0=0.5)[1],S=mys)
}