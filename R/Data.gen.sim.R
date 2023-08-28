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
