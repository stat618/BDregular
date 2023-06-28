dmean <- function(A){
  a = colMeans(A)
  A = t(t(A) - a)
  return(A)
}

# 标准化
standard=function(A){
  a=colMeans(A)
  b=apply(A,2,sd)
  B=t(t(A)-a)
  B=t(t(B)/b)
  return(B)
}
# 将两个矩阵合并为一个块对角矩阵
blockdiag=function(A,B){
  if(length(A)>1){
    p1=dim(A)[2];p2=dim(B)[2]
    C=diag(p1+p2)
    C[1:p1,1:p1]=A
    C[c((p1+1):(p1+p2)),c((p1+1):(p1+p2))]=B
  }else{
    p1=1;
    p2=dim(B)[2]
    C=diag(p1+p2)
    C[1,1]=A
    C[c(2:(1+p2)),c(2:(1+p2))]=B
  }
  return(C)
}

# SCAD的阈值函数
dSCAD=function(a,lam,gamma=3.7){
  a=abs(a)
  z=a
  z[a<lam]=lam
  z[a>lam]=(gamma*lam-z[a>lam])/(gamma-1)
  z[a>(gamma*lam)]=0
  return(z)
}

dMCP=function(a,lam,gamma=3){
  a=abs(a)
  z=lam-a/gamma
  z[a>(gamma*lam)]=0
  return(z)
}

soft=function(a,b){
  c=abs(a)-b
  c[c<0]=0
  c=c*sign(a)
  return(c)
}

shrinking=function(A,eps){
  d=eigen(A)$values
  d1=min(d)
  dp=max(d)
  us=max(eps,(d1+dp)/2)
  uf=sum((d-d1)^2)/sum(d-d1)
  u=max(us,uf)
  alpha=1-(eps-d1)/(u-d1)
  B=alpha*as.matrix(A)+(1-alpha)*u*diag(ncol(A))
  B=cov2cor(B)
  S=list(est=B,alpha=alpha,u=u)
  return(S)
}

SCADthreshold=function(S,lam,k=3){
  a=as.vector(S)
  b=abs(a)
  for(i in 1:k){
    c=dSCAD(b,lam)
    b=soft(b,c)
  }
  d=sign(a)*b
  d=matrix(d,ncol(S),ncol(S))
  diag(d)=1
  d=shrinking(d,eps=0.001)
  d=cov2cor(d$est)
  return(d)
}

MCPthreshold=function(S,lam,k=3){
  a=as.vector(S)
  b=abs(a)
  for(i in 1:k){
    c=dMCP(b,lam)
    b=soft(b,c)
  }
  d=sign(a)*b
  d=matrix(d,ncol(S),ncol(S))
  diag(d)=1
  d=shrinking(d,eps=0.001)
  d=cov2cor(d$est)
  return(d)
}

# library('Rcpp')
# sourceCpp(file = 'R/a.cpp')
