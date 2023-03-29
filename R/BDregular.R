
#' @title Accuracy matrix block diagonal regularization estimation
#' @description description demo
#' @details details demo
#' @param x param 1
#'
#' @return Precision matrix estimate
#' @export
#'
#' @examples BDregular(X)
BDregular = function(X,dim=nrow(X),lam="ada",shrink="MCP",min.eig=0.001){

  S=cov(X)
  D=sqrt(diag(S))
  R=cov2cor(S)
  if(lam=="ada"){
    lam=sqrt(log(ncol(X))/nrow(X))
  }
  if(shrink=="MCP"){
    R1=MCPthreshold(R,lam)
  }
  if(shrink=="SCAD"){
    R1=SCADthreshold(R,lam)
  }
  diag(R1)=0
  s=colSums(abs(R1))
  ind=order(s,decreasing=T)[1:dim]
  R2=R1[ind,ind]
  S2=t(R2*D)*D
  return(S2)

}
