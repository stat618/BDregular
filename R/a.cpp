
#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cov2cor(arma::mat cov) {

  arma::vec sds = sqrt(cov.diag());
  arma::mat outer_sds = sds * sds.t();
  arma::mat cor = cov / outer_sds;
  cor.diag().ones();
  return cor;
}

// [[Rcpp::export]]
arma::mat dmean(arma::mat A){
  arma::vec a = arma::mean(A,0).t();
  A = A.each_row() - a.t();
  return(A);
}
// [[Rcpp::export]]
arma::mat standard(arma::mat A){
  arma::vec a = arma::mean(A,0).t();
  arma::vec b = arma::stddev(A,0,0).t();
  A = A.each_row() - a.t();
  A = A.each_row() / b.t();
  return(A);
}
// [[Rcpp::export]]
arma::mat blockdiag(arma::mat A,arma::mat B){
  int p1 = A.n_cols;
  int p2 = B.n_cols;
  arma::mat C(p1+p2,p1+p2,arma::fill::zeros);
  C.submat(0,0,p1-1,p1-1) = A;
  C.submat(p1,p1,p1+p2-1,p1+p2-1) = B;
  return(C);
}
// [[Rcpp::export]]
arma::mat dSCAD(arma::mat a,double lam,double gamma=3.7){
  a = arma::abs(a);
  arma::mat z = a;
  z.elem(arma::find(a<lam)).fill(lam);
  z.elem(arma::find(a>lam)) = (gamma*lam-z.elem(arma::find(a>lam)))/(gamma-1);
  z.elem(arma::find(a>gamma*lam)).fill(0);
  return(z);
}
// [[Rcpp::export]]
arma::mat dMCP(arma::mat a,double lam,double gamma=3){
  a = arma::abs(a);
  arma::mat z = lam-a/gamma;
  z.elem(arma::find(a>gamma*lam)).fill(0);
  return(z);
}
// [[Rcpp::export]]
arma::mat soft(arma::mat a,double b){
  arma::mat c = arma::abs(a)-b;
  c.elem(arma::find(c<0)).fill(0);
  c = c % arma::sign(a);
  return(c);
}
// [[Rcpp::export]]
arma::mat shrinking(arma::mat A,double eps){
  arma::vec d = arma::eig_sym(A);
  double d1 = d.min();
  double dp = d.max();
  double us = std::max(eps,(d1+dp)/2);
  double uf = arma::sum(arma::pow(d-d1,2))/arma::sum(d-d1);
  double u = std::max(us,uf);
  double alpha = 1-(eps-d1)/(u-d1);
  arma::mat B = alpha*A+(1-alpha)*u*arma::eye(A.n_cols,A.n_cols);
  B = cov2cor(B);
  arma::mat S = arma::join_rows(B,arma::join_rows(arma::join_cols(alpha),arma::join_cols(u)));
  return(S);
}
// [[Rcpp::export]]
arma::mat SCADthreshold(arma::mat S,double lam,int k=3){
  arma::vec a = arma::vectorise(S);
  arma::vec b = arma::abs(a);
  for(int i=0;i<k;i++){
    arma::mat c = dSCAD(b,lam);
    b = soft(b,c);
  }
  arma::vec d = arma::sign(a) % b;
  arma::mat D = arma::diagmat(d);
  D = shrinking(D,0.001);
  D = cov2cor(D);
  return(D);
}
