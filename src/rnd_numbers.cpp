// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

//' @rdname dnorminvwish
//' @keywords internal
// [[Rcpp::export]]
arma::mat rmatn(arma::mat M, arma::mat Q, arma::mat P){
/*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  int p = P.n_rows;
  int q = Q.n_rows;
  arma::mat L = arma::chol(Q, "upper");
  arma::mat C = arma::chol(P, "lower");
  arma::mat X = arma::reshape(arma::vec(Rcpp::rnorm(p * q)), p, q);
  X = M + C * X * L;
  return(X);
}


//' @rdname dnorminvwish
//' @keywords internal
// [[Rcpp::export]]
arma::mat rinvwish(int v, arma::mat S){
  int p = S.n_rows;
  arma::mat L = arma::chol(arma::inv_sympd(S), "lower");
  arma::mat A(p,p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df));
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  arma::mat LA_inv = arma::inv(arma::trimatl(arma::trimatl(L) * arma::trimatl(A)));
  arma::mat X = LA_inv.t() * LA_inv;

  return(X);
}


//' @rdname dmultn
//' @keywords internal
// [[Rcpp::export]]
arma::vec rmultn(arma::vec m, arma::mat Sigma){
  /*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  int p = Sigma.n_rows;
  arma::vec X = Rcpp::rnorm(p);
  arma::mat L = arma::chol(Sigma, "lower");
  X = m + L * X;
  return(X);
}

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmultn(arma::mat x,
                      arma::rowvec m,
                      arma::mat Sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - m) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double dmatn2(arma::mat X, arma::mat M, arma::mat Q, arma::mat P,
             bool logd = false){
  /*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  int n = X.n_rows;
  int k = X.n_cols;
  double Q_logdet;
  double P_logdet;
  double Q_sign;
  double P_sign;
  arma::log_det(Q_logdet, Q_sign, Q);
  arma::log_det(P_logdet, P_sign, P);
  arma::mat ss = X - M;
  double exponent = -0.5 * arma::trace(arma::inv_sympd(Q) * ss.t() * arma::inv_sympd(P) * ss);
  double con = -((log(2 * M_PI)/2) * n * k + (Q_logdet/2) * n + (P_logdet/2) * k);
  double out = exponent + con;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double dmatn(arma::mat X, arma::mat M, arma::mat Q, arma::mat P,
             bool logd = false){
  /*-------------------------------------------------------
# Generate draws from a matricvariate normal distribution
#-------------------------------------------------------*/
  int n = X.n_rows;
  int k = X.n_cols;
  arma::mat Q_C = arma::chol(Q);
  arma::mat Q_C_inv = arma::inv(arma::trimatu(Q_C));
  arma::mat P_C = arma::chol(P);
  arma::mat P_C_inv = arma::inv(arma::trimatu(P_C));
  arma::mat ss = X - M;
  double Q_logdet = sum(log(arma::diagvec(Q_C)));
  double P_logdet = sum(log(arma::diagvec(P_C)));
  double exponent = -0.5 * arma::trace(Q_C_inv * Q_C_inv.t() * ss.t() * P_C_inv * P_C_inv.t() * ss);
  double con = -((log(2 * M_PI)/2) * n * k + (Q_logdet) * n + (P_logdet) * k);
  double out = exponent + con;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//' @rdname dnorminvwish
//' @keywords internal
// [[Rcpp::export]]
double dinvwish(arma::mat Sigma, int v, arma::mat S, bool logd = false){
  int q = Sigma.n_cols;
  double k = 0.5*v*q*log(2) + 0.25*q*(q-1)*log(M_PI);
  double lgam = 0;
  for (int i = 1; i <= q; i++) {
    lgam = lgam + lgamma((v+1-i)*0.5);
  }
  k = k + lgam;
  double A_logdet;
  double B_logdet;
  double A_sign;
  double B_sign;
  arma::log_det(A_logdet, A_sign, Sigma);
  arma::log_det(B_logdet, B_sign, S);

  double kern = 0.5*v*B_logdet -0.5*(v+q+1)*A_logdet-0.5*arma::trace(arma::inv_sympd(Sigma) * S);
  double out = -k + kern;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}
