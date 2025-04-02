#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List log_mh_cpp(int n, NumericVector par0, NumericVector sigmas,
            Function ll, Function lprior, Function lsampler,
            Function lsample_dens) {

  NumericVector acc(n);
  int m = par0.size();
  NumericMatrix pars(n, m);
  NumericVector par = par0;

  for(int i = 0; i < n; i++) {
    NumericVector park = lsampler(par0, par, sigmas);
    double a1 =  as<double>(ll(par)) - as<double>(ll(par));
    double a21 = as<double>(lprior(park, par0, sigmas));
    double a22 = as<double>(lprior(par, par0, sigmas));
    double a2 = a21 - a22;
    double a31 = as<double>(lsample_dens(par, par0, park, sigmas));
    double a32 = as<double>(lsample_dens(park, par0, par, sigmas));
    double a3 = a31 - a32;

    double a = a1 + a2 + a3;
    double alpha = exp(a);

    double u = unif_rand();

    if (u <= alpha) {
      acc[i] = 1;
      for(int j = 0; j < m; j++) {
        pars(i,j) = park[j];
      }
      par = park;
    }
    else {
      acc[i] = 0;
      for(int j = 0; j < m; j++) {
        pars(i,j) = par[j];
      }

    }
  }

  return List::create(
    Named("acceptance") = acc,
    Named("paramters") = pars
  );
}

// [[Rcpp::export]]
List mh_cpp(int n, NumericVector par0, NumericVector sigmas,
            Function likelihood, Function prior, Function sampler,
            Function sample_dens) {

  NumericVector acc(n);
  int m = par0.size();
  NumericMatrix pars(n, m);
  NumericVector par = par0;

  for(int i = 0; i < n; i++) {
    NumericVector park = sampler(par0, par, sigmas);
    double a1 =  as<double>(likelihood(par)) / as<double>(likelihood(par));
    double a21 = as<double>(prior(park, par0, sigmas));
    double a22 = as<double>(prior(par, par0, sigmas));
    double a2 = a21 / a22;
    double a31 = as<double>(sample_dens(par, par0, park, sigmas));
    double a32 = as<double>(sample_dens(park, par0, par, sigmas));
    double a3 = a31 / a32;

    double alpha = a1 * a2 * a3;

    double u = unif_rand();

    if (u <= alpha) {
      acc[i] = 1;
      for(int j = 0; j < m; j++) {
        pars(i,j) = park[j];
      }
      par = park;
    }
    else {
      acc[i] = 0;
      for(int j = 0; j < m; j++) {
        pars(i,j) = par[j];
      }

    }
  }

  return List::create(
    Named("acceptance") = acc,
    Named("paramters") = pars
  );
}
