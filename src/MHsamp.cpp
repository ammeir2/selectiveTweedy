#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeFz(NumericVector z, NumericVector beta, bool log) {
  NumericVector result(z.length()) ;
  for(int i = 0 ; i < z.length() ; i ++) {
    double linearTerm = 0 ;
    double zi = z[i] ;
    for(int j = 0 ; j < beta.length() ; j ++) {
      linearTerm += beta[j] * std::pow(zi, j + 1) ;
    }

    if(log) {
      result[i] = linearTerm ;
    } else {
      result[i] = std::exp(linearTerm) ;
    }
  }

  return result ;
}

// [[Rcpp::export]]
NumericVector rejectSampCpp(int maxsamp,
                         int maxtries,
                         double lthreshold,
                         double uthreshold,
                         NumericVector M,
                         NumericVector theta,
                         NumericVector beta) {
  double threshold ;
  double sign ;
  NumericVector proposal(1) ;
  double fdens, lapDens ;
  double ratio ;
  NumericVector result(maxsamp) ;
  int slot = 0 ;

  int k = 0 ;

  for(int i = 0 ; i <  maxtries ; i++) {
    sign = 1.0 * (1 - 2*rbinom(1, 1, 0.5)[0]) ;
    if(sign > 0) {
      threshold = uthreshold ;
    } else {
      threshold = -lthreshold ;
    }
    proposal[0] = (threshold + rexp(1, theta[k])[0]) ;
    proposal[0] = proposal[0] * sign ;
    fdens = computeFz(proposal, beta, true)[0] ;
    lapDens = std::log(0.5) + std::log(theta[k]) - theta[k] * std::abs(proposal[0]) ;
    if(runif(1)[0] < std::exp(fdens - lapDens - std::log(M[k]))) {
      result[slot++] = proposal[0] ;
    }

    if(++k == theta.length()) k = 0 ;

    if(slot == maxsamp) break ;
  }

  return result ;
}


double sampleUnivTruncNorm(double mu, double sd, double threshold) {
  double u = runif(1)[0] ;
  double phiThreshold, sample ;

  phiThreshold = R::pnorm5(threshold, mu, sd, 1, 0) ;
  sample = R::qnorm5(u * phiThreshold, mu, sd, 1, 0) ;
  return sample ;
}


// [[Rcpp::export]]
NumericVector mhSampler(double init, NumericVector beta,
                        double lthreshold, double uthreshold,
                        double sampsd,
                        int burnin, int trim, int nsamp,
                        NumericVector tries) {
  NumericVector sample(nsamp) ;
  NumericVector samp(1) ;
  samp[0] = init ;
  NumericVector proposal(1) ;
  int iteration = 0 ;
  int sampNum = 0 ;
  int sign ;
  double threshold, MHratio ;
  double ppos, pneg ;
  double rate = 1 / (uthreshold - lthreshold) ;

  while(true) {
    if(runif(1)[0] < 0.2) {
      sign = 1 - 2 * rbinom(1, 1, 0.5)[0] ;
      if(sign > 0) {
        threshold = uthreshold ;
      } else {
        threshold = -lthreshold ;
      }
      proposal[0] = sign * (threshold + R::rexp(rate)) ;
      MHratio = computeFz(proposal, beta, true)[0] - computeFz(samp, beta, true)[0] ;
      if(samp[0] > uthreshold) {
        MHratio += R::dexp(samp[0] - uthreshold, rate, true) ;
      } else {
        MHratio += R::dexp(lthreshold - samp[0], rate, true) ;
      }
      if(proposal[0] > uthreshold) {
        MHratio -= R::dexp(proposal[0] - uthreshold, rate, true) ;
      } else {
        MHratio -= R::dexp(lthreshold - proposal[0], rate, true) ;
      }
      MHratio = std::exp(MHratio) ;
      if(runif(1)[0] < MHratio) {
        samp[0] = proposal[0] ;
      }
    } else {
      tries[0] += 1 ;
      ppos = R::pnorm(uthreshold, samp[0], sampsd, 0, 1) ;
      pneg = R::pnorm(lthreshold, samp[0], sampsd, 1, 1) ;
      ppos = 1 / (1 + std::exp(pneg - ppos)) ;
      if(runif(1)[0] < ppos) {
        proposal[0] = -sampleUnivTruncNorm(-samp[0], sampsd, -uthreshold) ;
      } else {
        proposal[0] = sampleUnivTruncNorm(samp[0], sampsd, lthreshold) ;
      }
      MHratio = computeFz(proposal, beta, true)[0] - computeFz(samp, beta, true)[0] ;
      MHratio = std::exp(MHratio) ;
      if(runif(1)[0] < MHratio) {
        tries[1] += 1 ;
        samp[0] = proposal[0] ;
      }
    }

    if(iteration > (burnin - 1) &
       ((iteration - (burnin - 1)) % trim == 0)) {
      sample[sampNum++] = samp[0] ;
      if(sampNum == nsamp) break ;
    }
    Rcpp::Rcout<<samp<<"\n" ;
    iteration++ ;
  }

  return sample ;
}



