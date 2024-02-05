// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// loglikcpp
double loglikcpp(arma::vec rho, double sigerr, double sigdelta, arma::mat Rexp, Rcpp::List tensD, double alpha);
RcppExport SEXP _PIPScreening_loglikcpp(SEXP rhoSEXP, SEXP sigerrSEXP, SEXP sigdeltaSEXP, SEXP RexpSEXP, SEXP tensDSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type sigerr(sigerrSEXP);
    Rcpp::traits::input_parameter< double >::type sigdelta(sigdeltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Rexp(RexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tensD(tensDSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikcpp(rho, sigerr, sigdelta, Rexp, tensD, alpha));
    return rcpp_result_gen;
END_RCPP
}
// MetropoliswGibbs
Rcpp::List MetropoliswGibbs(int niter, arma::vec parwalk, arma::vec parinit, arma::vec Rexp, Rcpp::List tdensD, double alpha, arma::mat parprior, bool adaptive, Rcpp::List calibration);
RcppExport SEXP _PIPScreening_MetropoliswGibbs(SEXP niterSEXP, SEXP parwalkSEXP, SEXP parinitSEXP, SEXP RexpSEXP, SEXP tdensDSEXP, SEXP alphaSEXP, SEXP parpriorSEXP, SEXP adaptiveSEXP, SEXP calibrationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parwalk(parwalkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parinit(parinitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Rexp(RexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tdensD(tdensDSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type parprior(parpriorSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type calibration(calibrationSEXP);
    rcpp_result_gen = Rcpp::wrap(MetropoliswGibbs(niter, parwalk, parinit, Rexp, tdensD, alpha, parprior, adaptive, calibration));
    return rcpp_result_gen;
END_RCPP
}
// Metropolis
Rcpp::List Metropolis(int niter, arma::mat covwalk, arma::vec parinit, arma::vec Rexp, Rcpp::List tdensD, double alpha, arma::mat parprior, bool adaptive, Rcpp::List calibration);
RcppExport SEXP _PIPScreening_Metropolis(SEXP niterSEXP, SEXP covwalkSEXP, SEXP parinitSEXP, SEXP RexpSEXP, SEXP tdensDSEXP, SEXP alphaSEXP, SEXP parpriorSEXP, SEXP adaptiveSEXP, SEXP calibrationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covwalk(covwalkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parinit(parinitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Rexp(RexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tdensD(tdensDSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type parprior(parpriorSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type calibration(calibrationSEXP);
    rcpp_result_gen = Rcpp::wrap(Metropolis(niter, covwalk, parinit, Rexp, tdensD, alpha, parprior, adaptive, calibration));
    return rcpp_result_gen;
END_RCPP
}
// MCMC
Rcpp::List MCMC(int niterMwG, int niterMH, arma::vec parwalk, arma::vec parinit, arma::vec Rexp, Rcpp::List tdensD, double alpha, arma::mat parprior, bool adaptive, Rcpp::List calibration);
RcppExport SEXP _PIPScreening_MCMC(SEXP niterMwGSEXP, SEXP niterMHSEXP, SEXP parwalkSEXP, SEXP parinitSEXP, SEXP RexpSEXP, SEXP tdensDSEXP, SEXP alphaSEXP, SEXP parpriorSEXP, SEXP adaptiveSEXP, SEXP calibrationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type niterMwG(niterMwGSEXP);
    Rcpp::traits::input_parameter< int >::type niterMH(niterMHSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parwalk(parwalkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parinit(parinitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Rexp(RexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tdensD(tdensDSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type parprior(parpriorSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type calibration(calibrationSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC(niterMwG, niterMH, parwalk, parinit, Rexp, tdensD, alpha, parprior, adaptive, calibration));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PIPScreening_loglikcpp", (DL_FUNC) &_PIPScreening_loglikcpp, 6},
    {"_PIPScreening_MetropoliswGibbs", (DL_FUNC) &_PIPScreening_MetropoliswGibbs, 9},
    {"_PIPScreening_Metropolis", (DL_FUNC) &_PIPScreening_Metropolis, 9},
    {"_PIPScreening_MCMC", (DL_FUNC) &_PIPScreening_MCMC, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_PIPScreening(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
