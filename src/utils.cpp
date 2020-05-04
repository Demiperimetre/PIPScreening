#include <RcppArmadillo.h>
// [[Rcpp::depends(Matrix,RcppArmadillo)]]
#include <math.h>
#include <armadillo>


using namespace std;
using namespace arma;
using namespace Rcpp;



const double log2pi = std::log(2.0 * M_PI);

arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}




// scalaire puissance une matrice
arma::mat power(double v,arma::mat A)
{
  return exp(std::log(v)*mat(A));
}

// change les coordonnees de l'echelle logit vers 0,1
arma::vec expit(arma::vec a)
{
  return 1/(1+exp(-a));
}



// multivariate normal
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// correlation linkletter
arma::mat corfunctcpp(arma::vec rho,double alpha,Rcpp::List tensD)
{

  //Rcpp::Rcout << "res is now "  << tensD[0] << std::endl;
  int dx = rho.n_elem;
  arma::mat test = Rcpp::as<arma::mat>(tensD[0]);
  int n = test.n_rows;
  arma::mat res(n,n);
  res.ones();
  // Rcpp::Rcout << "res is now "  << test << std::endl;

  for (int i=0;i<dx;i++)
  {
    arma::mat prov = pow(2*Rcpp::as<arma::mat>(tensD[i]),alpha);
    res = res % power(rho(i),prov);
  }

  return res;

}


// prior dbeta
double priorrho(arma::vec gamma, arma::vec rho)
{
  double res=1;
  int dx = rho.n_elem;

  for (int i;i<dx;i++)
  {
    res*= R::dbeta(rho(i),gamma(i),1,false);
  }
  return res;
}


// prior dbeta
double loglikcpp(arma::vec rho,double sigerr, double sigdelta,arma::mat Rexp,Rcpp::List tensD,double alpha)
{
  int n = Rexp.n_rows;
  // Rcpp::Rcout << "res is now "  << n << std::endl;

  arma::mat id(n,n,fill::eye);
  arma::mat matvar = sigerr * id + sigdelta * corfunctcpp(rho,alpha,tensD);
  // attention calcul du determinant peut mener a des inf
  //double res = - .5 * log(det(matvar)) - .5 * as_scalar(Rexp.t() * solve(matvar,Rexp));
  double res =  - .5 * as_scalar(Rexp.t() * solve(matvar,Rexp));
  arma::vec vp = arma::eig_sym(matvar);
  double logdeterminant = sum(log(vp));
  res = res - 0.5 * logdeterminant;
  return res;
}


// Posterior computation
// calcule toujours en echelle transformee
// renvoie vercteur taille 2 log prior et loglik
// parprior matrice nbpar*2
// pour l instant on ne calibre pas, a faire en dehors ?
// on donne la longueur de rho afin d'eviter les probleme si on est en calibration et que les priors ne sont pas unif(0,1) pour theta
arma::vec logposteriorcpp(arma::vec parameters,arma::vec Rexp,Rcpp::List tdensD,double alpha,arma::mat parprior,int drho)
{
  int n = Rexp.n_elem;
  //Rcpp::Rcout << "res is now "  << n << std::endl;
  int dpar = parameters.n_elem;// = drho + 2 (var) + nb theta

  arma::vec retour(2);

  arma::vec rho(drho);
  arma::vec rhobrut(drho);

  // for (int i=0; i<(npar-2);i++)  {rho(i) = 1/(1+exp(-parameters(i)));}
  rhobrut = parameters(span(0,drho-1));
  rho = expit(rhobrut);
  double sigmaerr = exp(parameters(drho));
  double sigmadelta = exp(parameters(drho+1));


  // log prior of inverse gamma
  double logpriorsigerr = R::dgamma(1/sigmaerr,parprior(drho,0) ,1/parprior(drho,1),true) - log(pow(sigmaerr,2));// calcul inv gamma a partir gamma
  double logpriorsigdelta = R::dgamma(1/sigmadelta,parprior(drho+1,0) ,1/parprior(drho+1,1),true) - log(pow(sigmadelta,2));

  //log prior of rho
  double logpriorrho = log(priorrho(parprior(span(0,drho-1),0),rho));



  //log lik
  double loglik = loglikcpp(rho,sigmaerr,sigmadelta,Rexp,tdensD,alpha);


  // check if cal, if theta are out of bounds
  bool caloutbounds=false;
  if (dpar>drho+2)
  {
    int dtheta = dpar-drho-2;
    if (any(parameters(span(dpar-dtheta,dpar-1))<-15) or any(parameters(span(dpar-dtheta,dpar-1))>15))
    {caloutbounds=true;}
  }

  //  Rcpp::Rcout << "res is now "  << rho << std::endl;
  if (any(rhobrut>15) or any(rhobrut< -15) or caloutbounds) {retour(0) = - std::numeric_limits<double>::infinity();
    //  Rcpp::Rcout << "res is now "  << rho << std::endl;
  }
  else {retour(0)=logpriorsigerr+logpriorsigdelta+logpriorrho;}
  retour(1)=loglik;
  //double retour=logpriorsigerr+logpriorsigdelta+logpriorrho+loglik;
  return retour;
}



double logposteriorscalcpp(arma::vec parameters,arma::mat Rexp,Rcpp::List tdensD,double alpha,arma::mat parprior,int drho)
{
  arma::vec res=logposteriorcpp(parameters,Rexp,tdensD,alpha,parprior,drho);
  return res(0)+res(1);
}


// Metropolis within gibss,
//burnin a enlever a posterior
// les paramètres sont dans l'espace transformés
//' Metropolis within Gibbs Sampler
//'
//' @param niter number of iterations
//' @return Posterior sample in the transformed space
//' @export
// [[Rcpp::export]]
Rcpp::List MetropoliswGibbs(int niter,arma::vec parwalk,arma::vec parinit,arma::vec Rexp,Rcpp::List tdensD,double alpha,arma::mat parprior,bool adaptive,Rcpp::List calibration)
{
  bool cal=calibration[3];
  //if (cal)
  //{
  Rcpp::Function model = calibration[0];
  arma::vec Y = calibration[1];
  arma::mat X = calibration[2];
  int dx = X.n_cols; // dx=drho
  int dpar = parinit.n_elem;
  int dtheta = dpar - dx - 2;// when nb of theta notequal number of x
  //}

  //Rcpp::Rcout << "res is now "  << cal << std::endl;

  int iteradap = 300;
  int maxadap = 2000;


  arma::mat chain(niter,dpar);
  arma::vec courant=parinit;
  arma::vec prop=courant;
  arma::vec accep(dpar); //1 par dim acceptation globale
  accep.zeros();
  arma::vec accepprov(dpar);
  accepprov.zeros();


  //vec test(3);
  //test(0)=1;
  //test(1)=2;
  //test(2)=3;
  //Rcpp::Rcout << "res is now "  << Rcpp::as<vec>(model(X,test))<< std::endl;

  arma::vec Rexp1;
  if (cal)
  {


    // dernier parametres du vecteur dans la bonne echelle
    Rexp = Y- Rcpp::as<vec>(model(X,expit(courant(span(dpar-dtheta,dpar-1)))));
    //Rcpp::Rcout << "res is now "  << expit(courant(span(dpar-dx,dpar-1))) << std::endl;
    Rexp1=Rexp;
  }



  double logpostcourante = logposteriorscalcpp(parinit,Rexp,tdensD,alpha,parprior,dx);
  double logpostprop;
  double logratio;

  for (int i=0;i<niter;i++)
  {
    for (int k=0;k<dpar;k++)
    {
      // RW
      prop(k) = courant(k) + randn() * sqrt(parwalk(k)) ;
      // calcul de Rexp necessaire en cas de calib pour la logpost
      if (cal)
      {
        // dernier parameters du vecteur dans la bonne echelle
        Rexp = Y- Rcpp::as<vec>(model(X,expit(prop(span(dpar-dtheta,dpar-1)))));
      }
      logpostprop = logposteriorscalcpp(prop,Rexp,tdensD,alpha,parprior,dx);

      // Accept ou reject
      logratio = logpostprop - logpostcourante;
      if (log(randu()) < logratio and std::isinf(logpostprop)==false) // ajout pour empecher d accepter loglik infinie
      {
        courant(k) = prop(k);
        logpostcourante = logpostprop;
        accep(k)++;
        accepprov(k)++;
      }
      chain(i,k) = courant(k);
    }

    // adaptation
    if ((i+1)%iteradap==0 and  adaptive==true and i<maxadap)
    {
      for (int k=0;k<dpar;k++)
      {
        if (accepprov(k)/iteradap<.2) {parwalk(k) *= (1-.1); }
        if (accepprov(k)/iteradap>.5) {parwalk(k) *= (1+.1);  }
      }
      accepprov.zeros();
    }

  }

  return Rcpp::List::create(Rcpp::Named("chain")=chain,Rcpp::Named("AccepRate")=accep/niter,Rcpp::Named("parwalk")=parwalk,Rcpp::Named("test")=Rexp1);
}


//MH classique avec matrice de cov a estimer prealablement ?
// [[Rcpp::export]]
Rcpp::List Metropolis(int niter,arma::mat covwalk,arma::vec parinit,arma::vec Rexp,Rcpp::List tdensD,double alpha,arma::mat parprior,bool adaptive,Rcpp::List calibration)
{
  bool cal=calibration[3];
  //if (cal)
  //{
  Rcpp::Function model = calibration[0];
  arma::vec Y = calibration[1];
  arma::mat X = calibration[2];
  int dx = X.n_cols;
  //}
  int dpar = parinit.n_elem;
  int dtheta = dpar - dx - 2;// when nb of theta notequal number of x
  //Rcpp::Rcout << "res is now "  << cal << std::endl;

  // par adaptation
  int iteradap = 300;
  int maxadap = 2000;

  int n = Rexp.n_elem;

  arma::mat chain(niter,dpar);
  arma::vec vlogpost(niter); // stocker le calcul de la log vraisemblance
  arma::vec courant=parinit;
  arma::vec prop=courant;
  double accep=0;
  double accepprov=0;

  if (cal)
  {
    //Rcpp::Rcout << "res is now "  << Rcpp::as<vec>(model(X,expit(courant(span(dpar-dtheta,dpar-1)))))<< std::endl;
    // dernier parametres du vecteur dans la bonne echelle
    Rexp = Y- Rcpp::as<vec>(model(X,expit(courant(span(dpar-dtheta,dpar-1)))));

  }


  //Rcpp::Rcout << "res is now "  << expit(courant(span(dpar-dx,dpar-1))) << std::endl;

  double logpostcourante = logposteriorscalcpp(parinit,Rexp,tdensD,alpha,parprior,dx);
  double logpostprop;
  double logratio;
  arma::vec muRW(dpar);
  muRW.zeros();

  for (int i=0;i<niter;i++)
  {
    // RW
    prop = courant + mvrnormArma(1,muRW , covwalk).t();// mvrnormArma(1,muRW , covwalk).row(0);

    //Rcpp::Rcout << "res is now prop"  << prop << std::endl;

    // calcul de Rexp necessaire en cas de calib pour la logpost
    if (cal)
    {

      //Rcpp::Rcout << "res is now "  << Rcpp::as<vec>(model(X,expit(prop(span(dpar-dtheta,dpar-1)))))<< std::endl;
      // dernier parameters du vecteur dans la bonne echelle
      Rexp = Y- Rcpp::as<vec>(model(X,expit(prop(span(dpar-dtheta,dpar-1)))));
    }
    logpostprop = logposteriorscalcpp(prop,Rexp,tdensD,alpha,parprior,dx);

    // Accept ou reject
    logratio = logpostprop - logpostcourante;
    //Rcpp::Rcout << "res is now "  << logratio << std::endl;
    if (log(randu()) < logratio and std::isinf(logpostprop)==false) // ajout pour empecher valeur infinies
    {
      courant = prop;
      logpostcourante = logpostprop;
      accep++;
      accepprov++;
      //  Rcpp::Rcout << "accepte" << accep << std::endl;
    }
    chain.row(i) = courant.t();
    vlogpost(i) = logpostcourante;

    // adaptation
    if ((i+1)%iteradap==0 and  adaptive==true and i<maxadap)
    {

      if (accepprov/iteradap<.2) {covwalk *= (1-.1); }
      if (accepprov/iteradap>.5) {covwalk *= (1+.1);  }
      //    Rcpp::Rcout << "accepte" << accepprov << std::endl;
      //  Rcpp::Rcout << "cov" << covwalk << std::endl;
      accepprov=0;
    }


  }

  //double test= accep/niter;

  //Rcpp::Rcout << "res is now "  << test << std::endl;
  return Rcpp::List::create(Rcpp::Named("chain")=chain,Rcpp::Named("AccepRate")=accep/niter,Rcpp::Named("parwalk")=covwalk,Rcpp::Named("vlogpost")=vlogpost);
}


// gal mcmc
// [[Rcpp::export]]
Rcpp::List MCMC(int niterMwG,int niterMH,arma::vec parwalk,arma::vec parinit,arma::vec Rexp,Rcpp::List tdensD,double alpha,arma::mat parprior,bool adaptive,Rcpp::List calibration)
{
  int burnin = 100; // en accord avec le max d'adaptation du MwG OPTION a ajouter
  int dpar = parinit.n_elem;

  //Rcpp::Rcout << "res is now "  << parinit << std::endl;

  // on lance d a bord  MwG
  Rcpp::List resMwG = MetropoliswGibbs(niterMwG,parwalk,parinit,Rexp,tdensD,alpha,parprior,adaptive,calibration);


  //Rcpp::Rcout << "res is now "<< niterMwG << std::endl;

  arma::mat chain = resMwG["chain"];
  arma::mat subchain = chain(span(burnin,niterMwG-1),span(0,dpar-1));

  // estimation de la covariance
  arma::mat covwalk = cov(subchain);

  //Rcpp::Rcout << "res is now "  << subchain(0,span(0,dpar-1)) << std::endl;
  //Rcpp::Rcout << "res is now "  << chain.row(niterMwG-1) << std::endl;
  arma::vec accep=resMwG["AccepRate"];
  //  Rcpp::Rcout << "accep is now "  << accep << std::endl;

  arma::vec parinit2 = chain.row(niterMwG-1).t();
  //Rcpp::Rcout << "res is now "  << parinit2 << std::endl;

  Rcpp::List resMH = Metropolis(niterMH,covwalk/2,parinit2,Rexp,tdensD,alpha,parprior,adaptive,calibration);




  return Rcpp::List::create(Rcpp::Named("matcovMwG")=covwalk,Rcpp::Named("MH")=resMH);
}

/*

// gal mcmc
// [[Rcpp::export]]
vec Chibcomputationscpp(vec parstar,Rcpp::List resMH,vec Rexp,Rcpp::List tdensD,double alpha,mat parprior,Rcpp::List calibration)
{


  int dpar = parstar.n_elem;
  mat sample = resMH["chain"];
  vec vlogpostsample = resMH["vlogpost"];
  mat covwalk = resMH["parwalk"];
  int nsample = sample.n_rows;


  // calibration part
  bool cal=calibration[3];
  Rcpp::Function model = calibration[0];
  vec Y = calibration[1];
  mat X = calibration[2];
  int dx = X.n_cols;

  int dtheta = dpar - dx - 2;

  // posterior for parstar
  if (cal)
  {
    // dernier parameters du vecteur dans la bonne echelle
    Rexp = Y- Rcpp::as<vec>(model(X,expit(parstar(span(dpar-dtheta,dpar-1)))));
  }
  vec vlogpoststar = logposteriorcpp(parstar,Rexp,tdensD,alpha,parprior,dx);
  double logpoststar = sum(vlogpoststar);



  // numerator
  mat MATmin(nsample,2);
  MATmin.col(0) = exp(logpoststar-vlogpostsample);
  MATmin.col(1) = ones(nsample) ;
  vec Accepalpha = min(MATmin,1);
  vec qtrans(nsample); // proba de transition de parsample vers parstar
  qtrans = dmvnrm_arma(sample,parstar.t(),covwalk,false);
  double numerator = mean(qtrans%Accepalpha);
  // Rcpp::Rcout << "res is now "  << numerator << std::endl;

  //Rcpp::Rcout << "res is now "  << qtrans << std::endl;


  //Rcpp::Rcout << "tudo bem "   << std::endl;

  // denominator
  mat sampfromstar =  mvrnormArma(nsample,parstar , covwalk);
  vec logpostsampfromstar(nsample);
  vec vsamp(dpar);
  //Rcpp::Rcout << "res is now "  << size(vsamp) << std::endl;
  //Rcpp::Rcout << "res is now "  << size(sampfromstar) << std::endl;
  for (int i;i<nsample;i++)
  {
    vsamp = sampfromstar.row(i).t();
    if (cal)
    {
      // dernier parameters du vecteur dans la bonne echelle
      Rexp = Y- Rcpp::as<vec>(model(X,expit(vsamp(span(dpar-dtheta,dpar-1)))));
    }
    logpostsampfromstar(i) = logposteriorscalcpp(vsamp,Rexp,tdensD,alpha,parprior,dx);
  }
  mat MATminden(nsample,2);
  MATminden.col(0) = exp(logpostsampfromstar-logpoststar);
  MATminden.col(1) = ones(nsample) ;
  vec Accepalphaden = min(MATminden,1);
  double denominator = mean(Accepalphaden);

  vec retour(4);
  retour(0) = vlogpoststar(0);
  retour(1) = vlogpoststar(1);
  retour(2) = numerator;
  retour(3) = denominator;

  return retour;
}
*/
