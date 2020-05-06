#' Distances per coordinates
#'
#' @param X1 design matrix.
#'
#' @return List of distances per coordinates
#' @export
#'
#' @examples
#' @importFrom stats dist
tensordist  <-  function(X1)
{
  lapply(1:ncol(X1), function(i) as.matrix(dist(X1[,i],"manhattan")))
}

#' Compute Bayes Factor from Importance Sampling
#'
#' @param ech posterior sample in the logit scale
#' @param gamma values given the alpha parameter of the prior distribution on the range parameters
#' @param gammaref values of the reference model, by default 1 for all range parameters
#'
#' @return Bayes factor between model with prior given in gamma and reference model (gamma=1)
#'
#' @export
#' @examples
#' @importFrom stats dbeta
BFbridgeIS <-  function(ech,gamma,gammaref=rep(1,length(gamma)))
{
  n <-  nrow(ech)
  d <-  length(gamma)
  # gamma= as.numeric(gamma)
  rap <- numeric(n)
  for (i in 1:n)
  {
    echtransf <- 1/(1+exp(-ech[i,1:d]))
    rap[i] <- exp(sum(dbeta(echtransf,gamma,1,log = T)) - sum(dbeta(echtransf,gammaref,1,log = T)) )
  }
  return(mean(rap))
}



#' Compute the Probabilities that Input Variables are Active
#'
#' @param ech posterior sample in the logit scale
#' @param valmax value for alpha parameter for the inert prior distribution
#'
#' @return vector of activeness probabilities
#' @export
#'
#' @examples
computeProbActive <- function(ech,valmax=100)
{
  p <- ncol(ech)
  GAMMA  <-  expand.grid(replicate(p, c(1,valmax), simplify=FALSE))
  BFall <-  sapply(1:nrow(GAMMA),function(k) {
    BFbridgeIS(ech,gamma  <-  as.numeric(GAMMA[k,]))})

  IndicVar  <-  matrix(NA,nrow = p, ncol = 2^p)
  for (i in 1:p)
  {
    IndicVar[i,] <-  (GAMMA[,i]==1)*1
  }

  prob <-  IndicVar %*% BFall / sum(BFall)
  return(as.vector(prob))
}
