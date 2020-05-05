#' Distances per coordinates
#'
#' @param X1 design matrix.
#'
#' @return List of distances per coordinates
#' @export
#'
#' @examples
#' @importFrom stats dist
tensordist = function(X1)
{
  lapply(1:ncol(X1), function(i) as.matrix(dist(X1[,i],"manhattan")))
}

#' Distances per coordinates
#'
#' @param ech posterior sample
#'
#' @return List of distances per coordinates
#'
#' @examples
#' @importFrom stats dbeta
BFbridge = function(ech,gamma,gammaref=rep(1,length(gamma)))
{
  n = nrow(ech)
  d = length(gamma)
  # gamma= as.numeric(gamma)
  rap = numeric(n)
  for (i in 1:n)
  {
    rap[i] = exp(sum(dbeta(expit2(ech[i,1:d]),gamma,1,log = T)) - sum(dbeta(expit2(ech[i,1:d]),gammaref,1,log = T)) )
  }
  return(mean(rap))
}



