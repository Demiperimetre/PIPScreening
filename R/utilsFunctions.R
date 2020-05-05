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


