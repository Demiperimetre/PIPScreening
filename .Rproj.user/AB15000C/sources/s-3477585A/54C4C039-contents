#' Non Linear Simulator with Calibration Parameter in 3 dimension
#'
#' @param x Matrix of field design in dimension n \eqn{\times} p (with p>3).
#' Only the 3 first rows of x are taken into account
#' @param theta Calibration parameter of dimension 3
#'
#' @return vector of responses
#' @export
#'
#' @examples
sim3 <- function(x,theta)
{
  x <- x[,1:3]
  return((abs(4*x[,1]-2) + theta[1]) /(1+theta[1]) +
           (abs(4*x[,2]-2) + theta[2]) /(1+theta[2]) + (abs(4*x[,3]-2) + theta[3]) /(1+theta[3]))
}



#' Non Linear Simulator with Calibration Parameter in any Dimension
#'
#' @param x Matrix of field design in dimension n \eqn{\times} p.
#' @param theta Calibration parameter of dimension 3
#' @param beta binary vector. Same size as \code{theta} indicating
#'  with 1 which dimensions are active
#'
#' @return vector of responses corresponding to a linearized Sobol function
#' @export
#'
#' @examples
simGal <- function(x,theta,beta=rep(1,ncol(x)))
{
  v <- matrix(NA,nrow(x),length(beta))
  for (i in 1:length(beta))
    v[,i] <- (abs(4*x[,i]-2) + theta[i]) /(1+theta[i])
 return(apply(v[,beta==1],1,sum))
}


