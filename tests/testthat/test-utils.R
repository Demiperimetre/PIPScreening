test_that("MCMC chains", {
  # test comme scenarii papier
  N=50
  beta=c(1,0,1,0,0)
  theta=4:8/10
  sd=0.05
  non.linear = c(2,1,1,1,1)
  x1 <- runif(N,0,1)
  x2 <- runif(N,0,1)
  x3 <- runif(N,0,1)
  x4 <- runif(N,0,1)
  x5 <- x3 + rnorm(N,sd=.2)
  p <- non.linear
  covargal <- cbind(x1^p[1],x2^p[2],x3^p[3],x4^p[4],x5^p[5]) # reality cov
  x <- cbind(x1,x2,x3,x4,x5)
  covarmod <- cbind(x1,x2,x3) # model covariates
  mu <- simGal(covargal,theta,beta)
  y <- mu + rnorm(N,0,sd) # field exp
  mod <- sim3(x, theta[1:3])
  xnorm <- (x - matrix(apply(x,2,min),nrow=nrow(x),ncol=ncol(x),byrow=T))/matrix(apply(x,2,max)-apply(x,2,min),nrow=nrow(x),ncol=ncol(x),byrow=T)
  Yexp <- y; Xexpnorm <- xnorm; Xexp <- x;  Rexp <- (y-mod)
  calibration1 <- list(computermodel=sim3,Yexp=Yexp,Xexp=Xexp,FALSE)
  calibration2 <- list(computermodel=sim3,Yexp=Yexp,Xexp=Xexp,TRUE)
  tdistFULL <- tensordist(xnorm)
  pgamma <-  5
  parwalkinit  <-  c(rep(.1,pgamma),.1,.1)
  init  <-  c(rep(0,pgamma),.004,.2)
  a <- TRUE
  cpt <- 0
  parprior <- rbind(matrix(1,nrow=pgamma,ncol=2),c(4,.02),c(3,1))
  nMWG <- 200
  nMet <- 200
  resmcmc <- (MCMC(nMWG,nMet,parwalkinit,init,Rexp,tdistFULL,1.9,parprior,TRUE,calibration1))
  expect_length(resmcmc$MH$chain[,1],nMet)
  expect_length(resmcmc$MH$chain[1,],pgamma+2)
})
