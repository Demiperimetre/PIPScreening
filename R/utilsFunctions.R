# list of distance by coordinates
tensordist = function(X1)
{
  lapply(1:ncol(X1), function(i) as.matrix(dist(X1[,i],"manhattan")))
}


