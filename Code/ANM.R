require(reticulate)
source_python("GauProRegPy.py")

# Gaussian process regressor using reticulate and GauProRegPy
GaussianProcPredict = function(x, y){
  
  # call gpr from Python
  pred_df = gpr(x, y)
  
  # return predicted y
  return(pred_df$ypr)
}


rbf_dot = function(x, deg){
  
  # Set kernel size to median distance between points, if no kernel specified
  m = length(x)
  G = x*x
  Q = matrix(rep(G, each = m), nrow = m, ncol = m, byrow = T)
  H = Q + t(Q) - 2 * x %*% t(x)
  
  if(deg == -1){
    Hlt = as.vector(t(H))
    Hlt[which(!as.vector(t(lower.tri(H))))] = 0
    dists = as.vector(t(H)) - Hlt
    deg = sqrt(0.5*median(dists[dists>0]))
  }
  
  H = exp(-H / 2 / (deg**2))
  
  return(H)
}


FastHsicTest = function(x, y, sig=c(-1,-1), maxpts = 200){
  
  m = length(x)
  
  if(m > maxpts){
    
    indx = as.integer(floor(seq(0, m, ((m-1)/(maxpts-1)))))
    
    Xm = as.numeric(x[indx])
    Ym = as.numeric(y[indx])
  } else {
    
    Xm = as.numeric(x)
    Ym = as.numeric(y)  
  
  }
  
  m = length(Xm)
  
  id = matrix(0, nrow=m, ncol=m)
  for(i in 1:m)
    id[i,i] = 1
  
  H = id - 1 / m * matrix(1, nrow=m, ncol=m) 
  K = rbf_dot(Xm, sig[1])
  L = rbf_dot(Ym, sig[2])
  
  Kc = H %*% (K %*% H)
  Lc = H %*% (L %*% H)
  
  teststat = (1/m) * sum(t(Kc) * Lc)
  if(is.infinite(teststat) || is.na(teststat)){
    teststat = 0
  }
  
  return(teststat)
  
}

sdN = function(x){
  z = sd(x) * sqrt((length(x)-1)/length(x))
  return(z)
}

anm_score = function(x, y) {
  
  m = length(x)
  
  y = (y - mean(y))/sdN(y)
  x = (x - mean(x))/sdN(x)
  
  gp = GaussianProcPredict(x, y)
  yhat = gp - y
  
  yhat = (yhat - mean(yhat))/sdN(yhat)
  x = (x - mean(x))/sdN(x)
  
  stat = FastHsicTest(yhat, x)
  
  return(stat)
}