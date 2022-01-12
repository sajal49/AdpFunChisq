# If a method returned scores that are smaller the better:
# a negative sign is multiplied to make them larger the better.

# If a method returned scores that are larger the better:
# no changes are made to that method

# Causal Inference by Stochastic Complexity -- Originally smaller the better
CISC_M = function(tbl){
  
  t = Untable(tbl)
  stats = -Cisc(as.numeric(t[,1]), as.numeric(t[,2]))$x.to.y
  
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}


FI_M = function(tbl){
  
  estimate = fun.chisq.test(tbl)$estimate
  
  if(is.na(estimate)){
    estimate = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(estimate))
  
}

GKT_M = function(tbl){
  
  # Wei2017.et.al showed that Function index is equivalent to GKT
  estimate = fun.chisq.test(tbl)$estimate
  
  if(is.na(estimate)){
    estimate = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(estimate))
  
}

# Conditional Entropy -- Originally smaller the better
CE_M = function(tbl){
  
  untbl = Untable(tbl)
  stats = -condentropy( as.numeric(untbl[,2]), as.numeric(untbl[,1]) )
  
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}

# Fraction of Information -- Originally larger the better
FOI_M = function(tbl){
  
  untbl = Untable(tbl)
  H_y = entropy::entropy( as.numeric(untbl[,2]) )
  H_xy = condentropy( as.numeric(untbl[,2]), as.numeric(untbl[,1]) )
  fracinfo = (H_y - H_xy)/(H_y)
  
  if(is.nan(fracinfo)){
    fracinfo = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(fracinfo))
}

# Adapted Functional Chi-squared Test -- Originally smaller the better
AFC_M = function(tbl){
  
  pvalue = -AdpFunChisq(tbl, log.p = TRUE)$p.value
  
  if(is.na(pvalue)){
    pvalue = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(pvalue))
  
}

# Distance Correlation -- Originally smaller the better
DC_M = function(tbl){
  
  t = Untable(tbl)
  stats = -dc(as.numeric(t[,1]), as.numeric(t[,2]))$x.to.y
  
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}

# Digital Regression -- Originally larger the better
DR_M = function(tbl){
  
  t = Untable(tbl)
  stats = dr(as.numeric(t[,1]), as.numeric(t[,2]), level = -3)$x.to.y
  
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}

# Hidden Compact Representation -- Originally larger the better
HCR_M = function(tbl){
  
  t = Untable(tbl)
  stats = HCR.fast(X = as.numeric(t[,1]), Y = as.numeric(t[,2]))$score
  
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}

# IGCI -- Originally smaller the better
IGCI_M = function(tbl, mode="table"){
  
  if(mode == "table"){
    tbl = Untable(tbl)
    x = jitter(as.numeric(tbl[,1]), amount=0.1)
    y = jitter(as.numeric(tbl[,2]), amount=0.1)
  } else {
    x = as.numeric(tbl[,1])
    y = as.numeric(tbl[,2])
  }
  
  stats = -igci(x = x, y = y, refMeasure = 2, 
                estimator = 3)$C.x.to.y
  
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}

# Additive noise model -- Originally smaller the better

ANM_M = function(tbl, mode="table"){
 
  if(mode == "table"){
    tbl = Untable(tbl)
    x = jitter(as.numeric(tbl[,1]), amount=0.1)
    y = jitter(as.numeric(tbl[,2]), amount=0.1)
  } else {
    x = as.numeric(tbl[,1])
    y = as.numeric(tbl[,2])
  }
  
  stats = -anm_score(x=x, y=y)
  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
}

# 1-NN -- 1 nearest neigbhor regression, higher the better
KNN_M = function(tbl, mode="table"){
  
  if(mode == "table"){
    tbl = Untable(tbl)
    x = jitter(as.numeric(tbl[,1]), amount=0.1)
    y = jitter(as.numeric(tbl[,2]), amount=0.1)
  } else {
    x = as.numeric(tbl[,1])
    y = as.numeric(tbl[,2])
  }
  
  stats = knn.reg(train = x, y = y, k = 1)$R2Pred

  if(is.na(stats)){
    stats = -.Machine$integer.max # assign smallest possible value
  }
  
  return(as.numeric(stats))
  
}

# SCR -- Subcopula based regression, higher the better
SCR_M = function(tbl){
  
  # find row and column sums of tbl
  rsm = rowSums(tbl)
  csm = colSums(tbl)
  
  # add a very small value to tables to avoid 0 rows / columns
  if(any(rsm == 0) || any(csm == 0)){
    tbl = tbl + 0.1
  }
  
  # find total sum
  n = sum(tbl)
    
  # find pi and pj
  pi = rowSums(tbl)/n
  pj = colSums(tbl)/n
    
  # find cumulative marginal distributions
  # ui and vj
  ui = cumsum(pi)
  vj = cumsum(pj)
    
  # observed joint distribution
  tbl = tbl/n
    
  # construct subcopulas
  txy = tbl / pi
  tyx = t(t(tbl)/pj)
    
  # calculate r(v|u)
  rvu = vj %*% t(txy)
    
  # calculate r(u|v)
  ruv = ruv = t(ui) %*% tyx
    
  # calculate pxy**2
  vjpj = sum(vj*pj)
  pxy2 = sum(((rvu - vjpj)**2)*pi)/sum(((vj - vjpj)**2)*pj)
    
  # calculate pyx**2
  uipi = sum(ui*pi)
  pyx2 = sum(((ruv - uipi)**2)*pj) / sum(((ui - uipi)**2)*pi)
  
  if(is.na(pxy2)){
    pxy2 = -.Machine$integer.max # assign smallest possible value
  }
    
  return(pxy2)
}
