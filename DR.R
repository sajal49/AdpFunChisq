# The original codes for DC was obtained from:
# http://eda.mmci.uni-saarland.de/prj/cisc/

require(Rfast)
# for testing
require(assertthat)

unique.matlab <- function(Z){
  unz <- unique(Z)
  foccur <- sapply(unz, function(i){ min(which(Z == i)) })
  Z_new = match(Z, unz)
  return(list(unz, foccur, Z_new))
}

hist3 <- function(X, Y, X_values, Y_values){
  X.factor <- factor(X, levels = X_values)
  Y.factor <- factor(Y, levels = Y_values)
  p <- table(X.factor, Y.factor)
  return(p)
}

chi_sq_quant <- function(x, y, num_states_x, num_states_y){
  
  if(num_states_x==1||num_states_y==1){
    result <- 1
    TT <- 0
  }else{
    stats = Rfast::gchi2Test(x=x, y=y, logged = TRUE)[1,]
    TT = as.numeric(stats[1])
    result = as.numeric(stats[2])
  }
  return(c(result, TT))
}

fit_discrete <- function(X, Y, level = -3){
  
  ##########
  # parameter
  ##########
  num_iter <- 10
  num_pos_fct <- min(c(max(Y)-min(Y),20))
  
  fct <- NULL
  p_val <- NULL
  
  # rescaling:
  # X_new takes values from 1...X_new_max
  # Y_values are everything between Y_min and Y_max
  X.unique.res <- unique.matlab(X)
  X_values <- X.unique.res[[1]]
  aa <- X.unique.res[[2]]
  X_new <- X.unique.res[[3]]
  
  Y_values <- min(Y):max(Y)
  
  if(length(X_values) == 1 || length(Y_values) == 1){
    fct <- 1 * Y_values[1]
    p_val <- 1
  }else{
    p <- hist3(X, Y, X_values, Y_values)
    
    fct <- c()
    cand <- list()
    
    for (i in  1:length(X_values)) {
      a <- as.numeric(sort(p[i,]))
      b <- as.numeric(order(p[i,]))
      
      for (k in 1:ncol(p)) {
        if (k != b[length(b)]) {
          p[i,k] <- p[i,k]+1/(2 * abs(k - b[length(b)]))
        }else{
          p[i,k] <- p[i,k] + 1
        }
      }
      
      a <- as.numeric(sort(p[i,]))
      b <- as.numeric(order(p[i,]))
      
      cand[[i]] <- b
      fct <- c(fct, Y_values[b[length(b)]])
    }
    
    yhat <- fct[X_new]
    eps <- Y - yhat
    
    if(length(unique(eps))==1){
      #message("Warning!! there is a deterministic relation between X and Y.")
      p_val <- 0;
    }else{
      p_val <- chi_sq_quant(eps, X, length(unique(eps)),length(X_values))[1]
    }
    
    i <- 0
    pos_fct <- vector("list", num_pos_fct + 1)
    while (p_val < level && i < num_iter) {
      perx = sample(x = length(X_values), size = length(X_values), replace = FALSE)
      for (j_new in perx) {
        p_val_comp <- rep(NA, num_pos_fct + 1)
        p_val_comp2 <- rep(NA, num_pos_fct + 1)
        for (j in c(1:(num_pos_fct + 1))) {
          pos_fct[[j]] <- fct
          pos_fct[[j]][j_new] <- Y_values[cand[[j_new]][length(cand[[j_new]])-(j-1)]]
          yhat <- pos_fct[[j]][X_new]
          eps <- Y - yhat
          res.tmp <- chi_sq_quant(eps, X, length(unique(eps)),length(X_values))
          p_val_comp[j] <- res.tmp[1]
          p_val_comp2[j] <- res.tmp[2]
        }
        
        aa <- max(p_val_comp)
        j_max <- which.max(p_val_comp)
        if(aa < -6.907755) { # log of 0.001 
          aa <- min(p_val_comp2)
          j_max <- which.min(p_val_comp2)
        }
        
        fct <- pos_fct[[j_max]]
        yhat <- fct[X_new]
        eps <- Y - yhat
        p_val <- chi_sq_quant(eps, X, length(unique(eps)), length(X_values))[1]
      }
      i <- i + 1
    }
    fct <- fct + round(mean(eps))
  }
  return(p_val)
}

dr <- function(X, Y, level = -3){
  x.to.y <- fit_discrete(X, Y, level)
  y.to.x <- fit_discrete(Y, X, level)
  return(list(x.to.y = x.to.y,
              y.to.x = y.to.x))
}
test_dr <- function(){
  # Test 1
  X <- c(1, 1, 1, 3, 1, 1, 1, 1, 3, 3)
  Y <- c(3, 2, 3, 6, 3, 3, 3, 2, 5, 4)
  dr_res <- dr(X, Y, level = -3)
  assert_that(are_equal(exp(dr_res$x.to.y), 0.23965103644177588))
  assert_that(are_equal(exp(dr_res$y.to.x), 1.0))
  
  # Test 2
  X <- c(2, 2, 1, 4, 2, 2, 4, 4, 1, 2)
  Y <- c(4, 5, 2, 5, 4, 4, 5, 5, 3, 3)
  dr_res <- dr(X, Y, level = -3)
  assert_that(are_equal(exp(dr_res$x.to.y), 0.5459439747510072))
  assert_that(are_equal(exp(dr_res$y.to.x), 0.4302277231528573))
  
  # Test 3
  X <- c(1, 4, 1, 1, 3, 1, 2, 1, 2, 2)
  Y <- c(3, 7, 4, 3, 6, 3, 3, 4, 3, 3)
  dr_res <- dr(X, Y, level = -3)
  assert_that(are_equal(exp(dr_res$x.to.y), 0.4752910833430206))
  assert_that(are_equal(exp(dr_res$y.to.x), 0.4141793659679206))
  
  X = c(1,1,1,1, 2,2,2,2, 3,3,3,3)
  Y = c(1,1,1,1, 2,2,2,2, 1,1,1,1)
  dr_res = dr(X, Y, level=-3)
  assert_that(are_equal(exp(dr_res$x.to.y), 1.0))
  assert_that(are_equal(exp(dr_res$y.to.x), 0.08326451666355039))
  
  X = c(1,2,3, 3,2,1, 1,1,1, 2,2,2, 3,3,3)
  Y = c(1,1,1, 2,2,2, 1,2,3, 3,3,3, 3,2,1)
  dr_res = dr(X, Y, level=-3)
  assert_that(are_equal(exp(dr_res$x.to.y), 0.7191194024209784))
  assert_that(are_equal(exp(dr_res$y.to.x), 0.6380726074309759))
}