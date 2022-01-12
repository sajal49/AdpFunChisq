require(FunChisq)

# Parameters:

# x: (table, matrix or data.frame) a contingency table with 
# 'cause' on the row variable and 'effect' on the column variable.

# log.p: (logical) if TRUE, logged p-value will be returned.

AdpFunChisq = function(x, log.p=FALSE){
  
  if(!(any(class(x) %in% c("matrix", "data.frame", "table")))){
    stop("Adaptive functional chi-squared test only works with a 2D matrix, table or data.frame")
  }
  
  DNAME = deparse(substitute(x))
  
  # if x is portrait
  if(ncol(x) <= nrow(x)){
    stats = fun.chisq.test(x, log.p = log.p)
    hobj = structure(list(statistic = stats$statistic, p.value = stats$p.value,
                          data.name = DNAME, parameters = stats$parameter,
                          method = "Adapted Functional Chi-Squared Test"),
                     class = "htest")
  } else { # if x is landscape
    col_sum = colSums(x) # get column sums
    
    # store original dimensions
    og_dims = dim(x)
    
    # remove s-r most sparsely populated columns
    rem_indices = order(col_sum)[1:(og_dims[2] - og_dims[1])]
    F_x = x[,-rem_indices]
    
    # obtain residue
    residue = sum(col_sum[rem_indices])
    
    # if residue is 0
    if(residue == 0){
      stats = fun.chisq.test(F_x, log.p = log.p) # only compute FunChisq of F_x
      df = stats$parameter +  ((og_dims[2] - og_dims[1]) * (ncol(F_x)-1)) # adjusted degrees of freedom
      p.value = pchisq(q = stats$statistic, df = df, lower.tail = FALSE, log.p = log.p) # compute pvalue with new df
      hobj = structure(list(statistic = stats$statistic, p.value = p.value,
                            data.name = DNAME, parameters = df,
                            method = "Adapted Functional Chi-Squared Test"),
                       class = "htest")
    } else {
      E_x = matrix(0, nrow=(og_dims[2] - og_dims[1] + 1), ncol = ncol(F_x))
      
      # redistribute residue among s-r+1 rows with equal likelihood
      rrsum = rmultinom(n=1, size=residue, prob = rep(1/nrow(E_x), nrow(E_x)))
      
      # construct random table E_x
      for(i in 1:nrow(E_x)){
        E_x[i, ] = t(rmultinom(n=1, size=rrsum[i], prob = rep(1/ncol(E_x), ncol(E_x))))
      }
      
      # compute two FunChisq
      Chi_F_x = fun.chisq.test(F_x, log.p = TRUE)
      Chi_E_x = fun.chisq.test(E_x, log.p = TRUE)
      
      stats = Chi_F_x$statistic + Chi_E_x$statistic # add both statistics
      df = Chi_F_x$parameter + Chi_E_x$parameter # add both degrees of freedom
      p.value = pchisq(q = stats, df = df, lower.tail = FALSE, log.p = log.p) # compute pvalue
      hobj = structure(list(statistic = stats, p.value = p.value,
                            data.name = DNAME, parameters = df,
                            method = "Adapted Functional Chi-Squared Test"),
                       class = "htest")
    }
  }
  
  return(hobj)
}