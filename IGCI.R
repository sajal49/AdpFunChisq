#################################################################
# Author : Joe Song                                             #
# NMSU Song Lab, Department of Computer Science,                #
# New Mexico State University                                   #
# Migrated from "igci.m"                                        #
#################################################################

improved_slope_estimator <- function(x,y) 
{
  # see eq. (22) in http://arxiv.org/abs/1412.3773v2
  m = length(x)
  
  # [sx,ind] = sort(x)
  ind <- order(x)
  sx <- x[ind]
  
  Xs = c() 
  Ys = c() 
  Ns = c()
  last_index = 1
  
  for (i in 2:m) {
    if (x[ind[i]] != x[ind[last_index]]) {
      Xs = c(Xs, x[ind[last_index]])
      Ys = c(Ys, y[ind[last_index]])
      Ns = c(Ns, i - last_index)
      last_index = i
    }
  }
  
  Xs = c(Xs, x[ind[last_index]])
  Ys = c(Ys, y[ind[last_index]])
  Ns = c(Ns, m + 1 - last_index)
  
  m = length(Xs)
  a = 0
  Z = 0
  for (i in 1:(m-1)) {
    X1 = Xs[i]
    X2 = Xs[i+1]
    Y1 = Ys[i]  
    Y2 = Ys[i+1]
    
    if ((X2 != X1) && (Y2 != Y1)) {
      a = a + log(abs((Y2 - Y1) / (X2 - X1))) * Ns[i]
      Z = Z + Ns[i]
    }
  }
  a = a / Z
  return(a)
}

igci <- function(x,y,refMeasure,estimator)
{
  # Performs causal inference in a deterministic scenario (see [1] and [2] for details)
  # 
  # USAGE:
  #   f = igci(x,y,refMeasure,estimator)
  # 
  # INPUT:
  #   x          - m x 1 observations of x
  #   y          - m x 1 observations of y
  #   refMeasure - reference measure to use:
  #                  1: uniform
  #                  2: Gaussian
  #   estimator -  estimator to use:
  #                  1: entropy (eq. (12) in [1]),
  #                  2: integral approximation (eq. (13) in [1]).
  #                  3: new integral approximation (eq. (22) in [2]) that should
  #                     deal better with repeated values
  # 
  # OUTPUT: 
  #   f < 0:       the method prefers the causal direction x -> y
  #   f > 0:       the method prefers the causal direction y -> x
  # 
  # EXAMPLE: 
  #   x = randn(100,1); y = exp(x); igci(x,y,2,1) < 0
  #
  #
  # Copyright (c) 2010-2015  Povilas Daniušis, Joris Mooij
  # All rights reserved.  See the file LICENSE for license terms.
  # ----------------------------------------------------------------------------
  #
  # [1]  P. Daniusis, D. Janzing, J. M. Mooij, J. Zscheischler, B. Steudel,
  #      K. Zhang, B. Schoelkopf:  Inferring deterministic causal relations.
  #      Proceedings of the 26th Annual Conference on Uncertainty in Artificial 
  #      Intelligence (UAI-2010).  
  #      http://event.cwi.nl/uai2010/papers/UAI2010_0121.pdf
  # [2]  J. M. Mooij, J. Peters, D. Janzing, J. Zscheischler, B. Schoelkopf
  #      Distinguishing cause from effect using observational data: methods and benchmarks
  #      arXiv:1412.3773v2, submitted to Journal of Machine Learning Research
  #      http://arxiv.org/abs/1412.3773v2
  
  # ignore complex parts
  x <- Re(x)
  y <- Re(y)
  
  m <- length(x)
  
  # check input arguments
  if (! is.vector(x)) {
    stop('Dimensionality of x must be 1')
  }
  if (length(x) < 2) {
    stop('Not enough observations in x (must be >= 2)')
  }
  
  if (! is.vector(y)) {
    stop('Dimensionality of y must be 1')
  }
  if (length(x) < 2) {
    stop('Not enough observations in y (must be >= 2)')
  }
  
  if (length(x) != length(y)) {
    stop('Lenghts of x and y must be equal')
  }
  
  if(refMeasure == 1) {
    # uniform reference measure
    x = (x - min(x)) / (max(x) - min(x))
    y = (y - min(y)) / (max(y) - min(y))
  } else if(refMeasure == 2) {
    # Gaussian reference measure
    x = (x - mean(x)) / sd(x)
    y = (y - mean(y)) / sd(y)
  } else {
    warning('Warning: unknown reference measure - no scaling applied')
  }
  
  if(estimator == 1) {
    # difference of entropies
    
    # [x1,indXs] = sort(x);
    indXs <- order(x)
    x1 <- x[indXs]
    
    #[y1,indYs] = sort(y);
    indYs <- order(y)
    y1 <- y[indYs]
    
    n1 = length(x1)
    hx = 0.0
    for (i in 1:(n1-1)) {
      delta <- x1[i+1]-x1[i]
      if (delta) {
        hx <- hx + log(abs(delta))
      }
    }
    # hx = hx / (n1 - 1) + psi(n1) - psi(1)
    hx = hx / (n1 - 1) + digamma(n1) - digamma(1)
    
    n2 = length(y1)
    hy = 0.0
    for (i in 1:(n2-1)) {
      delta = y1[i+1]-y1[i]
      if (delta) {
        hy = hy + log(abs(delta))
      }
    }
    
    # hy = hy / (n2 - 1) + psi(n2) - psi(1)
    hy = hy / (n2 - 1) + digamma(n2) - digamma(1)
    
    f = hy - hx
    C.x.to.y <- hy - hx
    C.y.to.x <- hx - hy
    
  } else if(estimator == 2) {
    # integral-approximation based estimator
    a = 0
    b = 0
    
    #[sx,ind1] = sort(x);
    ind1 <- order(x)
    sx <- x[ind1]
    
    #[sy,ind2] = sort(y);
    ind2 <- order(y)
    sy <- y[ind2]
    
    for (i in 1:(m-1)) {
      X1 = x[ind1[i]]
      X2 = x[ind1[i+1]]
      Y1 = y[ind1[i]]
      Y2 = y[ind1[i+1]]
      if ( (X2 != X1) && (Y2 != Y1) ) {
        a = a + log(abs((Y2 - Y1) / (X2 - X1)))
      }
      X1 = x[ind2[i]]  
      X2 = x[ind2[i+1]]
      Y1 = y[ind2[i]]
      Y2 = y[ind2[i+1]]
      if ((Y2 != Y1) && (X2 != X1)) {
        b = b + log(abs((X2 - X1) / (Y2 - Y1)))
      }
    }
    
    f = (a - b) / m
    C.x.to.y <- a / (m-1)
    C.y.to.x <- b / (m-1)
    
  } else if(estimator == 3) {
    # integral-approximation based estimator
    # improved handling of values that occur multiple times
    # f = (improved_slope_estimator(x,y) - improved_slope_estimator(y,x)) / m
    C.x.to.y <- improved_slope_estimator(x,y) / (m-1)
    C.y.to.x <- improved_slope_estimator(y,x) / (m-1)
    f <- (C.x.to.y - C.y.to.x) * (m-1) / m
  } else {
    stop('Unknown estimator')
  }
  return(list(f=f, C.x.to.y=C.x.to.y, C.y.to.x=C.y.to.x))
}

test_igci <- function()
{
  examples <- list()
  
  examples[[1]] <- list(
    tag = "X -> Y: toy",
    x = c(0, -5, 3, 2, 1.3, 4.5),
    y = c(8, -10, 1, 0, 5, -6)
    
    # Truth (from octave):
    #1, 1, f=0.109394
    #1, 2, f=0.155888
    #1, 3, f=0.031178
    #2, 1, f=0.037252
    #2, 2, f=0.035652
    #2, 3, f=0.007130
  )
 
  N <- 100
  xr <- rnorm(N)
  examples[[2]] <- list(
    tag = "X -> Y: linear",
    x = xr,
    y = xr * 2 + rnorm(N) * .05
  )
  
  examples[[3]] <- list(
    tag = "X -> Y: quadratic",
    x = xr,
    y = (xr ^ 2 + rnorm(N) * .10)
  )

  xr <- c(-0.0473377, 1.4006309, -0.9017835, -0.0054639, 0.4842724, 
          0.7158695,  -0.9371985,  0.1532856,  -0.7102305, -0.3690689)

  examples[[4]] <- list(
    tag = "X -> Y: exponential",
    x = xr,
    y = exp(xr)
    
    # Truth from Octave:
    # 1, 1, f=-0.495513
    # 1, 2, f=-0.891924
    # 1, 3, f=-0.099103
    # 2, 1, f=-0.440130
    # 2, 2, f=-0.792235
    # 2, 3, f=-0.088026
  )
  
  for(ex in examples) {
    within(ex, 
           {
             cat(tag, ":\n")
             plot(x, y, main=tag)
             for(refMeasure in 1:2) {
               for(estimator in 1:3) {
                 cat("    ", c("uniform", "Gaussian")[refMeasure], "+", 
                     c("Diff Entropy", "Slope", "Slope Improved")[estimator], ":\n")
                 
                 within(igci(x, y, refMeasure, estimator),
                        {
                          cat("        f(X->Y) =", f, 
                              "C(X->Y)=", C.x.to.y, "C(Y->X)=", C.y.to.x, 
                              "\n")
                        } )
               }
             }
           })    
  }
}

# test_igci()
