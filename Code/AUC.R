# AUC.R -- Compute the area under the curve of function y=f(x)

trapezint <- function (x, y, a, b) {
  ## This function is copied exactly from the ROC package version 1.32.0
  ## It assumes that x is sorted.

  if (length(x) != length(y))
    stop("length x must equal length y")
  y <- y[x >= a & x <= b]
  x <- x[x >= a & x <= b]
  if (length(unique(x)) < 2)
    return(NA)
  ya <- approx(x, y, a, ties = max, rule = 2)$y
  yb <- approx(x, y, b, ties = max, rule = 2)$y
  x <- c(a, x, b)
  y <- c(ya, y, yb)
  h <- diff(x)
  lx <- length(x)
  0.5 * sum(h * (y[-1] + y[-lx]))
}

AUC <- function(FPR, TPR)
{ # This function can only be used to calculate AUROC, not AUPR
  if(length(unique(FPR)) < 2) {
    FPR <- c(FPR, 0, 1)
    TPR <- c(TPR, 0, 1)
  }
  o <- order(FPR, TPR)
  trapezint(FPR[o], TPR[o], 0, 1)
}

cal.AUC.PR <- function(recall, precision)
{
  if(length(unique(recall)) < 2) {
    return(0)
  }
  o <- order(recall, -precision)
  trapezint(recall[o], precision[o], 0, 1)
}

test.AUC <- function()
{
  examples <- list()

  examples[[1]] <- list(FPR=c(0, 0.5, 1), TPR=c(0, 0.5, 1), AUC=0.5)
  examples[[2]] <- list(FPR=c(0.5, 1, 0), TPR=c(0.5, 1, 0), AUC=0.5)
  examples[[3]] <- list(FPR=seq(0, 1, 0.1), TPR=seq(0, 1, 0.1), AUC=0.5)
  examples[[4]] <- list(FPR=c(0, 0, 1), TPR=c(0, 1, 1), AUC=1)
  examples[[5]] <- list(FPR=c(0, 1, 1), TPR=c(0, 0, 1), AUC=0)

  for(i in 1:length(examples) ) {
    ex <- examples[[i]]
    if(AUC(ex$FPR, ex$TPR) != ex$AUC) {
      cat("Testing function AUC (area under curve) ...\n")
      cat("Example", i)
      stop(" failed!\n")
    } else {
      # cat(" passed.\n")
    }
  }

  recall <- c(0.0, 0.5, 1.0, 1.0, 1.0, 1.0)
  precision <- c(1, 1, 1, 2/3, 0.5, 0.4)
  if(cal.AUC.PR(recall, precision) != 1) {
    cat("Testing function AUC (area under curve) ...\n")
    stop("Precision-recall example failed!\n")
  }
}

test.AUC()

for(pkg in c("ROCR")) {
  if(! pkg %in% installed.packages()) {
    install.packages(pkg, method="curl")
  }
  require(pkg, character.only = TRUE)
}

calculate.ROC.PR <- function(stats, true.classes)
{
  pred <- prediction(stats, true.classes)
  
  #ROC
  perf.ROC <- performance(pred, measure = "tpr", x.measure = "fpr")
  FPR <- perf.ROC@x.values[[1]]
  TPR <- perf.ROC@y.values[[1]]
  
  if(any(is.na(FPR))) browser()
  if(any(is.na(TPR))) browser()
  
  auc.ROC <- AUC(FPR, TPR)
  
  #PR
  perf.PR <- performance(pred, measure = "prec", x.measure = "rec")
  recall <- perf.PR@x.values[[1]]
  precision <- perf.PR@y.values[[1]]
  if(is.na(precision[1])) {
    precision[1] <- precision[2]
  }
  
  if(any(is.na(recall))) browser()
  if(any(is.na(precision))) browser()
  
  auc.PR <- cal.AUC.PR(recall, precision)
  
  if(0) {
    xvalues <- perf@x.values[[1]]
    
    if(min(xvalues) > 0) browser()
    
    xvalues[is.na(xvalues)] <- 0
    yvalues <- perf@y.values[[1]]
    yvalues[is.na(yvalues)] <- 0
    auc.PR <- AUC(xvalues, yvalues)
  }

  return(list(perf.ROC=perf.ROC, auc.ROC=auc.ROC, 
              perf.PR=perf.PR, auc.PR=auc.PR,
              FPR=FPR, TPR=TPR, 
              Precision=precision, Recall=recall))  
}

plot.ROC.with.AUC <- function(stats, true.classes, plot=FALSE, ...)
{
  # larger stats value are considered to correspond to larger true.classes labels
  # When using pvalues, if true is represented by true.classes = 1,
  #   and false by true.classes = 0, stats must be negative p-value
  # plot.ROC.with.AUC(stats = - p.val, true.classes, ... )

  res <- calculate.ROC.PR(stats, true.classes)
  within(res,
       {
         if(plot) {
           plot(perf.ROC, ylim=c(0,1), xlim=c(0,1), col="blue", ...)
           legend("bottomright", paste("AUROC =", format(auc.ROC, digits=4)))
           
           plot(perf.PR, ylim=c(0,1), xlim=c(0,1), col="chocolate", ...)
           #mtext(paste("AUPR =", format(auc, digits=5)), 1, -1, adj=1)
           legend("topright", paste("AUPR =", format(auc.PR, digits=4)))
         }
       }
  )

  return(c(AUROC=res$auc.ROC, AUPR=res$auc.PR))
}

test.plot.ROC.with.AUC <- function()
{
  stats <- c(5,1,2,3,4)
  true.classes <- c(0,0,0,1,1)

  plot.ROC.with.AUC(stats, true.classes, plot=TRUE)

  stats <- c(-5,1,2,3,4)
  true.classes <- c(0,0,0,1,1)

  plot.ROC.with.AUC(stats, true.classes, plot=TRUE)
}
