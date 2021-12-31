# AUC-multiple.R -- plot multiple ROC and PR curves in one plot

require(ggplot2)

ggplot.ROC.PR.curves <- function(
  list.stats, true.classes, colr, ltyps, plot=FALSE, title=NULL, digits=2, ...
)
{
  # larger stats value are considered to correspond to larger true.classes labels
  # When using pvalues, if true is represented by true.classes = 1,
  #   and false by true.classes = 0, stats must be negative p-value
  # ggplot.ROC.PR.curves(stats = - p.val, true.classes, ... )
  
  l <- length(list.stats)
  res <- vector("list", l)
  
  auc.PRs <- vector("double", l)
  auc.ROCs <- auc.PRs
  
  for(i in seq(l)) {
    
    res[[i]] <- calculate.ROC.PR(list.stats[[i]], true.classes)
    
    auc.PRs[i] <- res[[i]]$auc.PR
    auc.ROCs[i] <- res[[i]]$auc.ROC
  }
  
  if(plot) {
    
    ROC <- NULL
    
    for(i in seq(l)) {
      
      d <- data.frame(
        FPR = res[[i]]$FPR,
        TPR = res[[i]]$TPR,
        Method=paste0(names(list.stats)[i], ' (',
                      format(auc.ROCs[i], digits=digits), ')')
      )
      
      ROC <- rbind(ROC, d)
      
    }
    
    colra = c()
    meth_tb = table(ROC$Method)
    
    for(i in 1:length(colr)){
      
      colra = c(colra, rep(colr[i], meth_tb[i]))
      
    }
    
    p <- ggplot(ROC, aes(x=FPR, y=TPR, group=Method, 
                         color=Method, linetype=Method)) +
      labs(x="False positive rate", y="True positive rate", 
           title = title, color="Method (AUROC)", 
           linetype="Method (AUROC)") + 
      geom_line(size=1) + xlim(0, 1) + ylim(0, 1) +
      scale_color_manual(values=colr) +
      scale_linetype_manual(values=ltyps)+
      theme_bw(base_size = 13) +
      theme(plot.background = element_rect(fill = "transparent"), 
            axis.title = element_text(face = "bold", size = 20),
            legend.title = element_text(face = "bold", size = 16),
            title = element_text(face = "bold", size = 18),
            legend.text = element_text(size = 16))
    
    print(p)
    
    PRC <- NULL
    
    for(i in seq(l)) {
      
      d <- data.frame(
        Recall=res[[i]]$Recall,
        Precision=res[[i]]$Precision,
        Method=paste0(names(list.stats)[i], ' (',
                      format(auc.PRs[i], digits=digits), ')')
      )
      
      PRC <- rbind(PRC, d)
      
    }
    
    p <- ggplot(PRC, aes(x=Recall, y=Precision, 
                         group=Method, color=Method, 
                         linetype=Method)) +
      labs(title = title, color="Method (AUPR)", 
           linetype="Method (AUPR)") +
      scale_color_manual(values=colr) +
      geom_line(size=1) + xlim(0, 1) + ylim(0, 1) +
      scale_linetype_manual(values=ltyps)+
      theme_bw(base_size = 13) +
      theme(plot.background = element_rect(fill = "transparent"), 
            axis.title = element_text(face = "bold", size = 20),
            legend.title = element_text(face = "bold", size = 16),
            title = element_text(face = "bold", size = 18),
            legend.text = element_text(size = 16))
    
    print(p)
    
  }
  
  return(list(AUROC=auc.ROCs, AUPR=auc.PRs, 
              res=res))
}

plot.multiple.ROC.PR.curves <- function(
  list.stats, true.classes, plot=FALSE, ...
)
{
  # larger stats value are considered to correspond to larger true.classes labels
  # When using pvalues, if true is represented by true.classes = 1,
  #   and false by true.classes = 0, stats must be negative p-value
  # plot.ROC.with.AUC(stats = - p.val, true.classes, ... )
  
  l <- length(list.stats)
  res <- vector("list", l)
  
  auc.PRs <- vector("double", l)
  auc.ROCs <- auc.PRs
  
  for(i in seq(l)) {
    
    res[[i]] <- calculate.ROC.PR(list.stats[[i]], true.classes)
    
    auc.PRs[i] <- res[[i]]$auc.PR
    auc.ROCs[i] <- res[[i]]$auc.ROC
  }
  
  if(plot) {
    
    for(i in seq(l)) {
      within(res[[i]],
             {
               if(i == 1) {
                 plot(perf.ROC, ylim=c(0,1), xlim=c(0,1), col=i, lwd = l-i+1, ...)
               } else {
                 plot(perf.ROC, col=i, add=TRUE, lwd = l-i+1, ...)
               }
             }
      )
    }    
    
    legend("bottomright", bty="n", 
           lwd = rev(seq(l)), col=seq(l), 
           paste(names(list.stats), "AUROC =", 
                 format(auc.ROCs, digits=4)))
    
    for(i in seq_along(res)) {
      within(res[[i]],
             {
               if(i == 1) {
                 plot(perf.PR, ylim=c(0,1), xlim=c(0,1), col=i, 
                      lwd = l - i + 1, ...)
               } else {
                 plot(perf.PR, col=i, lwd = l - i + 1, add=TRUE, ...)
               }
             }
      )
    }
    #mtext(paste("AUPR =", format(auc, digits=5)), 1, -1, adj=1)
    legend("topright", bty = "n", lwd = rev(seq(l)), 
           col=seq(l), 
           paste(names(list.stats), "AUPR =", 
                 format(auc.PRs, digits=4)))
    
  }
  
  return(c(AUROC=auc.ROCs, AUPR=auc.PRs))
}

test.plot.multiple.ROC.PR.curves <- function()
{
  list.stats <- list(c(5,1,2,3,4), c(-5,1,2,3,4), c(1, -1, 7, 6, 4), c(0, 0, 1, 2, 2))
  names(list.stats) <- paste("Method", seq_along(list.stats))
  
  true.classes <- c(0,0,0,1,1)
  
  plot.multiple.ROC.PR.curves(list.stats, true.classes, plot=TRUE)
}
