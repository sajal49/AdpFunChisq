# Processing and Evaluation on:
# Abalone, for:
# Overcoming biases in causal inference of molecular interactions

# Manuscript:
# "A Quantitative Comparison of Dystal and Backpropagation"
# Clark.et.al 1996, Australian Conference on Neural Networks (ACNN'96).
# Original data obtained from: https://archive.ics.uci.edu/ml/datasets/Abalone

# The data (abalone.data) is available in Data/Abalone/ directory.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
require(DescTools, quietly = TRUE)
require(infotheo, quietly = TRUE)
require(HCR, quietly = TRUE)
require(FunChisq, quietly = TRUE)
require(GridOnClusters, quietly = TRUE)
require(entropy, quietly = TRUE)
require(pbapply, quietly = TRUE)
require(doParallel, quietly = TRUE)

# source R files
source("AdpFunChisq.R")
source("AUC-multiple.R")
source("IGCI.R")
source("AUC.R")
source("Methods.R")
source("Utility.R")
source("Cisc.R")
source("DC.R")
source("ANM.R")
source("DR.R")

# Global parameters
colr = c(ggplot2::alpha("forestgreen",0.7), # CISC
         ggplot2::alpha("chartreuse2",0.7), # CE
         ggplot2::alpha("magenta",0.7), # FOI
         ggplot2::alpha("firebrick2",0.7), # AFC
         ggplot2::alpha("darkorange",0.7), # DC
         ggplot2::alpha("brown",0.7), # DR
         ggplot2::alpha("dodgerblue2",0.7), # HCR
         ggplot2::alpha("hotpink3"), # SCR
         ggplot2::alpha("midnightblue",0.7), # GKT
         ggplot2::alpha("darkcyan",0.7), # IGCI
         ggplot2::alpha("darkgoldenrod",0.7)) # ANM

ltyps = rep("solid",11)

methodnames = c("CISC","CE","FOI","AFC","DC","DR","HCR","SCR","GKT","IGCI","ANM")

# run experiment
run_experiment = function(tbl){
  
  stats = c(CISC_M(tbl), CE_M(tbl), FOI_M(tbl), AFC_M(tbl), DC_M(tbl), DR_M(tbl), HCR_M(tbl), 
            SCR_M(tbl), GKT_M(tbl))
  return(stats)
}

run_cont_experiment = function(tbl){
  
  stats = c(IGCI_M(tbl, "data"), ANM_M(tbl, "data"))
  return(stats)
}

# causal versus non-causal experiment over increasing sample size
abalone_fvnf_test = function(){
  
  abalone = read.table("../Data/Abalone/abalone.data", sep=",")
  abalone = abalone[,c(1:4)]
  colnames(abalone) = c("Sex","Length","Diameter","Height")
  
  # sample size points
  sm = round(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.25, 0.5, 0.75, 0.99)*nrow(abalone))
  repeats = 100 # number of repeats
  
  # ground truth: Sex -> {Length, Diameter, Height}
  gt = rep(c(1, 0, 1, 0, 1, 0), repeats) 
  
  # result containers
  expscores = vector("list", length = length(sm))
  expstats = vector("list", length = length(sm))
  for(i in 1:length(sm)){
    
    # choose sample size
    expscore = as.data.frame(matrix(0, nrow=(ncol(abalone)-1)*2*repeats, ncol=length(methodnames)))
    colnames(expscore) = methodnames
    count = 1
    for(r in 1:repeats){
      
      # repeats
      abalone_s = abalone[sample(1:nrow(abalone), sm[i]),]
      x = as.factor(abalone_s[,1]) # x is fixed
      
      # discretize data jointly (Length, Diameter and Height are discretized together)
      d = discretize.jointly(log1p(apply(abalone_s[,-1], c(1,2), as.numeric)), k=c(2:10))
      
      # study the subset of data with repeats
      for(k in 2:ncol(abalone_s)){
        
        # make table
        tbl = table(x, d$D[,k-1])
        tblx = cbind(x, log1p(as.numeric(abalone_s[,k])))
        tbly = cbind(log1p(as.numeric(abalone_s[,k])), x)
        
        
        # run experiment
        expscore[count, 1:(length(methodnames)-2)] = run_experiment(tbl)
        expscore[count, c(length(methodnames)-1, length(methodnames))] = run_cont_experiment(tblx)
        expscore[count+1, 1:(length(methodnames)-2)] = run_experiment(t(tbl))
        expscore[count+1, c(length(methodnames)-1, length(methodnames))] = run_cont_experiment(tbly)
        count = count+2
      }
    }
      
    # collect and plot results
    expscores[[i]] = expscore
    ttle = paste0("Abalone FVNF","\n",
                  "Sample Size=",sm[i])
    expstats[[i]] = ggplot.ROC.PR.curves(expscore, gt, colr, ltyps, plot=TRUE, ttle)  
  }
  return(list(expscores = expscores, expstats = expstats))
}

# causal versus independent experiment over increasing sample size
abalone_indep_test = function(){
  
  abalone = read.table("../Data/Abalone/abalone.data", sep=",")
  abalone = abalone[,c(1:4)]
  colnames(abalone) = c("Sex","Length","Diameter","Height")
  
  # sample size points
  sm = round(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.25, 0.5, 0.75, 0.99)*nrow(abalone))
  repeats = 100 # sample the 3 cause->effect 100 times
  
  # Sex -> {Length, Diameter, Height}
  gt = rep(c(1, 0, 1, 0, 1, 0), repeats) 
  
  # result containers
  expscores = vector("list", length = length(sm))
  expstats = vector("list", length = length(sm))
  for(i in 1:length(sm)){
   
    # choose sample size
    expscore = as.data.frame(matrix(0, nrow=(ncol(abalone)-1)*2*repeats, ncol=length(methodnames)))
    colnames(expscore) = methodnames
    count = 1
    for(r in 1:repeats){
      
      # repeats
      abalone_s = abalone[sample(1:nrow(abalone), sm[i]),]
      x = as.factor(abalone_s[,1]) # x is fixed
      x_i = sample(x, length(x)) # sample x to create an independent pattern
      
      # discretize data jointly (Length, Diameter and Height are discretized together)
      d = discretize.jointly(log1p(apply(abalone_s[,-1], c(1,2), as.numeric)), k=c(2:10))
      for(k in 2:ncol(abalone_s)){
        
        # get y
        y = d$D[,k-1]
        y_i = sample(y, length(y)) # sample y to create an independent pattern
	
	      yc = log1p(as.numeric(abalone_s[,k]))
	      yc_i = sample(yc, length(yc))        

        expscore[count, 1:(length(methodnames)-2)] = run_experiment(table(x, y))
        expscore[count, c(length(methodnames)-1, length(methodnames))] = run_cont_experiment(cbind(x, yc))
        expscore[count+1, 1:(length(methodnames)-2)] = run_experiment(table(x_i, y_i))
        expscore[count+1, c(length(methodnames)-1, length(methodnames))] = run_cont_experiment(cbind(x_i, yc_i))
        count = count+2
      }
    }
    
    # collect and plot results
    expscores[[i]] = expscore
    ttle = paste0("Abalone FVI","\n",
                  "Sample Size=",sm[i])
    expstats[[i]] = ggplot.ROC.PR.curves(expscore, gt, colr, ltyps, plot=TRUE, ttle)  
    
  }
  return(list(expscores = expscores, expstats = expstats))
}

# summarize causal vs non-causal results over increase sample size
abalone_accuracy_test_lines = function(expscores){
  
  # sample sizes
  # total number of samples in Abalone is 4177
  sample_per = round(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.25, 0.5, 0.75, 0.99)*4177)
  repeats = 100 # number of repeats
  
  ttle = paste0("Abalone","\n",
                "Sex -> {Length,Diameter,Height}","\n",
                "versus Reverse")
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  # score collector
  dir_scores = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(expscores)))
  colnames(dir_scores) = methodnames
  
  # get the number of correct directions
  for(i in 1:length(expscores)) {
    
    # true and transpose results are stored alternatively in abalone
    score = unlist(lapply(expscores[[i]],function(x){
      return(sum(x[c(TRUE,FALSE)] > x[c(FALSE,TRUE)]))
    }), use.names = FALSE)
    dir_scores[i,] =  (score/(3*repeats))*100
  }
  
  # plot results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="Accuracy (%)", main=ttle, xlim=c(1,length(sample_per)), 
       xaxt="n", ylim = c(0, 100), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  
  for(k in 1:length(methodnames)){
    lines(dir_scores[,k]~c(1:length(sample_per)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(dir_scores[,k]~c(1:length(sample_per)), col=colr[k], pch=pchs[k], cex=1.5)
  }  
  
  axis(side = 1, at = c(1:length(sample_per)), labels = sample_per, las=2, cex.axis=1.2)
  mtext("Sample Size (n)", side = 1, line=4.5, cex=2)
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# summarize causal vs independent results over increase sample size
abalone_indep_test_lines = function(expstats){
  
  # sample sizes
  # total number of samples in Abalone is 4177
  sample_per = round(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.25, 0.5, 0.75, 0.99)*4177)
  repeats = 100 # number of repeats
  
  ttle = paste0("Abalone","\n",
                "Sex -> {Length,Diameter,Height}","\n",
                "versus Independent")
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  # auroc collector
  aurocs = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(expstats)))
  colnames(aurocs) = methodnames
  
  for(i in 1:length(expstats)) {
    aurocs[i,] =  expstats[[i]]$AUROC
  }
  
  # plot results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="AUROC", main=ttle, xlim=c(1,length(sample_per)), 
       xaxt="n", ylim = c(0, 1), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  
  for(k in 1:length(methodnames)){
    lines((aurocs[,k])~c(1:length(sample_per)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points((aurocs[,k])~c(1:length(sample_per)), col=colr[k], pch=pchs[k], cex=1.5)
  }  
  
  axis(side = 1, at = c(1:length(sample_per)), labels = sample_per, las=2, cex.axis=1.5)
  mtext("Sample Size (n)", side = 1, line=5, cex=2)
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5, horiz = FALSE)
  par(xpd=FALSE)
  
  # aupr collector
  auprs = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(expstats)))
  colnames(auprs) = methodnames
  
  for(i in 1:length(expstats)) {
    auprs[i,] =  expstats[[i]]$AUPR
  }
  
  # plot results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="AUPR", main=ttle, xlim=c(1,length(sample_per)), 
       xaxt="n", ylim = c(0, 1), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  
  for(k in 1:length(methodnames)){
    lines(auprs[,k]~c(1:length(sample_per)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(auprs[,k]~c(1:length(sample_per)), col=colr[k], pch=pchs[k], cex=1.5)
  }  
  
  axis(side = 1, at = c(1:length(sample_per)), labels = sample_per, las=2, cex.axis=1.5)
  mtext("Sample Size (n)", side = 1, line=4.5, cex=2)
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# Main procedure
AbaloneEval = function(){
  
  # causal versus non-functional evaluation
  pdf("../Results/AbaloneStudy/abalone_fvnf_auroc_aupr.pdf")
  abalone_fvnf = abalone_fvnf_test()
  dev.off()
  save(abalone_fvnf, file="../Results/AbaloneStudy/abalone_fvnf.Rds")
  
  # causal versus independent evaluation
  pdf("../Results/AbaloneStudy/abalone_fvi_auroc_aupr.pdf")
  abalone_fvi = abalone_indep_test()
  dev.off()
  save(abalone_fvi, file="../Results/AbaloneStudy/abalone_fvi.Rds")
  
  # plot directional accuracy over sample size
  pdf("../Results/AbaloneStudy/abalone_dir_test_over_sm.pdf")
  abalone_accuracy_test_lines(abalone_fvnf$expscores)
  dev.off()
  
  # plot AUROC/AUPR (functional versus independent) over sample size
  pdf("../Results/AbaloneStudy/abalone_fvi_test_over_sm.pdf")
  abalone_indep_test_lines(abalone_fvi$expstats)
  dev.off()
  
}
