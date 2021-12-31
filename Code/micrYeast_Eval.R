# Processing and Evaluation on:
# Perturbed Yeast Microarray for:
# Overcoming biases in causal inference of molecular interactions

# Manuscript:
# "Learning causal networks using inducible transcription factors and transcriptome‐wide time series"
# Hackett.et.al 2020, Molecular Systems Biology.
# Original data obtained from: https://idea.research.calicolabs.com/data

# Ground Truth:
# "YEASTRACT+: a portal for cross-species comparative genomics of transcription regulation in yeasts"
# Interactions obtained from http://www.yeastract.com/formfindregulators.php

# This code requires idea_tall_expression_data.tsv
# from https://idea.research.calicolabs.com/data to be put in Data/Yeast/ directory.
# Ground truth interactions (Yeast_TFTG.csv) are available in Data/Yeast/ directory.

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


# Construct a matrix using the selected metric ('colnum' in the yeast_microarray) 
Construct_dataset = function(yeast_microarray, gindex, colnum){
  
  # sample names for across all genes were found to be identical
  # create a custom sample name combining TF and time.
  sample_names = paste0("TF:",yeast_microarray$TF[gindex[,1]],";", 
                        "Time:",yeast_microarray$time[gindex[,1]])
  
  # prepare an empty data.frame
  yeast_microarray_matrix = as.data.frame(matrix(0, ncol=nrow(gindex), nrow=ncol(gindex)))
  colnames(yeast_microarray_matrix) = sample_names
  rownames(yeast_microarray_matrix) = genes
  
  for(i in 1:ncol(gindex)){
    yeast_microarray_matrix[i,] = yeast_microarray[gindex[,i], colnum]
  }
  
  # return data
  return(yeast_microarray_matrix)
}

# Find the corresponding TF:Time entry for each gene in the tall expression matrix
FindGIndex = function(yeast_microarray){
  
  genes = unique(yeast_microarray$GeneName)
  unique_samples = nrow(yeast_microarray)/length(genes)
  
  gindex = as.data.frame(matrix(0,ncol=length(genes), nrow=unique_samples))
  pb = txtProgressBar(min=1, max=length(genes), style = 3)
  for(i in 1:ncol(gindex)){
    gindex[,i] = which(!is.na(match(yeast_microarray$GeneName, genes[i])))
    setTxtProgressBar(pb, i)
  }
  
  # return indices
  return(gindex)
}

# Find TF and TG indices in the yeast_microarray data
FindTFTGIndex = function(yeast_microarray, TFTG){
  
  genes = rownames(yeast_microarray)
  
  TFind = match(toupper(TFTG[,1]), paste(genes,"P",sep=""))
  TGind = match(toupper(TFTG[,2]), genes)
  
  # remove entries with no match in the given data
  if(any(is.na(TFind))){
    rem_ind = which(is.na(TFind))
  }
  
  if(any(is.na(TGind))){
    rem_ind = which(is.na(TGind))
  }
  
  TFTG = TFTG[-unique(rem_ind),]
  
  TFind = TFind[-unique(rem_ind)]
  TGind = TGind[-unique(rem_ind)]
  
  TFTG = data.frame(TF=TFTG[,1], TG=TFTG[,2], TFind=TFind, TGind=TGind)
  
  # return TFTG with indices
  return(TFTG)
}

# causal versus non-causal experiment over increasing sample size
yeast_data_func_vs_indep = function(){
  
  load("../Data/Yeast/yeast_microarray_cleaned_log_ratio.RData")
  
  # save original
  yeast_microarray_og = yeast_microarray
  
  # sample sizes
  sample_per = c(0.05, 0.25, 0.5, 0.75, 1)
  
  # result containe
  expscores = vector("list", length = 5)
  
  # obtain parent and child indices in yeast_microarray
  parents = TFTG$TFind
  children = TFTG$TGind
  
  # repeat experiment at different sample sizes
  
  for(j in 1:length(sample_per))
  {
    sm_size = round(sample_per[j] * ncol(yeast_microarray_og)) # calculate sample size
    
    # sample a subset
    yeast_microarray = yeast_microarray_og[,sample(1:ncol(yeast_microarray_og), sm_size)]
    
    # obtain scores for all methods
    results = pblapply(c(1:length(parents)), function(i){
      
      p = as.numeric(yeast_microarray[parents[i],])
      c = as.numeric(yeast_microarray[children[i],])
      tblx = cbind(p, c)
      
      p_i = sample(p, length(p))
      c_i = sample(c, length(c))
      tblxi = cbind(p_i, c_i)
      
      # Calculate number of unique observations
      
      unq_dp = sum(!duplicated(cbind(p,c)))
      if(unq_dp < 50){
        # if less than 50 unique observations, then set unq_dp such that each cluster has atleast 5 points.
        unq_dp = round(unq_dp/5)
      }
      
      if(unq_dp < 3){ # if less than 3 unique points
        k = 2 # set k to 2
      } else {
        k = c(2:min(unq_dp, 10)) # set k in range of 2 to min(unq_dp, 10)
      }
      
      # Discretize jointly using GridOnCluster
      d = discretize.jointly(cbind(p,c), k = k)  
      
      # initiate garbage collection  
      gc()
      
      # prepare contingency table (causal)
      p = d$D[,1]
      c = d$D[,2]
      tbl = table(p, c)  
      
      # prepare contingency table (independent)
      p_i = sample(p, length(p))
      c_i = sample(c, length(c))
      tbl_i = table(p_i, c_i)
      
      # if the table has only 1 row or column, then return the worst possible score for all methods  
      if(ncol(tbl) == 1 || nrow(tbl) == 1){
        return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                     rep(-.Machine$integer.max, length(methodnames))))
      } else { # run experiment
        return(rbind(c(run_experiment(tbl), run_cont_experiment(tblx)), 
                     c(run_experiment(tbl_i), run_cont_experiment(tblxi))))
      }
      
    })
    
    # Convert result to data.frame
    results = ListToDataMatrixMethods(results)
    results = as.data.frame(results)
    colnames(results) = methodnames
    
    # ground truth
    gt = rep(c(1,0), each=(nrow(results)/2))
    
    # plot results
    ttle = paste0("Perturbed Microarray Yeast","\n",
                  "TF->TG vs Shuff. TF->TG","\n",
                  "Sample Size = ",sm_size)
    stats = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
    expscores[[j]] = list(stats = stats, scores = results)
  }
  
  # return result
  return(expscores)
}

# causal versus independent experiment over increasing sample size
yeast_data_func_vs_nfunc = function(){
  
  load("../Data/Yeast/yeast_microarray_cleaned_log_ratio.RData")
  
  # save original
  yeast_microarray_og = yeast_microarray
  
  # sample sizes
  sample_per = c(0.05, 0.25, 0.5, 0.75, 1)
  
  # result container
  expscores = vector("list", length = 5)
  
  # obtain parent and child indices in yeast_microarray
  parents = TFTG$TFind
  children = TFTG$TGind
  
  # repeat experiment at different sample sizes
  
  for(j in 1:length(sample_per))
  {
    sm_size = round(sample_per[j] * ncol(yeast_microarray_og)) # calculate sample size
    
    # sample a subset
    yeast_microarray = yeast_microarray_og[,sample(1:ncol(yeast_microarray_og), sm_size)]
    
    # obtain scores for all methods
    results = pblapply(c(1:length(parents)), function(i){
      
      p = as.numeric(yeast_microarray[parents[i],])
      c = as.numeric(yeast_microarray[children[i],])
      tblx = cbind(p, c)
      tbly = cbind(c, p)
      
      # Calculate number of unique observations
      
      unq_dp = sum(!duplicated(cbind(p,c)))
      if(unq_dp < 50){
        # if less than 50 unique observations, then set unq_dp such that each cluster has atleast 5 points.
        unq_dp = round(unq_dp/5)
      }
      
      if(unq_dp < 3){ # if less than 3 unique points
        k = 2 # set k to 2
      } else {
        k = c(2:min(unq_dp, 10)) # set k in range of 2 to min(unq_dp, 10)
      }
      
      # Discretize jointly using GridOnCluster
      d = discretize.jointly(cbind(p,c), k = k)  
      
      # initiate garbage collection  
      gc()
      
      # prepare contingency table
      p = d$D[,1]
      c = d$D[,2]
      tbl = table(p, c)  
      
      # if the table has only 1 row or column, then return the worst possible score for all methods  
      if(ncol(tbl) == 1 || nrow(tbl) == 1){
        return(rbind(rbind(rep(-.Machine$integer.max, length(methodnames)),
                           rep(-.Machine$integer.max, length(methodnames)))))
      } else { # run experiment
        return(rbind(c(run_experiment(tbl), run_cont_experiment(tblx)), 
                     c(run_experiment(t(tbl)), run_cont_experiment(tbly))))
      }
      
    })
    
    # Convert result to data.frame
    results = ListToDataMatrixMethods(results)
    results = as.data.frame(results)
    colnames(results) = methodnames
    
    # ground truth
    gt = rep(c(1,0), each=(nrow(results)/2))
    
    # plot results
    ttle = paste0("Perturbed Microarray Yeast","\n",
                  "TF->TG vs TG->TF","\n",
                  "Sample Size = ",sm_size)
    stats = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
    expscores[[j]] = list(stats = stats, scores = results)
  }
  
  # return result
  return(expscores)
}

# summarize causal vs non-causal results over increase sample size
yeast_accuracy_lines = function(expscores){
  
  # sample sizes
  # total number of samples in yeast is 1693
  sample_per = round(c(0.05, 0.25, 0.5, 0.75, 1)*1693)
  
  ttle = paste0("Pertubed Yeast Microarray","\n",
                "TF->TG versus TG->TF")
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  # score collector
  dir_scores = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(expscores)))
  colnames(dir_scores) = methodnames
  
  # get the number of correct directions over sample size
  for(j in 1:length(expscores)){
    
    scores = expscores[[j]]$scores
    ntest = nrow(scores)/2
    
    dir_test = (unlist(lapply(scores, function(i){
      return(sum(i[1:ntest] > i[(ntest+1):(ntest*2)]))
    }), use.names = FALSE) * 100)/ntest
    
    dir_scores[j,] = dir_test
  }
  
  # plot results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="Accuracy (%)", main=ttle, xaxt="n", xlim=c(1,length(sample_per)),
       ylim = c(0, 100), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  
  for(k in 1:length(methodnames)){
    lines(dir_scores[,k]~c(1:length(sample_per)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(dir_scores[,k]~c(1:length(sample_per)), col=colr[k], pch=pchs[k], cex=2)
  }  
  
  axis(side=1, at=c(1:length(sample_per)), labels = sample_per, las=2, cex.axis=1.5)
  mtext("Samples (n)", side=1, line=5, cex=2)
  
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# summarize causal vs independent results over increase sample size
yeast_indep_test_lines = function(expstats){
  
  # sample sizes
  # total number of samples in yeast is 1693
  sample_per = round(c(0.05, 0.25, 0.5, 0.75, 1)*1693)
  
  ttle = paste0("Pertubed Yeast Microarray","\n",
                "TF->TG versus Independent")
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  # auroc collector
  aurocs = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(expstats)))
  colnames(aurocs) = methodnames
  for(j in 1:length(expstats)){
    aurocs[j,] = expstats[[j]]$stats$AUROC
  }
  
  # plot auroc results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="AUROC", main=ttle, xaxt="n", xlim=c(1,length(sample_per)),
       ylim = c(0, 1), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  for(k in 1:length(methodnames)){
    lines(aurocs[,k]~c(1:length(sample_per)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(aurocs[,k]~c(1:length(sample_per)), col=colr[k], pch=pchs[k], cex=2)
  }  
  
  axis(side=1, at=c(1:length(sample_per)), labels = sample_per, las=2, cex.axis=1.5)
  mtext("Samples (n)", side=1, line=5, cex=2)
  
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
  # aupr collector
  auprs = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(expstats)))
  colnames(auprs) = methodnames
  for(j in 1:length(expstats)){
    auprs[j,] = expstats[[j]]$stats$AUPR
  }
  
  # plot aupr results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="AUPR", main=ttle, xaxt="n", xlim=c(1,length(sample_per)),
       ylim = c(0, 1), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  for(k in 1:length(methodnames)){
    lines(auprs[,k]~c(1:length(sample_per)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(auprs[,k]~c(1:length(sample_per)), col=colr[k], pch=pchs[k], cex=2)
  }  
  
  axis(side=1, at=c(1:length(sample_per)), labels = sample_per, las=2, cex.axis=1.5)
  mtext("Samples (n)", side=1, line=5, cex=2)
  
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# prepare data for analysis
Prepare_data = function(){
  
  # read all processing levels of data
  yeast_microarray = read.table("../Data/Yeast/idea_tall_expression_data.tsv", quote = "",
                                stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
  # get gindex
  gindex = FindGIndex(yeast_microarray)
  
  # log2 cleaned ratio
  yeast_microarray = Construct_dataset(yeast_microarray, gindex, 11)
  yeast_microarray = apply(yeast_microarray_log2ratiocl, c(1,2), as.numeric)
  
  # Get TFTG
  TFTG = read.csv("../Data/Yeast/Yeast_TFTG.csv")[,c(2,3)]
  
  # Find TFTG index in data
  TFTG = FindTFTGIndex(yeast_microarray, TFTG)
  
  # save data
  rm(list=setdiff(ls(), c("yeast_microarray","TFTG")))
  save.image("../Data/Yeast/yeast_microarray_cleaned_log_ratio.RData")
  
}

# Main procedure
YeastEval = function(){
  
  # Prepare data
  if(!file.exists("../Data/Yeast/yeast_microarray_cleaned_log_ratio.RData")){
    Prepare_data()
  }
  # causal versus independent evaluation
  pdf("../Results/YeastStudy/yeast_fvi_roc_pr.pdf")
  yeast_fvi = yeast_data_func_vs_indep()
  dev.off()
  save(yeast_fvi, file="../Results/YeastStudy/yeast_fvi.rds") 
  
  # causal versus non-causal evaluation
  pdf("../Results/YeastStudy/yeast_fvnf_roc_pr.pdf")
  yeast_fvnf = yeast_data_func_vs_nfunc()
  dev.off()
  save(yeast_fvnf, file="../Results/YeastStudy/yeast_fvnf.rds")
  
  # plot directional accuracy over sample size
  pdf("../Results/YeastStudy/yeast_dir_test_over_sm.pdf")
  yeast_accuracy_lines(yeast_fvnf)
  dev.off()
  
  # plot AUROC/AUPR (functional versus independent) over sample size
  pdf("../Results/YeastStudy/yeast_fvi_test_over_sm.pdf")
  yeast_indep_test_lines(yeast_fvi)
  dev.off()
}
