# Simulation study setup for:
# Overcoming biases in causal inference of molecular interactions

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
            SCR_M(tbl), GKT_M(tbl), IGCI_M(tbl), ANM_M(tbl))
  return(stats)
}
# show functional tables
show_tables = function(){
  
  pdf("../Results/Supplementary/portrait_table.pdf", height = 7, width = 5)
  portrait = matrix(c(10,0,0, 0,10,0, 0,0,10, 0,0,10, 0,10,0, 10,0,0), 6, 3, T)
  plot_table(portrait, main="Portrait Table (|X| > |Y|)", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
             col="cornflowerblue",value.cex = 3)
  mtext("Y", side = 1, line=3, cex=2)
  dev.off()
  
  pdf("../Results/Supplementary/square_table.pdf", height = 5, width = 5)
  square = matrix(c(10,0,0, 0,0,10, 0,10,0), 3, 3, T)
  plot_table(square, main="Square Table (|X| = |Y|)", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
             col="cornflowerblue",value.cex = 3)
  mtext("Y", side = 1, line=3, cex=2)
  dev.off()
  
  pdf("../Results/Supplementary/landscape_table.pdf", height = 5, width = 7)
  landscape = matrix(c(10,0,0,0,0,0, 0,10,0,0,0,0, 0,0,10,0,0,0), 3, 6, T)
  plot_table(landscape, main="Landscape Table (|X| < |Y|)", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
             col="cornflowerblue",value.cex = 3)
  mtext("Y", side = 1, line=3, cex=2)
  dev.off()
  
}

# Support function for computing and plotting AUROC over noise
PlotAUROC_over_Noise = function(ttle, expstats, methodnames, noise_levels, ltyps2, lwds, pchs){
  
  expauroc = as.data.frame(matrix(0, nrow=length(expstats), ncol=length(methodnames)))
  for(i in 1:length(expstats)){
    expauroc[i,] = expstats[[i]]$AUROC
  }
  colnames(expauroc) = methodnames
  
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="AUROC", main=ttle, xlim = c(0.01,1), 
       ylim = c(0, 1), cex.lab=2, cex.main=2, type="n", cex.axis=1.5, xaxt="n")
  grid()
  
  for(k in 1:length(methodnames)){
    lines(expauroc[,k]~noise_levels, col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(expauroc[,k]~noise_levels, col=colr[k], pch=pchs[k], cex=2)
  } 
  
  axis(side = 1, at = noise_levels, labels = noise_levels, las=2, cex.axis=1.5)
  mtext("Noise Levels", side = 1, line=5, cex=2)
  
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# Support function for computing and plotting directional accuracy over noise
PlotDIR_over_Noise = function(ttle, expscores, methodnames, noise_levels, ltyps2, lwds, pchs){
  
  dir_scores = as.data.frame(matrix(0, nrow=length(expscores), ncol=length(methodnames)))
  for(i in 1:length(expscores)){
    score = lapply(expscores[[i]],function(x){
      return(x[1,] > x[2,])
    })
    dir_scores[i,] =  apply(ListToDataMatrix(score), 2, sum)/200 
  }
  colnames(dir_scores) = methodnames
  
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="Accuracy (%)", main=ttle, xlim = c(0.01,1), xaxt="n", 
       ylim = c(0, 100), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  for(k in 1:length(methodnames)){
    lines((dir_scores[,k]*100)~noise_levels, col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points((dir_scores[,k]*100)~noise_levels, col=colr[k], pch=pchs[k], cex=2)
  }  
  axis(side = 1, at = noise_levels, labels = noise_levels, las=2, cex.axis=1.5)
  mtext("Noise Levels", side = 1, line=5, cex=2)
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# Summarize functional versus non-functional directional accuracy scores over an increasing decsision rate
# by means of violin plot.
# Requires a list of scores obtained from func_vs_nfunc() in SimStudy.R
SimStudy_FVNF_dir_scores_drate = function(scores){
  
  # decision rates
  drate = c(0.05, 0.20, 0.40, 0.60, 0.80, 1)
  
  nsm = length(scores) # samples
  nexp = 2 # marginals
  ntab = length(scores$fvnf2000$exp1stats) # table sizes
  nnoi = length(scores$fvnf2000$exp1stats[[1]]) # noise levels
  
  ########################################################################################################
  
  for(d in drate){
    
    # Directional accuracy score distribution
    dir_scores = as.data.frame(matrix(0, ncol=length(methodnames), nrow=nsm*nexp*ntab*nnoi))
    colnames(dir_scores) = methodnames
    count = 1
    
    # sample size 2000
    for(i in 1:ntab) { # tables
      for(j in 1:nnoi) { # noise levels
        
        # obtain the difference in scores between functional and non-functional patterns
        score = lapply(scores$fvnf2000$exp1scores[[i]][[j]],function(x){
          return(x[1,] - x[2,])
        })
        
        # convert scores to a matrix
        score = ListToDataMatrix(score)
        
        # determine the percentage of top differences to consider (based on drate, d)
        len = round(d*nrow(score))
        
        # compute effective score which is the number of top differences that are positive
        # Note: methods in Methods.R return scores that are higher the better
        eff_score = rep(0, length(methodnames))
        for(m in 1:length(methodnames)){
          x = sort(score[,m], decreasing = TRUE)[1:len]
          eff_score[m] = sum(x>0)
        }
        dir_scores[count,] =  eff_score/len # average over number of pairs considered
        
        # repeat for 2nd marginal scheme
        score = lapply(scores$fvnf2000$exp2scores[[i]][[j]],function(x){
          return(x[1,] - x[2,])
        })
        score = ListToDataMatrix(score)
        
        len = round(d*nrow(score))
        
        eff_score = rep(0, length(methodnames))
        for(m in 1:length(methodnames)){
          x = sort(score[,m], decreasing = TRUE)[1:len]
          eff_score[m] = sum(x>0)
        }
        dir_scores[count+1,] =  eff_score/len
        
        count = count+2
      }  
    }
    
    # sample size 200
    for(i in 1:ntab) { # tables
      for(j in 1:nnoi) { # noise levels
        
        # obtain the difference in scores between functional and non-functional patterns
        score = lapply(scores$fvnf200$exp1scores[[i]][[j]],function(x){
          return(x[1,] - x[2,])
        })
        
        # convert scores to a matrix
        score = ListToDataMatrix(score)
        
        # determine the percentage of top differences to consider (based on drate, d)
        len = round(d*nrow(score))
        
        # compute effective score which is the number of top differences that are positive
        # Note: methods in Methods.R return scores that are higher the better
        eff_score = rep(0, length(methodnames))
        for(m in 1:length(methodnames)){
          x = sort(score[,m], decreasing = TRUE)[1:len]
          eff_score[m] = sum(x>0)
        }
        dir_scores[count,] =  eff_score/len # average over number of pairs considered
        
        # repeat for 2nd marginal scheme
        score = lapply(scores$fvnf200$exp2scores[[i]][[j]],function(x){
          return(x[1,] - x[2,])
        })
        
        score = ListToDataMatrix(score)
        
        len = round(d*nrow(score))
        
        eff_score = rep(0, length(methodnames))
        for(m in 1:length(methodnames)){
          x = sort(score[,m], decreasing = TRUE)[1:len]
          eff_score[m] = sum(x>0)
        }
        dir_scores[count+1,] =  eff_score/len
        
        count = count+2
      }
    }
    
    # rank by mean directional accuracy distribution
    mean_dir = unlist(lapply(dir_scores, mean), use.names = FALSE)
    dir_scores = dir_scores[,order(mean_dir, decreasing = TRUE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(dir_scores+jitter(rep(0, length(methodnames)),amount=0.01),
                     col=colr[order(mean_dir, decreasing = TRUE)], 
                     main=paste0("Functional vs Non-Functional","\n",
                                 "D=",d*100,"%, Scores over 120 setups"), 
                     cex.axis=2.5, cex.main=2.6, las=2,
                     plotCentre="line", lwd=3, ylim=c(0,1))
    par(lwd=0.3)
    stripchart(dir_scores, add=TRUE, vertical = TRUE, jitter = 0.1, pch=23, cex=0.8,
               col="gray", bg="white", method = "jitter")
    abline(h=0.5, lty=2, lwd=3, col="blue")
    mtext(text="Accuracy", side = 2, line=4, cex=2.5)
  }
  
  ########################################################################################################
}

# summarize simulation study -- Functional versus Independent (AUROC)
# Requires 'Data/Simulation/SimResults.RData'
Summarize_simulation_study_FVI = function(){
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  noise_levels = c(0.01, 0.25, 0.5, 0.75, 1)
  
  load("../Results/SimStudy/fvi200.RData")
  load("../Results/SimStudy/fvi2000.RData")
  load("../Results/SimStudy/fvnf200.RData")
  load("../Results/SimStudy/fvnf2000.RData")
  
  ######################################################################################
  # Summarize over sample size
  # two table sizes 10x8 and 8x10
  # Marginal fixed at non-uniform column and row
  ######################################################################################
  
  # subset
  sm_fvi2000 = list(expscores = fvi2000$exp2scores[c(1,6)], expstats = fvi2000$exp2stats[c(1,6)])
  
  # sample size 2000
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/SM2000_10x8_NUC_NUR_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 10x8 table, NUR, NUC")
  PlotAUROC_over_Noise(ttle, sm_fvi2000$expstats[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # sample size 2000
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/SM2000_8x10_NUC_NUR_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 8x10 table, NUR, NUC")
  PlotAUROC_over_Noise(ttle, sm_fvi2000$expstats[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  # subset
  sm_fvi200 = list(expscores = fvi200$exp2scores[c(1,6)], expstats = fvi200$exp2stats[c(1,6)])
  
  # sample size 200
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/SM200_10x8_NUC_NUR_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "200 Samples, 10x8 table, NUR, NUC")
  PlotAUROC_over_Noise(ttle, sm_fvi200$expstats[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # sample size 200
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/SM200_8x10_NUC_NUR_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "200 Samples, 8x10 table, NUR, NUC")
  PlotAUROC_over_Noise(ttle, sm_fvi200$expstats[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  
  ######################################################################################
  # Summarize over marginal distributions
  # two table sizes 10x8 and 8x10
  # Sample size fixed at 2000
  ######################################################################################
  
  # subset
  mar_fvi2000 = list(expstats1 = fvi2000$exp1stats[c(1,6)], expstats2 = fvi2000$exp2stats[c(1,6)],
                     expstats3 = fvi2000$exp3stats[c(1,6)])
  
  # marginal UR NUC
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/UR_NUC_SM2000_10x8_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 10x8 table, UR, NUC")
  PlotAUROC_over_Noise(ttle, mar_fvi2000$expstats1[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # marginal UR NUC
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/UR_NUC_SM2000_8x10_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 8x10 table, UR, NUC")
  PlotAUROC_over_Noise(ttle, mar_fvi2000$expstats1[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  
  ######################################################################################
  
  # marginal NUR NUC -- Already have these
  
  ######################################################################################
  
  # marginal NUR UC
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/NUR_UC_SM2000_10x8_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 10x8 table, NUR, UC")
  PlotAUROC_over_Noise(ttle, mar_fvi2000$expstats3[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # marginal NUR UC
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/NUR_UC_SM2000_8x10_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 8x10 table, NUR, UC")
  PlotAUROC_over_Noise(ttle, mar_fvi2000$expstats3[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  ######################################################################################
  # Summarize over table size
  # Sample size fixed at 2000
  # Marginal fixed at non-uniform column and row
  ######################################################################################
  
  # subset
  tb_fvi2000 = fvi2000$exp2stats[c(1,3,6)]
  
  # table size 10x8 only AUROC -- Already have this
  
  ######################################################################################
  
  # table size 9x9 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVI/9x9_SM2000_NUC_NUR_FVI.pdf")
  ttle = paste0("Functional versus Independent","\n",
                "2000 Samples, 9x9 table, NUR, NUC")
  PlotAUROC_over_Noise(ttle, tb_fvi2000[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  # table size 8x10 only AUROC -- Already have this
  
  ######################################################################################
  
}

# summarize simulation study -- Functional versus non-functional (AUROC)
# Requires 'Data/Simulation/SimResults.RData'
Summarize_simulation_study_FVNF = function(){
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  noise_levels = c(0.01, 0.25, 0.5, 0.75, 1)
  
  load("../Results/SimStudy/fvi200.RData")
  load("../Results/SimStudy/fvi2000.RData")
  load("../Results/SimStudy/fvnf200.RData")
  load("../Results/SimStudy/fvnf2000.RData")
  
  ######################################################################################
  # Summarize over sample size
  # two table sizes 10x8 and 8x10
  # Marginal fixed at uniform column
  ######################################################################################
  
  # subset
  sm_fvnf2000 = fvnf2000$exp2stats[c(1,6)]
  
  # sample size 2000
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/SM2000_10x8_UC_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 10x8 table, UC")
  PlotAUROC_over_Noise(ttle, sm_fvnf2000[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # sample size 2000
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/SM2000_8x10_UC_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 8x10 table, UC")
  PlotAUROC_over_Noise(ttle, sm_fvnf2000[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
 
  # subset
  sm_fvnf200 = fvnf200$exp2stats[c(1,6)]
  
  # sample size 200
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/SM200_10x8_UC_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "200 Samples, 10x8 table, UC")
  PlotAUROC_over_Noise(ttle, sm_fvnf200[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # sample size 200
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/SM200_8x10_UC_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "200 Samples, 8x10 table, UC")
  PlotAUROC_over_Noise(ttle, sm_fvnf200[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  
  ######################################################################################
  # Summarize over marginal distributions
  # two table sizes 10x8 and 8x10
  # Sample size fixed at 2000
  ######################################################################################

  # subset
  mar_fvnf2000 = list(expstats1 = fvnf2000$exp1stats[c(1,6)], expstats2 = fvnf2000$exp2stats[c(1,6)])
  
  # marginal UR
  # table size 10x8 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/UR_SM2000_10x8_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 10x8 table, UR")
  PlotAUROC_over_Noise(ttle, mar_fvnf2000$expstats1[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # marginal UR
  # table size 8x10 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/UR_SM2000_8x10_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 8x10 table, UR")
  PlotAUROC_over_Noise(ttle, mar_fvnf2000$expstats1[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  
  ######################################################################################
  
  # marginal UC -- Already have these
  
  ######################################################################################
  
  ######################################################################################
  # Summarize over table size
  # Sample size fixed at 2000
  # Marginal fixed at uniform column
  ######################################################################################
  
  # subset
  tb_fvnf2000 = fvnf2000$exp2stats[c(1,3,6)]
  
  # table size 10x8 only AUROC -- Already have this
  
  ######################################################################################
  
  # table size 9x9 only AUROC
  pdf("../Results/Supplementary/SimStudy_FVNF/9x9_SM2000_UC_FVNF.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 9x9 table, UC")
  PlotAUROC_over_Noise(ttle, tb_fvnf2000[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  # table size 8x10 only AUROC -- Already have this
  
  ######################################################################################
  
}

# summarize simulation study -- Functional versus Non-functional (DIR)
# Requires 'Data/Simulation/SimResults.RData'
Summarize_simulation_study_FVNF_DIR = function(){
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  noise_levels = c(0.01, 0.25, 0.5, 0.75, 1)
  
  load("../Results/SimStudy/fvi200.RData")
  load("../Results/SimStudy/fvi2000.RData")
  load("../Results/SimStudy/fvnf200.RData")
  load("../Results/SimStudy/fvnf2000.RData")
  
  ######################################################################################
  # Summarize over sample size
  # two table sizes 10x8 and 8x10
  # Marginal fixed at uniform column
  ######################################################################################
  
  # subset
  sm_fvnf2000 = fvnf2000$exp2scores[c(1,6)]
  
  # sample size 2000
  # table size 10x8 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/SM2000_10x8_UC_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 10x8 table, UC")
  PlotDIR_over_Noise(ttle, sm_fvnf2000[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # sample size 2000
  # table size 8x10 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/SM2000_8x10_UC_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 8x10 table, UC")
  PlotDIR_over_Noise(ttle, sm_fvnf2000[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  # subset
  sm_fvnf200 = fvnf200$exp2scores[c(1,6)]
  
  # sample size 200
  # table size 10x8 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/SM200_10x8_UC_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "200 Samples, 10x8 table, UC")
  PlotDIR_over_Noise(ttle, sm_fvnf200[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # sample size 200
  # table size 8x10 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/SM200_8x10_UC_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "200 Samples, 8x10 table, UC")
  PlotDIR_over_Noise(ttle, sm_fvnf200[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  
  ######################################################################################
  # Summarize over marginal distributions
  # two table sizes 10x8 and 8x10
  # Sample size fixed at 2000
  ######################################################################################
  
  # subset
  mar_fvnf2000 = list(expscores1 = fvnf2000$exp1scores[c(1,6)], expscores2 = fvnf2000$exp2scores[c(1,6)])
  
  # marginal UR
  # table size 10x8 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/UR_SM2000_10x8_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 10x8 table, UR")
  PlotDIR_over_Noise(ttle, mar_fvnf2000$expscores1[[1]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  # marginal UR
  # table size 8x10 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/UR_SM2000_8x10_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 8x10 table, UR")
  PlotDIR_over_Noise(ttle, mar_fvnf2000$expscores1[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  
  ######################################################################################
  
  # marginal UC -- Already have these
  
  ######################################################################################
  
  ######################################################################################
  # Summarize over table size
  # Sample size fixed at 2000
  # Marginal fixed at uniform column
  ######################################################################################
  
  # subset
  tb_fvnf2000 = fvnf2000$exp2scores[c(1,3,6)]
  
  # table size 10x8 DIR -- Already have this
  
  ######################################################################################
  
  # table size 9x9 DIR
  pdf("../Results/Supplementary/SimStudy_FVNF/9x9_SM2000_UC_FVNF_DIR.pdf")
  ttle = paste0("Functional versus Non-Functional","\n",
                "2000 Samples, 9x9 table, UC")
  PlotDIR_over_Noise(ttle, tb_fvnf2000[[2]], methodnames, noise_levels, ltyps2, lwds, pchs)
  dev.off()
  
  ######################################################################################
  
  # table size 8x10 only AUROC -- Already have this
  
  ######################################################################################
  
}

# Plot a snapshot of the simulated data for both
# functional versus non-functional and functional versus independent
# Requires 'Data/Simulation/SimResults.RData'
PlotSimData = function(){
  
  load("../Results/SimStudy/fvi200.RData")
  load("../Results/SimStudy/fvi2000.RData")
  load("../Results/SimStudy/fvnf200.RData")
  load("../Results/SimStudy/fvnf2000.RData")
  
  # For function versus non-function
  # sample size=2000
  
  # noise levels 0.01, 0.5 and 1
  # table size 7x5 and 5x7
  
  # randomly choose a table from exp1, 7x5
  tnum = sample(1:200, 1)
  sub_pop = fvnf2000$populations[c(2,5)]
  pdf("../Results/Supplementary/SimStudy_FVNF/Population_exp1_SM2000_7x5_FVNF.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[1]]$exp1fun[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
    
    t2 = sub_pop[[1]]$exp1nfun[[i]][[tnum]]
    plot_table(t2, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  # randomly choose a table from exp1, 5x7
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVNF/Population_exp1_SM2000_5x7_FVNF.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[2]]$exp1fun[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
    
    t2 = sub_pop[[2]]$exp1nfun[[i]][[tnum]]
    plot_table(t2, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  ######################################################################################
  
  # randomly choose a table from exp2, 7x5
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVNF/Population_exp2_SM2000_7x5_FVNF.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[1]]$exp2fun[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
    
    t2 = sub_pop[[1]]$exp2nfun[[i]][[tnum]]
    plot_table(t2, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  # randomly choose a table from exp2, 5x7
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVNF/Population_exp2_SM2000_5x7_FVNF.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[2]]$exp2fun[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
    
    t2 = sub_pop[[2]]$exp2nfun[[i]][[tnum]]
    plot_table(t2, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  ######################################################################################

  # For function versus independent
  # sample size=2000
  
  # noise levels 0.01, 0.5 and 1
  # table size 7x5 and 5x7
  
  # randomly choose a functional table, 7x5
  tnum = sample(1:200, 1)
  sub_pop = fvi2000$populations[c(2,5)]
  pdf("../Results/Supplementary/SimStudy_FVI/Population_fun_SM2000_7x5_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[1]]$exp1fun[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  # randomly choose a functional table, 5x7
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_fun_SM2000_5x7_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[2]]$exp1fun[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  ######################################################################################
  
  # randomly choose a independent table from exp1, 7x5
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_exp1_SM2000_7x5_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[1]]$exp1ind[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  # randomly choose a independent table from exp1, 5x7
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_exp1_SM2000_5x7_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[2]]$exp1ind[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  ######################################################################################
  
  # randomly choose a independent table from exp2, 7x5
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_exp2_SM2000_7x5_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[1]]$exp2ind[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  # randomly choose a independent table from exp2, 5x7
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_exp2_SM2000_5x7_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[2]]$exp2ind[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  ######################################################################################
  
  # randomly choose a independent table from exp3, 7x5
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_exp3_SM2000_7x5_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[1]]$exp3ind[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
  # randomly choose a independent table from exp3, 5x7
  tnum = sample(1:200, 1)
  pdf("../Results/Supplementary/SimStudy_FVI/Population_exp3_SM2000_5x7_FVI.pdf")
  for(i in 1:3){
    
    t1 = sub_pop[[2]]$exp3ind[[i]][[tnum]]
    plot_table(t1, main="", xlab="", ylab="X", cex.lab=2, mar = c(3,3,3,3), cex.main=2, 
               col="cornflowerblue")
    mtext("Y", side = 1, line=3, cex=2)
  }
  dev.off()
  
}

# Evaluating directional inference on a subset of yeast data using finer discretization
yeast_data_func_vs_nfunc_finer = function(ncores){
  
  load("../Data/Yeast/yeast_microarray_cleaned_log_ratio.RData")
  
  # save original
  yeast_microarray_og = yeast_microarray
  
  # result collector
  expscores = vector("list", length = 10)
  
  # obtain parent and child indices in yeast_microarray
  parents = TFTG$TFind
  children = TFTG$TGind
  
  # randomly choose 1000 non self-cycle causal pairs
  self.cycle = which(parents == children)
  random_pairs = sample(c(1:length(parents))[-self.cycle], 1000)
  parents = parents[random_pairs]
  children = children[random_pairs]
  
  # repeat experiment 10 times
  for(j in 1:10)
  {
    # get random subset of 100 samples
    sm = sample(ncol(yeast_microarray_og), 100)
    yeast_microarray = yeast_microarray_og[,sm]
    
    # prepare cluster for parallel
    cl = makeForkCluster(nnodes = ncores)
    
    # obtain scores for all methods
    results = pblapply(c(1:length(parents)), cl=cl, function(i){
      
      print(i)
      p = as.numeric(yeast_microarray[parents[i],])
      c = as.numeric(yeast_microarray[children[i],])
      
      tblc = cbind(p, c)
      ttblc = cbind(c, p)
      
      # finer discretization; Use unique values as distinct levels (as in Peters.et.al 2010)
      pn = as.factor(p)
      cn = as.factor(c)
      tbl = table(pn, cn)
      
      # Discretize jointly using GridOnCluster for AFC, CE and FOI
      d = discretize.jointly(cbind(p,c), k = c(2:10)) 
      gc()
      
      # prepare contingency table
      p = d$D[,1]
      c = d$D[,2]
      tbljd = table(p, c)
      
      # grab causal versus reverse results
      stats = c(CISC_M(tbl), CE_M(tbljd), FOI_M(tbljd), AFC_M(tbljd), 
                DC_M(tbl), DR_M(tbl), HCR_M(tbl), SCR_M(tbljd), GKT_M(tbljd),
                IGCI_M(tblc, "data"), ANM_M(tblc, "data"))
      stats_t = c(CISC_M(t(tbl)), CE_M(t(tbljd)), FOI_M(t(tbljd)), AFC_M(t(tbljd)), 
                  DC_M(t(tbl)), DR_M(t(tbl)), HCR_M(t(tbl)), SCR_M(t(tbljd)), GKT_M(t(tbljd)),
                  IGCI_M(ttblc, "data"), ANM_M(ttblc, "data"))
      
      return(rbind(stats, stats_t))
    })
    
    # stop parallel clusters
    stopCluster(cl)
    
    # Convert result to data.frame
    results = ListToDataMatrixMethods(results)
    results = as.data.frame(results)
    colnames(results) = methodnames 
    
    # ground truth
    gt = rep(c(1,0), each=(nrow(results)/2))
    ttle = paste0("Run ",j,", ","Pertubed Microarray Yeast","\n",
                  "TF->TG vs TG->TF")
    
    # plot results
    stats = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
    expscores[[j]] = list(stats = stats, scores = results)
  }
  
  return(expscores)
}

# Yeast subset, mean directional accuracy over 10 runs using finer discretization
Yeast_dir_scores_finer = function(expscores){
  
  # accuracy collector
  dir_scores = as.data.frame(matrix(0, nrow=length(expscores), ncol=length(methodnames)))
  ntest = nrow(expscores[[1]]$scores)/2
  
  for(j in 1:nrow(dir_scores)){
    
    dir_test = (unlist(lapply(expscores[[j]]$scores, function(i){
      return(sum(i[1:ntest] > i[(ntest+1):(ntest*2)]))
    }), use.names = FALSE)*100)/ntest
    
    dir_scores[j, ] = dir_test 
  }
  
  ttle = paste0("Directional Accuracy","\n",
                "TF->TG versus TG->TF","\n", 
                "(Finer Discretization)")
  
  # calculate mean accuracy with standard deviation for error bars
  dir_mean = unlist(lapply(dir_scores, mean), use.names = FALSE)
  dir_sd = unlist(lapply(dir_scores, sd), use.names = FALSE)
  
  # order by mean accuracy from best to worst
  ord = order(dir_mean, decreasing = TRUE)
  meths = ordered(methodnames[ord], levels=methodnames[ord])
  dir_scores = data.frame(accu=dir_mean[ord], 
                          meths=meths,
                          sd=dir_sd[ord])
  # plot results
  plt = ggplot(data=dir_scores, aes(x=meths, y=accu)) +
    geom_bar(stat="identity", fill=colr[ord], color="black")+
    ylim(c(0,100)) +
    geom_errorbar(aes(ymin=accu-sd, ymax=accu+sd), width=.2,
                  position=position_dodge(.9)) +
    geom_text(aes(label=paste0(round(accu, digits = 2),"%")), vjust=0.3, 
              hjust=-0.3, color="black", size=7.5, angle=90)+
    theme_bw(base_size = 22) + xlab("Methods") + ylab("Accuracy") +
    ggtitle(ttle) +
    theme(title = element_text(face = "bold", size=21),
          legend.position = "none", 
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))
  print(plt)
}

# Evaluatiing directional inference on a subset of MPAL data using finer discretization
MPAL_func_vs_nfunc_finer = function(ncores){
  
  # Heavy memory usage, thus no parallel implementation
  
  load("../Data/MPAL/normalized_filtered_data.RData")
  
  causal_interactions = interactions$causal_interactions
  
  # protein as parent
  p_candidates = match(causal_interactions[,1], gsub("-.*","",alt_protein_names))
  p_candidates = which(!is.na(p_candidates))
  
  # find protein and rna index
  parent_index = match(causal_interactions[p_candidates,1], 
                       gsub("-.*","",alt_protein_names))
  child_index = match(causal_interactions[p_candidates,2],
                      gene_names)
  
  parent_index = parent_index[!is.na(child_index)]
  child_index = child_index[!is.na(child_index)]
  
  # result collector
  mpal_results_finer = vector("list", 10)
  
  # repeat experiment 10 times
  for(k in 1:10){
    
    # prepare cluster for parallel
    # cl = makeForkCluster(nnodes = ncores)
    
    # collect all scores
    results = pblapply(c(1:length(parent_index)), function(i){
      
      p = as.numeric(norm_proteins[parent_index[i],])
      c = as.numeric(norm_scdata[child_index[i],])
      z = is.na(p) | is.infinite(p) | c==0
      p = p[!z]
      c = c[!z]
      
      # randomly choose 100 single cells for protein/rna pair with more than 100 single cells
      if(length(p) >= 100){
        
        random_sample = sample(1:length(p), 100)
        p = p[random_sample]
        c = c[random_sample]
        
        tblc = cbind(p, c)
        ttblc = cbind(c, p)
        
        # finer discretization; Use unique values as distinct levels (as in Peters.et.al 2010)
        pn = as.factor(p)
        cn = as.factor(c)
        tbl = table(pn, cn)
        
        # Discretize jointly using GridOnCluster for AFC, CE and FOI
        d = discretize.jointly(cbind(p,c), k = c(2:10))
        gc()
        
        # prepare contingency table
        p = d$D[,1]
        c = d$D[,2]
        tbljd = table(p, c)
        
        # grab causal versus reverse results
        stats = c(CISC_M(tbl), CE_M(tbljd), FOI_M(tbljd), AFC_M(tbljd), 
                  DC_M(tbl), DR_M(tbl), HCR_M(tbl), SCR_M(tbljd), GKT_M(tbljd),
                  IGCI_M(tblc, "data"), ANM_M(tblc,  "data"))
        stats_t = c(CISC_M(t(tbl)), CE_M(t(tbljd)), FOI_M(t(tbljd)), AFC_M(t(tbljd)), 
                    DC_M(t(tbl)), DR_M(t(tbl)), HCR_M(t(tbl)), SCR_M(t(tbljd)), GKT_M(t(tbljd)),
                    IGCI_M(ttblc, "data"), ANM_M(ttblc, "data"))
        
        return(rbind(stats, stats_t))
        
      } else { # return worst scores for both directions
        return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                     rep(-.Machine$integer.max, length(methodnames))))
      }
      
    })
    
    # stop parallel clusters
    # stopCluster(cl)
    
    # Convert result to data.frame
    results = ListToDataMatrixMethods(results)
    results = as.data.frame(results)
    colnames(results) = methodnames 
    
    # ground truth
    gt = rep(c(1,0), each=(nrow(results)/2))
    ttle = paste0("Run ",k,", Single Cell Leukemia Data","\n",
                  "Causal Interactions vs Reverse")
    # plot results
    stats = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
    mpal_results_finer[[k]] = list(stats = stats, scores = results)
  }
  
  # return results
  return(mpal_results_finer)  
  
}

# MPAL subset, mean directional accuracy over 10 runs using finer discretization
MPAL_dir_scores_finer = function(expscores){
  
  # accuracy collector
  dir_scores = as.data.frame(matrix(0, nrow=length(expscores), ncol=length(methodnames)))
  ntest = nrow(expscores[[1]]$scores)/2
  
  for(j in 1:nrow(dir_scores)){
    
    dir_test = (unlist(lapply(expscores[[j]]$scores, function(i){
      return(sum(i[1:ntest] > i[(ntest+1):(ntest*2)]))
    }), use.names = FALSE)*100)/ntest
    
    dir_scores[j, ] = dir_test 
  }
  
  ttle = paste0("Directional Accuracy","\n",
                "Protein->RNA versus","\n", 
                "RNA->Protein (Finer Discretization)")
  
  # calculate mean accuracy with standard deviation for error bars
  dir_mean = unlist(lapply(dir_scores, mean), use.names = FALSE)
  dir_sd = unlist(lapply(dir_scores, sd), use.names = FALSE)
  
  # order by mean accuracy from best to worst
  ord = order(dir_mean, decreasing = TRUE)
  meths = ordered(methodnames[ord], levels=methodnames[ord])
  dir_scores = data.frame(accu=dir_mean[ord], 
                          meths=meths,
                          sd=dir_sd[ord])
  # plot results
  plt = ggplot(data=dir_scores, aes(x=meths, y=accu)) +
    geom_bar(stat="identity", fill=colr[ord], color="black")+
    ylim(c(0,100)) +
    geom_errorbar(aes(ymin=accu-sd, ymax=accu+sd), width=.2,
                  position=position_dodge(.9)) +
    geom_text(aes(label=paste0(round(accu, digits = 2),"%")), vjust=0.3, 
              hjust=-0.1, color="black", size=7.5, angle=90)+
    theme_bw(base_size = 22) + xlab("Methods") + ylab("Accuracy") +
    ggtitle(ttle) +
    theme(title = element_text(face = "bold", size=21),
          legend.position = "none", 
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))
  print(plt)
}

# MPAL subset, top known patterns reported by all methods using finer discretization
PlotTopKnownPatternsFiner= function(expscores){
  
  load("../Data/MPAL/normalized_filtered_data.RData")
  
  causal_interactions = interactions$causal_interactions
  
  # protein as parent
  p_candidates = match(causal_interactions[,1], gsub("-.*","",alt_protein_names))
  p_candidates = which(!is.na(p_candidates))
  
  # find protein and rna index
  parent_index = match(causal_interactions[p_candidates,1], 
                       gsub("-.*","",alt_protein_names))
  child_index = match(causal_interactions[p_candidates,2],
                      gene_names)
  
  parent_index = parent_index[!is.na(child_index)]
  child_index = child_index[!is.na(child_index)]
  
  # number of pairs
  ntest = nrow(expscores[[1]]$scores)/2
  
  # Plot top known causal patterns for each method
  for(i in 1:length(methodnames)){
    
    # get method scores for a random run
    meth_scores = expscores[[sample(1:10,1)]]$scores[,i]
    
    # pick one direction
    eff_scores = sapply(c(1:ntest), function(k){
      return(ifelse(meth_scores[k] > meth_scores[k+ntest], meth_scores[k], meth_scores[k+ntest]))
    })
    
    # is the picked direction correct?
    direction = ifelse(meth_scores[1:ntest] > meth_scores[(ntest+1):(2*ntest)], 1, -1)
    
    # rank the picked direction (Note: methods in Methods.R return scores that are higher the better)
    meth_score_ord = order(eff_scores, decreasing = TRUE)
    
    # plot top 3
    for(j in 1:3){
      
      png(paste0("../Results/Supplementary/MPAL_Finer/TopKnonwPatternsFiner/",methodnames[i],"_",j,".png"))
      
      # get effective index
      ind = meth_score_ord[j]
      dir = direction[ind]
      
      if(dir == -1){
        
        p = as.numeric(norm_proteins[parent_index[ind],])
        c = as.numeric(norm_scdata[child_index[ind],])
        
        z = is.na(p) | is.infinite(p) | c==0
        p = p[!z]
        c = c[!z]
        
        # randomly select 100 single cells, if protein/rna have more than 100 single cells
        if(length(p) >= 100){
          random_sample = sample(1:length(p), 100)
          p = p[random_sample]
          c = c[random_sample]
        }
        
        clr = RColorBrewer::brewer.pal("RdYlGn", n=6)[6]
        
        # plot pattern
        data = as.data.frame(cbind(c,p))
        plt = 
          ggplot(data, aes(x=c, y=p) ) +
          geom_point(color=clr, pch=21, lwd=3, fill="black", stroke=3) +
          coord_cartesian(ylim=c(min(p),max(p))) +
          geom_smooth(linetype="solid", color="firebrick4", se=FALSE,
                      method="glm", lwd=2, formula = y ~ x + I(x^2) + I(x^3))+
          scale_fill_distiller(palette = "RdYlGn", direction = -1) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_bw(base_size = 24) + ylab(paste0("PROT:",alt_protein_names[parent_index[ind]])) +
          xlab(paste0("RNA:",gene_names[child_index[ind]])) +
          theme(panel.background = element_rect(fill = 'white'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.title = element_text(face = "bold", size = 40))#,
        print(plt)
        
      } else {
        
        p = as.numeric(norm_proteins[parent_index[ind],])
        c = as.numeric(norm_scdata[child_index[ind],])
        
        z = is.na(p) | is.infinite(p) | c==0
        p = p[!z]
        c = c[!z]
        
        # randomly select 100 single cells, if protein/rna have more than 100 single cells
        if(length(p) >= 100){
          random_sample = sample(1:length(p), 100)
          p = p[random_sample]
          c = c[random_sample]
        }
        
        clr = RColorBrewer::brewer.pal("RdYlGn", n=6)[6]
        
        # plot patterns
        data = as.data.frame(cbind(p,c))
        plt = 
          ggplot(data, aes(x=p, y=c) ) +
          geom_point(color=clr, pch=21, lwd=3, fill="black", stroke=3) +
          coord_cartesian(ylim=c(min(c),max(c))) +
          geom_smooth(linetype="solid", color="firebrick4", se=FALSE,
                      method="glm", lwd=2, formula = y ~ x + I(x^2) + I(x^3))+
          scale_fill_distiller(palette = "RdYlGn", direction = -1) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_bw(base_size = 24) + xlab(paste0("PROT:",alt_protein_names[parent_index[ind]])) +
          ylab(paste0("RNA:",gene_names[child_index[ind]])) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.title = element_text(face = "bold", size = 40))#,
        print(plt)
        
      }
      
      dev.off()
    }
  }
}

# Evaluating directional inference on a subset of Abalone data using finer discretization
abalone_func_vs_nfunc_finer = function(ncores){
  
  
  abalone = read.table("../Data/Abalone/abalone.data", sep=",")
  abalone = abalone[,c(1:4)]
  colnames(abalone) = c("Sex","Length","Diameter","Height")
  
  # Sex -> {Length, Diameter, Height}
  sm = 100 # we will be choosing 100 samples randomly
  repeats = 100 # sample the 3 cause->effect pair 100 times
  
  # Sex -> {Length, Diameter, Height}
  gt = rep(c(1, 0, 1, 0, 1, 0), repeats)
  
  # result collector
  expscores = vector("list", length = 10)
  expstats = vector("list", length = 10)
  
  # repeat the experiment 10 times
  for(i in 1:10){
    
    # prepare parallel clusters
    cl = makeForkCluster(nnodes = ncores)
    
    # compute results for all methods
    results = pblapply(1:(3*repeats), cl=cl, function(r){
      
      if(r%%3==0){ # 100 instances of the 3 cause->effect pairs
        c=3
      } else {
        c=r%%3
      }
      
      set.seed(123*r*i)
      
      # randomly select 100 samples
      abalone_s = abalone[sample(1:nrow(abalone), sm),]
      
      # discretize data jointly (Length, Diameter and Height are discretized together)
      # for AFC, CE and FOI
      d = discretize.jointly(log1p(apply(abalone_s[,-1], c(1,2), as.numeric)), k=c(2:10))
      
      x = as.factor(abalone_s[,1])
      y = log1p(as.numeric(abalone_s[,c+1])) # ignoring Sex -> Sex pair
      
      tblc = cbind(as.numeric(x), y)
      ttblc = cbind(y, as.numeric(x))
      
      # finer discretization; Use unique values as distinct levels (as in Peters.et.al 2010)
      pn = x
      cn = as.factor(y)
      tbl = table(pn, cn)
      
      # prepare contingency table
      y = d$D[,c]
      tbljd = table(x, y)
      
      # grab causal versus reverse results
      stats = c(CISC_M(tbl), CE_M(tbljd), FOI_M(tbljd), AFC_M(tbljd), 
                DC_M(tbl), DR_M(tbl), HCR_M(tbl), SCR_M(tbljd), GKT_M(tbljd),
                IGCI_M(tblc, "data"), ANM_M(tblc, "data"))
      stats_t = c(CISC_M(t(tbl)), CE_M(t(tbljd)), FOI_M(t(tbljd)), AFC_M(t(tbljd)), 
                  DC_M(t(tbl)), DR_M(t(tbl)), HCR_M(t(tbl)), SCR_M(t(tbljd)), GKT_M(t(tbljd)),
                  IGCI_M(ttblc, "data"), ANM_M(ttblc, "data"))
      
      return(rbind(stats, stats_t))
    
    })
    
    # stop parallel cluster
    stopCluster(cl)
    
    # convert results to dara.frame
    results = ListToDataMatrixMethods(results)
    results = as.data.frame(results)
    colnames(results) = methodnames
    
    ttle = paste0("Run ",i," Abalone FVNF")
    
    # collect and plot results
    expscores[[i]] = results
    expstats[[i]] = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
  }
  
  return(list(expscores = expscores, expstats = expstats))
}

# Abalone subset, mean directional accuracy over 10 runs using finer discretization
abalone_dir_scores_finer = function(expscores){
  
  dir_scores = as.data.frame(matrix(0, nrow=length(expscores), ncol=length(methodnames)))
  ntest = nrow(expscores[[1]])/2
  
  for(j in 1:nrow(dir_scores)){
    
    dir_test = (unlist(lapply(expscores[[j]], function(i){
      return(sum(i[1:ntest] > i[(ntest+1):(ntest*2)]))
    }), use.names = FALSE)*100)/ntest
    
    dir_scores[j, ] = dir_test 
  }
  
  ttle = paste0("Directional Accuracy","\n",
                "Abalone","\n", 
                "(Finer Discretization)")
  
  # calculate mean accuracy with standard deviation for error bars
  dir_mean = unlist(lapply(dir_scores, mean), use.names = FALSE)
  dir_sd = unlist(lapply(dir_scores, sd), use.names = FALSE)
  
  # order by mean accuracy from best to worst
  ord = order(dir_mean, decreasing = TRUE)
  meths = ordered(methodnames[ord], levels=methodnames[ord])
  dir_scores = data.frame(accu=dir_mean[ord], 
                          meths=meths,
                          sd=dir_sd[ord])
  # plot results
  plt = ggplot(data=dir_scores, aes(x=meths, y=accu)) +
    geom_bar(stat="identity", fill=colr[ord], color="black")+
    ylim(c(0,150)) +
    geom_errorbar(aes(ymin=accu-sd, ymax=accu+sd), width=.2,
                  position=position_dodge(.9)) +
    geom_text(aes(label=paste0(round(accu, digits = 2),"%")), vjust=0.3, 
              hjust=-0.3, color="black", size=7.5, angle=90)+
    theme_bw(base_size = 22) + xlab("Methods") + ylab("Accuracy") +
    ggtitle(ttle) +
    theme(title = element_text(face = "bold", size=21),
          legend.position = "none", 
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))
  print(plt)
}

# Get runtime for each method
run_experiment_runtime = function(tbl){
  
  rt = rep(0, length(methodnames))
  
  rt[1] = system.time(CISC_M(tbl))[3]/2 # calculated twice
  rt[2] = system.time(CE_M(tbl))[3]
  rt[3] = system.time(FOI_M(tbl))[3]
  rt[4] = system.time(AFC_M(tbl))[3]
  rt[5] = system.time(DC_M(tbl))[3]/2 # calculated twice
  rt[6] = system.time(DR_M(tbl))[3]/2 # calculated twice
  rt[7] = system.time(HCR_M(tbl))[3]
  rt[8] = system.time(SCR_M(tbl))[3]/2 # calculated twice
  rt[9] = system.time(GKT_M(tbl))[3]
  rt[10] = system.time(IGCI_M(tbl))[3]/2 # calculated twice
  rt[11] = system.time(ANM_M(tbl))[3]
  
  return(rt)
}

# Runtime over sample and table size
Runtime_eval = function(){
  
  # table sizes
  tsize = seq(10, 100, 10)
  
  # simulate independent patterns -- by table size
  tbls = vector("list", length = length(tsize))
  for(i in 1:length(tsize)){
    tbls[[i]] = simulate_tables(n=10000, nrow=tsize[i], ncol=(tsize[i]+i), type="independent")$sample.list[[1]]
  }
  
  # collect time over table size
  tsize_time = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(tsize)))
  for(i in 1:length(tsize)){
    tsize_time[i,] = run_experiment_runtime(tbls[[i]])  
  }
  
  ########################################################################################################
  
  # sample size
  ssize = c(10, 50, 100, 500, 1000, 5000, 10000, 20000)
  
  # simulate independent patterns -- by sample size
  tbls = vector("list", length = length(ssize))
  for(i in 1:length(ssize)){
    tbls[[i]] = simulate_tables(n=ssize[i], nrow=10, ncol=15, type="independent")$sample.list[[1]]
  }
  
  # collect time over table size
  ssize_time = as.data.frame(matrix(0, ncol=length(methodnames), nrow=length(ssize)))
  for(i in 1:length(ssize)){
    ssize_time[i,] = run_experiment_runtime(tbls[[i]])  
  }
  
  ########################################################################################################
  return(list(tsize_time, ssize_time))
}

PlotRuntime = function(tsize_time, ssize_time){
  
  # table sizes
  tsize = seq(10, 100, 10)
  
  # sample size
  ssize = c(10, 50, 100, 500, 1000, 5000, 10000, 20000)
  
  ltyps2 = rep(2, length(methodnames))
  ltyps2[which(methodnames == "AFC")] = 1
  
  lwds = rep(3, length(methodnames))
  lwds[which(methodnames == "AFC")] = 4
  
  pchs = c(2, 3, 4, 19, 5, 6, 8, 9, 11, 13, 18)
  
  ttle = paste0("Runtime over increasing table size","\n",
                "Independent patterns","\n",
                "Sample Size = 10,000")
  # plot results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="Log+1 Runtime", main=ttle, xaxt="n", xlim=c(min(tsize), max(tsize)),
       ylim = c(log1p(min(tsize_time)), log1p(max(tsize_time))), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  
  for(k in 1:length(methodnames)){
    lines(log1p(tsize_time[,k])~tsize, col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(log1p(tsize_time[,k])~tsize, col=colr[k], pch=pchs[k], cex=2)
  }  
  
  axis(side=1, at=tsize, labels = paste0((tsize),"x",(tsize+c(1:length(tsize)))), las=2, cex.axis=1.5)
  mtext("Table size", side=1, line=5.5, cex=2)
  
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
  ########################################################################################################
  
  ttle = paste0("Runtime over increasing sample size","\n",
                "Independent patterns","\n",
                "Table Size = 10 x 15")
  # plot results
  par(mar=c(8,5,8,8), lwd=2)
  plot(0, xlab="", ylab="Log+1 Runtime", main=ttle, xaxt="n", xlim=c(1, length(ssize)),
       ylim = c(log1p(min(ssize_time)), log1p(max(ssize_time))), cex.lab=2, cex.main=2, type="n", cex.axis=1.5)
  grid()
  
  for(k in 1:length(methodnames)){
    lines(log1p(ssize_time[,k])~c(1:length(ssize)), col=colr[k], lwd=lwds[k], lty=ltyps2[k])
    points(log1p(ssize_time[,k])~c(1:length(ssize)), col=colr[k], pch=pchs[k], cex=2)
  }  
  
  axis(side=1, at=c(1:length(ssize)), labels = ssize, las=2, cex.axis=1.5)
  mtext("Sample size", side=1, line=5.5, cex=2)
  
  par(xpd=TRUE)
  legend("topright", legend = methodnames, col = colr, lty=ltyps2, lwd=lwds,
         pch = pchs, inset=c(-0.37,0), bty = "n", cex=1.5)
  par(xpd=FALSE)
  
}

# Main procedure
Suppl_Plots = function(){
  
  # plot functional tables
  show_tables()
  
  # Detailed summary to study biases -- functional versus independent (AUROC)
  Summarize_simulation_study_FVI()
  
  # Detailed summary to study biases -- functional versus non-functional (AUROC)
  Summarize_simulation_study_FVNF()
  
  # Detailed summary to study biases -- functional versus non-functional (DIR)
  Summarize_simulation_study_FVNF_DIR()
  
  # Plot functional versus non-functional directional scores over an increasing decsision rate
  pdf("../Results/Supplementary/SimStudy_FVNF/FVNF_DIR_ACCU_DIST_DRATE.pdf", width=10)
  load("../Results/SimStudy/fvi200.RData")
  load("../Results/SimStudy/fvi2000.RData")
  load("../Results/SimStudy/fvnf200.RData")
  load("../Results/SimStudy/fvnf2000.RData")
  SimStudy_FVNF_dir_scores_drate(list(fvnf2000=fvnf2000, fvnf200=fvnf200))
  dev.off()
  
  # A snapshot of the simulated data
  PlotSimData()
  
  # get runtime over table and sample size
  timeres = Runtime_eval()
  save(timeres, file="../Results/Supplementary/runtime_res.Rds")
  
  # plot runtime
  pdf("../Results/Supplementary/runtime_eval.pdf")
  PlotRuntime(timeres[[1]], timeres[[2]])
  dev.off()
  
  # Evaluate causal versus non-causal on Yeast subset using finer discretization
  if(!file.exists("../Data/Yeast/yeast_microarray_cleaned_log_ratio.RData")){
    Prepare_data()
  }
  pdf("../Results/Supplementary/Yeast_Finer/yeast_fvnf_roc_pr_finer.pdf")
  yeast_fvnf_finer = yeast_data_func_vs_nfunc_finer(ncores=10)
  dev.off()
  
  # Plot mean directional accuracy on Yeast subset using finer discretization
  pdf("../Results/Supplementary/Yeast_Finer/yeast_dir_test_finer.pdf")
  Yeast_dir_scores_finer(yeast_fvnf_finer)
  dev.off()
  save(yeast_fvnf_finer, 
       file="../Results/Supplementary/Yeast_Finer/yeast_fvnf_finer.RData")
  
  # Evaluate causal versus non-causal on MPAL subset using finer discretization
  if(!file.exists("../Data/MPAL/normalized_filtered_data.RData")){
    preprocess_data(PLOT=FALSE)
  }
  pdf("../Results/Supplementary/MPAL_Finer/mpal_fvnf_roc_pr_finer.pdf")
  mpal_fvnf_finer = MPAL_func_vs_nfunc_finer()
  dev.off()
  
  # Plot mean directional accuracy on MPAL subset using finer discretization
  pdf("../Results/Supplementary/MPAL_Finer/mpal_dir_test_finer.pdf")
  MPAL_dir_scores_finer(mpal_fvnf_finer)
  dev.off()
  save(mpal_fvnf_finer, 
       file="../Results/Supplementary/MPAL_Finer/mpal_fvnf_finer.RData")
  
  # Plot top known pattern using finer discretization on MPAL subset
  PlotTopKnownPatternsFiner(mpal_fvnf_finer)
  
  # Evaluate causal versus non-causal on Abalone subset using finer discretization
  pdf("../Results/Supplementary/Abalone_Finer/abalone_fvnf_roc_pr_finer.pdf")
  abalone_fvnf_finer = abalone_func_vs_nfunc_finer(ncores=10)
  dev.off()
  
  # Plot mean directional accuracy on Abalone subset using finer discretization
  pdf("../Results/Supplementary/Abalone_Finer/abalone_dir_test_finer.pdf")
  abalone_dir_scores_finer(abalone_fvnf_finer$expscores)
  dev.off()
  save(abalone_fvnf_finer, 
       file="../Results/Supplementary/Abalone_Finer/abalone_fvnf_finer.RData")
  
  
}