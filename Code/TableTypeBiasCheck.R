# Table Type Bias Check setup for:
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



Simulate_noisy_functional_tables = function(n, r, s, noise, ntables, ncores){
  
  tbs = vector("list", ntables)
  
  cl = makeForkCluster(nnodes = ncores)
  tbs = pblapply(cl = cl, X = 1:ntables, FUN = function(i){
    
    set.seed(round(as.numeric(Sys.time())/i))
    tb = simulate_tables(n = n, nrow = r, ncol = s, type="functional")$sample.list[[1]]
    tb = add.candle.noise(tb, u=noise, margin = 0)
    
  })
  stopCluster(cl)
  
  return(tbs)
}



TableTypeCheck = function(N, r, s, ncores){
  
  noise = c(0.01,seq(0.1,1,0.1))
  stat_collect = vector("list", length = length(noise))
  
  for(n in 1:length(noise)){
    
    print(paste0("Noise: ",noise[n]))
    tbs = Simulate_noisy_functional_tables(n = N, r = r, s = s, noise = noise[n],
                                           ntables = 1000, ncores = 10)
    
    cl = makeForkCluster(nnodes = ncores)
    meth_scores1 = pblapply(cl = cl, X=1:length(tbs), FUN = function(tb){
      run_experiment(tbs[[tb]])
    })
    stopCluster(cl)
    
    meth_scores1 = ListToDataMatrixR(meth_scores1)
    
    tbs = Simulate_noisy_functional_tables(n = N, r = s, s = r, noise = noise[n],
                                           ntables = 1000, ncores = 10)
    
    cl = makeForkCluster(nnodes = ncores)
    meth_scores2 = pblapply(cl = cl, X=1:length(tbs), FUN = function(tb){
      run_experiment(tbs[[tb]])
    })
    stopCluster(cl)
    
    meth_scores2 = ListToDataMatrixR(meth_scores2)
    
    colnames(meth_scores1) = methodnames
    colnames(meth_scores2) = methodnames
    
    stat_collect[[n]] = list(meth_scores1, meth_scores2)
  }
  
  return(stat_collect)
}

PlotPatterns = function(stat_collect, n, r, s){
  
  noise = c(0.01,seq(0.1,1,0.1))
  
  methform = c("-CISC", "-CE", "FOI", "-ln(AFC p-value)", "-DC", "ln (DR p-value)", "HCR",
               "SCR", "GKT", "-IGCI", "-ANM")
  
  par(mfrow=c(2, 6))
  par(oma = c(0, 1, 0, 0))
  par(mar = c(3.5, 4.5, 3, 1))
  
  for(m in 1:length(methodnames)){
    
    res_list = vector("list", length = length(noise) * 3)
    prtrt_med = rep(0, length(noise))
    lndscp_med = rep(0, length(noise))
    
    for(i in 1:length(noise)){
      
      if(m == 8){
        
        val1 = stat_collect[[i]][[1]][,m]
        if(any(val1 < 0)){
          val1[val1 < 0] = 0
        }
        res_list[[(i*3)-2]] = val1
        val2 = stat_collect[[i]][[2]][,m]
        if(any(val2 < 0)){
          val2[val2 < 0] = 0
        }
        res_list[[(i*3)-1]] = val2
        res_list[[i*3]] = NA
        
        prtrt_med[i] = median(val1)
        lndscp_med[i] = median(val2)
        
      } else {
        
        res_list[[(i*3)-2]] = stat_collect[[i]][[1]][,m]
        res_list[[(i*3)-1]] = stat_collect[[i]][[2]][,m]
        res_list[[i*3]] = NA
        
        prtrt_med[i] = median(stat_collect[[i]][[1]][,m])
        lndscp_med[i] = median(stat_collect[[i]][[2]][,m])
      }
    
    }  
    
    boxplot(res_list, col = c("blue","red","white"), lwd=2, pch=5, xaxt="n", 
            ylab=methform[m], cex.main=2, main=methodnames[m], 
            outcol=alpha("black",0.2), cex.lab=1.5, font.lab=2)
    grid()
    lines(prtrt_med~seq(1,31,3), col="blue", lwd=2, lty=1)
    lines(lndscp_med~seq(2,32,3), col="red", lwd=2, lty=1)
    axis(side = 1, at = seq(1, 31, 3)+0.5, labels = noise, las=2)
    
  }
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("center", legend =c("Table Type Bias", 
                             "X: Noise levels", 
                             "Y: Method Statistics", 
                             paste0("N: ",n),
                             "",
                             paste0("Portrait table (",r,"x",s,")"), 
                             paste0("Landscape table (",s,"x",r,")")), 
         pch=c(NA, NA, NA, NA, NA, 19, 19), pt.cex=3, cex=1.2, 
         bty='n', col = c(NA, NA, NA, NA, NA, 'blue', 'red'))
  
}


ListToDataMatrixR = function(data_list){
  
  new_data = matrix(0, ncol=length(data_list[[1]]), nrow=length(data_list))
  data = unlist(data_list, use.names = FALSE)
  
  for(i in 1:ncol(new_data)){
    access = rep(FALSE, ncol(new_data))
    access[i] = TRUE 
    new_data[,i] = data[access]
  }
  
  return(new_data)
}

TtyBiasCheck = function(){
  
  # 10 x 3 versus 3 x 10 at n = 100
  res_10x3_3x10_100 = TableTypeCheck(N = 100, r = 10, s = 3, ncores = 10)
  save(res_10x3_3x10_100, file="../Results/TableTypeStudy/res_10x3_3x10_100.RData")
  
  pdf("../Results/TableTypeStudy/tty_10x3_3x10_100.pdf", width = 15)
  PlotPatterns(stat_collect = res_10x3_3x10_100, n = 100, r = 10, s = 3)
  dev.off()
  
}