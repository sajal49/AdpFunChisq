# Simulation study setup for:
# Causal Inference by Functional and Statistical Dependency

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


# Generate synthetic population of a particular type / marginal structure
gen_population = function(nrow, ncol, n, ntables, type, rmar=1, cmar=1){
  
  # adjust row / column marginal according to parameters
  # 1 means uniform
  # NULL means no control
  # Any other number denotes the "amount" of non-uniformity
  
  if(is.null(rmar)){
    rowmar = NULL
  } else if(rmar == 1){
    rowmar = rep(1/nrow, nrow) # uniform
  } else if(rmar > 1) {
    rowmar = (c(1:nrow)**rmar) # amount of non-uniformity in rows
    rowmar = rowmar/sum(rowmar)
  } else {
    stop("Bad input for rmar!!")
  }
  
  if(is.null(cmar)){
    colmar = NULL
  } else if(cmar == 1){
    colmar = rep(1/ncol, ncol) # uniform
  } else if(cmar > 1){
    colmar = (c(1:ncol)**cmar) # amount of non-uniformity in columns
    colmar = colmar/sum(colmar)
  } else {
    stop("Bad input for cmar!!")
  }
  
  # base population with no noise
  tbls = simulate_tables(n=n, nrow=nrow, ncol=ncol, type=type, 
                         row.marginal = rowmar, col.marginal = colmar,
                         n.tables = ntables, noise = 0)$sample.list  
 
  # generate population
  return(tbls)
  
}

# Functional versus independent population
func_vs_indep_pop = function(nrow, ncol, n, ntables, noise_levels, margin){
  
  # functional is fixed at non-uniform row with no control on column marginal,
  # however, with non-uniform row, most tables will have non-uniform column marginals.
  
  # experiment 1
  # uniform row marginal and non-uniform column marginal independent
  exp1fun = vector("list", length = length(noise_levels))
  exp1ind = vector("list", length = length(noise_levels))
  
  # generate base tables
  base1fun = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "functional", rmar = 2, cmar = NULL)
  base1ind = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "independent", rmar = 1, cmar = 2)
  
  # add house noise to base tables
  for(i in 1:length(noise_levels)){
    exp1fun[[i]] = add.house.noise(base1fun, u=noise_levels[i], margin = margin)
    exp1ind[[i]] = base1ind
  }
  
  # ground truth
  gt1 = c(rep(1, ntables), rep(0, ntables)) # higher the better
  
  # experiment 2
  # non-uniform row marginal and non-uniform column marginal independent
  exp2fun = vector("list", length = length(noise_levels))
  exp2ind = vector("list", length = length(noise_levels))
  
  # generate base tables
  base2fun = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "functional", rmar = 2, cmar = NULL)
  base2ind = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "independent", rmar = 2, cmar = 2)
  
  # add house noise to base tables
  for(i in 1:length(noise_levels)){
    exp2fun[[i]] = add.house.noise(base2fun, u=noise_levels[i], margin = margin)
    exp2ind[[i]] = base2ind
  }
  
  # ground truth
  gt2 = c(rep(1, ntables), rep(0, ntables))
  
  # experiment 3
  # non-uniform row marginal and uniform column marginal independent
  exp3fun = vector("list", length = length(noise_levels))
  exp3ind = vector("list", length = length(noise_levels))
  
  # generate base tables
  base3fun = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "functional", rmar = 2, cmar = NULL)
  base3ind = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "independent", rmar = 2, cmar = 1)
  
  # add house noise to base tables
  for(i in 1:length(noise_levels)){
    exp3fun[[i]] = add.house.noise(base3fun, u=noise_levels[i], margin = margin)
    exp3ind[[i]] = base3ind
  }
  
  # ground truth
  gt3 = c(rep(1, ntables), rep(0, ntables))
  
  # return setup
  return(list(exp1fun = exp1fun, exp2fun = exp2fun, exp3fun = exp3fun, exp1ind=exp1ind, exp2ind=exp2ind,
              exp3ind = exp3ind, gt1=gt1, gt2=gt2, gt3=gt3))
}

# Functional versus non-functional population
func_vs_nfunc_pop = function(nrow, ncol, n, ntables, noise_levels, margin){
  
  # experiment 1
  # uniform row functional vs transpose (no control on column)
  exp1fun = vector("list", length = length(noise_levels))
  exp1nfun = vector("list", length = length(noise_levels))
  
  # generate base tables
  base1fun = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "many.to.one", rmar = 1, cmar = NULL)
  
  # add house noise to base tables
  for(i in 1:length(noise_levels)){
    exp1fun[[i]] = add.house.noise(base1fun, u=noise_levels[i], margin = margin)
    exp1nfun[[i]] = lapply(exp1fun[[i]], t)
  }
  
  # ground truth
  gt1 = c(rep(1, ntables), rep(0, ntables))
  
  # experiment 2
  # uniform column functional vs transpose (no control on row)
  exp2fun = vector("list", length = length(noise_levels))
  exp2nfun = vector("list", length = length(noise_levels))
  
  # generate base tables
  base2fun = gen_population(nrow = nrow, ncol = ncol, n = n, ntables = ntables, 
                            type = "many.to.one", rmar = NULL, cmar = 1)
  
  # add noise to base tables
  for(i in 1:length(noise_levels)){
    exp2fun[[i]] = add.house.noise(base2fun, u=noise_levels[i], margin = margin)
    exp2nfun[[i]] = lapply(exp2fun[[i]], t)
  }
  
  # ground truth
  gt2 = c(rep(1, ntables), rep(0, ntables))
  
  # return setup
  return(list(exp1fun = exp1fun, exp2fun = exp2fun, exp1nfun=exp1nfun, exp2nfun=exp2nfun,
              gt1=gt1, gt2=gt2))
}

# Functional vs independent experiment
func_vs_indep = function(n, ncores){
  
  # 6 table sizes x 5 noise levels x 3 marginal schemes 
  # noise levels
  noise_levels = c(0.01, 0.25, 0.5, 0.75, 1)
  
  # table sizes
  r = c(10, 7, 9, 3, 5, 8)
  s = c(8, 5, 9, 5, 7, 10)
  
  # result containers
  exp1scores = vector("list", length = length(r))
  exp1stats = vector("list", length = length(r))
  exp2scores = vector("list", length = length(r))
  exp2stats = vector("list", length = length(r))
  exp3scores = vector("list", length = length(r))
  exp3stats = vector("list", length = length(r))
  populations = vector("list", length = length(r))
  
  # begin experiment
  for(i in 1:length(r)){
    
    cat("Functional vs Independent","\n",
        "Table size=",r[i],"x",s[i],"\n",
        "Sample size=",n,"\n", sep = "")
    population = func_vs_indep_pop(nrow = r[i], ncol = s[i], n = n, ntables = 200, 
                                   noise_levels = noise_levels, margin = 0) # generate population
    populations[[i]] = population
    
    ########################################################################################################
    
    # experiment 1
    cat("Experiment 1","\n",
        "non-unif row/col functional vs unif row and non-unif col independent", "\n",
        sep = "")
    exp1score = vector("list", length = length(noise_levels))
    exp1stat = vector("list", length = length(noise_levels))
    
    for(nc in 1:length(noise_levels)){
      cat("At noise level=",noise_levels[nc],"\n")
      
      cl = makeForkCluster(nnodes = ncores)
      exp1score[[nc]] = pblapply(cl = cl, X = 1:length(population$exp1fun[[nc]]), FUN = function(j){
        return(rbind(run_experiment(population$exp1fun[[nc]][[j]]),
                     run_experiment(population$exp1ind[[nc]][[j]])))
      }) # grab results
      stopCluster(cl)
      
      scores = as.data.frame(ListToDataMatrixMethods(exp1score[[nc]]))
      colnames(scores) = methodnames
      ttle = paste0("Functional vs Independent","\n",
                    "NU Mar Func vs UR NUC Indep","\n",
                    "S=",n,", T=",r[i],"x",s[i],"\n",
                    "Noise=",noise_levels[nc])
      # plot ROC/PR curves
      exp1stat[[nc]] = ggplot.ROC.PR.curves(scores, population$gt1, colr, ltyps, plot=TRUE, ttle)
    }
    
    # store results
    exp1scores[[i]] = exp1score
    exp1stats[[i]] = exp1stat
    
    ########################################################################################################
    
    # experiment 2
    cat("Experiment 2","\n",
        "non-unif row/col functional vs non-unif row/col independent", "\n",
        sep = "")
    exp2score = vector("list", length = length(noise_levels))
    exp2stat = vector("list", length = length(noise_levels))
    
    for(nc in 1:length(noise_levels)){
      cat("At noise level=",noise_levels[nc],"\n")
      
      cl = makeForkCluster(nnodes = ncores)
      exp2score[[nc]] = pblapply(cl = cl, X = 1:length(population$exp2fun[[nc]]), FUN = function(j){
        return(rbind(run_experiment(population$exp2fun[[nc]][[j]]),
                     run_experiment(population$exp2ind[[nc]][[j]])))
      }) # grab results
      stopCluster(cl)
      
      scores = as.data.frame(ListToDataMatrixMethods(exp2score[[nc]]))
      colnames(scores) = methodnames
      ttle = paste0("Functional vs Independent","\n",
                    "NU Mar Func vs NU Mar Indep","\n",
                    "S=",n,", T=",r[i],"x",s[i],"\n",
                    "Noise=",noise_levels[nc])
      # plot ROC/PR
      exp2stat[[nc]] = ggplot.ROC.PR.curves(scores, population$gt2, colr, ltyps, plot=TRUE, ttle)
    }
    
    # store results
    exp2scores[[i]] = exp2score
    exp2stats[[i]] = exp2stat
    
    ########################################################################################################
    
    # experiment 3
    cat("Experiment 3","\n",
        "non-unif row/col functional vs non-unif row and unif col independent", "\n",
        sep = "")
    exp3score = vector("list", length = length(noise_levels))
    exp3stat = vector("list", length = length(noise_levels))
    
    for(nc in 1:length(noise_levels)){
      cat("At noise level=",noise_levels[nc],"\n")
      
      cl = makeForkCluster(nnodes = ncores)
      exp3score[[nc]] = pblapply(cl = cl, X = 1:length(population$exp3fun[[nc]]), FUN = function(j){
        return(rbind(run_experiment(population$exp3fun[[nc]][[j]]),
                     run_experiment(population$exp3ind[[nc]][[j]])))
      }) # grab results
      stopCluster(cl)
      
      scores = as.data.frame(ListToDataMatrixMethods(exp3score[[nc]]))
      colnames(scores) = methodnames
      ttle = paste0("Functional vs Independent","\n",
                    "NU Mar Func vs NUR UC Indep","\n",
                    "S=",n,", T=",r[i],"x",s[i],"\n",
                    "Noise=",noise_levels[nc])
      # plot results
      exp3stat[[nc]] = ggplot.ROC.PR.curves(scores, population$gt3, colr, ltyps, plot=TRUE, ttle)
    }
    
    # store results
    exp3scores[[i]] = exp3score
    exp3stats[[i]] = exp3stat
    
    ########################################################################################################
  }
  
  # return population / result containers
  return(list(exp1scores = exp1scores, exp2scores = exp2scores, exp3scores = exp3scores,
              exp1stats = exp1stats, exp2stats = exp2stats, exp3stats = exp3stats,
              populations = populations))
}

# Functional vs non-functional experiment
func_vs_nfunc = function(n, ncores){
  
  # 6 table sizes x 5 noise levels x 2 marginal schemes
  noise_levels = c(0.01, 0.25, 0.5, 0.75, 1)
  
  # table sizes
  r = c(10, 7, 9, 3, 5, 8)
  s = c(8, 5, 9, 5, 7, 10)
  
  # result containers
  exp1scores = vector("list", length = length(r))
  exp1stats = vector("list", length = length(r))
  exp2scores = vector("list", length = length(r))
  exp2stats = vector("list", length = length(r))
  populations = vector("list", length = length(r))
  
  # begin experiment
  for(i in 1:length(r)){
    cat("Functional vs Non Functional","\n",
        "Table size=",r[i],"x",s[i],"\n",
        "Sample size=",n,"\n", sep = "")
    population = func_vs_nfunc_pop(nrow = r[i], ncol = s[i], n = n, ntables = 200, 
                                   noise_levels = noise_levels, margin = 0) # generate population
    populations[[i]] = population
    
    ########################################################################################################
    
    # experiment 1
    cat("Experiment 1","\n",
        "uniform row functional vs transpose", "\n",
        sep = "")
    exp1score = vector("list", length = length(noise_levels))
    exp1stat = vector("list", length = length(noise_levels))
    
    for(nc in 1:length(noise_levels)){
      cat("At noise level=",noise_levels[nc],"\n")
      
      cl = makeForkCluster(nnodes = ncores)
      exp1score[[nc]] = pblapply(cl = cl, X = 1:length(population$exp1fun[[nc]]), FUN = function(j){
        return(rbind(run_experiment(population$exp1fun[[nc]][[j]]),
                     run_experiment(population$exp1nfun[[nc]][[j]])))
      }) # grab results
      stopCluster(cl)
      
      scores = as.data.frame(ListToDataMatrixMethods(exp1score[[nc]]))
      colnames(scores) = methodnames
      ttle = paste0("Functional vs Non-Functional","\n",
                    "UR Func vs NFunc","\n",
                    "S=",n,", T=",r[i],"x",s[i],"\n",
                    "Noise=",noise_levels[nc])
      # plot results
      exp1stat[[nc]] = ggplot.ROC.PR.curves(scores, population$gt1, colr, ltyps, plot=TRUE, ttle)
    }
    
    # store results
    exp1scores[[i]] = exp1score
    exp1stats[[i]] = exp1stat
    
    ########################################################################################################
    
    # experiment 2
    cat("Experiment 2","\n",
        "uniform column functional vs transpose", "\n",
        sep = "")
    exp2score = vector("list", length = length(noise_levels))
    exp2stat = vector("list", length = length(noise_levels))
    
    for(nc in 1:length(noise_levels)){
      cat("At noise level=",noise_levels[nc],"\n")
      
      cl = makeForkCluster(nnodes = ncores)
      exp2score[[nc]] = pblapply(cl = cl, X = 1:length(population$exp2fun[[nc]]), FUN = function(j){
        return(rbind(run_experiment(population$exp2fun[[nc]][[j]]),
                     run_experiment(population$exp2nfun[[nc]][[j]])))
      }) # grab results
      stopCluster(cl)
      
      scores = as.data.frame(ListToDataMatrixMethods(exp2score[[nc]]))
      colnames(scores) = methodnames
      ttle = paste0("Functional vs Non-Functional","\n",
                    "UC Func vs NFunc","\n",
                    "S=",n,", T=",r[i],"x",s[i],"\n",
                    "Noise=",noise_levels[nc])
      # plot results
      exp2stat[[nc]] = ggplot.ROC.PR.curves(scores, population$gt2, colr, ltyps, plot=TRUE, ttle)
    }
    
    # store results
    exp2scores[[i]] = exp2score
    exp2stats[[i]] = exp2stat
    
    ########################################################################################################
  }
  
  # return result containers
  return(list(exp1scores = exp1scores, exp2scores = exp2scores, 
              exp1stats = exp1stats, exp2stats = exp2stats,
              populations = populations))
}

# Summarize functional versus independent AUROC and AUPR scores
# by means of violin plot (if PLOT=TRUE). 
# Requires a list of scores obtained from func_vs_indep()
func_vs_indep_summary = function(scores, PLOT=TRUE){
  
  nsm = length(scores) # samples
  nexp = 3 # marginals
  ntab = length(scores$fvi2000$exp1stats) # table sizes
  nnoi = length(scores$fvi2000$exp1stats[[1]]) # noise levels
  
  ########################################################################################################
  
  # AUROC score distribution
  expauroc = as.data.frame(matrix(0, ncol = length(methodnames), nrow=nsm*nexp*ntab*nnoi))
  count = 1
  
  # sample size 2000
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expauroc[count, ] = scores$fvi2000$exp1stats[[i]][[j]]$AUROC
      expauroc[count+1, ] = scores$fvi2000$exp2stats[[i]][[j]]$AUROC
      expauroc[count+2, ] = scores$fvi2000$exp3stats[[i]][[j]]$AUROC
      count = count + 3
    }  
  }
  
  # sample size 200
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expauroc[count, ] = scores$fvi200$exp1stats[[i]][[j]]$AUROC
      expauroc[count+1, ] = scores$fvi200$exp2stats[[i]][[j]]$AUROC
      expauroc[count+2, ] = scores$fvi200$exp3stats[[i]][[j]]$AUROC
      count = count + 3
    }  
  }
  
  # set column names
  colnames(expauroc) = methodnames
  
  if(PLOT){
    # rank by mean AUROC distribution
    mean_expauroc = unlist(lapply(expauroc,mean), use.names = FALSE)
    expauroc_o = expauroc[,order(mean_expauroc, decreasing = TRUE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(expauroc_o + jitter(rep(0, nrow(expauroc_o))),
                     col=colr[order(mean_expauroc, decreasing = TRUE)], 
                     main=paste0("Functional vs Independent","\n",
                                 "Scores over 180 setups"), 
                     cex.axis=2.5, cex.main=2.5, las=2,
                     plotCentre="line", lwd=3, ylim=c(0,1))
    par(lwd=0.3)
    stripchart(expauroc_o, add=TRUE, vertical = TRUE, jitter = 0.1, pch=23, cex=0.8,
               col="gray", bg="white", method = "jitter")
    abline(h = 0.5, lty=2, lwd=3, col="blue")
    mtext(text="AUROC", side = 2, line=4, cex=2.5)
  }
  
  
  ########################################################################################################
  
  # AUPR score distribution
  expaupr = as.data.frame(matrix(0, ncol = length(methodnames), nrow=nsm*nexp*ntab*nnoi))
  count = 1
  
  # sample size 2000
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expaupr[count, ] = scores$fvi2000$exp1stats[[i]][[j]]$AUPR
      expaupr[count+1, ] = scores$fvi2000$exp2stats[[i]][[j]]$AUPR
      expaupr[count+2, ] = scores$fvi2000$exp3stats[[i]][[j]]$AUPR
      count = count + 3
    }  
  }
  
  # sample size 200
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expaupr[count, ] = scores$fvi200$exp1stats[[i]][[j]]$AUPR
      expaupr[count+1, ] = scores$fvi200$exp2stats[[i]][[j]]$AUPR
      expaupr[count+2, ] = scores$fvi200$exp3stats[[i]][[j]]$AUPR
      count = count + 3
    } 
  }  
  
  # set column names
  colnames(expaupr) = methodnames
  
  if(PLOT){
    # rank by mean AUPR distribution
    mean_expaupr = unlist(lapply(expaupr,mean), use.names = FALSE)
    expaupr_o = expaupr[,order(mean_expaupr, decreasing = TRUE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(expaupr_o + jitter(rep(0, nrow(expaupr_o))),
                     col=colr[order(mean_expaupr, decreasing = TRUE)], 
                     main=paste0("Functional vs Independent","\n",
                                 "Scores over 180 setups"),
                     cex.axis=2.5, cex.main=2.5, las=2,
                     plotCentre="line", lwd=3, ylim=c(0,1))
    par(lwd=0.3)
    stripchart(expaupr_o, add=TRUE, vertical = TRUE, jitter = 0.1, pch=23, cex=0.8,
               col="gray", bg="white", method = "jitter")
    abline(h = 0.5, lty=2, lwd=3, col="blue")
    mtext(text="AUPR", side = 2, line=4, cex=2.5)
  }
  
  ########################################################################################################
  
  return(list(AUROC_C = expauroc, AUPR_C = expaupr))
}

# Summarize functional versus non-functional AUROC and AUPR scores
# by means of violin plot (if PLOT=TRUE).
# Requires a list of scores obtained from func_vs_nfunc()
func_vs_nfunc_summary = function(scores, PLOT=TRUE){
  
  nsm = length(scores) # samples
  nexp = 2 # marginals
  ntab = length(scores$fvnf2000$exp1stats) # table sizes
  nnoi = length(scores$fvnf2000$exp1stats[[1]]) # noise levels
  
  ########################################################################################################
  
  # AUROC score distribution
  expauroc = as.data.frame(matrix(0, ncol = length(methodnames), nrow=nsm*nexp*ntab*nnoi))
  count = 1
  
  # sample size 2000
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expauroc[count, ] = scores$fvnf2000$exp1stats[[i]][[j]]$AUROC
      expauroc[count+1, ] = scores$fvnf2000$exp2stats[[i]][[j]]$AUROC
      count = count + 2
    }
  }  
  
  # sample size 200
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expauroc[count, ] = scores$fvnf200$exp1stats[[i]][[j]]$AUROC
      expauroc[count+1, ] = scores$fvnf200$exp2stats[[i]][[j]]$AUROC
      count = count + 2
    }  
  }
  
  # set column names
  colnames(expauroc) = methodnames
  
  if(PLOT){
    # rank by mean AUROC distribution
    mean_expauroc = unlist(lapply(expauroc,mean), use.names = FALSE)
    expauroc_o = expauroc[,order(mean_expauroc, decreasing = TRUE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(expauroc_o + jitter(rep(0, nrow(expauroc_o))),
                     col=colr[order(mean_expauroc, decreasing = TRUE)], 
                     main=paste0("Functional vs Non-Functional","\n",
                                 "Scores over 120 setups"), 
                     cex.axis=2.5, cex.main=2.5, las=2,
                     plotCentre="line", lwd=3, ylim=c(0,1))
    par(lwd=0.3)
    stripchart(expauroc_o, add=TRUE, vertical = TRUE, jitter = 0.1, pch=23, cex=0.8,
               col="gray", bg="white", method = "jitter")
    abline(h=0.5, lty=2, lwd=3, col="blue")
    mtext(text="AUROC", side = 2, line=4, cex=2.5)
    
  }
  
  ########################################################################################################
  
  # AUPR score distribution
  expaupr = as.data.frame(matrix(0, ncol = length(methodnames), nrow=nsm*nexp*ntab*nnoi))
  count = 1
  
  # sample size 2000
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expaupr[count, ] = scores$fvnf2000$exp1stats[[i]][[j]]$AUPR
      expaupr[count+1, ] = scores$fvnf2000$exp2stats[[i]][[j]]$AUPR
      count = count + 2
    }  
  }
  
  # sample size 200
  for(i in 1:ntab){ # table sizes
    for(j in 1:nnoi){ # noise levels
      expaupr[count, ] = scores$fvnf200$exp1stats[[i]][[j]]$AUPR
      expaupr[count+1, ] = scores$fvnf200$exp2stats[[i]][[j]]$AUPR
      count = count + 2
    }  
  }
  
  # set column names
  colnames(expaupr) = methodnames
  
  if(PLOT){
    # rank by mean AUPR distribution
    mean_expaupr = unlist(lapply(expaupr,mean), use.names = FALSE)
    expaupr_o = expaupr[,order(mean_expaupr, decreasing = TRUE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(expaupr_o + jitter(rep(0, nrow(expaupr_o))),
                     col=colr[order(mean_expaupr, decreasing = TRUE)], 
                     main=paste0("Functional vs Non-Functional","\n",
                                 "Scores over 120 setups"),
                     cex.axis=2.5, cex.main=2.5, las=2,
                     plotCentre="line", lwd=3, ylim=c(0,1))
    par(lwd=0.3)
    stripchart(expaupr_o, add=TRUE, vertical = TRUE, jitter = 0.1, pch=23, cex=0.8,
               col="gray", bg="white", method = "jitter")
    abline(h=0.5, lty=2, lwd=3, col="blue")
    mtext(text="AUPR", side = 2, line=4, cex=2.5)
    
  }
  
  ########################################################################################################
  
  return(list(AUROC_C = expauroc, AUPR_C = expaupr))
}

# Summarize functional versus non-functional directional accuracy scores
# by means of violin plot (if PLOT=TRUE).
# Requires a list of scores obtained from func_vs_nfunc()
func_vs_nfunc_dir_summary = function(scores, PLOT=TRUE){
  
  nsm = length(scores) # samples
  nexp = 2 # marginals
  ntab = length(scores$fvnf2000$exp1stats) # table sizes
  nnoi = length(scores$fvnf2000$exp1stats[[1]]) # noise levels
  
  ########################################################################################################
  
  # Directional accuracy score distribution
  dir_scores = as.data.frame(matrix(0, ncol=length(methodnames), nrow=nsm*nexp*ntab*nnoi))
  colnames(dir_scores) = methodnames
  count = 1
  
  # sample size 2000
  for(i in 1:ntab) { # tables
    for(j in 1:nnoi) { # noise levels
      
      # obtain scores for both direction and return the boolean equivalent of
      # whether the functional values are 'higher' than non-functional values
      
      # experiment 1
      score = lapply(scores$fvnf2000$exp1scores[[i]][[j]],function(x){
        return(x[1,] > x[2,])
      })
      # average over 200 tables in each setup
      dir_scores[count,] =  apply(ListToDataMatrix(score), 2, sum)/200
      
      # experiment 2
      score = lapply(scores$fvnf2000$exp2scores[[i]][[j]],function(x){
        return(x[1,] > x[2,])
      })
      # average over 200 tables in each setup
      dir_scores[count+1,] =  apply(ListToDataMatrix(score), 2, sum)/200
      count = count+2
    }  
  }
  
  # sample size 200
  for(i in 1:ntab) { # tables
    for(j in 1:nnoi) { # noise levels
      
      # obtain scores for both direction and return the boolean equivalent of
      # whether the functional values are 'higher' than non-functional values
      
      # experiment 1
      score = lapply(scores$fvnf200$exp1scores[[i]][[j]],function(x){
        return(x[1,] > x[2,])
      })
      # average over 200 tables in each setup
      dir_scores[count,] =  apply(ListToDataMatrix(score), 2, sum)/200
      
      # experiment 2
      score = lapply(scores$fvnf200$exp2scores[[i]][[j]],function(x){
        return(x[1,] > x[2,])
      })
      # average over 200 tables in each setup
      dir_scores[count+1,] =  apply(ListToDataMatrix(score), 2, sum)/200
      count = count+2
    }
  }
  
  if(PLOT){
    
    # rank by mean directional accuracy distribution
    mean_dir = unlist(lapply(dir_scores, mean), use.names = FALSE)
    dir_scores_o = dir_scores[,order(mean_dir, decreasing = TRUE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(dir_scores_o + jitter(rep(0, nrow(dir_scores_o))),
                     col=colr[order(mean_dir, decreasing = TRUE)], 
                     main=paste0("Functional vs Non-Functional","\n",
                                 "Scores over 120 setups"), 
                     cex.axis=2.5, cex.main=2.6, las=2,
                     plotCentre="line", lwd=3, ylim=c(0,1))
    par(lwd=0.3)
    stripchart(dir_scores_o, add=TRUE, vertical = TRUE, jitter = 0.1, pch=23, cex=0.8,
               col="gray", bg="white", method = "jitter")
    abline(h=0.5, lty=2, lwd=3, col="blue")
    mtext(text="Accuracy", side = 2, line=4, cex=2.5)
  }
  
  ########################################################################################################
  
  return(list(DIR_SCORES = dir_scores))
}

# Summarize overall ranking using all AUROC, AUPR and 
# directional accuracy scores by means of violin plot (if PLOT=TRUE).
# Requires a list of scores obtained from func_vs_nfunc() and func_vs_indep()
Overall_ranking = function(scores, PLOT=TRUE){
  
  # get AUROC / AUPR scores from functional versus independent
  fvi = func_vs_indep_summary(list(fvi2000=scores$fvi2000, fvi200=scores$fvi200), PLOT = FALSE)
  
  # get AUROC / AUPR scores from functional versus non-functional
  fvnf = func_vs_nfunc_summary(list(fvnf2000=scores$fvnf2000, fvnf200=scores$fvnf200), PLOT = FALSE)
  
  # get directional accuracy scores from functional versus non-functional
  fvnf_dir = func_vs_nfunc_dir_summary(list(fvnf2000=scores$fvnf2000, fvnf200=scores$fvnf200), PLOT = FALSE)
  
  # merge scores
  all_scores = rbind(fvnf_dir$DIR_SCORES, fvnf$AUROC_C, fvnf$AUPR_C, fvi$AUROC_C, fvi$AUPR_C)
  
  # rank each score setup (using minimum for resolving ties)
  # rank works with smaller the better, thus we multiply with a negative sign
  for(i in 1:nrow(all_scores)){
    all_scores[i,] = rank(-all_scores[i,], ties.method = "min")
  }
  
  if(PLOT){
  
    # order by mean rank
    mean_scores = unlist(lapply(all_scores, mean), use.names = FALSE)
    all_scores_o = all_scores[,order(mean_scores, decreasing = FALSE)]
    
    # plot
    par(mar=c(7,6,6,6), lwd=5)
    vioplot::vioplot(all_scores_o + jitter(rep(0, nrow(all_scores_o))),
                     col=colr[order(mean_scores, decreasing = FALSE)], 
                     main="Overall ranking over 720 setups", 
                     cex.axis=2.5, cex.main=2.5, las=2,
                     plotCentre="line", ylim=rev(c(1,length(methodnames))), lwd=3)
    mtext(text="Overall Ranking", side = 2, line=4, cex=2.5)
  }
  
  return(list(ALL_SCORES = all_scores))
}

# Main procedure
SimStudy = function(ncores){
  
  # Generates all simulation results / simulation plots.
  
  # Functional versus independent at sample size 2000
  pdf("../Results/SimStudy/Functional_vs_Independent_SampleSize_2000.pdf")
  fvi2000 = func_vs_indep(n=2000, ncores=ncores)
  dev.off()
  save(fvi2000, file="../Results/SimStudy/fvi2000.RData")
  
  # Functional versus independent at sample size 200
  pdf("../Results/SimStudy/Functional_vs_Independent_SampleSize_200.pdf")
  fvi200 = func_vs_indep(n=200, ncores=ncores)
  dev.off()
  save(fvi200, file="../Results/SimStudy/fvi200.RData")
  
  # Functional verus non-functional at sample size 2000
  pdf("../Results/SimStudy/Functional_vs_Nonfunctional_SampleSize_2000.pdf")
  fvnf2000 = func_vs_nfunc(n=2000, ncores=ncores)
  dev.off()
  save(fvnf2000, file="../Results/SimStudy/fvnf2000.RData")
  
  # Functionoal versus non-functional at sample size 200
  pdf("../Results/SimStudy/Functional_vs_Nonfunctional_SampleSize_200.pdf")
  fvnf200 = func_vs_nfunc(n=200, ncores=ncores)
  dev.off()
  save(fvnf200, file="../Results/SimStudy/fvnf200.RData")
  
  load("../Results/SimStudy/fvi2000.RData")
  load("../Results/SimStudy/fvi200.RData")
  load("../Results/SimStudy/fvnf2000.RData")
  load("../Results/SimStudy/fvnf200.RData")
  
  # Summary of functional versus independent
  pdf("../Results/SimStudy/Functional_vs_Independent_AUROC_AUPR_DIST.pdf", width = 10)
  fvisumm = func_vs_indep_summary(list(fvi2000=fvi2000, fvi200=fvi200))
  dev.off()
  
  # Summary of functional versus non-functional
  pdf("../Results/SimStudy/Functional_vs_Nonfunctional_AUROC_AUPR_DIST.pdf", width = 10)
  fvnfsumm = func_vs_nfunc_summary(list(fvnf2000=fvnf2000, fvnf200=fvnf200))
  dev.off()
  pdf("../Results/SimStudy/Functional_vs_Nonfunctional_DIR_ACCU_DIST.pdf", width = 10)
  fvnf_dirsumm = func_vs_nfunc_dir_summary(list(fvnf2000=fvnf2000, fvnf200=fvnf200))
  dev.off()
  
  # Overall ranking in terms of AUROC/AUPR/Directional accuracy distribution
  pdf("../Results/SimStudy/Overall_ranking_sim_study.pdf", width = 10)
  ovsumm = Overall_ranking(list(fvi2000=fvi2000, fvi200=fvi200, fvnf2000=fvnf2000, 
                                fvnf200=fvnf200))
  dev.off()
}