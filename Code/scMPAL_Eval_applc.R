# Processing, Evaluation and Application on:
# Single-Cell Multiomics of Mutliple-Phenotype Acute Leukemia, for:
# Causal Inference by Functional and Statistical Dependency

# Manuscript:
# "Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia"
# Granja.et.al 2019, Nature Biotechnology.
# Original data obtained from: https://github.com/GreenleafLab/MPAL-Single-Cell-2019/

# Ground Truth:
# "Pathway Commons, a web resource for biological pathway data"
# Cerami.et.al 2010, Nucleic Acids Research.
# Interactions obtained from https://www.pathwaycommons.org/archives/PC2/v11/

# This code requires:
# scADT-All-Hematopoiesis-MPAL-191120.rds and scRNA-All-Hematopoiesis-MPAL-191120.rds
# from https://github.com/GreenleafLab/MPAL-Single-Cell-2019/ to be put in Data/MPAL/ directory.
# PathwayCommons11.All.hgnc.sif
# from https://www.pathwaycommons.org/archives/PC2/v11/ to be put in Data/MPAL/ directory.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
require(DescTools, quietly = TRUE)
require(infotheo, quietly = TRUE)
require(HCR, quietly = TRUE)
require(FunChisq, quietly = TRUE)
require(GridOnClusters, quietly = TRUE)
require(entropy, quietly = TRUE)
require(pbapply, quietly = TRUE)
require(scater, quietly = TRUE)
require(scran, quietly = TRUE)
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
            SCR_M(tbl), FI_M(tbl))
  return(stats)
}

run_cont_experiment = function(tbl){
  
  stats = c(IGCI_M(tbl, "data"), ANM_M(tbl, "data"))
  return(stats)
}

# read data (assumes the .rds files are available under Data/MPAL/)
read_data = function(){
  
  proteins = readRDS("../Data/MPAL/scADT-All-Hematopoiesis-MPAL-191120.rds")
  scData = readRDS("../Data/MPAL/scRNA-All-Hematopoiesis-MPAL-191120.rds")
  
  # match single cells to only keep corresponding entries
  matched_colnm = match(colnames(proteins), colnames(scData))
  return(list(proteins = proteins, scData = scData[,matched_colnm]))
  
}

# pre-process data according to Granja.et.al 2019
preprocess_data = function(PLOT=FALSE){
  
  # read data
  dt = read_data()
  proteins = dt$proteins
  scData = dt$scData
  rm(dt)
  
  # plot UMAP if required
  if(PLOT){
    
    # Visualize raw data
    Visualize_data(proteins, "Immune Surface Proteins", "../Results/MpalStudy/protein_UMAP_projection.pdf")
    Visualize_data(scData, "Single Cell RNAs", "../Results/MpalStudy/scRNAs_UMPA_projection.pdf")  
  
  }
  
  # filter zero genes
  prot_zerog = rowSums(assay(proteins), na.rm = TRUE)
  rna_zerog = rowSums(assay(scData), na.rm = TRUE)
  if(any(prot_zerog == 0)){
    proteins = proteins[prot_zerog != 0,]
  }
  if(any(rna_zerog == 0)){
    scData = scData[rna_zerog != 0,]
  }
  
  # filter zero cells
  prot_zeroc = colSums(assay(proteins), na.rm=TRUE)
  rna_zeroc = colSums(assay(scData), na.rm = TRUE)
  if(any(prot_zeroc == 0 | rna_zeroc == 0)){
    rem_ind = which(prot_zeroc == 0 | rna_zeroc == 0)
    proteins = proteins[,-rem_ind]
    scData = scData[,-rem_ind]
  }
  
  # Normalize protein
  norm_proteins = Normalize_proteins(assay(proteins), ncores=6)
  
  # Normalize scRNA
  norm_scdata = Normalize_scRNA(assay(scData), ncores=6)

  # prepare meta-info
  gene_names = rownames(scData)
  protein_names = rownames(proteins)
  sample_names = colnames(scData)
  bioclassification = proteins$ProjectClassification
  
  # convert protein CD names to gene names
  alt_protein_names = c("MME", "IL3RA", "CD14", "CD19", "CD3D", "CD33", "CD34", "CD38", "CD4", 
                        "PTPRC-A", "CD7", "CD8A", "THY1", "FUT4", "FCGR3A", "NCAM1", "IL2RA", "PTPRC-O",
                        "PDCD1", "TIGIT", "IL7R")
  
  # Get known interactions from Pathway Commons
  interactions = Get_known_protein_interactions(alt_protein_names)
  
  # Get gene variance by means of number of non-zero single cells * (MAD of non-zero single cells)
  gene_var = Get_gene_var(norm_scdata, 6)
  
  # Save pre-processed data under Data/MPAL/
  rm(list = setdiff(ls(), c("norm_proteins", "norm_scdata","gene_names","protein_names",
                            "sample_names", "bioclassification", "alt_protein_names",
                            "interactions", "gene_var")))
  save.image("../Data/MPAL/normalized_filtered_data.RData")
}

# Calculate gene variance
Get_gene_var = function(norm_scdata, ncores){
  
  cl = makeForkCluster(nnodes = ncores)
  
  # gene variance = number of non-zero single cells * (MAD of non-zero single cells)
  gene_var = parLapply(cl, c(1:nrow(norm_scdata)), function(i){
    x = norm_scdata[i,]
    x = mad(x[x!=0])*sum(x!=0)
  })
  
  gene_var = unlist(gene_var, use.names = FALSE)
  stopCluster(cl)
  return(gene_var)
  
}

# get known protein interactions from Pathwaycommons.org
Get_known_protein_interactions = function(alt_protein_names){
  
  # read interactions
  mod_prot_names = unique(gsub("-.*","",alt_protein_names))
  all_interactions = read.table("../Data/MPAL/PathwayCommons11.All.hgnc.sif", sep="\t")
  
  # find interactions with alt_protein_names
  int_ind1 = match(all_interactions[,1], mod_prot_names)
  int_ind2 = match(all_interactions[,3], mod_prot_names)
  keep_ind = which(!is.na(int_ind1) | !is.na(int_ind2))
  
  relevant_interactions = all_interactions[keep_ind,]
  
  # causal interactions
  causal_ind = grep("controls", relevant_interactions[,2])
  causal_interactions = relevant_interactions[causal_ind, ]
  causal_interactions = causal_interactions[,c(1,3)]
  causal_interactions = causal_interactions[!duplicated(causal_interactions),]
  
  return(list(relevant_interactions = relevant_interactions, causal_interactions = causal_interactions))
}

# Normalize single-cell RNA data
Normalize_scRNA = function(scdata_exp, ncores){
  
  cl = makeForkCluster(nnodes = ncores)
  
  # Transform using CPT
  norm_scdata = parLapply(cl = cl, c(1:ncol(scdata_exp)), function(i){
    x = scdata_exp[,i]
    x = log2(((x * 10000) / sum(x))+1)
    return(x)
  })
  stopCluster(cl)
  
  norm_scdata = t(ListToDataMatrix(norm_scdata))
  return(norm_scdata)
}

# Normalize single-cell protein data
Normalize_proteins = function(prot_exp, ncores){
  
  cl = makeForkCluster(nnodes = ncores)
  
  # Transform using centered log ratio
  norm_prot = parLapply(cl = cl, c(1:ncol(prot_exp)), function(i){
    x = prot_exp[,i]
    logx = log(x)
    geo_mean_x = exp(mean(logx[!is.infinite(logx)], na.rm=TRUE))
    clr_x = log(x / geo_mean_x)
    return(clr_x)
  })
  stopCluster(cl)
  
  norm_prot = t(ListToDataMatrix(norm_prot))
  return(norm_prot)
}

# Visualize UMAP of a singlecellexperiment object
Visualize_data = function(scexp, ttle, pdfpath){
  
  bioclass = as.factor(scexp$ProjectClassification)
  
  colr = rainbow(length(levels(bioclass)))
  colrs = colr[bioclass]
  
  pdf(pdfpath)
  par(mar=c(5,5,5,9))
  plot(scexp$ProjectUMAP2~proteins$ProjectUMAP1, col=colrs, pch=19, xlab="UMAP Dimension 1",
       ylab = "UMAP Dimension 2", main=ttle, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
  par(xpd=TRUE)
  legend("topright", legend = levels(bioclass), fill = colr, inset=c(-0.43,0), bty = "n", cex=1.1)
  par(xpd=FALSE)
  dev.off()
}

# Evaluate known interaction patterns -- causal versus non-causal
MPAL_Known_Interaction_FVNF_Evaluation = function(){
  
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
  
  # Compute all scores
  results = pblapply(c(1:length(parent_index)), function(i){
    
    # Only consider non-zero single cell in both protein and RNA
    p = as.numeric(norm_proteins[parent_index[i],])
    c = as.numeric(norm_scdata[child_index[i],])
    z = is.na(p) | is.infinite(p) | c==0
    p = p[!z]
    c = c[!z]
    tblx = cbind(p, c)
    tbly = cbind(c, p)
    
    # if <10 elements then return the worst possible score for all methods    
    if(length(p) < 10){
      return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                   rep(-.Machine$integer.max, length(methodnames))))
    } else {
      
      # Discretize jointly using GridOnCluster
      d = discretize.jointly(data = cbind(p,c),
                             k=Inf,
                             cluster_method = 'kmeans+silhouette',
                             min_level = 2,
                             grid_method = 'Sort+split')  
      
      # initiate garbage collection
      gc()
      
      # prepare contingency table
      p = d$D[,1]
      c = d$D[,2]
      tbl = table(p, c)
      
      # if the table has only 1 row or column, then return the worst possible score for all methods  
      if(ncol(tbl) == 1 || nrow(tbl) == 1){
        return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                     rep(-.Machine$integer.max, length(methodnames))))
      } else { # run experiment
        return(rbind(c(run_experiment(tbl), run_cont_experiment(tblx)), 
                     c(run_experiment(t(tbl)), run_cont_experiment(tbly))))
      }
    }
    
  })
  
  # Convert result to data.frame
  results = ListToDataMatrixMethods(results)
  results = as.data.frame(results)
  colnames(results) = methodnames
  
  # ground truth
  gt = rep(c(1,0), each=(nrow(results)/2))
  
  # plot results
  ttle = paste0("Single Cell Leukemia Data","\n",
                "Causal Interactions vs Reverse")
  stats = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
  return(list(stats = stats, scores = results))  
  
}

# Evaluate known interaction patterns -- causal versus independent
MPAL_Known_Interaction_FVI_Evaluation = function(){
  
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
  
  # Compute all scores
  results = pblapply(c(1:length(parent_index)), function(i){
    
    # Only consider non-zero single cell in both protein and RNA
    p = as.numeric(norm_proteins[parent_index[i],])
    c = as.numeric(norm_scdata[child_index[i],])
    z = is.na(p) | is.infinite(p) | c==0
    p = p[!z]
    c = c[!z]
    tblx = cbind(p, c)
    
    p_i = sample(p, length(p))
    c_i = sample(c, length(c))
    tblxi = cbind(p_i, c_i)
    
    # if <10 elements then return the worst possible score for all methods    
    if(length(p) < 10){
      return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                   rep(-.Machine$integer.max, length(methodnames))))
    } else {
      
      # Discretize jointly using GridOnCluster
      d = discretize.jointly(data = cbind(p,c),
                             k=Inf,
                             cluster_method = 'kmeans+silhouette',
                             min_level = 2,
                             grid_method = 'Sort+split')  
      
      # initiate garbage collection
      gc()
      
      # prepare contingency tables
      
      # causal interaction
      p = d$D[,1]
      c = d$D[,2]
      
      # independent pattern
      p_i = sample(p, length(p))
      c_i = sample(c, length(c))
      
      # create tables
      tbl = table(p, c)
      tbl_i = table(p_i, c_i)
      
      # if the table has only 1 row or column, then return the worst possible score for all methods  
      if(ncol(tbl) == 1 || nrow(tbl) == 1){
        return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                     rep(-.Machine$integer.max, length(methodnames))))
      } else { # run experiment
        return(rbind(c(run_experiment(tbl), run_cont_experiment(tblx)), 
                     c(run_experiment(tbl_i), run_cont_experiment(tblxi))))
      }
    }
    
  })
  
  # Convert result to data.frame
  results = ListToDataMatrixMethods(results)
  results = as.data.frame(results)
  colnames(results) = methodnames
  
  # ground truth
  gt = rep(c(1,0), each=(nrow(results)/2))
  
  # plot results
  ttle = paste0("Single Cell Leukemia Data","\n",
                "Causal Interactions vs Shuffled Interactions")
  stats = ggplot.ROC.PR.curves(results, gt, colr, ltyps, plot=TRUE, ttle)
  return(list(stats = stats, scores = results))  
  
}

# Plot directional accuracy scores
# Requires scores from MPAL_Known_Interaction_FVNF_Evaluation()
Dir_accu_scores = function(expscores){
  
  # number of experiments
  scores = expscores$scores
  ntest = nrow(scores)/2
  
  # get the right direction
  dir_test = unlist(lapply(scores, function(i){
    return(sum(i[1:ntest] > i[(ntest+1):(ntest*2)]))
  }), use.names = FALSE)
  
  # title of the plot
  ttle = paste0("Directional Accuracy","\n",
                "Protein->RNA versus","\n", 
                "RNA->Protein")
  
  # order by directional accuracy
  ord = order(dir_test, decreasing = TRUE)
  meths = ordered(methodnames[ord], levels=methodnames[ord])
  dir_scores = data.frame(accu=((dir_test[ord]*100)/ntest), 
                          meths = meths)
  
  # plot accuracy
  plt = ggplot(data=dir_scores, aes(x=meths, y=accu)) +
    geom_bar(stat="identity", fill=colr[ord], color="black")+
    ylim(c(0,110)) +
    geom_text(aes(label=paste0(round(accu, digits = 2),"%")), vjust=0.3, 
              hjust=-0.1, color="black", size=6.5, angle=90)+
    theme_bw(base_size = 22) + xlab("Methods") + ylab("Accuracy") +
    ggtitle(ttle) +
    theme(title = element_text(face = "bold", size=21),
          legend.position = "none", 
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1))
  print(plt)
  
}

# Plot causal versus independent pattern
# Requires scores from MPAL_Known_Interaction_FVI_Evaluation()
Causal_vs_independent_scores = function(expscores){
  
  # number of experiments
  scores = expscores$scores
  ntest = nrow(scores)/2
  
  # get stats for plotting
  gt = rep(c(1,0),each=ntest)
  stats = ggplot.ROC.PR.curves(list.stats = scores, true.classes=gt, 
                               colr = colr, ltyps = ltyps, plot=FALSE, title = NULL)
  auc.ROCs = stats$AUROC
  auc.PRs = stats$AUPR
  res = stats$res
  
  # title for ROC
  ttle1 = paste0("Receiver Operator Characteristic","\n",
                 "Protein->RNA versus","\n",
                 "Shuffled Protein->RNA")
  
  ROC <- NULL
  
  for(i in seq(ncol(scores))) {
    d <- data.frame(
      FPR = res[[i]]$FPR,
      TPR = res[[i]]$TPR,
      Method=paste0(names(scores)[i], ' (',
                    format(auc.ROCs[i], digits=2), ')')
    )
    ROC <- rbind(ROC, d)
  }
  
  # plot ROC
  p <- ggplot(ROC, aes(x=FPR, y=TPR, group=Method, 
                       color=Method, linetype=Method)) +
    scale_color_manual(values=colr[order(methodnames)]) +
    scale_linetype_manual(values=ltyps[order(methodnames)])+
    labs(x="False positive rate", y="True positive rate", 
         title = ttle1, color="Method (AUROC)", 
         linetype="Method (AUROC)") + 
    geom_line(size=1) + xlim(0, 1) + ylim(0, 1) +
    theme_bw(base_size = 22) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          title = element_text(face = "bold", size=21),
          legend.text = element_text(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  print(p)
  
  # title for PR
  ttle2 = paste0("Precision Recall","\n",
                 "Protein->RNA versus","\n",
                 "Shuffled Protein->RNA")
  PRC <- NULL
  for(i in seq(ncol(scores))) {
    
    d <- data.frame(
      Recall=res[[i]]$Recall,
      Precision=res[[i]]$Precision,
      Method=paste0(names(scores)[i], ' (',
                    format(auc.PRs[i], digits=2), ')')
    )
    
    PRC <- rbind(PRC, d)
    
  }
  
  # plot PR
  p <- ggplot(PRC, aes(x=Recall, y=Precision, 
                       group=Method, color=Method, 
                       linetype=Method)) +
    labs(title = ttle2, color="Method (AUPR)", 
         linetype="Method (AUPR)") +
    scale_color_manual(values=colr[order(methodnames)]) +
    geom_line(size=1) + xlim(0, 1) + ylim(0, 1) +
    scale_linetype_manual(values=ltyps[order(methodnames)])+
    theme_bw(base_size = 22) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          title = element_text(face = "bold", size=21),
          legend.text = element_text(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  print(p)
}

# Plot top known patterns
# Requires scores from MPAL_Known_Interaction_FVNF_Evaluation()
PlotTopKnownPatterns = function(scores){
  
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
  
  # number of experiments
  ntest = nrow(scores)/2
  
  # Plot top known causal patterns for each method
  for(i in 1:length(methodnames)){
    
    # get method scores
    meth_scores = scores[,i]
    
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
      
      png(paste0("../Results/MpalStudy/TopPatterns/Known/",methodnames[i],"_",j,".png"))
      
      # get effective index
      ind = meth_score_ord[j]
      dir = direction[ind]
      
      # plot pattern with a cross
      if(dir == -1){
        
        p = as.numeric(norm_proteins[parent_index[ind],])
        c = as.numeric(norm_scdata[child_index[ind],])
        
        z = is.na(p) | is.infinite(p) | c==0
        p = p[!z]
        c = c[!z]
        
        clr = RColorBrewer::brewer.pal("RdYlGn", n=6)[6]
        
        bins = ceiling(length(p)/1000) # 1000 points in 1 bin
        
        data = as.data.frame(cbind(c,p))
        plt = 
          ggplot(data, aes(x=c, y=p) ) +
          geom_point(color=clr, pch=19, lwd=3) +
          stat_density_2d(aes(fill = (..level..)), geom = "polygon", contour = TRUE, bins=bins) +
          coord_cartesian(ylim=c(min(p),max(p))) +
          geom_text(x=min(c), y=max(p), label=expression(symbol("\264")), hjust=-0.1, vjust=1.1,
                    size=50, col="red") +
          geom_smooth(linetype="solid", color="firebrick4", se=FALSE, 
                      method="gam", lwd=2, formula = y ~ s(x, bs = "cs"))+
          scale_fill_distiller(palette = "RdYlGn", direction = -1) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_bw(base_size = 24) + ylab(paste0("PROT:",alt_protein_names[parent_index[ind]])) +
          xlab(paste0("RNA:",gene_names[child_index[ind]])) +
          theme(panel.background = element_rect(fill = 'white'),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.title = element_text(face = "bold", size = 40))#,
          print(plt)
        
      } else { # plot pattern with a check
        
        p = as.numeric(norm_proteins[parent_index[ind],])
        c = as.numeric(norm_scdata[child_index[ind],])
        
        z = is.na(p) | is.infinite(p) | c==0
        p = p[!z]
        c = c[!z]
        
        clr = RColorBrewer::brewer.pal("RdYlGn", n=6)[6]
        
        bins = ceiling(length(p)/1000) # 1000 points in 1 bin
        
        data = as.data.frame(cbind(p,c))
        plt = 
          ggplot(data, aes(x=p, y=c) ) +
          geom_point(color=clr, pch=19, lwd=3) +
          stat_density_2d(aes(fill = (..level..)), geom = "polygon", contour = TRUE, bins=bins) +
          coord_cartesian(ylim=c(min(c),max(c))) +
          geom_text(x=min(p), y=max(c), label=expression(symbol("\326")), hjust=-0.1, vjust=1.1,
                    size=50, col="dodgerblue3") +
          geom_smooth(linetype="solid", color="firebrick4", se=FALSE, 
                      method="gam", lwd=2, formula = y ~ s(x, bs = "cs"))+
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


# Apply all methods to a subset including 4 surface proteins and genes in viral carcinogen
MPAL_ViralCarcinogen_Interactions_FVNF = function(){
  
  load("../Data/MPAL/sub_data_for_viral_carcinogens_4_surface_protein.RData")
  
  # each protein gets paired with RNAGenes
  parent_index = rep(c(1:4), each=nrow(rnaGenes))
  child_index = rep(c(1:nrow(rnaGenes)), 4)
  
  # Compute all scores
  results = pblapply(c(1:length(parent_index)), function(i){
    
    # Only consider non-zero single cell in both protein and RNA
    p = as.numeric(surface_proteins[parent_index[i],])
    c = as.numeric(rnaGenes[child_index[i],])
    z = is.na(p) | is.infinite(p) | c==0
    p = p[!z]
    c = c[!z]
    
    if(length(p) > 20000){
      indx = sample(length(p), 20000)
      p = p[indx]
      c = c[indx]
    }
    
    tblx = cbind(p, c)
    tbly = cbind(c, p)
    
    # if <10 elements then return the worst possible score for all methods  
    if(length(p) < 10){
      return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                   rep(-.Machine$integer.max, length(methodnames))))
    } else {
      
      # Discretize jointly using GridOnCluster
      tbl = discretize.jointly(data = cbind(p,c),
                               k=Inf,
                               cluster_method = 'kmeans+silhouette',
                               min_level = 2,
                               grid_method = 'Sort+split')  
      
      # initiate garbage collection
      gc()
      
      # prepare contingency table
      p = tbl$D[,1]
      c = tbl$D[,2]
      tbl = table(p, c)
      
      # if the table has only 1 row or column, then return the worst possible score for all methods  
      if(ncol(tbl) == 1 || nrow(tbl) == 1){
        return(rbind(rep(-.Machine$integer.max, length(methodnames)),
                     rep(-.Machine$integer.max, length(methodnames))))
      } else { # run experiment
        return(rbind(c(run_experiment(tbl), run_cont_experiment(tblx)), 
                     c(run_experiment(t(tbl)), run_cont_experiment(tbly))))
      }
    }
    
  })
  
  # Convert result to data.frame
  results = ListToDataMatrixMethods(results)
  results = as.data.frame(results)
  colnames(results) = methodnames
  
  return(results)
  
}

# Plot top patterns for pathway -- Viral Carcinogenesis
# Requires scores from MPAL_ViralCarcinogen_Interactions_FVNF()
PlotTopPathwayPatterns = function(scores){
  
  load("../Data/MPAL/sub_data_for_viral_carcinogens_4_surface_protein.RData")
  
  # each protein gets paired with 1000 RNAs
  parent_index = rep(c(1:4), each=nrow(rnaGenes))
  child_index = rep(c(1:nrow(rnaGenes)), 4)
  
  # number of experiments
  ntest = nrow(scores)/2
  
  # Plot top putative causal patterns for each method
  for(i in 1:length(methodnames)){
    
    # get method scores
    meth_scores = scores[,i]
    
    # pick one direction
    eff_scores = sapply(c(1:ntest), function(k){
      return(ifelse(meth_scores[k] > meth_scores[k+ntest], meth_scores[k], meth_scores[k+ntest]))
    })
    
    # is the picked direction protein to RNA (1) or RNA to protein (-1)
    direction = ifelse(meth_scores[1:ntest] > meth_scores[(ntest+1):(2*ntest)], 1, -1)
    
    # rank the picked direction (Note: methods in Methods.R return scores that are higher the better)
    meth_score_ord = order(eff_scores, decreasing = TRUE)
    
    # plot top 5
    for(j in 1:5){
      
      png(paste0("../Results/MpalStudy/TopPatterns/Discovery/",methodnames[i],"_",j,".png"))
      
      # effective index
      ind = meth_score_ord[j]
      dir = direction[ind]
      
      if(dir == -1){
        
        p = as.numeric(surface_proteins[parent_index[ind],])
        c = as.numeric(rnaGenes[child_index[ind],])
        
        z = is.na(p) | is.infinite(p) | c==0
        p = p[!z]
        c = c[!z]
        
        if(length(p) > 20000){
          indx = sample(length(p), 20000)
          p = p[indx]
          c = c[indx]
        }
        
        bins = round(length(p)/1000) # 1000 points in 1 bin
        clr = RColorBrewer::brewer.pal("Spectral", n=6)[6]
        
        data = as.data.frame(cbind(c,p))
        plt = 
          ggplot(data, aes(x=c, y=p) ) +
          geom_point(color=clr) +
          stat_density_2d(aes(fill = (..level..)), geom = "polygon", contour = TRUE, bins=bins) +
          coord_cartesian(ylim=c(min(p),max(p))) +
          geom_smooth(linetype="solid", color="firebrick4", se=FALSE, 
                      method="gam", lwd=2, formula = y ~ s(x, bs = "cs"))+
          scale_fill_distiller(palette = "Spectral", direction = -1) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_bw(base_size = 24) + ylab(paste0("PROT:",surface_markers[parent_index[ind]])) +
          xlab(paste0("RNA:",pthwGenes[child_index[ind]])) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.title = element_text(face = "bold", size = 40))
        print(plt)
        
      } else {
        
        p = as.numeric(surface_proteins[parent_index[ind],])
        c = as.numeric(rnaGenes[child_index[ind],])
        
        z = is.na(p) | is.infinite(p) | c==0
        p = p[!z]
        c = c[!z]
        
        if(length(p) > 20000){
          indx = sample(length(p), 20000)
          p = p[indx]
          c = c[indx]
        }
        
        clr = RColorBrewer::brewer.pal("Spectral", n=6)[6]
        
        bins = ceiling(length(p)/1000) # 1000 points in 1 bin
        
        data = as.data.frame(cbind(p,c))
        plt = 
          ggplot(data, aes(x=p, y=c) ) +
          geom_point(color=clr) +
          stat_density_2d(aes(fill = (..level..)), geom = "polygon", contour = TRUE, bins=bins) +
          coord_cartesian(ylim=c(min(c),max(c))) +
          geom_smooth(linetype="solid", color="firebrick4", se=FALSE, 
                      method="gam", lwd=2, formula = y ~ s(x, bs = "cs"))+
          scale_fill_distiller(palette = "Spectral", direction = -1) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_bw(base_size = 24) + xlab(paste0("PROT:",surface_markers[parent_index[ind]])) +
          ylab(paste0("RNA:",pthwGenes[child_index[ind]])) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                legend.position = "none", axis.title = element_text(face = "bold", size = 40))#,
        print(plt)
        
      }
      dev.off()
    }
    
  }
}

DiscretizationLvlCheck = function(){
  
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
  
  # Compute all scores
  results = pblapply(c(1:length(parent_index)), function(i){
    
    # Only consider non-zero single cell in both protein and RNA
    p = as.numeric(norm_proteins[parent_index[i],])
    c = as.numeric(norm_scdata[child_index[i],])
    z = is.na(p) | is.infinite(p) | c==0
    p = p[!z]
    c = c[!z]
    
    # if <10 elements then return the worst possible score for all methods    
    if(length(p) < 10){
      return(c(0,0))
    } else {
      
      # Discretize jointly using GridOnCluster
      d = discretize.jointly(data = cbind(p,c),
                             k=Inf,
                             cluster_method = 'kmeans+silhouette',
                             min_level = 2,
                             grid_method = 'Sort+split')  
      
      # initiate garbage collection
      gc()
      
      return(as.numeric(apply(d$D, 2, max)))
    }
    
  })
  
  dis_results = unlist(results, use.names = FALSE)
  
  cat("Minimum Discretization Levels: ", min(dis_results[dis_results!=0]), "\n")
  cat("Mean Discretization Levels: ",mean(dis_results[dis_results!=0]),"\n")
  cat("Maximum Discretization Levels: ", max(dis_results[dis_results!=0]), "\n")
  
  # table type stats
  ttype = rep("", length(results))
  for(i in 1:length(results)){
    if(results[[i]][1] == results[[i]][2]){
      ttype[i] = "square"
    } else if(results[[i]][1] > results[[i]][2]){
      ttype[i] = "portrait"
    } else {
      ttype[i] = "landscape"
    }
  }
  
  MPAL_sc_table_types = ttype
  tab_ttype = table(MPAL_sc_table_types)
  
  cat("\nTable type summary counts: \n")
  print(tab_ttype)
  cat("\nTable type summary proportions: \n")
  print(round((tab_ttype/sum(tab_ttype)), digits = 2))
  
  save(results,file="../Data/MPAL/Table_Stats.RData")
  
}

DiscretizationLvlCheck_ViralCarcinogen = function(){
  
  load("../Data/MPAL/sub_data_for_viral_carcinogens_4_surface_protein.RData")
  
  # each protein gets paired with RNAGenes
  parent_index = rep(c(1:4), each=nrow(rnaGenes))
  child_index = rep(c(1:nrow(rnaGenes)), 4)
  
  # Compute all scores
  results = pblapply(c(1:length(parent_index)), function(i){
    
    # Only consider non-zero single cell in both protein and RNA
    p = as.numeric(surface_proteins[parent_index[i],])
    c = as.numeric(rnaGenes[child_index[i],])
    z = is.na(p) | is.infinite(p) | c==0
    p = p[!z]
    c = c[!z]
    
    if(length(p) > 20000){
      indx = sample(length(p), 20000)
      p = p[indx]
      c = c[indx]
    }
    
    # if <10 elements then return the worst possible score for all methods  
    if(length(p) < 10){
      return(c(0,0))
    } else {
      
      # Discretize jointly using GridOnCluster
      d = discretize.jointly(data = cbind(p,c),
                             k=Inf,
                             cluster_method = 'kmeans+silhouette',
                             min_level = 2,
                             grid_method = 'Sort+split')  
      
      # initiate garbage collection
      gc()
      
      return(as.numeric(apply(d$D, 2, max)))
    }
    
  })
  
  dis_results = unlist(results, use.names = FALSE)
  
  cat("Minimum Discretization Levels: ", min(dis_results[dis_results!=0]), "\n")
  cat("Mean Discretization Levels: ",mean(dis_results[dis_results!=0]),"\n")
  cat("Maximum Discretization Levels: ", max(dis_results[dis_results!=0]), "\n")
  
  # table type stats
  ttype = rep("", length(results))
  for(i in 1:length(results)){
    if(results[[i]][1] == results[[i]][2]){
      ttype[i] = "square"
    } else if(results[[i]][1] > results[[i]][2]){
      ttype[i] = "portrait"
    } else {
      ttype[i] = "landscape"
    }
  }
  
  MPAL_sc_table_types = ttype
  tab_ttype = table(MPAL_sc_table_types)
  
  cat("\nTable type summary counts: \n")
  print(tab_ttype)
  cat("\nTable type summary proportions: \n")
  print(round((tab_ttype/sum(tab_ttype)), digits = 2))
  
  save(results,file="../Data/MPAL/Table_Stats_ViralCarcinogen.RData")
  
}


# Main procedure
# Call MPAL_study() for evaluations only
# Call MPAL_study(ttype_check=TRUE) for evaluations along with table type summary.
# Call MPAL_study(prep_data=TRUE, ttype_check=TRUE) for pre-processing and all evaluations.
MPAL_study = function(prep_data = FALSE, ttype_check=FALSE){
  
  # Pre-process data if desired. 
  if(prep_data){
    preprocess_data(PLOT=TRUE)  
  }
  
  ########################################################################################################
  
  # evaluate all methods on known causal interactions versus their reverse
  pdf("../Results/MpalStudy/MPAL_known_causal_versus_reverse.pdf")
  known_eval_fvnf_scores = MPAL_Known_Interaction_FVNF_Evaluation()
  dev.off()
  save(known_eval_fvnf_scores, file="../Data/MPAL/known_eval_fvnf_scores.RData")

  # evaluate all methods on known causal interactions versus independent
  pdf("../Results/MpalStudy/MPAL_known_causal_versus_shuffled.pdf")
  known_eval_fvi_scores = MPAL_Known_Interaction_FVI_Evaluation()
  dev.off()
  save(known_eval_fvi_scores, file="../Data/MPAL/known_eval_fvi_scores.RData")

  # Plot directional accuracy scores
  pdf("../Results/MpalStudy/MPAL_directional_accuracy.pdf")
  Dir_accu_scores(known_eval_fvnf_scores)
  dev.off()

  # Plot ROC / PR (formatted for manuscript) for causal versus independent
  pdf("../Results/MpalStudy/MPAL_ROC_PR_causal_versus_independent.pdf")
  Causal_vs_independent_scores(known_eval_fvi_scores)
  dev.off()
  
  # Plot top 3 patterns for known
  load("../Data/MPAL/known_eval_fvnf_scores.RData")
  PlotTopKnownPatterns(known_eval_fvnf_scores$scores)
  
  if(ttype_check){
    # Table size and type stats
    DiscretizationLvlCheck()  
  }
  
  ########################################################################################################
  
  # find pairs in Viral carcinogens pathway
  pathway_fvnf_scores = MPAL_ViralCarcinogen_Interactions_FVNF()
  save(pathway_fvnf_scores, file="../Data/MPAL/pathway_fvnf_scores.RData")
  
  # Plot top 5 patterns for Viral carcinogens pathway
  load("../Data/MPAL/pathway_fvnf_scores.RData")
  PlotTopPathwayPatterns(pathway_fvnf_scores)
  
  if(ttype_check){
    # Table size and type stats
    DiscretizationLvlCheck_ViralCarcinogen()  
  }
  
  ########################################################################################################
}
