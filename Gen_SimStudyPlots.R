# Assuming that 'SimResults' is already available under Data/Simulation/
# Regenerate all simulation plots presented in the main manuscript

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("SimStudy.R")
load("../Data/Simulation/SimResults.RData")

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