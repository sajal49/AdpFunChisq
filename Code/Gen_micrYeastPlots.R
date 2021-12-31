# Regenerate all plots presented in the main manuscript

# Manuscript:
# "Learning causal networks using inducible transcription factors and transcriptome‐wide time series"
# Hackett.et.al 2020, Molecular Systems Biology.
# Original data obtained from: https://idea.research.calicolabs.com/data

# Ground Truth:
# "YEASTRACT+: a portal for cross-species comparative genomics of transcription regulation in yeasts"
# Interactions obtained from http://www.yeastract.com/formfindregulators.php

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("micrYeast_Eval.R")

# plot directional accuracy over sample size
load("../Results/YeastStudy/yeast_fvnf.Rds")
pdf("../Results/YeastStudy/yeast_dir_test_over_sm.pdf")
yeast_accuracy_lines(yeast_fvnf)
dev.off()

# plot AUROC/AUPR (functional versus independent) over sample size
load("../Results/YeastStudy/yeast_fvi.Rds")
pdf("../Results/YeastStudy/yeast_fvi_test_over_sm.pdf")
yeast_indep_test_lines(yeast_fvi)
dev.off() 