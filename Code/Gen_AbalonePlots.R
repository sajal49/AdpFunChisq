# Regenerate all plots presented in the main manuscript

# Manuscript:
# "A Quantitative Comparison of Dystal and Backpropagation"
# Clark.et.al 1996, Australian Conference on Neural Networks (ACNN'96).
# Original data obtained from: https://archive.ics.uci.edu/ml/datasets/Abalone

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Abalone_Eval.R")

# plot directional accuracy over sample size
load("../Results/AbaloneStudy/abalone_fvnf.Rds")
pdf("../Results/AbaloneStudy/abalone_dir_test_over_sm.pdf")
abalone_accuracy_test_lines(abalone_fvnf$expscores)
dev.off()

# plot AUROC/AUPR (functional versus independent) over sample size
load("../Results/AbaloneStudy/abalone_fvi.Rds")
pdf("../Results/AbaloneStudy/abalone_fvi_test_over_sm.pdf")
abalone_indep_test_lines(abalone_fvi$expstats)
dev.off()