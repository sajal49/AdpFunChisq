# Assuming that all data subsets / precomputed scores are available under Data/MPAL/
# Regenerate all plots presented in the main manuscript

# Manuscript:
# "Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia"
# Granja.et.al 2019, Nature Biotechnology.
# Original data obtained from: https://github.com/GreenleafLab/MPAL-Single-Cell-2019/

# Ground Truth:
# "Pathway Commons, a web resource for biological pathway data"
# Cerami.et.al 2010, Nucleic Acids Research.
# Interactions obtained from https://www.pathwaycommons.org/archives/PC2/v11/

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("scMPAL_Eval_applc.R")

# Plot directional accuracy scores
load("../Results/MpalStudy/known_eval_fvnf_scores.RData")
pdf("../Results/MpalStudy/MPAL_directional_accuracy.pdf")
Dir_accu_scores(known_eval_fvnf_scores)
dev.off()

# Plot ROC / PR (formatted for manuscript) for causal versus independent
load("../Results/MpalStudy/known_eval_fvi_scores.RData")
pdf("../Results/MpalStudy/MPAL_ROC_PR_causal_versus_independent.pdf")
Causal_vs_independent_scores(known_eval_fvi_scores)
dev.off()

# Plot top 3 patterns for known
PlotTopKnownPatterns(known_eval_fvnf_scores$scores)

# Plot top 5 patterns for pathway
load("../Results/MpalStudy/pathway_fvnf_scores.RData")
PlotTopViralCarcinogenPatterns(pathway_fvnf_scores)
