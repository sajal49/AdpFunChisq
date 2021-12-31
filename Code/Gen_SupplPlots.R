# Assuming that results are already available under Results/Supplementary and Results/SimStudy
# Regenerate all plots presented in the supplementary

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SupplementaryPlots.R")

load("../Results/SimStudy/fvi200.RData")
load("../Results/SimStudy/fvi2000.RData")
load("../Results/SimStudy/fvnf200.RData")
load("../Results/SimStudy/fvnf2000.RData")
load("../Results/Supplementary/Abalone_Finer/abalone_fvnf_finer.RData")
load("../Results/Supplementary/MPAL_Finer/mpal_fvnf_finer.RData")
load("../Results/Supplementary/Yeast_Finer/yeast_fvnf_finer.RData")
load("../Results/Supplementary/runtime_res.Rds")

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
SimStudy_FVNF_dir_scores_drate(list(fvnf2000=fvnf2000, fvnf200=fvnf200))
dev.off()

# A snapshot of the simulated data
PlotSimData()

# plot runtime over table and sample size
pdf("../Results/Supplementary/runtime_eval.pdf")
PlotRuntime(timeres[[1]], timeres[[2]])
dev.off()

# Plot mean directional accuracy on Yeast subset using finer discretization
pdf("../Results/Supplementary/Yeast_Finer/yeast_dir_test_finer.pdf")
Yeast_dir_scores_finer(yeast_fvnf_finer)
dev.off()

# Plot mean directional accuracy on MPAL subset using finer discretization
pdf("../Results/Supplementary/MPAL_Finer/mpal_dir_test_finer.pdf")
MPAL_dir_scores_finer(mpal_fvnf_finer)
dev.off()

# Plot top known pattern using finer discretization on MPAL subset
PlotTopKnownPatternsFiner(mpal_fvnf_finer)

# Plot mean directional accuracy on Abalone subset using finer discretization
pdf("../Results/Supplementary/Abalone_Finer/abalone_dir_test_finer.pdf")
abalone_dir_scores_finer(abalone_fvnf_finer$expscores)
dev.off()

