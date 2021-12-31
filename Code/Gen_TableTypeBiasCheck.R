# Regenerate all plots presented in the main manuscript

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("TableTypeBiasCheck.R")

# plot table type bias between 3x10 versus 10x3 at n = 100 and increasing noise.
load("../Results/TableTypeStudy/res_10x3_3x10_100.RData")
pdf("../Results/TableTypeStudy/tty_10x3_3x10_100.pdf", width = 15)
PlotPatterns(stat_collect = res_10x3_3x10_100, n = 100, r = 10, s = 3)
dev.off()