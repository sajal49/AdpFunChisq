# AdaptiveFunChisq

The source code can be used to repeat all pre-processing steps and evaluations presented in the manuscript as well as in supplementary; however, some data files are required to be download and placed in the correct directory. Please check out the README files in each Data/ directory. Some evaluations (such as those on yeast and leukemia datasets) may take multiple days to run.

Alternatively, all pre-processed data and some precomputed result libraries (.RData/.RDS) are available to download at http://www.cs.nmsu.edu/~joemsong/AFC/Supplementary_Code.zip to skip the pre-processing steps and only run evaluations. Additionally, We also provide means to only regenerate plots presented in the manuscript and supplementary without re-running all evaluations.

1. To repeat all experiments and plots for Abalone, in the manuscript: 
    * source Abalone_Eval.R
    * Run AbaloneEval()
2. To only regenerate all plots for Abalone, in the main manuscript:
    * source Gen_AbalonePlots.R
3. To repeat all experiments and plots for MPAL, in the manuscript:
    * Download scADT-All-Hematopoiesis-MPAL-191120.rds and scRNA-All-Hematopoiesis-MPAL-191120.rds from https://github.com/GreenleafLab/MPAL-Single-Cell-2019/ and place them in Data/MPAL/.
    * Download PathwayCommons11.All.hgnc.sif from https://www.pathwaycommons.org/archives/PC2/v11/ and place it in Data/MPAL/.
    * source scMPAL_Eval_applc.R
    * Run MPAL_study()
4. To only regenerate all plots for MPAL, in the main manuscript: 
    * source Gen_scMPALPlots.R
5. To repeat all experiments and plots for perturbed yeast microarray, in the manuscript:
    * Download idea_tall_expression_data.tsv from https://idea.research.calicolabs.com/data and place it in Data/Yeast/.
    * source micrYeast_Eval.R
    * Run YeastEval()
6. To only regenerate all plots for yeast, in the main manuscript:
    * source Gen_micrYeastPlots.R
7. To repeat all experiments and plot for simulation study, in the main manuscript:
    * source SimStudy.R 
    * Run SimStudy()
8. To only regenerate all plots for simulation study, in the main manuscript: 
    * source Gen_SimStudyPlots.R
9. To repeat the table type bias check, in the main manuscript: 
    * source TableTypeBiasCheck.R
    * Run TtyBiasCheck()
10. To regenerate the plot for table type bias check, in the manuscript:
    * source Gen_TableTypeBiasCheck.R
11. To repeat all experiments and plots, in the supplementary:
    * Download scADT-All-Hematopoiesis-MPAL-191120.rds and scRNA-All-Hematopoiesis-MPAL-191120.rds from https://github.com/GreenleafLab/MPAL-Single-Cell-2019/ and place them in Data/MPAL/, if not done already.
    * Download PathwayCommons11.All.hgnc.sif from https://www.pathwaycommons.org/archives/PC2/v11/ and place it in Data/MPAL/, if not done already.
    * Download idea_tall_expression_data.tsv from https://idea.research.calicolabs.com/data and place it in Data/Yeast/, if not done already.
    * source SupplementaryPlots.R
    * Run Suppl_Plots()
12. To only regenerate all plots for supplementary 
    * source Gen_SupplPlots.R
13. Implementations and Ports:
    * Adaptive functional chi-squared test in AdpFunChisq.R 
    * Digital Regression in DR.R
    * Distance Correlation in DC.R
    * Causal Inference by Stochastic Complexity in Cisc.R 
    * Fraction of Information in Methods.R
    * Sub-Copula Regression in Methods.R
    * Goodman-Kruskal Tau in Methods.R
    * Additive Noise Model in ANM.R and GauProRegPy.py
    * Information-Geometric Causal Inference in IGCI.R
