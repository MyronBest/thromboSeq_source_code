# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Summary       : Script to load thromboSeq processed HTSeq, perform data QCs, differential expression 
#                 of splice junction analysis, and PSO-enhanced development of thromboSeq classification
#                 algorithms
# Date          : 1st of September 2018
# Revised       : None
# req. packages : methods       [general functions]
#                 reshape       [visualize data in a heatmap]
#                 affy          [visualize data in a heatmap]
#                 ggplot2       [visualize data in a heatmap]
#                 dendroextras  [visualize data in a heatmap]
#                 e1071         [SVM-algorithm and classification]
#                 ROCR          [generate ROC curves and calculate AUC-values]
#                 caret         [calculate highly correlated transcripts during spliced RNA biomarker panel selection]
#                 pROC          [calculate ROC 95% confidence intervals]
#                 ppso          [run particle-swarm optimization algorithms, 
#                                included as ppsoThromboSeq.tar in the Supplemental Code (/bin/ppsoThromboSeq.tar); 
#                                install via: install.packages("ppsoThromboSeq.tar", type="source", repos=NULL)]
#                 edgeR         [store data in a data object called a DGEList, calculate normalization factor, 
#                                perform ANOVA Likelihood-ratio analysis]
#                 doSNOW/doMC   [parallel computing in order to enable efficient data analysis]
#                 foreach       [enable for repeated executions]
#                 RUVSeq        [perform RUVg data correction]
#
# note: If you do not have one of these packages installed, install by source("https://bioconductor.org/biocLite.R") and then biocLite("packageName")

# open R by running in terminal the following command: R

# continue in R
# load the required functions by running the following commands: 
source('bin/thromboSeqTools_PreProcessing_2.R')
source('bin/thromboSeqTools_ANOVA.R')
source('bin/thromboSeqTools_PSO.R')

# Collect the intron-spanning read counts outputted by the mapping pipeline 
# by running the following command:
counts <- collect.read.counts(inputDir = "symlinksHTSEQ/")

# Prepare a DGE object for subsequent data analysis.
# For this, run the following command:
dge <- prepare.dge.object(sample.info = "sampleInfo.csv",
                          gene.info = "bin/dgeGenesEnsembl75.RData",
                          read.counts = counts)

# Filter the dataset for low abundant RNAs.
# For this, run the following command:
dge <- filter.for.platelet.transcriptome(dge)

# Perform thromboSeqQC pipeline.
# For this, run the following command:
dgeIncludedSamples <- thromboSeqQC(dge)

# Perform thromboSeq ANOVA differential expression analysis of splice junctions
# For this, run the following command:
thromboSeq.anova <- thromboSeqANOVA(dgeIncludedSamples)

# Perform PSO-enhanced thromboSeq classifier development.
# For this, run the following command:
thromboPSO <- thromboSeqPSO(dge = dgeIncludedSamples,
                            n.particles = 100,
                            n.iterations = 10,
                            number.cores = 8)

# Validate the developed thromboSeq algorithm.
# For this, run the following command:
thromboPSOreadout <- thromboSeqPSO.readout(replace.counts.validation = 0,  
                                           number.cores = 8)

# Perform control experiments for thromboSeq classifier development.
# For this, run the following command:
control <- thromboSeqPSO.controls(thromboSeqPSO.shuffled = T,
                                  thromboSeqPSO.iterations = T,
                                  n.shuffled = 100,
                                  n.iterations = 100,
                                  number.cores = 8)

# Optionally save and store the R-session
sessionInfo <- sessionInfo()
save.image(paste(timestamp(),".RData",sep=""))
# End