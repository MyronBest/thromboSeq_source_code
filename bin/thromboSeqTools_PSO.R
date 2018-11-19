# R-tools for thromboSeq
# Functions included for PSO-enhanced classification algorithm development
# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Date          : 1st of September 2018
# Revision      : 19th of November 2018

thromboSeqPSO <- function(dge = dgeIncludedSamples,
                          percentage.for.training = 40,
                          percentage.for.evaluation = 30,
                          training.samples.provided = NA,
                          evaluation.samples.provided = NA,
                          swarm.parameters = c("lib.size","fdr","correlatedTranscripts","rankedTranscripts"),
                          swarm.boundaries = c(-0.1, 1.0, 0.1, 1.0, 0.5, 1.0, 200, "all"),
                          k.variables = 3,
                          variable.to.assess = c("Age","lib.size"),
                          variable.threshold = c(0.2, 0.8),
                          ruvg.pvalue.threshold.group = 1e-2,
                          ruvg.pvalue.threshold.strongest.variable = 1e-2,
                          select.biomarker.FDR = FALSE,
                          minimum.n.transcripts.biomarkerpanel = 2,
                          svm.gamma.range = 2^(-20:0),
                          svm.cost.range = 2^(0:20),
                          number.cross.splits = 2,
                          n.particles.gamma.cost.optimization = 50,
                          n.iterations.gamma.cost.optimization = 4,
                          n.particles = 100,
                          n.iterations = 10,
                          figureDir = "figureOutputFolder",
                          number.cores = 2, 
                          verbose = TRUE){
  # Perform PSO-enhanced thromboSeq classification algorithm development. For this, first
  # select the training and evaluation series. Next, perform PSO-optimization. Finally, summarize data
  # and output the PSO progression plot.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   percentage.for.training: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the training series.
  #   percentage.for.evaluation: Numeric value indicating the percentage of samples per group to be
  #                            assigned to the evaluation series.
  #   training.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           training series.
  #   evaluation.samples.provided: Vector with specified column names of samples that have to be assigned to the 
  #                           evaluation series.
  #   swarm.parameters: Vector with parameters to be optimized by swarm intelligence. By default the parameters FDR,
  #                     stable genes based on lib size correlation, correlated genes, and ranked genes are included.
  #                     Additional clinical parameters may be added by adding the sample info column names to this vector
  #                     plus the 'variable.to.assess'-vector and a default value to the 'variable.threshold'-vector.
  #   swarm.boundaries: Vector with lower and upper boundaries per variable in swarm parameters, separated
  #                     by a comma. Number of input values should match with the number of swarm parameters provided.
  #   k.variables: Number of (expected) variables/axes in the dataset.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   ruvg.pvalue.threshold.group: Numeric-value representing the lowest p-value that should be 
  #                                be reached by the correlation between the counts and  any variable
  #                                in order to bypass the wanted variable 'group' to be selected.
  #   ruvg.pvalue.threshold.strongest.variable: Numeric-value representing the lowest p-value that should 
  #                                be reached by the correlation between the counts and the specific 
  #                                variable in order to this variable to be assigned to the RUVg axis.
  #   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
  #                         or PValue (FALSE) statistics.
  #   minimum.n.transcripts.biomarkerpanel: Numeric value with minimum number of RNAs to be included in the 
  #                                         biomarker panel.
  #   svm.gamma.range: Numeric value for the range of the grid search for the best gamma parameter in SVM.
  #   svm.cost.range: Numeric value for the range of the grid search for the best cost parameter in SVM.
  #   number.cross.splits: Numeric value with the number of subseriess employed by SVM algorithm for internal tuning.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   n.particles: Numeric-value with number of PSO particles per iteration for classifier development.
  #   n.iterations: Numeric-value with number of PSO iterations in total for classifier development.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #  Returns the best settings selected by PSO.
  #  Also produces output files of the training and optimization process, including the ppso.log, ppso.log, and individual RData
  #  files in the folder outputPSO that contain the trained support vectors.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!percentage.for.training %in% seq(1,100,by=1e-3)){
    stop("percentage.for.training should be within 1 and 100")
  }
  
  if (!percentage.for.evaluation %in% seq(1,100,by=1e-3)){
    stop("percentage.for.evaluation should be within 1 and 100")
  }
  
  if (is.numeric(percentage.for.training) & length(training.samples.provided) > 1){
    print("Both percentage for training and a separate training list provided. The provided training list will be used.")
  }
  
  if (length(training.samples.provided) > 1 & length(evaluation.samples.provided) == 1){
    stop("list of training samples provided but no evaluation samples specified. Please specify evaluation samples.")
  }
  
  if (length(swarm.parameters)*2 != length(swarm.boundaries)){
    stop("number of swarm.parameters and provided swarm.boundaries do not match.")
  }
  
  if (!is.numeric(k.variables)){
    stop("Provide numeric value for k.variables")
  }
  
  if (k.variables < 1){
    stop("k.variables should be 1 or more")
  }
  
  if (all(variable.to.assess %in% colnames(dge$samples)) != TRUE){
    stop("Inputted variables do not match column names of the sample info table.")
  }
  
  if (length(variable.to.assess) != length(variable.threshold)){
    stop("Number of variables in variable.to.assess should match number of thresholds")
  }
  
  if (!is.numeric(ruvg.pvalue.threshold.group)){
    stop("Provide numeric value for ruvg.pvalue.threshold.group")
  }
  
  if (!is.numeric(ruvg.pvalue.threshold.strongest.variable)){
    stop("Provide numeric value for ruvg.pvalue.threshold.strongest.variable")
  }
  
  if (!is.numeric(minimum.n.transcripts.biomarkerpanel)){
    stop("Provide numeric value for minimum.n.transcripts.biomarkerpanel")
  }
  
  if (!is.numeric(number.cross.splits)){
    stop("Provide numeric value for number.cross.splits")
  }
  
  if (!is.numeric(n.particles.gamma.cost.optimization)){
    stop("Provide numeric value for n.particles.gamma.cost.optimization")
  }
  
  if (!is.numeric(n.iterations.gamma.cost.optimization)){
    stop("Provide numeric value for n.iterations.gamma.cost.optimization")
  }
  
  if (!is.numeric(n.particles)){
    stop("Provide numeric value for n.particles")
  }
  
  if (!is.numeric(n.iterations)){
    stop("Provide numeric value for n.iterations")
  }
  
  library(parallel)
  if (number.cores > detectCores()){
    print(paste("Please note that more computer threads have been selected than available with this machine. This may cause problems with start of the PSO process. It is advised to reduce the number of cores to ", detectCores(), sep = ""))
  }
    
  # load required packages
  if (verbose == TRUE){
    print("Load required packages ppso, edgeR, e1071, and RUVSeq")
  } 
  suppressMessages(library(ppso))
  suppressMessages(library(edgeR))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(e1071))
  
  # create subdirectory for PSO-process
  if (!file.exists("pso-enhanced-thromboSeq")){
    dir.create("pso-enhanced-thromboSeq", recursive = T)
  }
  # store current directory and subdirectory
  workDir_main <- getwd()
  setwd("pso-enhanced-thromboSeq/") # transfer to subdirectory
  workDir <- getwd()
  
  if (is.numeric(percentage.for.training) & length(training.samples.provided) == 1){
  # randomly select samples for the training and evaluation series
  # here it is assumed the group size and potential confounding factors 
  # (e.g. age of the individuals) are similar among both groups
  set.seed(1000) # lock random number generator
  series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
    n.samples.training <- round(length(which(
      dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
    ) 
    
    training.samples.subset <- sample(
      colnames(dge)[dge$samples$group == levels(dge$samples$group)[i]],
      size = n.samples.training,
      replace = F
    )
    
    # container
    series <- list()
    series[["training.samples.subset"]] <- training.samples.subset
    series
  }
  training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))  
  write.csv(training.samples, "trainingSamples.csv")
  
  series.evaluation <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
    n.samples.evaluation <- round(length(which(
      dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.evaluation / 100)
    ) 
    
    evaluation.samples.subset <- sample(
      colnames(dge[, dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% training.samples]),
      size = n.samples.evaluation,
      replace = F
    )
    
    # container
    series <- list()
    series[["evaluation.samples.subset"]] <- evaluation.samples.subset
    series
  }
  evaluation.samples <- unlist(lapply(series.evaluation,  function(x){x[["evaluation.samples.subset"]]}))  
  write.csv(evaluation.samples,"evaluationSamples.csv")
  } else {
    training.samples <- training.samples.provided
    evaluation.samples <- evaluation.samples.provided
  }
  
  if (verbose == TRUE){
    print("Samples selected for training series:")
    print(summary((dge[, training.samples])$samples$group))
    print("Samples selected for evaluation series:")
    print(summary((dge[, evaluation.samples])$samples$group))
  } 
  
  # store these particular variables
  k.variables = k.variables
  variable.to.assess = variable.to.assess
  ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group
  ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable
  
  # input parameter boundaries for ppso:
  if ('all' %in% swarm.boundaries){
    # in case rankedTranscripts supplied with 'all', replace by numeric value of all RNAs detected
    swarm.boundaries[swarm.boundaries == 'all'] <- nrow(dge$counts)
  }
  
  # prepare matrix with per variable to optimize the lower and upper boundary value
  parameter.bounds <- matrix(ncol = 2, nrow = length(swarm.parameters))
  rownames(parameter.bounds) <- swarm.parameters
  parameter.bounds[,1] <- swarm.boundaries[seq(1, length(swarm.boundaries), by = 2)]
  parameter.bounds[,2] <- swarm.boundaries[seq(2, length(swarm.boundaries), by = 2)]
  parameter.bounds <- apply(parameter.bounds, 2, as.numeric) # ensure the matrix contains numeric values
  
  # create an RData file with all data in this R session, 
  # employed by ppso for analyzing all individual particles
  save(list = ls(envir = environment(), all.names = TRUE), 
       file = "Pre-PSO-snapshot.RData", envir = environment()) 
  if (!file.exists("outputPSO")){
    dir.create("outputPSO", recursive = T)
  } # separate output folder for individual particles
  
  # run ppso
  # Of note, this ppso-function may require days to complete. The algorithm will read each 
  # second whether the 'slaves' have produced new results and updates automatically the swarming process.
  # Continuous output can be seen in the 'outputPSO'-folder, the individual slave-log-files, 
  # and following each iteration of n particles the ppso.log-file.
  set.seed(1000) # lock random number generator
  resultPSO <- optim_ppso_robust(objective_function        = thrombo.algo,
                                 nslaves                   = number.cores - 1,
                                 number_of_parameters      = nrow(parameter.bounds),
                                 plot_progress             = FALSE,
                                 number_of_particles       = n.particles,
                                 max_number_of_iterations  = n.iterations,
                                 max_number_function_calls = n.particles * n.iterations,
                                 parameter_bounds          = parameter.bounds,
                                 tryCall                   = TRUE, 
                                 verbose                   = if (verbose == TRUE){TRUE} else {FALSE},
                                 lhc_init                  = TRUE, 
                                 wait_complete_iteration   = TRUE,
                                 logfile                   = paste(workDir, "/ppso.log", sep = ""),
                                 projectfile               = paste(workDir, "/ppso.pro", sep = ""),
                                 break_file                = "stopPPSO.txt"
  )
  # end of ppso
  
  # Select the best parameter setting
  logged.PSO.distribution <- read.csv("ppso.log", sep = "\t")
  logged.PSO.distribution <- logged.PSO.distribution[
    order(logged.PSO.distribution$objective_function), 
    ]
  logged.PSO.distribution.Index <- logged.PSO.distribution[
    which(logged.PSO.distribution$objective_function == min(logged.PSO.distribution$objective_function)), 
    ]
  
  set.seed(1000)
  # in case more particles have the same AUC output value, select randomly one as the particle for readout
  if (nrow(logged.PSO.distribution.Index) > 1){
    logged.PSO.distribution.Index <- logged.PSO.distribution.Index[
      sample(1 : nrow(logged.PSO.distribution.Index), size = 1),
      ]
  }
  
  # plot the pso-progress
  # add swarm.parameters to columns of ppso-overview log-file
  colnames(logged.PSO.distribution) <- c("time",swarm.parameters,"objective_function","worker")
  colnames(logged.PSO.distribution.Index) <- c("time",swarm.parameters,"objective_function","worker")
  # for each PSO-optimization parameter plot the input value to 1-AUC-result
  # highlight the ultimate best particle in red
  # generate and store plot in output folder
  if (!file.exists(paste(workDir_main,"/",figureDir,sep=""))){
    dir.create(paste(workDir_main,"/",figureDir,sep=""), recursive = T)
  }
  
  pdf(paste(workDir_main,"/",figureDir,"/PSO_optimizationPlots.pdf",sep=""), paper="a4")
  par(mfrow=c(2,2))
  for (parameter in swarm.parameters){
    plot(logged.PSO.distribution[,parameter],logged.PSO.distribution$objective_function,
         xlim = c(min(logged.PSO.distribution[,parameter]), max(logged.PSO.distribution[,parameter])), 
         ylim = c(0, max(logged.PSO.distribution$objective_function)),
         pch = 20, 
         ylab = "1-AUC", 
         xlab = parameter)
    points(logged.PSO.distribution.Index[,parameter], logged.PSO.distribution.Index$objective_function, 
           col = "red", 
           pch = 19)
  }
  dev.off()
  
  # return the string with the best settings selected by PSO.
  return(logged.PSO.distribution.Index)
}

thrombo.algo <- function(x){
  # thromboSeq classification algorithm to be optimized by PSO. In brief, the classification 
  # algorithm performs data correction and normalization employing only the training series. 
  # Next, it performs ANOVA differential expression analysis of splice junctions. Then, it 
  # identifies highly correlated RNAs and removes those according to the PSO-proposed threshold. 
  # Following, it trains the first SVM algorithm of which the most contributing RNAs are 
  # maintained, as suggested by PSO. Then, an updated SVM-model is trained of which the gamma 
  # and cost setting are also PSO-optimized. Finally, the classification model is stored and the 
  # evaluation series are classified.
  # 
  # Args
  #   x: Vector with values proposed by PSO (ppso-function)
  #
  # Returns:
  #   Inverse AUC-value of the algorithm when evaluation series were classified using the proposed
  #   algorithm threshold settings
  
  # load R environment data
  load(paste(getwd(), "/Pre-PSO-snapshot.RData", sep = ""))
  # load packages
  suppressMessages(library(edgeR, warn.conflicts = F, quietly = T))
  suppressMessages(library(e1071, warn.conflicts = F, quietly = T))
  suppressMessages(library(foreach, warn.conflicts = F, quietly = T))
  suppressMessages(library(ROCR, warn.conflicts = F, quietly = T))
  suppressMessages(library(RUVSeq, warn.conflicts = F, quietly = T))
  suppressMessages(library(caret, warn.conflicts = F, quietly = T))
  suppressMessages(library(ppso, warn.conflicts = F, quietly = T))
  # load functions
  source(paste(workDir_main, "/bin/thromboSeqTools_PreProcessing_2.R", sep = ""))
  source(paste(workDir_main, "/bin/thromboSeqTools_ANOVA.R", sep = ""))
  source(paste(workDir_main, "/bin/thromboSeqTools_PSO.R", sep = ""))
  
  # prepare data according to training and evaluation series
  dgeTraining <- dge[, training.samples]
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.evaluation <- dge$samples[evaluation.samples, "group"]
  
  # ruv.correction
  # if variable.to.assess in PSO optimization, replace default value by PSO-proposed value
  for (variable in variable.to.assess){
    if (variable %in% swarm.parameters){
      variable.threshold[which(variable.to.assess==variable)] <- x[which(variable==swarm.parameters)]
    }
  }
  
  dgeTraining <- perform.RUVg.correction(dge = dge, 
                                         k.variables = k.variables, 
                                         variable.to.assess = variable.to.assess,
                                         variable.threshold = variable.threshold, 
                                         ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                         ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                         training.series.only = T,
                                         training.series.only.samples = training.samples,
                                         verbose = verbose)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$raw.lib.size <- dgeTraining$samples$lib.size
  dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
  
  # perform TMM normalization
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
  
  # Likelihood-ratio test modified for thromboSeq (ANOVA)
  dgeTraining$ruv.counts <- dgeTraining$counts
  dgeTraining$counts <- dgeTraining$raw.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$raw.lib.size
  thromboSeq.anova <- anovaLikeEdgeRthromboSeq(dgeTraining[, training.samples], 
                                               method = "TMM",
                                               normalize.on.training.series = TRUE, 
                                               samples.for.training = training.samples,
                                               ref.column = dgeTraining$refSample)
  
  # select FDR threshold. When FDR/fdr was included as a PSO optimization variable, select PSO-proposed value
  if (any(c("fdr","FDR") %in% swarm.parameters)){
    fdr <- x[which(swarm.parameters %in% c("fdr","FDR"))]
  } else {
    fdr <- 0.05 # in case fdr/FDR was not provided by user
  }
  
  # select biomarker panel using either FDR or p-value statistics
  if (select.biomarker.FDR == TRUE){
    selected.transcripts <- rownames(thromboSeq.anova)[
      thromboSeq.anova$FDR < fdr & 
        thromboSeq.anova$logCPM > 3 & 
        thromboSeq.anova$chromosome_name %in% c(1:22, "X")
      ]  
  } else {
    selected.transcripts <- rownames(thromboSeq.anova)[
      thromboSeq.anova$PValue < fdr & 
        thromboSeq.anova$logCPM > 3 & 
        thromboSeq.anova$chromosome_name %in% c(1:22, "X")
      ]  
  }
  
  # only continue when biomarker panel size is more than provided threshold
  if (length(selected.transcripts) > minimum.n.transcripts.biomarkerpanel){
    # remove transcripts with a Pearson's correlation to any other transcripts over the PSO-proposed threshold
    correlation.matrix <- cor(t(normalized.counts[selected.transcripts, training.samples]))
    if ("correlatedTranscripts" %in% swarm.parameters){
      correlated.transcripts <- x[which(swarm.parameters == "correlatedTranscripts")]
    } else {
      correlated.transcripts <- 1.0 # in case this variable was not provided by user
    }
    highly.correlated <- colnames(correlation.matrix)[findCorrelation(correlation.matrix, cutoff = correlated.transcripts)]
    
    # select biomarker panel using either FDR or p-value statistics
    if (select.biomarker.FDR == TRUE){
      selected.transcripts <- rownames(thromboSeq.anova)[
        thromboSeq.anova$FDR < fdr & thromboSeq.anova$logCPM > 3 & thromboSeq.anova$chromosome_name %in% c(1:22, "X") &
          !thromboSeq.anova$ensembl_gene_id %in% highly.correlated
        ]  
    } else {
      selected.transcripts <- rownames(thromboSeq.anova)[
        thromboSeq.anova$PValue < fdr & thromboSeq.anova$logCPM > 3 & thromboSeq.anova$chromosome_name %in% c(1:22, "X") &
          !thromboSeq.anova$ensembl_gene_id %in% highly.correlated
        ]  
    }
    
    # only continue when biomarker panel size is more than provided threshold
    if (length(selected.transcripts) > minimum.n.transcripts.biomarkerpanel){
      # first SVM-model, with grid search for gamma and cost
      tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                            y           = real.groups.training,
                            gamma       = svm.gamma.range,
                            cost        = svm.cost.range,
                            tunecontrol = tune.control(cross = number.cross.splits),
                            probability = TRUE
      )
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      # rank transcripts (features) of the biomarker panel according to their relative contribution to the SVM model
      if (nlevels(dgeTraining$samples$group) == 2){
        # binary feature ranking algorithm
        svm.ranking <- svmrfeFeatureRanking(
          t(normalized.counts[selected.transcripts, training.samples]), 
          real.groups.training, tuned.svm.model$gamma, tuned.svm.model$cost)
      } else {
        # multiclass feature ranking algorithm
        svm.ranking <- svmrfeFeatureRankingForMulticlass(
          t(normalized.counts[selected.transcripts, training.samples]), 
          real.groups.training, tuned.svm.model$gamma, tuned.svm.model$cost)
      }
      
      # employ PSO-proposed threshold and set final spliced RNA biomarker panel of this PSO particle
      if ("rankedTranscripts" %in% swarm.parameters){
        ranked.transcripts <- round(x[which(swarm.parameters == "rankedTranscripts")], digits = 0)
      } else {
        ranked.transcripts <- nrow(dgeTraining$counts) # in case rankedTranscripts was not provided by user
      }
      selected.transcripts <- (cbind(colnames(t(normalized.counts[selected.transcripts, training.samples])), svm.ranking)[
        which(as.numeric(cbind(colnames(t(normalized.counts[selected.transcripts, training.samples])), svm.ranking)[,2]) <= as.numeric(ranked.transcripts)),
        ])[,1]
      
      # second SVM-model, re-train with adjusted biomarker panel and employ a grid search for gamma and cost
      tuned.svm <- tune.svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                            y           = real.groups.training,
                            gamma       = svm.gamma.range,
                            cost        = svm.cost.range,
                            tunecontrol = tune.control(cross = number.cross.splits),
                            probability = TRUE
      )
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      ## optimize gamma and cost by a second PSO-optimization algorithm
      # employ the training series for SVM algorithm training, and the evaluation series for optimization
      # make sure the values for cost and gamma are not infinite
      if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
        # if necessary - broaden the range of cost and gamma
        if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
          svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
        }
        if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
          svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
        }
        svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                             svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
        svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                            svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])

        # input parameters:
        parameter.bounds.gamma.cost <- matrix(ncol = 2, nrow = 2)
        rownames(parameter.bounds.gamma.cost) <- c("gamma", "cost")
        parameter.bounds.gamma.cost[, 1] <- c(svm.gamma.range[1], svm.cost.range[1])
        parameter.bounds.gamma.cost[, 2] <- c(svm.gamma.range[2], svm.cost.range[2])

        # PSO
        set.seed(2000)

        result.internal.gamma.cost <- optim_pso(objective_function        = if (nlevels(dgeTraining$samples$group) == 2){
                                                                                  thrombo.svm.gamma.cost
                                                                            } else {
                                                                                  thrombo.svm.gamma.cost.multiclass
                                                                            },
                                                number_of_parameters      = nrow(parameter.bounds.gamma.cost),
                                                plot_progress             = FALSE,
                                                number_of_particles       = n.particles.gamma.cost.optimization,
                                                max_number_of_iterations  = n.iterations.gamma.cost.optimization,
                                                max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
                                                parameter_bounds          = parameter.bounds.gamma.cost,
                                                tryCall                   = TRUE,
                                                verbose                   = FALSE,
                                                lhc_init                  = FALSE,
                                                wait_complete_iteration   = TRUE,
                                                logfile                   = NULL,
                                                projectfile               = NULL,
                                                break_file                = "stopPPSO.txt"
        )

        # employ PSO proposed gamma and cost parameters
        tuned.svm <- svm(x           = t(normalized.counts[selected.transcripts, training.samples]),
                         y           = real.groups.training,
                         gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                         cost        = as.numeric(result.internal.gamma.cost$par[2]),
                         tunecontrol = tune.control(cross = number.cross.splits),
                         probability = TRUE
        )

        # extract best model
        tuned.svm.model <- tuned.svm
      }
      
      # prepare counts for sample prediction
      normalized.counts.prediction <- normalized.counts[selected.transcripts, evaluation.samples]
      
      # store data
      dgeTraining$biomarker.transcripts <- selected.transcripts
      dgeTraining$tuned.svm.model <- tuned.svm.model
      
      # this file is stored in the outputPSO-folder and can be employed for evaluation and validation 
      # if this particle is selected as the 'best option'
      save(dgeTraining, 
           file = paste("outputPSO/", paste(x, collapse = "-"), ".RData", sep = ""))
      
      # predict evaluation series
      prediction.class <- predict(tuned.svm.model,
                                  newdata = t(normalized.counts.prediction), 
                                  probability = TRUE)
      confusion.matrix <- table(prediction.class, real.groups.evaluation)
      confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
      
      # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
      if (nlevels(dgeTraining$samples$group) == 2){
        # create classification overview
        svm.summary <- data.frame(
          sampleName = attributes(prediction.class)$names,
          predicted = as.character((prediction.class)[1:length(prediction.class)]),
          real = real.groups.evaluation
        )
        svm.summary <- cbind(svm.summary, 
                             data.frame(attributes(prediction.class)$probabilities[,
                               which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                             data.frame(attributes(prediction.class)$probabilities[,
                               which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
        )
        colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
        # ROC
        rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                            svm.summary[, 3], 
                            label.ordering = levels(dgeTraining$samples$group)
                            )
        perfa <- performance(rocra, "tpr", "fpr")
        if (verbose == TRUE){
          print(paste("AUC Evaluation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                      sep = ""))
        }
        roc.summary <- data.frame(
          cutOffs = unlist(attributes(rocra)$cutoffs),
          tp = unlist(attributes(rocra)$tp),
          tn = unlist(attributes(rocra)$tn),
          fp = unlist(attributes(rocra)$fp),
          fn = unlist(attributes(rocra)$fn),
          accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
            (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
          xValues = unlist(attributes(perfa)$x.values),
          yValues = unlist(attributes(perfa)$y.values)
        )
        roc.optimal.accuracy <- max(roc.summary$accuracy)
        
        # metrics for readout
        AUC <- 1 - attributes(performance(rocra, 'auc'))$y.values[[1]]
        return(AUC)
      } else {
        AUC <- 1 - confusion.matrix.evaluated$diag
        return(AUC)
      } 
    } else {
      AUC <- Inf
      return(AUC)
    }
  } else {
    AUC <- Inf
    return(AUC)
  }
}

# Functions to identify the relative contribution of features (transcripts)
# to the SVM model for binary (svmrfeFeatureRanking-function) and multiclass
# (svmrfeFeatureRankingForMulticlass-function) algorithm development.
# the svm.weights function is employed in the other two individual functions
# adapted from http://www.bdmg.com.ar/wp-content/uploads/2011/11/SVM_RFE_R_implementation.pdf
svmrfeFeatureRanking <- function(x, 
                                 y,
                                 bestGamma,
                                 bestCost){
  n = ncol(x)
  survivingFeaturesIndexes = seq(1 : n)
  featureRankedList = vector(length = n)
  rankedFeatureIndex = n
  while(length(survivingFeaturesIndexes) > 0){
    # train the support vector machine
    svmModel <- svm(x[, survivingFeaturesIndexes], y, cost = bestCost, gamma = bestGamma,
                    cachesize = 1000, scale = T, type = "C-classification")
    # compute the weight vector
    w = t(svmModel$coefs)%*%svmModel$SV
    #compute ranking criteria
    rankingCriteria = w * w
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix
    # update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
  }
  return (featureRankedList)
}

svmrfeFeatureRankingForMulticlass <- function(x,
                                              y,
                                              bestCost,
                                              bestGamma){ 
  n = ncol(x) 
  survivingFeaturesIndexes = seq(1 : n) 
  featureRankedList = vector(length = n) 
  rankedFeatureIndex = n 
  while(length(survivingFeaturesIndexes) > 0){ 
    #train the support vector machine 
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = bestCost, gamma = bestGamma, cachesize = 1000,  
                   scale = T, type = "C-classification")
    #compute the weight vector 
    multiclassWeights = svm.weights(svmModel) 
    #compute ranking criteria 
    multiclassWeights = multiclassWeights * multiclassWeights 
    rankingCriteria = 0 
    for(i in 1:ncol(multiclassWeights)) 
      rankingCriteria[i] = mean(multiclassWeights[, i]) 
    #rank the features 
    (ranking = sort(rankingCriteria, index.return = TRUE)$ix) 
    # update feature ranked list 
    (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]) 
    rankedFeatureIndex = rankedFeatureIndex - 1
    # eliminate the feature with smallest ranking criterion 
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]]) 
  } 
  return (featureRankedList)
} 

svm.weights <- function(model){ 
  w = 0 
  if(model$nclasses == 2){ 
    w = t(model$coefs)%*%model$SV 
  } else {    
    #when we deal with OVO svm classification 
    ## compute start-index 
    start <- c(1, cumsum(model$nSV) + 1) 
    start <- start[-length(start)] 
    calcw <- function (i,j) { 
      ## ranges for class i and j: 
      ri <- start[i] : (start[i] + model$nSV[i] - 1) 
      rj <- start[j] : (start[j] + model$nSV[j] - 1) 
      ## coefs for (i,j): 
      coef1 <- model$coefs[ri, j-1] 
      coef2 <- model$coefs[rj, i] 
      ## return w values: 
      w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,] 
      return(w) 
    } 
    W=NULL 
    for (i in 1 : (model$nclasses - 1)){ 
      for (j in (i + 1) : model$nclasses){ 
        wi=calcw(i,j) 
        W=rbind(W,wi) 
      } 
    } 
    w=W 
  } 
  return(w) 
} 

# Svm functions, binary and multiclass.
# Optimize towards AUC in binary classifier, and accuracy in multiclass classifier

# Function do, required by thrombo.svm.gamma.cost and thrombo.svm.gamma.cost.multiclass.
# determines whether the value outputted by do can be handled by R.
do <- function(x){
  2^x
}

thrombo.svm.gamma.cost <- function(x) {
  # Function to improve gamma and cost settings in trained SVM model, binary comparison.
  # 
  # Args
  #   x: Vector with PSO-proposed gamma and cost values.
  #
  # Returns:
  #   Inverse AUC-value of the algorithm when evaluation series were classified using the proposed
  #   algorithm threshold settings.
  
  set.seed(1000) # lock randomness
  gammaRange <- x[1] # first value in x-vector is gamma value
  costRange <- x[2] # second value in x-vector is gamma value
  internalTraining <- colnames(normalized.counts) # training series
  internalTesting <- evaluation.samples # evaluation series
  realGroupsTesting <- droplevels(dgeTraining$samples[internalTesting, "group"])
  if(!do(costRange) == Inf){ # ensure the value for cost does not result in non-handible values
    if(!do(gammaRange) == Inf){ # ensure the value for gamma does not result in non-handible values
      # train SVM model with the proposed gamma and cost values
      m <- svm(t(normalized.counts)[internalTraining, selected.transcripts],
               droplevels(dgeTraining$samples[internalTraining, "group"]), 
               kernel = "radial", 
               cost = do(costRange), 
               gamma = do(gammaRange),
               type = "C-classification", 
               probability = TRUE, 
               cachesize = 500, 
               scale = T)
      # predict the evaluation series in the trained SVM model
      yscore <- attr(predict(m, t(normalized.counts)[internalTesting, selected.transcripts], probability= T), "probabilities")[, levels(realGroupsTesting)[1]]
      # construct AUC-value
      rocra <- prediction(as.numeric(as.character(yscore)), realGroupsTesting)
      perfa <- performance(rocra, "tpr", "fpr")
      AUC <- 1 - attributes(performance(rocra, 'auc'))$y.values[[1]]
      return(AUC)
    } else { 
      AUC <- Inf
      return(AUC)
    }
  } else { 
    AUC <- Inf
    return(AUC)
  }
}
thrombo.svm.gamma.cost.multiclass <- function(x) {
  # Function to improve gamma and cost settings in trained SVM model, multiclass comparison.
  # 
  # Args
  #   x: Vector with PSO-proposed gamma and cost values.
  #
  # Returns:
  #   Inverse accuracy of the algorithm when evaluation series were classified using the proposed
  #   algorithm threshold settings.
  
  set.seed(1000) # lock randomness
  gammaRange <- x[1] # first value in x-vector is gamma value
  costRange <- x[2] # second value in x-vector is gamma value
  internalTraining <- colnames(normalized.counts) # training series
  internalTesting <- evaluation.samples # evaluation series
  realGroupsTesting <- droplevels(dgeTraining$samples[internalTesting, "group"])
  if(!do(costRange) == Inf){# ensure the value for cost does not result in non-handible values
    if(!do(gammaRange) == Inf){ # ensure the value for gamma does not result in non-handible values
      # train SVM model with the proposed gamma and cost values
      m <- svm(t(normalizedCounts)[internalTraining, selected.transcripts],
               droplevels(dgeTraining$samples[internalTraining,"group"]), 
               kernel = "radial", 
               cost = do(costRange), 
               gamma = do(gammaRange),
               type = "C-classification", 
               probability = TRUE, 
               cachesize = 500, 
               scale = T)
      # predict the evaluation series in the trained SVM model and calculate inverse accuracy
      Acc <- 1 - classAgreement(table(predict(m, t(normalizedCounts)[
        internalTesting,selectedTranscripts_svmRanking], probability = T), realGroupsTesting), match.names = T)$diag
      return(Acc)
    } else { 
      Acc <- Inf
      return(Acc)
    }
  } else { 
    Acc <- Inf
    return(Acc)
  }
}

thromboSeqPSO.readout <- function(dge = dgeIncludedSamples,
                                  pso.output.log.file = "ppso.log",
                                  pso.output.pro.file = "ppso.pro",
                                  readout.training = TRUE,
                                  readout.evaluation = TRUE,
                                  readout.validation = TRUE,
                                  replace.counts.validation = 0,
                                  R.snapshot = "Pre-PSO-snapshot.RData",
                                  figureDir = "figureOutputFolder",
                                  number.cores = 2,
                                  verbose = TRUE){
  # Perform readout of the PSO-enhanced thromboSeq classification algorithm. 
  #
  # Args:
  #   pso.output.log.file: String with log-file name of the PSO-optimization process.
  #   pso.output.pro.file: String with pro-file name of the PSO-optimization process.
  #   readout.training: Whether or not to perform LOOCV analysis of training series (TRUE/FALSE).
  #   readout.evaluation: Whether or not to classify samples of the evaluation series (TRUE/FALSE).
  #   readout.validation: Whether or not to classify samples of the validation series (TRUE/FALSE).
  #   replace.counts.validation: Numeric-value indicating the number of reads counts (0 - value)
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Vector with best parameter settings. Also returns in working directory table with algorithm metrics. 
  #   In addition, outputs RData output files with specific data regarding the training, evaluation and 
  #   validation series, and a pdf-plot of the ROC-curve.
  
  if (all(!is.numeric(replace.counts.validation), !is.nan(replace.counts.validation))){
    stop("provide numeric value or NaN for replace.counts.validation")
  }
  
  # load required packages
  suppressMessages(library(foreach))
  suppressMessages(library(RUVSeq))
  suppressMessages(library(e1071))
  suppressMessages(library(ppso))
  suppressMessages(library(ROCR))
  suppressMessages(library(pROC))
  suppressMessages(library(reshape))
  suppressMessages(library(ggplot2))
  
  # Select the best parameter setting
  logged.PSO.distribution <- read.csv(pso.output.log.file, sep = "\t")
  logged.PSO.distribution <- logged.PSO.distribution[
    order(logged.PSO.distribution$objective_function), 
    ]
  logged.PSO.distribution.Index <- logged.PSO.distribution[
    which(logged.PSO.distribution$objective_function == min(logged.PSO.distribution$objective_function)), 
    ]
  
  set.seed(1000)
  # in case more particles have the same AUC output value, select randomly one as the particle for readout
  if (nrow(logged.PSO.distribution.Index) > 1){
    logged.PSO.distribution.Index <- logged.PSO.distribution.Index[
      sample(1 : nrow(logged.PSO.distribution.Index), size = 1),
      ]
  }
  best.selection <- paste(logged.PSO.distribution.Index[, c(seq(2,ncol(logged.PSO.distribution.Index)-2, by=1))], collapse = "-") # collapse data
  
  if (verbose == TRUE){
    print(paste("Best selection: ", best.selection, sep = ""))
  }
  
  if (readout.training == TRUE){
    results.classification.training <- thrombo.algo.classify.training.set(dge = dgeIncludedSamples,
                                                                          best.particle = best.selection,
                                                                          R.snapshot = "Pre-PSO-snapshot.RData",
                                                                          number.cores = 2)
    save(results.classification.training, file = "results.classification.training.RData")
  }
  
  if (readout.evaluation == TRUE){
    results.classification.evaluation <- thrombo.algo.classify.evaluation.set(dge = dgeIncludedSamples,
                                                                              best.particle = best.selection,
                                                                              R.snapshot = "Pre-PSO-snapshot.RData")
    save(results.classification.evaluation, file = "results.classification.evaluation.RData")
  }
  
  if (readout.validation == TRUE){
    results.classification.validation <- thrombo.algo.classify.validation.set(dge = dgeIncludedSamples,
                                                                              best.particle = best.selection,
                                                                              R.snapshot = "Pre-PSO-snapshot.RData")
    save(results.classification.validation, file = "results.classification.validation.RData")
  }

  if (verbose == TRUE){
    print("summarize data")
  }
  
  # summarize data in a matrix
  matrix <- matrix(nrow = 3, ncol = 4)
  rownames(matrix) <- c("Training", "Evaluation", "Validation")
  colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy")
  if (readout.training == TRUE){
  matrix[1, ] <- c(results.classification.training$number.samples,
                   results.classification.training$AUCorDiagonal,
                  paste(round(results.classification.training$roc.95ci$ci[1], digits = 2), "-",
                        round(results.classification.training$roc.95ci$ci[3], digits = 2), sep = ""),
                  results.classification.training$roc.optimal.accuracy)
  }
  if (readout.evaluation == TRUE){
  matrix[2,] <- c(length(results.classification.evaluation$samples.for.evaluation),
                  results.classification.evaluation$AUCorDiagonal,
                  paste(round(results.classification.evaluation$roc.95ci$ci[1], digits = 2), "-",
                        round(results.classification.evaluation$roc.95ci$ci[3], digits = 2), sep = ""),
                  results.classification.evaluation$roc.optimal.accuracy)
  }
  if (readout.validation == TRUE){
  matrix[3,] <- c(length(results.classification.validation$samples.for.evaluation),
                  results.classification.validation$AUCorDiagonal,
                  paste(round(results.classification.validation$ci.roc$ci[1], digits = 2), "-",
                        round(results.classification.validation$ci.roc$ci[3], digits = 2), sep = ""),
                  results.classification.validation$roc.optimal.accuracy)
  }
  write.csv(matrix, file = "ROCcurve_Metrics.csv")
  
  # if binary classifier print ROC-curve. For multiclass comparisons print confusion matrices
  if (nlevels(dge$samples$group) == 2){
    # Plot ROC-curve
    # check for correct output-directory for ROC Curve
    currentDir <- getwd()
    if (file.exists(figureDir) == TRUE){
      pdf(paste(figureDir, "/ROCcurve.pdf", sep = ""))  
    } else {
      setwd('..') # move to (previous) mother directory
      if (file.exists(figureDir) == TRUE){
        pdf(paste(figureDir, "/ROCcurve.pdf", sep = ""))  
      } else {
        dir.create(figureDir)
        pdf(paste(figureDir, "/ROCcurve.pdf", sep = ""))  
      }
    }
    if (readout.training == TRUE){
      plot(results.classification.training$perfa, 
           lwd = 4, 
           col = "grey")
      par(new=T)
    }
    if (readout.evaluation == TRUE){
      plot(results.classification.evaluation$perfa, 
           lwd = 4, 
           col = "#B03B3D")
      par(new=T)
    }
    if (readout.validation == TRUE){
      plot(results.classification.validation$perfa, 
           lwd = 4, 
           col = "#3C66A6")
      par(new=T)
    }
    dev.off()
    setwd(currentDir)
  } else {
    # confusion matrices
    # check for correct output-directory for confusion matrices
    currentDir <- getwd()
    if (!file.exists(figureDir) == TRUE){
      setwd('..') # move to (previous) mother directory
      if (!file.exists(figureDir) == TRUE){
        dir.create(figureDir)
      }
    }
    if (readout.training == TRUE){
      confusion.matrix <- as.data.frame(
                              cast(results.classification.training$svm.summary, 
                                   predicted.group ~ real.group, 
                                   length, 
                                   value = "sample.ID"))
      rownames(confusion.matrix) <- confusion.matrix$predicted
      confusion.matrix <- confusion.matrix[,-1]
      lev <- sort(unique(c(colnames(confusion.matrix), rownames(confusion.matrix))))
      confusion.matrix <- confusion.matrix[lev, lev]
      colnames(confusion.matrix) <- lev
      rownames(confusion.matrix) <- lev
      confusion.matrix <- as.matrix(confusion.matrix)
      confusion.matrix[is.na(confusion.matrix)] <- 0
      melted.confusion.matrix <- as.data.frame(confusion.matrix)
      melted.confusion.matrix$predicted <- rownames(melted.confusion.matrix)
      melted.confusion.matrix <- melt(melted.confusion.matrix, id.vars = "predicted")
      colnames(melted.confusion.matrix) <- c("predicted", "real", "frequency")
      
      tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
        geom_tile(aes(fill = frequency)) +
        scale_fill_continuous(low = "white", high = "red") + 
        geom_text(aes(label = frequency)) +
        ggtitle("SVM classification results") + 
        labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
        theme_bw() + 
        theme(legend.position = "top")
      
      pdf(paste(figureDir, "TrainingSeriesConfusionMatrix.pdf", sep = ""), paper = "a4")
      print(tiles)
      dev.off()
    }
    if (readout.evaluation == TRUE){
      confusion.matrix <- as.data.frame(
        cast(results.classification.evaluation$svm.summary, 
             predicted.group ~ real.group, 
             length, 
             value = "sampleName"))
      rownames(confusion.matrix) <- confusion.matrix$predicted
      confusion.matrix[order(colnames(confusion.matrix)[2:ncol(confusion.matrix)]),]
      confusion.matrix <- confusion.matrix[,-1]
      lev <- sort(unique(c(colnames(confusion.matrix), rownames(confusion.matrix))))
      confusion.matrix <- confusion.matrix[lev, lev]
      colnames(confusion.matrix) <- lev
      rownames(confusion.matrix) <- lev
      confusion.matrix <- as.matrix(confusion.matrix)
      confusion.matrix[is.na(confusion.matrix)] <- 0
      melted.confusion.matrix <- as.data.frame(confusion.matrix)
      melted.confusion.matrix$predicted <- rownames(melted.confusion.matrix)
      melted.confusion.matrix <- melt(melted.confusion.matrix, id.vars = "predicted")
      colnames(melted.confusion.matrix) <- c("predicted", "real", "frequency")
      tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
        geom_tile(aes(fill = frequency)) +
        scale_fill_continuous(low = "white", high = "red") + 
        geom_text(aes(label = frequency)) +
        ggtitle("SVM classification results") + 
        labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
        theme_bw() + 
        theme(legend.position = "top")
      pdf(paste(figureDir, "EvaluationSeriesConfusionMatrix.pdf", sep = ""), paper = "a4")
      print(tiles)
      dev.off()
    }
    if (readout.validation == TRUE){
      confusion.matrix <- as.data.frame(
        cast(results.classification.validation$svm.summary, 
             predicted.group ~ real.group, 
             length, 
             value = "sampleName"))
      rownames(confusion.matrix) <- confusion.matrix$predicted
      confusion.matrix <- confusion.matrix[,-1]
      lev <- sort(unique(c(colnames(confusion.matrix), rownames(confusion.matrix))))
      confusion.matrix <- confusion.matrix[lev, lev]
      colnames(confusion.matrix) <- lev
      rownames(confusion.matrix) <- lev
      confusion.matrix <- as.matrix(confusion.matrix)
      confusion.matrix[is.na(confusion.matrix)] <- 0
      melted.confusion.matrix <- as.data.frame(confusion.matrix)
      melted.confusion.matrix$predicted <- rownames(melted.confusion.matrix)
      melted.confusion.matrix <- melt(melted.confusion.matrix, id.vars = "predicted")
      colnames(melted.confusion.matrix) <- c("predicted", "real", "frequency")
      
      tiles <- ggplot(melted.confusion.matrix, aes(x = real, y = predicted)) +
        geom_tile(aes(fill = frequency)) +
        scale_fill_continuous(low = "white", high = "red") + 
        geom_text(aes(label = frequency)) +
        ggtitle("SVM classification results") + 
        labs(x = "Real group", y = "Predicted group", fill = "Frequency") +
        theme_bw() + 
        theme(legend.position = "top")
      pdf(paste(figureDir, "ValidationSeriesConfusionMatrix.pdf", sep = ""), paper = "a4")
      print(tiles)
      dev.off()
    }
    
    setwd(currentDir)
  }
  return(best.selection)
}

thromboSeqPSO.controls <- function(thromboSeqPSO.shuffled = TRUE,
                                   thromboSeqPSO.iterations = TRUE,
                                   n.shuffled = 1000,
                                   n.iterations = 1000,
                                   best.particle.input = thromboPSOreadout,
                                   replace.counts.validation = 0,
                                   n.particles.gamma.cost.optimization = 50,
                                   n.iterations.gamma.cost.optimization = 4,
                                   R.snapshot = "Pre-PSO-snapshot.RData",
                                   number.cores = 2,
                                   verbose = TRUE) {
  # Perform shuffled labels and training series iterations for PSO-enhanced thromboSeq classifier.
  #
  # Args:
  #   thromboSeqPSO.shuffled: TRUE/FALSE whether or not to perform shuffled labels experiment in which the group 
  #                           labels of samples in the training series are randomly assigned. Indicates specificity 
  #                           of the developed classification algorithm.
  #   thromboSeqPSO.iterations: TRUE/FALSE whether or not to perform training series iterations experiment, in which 
  #                             the samplesassigned to the training series are randomly shuffled. Indicates the sensitivity 
  #                             of the developed classification algorithm
  #   n.shuffled: Vector with number of iterations for thromboSeqPSO.shuffled.
  #   n.iterations: Vector with number of iterations for thromboSeqPSO.iterations.
  #   best.particle.input: Vector with the merged best.selection values.
  #   replace.counts.validation: Numeric-value indicating the number of reads counts (0 - value)
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Table with algorithm metrics.
  
  if (thromboSeqPSO.shuffled == TRUE) {
    if (verbose == TRUE){
      print("start shuffled labels")
    }
    
    if (!file.exists("outputShuffled")) {
      dir.create("outputShuffled", recursive = T)
    }
    
    # loop the training process i times (total number of iterations) and perform in each iteration algorithm
    # training and validation. Results are stored in a container and summarized afterwards.
    suppressMessages(library(foreach))
    suppressMessages(library(doMC))
    registerDoMC(cores = number.cores)
    shuffled.loop <- foreach(i = 1 : n.shuffled) %dopar% {
      
      results.classification.evaluation.shuffled <- thrombo.algo.classify.evaluation.shuffled(dge = dgeIncludedSamples,
                                                                                              best.particle = best.particle.input,
                                                                                              shuffle.iteration = i)
      save(results.classification.evaluation.shuffled, 
           file = paste("outputShuffled/", i, "-evaluation.RData", sep = ""))
      results.classification.validation.shuffled <- thrombo.algo.classify.validation.shuffled(dge = dgeIncludedSamples,
                                                                                              n.to.impute = replace.counts.validation,
                                                                                              best.particle = best.particle.input,
                                                                                              shuffle.iteration = i)
      save(results.classification.validation.shuffled, 
           file = paste("outputShuffled/", i, "-validation.RData", sep = ""))
      
      results <- list()
      results[["i"]] <- i
      results[["evaluation.acc"]] <- results.classification.evaluation.shuffled$roc.optimal.accuracy
      results[["evaluation.AUC"]] <- results.classification.evaluation.shuffled$AUCorDiagonal
      results[["validation.acc"]] <- results.classification.validation.shuffled$roc.optimal.accuracy
      results[["validation.AUC"]] <- results.classification.validation.shuffled$AUCorDiagonal
      results
    }
    
    # summarize
    output.shuffled <- data.frame(
      i = unlist(lapply(shuffled.loop, function(x){x[["i"]]})),
      evaluation.acc = unlist(lapply(shuffled.loop, function(x){x[["evaluation.acc"]]})),
      evaluation.AUC = unlist(lapply(shuffled.loop, function(x){x[["evaluation.AUC"]]})),
      validation.acc = unlist(lapply(shuffled.loop, function(x){x[["validation.acc"]]})),
      validation.AUC = unlist(lapply(shuffled.loop, function(x){x[["validation.AUC"]]}))
    )
  }
  
  if (thromboSeqPSO.iterations == TRUE) {
    if (verbose == TRUE){
      print("start training set iterations")
    }
    
    if (!file.exists("outputIterations")) {
      dir.create("outputIterations", recursive = T)
    }
    
    n.iterations.datasets <- n.iterations
    
    # load R snapshot data
    load(R.snapshot)
   
     # loop the training process i times (total number of iterations) and perform in each iteration algorithm
    # training and validation. Results are stored in a container and summarized afterwards.
    suppressMessages(library(foreach))
    suppressMessages(library(doMC))
    registerDoMC(cores = number.cores)
    iterations.loop <- foreach(iter = 1 : n.iterations.datasets) %dopar% {
      # reset seed for random selection of training-evaluation series
      rm(list=".Random.seed", envir = globalenv())
      
      # load particle
      load(paste("outputPSO/", best.particle.input, ".RData", sep = ""))
      dgeParticle <- dgeTraining
      
      # randomly select samples for the training and evaluation series.
      # here it is assumed the group size and potential confounding factors 
      # (e.g. age of the individuals) are similar among both groups
      series.training <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
        n.samples.training <- round(length(which(
          dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.training / 100)
        ) 
        
        training.samples.subset <- sample(
          colnames(dge)[dge$samples$group == levels(dge$samples$group)[i]],
          size = n.samples.training,
          replace = F
        )
        
        # container
        series <- list()
        series[["training.samples.subset"]] <- training.samples.subset
        series
      }
      training.samples <- unlist(lapply(series.training, function(x){x[["training.samples.subset"]]}))  
      # store training series
      write.csv(training.samples, 
                paste("outputIterations/", iter, "-trainingSamples_subsampling.csv", sep = ""))
      
      series.evaluation <- foreach(i = 1 : nlevels(dge$samples$group)) %do% {
        n.samples.evaluation <- round(length(which(
          dge$samples$group == levels(dge$samples$group)[i])) * (percentage.for.evaluation / 100)
        ) 
        
        evaluation.samples.subset <- sample(
          colnames(dge[, dge$samples$group == levels(dge$samples$group)[i] & !colnames(dge) %in% training.samples]),
          size = n.samples.evaluation,
          replace = F
        )
        
        # container
        series <- list()
        series[["evaluation.samples.subset"]] <- evaluation.samples.subset
        series
      }
      evaluation.samples <- unlist(lapply(series.evaluation,  function(x){x[["evaluation.samples.subset"]]}))  
      # store evaluation series
      write.csv(evaluation.samples,
                paste("outputIterations/", iter, "-evaluationSamples_subsampling.csv", sep = ""))
      # save R environment
      save(list = ls(envir = environment(), all.names = TRUE), 
           file = paste("outputIterations/", iter, "-Pre-PSO-like.RData", sep = ""), 
           envir = environment()) 
      
      # select series
      real.groups.training <- dge$samples[training.samples, "group"]
      real.groups.prediction <- dge$samples[evaluation.samples, "group"]
      
      dgeTraining <- perform.RUVg.correction(dge = dge, 
                                             k.variables = k.variables, 
                                             variable.to.assess = variable.to.assess,
                                             variable.threshold = variable.threshold, 
                                             ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                             ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                             training.series.only = TRUE, 
                                             training.series.only.samples = training.samples,
                                             verbose = verbose
                                             )
      dgeTraining$counts <- dgeTraining$ruv.counts
      dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
      # TMM-normalisation
      dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                               normalize.on.training.series = TRUE, 
                                               samples.for.training = training.samples,
                                               ref.sample.readout = FALSE) # calculate normalization factors
      dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                               normalize.on.training.series = TRUE, 
                                               samples.for.training = training.samples,
                                               ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
      dgeTraining$samples <- droplevels(dgeTraining$samples)
      
      # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
      normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
      
      # first SVM-model, with grid search for gamma and cost
      tuned.svm <- tune.svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                            y           = real.groups.training,
                            gamma       = svm.gamma.range,
                            cost        = svm.cost.range,
                            tunecontrol = tune.control(cross = number.cross.splits),
                            probability = TRUE
      )
      # extract best model
      tuned.svm.model <- tuned.svm[["best.model"]]
      
      ## optimize gamma and cost by a second PSO-optimization algorithm
      # employ the training series for SVM algorithm training, and the evaluation series for optimization
      # make sure the values for cost and gamma are not infinite
      if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
        # if necessary - broaden the range of cost and gamma
        if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
          svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
        }
        if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
          svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
        }
        svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                             svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
        svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                            svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])
        
        # input parameters:
        parameter.bounds.gamma.cost = matrix(ncol = 2, nrow = 2)
        rownames(parameter.bounds.gamma.cost) <- c("gamma","cost")
        parameter.bounds.gamma.cost[,1] <- c(svm.gamma.range[1],svm.cost.range[1])
        parameter.bounds.gamma.cost[,2] <- c(svm.gamma.range[2],svm.cost.range[2])
        
        # PSO
        set.seed(2000)
        selected.transcripts <- dgeParticle$biomarker.transcripts
        result.internal.gamma.cost <- optim_pso(objective_function = if (nlevels(dgeTraining$samples$group) == 2){
          thrombo.svm.gamma.cost
        } else {
          thrombo.svm.gamma.cost.multiclass
        },
        number_of_parameters      = nrow(parameter.bounds.gamma.cost),
        plot_progress             = FALSE,
        number_of_particles       = n.particles.gamma.cost.optimization,
        max_number_of_iterations  = n.iterations.gamma.cost.optimization,
        max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
        parameter_bounds          = parameter.bounds.gamma.cost,
        tryCall                   = TRUE, 
        verbose                   = FALSE,
        lhc_init                  = FALSE, 
        wait_complete_iteration   = TRUE,
        logfile                   = NULL,
        projectfile               = NULL,
        break_file                = "stopPPSO.txt"
        )
        
        # employ PSO proposed gamma and cost parameters
        tuned.svm <- svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                         y           = real.groups.training,
                         gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                         cost        = as.numeric(result.internal.gamma.cost$par[2]),
                         tunecontrol = tune.control(cross = number.cross.splits),
                         probability = TRUE
        )
        
        # extract best model
        tuned.svm.model <- tuned.svm
      }
      
      # prepare counts for sample prediction
      normalized.counts.prediction <- normalized.counts[dgeParticle$biomarker.transcripts, evaluation.samples]
      
      # store data
      dgeParticle$tuned.svm.model <- tuned.svm.model
      save(dgeParticle, 
           file = paste("outputIterations/", iter, "-dgeParticle.RData", sep = ""))
      
      # perform classification of evaluation and validation series and store output
      thrombo.algo.classify.evaluation.iterations <- thrombo.algo.classify.evaluation.set(dge = dgeIncludedSamples,
                                                                                          iterations = TRUE,
                                                                                          n_iter = iter,
                                                                                          R.snapshot = "Pre-PSO-snapshot.RData"
      )
      save(thrombo.algo.classify.evaluation.iterations, 
           file = paste("outputIterations/", iter, "-results.evaluation.iterations.RData", sep = ""))
      thrombo.algo.classify.validation.iterations <- thrombo.algo.classify.validation.set(dge = dgeIncludedSamples,
                                                                                          iterations = TRUE,
                                                                                          n_iter = iter,
                                                                                          replace.counts.validation = replace.counts.validation,
                                                                                          R.snapshot = "Pre-PSO-snapshot.RData"
      )
      save(thrombo.algo.classify.validation.iterations, 
           file = paste("outputIterations/", iter, "-results.validation.iterations.RData", sep = ""))
      
      results <- list()
      results[["i"]] <- iter
      results[["evaluation.acc"]] <- thrombo.algo.classify.evaluation.iterations$roc.optimal.accuracy
      results[["evaluation.AUC"]] <- thrombo.algo.classify.evaluation.iterations$AUCorDiagonal
      results[["validation.acc"]] <- thrombo.algo.classify.validation.iterations$roc.optimal.accuracy
      results[["validation.AUC"]] <- thrombo.algo.classify.validation.iterations$AUCorDiagonal
      results
    }
    
    # summarize
    output.iterations <- data.frame(
      i = unlist(lapply(iterations.loop, function(x){x[["i"]]})),
      evaluation.acc = unlist(lapply(iterations.loop, function(x){x[["evaluation.acc"]]})),
      evaluation.AUC = unlist(lapply(iterations.loop, function(x){x[["evaluation.AUC"]]})),
      validation.acc = unlist(lapply(iterations.loop, function(x){x[["validation.acc"]]})),
      validation.AUC = unlist(lapply(iterations.loop, function(x){x[["validation.AUC"]]}))
    )
  }
  
  # summarize data
  # also check whether output files from the classifications are available
  matrix <- matrix(nrow = 6, ncol = 4)
  rownames(matrix) <- c("Training", "Evaluation", "Validation", "NA", "NA", "NA")
  colnames(matrix) <- c("n", "AUC", "95%-CI", "Accuracy")
  if (file.exists("results.classification.training.RData") == TRUE){
    load("results.classification.training.RData")
    matrix[1, ] <- c(results.classification.training$number.samples,
                     results.classification.training$AUCorDiagonal,
                     paste(round(results.classification.training$roc.95ci$ci[1], digits = 2), "-",
                           round(results.classification.training$roc.95ci$ci[3], digits = 2), sep = ""),
                     results.classification.training$roc.optimal.accuracy)
  }
  
  if (file.exists("results.classification.evaluation.RData") == TRUE){
    load("results.classification.evaluation.RData")
    matrix[2,] <- c(length(results.classification.evaluation$samples.for.evaluation),
                    results.classification.evaluation$AUCorDiagonal,
                    paste(round(results.classification.evaluation$roc.95ci$ci[1], digits = 2), "-",
                          round(results.classification.evaluation$roc.95ci$ci[3], digits = 2), sep = ""),
                    results.classification.evaluation$roc.optimal.accuracy)
  }
  
  if (file.exists("results.classification.validation.RData") == TRUE){
    load("results.classification.validation.RData")
    matrix[3,] <- c(length(results.classification.validation$samples.for.evaluation),
                    results.classification.validation$AUCorDiagonal,
                    paste(round(results.classification.validation$ci.roc$ci[1], digits = 2), "-",
                          round(results.classification.validation$ci.roc$ci[3], digits = 2), sep = ""),
                    results.classification.validation$roc.optimal.accuracy)
  }
  if (thromboSeqPSO.shuffled == TRUE){
    matrix[5, 1] <- paste("Shuffled class labels (n=", nrow(output.shuffled), "): Median AUC Validation series: ",
                         round(median(output.shuffled$validation.AUC), digits = 2), ", and IQR: ",
                         round(IQR(output.shuffled$validation.AUC), digits = 2), sep = "")
    
  }
  if (thromboSeqPSO.iterations == TRUE){
    matrix[6, 1] <- paste("Iterations (n=", nrow(output.iterations), "): Median AUC Validation series: ",
                         round(median(output.iterations$validation.AUC), digits = 2), ", and IQR: ",
                         round(IQR(output.iterations$validation.AUC), digits = 2), sep = "")
    
  }
  write.csv(matrix, file = "ROCcurve_Metrics.csv")
  return(matrix)
}

perform.RUVg.correction.validation <- function(dge = dge,
                                               output.particle = dgeParticle){
  # Performs the RUVSeq confounding variable correction specifically in a validation setting, i.e.
  # with known stable transcripts, confounding axes, and correction factors.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   output.particle: DGEList compiled during the PSO process with the setting specific
  #                    stable transcripts, ruv-axes, and correction factors.
  #
  # Returns:
  #   DGEList including the corrected raw read counts.
  
  # remove the factors that were identified as potential confounding variables from the dataset
  # collect all output from the specific particle for correction of the counts of the to-be classified sample series
  axis.group <- output.particle$axis.group
  axis.na <- output.particle$axis.na
  axis.confounding <- output.particle$axis.confounding
  axis.all <- output.particle$axis.all
  axis.drop <- 0
  axis.removed <- FALSE
  k.variables <- length(output.particle$axis.all)
  # loop the included k.variables and correct if necessary
  tmp.loop <- foreach(i = 1 : k.variables) %do% {
    if (i == 1) {
      if (i %in% axis.confounding) {
        RUVg.post.correction <- RUVg(as.matrix(dge$counts), output.particle$stable.transcripts, k = axis.drop + 1, drop = 0)
        axis.removed <- TRUE
      } else {
        axis.drop <- 1
      }
    } else {
      if (i %in% axis.confounding) {
        if(axis.removed == TRUE) {
          RUVg.post.correction <- RUVg(as.matrix(RUVg.post.correction$normalizedCounts), output.particle$stable.transcripts, k = axis.drop + 1, drop = axis.drop)
        } else {
          RUVg.post.correction <- RUVg(as.matrix(dge$counts), output.particle$stable.transcripts, k = axis.drop + 1, drop = axis.drop)
          axis.removed <- TRUE
        }
      } else {
        axis.drop <- axis.drop + 1
      }
    }
  }
  
  # prepare a new corrected countmatrix, update the total library size
  dge$raw.counts <- dge$counts
  if (length(axis.confounding) > 0){
    # if RUVg correction has been performed, add the updated countmatrix to the dge-object
    dge$ruv.counts <- as.data.frame(RUVg.post.correction$normalizedCounts)
  }
  dge$samples$ruv.lib.size <- colSums(dge$ruv.counts)
  
  # return corrected DGEList
  return(dge)
}

thrombo.algo.classify.training.set <- function(dge = dgeIncludedSamples,
                                               n.particles.gamma.cost.optimization = 50,
                                               n.iterations.gamma.cost.optimization = 4,
                                               best.particle = best.particle.input,
                                               R.snapshot = "Pre-PSO-snapshot.RData",
                                               number.cores = number.cores,
                                               verbose = TRUE){
  # Performs leave-one-out cross-validation (LOOCV) analysis of the training samples series.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   best.particle: Vector with the merged best.selection values.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #   
  # Returns:
  #  Result-container with classification details and metrics.
  
  # load R snapshot data
  load(R.snapshot)
  load(paste("outputPSO/", best.particle, ".RData", sep = ""))
  
  suppressMessages(library(foreach))
  suppressMessages(library(doMC))
  registerDoMC(cores = number.cores)
  loocv.loop <- foreach(leave.out.sample = training.samples) %dopar% {
    # assign samples
    training.samples.loocv <- training.samples[leave.out.sample != training.samples]
    real.groups.training <- dge$samples[training.samples.loocv, "group"]
    real.groups.prediction <- dge$samples[leave.out.sample, "group"]
    # load particle data
    load(paste("outputPSO/", best.particle, ".RData", sep = ""))
   
    dgeParticle <- dgeTraining
    # perform RUVg correction
    dgeTraining <- perform.RUVg.correction.validation(dge = dge, 
                                                      output.particle = dgeParticle)
    dgeTraining$counts <- dgeTraining$ruv.counts
    dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
    
    # TMM-normalization
    dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                             normalize.on.training.series = TRUE, 
                                             samples.for.training = training.samples.loocv,
                                             ref.sample.readout = FALSE) # calculate normalization factors
    dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                             normalize.on.training.series = TRUE, 
                                             samples.for.training = training.samples.loocv,
                                             ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
    dgeTraining$samples <- droplevels(dgeTraining$samples)
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
    
    # train SVM-model, with grid search for gamma and cost.
    tuned.svm <- tune.svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples.loocv]),
                          y           = real.groups.training,
                          gamma       = svm.gamma.range,
                          cost        = svm.cost.range,
                          tunecontrol = tune.control(cross = number.cross.splits),
                          probability = TRUE
    )
    # extract best model
    tuned.svm.model <- tuned.svm[["best.model"]]
    
    if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
      # if necessary - broaden the range of cost and gamma
      if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
        svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
      }
      if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
        svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
      }
      svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                           svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
      svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                          svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])
      
      # input parameters:
      parameter.bounds.gamma.cost = matrix(ncol = 2, nrow = 2)
      rownames(parameter.bounds.gamma.cost) <- c("gamma","cost")
      parameter.bounds.gamma.cost[,1] <- c(svm.gamma.range[1],svm.cost.range[1])
      parameter.bounds.gamma.cost[,2] <- c(svm.gamma.range[2],svm.cost.range[2])
      
      # PSO
      set.seed(2000)
      selected.transcripts <- dgeParticle$biomarker.transcripts
      result.internal.gamma.cost <- optim_pso(objective_function = if (nlevels(dgeTraining$samples$group) == 2){
        thrombo.svm.gamma.cost
      } else {
        thrombo.svm.gamma.cost.multiclass
      },
      number_of_parameters      = nrow(parameter.bounds.gamma.cost),
      plot_progress             = FALSE,
      number_of_particles       = n.particles.gamma.cost.optimization,
      max_number_of_iterations  = n.iterations.gamma.cost.optimization,
      max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
      parameter_bounds          = parameter.bounds.gamma.cost,
      tryCall                   = TRUE,
      verbose                   = FALSE,
      lhc_init                  = FALSE,
      wait_complete_iteration   = TRUE,
      logfile                   = NULL,
      projectfile               = NULL,
      break_file                = "stopPPSO.txt"
      )
      
      # employ PSO proposed gamma and cost parameters
      tuned.svm <- svm(x           = t(normalized.counts[selected.transcripts, training.samples.loocv]),
                       y           = real.groups.training,
                       gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                       cost        = as.numeric(result.internal.gamma.cost$par[2]),
                       tunecontrol = tune.control(cross = number.cross.splits),
                       probability = TRUE
      )
      
      # extract best model
      tuned.svm.model <- tuned.svm
    }
    
    # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
    normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, leave.out.sample]
    # perform prediction
    prediction.class <- predict(tuned.svm.model,
                                newdata = t(normalized.counts.prediction), 
                                probability = TRUE)
    
    # summarize results
    result <- list()
    tryName <- paste("Try", leave.out.sample, sep = "")
    result[["training.samples"]] <- training.samples.loocv
    result[["prediction.sample"]] <- leave.out.sample
    result[["predicted.group"]] <- as.character(prediction.class)
    result[["real.group"]] <- as.character(real.groups.prediction)
    result[["biomarker.panel"]] <- length(dgeParticle$biomarker.transcripts)
    result[["predictive.strength"]] <- attributes(prediction.class)$probabilities
    result
  }
  
  # load additional packages
  suppressMessages(library(reshape, warn.conflicts = F, quietly = T))
  suppressMessages(library(ROCR, warn.conflicts = F, quietly = T))
  suppressMessages(library(pROC, warn.conflicts = F, quietly = T))
  
  # summarize data into a data frame
  svm.summary <- data.frame(
    sample.ID = unlist(lapply(loocv.loop, function(x){x[["prediction.sample"]]})),
    predicted.group = unlist(lapply(loocv.loop, function(x){x[["predicted.group"]]})),
    real.group = unlist(lapply(loocv.loop, function(x){x[["real.group"]]})),
    biomarker.panel = unlist(lapply(loocv.loop, function(x){x[["biomarker.panel"]]}))
  )
  # add predictive values to the data frame
  svm.summary <- suppressWarnings(cbind(svm.summary, do.call("rbind", lapply(loocv.loop, function(x){x[["predictive.strength"]]}))))
  # warning suppressed, indicates: duplicated rownames, hence not used.
  
  # summarize in confusion matrix
  confusion.matrix <- cast(svm.summary, 
                           predicted.group ~ real.group, 
                           length, 
                           value = "sample.ID")
  confusion.matrix.evaluated <- classAgreement(as.matrix(confusion.matrix), match.names = T)
  
  if (nlevels(dgeTraining$samples$group) == 2){
    # in case two classes are included, generate ROC curve
    rocra <- prediction(as.numeric(as.character(svm.summary[, 6])), 
                        svm.summary[, 3], 
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Training Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, ncol(svm.summary)], 
                    ci = TRUE
    )
  } else {
    perfa <- NA
    AUC <- confusion.matrix.evaluated$diag
    roc.optimal.accuracy <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.summary <- NA
  }
  # summarize data
  result <- list()
  result[["number.samples"]] <- length(training.samples)
  result[["svm.summary"]] <- svm.summary
  result[["AUCorDiagonal"]] <- AUC
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["roc.95ci"]] <- roc.95ci
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  return(result)
}

thrombo.algo.classify.evaluation.set <- function(dge = dgeIncludedSamples, 
                                                 best.particle = best.particle.input, 
                                                 iterations = FALSE,
                                                 n_iter = NULL,
                                                 R.snapshot = "Pre-PSO-snapshot.RData",
                                                 verbose = TRUE){
  # Performs classification of samples included in the evaluation series in a PSO-optimized
  # classification algorithm.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   best.particle: Vector with the merged best.selection values.
  #   iterations: Whether or not (TRUE/FALSE) this function is called in the iterations control experiments.
  #   n_iter: Numeric value with the number of the iteration (only when iterations == TRUE).
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #  Result-container with classification details and metrics.
  
  # load R snapshot data
  load(R.snapshot)
  
  # load data according to whether the iterations are enabled or not
  if (iterations == TRUE){
    load(paste("outputIterations/", n_iter, "-Pre-PSO-like.RData", sep = ""))
    load(paste("outputIterations/", n_iter, "-dgeParticle.RData", sep=""))
    
    training.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-trainingSamples_subsampling.csv", sep = ""))[, 2])
    evaluation.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-evaluationSamples_subsampling.csv", sep = ""))[, 2])
  } else {
    # load particle-specific output file
    load(paste("outputPSO/", best.particle, ".RData", sep = ""))
    dgeParticle <- dgeTraining
  }
  
  # assign samples to training and evaluation
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.evaluation <- dge$samples[evaluation.samples, "group"]
  # perform RUV correction
  dgeTraining <- perform.RUVg.correction.validation(dge,
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
  # normalize using TMM normalization
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, evaluation.samples]
  # perform classification
  prediction.class <- predict(dgeParticle$tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  confusion.matrix <- table(prediction.class, real.groups.evaluation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.evaluation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
    # ROC
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3], 
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Evaluation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, ncol(svm.summary)], 
                    ci = TRUE
    )
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
  }
  
  # summarize data
  result <- list()
  result[["samples.for.training"]] <- training.samples
  result[["samples.for.evaluation"]] <- evaluation.samples
  result[["biomarker.panel.size"]] <- length(dgeParticle$biomarker.transcripts)
  result[["ruv.confounding.axes"]] <- dgeParticle$axis.confounding
  result[["svm.summary"]] <- svm.summary
  result[["confusion.matrix"]] <- confusion.matrix
  result[["confusion.matrix.evaluated"]] <- confusion.matrix.evaluated
  result[["AUCorDiagonal"]] <- AUC
  result[["roc.95ci"]] <- roc.95ci
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  return(result)
}

thrombo.algo.classify.validation.set <- function(dge = dgeIncludedSamples, 
                                                 best.particle = best.particle.input, 
                                                 iterations = FALSE,
                                                 n_iter = NULL,
                                                 replace.counts.validation = 0,
                                                 R.snapshot = "Pre-PSO-snapshot.RData",
                                                 verbose = TRUE){
  # Performs classification of samples included in the validation series in a PSO-optimized
  # classification algorithm.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   best.particle: Vector with the merged best.selection values.
  #   iterations: Whether or not (TRUE/FALSE) this function is called in the iterations control experiments.
  #   iter: Numeric value with the number of the iteration (only when iterations == TRUE).
  #   replace.counts.validation: Numeric-value indicating the number of reads counts (0 - value) 
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Result-container with classification details and metrics.
  
  # load R snapshot data
  load(R.snapshot)
  # load data according to whether the iterations are enabled or not
  if (iterations == TRUE){
    load(paste("outputIterations/", n_iter, "-Pre-PSO-like.RData", sep = ""))
    load(paste("outputIterations/", n_iter, "-dgeParticle.RData", sep = ""))
    
    training.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-trainingSamples_subsampling.csv", sep = ""))[, 2])
    evaluation.samples <- as.character(read.csv(
      paste("outputIterations/", n_iter, "-evaluationSamples_subsampling.csv", sep = ""))[, 2])
  } else {
    # load particle-specific output file
    load(paste("outputPSO/", best.particle, ".RData", sep = ""))
    dgeParticle <- dgeTraining
  }
  
  # assign samples to training and evaluation
  validation.samples <- colnames(dge)[!colnames(dge) %in% c(training.samples, evaluation.samples)]
  
  # narrow the dge to those samples relevant
  dge <- dge[, c(training.samples, validation.samples)]
  dge$samples <- droplevels(dge$samples)
  
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.validation <- dge$samples[validation.samples, "group"]
  
  # perform RUV correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge, 
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  
  # enable to replace counts with 0 to provided counts in the validation series
  # by the median of those in the training series.
  if (!replace.counts.validation %in% c("NaN", NaN)){ # omit this function when NaN is inputted
    for (assess.sample in validation.samples) { # for each sample in the validation series
      tmpA <- matrix(dgeTraining$counts[, assess.sample]) # select the counts
      sel <- which(tmpA %in% seq(0, replace.counts.validation)) # identify which counts have too little raw reads detected
      if (length(sel) > 1){ # if more than one selected, calculate median read counts of these genes in the training series and replace
        tmpB <- round(apply(dgeTraining$counts[sel, training.samples], 1, median))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, replace.counts.validation)), assess.sample] <- tmpB
      } else if (length(sel) == 1) { # if one selected, calculate median read counts of this gene in the training series and replace
        tmpB <- round(median(as.numeric(dgeTraining$counts[sel, ])))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, replace.counts.validation)), assess.sample] <- tmpB
      }
    }
  }
  
  # calculate newest total read counts (lib.size)
  dgeTraining$samples$lib.size <- colSums(dgeTraining$counts)
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, validation.samples]
  # perform classification
  prediction.class <- predict(dgeParticle$tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  confusion.matrix <- table(prediction.class, real.groups.validation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.validation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
    # ROC
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3], 
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Validation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, ncol(svm.summary)], 
                    ci = TRUE
    )
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
  }
  
  # summarize data
  result <- list()
  result[["samples.for.training"]] <- training.samples
  result[["samples.for.validation"]] <- validation.samples
  result[["biomarker.panel.size"]] <- length(dgeParticle$biomarker.transcripts)
  result[["ruv.confounding.axes"]] <- dgeParticle$axis.confounding
  result[["svm.summary"]] <- svm.summary
  result[["confusion.matrix"]] <- confusion.matrix
  result[["confusion.matrix.evaluated"]] <- confusion.matrix.evaluated
  result[["AUCorDiagonal"]] <- AUC
  result[["ci.roc"]] <- roc.95ci
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  return(result)
}


thrombo.algo.classify.evaluation.shuffled <- function(dge = dgeIncludedSamples, 
                                                      best.particle = best.particle.input, 
                                                      svm.gamma.range = 2^(-20:0),
                                                      svm.cost.range = 2^(0:20),
                                                      number.cross.splits = 2,
                                                      n.particles.gamma.cost.optimization = 50,
                                                      n.iterations.gamma.cost.optimization = 4,
                                                      shuffle.iteration = i,
                                                      R.snapshot = "Pre-PSO-snapshot.RData",
                                                      number.cores = number.cores,
                                                      verbose = TRUE){
  # Performs development of a classification algorithm based on the biomarker panel in the best PSO 
  # particle but with shuffled group labels of the training series.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   input: Vector with the merged best.selection values.
  #   svm.gamma.range: Numeric value for the range of the grid search for the best gamma parameter in SVM.
  #   svm.cost.range: Numeric value for the range of the grid search for the best cost parameter in SVM.
  #   number.cross.splits: Numeric value with the number of subseriess employed by SVM algorithm for internal tuning.
  #   n.particles.gamma.cost.optimization: Numeric-value with number of PSO particles to be employed for gamma/cost optimization.
  #   n.iterations.gamma.cost.optimization: Numeric-value with number of PSO iterations to be employed for gamma/cost optimization.
  #   shuffle.iteration: Numeric value to be provided by thromboSeqPSO.controls-function.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   Result-container with classification details and metrics.
  
  # load R environment
  load(R.snapshot)
  # reset seed for random selection of training-evaluation series
  rm(list=".Random.seed", envir = globalenv())
  
  # load particle data
  load(paste("outputPSO/", best.particle, ".RData", sep = ""))
  dgeParticle <- dgeTraining
  
  # collect group conditions
  real.groups.training <- dge$samples[training.samples, "group"]
  real.groups.evaluation <- dge$samples[evaluation.samples, "group"]
  # perform RUVg correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge, 
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
  
  # TMM-normalization
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dgeTraining, log = T, normalized.lib.sizes = T) 
  
  ### 
  ### shuffle here the group labels of the samples in the training series
  ### of note, the set.seed-function is not employed here to ensure each classification is random
  real.groups.training <- sample(dge$samples[training.samples, "group"])
  
  # first SVM-model, with grid search for gamma and cost
  tuned.svm <- tune.svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                        y           = real.groups.training,
                        gamma       = svm.gamma.range,
                        cost        = svm.cost.range,
                        tunecontrol = tune.control(cross = number.cross.splits),
                        probability = TRUE
  )
  # extract best model
  tuned.svm.model <- tuned.svm[["best.model"]]
  
  ## optimize gamma and cost by a second PSO-optimization algorithm
  # employ the training series for SVM algorithm training, and the evaluation series for optimization
  # make sure the values for cost and gamma are not infinite
  if (!do(tuned.svm.model$cost) == Inf & !do(tuned.svm.model$gamma) == Inf){
    # if necessary - broaden the range of cost and gamma
    if (tuned.svm.model$gamma == svm.gamma.range[length(svm.gamma.range)] | tuned.svm.model$gamma == svm.gamma.range[1]){
      svm.gamma.range <- 2 ^ ((log2(svm.gamma.range)[1] - 1) : svm.gamma.range[length(svm.gamma.range)])
    }
    if (tuned.svm.model$cost == svm.cost.range[length(svm.cost.range)] | tuned.svm.model$cost == svm.cost.range[1]){
      svm.cost.range <- 2 ^ ((log2(svm.cost.range)[1] - 1) : (log2(svm.cost.range)[length(svm.cost.range)] + 1))
    }
    svm.gamma.range <- c(svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) - 1],
                         svm.gamma.range[which(svm.gamma.range == tuned.svm.model$gamma) + 1])
    svm.cost.range <- c(svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) - 1],
                        svm.cost.range[which(svm.cost.range == tuned.svm.model$cost) + 1])
    
    # input parameters:
    parameter.bounds.gamma.cost = matrix(ncol = 2, nrow = 2)
    rownames(parameter.bounds.gamma.cost) <- c("gamma","cost")
    parameter.bounds.gamma.cost[,1] <- c(svm.gamma.range[1],svm.cost.range[1])
    parameter.bounds.gamma.cost[,2] <- c(svm.gamma.range[2],svm.cost.range[2])
    
    # PSO
    set.seed(2000)
    selected.transcripts <- dgeParticle$biomarker.transcripts
    result.internal.gamma.cost <- optim_pso(objective_function = if (nlevels(dgeTraining$samples$group) == 2){
                                                                    thrombo.svm.gamma.cost
                                                                  } else {
                                                                    thrombo.svm.gamma.cost.multiclass
                                                                  },
    number_of_parameters      = nrow(parameter.bounds.gamma.cost),
    plot_progress             = FALSE,
    number_of_particles       = n.particles.gamma.cost.optimization,
    max_number_of_iterations  = n.iterations.gamma.cost.optimization,
    max_number_function_calls = n.particles.gamma.cost.optimization * n.iterations.gamma.cost.optimization,
    parameter_bounds          = parameter.bounds.gamma.cost,
    tryCall                   = TRUE, 
    verbose                   = FALSE,
    lhc_init                  = FALSE, 
    wait_complete_iteration   = TRUE,
    logfile                   = NULL,
    projectfile               = NULL,
    break_file                = "stopPPSO.txt"
    )
    
    # employ PSO proposed gamma and cost parameters
    tuned.svm <- svm(x           = t(normalized.counts[dgeParticle$biomarker.transcripts, training.samples]),
                     y           = real.groups.training,
                     gamma       = as.numeric(result.internal.gamma.cost$par[1]),
                     cost        = as.numeric(result.internal.gamma.cost$par[2]),
                     tunecontrol = tune.control(cross = number.cross.splits),
                     probability = TRUE
    )
    
    # extract best model
    tuned.svm.model <- tuned.svm
  }
  
  # prepare counts for sample classification
  normalized.counts.prediction <- normalized.counts[dgeParticle$biomarker.transcripts, evaluation.samples]
  
  # store support vectors
  save(tuned.svm.model, 
       file = paste("outputShuffled/", shuffle.iteration, "-Testing-SupportVectors.RData", sep = "")
  )
  # perform classification
  prediction.class <- predict(tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  # summarize in confusion matrix
  confusion.matrix <- table(prediction.class, real.groups.evaluation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.evaluation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
     # ROC
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3],
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Evaluation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
  }
  # summarize results
  result <- list()
  result[["training.samples"]] <- training.samples
  result[["evaluation.samples"]] <- evaluation.samples
  result[["AUCorDiagonal"]] <- AUC
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  
  return(result)
}

thrombo.algo.classify.validation.shuffled <- function(dge = dgeIncludedSamples,
                                                      best.particle = best.particle.input,
                                                      shuffle.iteration = i,
                                                      n.to.impute = 0, 
                                                      R.snapshot = "Pre-PSO-snapshot.RData",
                                                      verbose = TRUE){
  # Performs validation of a classification algorithm trained using shuffled group labels.
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   best.particle: Vector with the merged best.selection values.
  #   shuffle.iteration: Numeric value to be provided by thromboSeqPSO.controls-function.
  #   n.to.impute: Numeric-value indicating the number of reads counts (0 - value) 
  #   at which the validation series samples will be supplemented by counts from the training series. In
  #   case of "NaN" this step will be omitted.
  #   R.snapshot: String with RData-file name in which the R environment has been stored.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #  Result-container with classification details and metrics.
  
  # load R environment
  load(R.snapshot)
  # load particle
  load(paste("outputPSO/", best.particle, ".RData", sep = ""))
  dgeParticle <- dgeTraining
  # assign groups and select samples
  real.groups.training <- dge$samples[training.samples, "group"]
  validation.samples <- colnames(dge)[!colnames(dge) %in% c(training.samples, evaluation.samples)]
  real.groups.validation <- dge$samples[validation.samples, "group"]
  
  # RUVg correction
  dgeTraining <- perform.RUVg.correction.validation(dge = dge, 
                                                    output.particle = dgeParticle)
  dgeTraining$counts <- dgeTraining$ruv.counts
  dgeTraining$samples$lib.size <- dgeTraining$samples$ruv.lib.size
  
  # enable to replace counts with 0 to provided counts in the validation series
  # by the median of those in the training series.
  if (!n.to.impute %in% c("NaN", NaN)){ # omit this function when NaN is inputted
    for (assess.sample in validation.samples) { # for each sample in the validation series
      tmpA <- matrix(dgeTraining$counts[, assess.sample]) # select the counts
      sel <- which(tmpA %in% seq(0, n.to.impute)) # identify which counts have too little raw reads detected
      if (length(sel) > 1){ # if more than one selected, calculate median read counts of these genes in the training series and replace
        tmpB <- round(apply(dgeTraining$counts[sel, training.samples], 1, median))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, n.to.impute)), assess.sample] <- tmpB
      } else if (length(sel) == 1) { # if one selected, calculate median read counts of this gene in the training series and replace
        tmpB <- round(median(dgeTraining$ruv.counts[sel, ]))
        dgeTraining$counts[which(dgeTraining$counts[, assess.sample] %in% seq(0, n.to.impute)), assess.sample] <- tmpB
      }
    }
  }
  
  dgeTraining <- calcNormFactorsThromboseq(dgeTraining,
                                           normalize.on.training.series = TRUE, 
                                           samples.for.training = training.samples,
                                           ref.sample.readout = FALSE) # calculate normalization factors
  dgeTraining$samples <- droplevels(dgeTraining$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts.prediction <- cpm(dgeTraining, log = T, normalized.lib.sizes = T)[dgeParticle$biomarker.transcripts, validation.samples]
  # load SVM model
  load(paste("outputShuffled/", shuffle.iteration, "-Testing-SupportVectors.RData", sep = ""))
  
  # perform classification
  prediction.class <- predict(tuned.svm.model,
                              newdata = t(normalized.counts.prediction), 
                              probability = TRUE)
  confusion.matrix <- table(prediction.class, real.groups.validation)
  confusion.matrix.evaluated <- classAgreement(confusion.matrix, match.names = T)
  
  # create classification overview
  svm.summary <- data.frame(
    sampleName = attributes(prediction.class)$names,
    predicted.group = as.character((prediction.class)[1:length(prediction.class)]),
    real.group = real.groups.validation
  )
  svm.summary <- cbind(svm.summary, 
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[1])]),
                       data.frame(attributes(prediction.class)$probabilities[,
                                                                             which(colnames(attributes(prediction.class)$probabilities) == levels(prediction.class)[2])])
  )
  colnames(svm.summary)[c(4:5)] <- levels(prediction.class)
  
  # prepare output, for binary comparison calculate AUC-value, and for multiclass comparison the overall accuracy
  if (nlevels(dgeTraining$samples$group) == 2){
    # ROC 
    rocra <- prediction(as.numeric(as.character(svm.summary[, 5])), 
                        svm.summary[, 3],
                        label.ordering = levels(dgeTraining$samples$group))
    perfa <- performance(rocra, "tpr", "fpr")
    AUC <- attributes(performance(rocra, 'auc'))$y.values[[1]]
    if (verbose == TRUE){
      print(paste("AUC Validation Series: ", attributes(performance(rocra, 'auc'))$y.values[[1]], 
                  sep = ""))
    }
    roc.summary <- data.frame(
      cutOffs = unlist(attributes(rocra)$cutoffs),
      tp = unlist(attributes(rocra)$tp),
      tn = unlist(attributes(rocra)$tn),
      fp = unlist(attributes(rocra)$fp),
      fn = unlist(attributes(rocra)$fn),
      accuracy = (unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$tn)) /
        (unlist(attributes(rocra)$fp) + unlist(attributes(rocra)$tn) + unlist(attributes(rocra)$tp) + unlist(attributes(rocra)$fn)),
      xValues = unlist(attributes(perfa)$x.values),
      yValues = unlist(attributes(perfa)$y.values)
    )
    roc.optimal.accuracy <- max(roc.summary$accuracy)
    
    # calculate confidence interval
    roc.95ci <- roc(svm.summary$real,
                    svm.summary[, ncol(svm.summary)], 
                    ci = TRUE
    )
    
  } else {
    AUC <- NA
    roc.95ci <- list(NA)
    roc.95ci$ci <- NA
    roc.optimal.accuracy <- confusion.matrix.evaluated$diag
    perfa <- NA
    roc.summary <- NA
  }
  # summarize data
  result <- list()
  result[["samples.for.training"]] <- training.samples
  result[["samples.for.evaluation"]] <- validation.samples
  result[["biomarker.panel.size"]] <- length(dgeParticle$biomarker.transcripts)
  result[["ruv.confounding.axes"]] <- dgeParticle$axis.confounding
  result[["svm.summary"]] <- svm.summary
  result[["confusion.matrix"]] <- confusion.matrix
  result[["confusion.matrix.evaluated"]] <- confusion.matrix.evaluated
  result[["AUCorDiagonal"]] <- AUC
  result[["ci.roc"]] <- roc.95ci
  result[["roc.optimal.accuracy"]] <- roc.optimal.accuracy
  result[["perfa"]] <- perfa
  result[["ROC"]] <- roc.summary
  result
  
  return(result)
}

# e1071 custom svm function
# tune function implemented from e1071 package and included in SVM training (tune.svm)
# function has been adjusted to lock randomness of the internal training series selection 
# (i.e. cross split). Adjustment highlighted by 'MB: ...'.
# additional functions were included in this file to ensure proper use while analyzing
tune.svm <- function(x, y = NULL, data = NULL, degree = NULL, gamma = NULL,
                     coef0 = NULL, cost = NULL, nu = NULL, class.weights = NULL,
                     epsilon = NULL, ...) {
  call <- match.call()
  call[[1]] <- as.symbol("best.svm")
  ranges <- list(degree = degree, gamma = gamma,
                 coef0 = coef0, cost = cost, nu = nu,
                 class.weights = class.weights, epsilon = epsilon)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    ranges = NULL
  modeltmp <- if (inherits(x, "formula"))
    tune("svm", train.x = x, data = data, ranges = ranges, ...)
  else
    tune("svm", train.x = x, train.y = y, ranges = ranges, ...)
  if (!is.null(modeltmp$best.model))
    modeltmp$best.model$call <- call
  modeltmp
}

tune.control <- function(random = FALSE,
                         nrepeat = 1,
                         repeat.aggregate = mean,
                         sampling = c("cross", "fix", "bootstrap"),
                         sampling.aggregate = mean,
                         sampling.dispersion = sd,
                         cross = 10,
                         fix = 2 / 3,
                         nboot = 10,
                         boot.size = 9 / 10,
                         best.model = TRUE,
                         performances = TRUE,
                         error.fun = NULL) {
  structure(list(random = random,
                 nrepeat = nrepeat,
                 repeat.aggregate = repeat.aggregate,
                 sampling = match.arg(sampling),
                 sampling.aggregate = sampling.aggregate,
                 sampling.dispersion = sampling.dispersion,
                 cross = cross,
                 fix = fix,
                 nboot = nboot,
                 boot.size = boot.size,
                 best.model = best.model,
                 performances = performances,
                 error.fun = error.fun
  ),
  class = "tune.control"
  )
}

tune <- function(method, train.x, train.y = NULL, data = list(),
                 validation.x = NULL, validation.y = NULL,
                 ranges = NULL, predict.func = predict,
                 tunecontrol = tune.control(),
                 ...
) {
  call <- match.call()
  
  ## internal helper functions
  resp <- function(formula, data) {
    
    model.response(model.frame(formula, data))
  }
  
  classAgreement <- function (tab) {
    n <- sum(tab)
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      p0 <- sum(diag(tab[lev, lev])) / n
    } else {
      m <- min(dim(tab))
      p0 <- sum(diag(tab[1:m, 1:m])) / n
    }
    p0
  }
  
  ## parameter handling
  if (tunecontrol$sampling == "cross")
    validation.x <- validation.y <- NULL
  useFormula <- is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0))
    data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x))
  if (is.data.frame(train.y))
    train.y <- as.matrix(train.y)
  
  ## prepare training indices
  if (!is.null(validation.x)) tunecontrol$fix <- 1
  n <- nrow(if (useFormula) data else train.x)
  
  # MB: fix set.seed here
  getOption("myseed")
  foo <- function() {
    if (!is.null(seed <- getOption("myseed")))
      set.seed(seed)
    sample(n)
  }
  options(myseed = 1000)
  perm.ind <- foo()
  
  # perm.ind <- sample(n) # sampling - this is random
  if (tunecontrol$sampling == "cross") {
    if (tunecontrol$cross > n)
      stop(sQuote("cross"), " must not exceed sampling size!")
    if (tunecontrol$cross == 1)
      stop(sQuote("cross"), " must be greater than 1!")
  }
  train.ind <- if (tunecontrol$sampling == "cross")
    tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
  else if (tunecontrol$sampling == "fix")
    list(perm.ind[1:trunc(n * tunecontrol$fix)])
  else ## bootstrap
    lapply(1:tunecontrol$nboot,
           function(x) sample(n, n * tunecontrol$boot.size, replace = TRUE))
  
  ## find best model
  parameters <- if (is.null(ranges))
    data.frame(dummyparameter = 0)
  else
    expand.grid(ranges)
  p <- nrow(parameters)
  if (!is.logical(tunecontrol$random)) {
    if (tunecontrol$random < 1)
      stop("random must be a strictly positive integer")
    if (tunecontrol$random > p) tunecontrol$random <- p
    parameters <- parameters[sample(1:p, tunecontrol$random),]
  }
  model.variances <- model.errors <- c()
  
  ## - loop over all models
  for (para.set in 1:p) {
    sampling.errors <- c()
    
    ## - loop over all training samples
    for (sample in 1:length(train.ind)) {
      repeat.errors <- c()
      
      ## - repeat training `nrepeat' times
      for (reps in 1:tunecontrol$nrepeat) {
        
        ## train one model
        pars <- if (is.null(ranges))
          NULL
        else
          lapply(parameters[para.set,,drop = FALSE], unlist)
        
        
        # set more decimals instead of round at 7 digits
        options(digits = 8)
        #if(para.set==1){
        # print(list(train.x[train.ind[[sample]],],
        #            y = train.y[train.ind[[sample]]]))
        # }
        model <- if (useFormula)
          do.call(method, c(list(train.x,
                                 data = data,
                                 subset = train.ind[[sample]]),
                            pars, list(...)
          )
          )
        else
          do.call(method, c(list(train.x[train.ind[[sample]],],
                                 y = train.y[train.ind[[sample]]]),
                            pars, list(...)
          )
          )
        
        ## predict validation set
        pred <- predict.func(model,
                             if (!is.null(validation.x))
                               validation.x
                             else if (useFormula)
                               data[-train.ind[[sample]],,drop = FALSE]
                             else if (inherits(train.x, "matrix.csr"))
                               train.x[-train.ind[[sample]],]
                             else
                               train.x[-train.ind[[sample]],,drop = FALSE]
        )
        # print(pred)
        ## compute performance measure
        true.y <- if (!is.null(validation.y))
          validation.y
        else if (useFormula) {
          if (!is.null(validation.x))
            resp(train.x, validation.x)
          else
            resp(train.x, data[-train.ind[[sample]],])
        } else
          train.y[-train.ind[[sample]]]
        # print(true.y)
        # print(validation.y)
        
        if (is.null(true.y)) true.y <- rep(TRUE, length(pred))
        
        repeat.errors[reps] <- if (!is.null(tunecontrol$error.fun))
          tunecontrol$error.fun(true.y, pred)
        else if ((is.logical(true.y) || is.factor(true.y)) && (is.logical(pred) || is.factor(pred) || is.character(pred))) ## classification error
          1 - classAgreement(table(pred, true.y))
        else if (is.numeric(true.y) && is.numeric(pred)) ## mean squared error
          crossprod(pred - true.y) / length(pred)
        else
          stop("Dependent variable has wrong type!")
      }
      # print(repeat.errors[reps])
      sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
    }
    model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
    model.variances[para.set] <- tunecontrol$sampling.dispersion(sampling.errors)
  }
  
  ## return results
  best <- which.min(model.errors)
  pars <- if (is.null(ranges))
    NULL
  else
    lapply(parameters[best,,drop = FALSE], unlist)
  structure(list(best.parameters  = parameters[best,,drop = FALSE],
                 best.performance = model.errors[best],
                 method           = if (!is.character(method))
                   deparse(substitute(method)) else method,
                 nparcomb         = nrow(parameters),
                 train.ind        = train.ind,
                 sampling         = switch(tunecontrol$sampling,
                                           fix = "fixed training/validation set",
                                           bootstrap = "bootstrapping",
                                           cross = if (tunecontrol$cross == n) "leave-one-out" else
                                             paste(tunecontrol$cross,"-fold cross validation", sep="")
                 ),
                 performances     = if (tunecontrol$performances) cbind(parameters, error = model.errors, dispersion = model.variances),
                 best.model       = if (tunecontrol$best.model) {
                   modeltmp <- if (useFormula)
                     do.call(method, c(list(train.x, data = data),
                                       pars, list(...)))
                   else
                     do.call(method, c(list(x = train.x,
                                            y = train.y),
                                       pars, list(...)))
                   call[[1]] <- as.symbol("best.tune")
                   modeltmp$call <- call
                   modeltmp
                 }
  ),
  class = "tune"
  )
}

best.tune <- function(...) {
  call <- match.call()
  modeltmp <- tune(...)$best.model
  modeltmp$call <- call
  modeltmp
}

print.tune <- function(x, ...) {
  if (x$nparcomb > 1) {
    cat("\nParameter tuning of ", sQuote(x$method), ":\n\n", sep="")
    cat("- sampling method:", x$sampling,"\n\n")
    cat("- best parameters:\n")
    tmp <- x$best.parameters
    rownames(tmp) <- ""
    print(tmp)
    cat("\n- best performance:", x$best.performance, "\n")
    cat("\n")
  } else {
    cat("\nError estimation of ", sQuote(x$method), " using ", x$sampling, ": ",
        x$best.performance, "\n\n", sep="")
  }
}

summary.tune <- function(object, ...)
  structure(object, class = "summary.tune")

print.summary.tune <- function(x, ...) {
  print.tune(x)
  if (!is.null(x$performances) && (x$nparcomb > 1)) {
    cat("- Detailed performance results:\n")
    print(x$performances)
    cat("\n")
  }
}

hsv_palette <- function(h = 2/3, from = 0.7, to = 0.2, v = 1)
  function(n) hsv(h = h, s = seq(from, to, length.out = n), v = v)

plot.tune <- function(x,
                      type=c("contour","perspective"),
                      theta=60,
                      col="lightblue",
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      swapxy = FALSE,
                      transform.x = NULL,
                      transform.y = NULL,
                      transform.z = NULL,
                      color.palette = hsv_palette(),
                      nlevels = 20,
                      ...)
{
  if (is.null(x$performances))
    stop("Object does not contain detailed performance measures!")
  k <- ncol(x$performances)
  if (k > 4) stop("Cannot visualize more than 2 parameters")
  type = match.arg(type)
  
  if (is.null(main))
    main <- paste("Performance of `", x$method, "'", sep="")
  
  if (k == 3)
    plot(x$performances[,1:2], type = "b", main = main)
  else  {
    if (!is.null(transform.x))
      x$performances[,1] <- transform.x(x$performances[,1])
    if (!is.null(transform.y))
      x$performances[,2] <- transform.y(x$performances[,2])
    if (!is.null(transform.z))
      x$performances[,3] <- transform.z(x$performances[,3])
    if (swapxy)
      x$performances[,1:2] <- x$performances[,2:1]
    x <- xtabs(error~., data = x$performances[,-k])
    if (is.null(xlab)) xlab <- names(dimnames(x))[1 + swapxy]
    if (is.null(ylab)) ylab <- names(dimnames(x))[2 - swapxy]
    if (type == "perspective")
      persp(x=as.double(rownames(x)),
            y=as.double(colnames(x)),
            z=x,
            xlab=xlab,
            ylab=ylab,
            zlab="accuracy",
            theta=theta,
            col=col,
            ticktype="detailed",
            main = main,
            ...
      )
    else
      filled.contour(x=as.double(rownames(x)),
                     y=as.double(colnames(x)),
                     xlab=xlab,
                     ylab=ylab,
                     nlevels=nlevels,
                     color.palette = color.palette,
                     main = main,
                     x, ...)
  }
}

predict.svm <-
  function (object, newdata,
            decision.values = FALSE,
            probability = FALSE,
            ...,
            na.action = na.omit)
  {
    if (missing(newdata))
      return(fitted(object))
    
    if (object$tot.nSV < 1)
      stop("Model is empty!")
    
    
    if(inherits(newdata, "Matrix")) {
      loadNamespace("SparseM")
      loadNamespace("Matrix")
      newdata <- as(newdata, "matrix.csr")
    }
    if(inherits(newdata, "simple_triplet_matrix")) {
      loadNamespace("SparseM")
      ind <- order(newdata$i, newdata$j)
      newdata <- new("matrix.csr",
                     ra = newdata$v[ind],
                     ja = newdata$j[ind],
                     ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))),
                     dimension = c(newdata$nrow, newdata$ncol))
    }
    
    sparse <- inherits(newdata, "matrix.csr")
    if (object$sparse || sparse)
      loadNamespace("SparseM")
    
    act <- NULL
    if ((is.vector(newdata) && is.atomic(newdata)))
      newdata <- t(t(newdata))
    if (sparse)
      newdata <- SparseM::t(SparseM::t(newdata))
    preprocessed <- !is.null(attr(newdata, "na.action"))
    rowns <- if (!is.null(rownames(newdata)))
      rownames(newdata)
    else
      1:nrow(newdata)
    if (!object$sparse) {
      if (inherits(object, "svm.formula")) {
        if(is.null(colnames(newdata)))
          colnames(newdata) <- colnames(object$SV)
        newdata <- model.matrix(delete.response(terms(object)),
                                as.data.frame(newdata))
        newdata <- na.action(newdata)
        act <- attr(newdata, "na.action")
      } else {
        newdata <- na.action(as.matrix(newdata))
        act <- attr(newdata, "na.action")
      }
    }
    
    if (!is.null(act) && !preprocessed)
      rowns <- rowns[-act]
    
    if (any(object$scaled))
      newdata[,object$scaled] <-
      scale(newdata[,object$scaled, drop = FALSE],
            center = object$x.scale$"scaled:center",
            scale  = object$x.scale$"scaled:scale"
      )
    
    if (ncol(object$SV) != ncol(newdata))
      stop ("test data does not match model !")
    
    ret <- .C ("svmpredict",
               as.integer (decision.values),
               as.integer (probability),
               
               ## model
               as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
               as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
               as.integer (if (object$sparse) object$SV@ia else 0),
               as.integer (if (object$sparse) object$SV@ja else 0),
               as.double  (as.vector(object$coefs)),
               as.double  (object$rho),
               as.integer (object$compprob),
               as.double  (if (object$compprob) object$probA else 0),
               as.double  (if (object$compprob) object$probB else 0),
               as.integer (object$nclasses),
               as.integer (object$tot.nSV),
               as.integer (object$labels),
               as.integer (object$nSV),
               as.integer (object$sparse),
               
               ## parameter
               as.integer (object$type),
               as.integer (object$kernel),
               as.integer (object$degree),
               as.double  (object$gamma),
               as.double  (object$coef0),
               
               ## test matrix
               as.double  (if (sparse) newdata@ra else t(newdata)),
               as.integer (nrow(newdata)),
               as.integer (if (sparse) newdata@ia else 0),
               as.integer (if (sparse) newdata@ja else 0),
               as.integer (sparse),
               
               ## decision-values
               ret = double(nrow(newdata)),
               dec = double(nrow(newdata) * object$nclasses * (object$nclasses - 1) / 2),
               prob = double(nrow(newdata) * object$nclasses),
               
               PACKAGE = "e1071"
    )
    
    ret2 <- if (is.character(object$levels)) # classification: return factors
      factor (object$levels[ret$ret], levels = object$levels)
    else if (object$type == 2) # one-class-classification: return TRUE/FALSE
      ret$ret == 1
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
      ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
    else
      ret$ret
    
    names(ret2) <- rowns
    ret2 <- napredict(act, ret2)
    
    if (decision.values) {
      colns = c()
      for (i in 1:(object$nclasses - 1))
        for (j in (i + 1):object$nclasses)
          colns <- c(colns,
                     paste(object$levels[object$labels[i]],
                           "/", object$levels[object$labels[j]],
                           sep = ""))
        attr(ret2, "decision.values") <-
          napredict(act,
                    matrix(ret$dec, nrow = nrow(newdata), byrow = TRUE,
                           dimnames = list(rowns, colns)
                    )
          )
    }
    
    if (probability && object$type < 2) {
      if (!object$compprob)
        warning("SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
      else
        attr(ret2, "probabilities") <-
          napredict(act,
                    matrix(ret$prob, nrow = nrow(newdata), byrow = TRUE,
                           dimnames = list(rowns, object$levels[object$labels])
                    )
          )
    }
    
    ret2
  }

# End