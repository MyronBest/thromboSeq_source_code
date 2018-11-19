# R-tools for thromboSeq
# Functions included in pre-processing of the data
# Author        : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Date          : 1st of September 2018
# Revision      : 19th of November 2018

collect.read.counts <- function(inputDir, 
                                verbose = TRUE){
  # Recognizes and loads all HTSeq files in the input directory and merges the individual 
  # HTSeq files into one read count table
  #
  # Args:
  #   inputDir: Directory in which the FASTQ-files have been mapped and output has been stored.
  #             Does search recursively.
  #   verbose:  Whether or not to show function output.
  #
  # Returns:
  #   Count table with merged HTSeq files.
  
  if (!file.exists(inputDir)){
    stop("Input directory does not exist. Provide valid input directory")
  }
  
  # list all HTSeq files recursively in the input directory
  hts <- list.files(path = inputDir, 
                    pattern = "\\.htseq.ssv$", 
                    full.names = TRUE, 
                    recursive = TRUE
                    )
  
  if (verbose == TRUE){
    print(paste("Number of HTSeq files detected: ", length(hts), sep = ""))  
  }
  
  # read individual HTSeq files
  # merge into a single count matrix
  counts <- list()
  for(sample in hts){
    counts[[sample]] <- read.table(
      file = sample, 
      header = FALSE,
      sep = "\t",
      row.names = 1
    )
    # set sample name
    name <- basename(sample)
    name <- gsub(name, pattern = "\\.htseq.ssv", replacement = "")
    colnames(counts[[sample]]) <- name
  }
  # merge read files into single count matrix
  counts <- do.call("cbind", counts)
  # remove QC output from HTSeq in the final rows of the count table
  counts <- counts[grep("^__", rownames(counts), invert = T), ]
  
  # save counts into a R-data file
  save(counts, file = "allSamples_readCounts_ISreads.RData")
  
  # return the count matrix
  return(counts)
}

prepare.dge.object <- function(sample.info, 
                               gene.info = 'bin/dgeGenesEnsembl75.RData', 
                               read.counts = counts,
                               verbose = TRUE){
  # Reads the sample info csv file as provided by the user, matches the sample IDs in the sample info 
  # with those in the count table, loads predefined gene info data, and summarizes all data into 
  # an edgeR DGE-object for downstream analyses.
  #
  # Args:
  #   sample.info: String referring to a csv-file-name that contains per sample 
  #                the patient characteristics.
  #   gene.info: String referring to the RData file that contains annotated gene info data.
  #   counts: Data frame count table as generated using the collect.read.counts-function.
  #   verbose:  Whether or not to show function output.
  #   
  # Returns:
  #   DGEList with count table, sample info table and gene info table.
  
  if (!file.exists(sample.info)){
    stop("Input csv-file does not exist.")
  }
  
  # read the csv-file containing the sample info
  sample.info.table <- data.frame(read.csv(sample.info, header = T, sep = ",")) # load the csv file, comma separated
  if (ncol(sample.info.table) == 1){
    sample.info.table <- data.frame(read.csv(sample.info, header = T, sep = ";")) # load the csv file, semicolumn separated
  }
  
  if (verbose == TRUE){
  print(paste("The sample info file contains data from ", nrow(sample.info.table), " samples", sep = ""))
  }
  
  # check if same number of samples are included in both the sample info file as count table
  if (nrow(sample.info.table) != ncol(counts)){
    stop("Mismatch in number of samples in the sample info file and count table.")
  }
  
  # sort sample.info according to colnames of read count table
  rownames(sample.info.table) <- as.character(sample.info.table$ID)
  sample.info.table <- sample.info.table[order(match(rownames(sample.info.table), colnames(counts))), ]
  
  # load gene.info file
  load(gene.info)
  
  # generate DGE-object, as implemented from the edgeR package
  library(edgeR)
  dge <- DGEList(counts = counts,
                 group = sample.info.table$Group,
                 genes = genes
  )
  dge$samples <- cbind(dge$samples, sample.info.table) # add the sample info to the object
  
  # return DGEList-object
  return(dge)
}

filter.for.platelet.transcriptome <- function(dge, 
                                              minimum.read.counts = 30,
                                              verbose = TRUE){
  # Filters the DGE-object containing the raw data for low-abundant RNAs.
  # 
  # Args:
  #   dge: DGEList outputted by the prepare.dge.object-function, contains raw count table,
  #       sample info and gene info.
  #   minimum.read.counts: Numeric-value containing the minimum number of gene counts to 
  #       be detected in at least 90% of the samples.
  #   verbose:  Whether or not to show function output.
  #
  # Returns:
  #   DGEList with filtered count table, sample info table and gene info table.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!is.numeric(minimum.read.counts)){
    stop("Provide numeric value for minimum.read.counts")
  }
  
  # filter for low-abundant transcripts, 
  # i.e. those transcripts with less than minimum.read.counts reads in more 
  # than 90% of all samples
  bad.genes <- names(which(apply(dge$counts, 1, function(x){sum(x < minimum.read.counts)}) > 0.9 * ncol(dge)))
  dgeFiltered <- dge[which(!rownames(dge) %in% bad.genes), ]
  
  if (verbose == TRUE){
  print(paste("Number of transcripts detected: ", nrow(dgeFiltered$counts), sep = ""))
  }
  
  # return DGEList-object
  return(dgeFiltered)
}

perform.RUVg.correction <- function(dge = dge,
                                    k.variables = 3,
                                    variable.to.assess = c("Age","lib.size"),
                                    variable.threshold = c(0.2, 0.8),
                                    ruvg.pvalue.threshold.group = 1e-2,
                                    ruvg.pvalue.threshold.strongest.variable = 1e-2,
                                    training.series.only = FALSE,
                                    training.series.only.samples = training.samples,
                                    number.cores = 2, 
                                    verbose = TRUE){
  # Performs the RUVSeq confounding variable correction. First, it correlates the provided potential
  # confounding factors to the raw read counts, then the stable transcripts are determined
  # according to the thresholds provided, following RUVg is executed and the potential 
  # confounding factors are correlated to the by RUVg presented axes, and finally the potential
  # confounding factors are removed from the data.
  # 
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   k.variables: Number of (expected) variables/axes in the dataset.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   ruvg.pvalue.threshold.group: Numeric-value representing the lowest p-value that should be 
  #                                reached by the correlation between the counts and any variable
  #                                in order to bypass the wanted variable 'group' to be selected.
  #   ruvg.pvalue.threshold.strongest.variable: Numeric-value representing the lowest p-value that should 
  #                                be reached by the correlation between the counts and the specific 
  #                                variable in order to this variable to be assigned to the RUVg axis.
  #   training.series.only: TRUE/FALSE whether training series should be selected only for stable gene panel
  #                         selection and RUV-correction.
  #   verbose:  Whether or not to show function output.
  #
  # Returns:
  #  DGEList including the corrected raw read counts.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
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
  
  # in case only the training series is employed for RUV, reduce allSamples here
  if (training.series.only == TRUE){
    allSamples <- training.series.only.samples
  } else {
    allSamples <- colnames(dge)
  }
  
  # add group to the variable to assess when group is not included
  if (!"group" %in% variable.to.assess){
    variable.to.assess <- c(variable.to.assess, "group")
  }
  
  # variables to assess; prepare matrix to fill with correlation coefficients during the following steps
  matrix.variables <- matrix(NA, ncol = length(variable.to.assess) + 1, nrow = nrow(dge$counts) + 1) 
  colnames(matrix.variables) <- c("geneID", variable.to.assess)
  matrix.variables[, 1] <- c(NA, rownames(dge$counts))
  for (variable in variable.to.assess){
    if (is.numeric(dge$samples[, variable])){
      matrix.variables[1, variable] <- "numeric"
    } else {
      matrix.variables[1, variable] <- "categorial"
    }
  }
  
  # Loop all provided potential confounding factors and correlate to the raw read counts per sample
  conf.loop <- data.frame(foreach(i = 1 : nrow(dge$counts), .combine = rbind) %do% {
    for (variable in variable.to.assess){
      if (matrix.variables[1, variable] == "numeric"){
        # for numeric value use cor-function
        matrix.variables[i + 1, variable ] <- cor(dge$samples[allSamples, variable], 
                                                  as.numeric(dge$counts[i, allSamples]), 
                                                  use = "complete.obs")
      } else {
        # for categorial data employ ANOVA
        matrix.variables[i + 1, variable] <- as.numeric(unlist(summary(
          aov(as.numeric(dge$counts[i, allSamples]) ~ dge$samples[allSamples, variable])))[9])
      }
    }
  })
  
  # summarize stable genes
  stable.transcripts = NA
  for (variable in variable.to.assess){
    if (variable %in% c("Group","group")){
      stable.transcripts = stable.transcripts
    } else if (variable %in% c("Age", "age")){
      stable.transcripts <- c(stable.transcripts,
                              (matrix.variables[2 : nrow(matrix.variables), 1])[
                                which((as.numeric(matrix.variables[2 : nrow(matrix.variables), variable]) > 
                                         variable.threshold[which(variable == variable.to.assess)])
                                      == TRUE)],
                              (matrix.variables[2 : nrow(matrix.variables), 1])[
                                which((as.numeric(matrix.variables[2 : nrow(matrix.variables), variable]) < 
                                         -variable.threshold[which(variable == variable.to.assess)])
                                      == TRUE)]
      )
    } else {
      stable.transcripts <- c(stable.transcripts,
                              (matrix.variables[2 : nrow(matrix.variables), 1])[
                                which((as.numeric(matrix.variables[2 : nrow(matrix.variables), variable]) > 
                                         variable.threshold[which(variable == variable.to.assess)])
                                      == TRUE)]  
      )
    }
  }
  
  # select unique stable genes and remove the first dummy 'NA'
  stable.transcripts <- unique(stable.transcripts[which(!stable.transcripts == "NA")])
  if (verbose == TRUE){
    print(paste("Number of stable transcripts selected for RUV-module correction:", 
                as.numeric(length(stable.transcripts))))
  }
  
  if (length(stable.transcripts) > 3){
  
  # perform RUVg data correction
  RUVg.pre.correction <- RUVg(x = as.matrix(dge[, allSamples]$counts), 
                              cIdx = stable.transcripts, 
                              k = k.variables)
  
  # variables to assess; prepare matrix to fill with correlation or ANOVA p-values during the following steps
  matrix.variables.ruvg <- matrix(NA, ncol = length(variable.to.assess) + 2, nrow = k.variables + 1) 
  colnames(matrix.variables.ruvg) <- c("Axis", variable.to.assess, "Verdict")
  for (variable in variable.to.assess){
    if (is.numeric(dge$samples[, variable])){
      matrix.variables.ruvg[1, variable] <- "numeric"
    } else {
      matrix.variables.ruvg[1, variable] <- "categorial"
    }
  }

  # loop through the number of k.variables and check for each variable whether it 
  # correlates to the potential confounding factors, 
  RUVg.loop <- data.frame(foreach(i = 1 : k.variables, .combine = rbind) %do% {
    matrix.variables.ruvg[i + 1, "Axis"] <- i
    for (variable in variable.to.assess){
      if (matrix.variables.ruvg[1, variable] == "numeric"){
        matrix.variables.ruvg[i + 1, variable ] <- as.numeric(cor.test(RUVg.pre.correction$W[, i], 
                                                                       dge$samples[allSamples, variable])$p.value)
      } else {
        matrix.variables.ruvg[i + 1, variable] <- min(summary(aov(
          RUVg.pre.correction$W[, i] ~ dge$samples[allSamples, variable]))[[1]][["Pr(>F)"]], na.rm = T)
      }
    }
    
    # Filter summarized results and ensure that a statistically significant correlation to the wanted biological variation
    # in the data between the groups is not erroneously removed from the dataset.
    if (as.numeric(matrix.variables.ruvg[i + 1, "group"]) > ruvg.pvalue.threshold.group) {
      strongest.variable <- colnames(matrix.variables.ruvg)[
        which(as.numeric(matrix.variables.ruvg[i + 1, 2 : ncol(matrix.variables.ruvg) - 1]) == 
                min(as.numeric(matrix.variables.ruvg[i + 1, 2 : ncol(matrix.variables.ruvg) - 1])))
        ]
    } else {
      strongest.variable <- "group"
    }
    
    # Filter for statistically significant correlation.
    if (as.numeric(matrix.variables.ruvg[i + 1, strongest.variable]) < ruvg.pvalue.threshold.strongest.variable) {
      matrix.variables.ruvg[i + 1, "Verdict"] <- strongest.variable
    } 
  })

  # remove the factors that were identified as potential confounding variables from the dataset
  axis.group <- (matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Axis"])[
    which(matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Verdict"] == "group")
    ]
  axis.na <- (matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Axis"])[
    is.na(matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Verdict"])
    ]
  axis.confounding <- (matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Axis"])[
    which(!(matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Axis"]) %in% 
            c(axis.group, axis.na))
    ]
  axis.all <- (matrix.variables.ruvg[2 : nrow(matrix.variables.ruvg), "Axis"])
  axis.drop <- 0
  axis.removed <- FALSE
  
  # loop all axes and decide whether to remove using RUVg or drop
  axis.loop <- foreach(i = 1 : k.variables) %do% {
    if (i == 1) {
      if (i %in% axis.confounding) {
        # remove 
        RUVg.post.correction <- RUVg(dge$counts, stable.transcripts, k = axis.drop + 1, drop = 0)
        axis.removed <- TRUE
      } else {
        # drop
        axis.drop <- 1
      }
    } else {
      if (i %in% axis.confounding) {
        # remove
        if(axis.removed == TRUE) {
          RUVg.post.correction <- RUVg(RUVg.post.correction$normalizedCounts, stable.transcripts, 
                                       k = axis.drop + 1, drop = axis.drop)
        } else {
          RUVg.post.correction <- RUVg(dge$counts, stable.transcripts, k = axis.drop + 1, drop = axis.drop)
          axis.removed <- TRUE
        }
      } else {
        # drop
        axis.drop <- axis.drop + 1
      }
    }
  }
  
  # prepare a new corrected DGEList countmatrix, update the total library size
  # set raw input counts aside
  dge$raw.counts <- dge$counts
  # if data has been adjusted, apply ruv.counts with corrected data
  # else replace ruv.counts with the original raw counts to prevent 
  # issues lateron in the pipeline
  if (length(axis.confounding) > 0){
    dge$ruv.counts <- as.data.frame(RUVg.post.correction$normalizedCounts)
  } else {
    dge$ruv.counts <- dge$raw.counts
  }
  dge$samples$ruv.lib.size <- colSums(dge$ruv.counts)
  
  # store output data in DGEList-object
  dge$stable.transcripts <- stable.transcripts
  dge$RUVg.loop <- RUVg.loop
  dge$axis.group <- axis.group
  dge$axis.na <- axis.na
  dge$axis.confounding <- axis.confounding
  dge$axis.all <- axis.all
  
  } else {
    # prepare a new corrected DGEList countmatrix, update the total library size
    # set raw input counts aside
    dge$raw.counts <- dge$counts
    # if data has been adjusted, apply ruv.counts with corrected data
    # else replace ruv.counts with the original raw counts to prevent 
    # issues lateron in the pipeline
    dge$ruv.counts <- dge$raw.counts
    dge$samples$ruv.lib.size <- colSums(dge$ruv.counts)
    
    dge$stable.transcripts <- NA
    dge$RUVg.loop <- NA
    dge$axis.group <- NA
    dge$axis.na <- NA
    dge$axis.confounding <- NA
    dge$axis.all <- NA
  }
  
  # return DGEList
  return(dge)
}

thromboSeqQC <- function(dge = dgeIncludedSamples, 
                         min.number.reads.per.RNAs.detected = 0, 
                         min.number.total.RNAs.detected = 750, 
                         k.variables = 3,
                         variable.to.assess = c("Age","lib.size"),
                         variable.threshold = c(0.2, 0.8),
                         ruvg.pvalue.threshold.group = 1e-2,
                         ruvg.pvalue.threshold.strongest.variable = 1e-2,
                         training.series.only = FALSE,
                         training.series.only.samples = NA,
                         leave.sample.out.threshold = 0.5, 
                         figureDir = "figureOutputFolder", 
                         number.cores = 2, 
                         verbose = FALSE){
  # Performs the thromboSeqQC steps, that includes 
  # 1) exclusion of samples with too little RNAs detected
  # 2) leave-one-sample-out cross-correlation, i.e. excluding samples with low correlation to each other 
  #
  # Args:
  #   dge: DGEList with dataset count table, sample info and gene info tables.
  #   min.number.reads.per.RNAs.detected: Minimum number of intron-spanning reads detected per transcript.
  #   min.number.total.RNAs.detected: Minimum number of total transcripts at least detected per sample.
  #   k.variables: Number of (expected) variables/axes in the dataset.
  #   variable.to.assess: Vector containing the column names of the sample info
  #                       that are potential confounding factors. Of note, for the columns
  #                       'age' or 'Age' the transcripts with a correlation coefficient below
  #                       and over the provided variable.thresholds are included (see Best et al.
  #                       Cancer Cell 2017).
  #   variable.threshold: Vector containing manually set thresholds for the potential
  #                       confounding factors in the same order as for variable.to.assess.
  #   ruvg.pvalue.threshold.group: Numeric-value representing the lowest p-value that should be 
  #                                reached by the correlation between the counts and  any variable
  #                                in order to bypass the wanted variable 'group' to be selected.
  #   ruvg.pvalue.threshold.strongest.variable: Numeric-value representing the lowest p-value that should 
  #                                be reached by the correlation between the counts and the specific 
  #                                variable in order to this variable to be assigned to the RUVg axis.
  #   training.series.only: TRUE/FALSE whether training series should be selected only for stable gene panel
  #                         selection and RUV-correction.
  #   leave.sample.out.threshold: Numeric-value that indicates lowest Pearson's correlation
  #                               coefficient for leave-on-sample-out cross-correlation.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational threads to be used in parallel.
  #   verbose: Whether or not to show function output.
  #
  # Returns:
  #  QC-filtered DGEList with count table, sample info table and gene info table.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!is.numeric(min.number.reads.per.RNAs.detected)){
    stop("Provide numeric value for min.number.reads.per.RNAs.detected")
  }
  
  if (!is.numeric(min.number.total.RNAs.detected)){
    stop("Provide numeric value for min.number.total.RNAs.detected")
  }
  
  if (!is.numeric(number.cores)){
    stop("Provide numeric value for number.cores")
  }
  
  # load required packages
  if (verbose == TRUE){
    print("Load required packages doMC, foreach, and RUVSeq")
  }
  suppressMessages(library(doMC))
  suppressMessages(library(foreach))
  suppressMessages(library(RUVSeq))
  
  if (verbose == TRUE){
    print("Determine number of RNAs detected per sample")
  }
  
  # exclude samples that have too little transcripts (intron-spanning reads) detected.
  # loop all samples individually, collect per sample the number of RNAs with more RNAs
  # detected than the provided threshold
  registerDoMC(cores = number.cores)
  detection.loop <- foreach(i = 1 : ncol(dge$counts)) %dopar% {
    # calculate number of genes with more than zero raw read counts
    # in case all genes are detected, summary would have only two output
    # options and is.na would result in FALSE, hence the if-statement
    sample.raw.counts <- dge$counts[, colnames(dge$counts)[i]]
    if (as.character(is.na(summary(as.numeric(sample.raw.counts) > 0)[3]))){
      sample.RNAs.detected <- as.numeric(summary(as.numeric(sample.raw.counts) > 
                                                   min.number.reads.per.RNAs.detected)[2])  
    } else {
      sample.RNAs.detected <- as.numeric(summary(as.numeric(sample.raw.counts) > 
                                                   min.number.reads.per.RNAs.detected)[3])  
    }
    
    # store data in container
    cont <- list()
    cont[["sample"]] <- colnames(dge$counts)[i] # collect sample ID
    cont[["sample.libsize"]] <- dge$samples[colnames(dge$counts)[i], "lib.size"] # collect lib.size
    cont[["sample.RNAs.detected"]] <- sample.RNAs.detected
    cont
  }
  
  # summarize data from loop into a data.frame
  detection.data.frame <- data.frame(
    sample = unlist(lapply(detection.loop, function(x){x[["sample"]]})),
    lib.size = unlist(lapply(detection.loop, function(x){x[["sample.libsize"]]})),
    RNAs.detected = unlist(lapply(detection.loop, function(x){x[["sample.RNAs.detected"]]}))
  )
  
  if (verbose == TRUE){
    print(paste("Median number of transcripts detected in dataset:", 
                round(median(detection.data.frame$RNAs.detected), digits = 0)))
  }
  
  # generate and store plot in output folder
  if (!file.exists(figureDir)){
    dir.create(figureDir, recursive = T)
  }
  
  upper.boundary.lib.size <- seq(from = 0, to = 10e8, by = 2e5)[
    which.min(abs(seq(from = 0, to = 10e8, by = 2e5) - max(detection.data.frame$lib.size))) + 1
    ] # automatically set best upper boundary for x-axis
  upper.boundary.RNAs.detected <- seq(from = 0, to = 10e6, by = 5e2)[
    which.min(abs(seq(from = 0, to = 10e6, by = 5e2) - max(detection.data.frame$RNAs.detected))) + 1
    ] # automatically set best upper boundary for y-axis
  
  # Dot plot of library size to RNAs detected
  # min.number.total.RNAs.detected as line included
  pdf(paste(figureDir, "/Plot-isReadsLibSize-vs-nTranscriptsDetected.pdf", sep = ""))
  plot(x =    detection.data.frame$lib.size, 
       y =    detection.data.frame$RNAs.detected,
       pch =  20, 
       col =  "#2b83ba", 
       xaxt = "n", 
       xlim = c(0, upper.boundary.lib.size),
       ylim = c(0, upper.boundary.RNAs.detected)
  )
  axis(1, at = seq(0, 9e6, by = 0.4e6), las=2)
  abline(h = min.number.total.RNAs.detected, lwd = 2)
  dev.off()
  
  # select samples that have to be excluded because of too little RNAs detected
  samples.excluded.RNAs.detected <- detection.data.frame[detection.data.frame$RNAs.detected < 
                                                           min.number.total.RNAs.detected, ]
  
  # update DGEList, and remove excluded samples
  dgeIncludedSamples <- dge[, !colnames(dge) %in% samples.excluded.RNAs.detected$sample]
  dgeIncludedSamples$samples <- droplevels(dgeIncludedSamples$samples)
  
  # perform RUV correction
  dgeIncludedSamples <- perform.RUVg.correction(dge = dgeIncludedSamples,
                                                k.variables = k.variables,
                                                variable.to.assess = variable.to.assess,
                                                variable.threshold = variable.threshold,
                                                ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                                ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                                training.series.only = training.series.only,
                                                training.series.only.samples = training.series.only.samples)
  
  # perform TMM normalization
  dgeIncludedSamples <- calcNormFactorsThromboseq(dgeIncludedSamples,
                                                  normalize.on.training.series = FALSE, 
                                                  ref.sample.readout = FALSE) # calculate normalization factors
  dgeIncludedSamples <- calcNormFactorsThromboseq(dgeIncludedSamples,
                                                  normalize.on.training.series = FALSE, 
                                                  ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
  dgeIncludedSamples$samples <- droplevels(dgeIncludedSamples$samples)
  
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dgeIncludedSamples, log = T, normalized.lib.sizes = T) 
  
  if (verbose == TRUE){
    print("Perform leave-one-sample-out cross-correlation QC")
  }
  
  # correlate all samples to all other samples in a leave-one-sample-out cross-correlation analysis
  # summarize data and filter those with Pearson's correlation below provided threshold
  all.samples <- colnames(normalized.counts)
  registerDoMC(cores = number.cores)
  tmp.loop <- data.frame(foreach(test.sample = all.samples, .combine = rbind) %dopar% {
    cor(apply(normalized.counts[, all.samples[all.samples != test.sample]], 1, median), normalized.counts[, test.sample])
  })
  leave.sample.out.output <- cbind(all.samples, tmp.loop$foreach.test.sample...all.samples...combine...rbind...dopar..)
  leave.sample.out.output <- leave.sample.out.output[order(as.numeric(leave.sample.out.output[, 2])), ]
  
  pdf(paste(figureDir,"/Plot-cross-correlation-sampleFilter.pdf", sep = ""))
  par(mai   = c(1, 1, 0.8, 1.9))
  plot(x    = as.numeric(as.character(leave.sample.out.output[, 2])), 
       pch  = 19, 
       xaxt = 'n', 
       ylim = c(0.4, 1.0)
  )
  abline(h = leave.sample.out.threshold, lwd = 2)
  dev.off()
  
  # exclude samples that did not pass QC filters
  samples.excluded.cross.corr <- leave.sample.out.output[, 1][
    which(as.numeric(leave.sample.out.output[, 2]) < leave.sample.out.threshold)]
  samples.excluded.logcpm <- rownames(data.frame(sort(apply(normalized.counts, 2, median), decreasing = T)))[
    which(as.numeric(as.character(data.frame(sort(apply(normalized.counts, 2, median), decreasing = T))$
                                    sort.apply.normalized.counts..2..median...decreasing...T.)) < 3)]
  samples.excluded.all <- c(samples.excluded.RNAs.detected$sample,
                            samples.excluded.cross.corr,
                            samples.excluded.logcpm
  )
  
  if (length(samples.excluded.all) >= 1){
    write.csv(samples.excluded.all, file = "samplesExcludedQC.csv")
  } else {
    print("No samples were excluded due to QC-failure")
  }
  
  dgeIncludedSamples <- dgeIncludedSamples[, colnames(dgeIncludedSamples)[
    which(!colnames(dgeIncludedSamples) %in% samples.excluded.all)]
    ]
  
  # plot and store boxplot with age distribution per group
  if('Age' %in% colnames(dgeIncludedSamples$samples)){
  pdf(paste(figureDir, "/BoxplotAgePerGroup.pdf", sep = ""))
  par(mai = c(1, 1, 0.8, 1.9))
  boxplot(dgeIncludedSamples$samples$Age~dgeIncludedSamples$samples$group,
          ylim   = c(0, 80),
          col    = c("#a8ddb5", "#43a2ca"),
          boxwex = 0.5
  )
  dev.off()
  }
  
  # store dataset in RData file
  if (verbose == TRUE){
    print("Save dataset")
  }
  save(dgeIncludedSamples, 
       file = paste("PSO-enhanced-dataset-", timestamp(), "-.RData", sep = ""))
  
  # return DGEList
  return(dgeIncludedSamples)
}

# calcNormFactors using edgeR
# custom calcNormFactorsThromboseq that includes TMM-reference sample selection from only
# the training series.
# the original coding was extracted from the edgeR library package
# adjustment has been highlighted with 'MB: ...'
calcNormFactorsThromboseq <- function(object, 
                                      method = c("TMM", "RLE", "upperquartile", "none"), 
                                      normalize.on.training.series = TRUE, 
                                      samples.for.training = NULL, 
                                      ref.sample.readout = FALSE, ...)
  UseMethod("calcNormFactorsThromboseq")

# function adapted from edgeR-package
# calculates the TMM-normalization factor, and selection of TMM reference sample can
# be narrowed down to samples in the training series only. Also, the function can 
# be enabled to select and store the reference sample (ref.sample.readout-option).

calcNormFactorsThromboseq.DGEList <- function(object, method = c("TMM", "RLE", "upperquartile", "none"), refColumn = NULL, 
                                              logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10, p = 0.75, 
                                              normalize.on.training.series = TRUE, samples.for.training = samples.for.training, 
                                              ref.sample.readout = ref.sample.readout, ...)
  #	Scale normalization of RNA-Seq data, for DGEList objects
  #	Created 2 October 2014.  Last modified 27 August 2015.
{
  if(ref.sample.readout == FALSE){
    object$samples$norm.factors <- calcNormFactorsThromboseq.default(object=object$counts, objectTwo=object$samples, 
                                                                     lib.size=object$samples$lib.size, method=method, 
                                                                     refColumn=refColumn, 
                                                                     logratioTrim=logratioTrim, sumTrim=sumTrim, 
                                                                     doWeighting=doWeighting, Acutoff=Acutoff, p=p,
                                                                     normalize.on.training.series = normalize.on.training.series, 
                                                                     samples.for.training = samples.for.training,
                                                                     ref.sample.readout = ref.sample.readout, ...)
    object 
  } else {
    object$refSample <- calcNormFactorsThromboseq.default(object=object$counts, objectTwo=object$samples, lib.size=object$samples$lib.size, method=method, refColumn=NULL,
                                                          logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff, p=p,
                                                          normalize.on.training.series = normalize.on.training.series, samples.for.training = samples.for.training,
                                                          ref.sample.readout = ref.sample.readout, ...) # ref.sample = object$ref.sample,
    object
  }
}

calcNormFactorsThromboseq.default <- function (object, objectTwo, lib.size = NULL, method = c("TMM", 
                                                                                              "RLE", "upperquartile", "none"), refColumn = NULL, logratioTrim = 0.3, 
                                               sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10, p = p, 
                                               normalize.on.training.series = TRUE, samples.for.training = NULL, 
                                               ref.sample.readout = ref.sample.readout, ...) 
{
  x <- as.matrix(object)
  if (any(is.na(x))) 
    stop("NA counts not permitted")
  if (is.null(lib.size)) 
    lib.size <- colSums(x)
  if (any(is.na(lib.size))) 
    stop("NA lib.sizes not permitted")
  method <- match.arg(method)
  allzero <- .rowSums(x > 0, nrow(x), ncol(x)) == 0
  if (any(allzero)) 
    x <- x[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 1) 
    method = "none"
  # MB: in case only training set has to be employed for TMM reference sample selection
  # normalize.on.training.series is TRUE, and only those samples are eligible for TMM selection.
  if (normalize.on.training.series == TRUE) {
    refGroupSamples <- rownames(objectTwo)[rownames(objectTwo) %in% samples.for.training]
  } else {
    refGroupSamples <- rownames(objectTwo)
  }
  # MB: select lib.size values of the selected samples
  lib.size.refGroup <- objectTwo$lib.size[which(rownames(objectTwo) %in% refGroupSamples)]
  # MB: in case one selection fo the TMM-reference samples is required,
  # here select TRUE, otherwise calculate normalization factors (FALSE)
  if (ref.sample.readout == FALSE) {
    f <- switch(method, TMM = {
      # MB: added filter for refGroupSamples
      f75 <- .calcFactorQuantile(data = x[,refGroupSamples], 
                                 lib.size = lib.size.refGroup, p = p)
      if (is.null(refColumn)) refColumn <- which.min(abs(f75 - mean(f75)))
      if (length(refColumn) == 0 | refColumn < 1 | refColumn > 
          ncol(x)) refColumn <- 1
      lib.size.refSample <- objectTwo$lib.size[which(rownames(objectTwo) == 
                                                       attributes(refColumn)$names)]
      f <- rep(NA, ncol(x))
      for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted(obs = x[,i], 
                                                       # MB: force to select the specific ref.column
                                                       ref = x[,which(colnames(x) == attributes(refColumn)$names)], 
                                                       libsize.obs = NULL, 
                                                       libsize.ref = NULL, 
                                                       logratioTrim = logratioTrim, 
                                                       sumTrim = sumTrim, 
                                                       doWeighting = doWeighting, 
                                                       Acutoff = Acutoff)
      f
    }, RLE = .calcFactorRLE(x)/lib.size, upperquartile = .calcFactorQuantile(x, lib.size, p = p), none = rep(1, ncol(x)))
  } else {
    f <- switch(method, TMM = {
      # MB: added filter for refGroupSamples
      f75 <- .calcFactorQuantile(data = x[, refGroupSamples], 
                                 lib.size = lib.size.refGroup, p = p)
      if (is.null(refColumn)) refColumn <- which.min(abs(f75 - mean(f75)))
      if (length(refColumn) == 0 | refColumn < 1 | refColumn > 
          ncol(x)) refColumn <- 1
      lib.size.refSample <- objectTwo$lib.size[which(rownames(objectTwo) == 
                                                       attributes(refColumn)$names)]
      f <- rep(NA, ncol(x))
      for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted(obs = x[, i],
                                                       # MB: force to select the specific ref.column
                                                       ref = x[, which(colnames(x) == attributes(refColumn)$names)], 
                                                       libsize.obs = NULL, 
                                                       libsize.ref = NULL, 
                                                       logratioTrim = logratioTrim, 
                                                       sumTrim = sumTrim, 
                                                       doWeighting = doWeighting, 
                                                       Acutoff = Acutoff)
      refColumn
    }, RLE = .calcFactorRLE(x)/lib.size, upperquartile = .calcFactorQuantile(x, lib.size, p = p), none = rep(1, ncol(x)))
  }
  # MB: in case only training set has to be employed for TMM reference sample selection
  # normalize.on.training.series is TRUE, and only those samples are eligible for TMM selection.
  # calculate the TMM-normalization factor for only the samples in the training series.
  if (ref.sample.readout == FALSE) {
    if (normalize.on.training.series == TRUE) {
      f <- f/exp(mean(log(f[which(rownames(objectTwo) %in% 
                                    samples.for.training)])))
    }
    else {
      f <- f/exp(mean(log(f)))
    }
  }
  f
}