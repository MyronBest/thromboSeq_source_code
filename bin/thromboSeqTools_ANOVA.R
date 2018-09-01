# R-tools for thromboSeq
# Functions included for ANOVA analysis
# Authors       : Myron G. Best & Sjors G.J.G. In 't Veld
# Email         : m.best@vumc.nl; g.intveld1@vumc.nl
# Date          : 1st of September 2018
# Revision      : None

thromboSeqANOVA <- function(dge = dgeIncludedSamples,
                            k.variables = 3,
                            variable.to.assess = c("Age","lib.size"),
                            variable.threshold = c(0.2, 0.8),
                            ruvg.pvalue.threshold.group = 1e-2,
                            ruvg.pvalue.threshold.strongest.variable = 1e-2,
                            training.series.only = FALSE,
                            select.biomarker.FDR = FALSE,
                            plot = TRUE,
                            clinical.info.heatmap = c("group"),
                            swarm.optimization = TRUE,
                            n.particles = 200,
                            n.iterations = 12,
                            figureDir = "figureOutputFolder",
                            number.cores = 2,
                            verbose = TRUE){
  # Compute ANOVA Likelihood ratio for thromboSeq. Optimizes the heatmap clustering
  # with particle swarm optimization.
  #
  # Args:
  #   dge: DGEList with count table, sample info table and gene info table.
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
  #   select.biomarker.FDR: TRUE/FALSE whether the ANOVA output should be filtered by FDR (TRUE) 
  #                         or PValue (FALSE) statistics.
  #   plot:                 Whether or not to produce heatmaps from ANOVA results.
  #   clinical.info.heatmap: Vector containing sample info column names of clinical info to be included
  #                          in the heatmaps color-coding dendrogram.
  #   swarm.optimization: TRUE/FALSE whether the heatmap clustering should be improved by PSO.
  #   n.particles: Numeric-value with number of PSO particles per iteration.
  #   n.iterations: Numeric-value with number of PSO iterations in total.
  #   figureDir: String with directory in which figures can be outputted.
  #   number.cores: Vector indicating number of computational cores to be used.
  #   verbose: Whether to print output or not (TRUE/FALSE).
  #
  # Returns:
  #   thromboSeq.anova data frame with output of the ANOVA Likelihood analysis. Prints
  #   heatmaps and csv file with ANOVA results.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!is.numeric(number.cores)){
    stop("Provide numeric value for number.cores")
  }
  
  if (!is.numeric(n.particles)){
    stop("Provide numeric value for n.particles")
  }
  
  if (!is.numeric(n.iterations)){
    stop("Provide numeric value for n.iterations")
  }
  
  if (all(!clinical.info.heatmap %in% colnames(dge$samples))){
    stop("Provide existing column names in clinical.info.heatmap")
  }
  
  # load required packages
  if (verbose == TRUE){
    print("Load required packages doMC, foreach, and RUVSeq")
  }
  suppressMessages(library(doMC))
  suppressMessages(library(foreach))
  suppressMessages(library(RUVSeq))
  
  # perform RUV correction
  dge <- perform.RUVg.correction(dge = dge,
                                 k.variables = k.variables, 
                                 variable.to.assess = variable.to.assess,
                                 variable.threshold = variable.threshold, 
                                 ruvg.pvalue.threshold.group = ruvg.pvalue.threshold.group,
                                 ruvg.pvalue.threshold.strongest.variable = ruvg.pvalue.threshold.strongest.variable,
                                 training.series.only = training.series.only)
  # replace raw counts by RUV-corrected counts, and recalculate lib-size
  dge$counts <- dge$ruv.counts
  dge$samples$lib.size <- colSums(dge$counts)
  # perform TMM normalisation
  dge <- calcNormFactorsThromboseq(dge,
                                   normalize.on.training.series = FALSE, 
                                   ref.sample.readout = FALSE) # calculate normalization factors
  dge <- calcNormFactorsThromboseq(dge,
                                   normalize.on.training.series = FALSE, 
                                   ref.sample.readout = TRUE) # store the reference sample employed for TMM-normalization
  dge$samples <- droplevels(dge$samples)
  # calculate counts-per-million matrix (log-transformed and normalized via the TMM-normalization factor)
  normalized.counts <- cpm(dge, log = T, normalized.lib.sizes = T) 
  
  # Likelihood-ratio test modified for thromboSeq (ANOVA)
  thromboSeq.anova <- anovaLikeEdgeRthromboSeq(dge, 
                                               method = "TMM",
                                               normalize.on.training.series = FALSE, 
                                               ref.column = dge$refSample
  )
  
  if (plot == TRUE){
    write.csv(thromboSeq.anova, 
              file = "thromboSeq.anova-topTable.csv", 
              row.names = F)
    # select transcripts with (False Discovery Rate) FDR or p-value < 0.05
    if (select.biomarker.FDR == TRUE){
      selected.transcripts <- rownames(thromboSeq.anova)[
        thromboSeq.anova$FDR < 0.05 & 
          thromboSeq.anova$logCPM > 3 & 
          thromboSeq.anova$chromosome_name %in% c(1:22, "X")
        ]  
    } else {
      selected.transcripts <- rownames(thromboSeq.anova)[
        thromboSeq.anova$PValue < 0.05 & 
          thromboSeq.anova$logCPM > 3 & 
          thromboSeq.anova$chromosome_name %in% c(1:22,"X")
        ]  
    }
    
    # Prepare clinical info (group) and color coding for heatmap
    relevant.info <- as.matrix(dge$samples[, clinical.info.heatmap])
    col.info.all.samples <- matrix(ncol = length(clinical.info.heatmap),
                                   nrow = nrow(relevant.info))
    # loop the provided clinical info variables and select a categorial
    # color code bar or numeric color code bar
    for (variable in clinical.info.heatmap){
      if (nlevels(dge$samples[, variable]) == 2){
        if (is.numeric(dge$samples[, variable])){
          col.info.all.samples[, which(variable==clinical.info.heatmap)] <- 
            numericToColors(relevant.info[,which(variable==clinical.info.heatmap)])
        } else {
          col.info.all.samples[, which(variable==clinical.info.heatmap)] <- 
            categoriesToColors(relevant.info[, which(variable==clinical.info.heatmap)], colors = c("#a8ddb5", "#43a2ca"))
        }
      } else {
        if (is.numeric(dge$samples[, variable])){
          col.info.all.samples[, which(variable==clinical.info.heatmap)] <- 
            numericToColors(as.numeric(relevant.info[, which(variable==clinical.info.heatmap)]))
        } else {
          col.info.all.samples[, which(variable==clinical.info.heatmap)] <- 
            categoriesToColors(relevant.info[, which(variable==clinical.info.heatmap)])
        }
      }
    }     
    
    # perform clustering and store heatmap in pdf
    pdf(paste(figureDir, "/Heatmap-ANOVAtest.pdf", sep = ""))
    h <- rna.clustering(x                = normalized.counts[selected.transcripts,], 
                        margins          = c(5,2), 
                        keysize          = 1, 
                        key              = T, 
                        density.info     = "histogram", 
                        column.colors    = col.info.all.samples, 
                        NumColSideColors = ncol(col.info.all.samples), 
                        minmaxCol        = 1.5)
    dev.off()
    
    # PSO-enhanced ANOVA heatmap
    # improve clustering significance by swarm optimization, 
    # i.e. select the best FDR or p-value threshold for dendrogram clustering
    if (swarm.optimization == TRUE){
    # calculate current p-value of sample dendrogram (clustering, Fisher's exact test)
    suppressMessages(library(dendroextras))
    heatmap.fisher <- fisher.test(table(data.frame(cbind(as.matrix(slice(
      h$colDendrogram, k = nlevels(dge$samples$group))),
      as.character(dge$samples[
        rownames(as.matrix(slice(h$colDendrogram, k = nlevels(dge$samples$group)))),
        "group"])))))
    
    # optimize the heatmap clustering using PSO. Improve based on the Fisher's exact test
    # clustering.
    pso.parameter.bounds <- matrix(ncol = 2, nrow = 1)
    rownames(pso.parameter.bounds) <- c("fdrTinput")
    pso.parameter.bounds[,1] <- c(1e-10)
    pso.parameter.bounds[,2] <- c(1.0)
    save(list = ls(envir = environment(), all.names = TRUE), 
         file = "Pre-PSO-snapshot.RData", envir = environment()) 
    
    suppressMessages(library(ppso))
    set.seed(1000) # lock random number generator, required for data reproducibility
    result <- optim_pso(objective_function        = thrombo.algo.heatmap,
                        number_of_parameters      = nrow(pso.parameter.bounds),
                        plot_progress             = FALSE,
                        number_of_particles       = n.particles, 
                        max_number_of_iterations  = n.iterations, 
                        max_number_function_calls = n.particles * n.iterations,
                        parameter_bounds          = pso.parameter.bounds,
                        tryCall                   = TRUE,
                        lhc_init                  = TRUE, 
                        wait_complete_iteration   = TRUE,
                        logfile                   = "ppso-heatmap.log", 
                        projectfile               = "ppso-heatmap.pro",
                        load_projectfile          = "no"
    )
    file.remove("Pre-PSO-snapshot.RData")
    
    # select PSO-enhanced biomarker panel for heatmap clustering
    if (select.biomarker.FDR == TRUE){
      selected.transcripts <- rownames(thromboSeq.anova)[thromboSeq.anova$FDR < as.numeric(result$par) & 
                                                           thromboSeq.anova$logCPM > 3 & 
                                                           thromboSeq.anova$chromosome_name %in% c(1:22, "X")]  
    } else {
      selected.transcripts <- rownames(thromboSeq.anova)[thromboSeq.anova$PValue < as.numeric(result$par) & 
                                                           thromboSeq.anova$logCPM > 3 & 
                                                           thromboSeq.anova$chromosome_name %in% c(1:22, "X")]  
    }
    # perform clustering and store heatmap in pdf
    pdf(paste(figureDir, "/Heatmap-ANOVAtest-PSOoptimized.pdf", sep = ""))
    h <- rna.clustering(x                = normalized.counts[selected.transcripts,], 
                        margins          = c(5,2), 
                        keysize          = 1, 
                        key              = T, 
                        density.info     = "histogram", 
                        column.colors    = col.info.all.samples, 
                        NumColSideColors = ncol(col.info.all.samples), 
                        minmaxCol        = 1.5)
    dev.off()
    
    # calculate optimized p-value of sample dendrogram (clustering, Fisher's exact test)
    heatmap.fisher <- fisher.test(table(data.frame(cbind(as.matrix(slice(
      h$colDendrogram, k = nlevels(dge$samples$group))),
      as.character(dge$samples[
        rownames(as.matrix(slice(h$colDendrogram, k = nlevels(dge$samples$group)))),
        "group"])))))
    print(paste("Fisher exact test PSO-enhanced heatmap:", round(as.numeric(heatmap.fisher$p.value), digits = 6)))
    }
  }
  
  # return ANOVA output
  return(thromboSeq.anova)
}

makeColSideColors <- function(x, colors){
  # Prepare data frame with distinct color codes for included variables printed with heatmap clustering.
  #
  # Args:
  #   x: Data frame with relevant info (rows samples, columns variables) to be color-coded.
  #   colors: Vector with colors to include.
  #
  # Returns:
  #   Data frame with color coding for heatmap clustering.
  
  require(affy)
  col <- data.frame(row.names = rownames(x))
  for (column in colnames(x)){
    if (is.numeric(x[, column])){
      col[, column] <- numericToColors(x[, column])
    } else {
      col[, column] <- categoriesToColors(x[, column], colors = colors)		
    }
  }
  return(as.matrix(col))
}

categoriesToColors <- function(x, colors, palette = "Set1"){
  # Select colors for categorial variables for data heatmap visualization.
  #
  # Args:
  #   x: Data frame with relevant info (rows samples, columns variables) to be color-coded.
  #   colors: Vector with colors to include.
  #   
  # Returns:
  #   Data frame with color coding for heatmap clustering.
  
  if (is.numeric(x)){ 
    stop("You should provide non-numeric vector")
  }
  
  # determine the amount of categories
  n <- length(unique(x))
  if(n < 3){n <- 3}
  require(RColorBrewer)
  if (missing(colors)){
    colors <- brewer.pal(n, name = palette)[factor(x)]
  } else {
    colors <- colors[factor(x)]
  }
  colors[which(is.na(colors))] <- "grey"
  return(colors)
}

numericToColors <- function(x, na.color = "grey"){
  # Select colors for numeric variables for data heatmap visualization.
  #
  # Args:
  #   x: Data frame with relevant info (rows samples, columns variables) to be color-coded.
  #   colors: Vector with colors to include.
  #
  # Returns:
  #   Data frame with color coding for heatmap clustering.
  
  require(gplots)
  if (!is.numeric(x)){ 
    stop("You should provide numeric vector")
  }
  # determine the amount of numbers	
  n <- length(x)
  colors <- rep(na.color, n)
  no.na <- which(!is.na(x))
  require(RColorBrewer)
  colors[no.na] <- colorpanel(length(no.na),
                              low="darkblue",
                              mid="white",
                              high="darkred")[rank(x[no.na])]
  return(colors)
}

thrombo.algo.heatmap <- function(x){
  # Function to calculate sample clustering in heatmap, to be optimized by PSO.
  #
  # Args:
  #   x: Input FDR or p-value provided by PSO.
  #
  # Returns:
  #   Fisher's exact test p-value of sample clustering.
  
  load("Pre-PSO-snapshot.RData")
  
  # calculate biomarker panel based on PSO-proposed threshold.
  if(select.biomarker.FDR == TRUE){
    selected.transcripts <- rownames(thromboSeq.anova)[thromboSeq.anova$FDR < x[1] &
                                                         thromboSeq.anova$logCPM > 3 &
                                                         thromboSeq.anova$chromosome_name %in% c(1:22, "X")]
  } else {
    selected.transcripts <- rownames(thromboSeq.anova)[thromboSeq.anova$PValue < x[1] &
                                                         thromboSeq.anova$logCPM > 3 &
                                                         thromboSeq.anova$chromosome_name %in% c(1:22, "X")]
  }
  
  # make sure at least 2 transcripts are employed for clustering (if 1, clustering is not possible)
  library(dendroextras)
  
  if (length(selected.transcripts) > 2){
    fisher <- fisher.test(table(data.frame(cbind(as.matrix(slice(as.dendrogram(
      clustWard(distPearson(t(normalized.counts[selected.transcripts, ])))), 
      k = nlevels(dge$samples$group))),
      as.character(dge$samples[
        rownames(as.matrix(slice(as.dendrogram(clustWard(distPearson(t(normalized.counts[selected.transcripts, ])))), 
                                 k=nlevels(dge$samples$group)))), "group"])))))
    # return Fisher's p-value
    return(fisher$p.value)
  } else {
    return(1)
  }
}

# ANOVA Likelihood ratio testing using edgeR
# custom calcNormFactorsThromboseq that includes TMM-reference sample selection from only
# the training series.
# the original coding was extracted from the edgeR library package
# adjustment has been highlighted with 'MB: ...'
anovaLikeEdgeRthromboSeq <- function(y, 
                                     model = "group", 
                                     normalization = "TMM", 
                                     method = c("TMM", "RLE", "upperquartile", "none"), 
                                     ref.column= y$refSample, 
                                     logratioTrim = .3, 
                                     sumTrim = 0.05, 
                                     doWeighting = TRUE, 
                                     Acutoff = -1e10, 
                                     p = 0.75, 
                                     BCV, 
                                     normalize.on.training.series = TRUE, 
                                     samples.for.training = NULL){
  if (missing(y)){
    stop("Provide DGElist object")
  }
  stopifnot(class(y) == "DGEList")
  # the DGEList should containt y$samples with one column
  # named groups which is the factor for the analysis
  model <- paste("~0 + ", model, sep = "")
  model <- as.formula(model)
  design <- model.matrix(model, data = y$samples)
  colnames(design)  <- gsub("^group(.*$)", "\\1", colnames(design))
  # get all pairwise combinations in vector
  combinationsMatrix <- combn(levels(y$samples$group), 2)
  comparisonsCharacter <- apply(combinationsMatrix, 2, function(x){paste(x[c(2, 1)], collapse = " - ")})
  # make all contrasts
  pairwise.contrasts <- makeContrasts(contrasts=eval(comparisonsCharacter), levels = design)
  # normalize counts and fit general model
  x <- as.matrix(y$counts)
  lib.size <- y$samples$lib.size
  # MB: in case only training set has to be employed for TMM reference sample selection
  # normalize.on.training.series is TRUE, and only those samples are eligible for TMM selection.
  if(normalize.on.training.series == TRUE){
    refGroupSamples <- colnames(y)[colnames(y) %in% samples.for.training]
  } else {
    refGroupSamples <- colnames(y)
  }
  # MB: select lib.size values of the selected samples
  lib.size.refGroup <- y$samples$lib.size[which(colnames(y) %in% refGroupSamples)]
  
  f <- switch(method,
              TMM = {
                # MB: added filter for refGroupSamples
                f75 <- .calcFactorQuantile(data = x[,refGroupSamples], lib.size = lib.size.refGroup, p=p)
                f <- rep(NA,ncol(x))
                # MB: force to select the specific ref.column
                lib.size.refSample <- y$samples$lib.size[which(rownames(y$samples)==attributes(ref.column)$names)]
                for(i in 1:ncol(x))
                  f[i] <- .calcFactorWeighted(obs=x[,i],
                                              # MB: force to select the specific ref.column
                                              ref=x[, which(colnames(x) == attributes(ref.column)$names)], 
                                              libsize.obs=NULL, 
                                              libsize.ref=NULL, 
                                              logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
                f
                
              },
              RLE = .calcFactorRLE(x)/lib.size,
              upperquartile = .calcFactorQuantile(x,lib.size,p=p),
              none = rep(1,ncol(x))
  )
  
  # MB: in case only training set has to be employed for TMM reference sample selection
  # normalize.on.training.series is TRUE, and only those samples are eligible for TMM selection.
  # calculate the TMM-normalization factor for only the samples in the training series.
  if(normalize.on.training.series == TRUE){
    f <- f / exp(mean(log(f[which(rownames(y$samples) %in% samples.for.training)])))
  } else {
    f <- f / exp(mean(log(f)))
  }
  
  #	Output
  y$samples$norm.factors <- f
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design)
  y <- estimateGLMTagwiseDisp(y,design)
  fit <- glmFit(y,design)
  
  if(!missing(BCV)){
    # plot BCV
    png(BCV,15,15,"cm",res=300)
    plotBCV(y)
    dev.off()    
  }
  
  # get results from the fit
  top.table <- topTags(glmLRT(fit,contrast=pairwise.contrasts),n=nrow(y))$table
  return(top.table)
}


.calcFactorRLE <- function (data)
{
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
{
  #	i <- apply(data<=0,1,all)
  #	if(any(i)) data <- data[!i,,drop=FALSE]
  y <- t(t(data)/lib.size)
  f <- apply(y,2,function(x) quantile(x,p=p))
  #	f/exp(mean(log(f)))
}

.calcFactorWeighted <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
  #	TMM between two libraries
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  
  if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
  if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref
  
  logR <- log2((obs/nO)/(ref/nR))			# log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2	# absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref	 # estimated asymptotic variance
  
  #	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  if(max(abs(logR)) < 1e-6) return(1)
  
  #	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  #	a fix from leonardo ivan almonacid cardenas, since rank() can return
  #	non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
  
  if(doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  else
    f <- mean(logR[keep], na.rm=TRUE)
  
  #	Results will be missing if the two libraries share no features with positive counts
  #	In this case, return unity
  if(is.na(f)) f <- 0
  2^f
}

# here below the functions were
# downloaded and implemented from https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R.
# these functions are commonly employed for heatmap visualization.
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.1,
                      cexRow = 0.05 + 0.7/log10(nr),
                      cexCol = 0.05 + 0.7/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(colnames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

# Define clustering functions
clustWard <- function(x) hclust(x, method = "ward.D")
distPearson <- function(x) as.dist(1 - cor(t(x), method = "pearson"))

# plot function together with heatmap.3
rna.clustering <- function(x,
                           scale = "row",
                           hc = clustWard,
                           dist = distPearson,
                           minmaxCol = 2,
                           palette = "RdBu",
                           col,
                           column.colors,
                           row.colors,
                           labRow = F, ...){
  require(affy)
  require(RColorBrewer)
  require(gplots)
  br <- seq(-1*minmaxCol,minmaxCol,length.out=41)
  if(missing(col)){
    col <- colorRampPalette(rev(brewer.pal(n=11,name=palette)))(length(br)-1)
  }
  if(missing(column.colors) & missing(row.colors)){
    h <- heatmap.3(
      x = x,
      col = col,
      breaks=br,
      labRow = labRow,
      scale=scale,
      hclustfun = hc,
      distfun = dist,
      ...
    )
  } else if(!missing(column.colors) & missing(row.colors)) {
    h <- heatmap.3(
      x = x,
      col = col,
      breaks=br,
      labRow = labRow,
      ColSideColors = column.colors,
      scale=scale,
      hclustfun = hc,
      distfun = dist,
      ...
    )
  } else if(missing(column.colors) & !missing(row.colors)) {
    h <- heatmap.3(
      x = x,
      col = col,
      breaks=br,
      labRow = labRow,
      RowSideColors = row.colors,
      scale=scale,
      hclustfun = hc,
      distfun = dist,
      ...
    )
  } 
  else if(!missing(column.colors) & !missing(row.colors)) { 
    h <- heatmap.3(
      x = x,
      col = col,
      breaks=br,
      labRow = labRow,
      ColSideColors = column.colors,
      RowSideColors = row.colors,
      scale=scale,
      hclustfun = hc,
      distfun = dist,
      ...
    )
  } else {
    stop("Something went wrong ...")
  }
  return(h)
}


