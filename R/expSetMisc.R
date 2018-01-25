################################################################################
# Functions related to the ExpressionSet with experimental values
################################################################################


################################################################################
# Build a table of data from files containing its columns. Filenames become
# column names. Rownames are conserved from the first field of each row in the
# file. Filtering can be done either by rownames or by minimum value for any
# data across the row.
expSetFromCounts <- function(samples, annotation) {

  inFile <- samples$file
  names(inFile) <- samples$name
  stopifnot(all(grepl("\\.tsv$|\\.count$", inFile)))

  tabs <- list()
  for (i in unique(inFile)) {
    tabs[[i]] <- utils::read.table(file = i, header = TRUE, check.names = FALSE,
                                   row.names = 1)
    tabs[[i]]$GeneID <- rownames(tabs[[i]])
  }
  mat <- data.frame()
  for(i in seq_along(inFile)) {
    mat <- if (i < 2) {
      mat <- data.frame(rownames(tabs[[inFile[i]]]),
                        tabs[[inFile[i]]][[names(inFile)[i]]])
      names(mat) <- c("GeneID", names(inFile)[i])
      mat
    } else {
      n <- names(mat)
      mat <- merge(mat, tabs[[inFile[i]]][, c("GeneID", names(inFile)[i])],
            all = TRUE, by = "GeneID")
      names(mat) <- c(n, names(inFile)[i])
      mat
    }
  }
  rownames(mat) <- mat$GeneID
  mat <- mat[grep("^__", mat$GeneID, invert = TRUE),
             !colnames(mat) %in% c("GeneID")]
  mat <- mat[rowSums(mat) > 0, ]

  expSet(mat, annotation, samples)
}


################################################################################
# Builds phenoData object
phenData <- function(samples, sf) {

  cN <- setdiff(colnames(samples), c("name", "file", "comment"))
  pd <- samples[, cN, drop = FALSE]
  pd$sf <- sf

  methods::new("AnnotatedDataFrame", data = pd,
      varMetadata = data.frame(labelDescription = colnames(pd),
                               row.names = colnames(pd)))
}


################################################################################
# Builds featureData object
featData <- function(annotation) {

  cn <- colnames(annotation)
  scols <- "GeneID" != cn

  methods::new("AnnotatedDataFrame",
      data = data.frame(annotation[, scols],
                        row.names = annotation[, !scols]),
      varMetadata = data.frame(labelDescription = cn[scols],
                               row.names = cn[scols]))
}


################################################################################
# Builds an ExpressionSet
expSet <- function(counts, annotation, samples) {

  genes <- annotation$gene$GeneID %in% rownames(counts)
  fD <- featData(annotation$gene[genes, ])
  counts <- as.matrix(counts[rownames(fD), ])
  sf <- sizeFactors(counts, samples)
  pD <- phenData(samples, sf)

  Biobase::ExpressionSet(assayData = Biobase::assayDataNew(exprs = counts,
                                         sizef = applySizeFactors(counts, sf)),
                phenoData = pD, featureData = fD)
}

