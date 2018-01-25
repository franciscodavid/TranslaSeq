################################################################################
# Count BAM/SAM hits
################################################################################


################################################################################
# Counts a set of input BAM files, saving table of counts as an ExpressionSet
# and returning the summary statistics for the counting process.
countAlignments <- function(samples, outdir, annotation, threads = 1,
                           isGTF = FALSE) {

  stopifnot(all(grepl("\\.bam$|\\.sam$", samples$file)))
  require(Rsubread)

  fileOut <- paste(outdir, "counts", sep = '/')
  dir.create(fileOut, showWarnings = FALSE)
  fileOut <- paste(fileOut,
                   gsub("[^/]*/", "",
                        sub("\\.[bs]am$", ".count", samples$file)),
                   sep = '/')

  # Do not align already aligned files
  cnd <- !file.exists(fileOut)
  for(i in fileOut[!cnd]) message(paste0("File ", i, " already exists."))
  fileOutp <- fileOut[cnd]
  fileInp <- samples$file[cnd]

  message("Counting gene hits in alignment files... ", appendLF = FALSE)
  sink(paste(outdir, "pipeline.log", sep="/"), append = TRUE)
  fc <- Rsubread::featureCounts(fileInp, annot.ext = annotation$exon,
                                nthreads = threads, isGTFAnnotationFile = isGTF)
  sink(NULL)
  message("DONE!")
  rownames(fc$stat) <- fc$stat$Status
  fc$stat <- fc$stat[-1]
  if(!is.null(samples$name))
    colnames(fc$counts) <- colnames(fc$stat) <- samples$name
  countTab <- fc$counts[rowSums(fc$counts) > 0, , drop = FALSE]
  s <- fc$stat
  rownames(s) <- paste0("__", rownames(s))
  tab <- rbind(countTab, s)
  for (i in colnames(tab)) {
    utils::write.table(tab[, i, drop = FALSE],
                file = fileOutp[1], col.names = TRUE)
    fileOutp <- fileOutp[-1]
  }

  samples$file <- fileOut
  expSetFromCounts(samples, annotation)
}

