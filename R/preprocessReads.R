################################################################################
# Preprocess FastQ records. Filter, unclip, trim and size-select
################################################################################


################################################################################
# Preprocess FastQ records. Filter, unclip, trim and size-select
preprocessReads <- function(samples, outdir, Radapt, threads = 1,
                            platform = "Illumina", verbose  = FALSE) {

  fqFN <- samples$file
  stopifnot(all(grepl("\\.fastq\\.gz$|\\.fastq$", fqFN)))

  path <- paste(outdir, "fastq", sep="/")
  fastq <- paste(path, basename(fqFN), sep = "/")
  dir.create(path, showWarnings = FALSE)
  v <- verbose
  doMC::registerDoMC(threads)
  `%dopar%` <- foreach::`%dopar%`
  i <- NULL
  ret <- foreach::foreach (i=1:length(fqFN)) %dopar% {
    if (!file.exists(fastq[i])) {
      # Load FastQ file
      if (v) message(paste0("[", fqFN[i], "] Reading FASTQ file..."))
      fastqFile <- ShortRead::readFastq(fqFN[i])
      totalReads <- length(fastqFile)
      if (v) message(paste0("[", fqFN[i], "] Total reads:", totalReads))

      # Discard filtered reads
      if (v) message(paste0("[", fqFN[i], "] Discarding filtered reads..."))
      if (platform == "Illumina") {
        filter <- grep('^.* [^:]*:N:[^:]*:',
                       as.character(ShortRead::id(fastqFile)))
      } else filter <- 1:length(fastqFile)
      fastqFile <- fastqFile[filter, ]
      postFiltR <- length(fastqFile)
      if (v) message(paste0("[", fqFN[i], "] Reads passing filter:", postFiltR))

      # Cut adaptor sequence and trim first base
      if (v) message(paste0("[", fqFN[i], "] Trimming first base..."))
      w <- BiocGenerics::width(fastqFile)
      Rp <- paste0(c(Radapt, rep('N', max(w) - nchar(Radapt))), collapse = "")
      fastqFile <- Biostrings::trimLRPatterns(Rpattern = Rp, Lpattern = 'N',
                                              Rfixed = FALSE, Lfixed = FALSE,
                                              subject = fastqFile)

      # Discard unclipped reads
      if (v) message(paste0("[", fqFN[i], "] Discarding unclipped reads..."))
      fastqFile <- fastqFile[w - BiocGenerics::width(fastqFile) > 1, ]
      clippedReads <- length(fastqFile)
      if (v) message(paste0("[", fqFN[i], "] Unclipped reads:", clippedReads))

      # Keep >= 25nt reads
      if (v) message(paste0("[", fqFN[i], "] Size-selecting >= 25nt reads..."))
      fastqFile <- fastqFile[BiocGenerics::width(fastqFile) >= 25, ]
      minLenReads <- length(fastqFile)
      if (v) message(paste0("[", fqFN[i], "] +24nt reads:", minLenReads))

      # Write FASTQ output
      if (v) message(paste0("[", fqFN[i], "] Writing FASTQ output..."))
      ShortRead::writeFastq(fastqFile, fastq[i], input_format = "FASTQ",
                            compress = TRUE)
      if (v) message(paste0("[", fqFN[i], "] Done."))

      # Release memory
      rm(fastqFile)
      gc()
    }

    # Return output file's pathname
    fastq[i]
  }

  unlist(ret)
}

