###############################################################################
# Align FASTQ against a prepared genome
###############################################################################

###############################################################################
# Align FASTQ against a prepared genome
alignReads <- function(samples, outdir, alnindex, annotation, threads = 1) {

  stopifnot(all(grepl("\\.fastq\\.gz$|\\.fastq$", samples$file)))
  require(Rsubread)

  # Prepare output directory and filename
  fileOut <- paste(outdir, "alignments", sep="/")
  dir.create(fileOut, showWarnings = FALSE)
  fileOut <- paste(fileOut, samples$name, sep="/")
  fileOut <- paste(fileOut, "bam", sep='.')

  # Do not align already aligned files
  cnd <- !file.exists(paste0(fileOut, ".indel"))
  for(i in fileOut[!cnd]) message(paste0("File ", i, " already exists."))
  fileOutp <- fileOut[cnd]
  fileInp <- samples$file[cnd]

  # Perform alignment step
  message("Alignment with reference genome... ", appendLF = FALSE)
  sink(paste(outdir, "pipeline.log", sep="/"), append = TRUE)
  Rsubread::align(index = alnindex, readfile1 = fileInp,
                  output_file = fileOutp, output_format = "BAM",
                  nthreads = threads, useAnnotation = TRUE,
                  annot.ext = annotation$exon, isGTF = FALSE)
  sink(NULL)
  message("DONE!")

  # Return output's filename
  fileOut
}

