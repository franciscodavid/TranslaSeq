################################################################################
# TranslaSeq pipeline
################################################################################

#' Translational control assessment from ribosome footprint and total RNA
#' libraries
#'
#' This function compares gene counts from ribosome footprint libraries with
#' gene counts from total RNA libraries for the same samples across different
#' experimental conditions to assess translational efficiency changes for
#' individual transcripts and their statistical significance.
#'
#' @usage TranslaSeq(metadata, refname, fafile, gtffile, ctrlabel, condition,
#'            outdir = 'TranslaSeq.out', preprocess = F, threads = 1,
#'            Radapt = 'CTGTAGGCACCATCAAT', platform = 'Illumina', verbose = F)
#'
#' @param metadata A dataframe with samples metadata (see Details section).
#' @param refname Name given to the genome annotation used in the analysis.
#' @param fafile Filepath or URL address of the genome sequence FASTA file.
#' @param gtffile Filepath or URL address of the genome annotation GTF file.
#' @param ctrlabel Text label for the control condition in metadata dataframe.
#' @param condition Column name containing the condition variable in metadata.
#' @param outdir Path in which output results will be stored.
#' @param preprocess Boolean indicating whether FASTQ input should be
#'        preprocessed to remove adaptors.
#' @param threads Number of threads to be used in processing the input data.
#' @param Radapt Adaptor sequence to be removed from the input FASTQ files.
#' @param platform Name of the platform which generated input data (default:
#'        Illumina).
#' @param verbose A boolean used for debugging the entire pipeline.
#'
#' @return A dataframe specifying mean ribosome footprint counts, total RNA
#'         counts and translation efficiency ratio for the control samples along
#'         with log2 fold changes for the case samples and their translation
#'         efficiency p-value and adjusted p-value.
#'
#' @details In a metadata dataframe each row represents an input file. This data
#'          structure has the following mandatory columns: \cr
#'          \enumerate{
#'            \item name: the sample name.
#'            \item file: the path to the input file.
#'            \item type: the library type, either 'rna' or 'rpf'.
#'            \item comment: string describing the input file.
#'          }
#'          All other wanted data fields must be inserted between type and
#'          comment columns. To sum up, the first column in metadata dataframe
#'          must be 'name'; the second one 'file'; the third one 'type'; then
#'          any number of arbitrary columns with other data fields and the last
#'          column must be 'comment'.
#'
#'          Valid file names should contain these suffixes: \cr
#'          \itemize{
#'            \item .fastq[.gz] files providing [compressed] FASTQ data.
#'            \item .bam|.sam files providing alignment data.
#'            \item .tsv|.count files providing count values.
#'          }
#'          All files from a metadata dataframe must be of the same type.
#'
#'          From the arbitrary columns, at least one should be named as the
#'          'condition' argument. In its field values, at least one sample
#'          must has the label 'ctrlabel', which will be the control condition.
#'          All other labels different from 'ctrlabel' will be treated as case
#'          conditions to be compared against the control one.
#'
#' @seealso A working example can be found at
#'          \url{https://franciscodavid.github.io/TranslaSeq/vignette.html}
#'
#' @author Francisco D. Morón-Duran
#'
#' @export
TranslaSeq <- function (metadata, refname, fafile, gtffile, ctrlabel, condition,
                        outdir = "TranslaSeq.out", preprocess = FALSE,
                        threads = 1, Radapt = "CTGTAGGCACCATCAAT",
                        platform = "Illumina", verbose = FALSE) {

  stopifnot(threads <= parallel::detectCores(), file.exists(metadata))

  get_os <- function() {
    # Credit for this function to
    # https://www.r-bloggers.com/identifying-the-os-from-r/
    sysinf <- Sys.info()
    if (!is.null(sysinf)) {
      os <- sysinf['sysname']
      if (os == 'Darwin') os <- "osx"
    } else {
      os <- .Platform$OS.type
      if (grepl("^darwin", R.version$os)) os <- "osx"
      if (grepl("linux-gnu", R.version$os)) os <- "linux"
    }
    tolower(os)
  }

  ty <- c('character', 'character', 'factor', 'factor', 'character')
  names(ty) <- c("name", "file", "type", condition, "comment")
  samples <- utils::read.table(metadata, header = TRUE, check.names = FALSE,
                        colClasses = ty)
  rownames(samples) <- samples$name
  dir.create(outdir, showWarnings = FALSE)
  genome <- importGenome(refname, fafile, gtffile, outdir, threads)
  BAM <- paste0(outdir, "/alignments/", samples$name, ".bam")
  SAM <- paste0(outdir, "/alignments/", samples$name, ".sam")
  CNT <- paste0(outdir, "/counts/", samples$name, ".count")
  DEF <- paste(dirname(metadata), samples$file, sep = '/')
  startingFiles <- function(files, default) {
    for(f in files) if(all(file.exists(f))) return(f)
    default
  }
  samples$file <- startingFiles(list(CNT, SAM, BAM), DEF)
  file.remove(paste(outdir, "pipeline.log", sep="/"))
  experiment <- if(all(grepl("\\.fastq\\.gz$|\\.fastq$", samples$file))) {
    if(get_os() %in% c("osx", "linux")) {
      if("Rsubread" %in% installed.packages()[,'Package']) {
        samples$file <- if (preprocess) {
          message("Starting from FASTQ files.")
          f <- preprocessReads(samples, outdir, Radapt, threads, platform, verbose)
          samples$file <- f
          alignReads(samples, outdir, genome$alnindex, genome$annotation, threads)
        } else {
          message("Starting from raw FASTQ files.")
          alignReads(samples, outdir, genome$alnindex, genome$annotation, threads)
        }
        countAlignments(samples, outdir, genome$annotation, threads)
      } else {
        warning("You should install Rsubread package for this option.")
        "fnotsupp"
      }
    } else {
      warning("Sequence reads are not supported in this enviroment.")
      "fnotsupp"
    }
  } else if (all(grepl("\\.bam$|\\.sam$", samples$file))) {
    if(get_os() %in% c("osx", "linux")) {
      if("Rsubread" %in% installed.packages()[,'Package']) {
        message("Starting from alignment files.")
        countAlignments(samples, outdir, genome$annotation, threads)
      } else {
        warning("You should install Rsubread package for this option.")
        "fnotsupp"
      }
    } else {
      warning("Sequence reads are not supported in this enviroment.")
      "fnotsupp"
    }
  } else if (all(grepl("\\.tsv$||\\.count$", samples$file))) {
    message("Starting from count values.")
    expSetFromCounts(samples, genome$annotation)
  } else "fnotsupp"
  if (is.character(experiment) && experiment == "fnotsupp") {
    stop("Filetype not supported.")
  }
  runTests(experiment, ctrlabel, condition, threads)
}

