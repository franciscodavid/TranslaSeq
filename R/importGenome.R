################################################################################
# Importing genome into the pipeline
################################################################################

################################################################################
# Import genomic data (FASTA, GTF) files into the environment
importGenome <- function (refName, faFile, gtfFile, outdir, threads) {

  require(Rsubread)
  path <- paste(outdir, "reference", sep = "/")
  dir.create(path, showWarnings = FALSE)
  pathRef <- paste(path, refName, sep = "/")
  dests <- c(fa = paste0(pathRef, ".fa"), gtf = paste0(pathRef, ".gtf"))
  origs <- c(fa = faFile, gtf = gtfFile)
  descs <- c(fa = "genome", gtf = "annotation")
  for (i in seq_along(dests)) {
    msg1 <- paste0("Downloading ", descs[i], " file from ", origs[i])
    msg2 <- paste0("Using ", descs[i], " file in ", origs[i])
    msg3 <- paste0("Using existent ", descs[i], " file in ", dests[i])
    msg4 <- paste0("Uncompressing ", descs[i], " file.")
    if(!file.exists(dests[i])) {
      if(grepl("\\.gz$", origs[i])) {
        destGz <- paste(dests[i], "gz", sep = ".")
        if(!file.exists(destGz)) {
          if(R.utils::isUrl(origs[i])) {
            message(msg1)
            utils::download.file(origs[i], destfile = destGz)
          } else {
            message(msg2)
            file.copy(from = origs[i], to = destGz)
          }
        } else message(msg3)
        message(msg4)
        R.utils::gunzip(destGz, overwrite = FALSE)
      } else {
        if(R.utils::isUrl(origs[i])) {
          message(msg1)
          utils::download.file(origs[i], destfile = dests[i])
        } else {
          message(msg2)
          file.copy(from = origs[i], to = dests[i])
        }
      }
    } else message(msg3)
  }
  if(!file.exists(paste0(pathRef, ".reads"))) {
    message(paste0("Building reference genome for ", dests["fa"], "... "),
            appendLF = FALSE)
    sink(paste(outdir, "pipeline.log", sep="/"), append = TRUE)
    Rsubread::buildindex(basename = pathRef, reference = dests["fa"])
    sink(NULL)
    message("DONE!")
  }

  list(alnindex = pathRef, annotation = getAnnotation(dests["gtf"], threads))
}

