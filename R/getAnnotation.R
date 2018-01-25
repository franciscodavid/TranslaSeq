################################################################################
# Wrapper to C++ function
################################################################################


getAnnotation <- function(annotation, threads = 1) {
  RcppParallel::setThreadOptions(numThreads = threads)
  ret <- gtf2SAF(annotation)
  RcppParallel::setThreadOptions()
  ret
}

