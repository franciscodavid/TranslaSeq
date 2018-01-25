################################################################################
# Get matrix of normalized counts
################################################################################

################################################################################
# Normalize by Size Factors
normalizeCounts <- function(tab, samples) {

  applySizeFactors(tab, sizeFactors(tab, samples))
}

################################################################################
# Apply given Size Factors
applySizeFactors <- function(tab, sf) {

  sweep(tab, 2, sf, "/")
}

################################################################################
# Calculate Size Factors
sizeFactors <- function(tab, samples) {

  rnaSamples <- samples[samples$type == "rna", "name"]
  rpfSamples <- samples[samples$type == "rpf", "name"]
  tab_rpf <- tab[, rpfSamples]
  tab_rna <- tab[, rnaSamples]
  lgm_rpf <- rowMeans(log(tab_rpf))
  lgm_rna <- rowMeans(log(tab_rna))
  ret <- c(apply(tab_rpf, 2, function(cnts) {
    exp(stats::median((log(cnts)-lgm_rpf)[is.finite(lgm_rpf) & cnts > 0]))
  }), apply(tab_rna, 2, function(cnts) {
    exp(stats::median((log(cnts)-lgm_rna)[is.finite(lgm_rna) & cnts > 0]))
  }))[colnames(tab)]

  ret
}

