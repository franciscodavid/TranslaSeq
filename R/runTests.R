################################################################################
# Perform differential expression analysis and detection of translational ctrl
################################################################################


################################################################################
# Wrapper to call tests
runTests <- function(expSet, ctrlLab, variable, threads = 1) {

  labs <- unique(as.character(Biobase::pData(expSet)[[variable]]))
  labs <- labs[! labs %in% ctrlLab]
  doMC::registerDoMC(threads)
  `%dopar%` <- foreach::`%dopar%`
  i <- NULL
  ret <- foreach::foreach (i = 1:length(labs)) %dopar% {
    message(paste0("Analyzing condition ", labs[i], " vs ", ctrlLab, "."))
    cond <- as.character(Biobase::pData(expSet)[[variable]]) %in% c(labs[i],
                                                                    ctrlLab)
    samples <- Biobase::sampleNames(expSet)[cond]
    ret2 <- runTest(expSet[, samples], ctrlLab, variable)
    message(paste0("Condition ", labs[i], " vs ", ctrlLab, " done."))
    ret2
  }
  names(ret) <- labs

  ret
}


################################################################################
# Wrapper to call tests
runTest <- function(expSet, ...) {

  testTC(expSet, ...)
}


################################################################################
# TC Analysis
testTC <- function(expSet, ctrlLab, variable) {

  rpf <- as.character(Biobase::pData(expSet)$type) == "rpf"
  ctr <- as.character(Biobase::pData(expSet)[[variable]]) == ctrlLab
  ctrlRpf <- Biobase::sampleNames(expSet[, rpf & ctr])
  caseRpf <- Biobase::sampleNames(expSet[, rpf & !ctr])
  ctrlRna <- Biobase::sampleNames(expSet[, !rpf & ctr])
  caseRna <- Biobase::sampleNames(expSet[, !rpf & !ctr])
  normCounts <- Biobase::assayData(expSet)$sizef
  ctrlRpfMean <- rowMeans(normCounts[, ctrlRpf, drop = FALSE])
  ctrlRnaMean <- rowMeans(normCounts[, ctrlRna, drop = FALSE])
  zeros <- ctrlRpfMean == 0 | ctrlRnaMean == 0
  ctrlRpfMean[zeros] <- ctrlRpfMean[zeros] + 1
  ctrlRnaMean[zeros] <- ctrlRnaMean[zeros] + 1
  ctrlTE <- ctrlRpfMean/ctrlRnaMean
  caseRpfMean <- rowMeans(normCounts[, caseRpf, drop = FALSE])
  caseRnaMean <- rowMeans(normCounts[, caseRna, drop = FALSE])
  zeros <- caseRpfMean == 0 | caseRnaMean == 0
  caseRpfMean[zeros] <- caseRpfMean[zeros] + 1
  caseRnaMean[zeros] <- caseRnaMean[zeros] + 1
  caseTE <- caseRpfMean/caseRnaMean
  log2FC_TE <- log2(caseTE) - log2(ctrlTE)
  log2FC_RPF <- log2(caseRpfMean) - log2(ctrlRpfMean)
  log2FC_RNA <- log2(caseRnaMean) - log2(ctrlRnaMean)
  design <- stats::as.formula(paste0('~', variable, " + type + ", variable,
                                     ":type"))
  samples <- c(ctrlRpf, ctrlRna, caseRpf, caseRna)
  eS <- expSet[, samples]
  pD <- Biobase::pData(eS)
  pD[[variable]] <- factor(as.character(pD[[variable]]))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = Biobase::exprs(eS),
                                        colData = pD, design = design)
  DESeq2::sizeFactors(dds) <- Biobase::pData(eS)$sf
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c(0,0,0,1))
  res <- res[!is.na(res$padj), ]
  res <- res[order(res$pvalue, res$padj), ]
  res$TElog2FC <- log2FC_TE[rownames(res)]
  res$caseTE <- caseTE[rownames(res)]
  res$ctrlTE <- ctrlTE[rownames(res)]
  res$caseRPF <- caseRpfMean[rownames(res)]
  res$caseRNA <- caseRnaMean[rownames(res)]
  res$ctrlRPF <- ctrlRpfMean[rownames(res)]
  res$ctrlRNA <- ctrlRnaMean[rownames(res)]
  res$RPFlog2FC <- log2FC_RPF[rownames(res)]
  res$RNAlog2FC <- log2FC_RNA[rownames(res)]

  #res[, c("ctrlRPF", "ctrlRNA", "ctrlTE", "caseRPF",
  #        "caseRNA", "caseTE", "TElog2FC", "RPFlog2FC", "RNAlog2FC",
  #        "stat", "pvalue", "padj")]
  res[, c("ctrlRPF", "RPFlog2FC", "ctrlRNA", "RNAlog2FC", "ctrlTE", "TElog2FC",
          "stat", "pvalue", "padj")]
}

