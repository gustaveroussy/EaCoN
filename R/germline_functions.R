EaCoN.Predict.Germline <- function(ASCATobj = NULL, prior = "G", bafbin.size = 1E+07, modelName = "E", toclustname = "BAF", mcMin = 2, mcMax = 4, nfactor = 4, BAF.filter = .9, BAF.cutter = 0, segmentLength = 5, genome = "hg19") {

  ## TEMP
  # ASCATobj = my.ascat.data
  # prior = "G"
  # bafbin.size = 1E+07
  # modelName = "E"
  # toclustname = "BAF"
  # mcMin = 2
  # mcMax = 4
  # nfactor = 4
  # segmentLength = 5
  # genome = "hg19"
  # BAF.filter <- .90
  # BAF.cutter <- .02

  Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
  rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)

  Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[, 1]), 1]
  names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[, 1])]
  Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA), 1]
  names(Tumor_LogR_noNA) <- names(Tumor_BAF_noNA)
  bafcdf <- data.frame(idx = seq_len(nrow(ASCATobj$Tumor_BAF)), ASCATobj$SNPpos, value = ASCATobj$Tumor_BAF[,1], stringsAsFactors = FALSE)
  tbsm <- bsm <- bafcdf$value
  tbsm[which(tbsm < 0)] <- -tbsm[which(tbsm < 0)]
  tbsm[which(tbsm > 1)] <- 2 - tbsm[which(tbsm > 1)]
  bsm <- ifelse(bsm < .5, bsm, 1 - bsm)
  if (toclustname == "BAF") bafcdf$value <- tbsm else if (toclustname == "mBAF") bafcdf$value <- bsm else stop("Unknown toclustname value !")
  bafcdf.nna <- bafcdf[!is.na(bafcdf$value),]

  ## BAF.cutter modif
  baf.Xin <- bafcdf.nna$value >= BAF.cutter & bafcdf.nna$value <= (1-BAF.cutter)

  ## Computing genome-wide prior (if needed)
  set.seed(123456)
  if (!is.null(prior)) {
    if (prior == "G") {
      message(tmsg("Computing whole genome BAF clustering prior ..."))
      # priorX <- mclust::defaultPrior(bafcdf.nna$value, G = mcMin:mcMax, modelName = modelName)
      priorX <- mclust::defaultPrior(bafcdf.nna$value[baf.Xin], G = mcMin:mcMax, modelName = modelName)
    }
  }

  ## Computing banded BAF clusters
  message(tmsg("BAF fragments clustering ..."))
  `%do%` <- foreach::"%do%"
  mcreslist <- foreach::foreach(k = unique(ASCATobj$SNPpos$chrs), .combine = "c") %do% {
    WESk <- bafcdf.nna[bafcdf.nna$chrs == k,]
    krange <- range(WESk$pos, na.rm = TRUE)
    binmax <- ceiling(diff(krange) / bafbin.size)

    if (!is.null(prior)) {
      if (prior == "K") {
        set.seed(123456)
        priorX <- mclust::defaultPrior(WESk$value, G = mcMin:mcMax, modelName = modelName)
      }
    }

    mcresk <- foreach::foreach(bk = seq_len(binmax)) %do% {

      bleft <- krange[1] + (bafbin.size * (bk - 1))
      bright <- krange[1] + (bafbin.size * bk)
      # inbk <- WESk$pos >= bleft & WESk$pos < bright
      inbk.full <- WESk$pos >= bleft & WESk$pos < bright
      # if(!any(inbk)) return()
      if(!any(inbk.full)) return()
      BAFk.full <- WESk$value[inbk.full]
      Xlo <- BAFk.full < BAF.cutter
      Xhi <- BAFk.full > (1 - BAF.cutter)
      Xin <- BAFk.full >= BAF.cutter & BAFk.full <= (1 - BAF.cutter)
      retvec <- rep(NA, length(BAFk.full))
      if(!any(Xin)) {
        retvec[Xlo] <- mcMin
        retvec[Xhi] <- mcMax
        return(retvec)
      }
      # BAFk <- WESk$value[inbk]
      BAFk <- BAFk.full[Xin]
      len.BAFk <- length(BAFk)
      # saveRDS(BAFk, "BAFk.RDS")
      # if (length(unique(BAFk)) == 1) return(rep(NA, len.BAFk))
      if (length(unique(BAFk)) == 1) return(retvec)
      set.seed(123456)
      if (is.null(prior)) {
        mcBAFk <- try(suppressWarnings(mclust::Mclust(BAFk, G = mcMin:mcMax, modelNames = modelName, verbose = FALSE)))
      } else  {
        mcBAFk <- try(suppressWarnings(mclust::Mclust(BAFk, G = mcMin:mcMax, modelNames = modelName, verbose = FALSE, prior = priorControl(scale = priorX$scale))))
      }
      ## Handling case of failed mclust
      # if (is.null(mcBAFk)) return(rep(NA, len.BAFk))
      if (is.null(mcBAFk)) return(retvec)
      # if (is.character(mcBAFk)) return(rep(NA, len.BAFk))
      if (is.character(mcBAFk)) return(retvec)

      kclassif <- mcBAFk$classification
      # if (all(is.na(kclassif))) return(rep(NA, len.BAFk))

      ### rescue method to restore symmetry in BAF mode
      if (max(kclassif, na.rm = TRUE) >= 3) {
        # c2range <- range(BAFk[kclassif == 2])
        cINrange <- range(BAFk[kclassif > min(kclassif, na.rm = TRUE) & kclassif < max(kclassif, na.rm = TRUE)])
        inv.cINrange <- 1 - cINrange
        objrange <- range(c(cINrange, inv.cINrange))
        inobj <- BAFk >= objrange[1] & BAFk <= objrange[2]
        kclassif[inobj] <- 2
      }

      # kuncertain <- mcBAFk$uncertainty
      # kclassif[kuncertain > 0.5] <- NA

      # if (length(unique(kclassif[!is.na(kclassif)])) == 1) {
      #   message(paste0(k, " ", bk, " ", length(BAFk), " ", unique(kclassif)))
      #   mcBAFk <- try(mclust::Mclust(BAFk, G = 2:mcMax, modelNames = modelName, verbose = FALSE), silent = TRUE)
      #   if (is.null(mcBAFk)) return(rep(NA, length(BAFk)))
      #   if (is.character(mcBAFk)) return(rep(NA, length(BAFk)))
      #   kclassif <- mcBAFk$classification
      #   kuncertain <- mcBAFk$uncertainty
      #   kclassif[kuncertain > 0.5] <- NA
      # }

      retvec[Xlo] <- min(kclassif)
      retvec[Xhi] <- max(kclassif)
      retvec[Xin] <- kclassif
      # return(kclassif)
      return(retvec)
    }
    return(mcresk)
  }

  unl.mcres <- unlist(mcreslist)
  # summary(unl.mcres)
  # table(unl.mcres)

  # message("MCLOGIC")
  mclogic <- foreach::foreach(b = mcreslist) %do%  {
    # if(!is.numeric(b)) message(b)
    if (length(b) == 0) return()
    ## Handling case of single class (NA by default)
    # if (length(unique(b[!is.na(b)])) == 1) return(rep(NA, length(b)))

    if(all(is.na(b))) return(rep(NA, length(b)))
    klogical <- rep(FALSE, length(b))
    klogical[b == max(b, na.rm = TRUE)] <- TRUE
    klogical[b == min(b, na.rm = TRUE)] <- TRUE
    # klogical[is.na(b)] <- NA
    # if (max(b, na.rm = TRUE) == 2) klogical[b == 2] <- NA
    return(klogical)
  }
  unl.logic <- unlist(mclogic)

  # unl.logic <- rep(TRUE, nrow(bafcdf.nna.complete))
  # unl.logic[baf.Xin] <- unl.logic.tmp

  ## Rate of hetero features
  hrate <- length(which(!unl.logic)) / length(unl.logic[!is.na(unl.logic)])
  # message(hrate)

  ## Filtering noisy features
  if (!is.null(BAF.filter)) {
    message(tmsg(paste0("Filtering BAF noise at ", BAF.filter * 100, "%...")))
    # allProbes = 1:length(Tumor_BAF_noNA)
    allProbes <- seq_along(unl.logic)
    # nonHomoProbes = allProbes[is.na(Hom) | Hom == FALSE]
    Hom <- unl.logic
    Hom[!Hom] <- NA
    names(Hom) <- names(Tumor_BAF_noNA)
    nonHomoProbes = allProbes[is.na(Hom) | Hom == FALSE]
    extraHetero <- round(length(nonHomoProbes) * BAF.filter)
    lowestDist <- NULL
    bsmHNA <- bsm[!is.na(bsm)]
    bsmHNA[!is.na(Hom) & Hom] = NA
    for (k in seq_along(ASCATobj$ch)) {
      # print(k)
      chrke <- ASCATobj$ch[[k]]
      chrNonHomoProbes = intersect(nonHomoProbes, chrke)
      if (length(chrNonHomoProbes) > 5) {
        segmentLength2 = min(length(chrNonHomoProbes) - 1, segmentLength)
        chrNonHomoProbesStartWindowLeft = c(rep(NA, segmentLength2), chrNonHomoProbes[1:(length(chrNonHomoProbes) - segmentLength2)])
        chrNonHomoProbesEndWindowLeft = c(NA, chrNonHomoProbes[1:(length(chrNonHomoProbes) - 1)])
        chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)], NA)
        chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2 + 1):length(chrNonHomoProbes)], rep(NA, segmentLength2))
        chrNonHomoProbesStartWindowMiddle = c(rep(NA, segmentLength2/2), chrNonHomoProbes[1:(length(chrNonHomoProbes) - segmentLength2/2)])
        chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2 + 1):length(chrNonHomoProbes)], rep(NA, segmentLength2/2))
        chrLowestDist = NULL
        for (probeNr in 1:length(chrNonHomoProbes)) {
          probe = chrNonHomoProbes[probeNr]
          if (!is.na(chrNonHomoProbesStartWindowLeft[probeNr]) & !is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
            medianLeft = stats::median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], na.rm = T)
          }
          else {
            medianLeft = NA
          }
          if (!is.na(chrNonHomoProbesStartWindowRight[probeNr]) & !is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
            medianRight = stats::median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]], na.rm = T)
          }
          else {
            medianRight = NA
          }
          if (!is.na(chrNonHomoProbesStartWindowMiddle[probeNr]) & !is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
            medianMiddle = stats::median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]), na.rm = T)
          }
          else {
            medianMiddle = NA
          }
          chrLowestDist[probeNr] = min(abs(medianLeft - bsmHNA[probe]), abs(medianRight - bsmHNA[probe]), abs(medianMiddle - bsmHNA[probe]), Inf, na.rm = T)
        }
      }
      else {
        chrLowestDist = NULL
        if (length(chrNonHomoProbes) > 0) {
          chrLowestDist[1:length(chrNonHomoProbes)] = 1
        }
      }
      lowestDist = c(lowestDist, chrLowestDist)
    }
    lowestDistUndecided <- lowestDist[is.na(Hom[nonHomoProbes])]
    names(lowestDistUndecided) = names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]
    sorted = sort(lowestDistUndecided)
    Hom[names(sorted[1:min(length(sorted), extraHetero)])] <-  FALSE
  }

  ## Rescuing homo regions
  message(tmsg("Rescuing putative homozygous BAF regions ..."))
  mcrescued <- foreach::foreach(b = mclogic) %do% {
    blen <- length(b)
    if(blen > 0) {
      if(all(!is.na(b))) {
        bhrate <- length(which(!b)) / length(b[!is.na(b)])
        if (bhrate < hrate / nfactor) b <- rep(FALSE, length(b))
      }
    }
    return(b)
  }
  unl.rescued <- unlist(mcrescued)

  ## Applying 95% BAF noise filtering
  unl.rescued[is.na(Hom)] <- TRUE

  # summary(unl.rescued)
  # table(unl.rescued)

  # message('PLOT')
  hoObj <- list(germlinegenotypes = matrix(NA, nrow = nrow(ASCATobj$Tumor_BAF), ncol = 1, dimnames = list(rownames(ASCATobj$SNPpos), ASCATobj$samples[1])), failedarrays = NULL)
  # message('chk1')
  hoObj$germlinegenotypes[bafcdf.nna$idx,1] <- unl.rescued
  # message('chk2')
  mcreslist.col <- unl.mcres
  mcreslist.col[is.na(mcreslist.col)] <- 4
  mcreslist.col <- c("black", "red", "green4", "purple", "blue", "pink")[mcreslist.col]
  # message('chk3')
  mclogic.col <- unl.logic
  mclogic.col <- as.numeric(mclogic.col)
  mclogic.col[is.na(mclogic.col)] <- 2
  mclogic.col <- c("red", "blue", "grey50")[(mclogic.col+1)]
  # message('chk4')
  mcres.col <- unl.rescued
  mcres.col <- as.numeric(mcres.col)
  mcres.col[is.na(mcres.col)] <- 2
  mcres.col <- c("red", "blue", "grey50")[(mcres.col+1)]
  # message('chk6')
  data(list = genome, package = "chromosomes", envir = environment())
  if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) {
    bafcdf.nna$genopos <- bafcdf.nna$pos + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[paste0("chr", bafcdf.nna$chrs)])]
  } else bafcdf.nna$genopos <- bafcdf.nna$pos + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[as.character(bafcdf.nna$chrs)])]
  # message(str(bafcdf.nna$genopos))
  # message('chk6')
  png(paste0(ASCATobj$samples[1], ".PredictedGermline.png"), width = 1700, height = 950)
  par(mfrow = c(3,1))
  # message('chk7')
  plot(bafcdf.nna$genopos, bafcdf.nna$value, col = mcreslist.col, pch = ".", yaxs = "i", xaxs = "i", ylim = c(0,1), cex = 5, xlab = "Genomic position", ylab = "BAF", main = paste0("Clustered BAF (", nrow(bafcdf.nna), ")"))
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lwd = 2)
  # message('chk8')
  plot(bafcdf.nna$genopos, bafcdf.nna$value, col = mclogic.col, pch = ".", yaxs = "i", xaxs = "i", ylim = c(0,1), cex = 5, xlab = "Genomic position", ylab = "BAF", main = paste0("Logical BAF (", length(which(unl.logic)), ")"))
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lwd = 2)
  # message('chk9')
  plot(bafcdf.nna$genopos, bafcdf.nna$value, col = mcres.col, pch = ".", yaxs = "i", xaxs = "i", ylim = c(0,1), cex = 5, xlab = "Genomic position", ylab = "BAF", main = paste0("Filtered, rescued BAF (", length(which(unl.rescued)), ")"))
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lwd = 2)
  # message('chk10')
  dev.off()

  # message("DONE!")
  return(hoObj)
}

