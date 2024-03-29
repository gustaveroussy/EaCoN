# EaCoN.BAF.Scaler.old <- function(ASCATobj = NULL, bafbin.size = 1E+07, toclustname = "BAF") {
#   
#   Homozygous = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
#   colnames(Homozygous) = colnames(ASCATobj$Tumor_LogR)
#   rownames(Homozygous) = rownames(ASCATobj$Tumor_LogR)
#   
#   Tumor_BAF_noNA = ASCATobj$Tumor_BAF[!is.na(ASCATobj$Tumor_BAF[, 1]), 1]
#   names(Tumor_BAF_noNA) = rownames(ASCATobj$Tumor_BAF)[!is.na(ASCATobj$Tumor_BAF[, 1])]
#   Tumor_LogR_noNA = ASCATobj$Tumor_LogR[names(Tumor_BAF_noNA), 1]
#   names(Tumor_LogR_noNA) <- names(Tumor_BAF_noNA)
#   bafcdf <- data.frame(idx = seq_len(nrow(ASCATobj$Tumor_BAF)), ASCATobj$SNPpos, value = ASCATobj$Tumor_BAF[,1], stringsAsFactors = FALSE)
#   tbsm <- bsm <- bafcdf$value
#   tbsm[which(tbsm < 0)] <- -tbsm[which(tbsm < 0)]
#   tbsm[which(tbsm > 1)] <- 2 - tbsm[which(tbsm > 1)]
#   bsm <- ifelse(bsm < .5, bsm, 1 - bsm)
#   if (toclustname == "BAF") bafcdf$value <- tbsm else if (toclustname == "mBAF") bafcdf$value <- bsm else stop(tmsg("Unknown toclustname value !"))
#   bafcdf.nna <- bafcdf[!is.na(bafcdf$value),]
#   
#   `%do%` <- foreach::"%do%"
#   scaledBAF <- foreach::foreach(k = unique(ASCATobj$SNPpos$chrs), .combine = "c") %do% {
#     WESk <- bafcdf.nna[bafcdf.nna$chrs == k,]
#     if(nrow(WESk) == 0) return(NULL)
#     krange <- range(WESk$pos, na.rm = TRUE)
#     binmax <- ceiling(diff(krange) / bafbin.size)
#     
#     ## Scaling
#     scaledBAF.k <- foreach::foreach(bk = seq_len(binmax), .combine = "c") %do% {
#       bleft <- krange[1] + (bafbin.size * (bk - 1))
#       bright <- krange[1] + (bafbin.size * bk)
#       inbk.full <- WESk$pos >= bleft & WESk$pos < bright
#       
#       if(!any(inbk.full)) return()
#       BAFk.full <- WESk$value[inbk.full]
#       
#       if (length(BAFk.full) >= 100) {
#         bafden <- density(BAFk.full, adjust = .5)
#         lowedge <- bafden$x[bafden$x < .25][which(bafden$y[bafden$x < .25] == max(bafden$y[bafden$x < .25]))]
#         highedge <- bafden$x[bafden$x > .75][which(bafden$y[bafden$x > .75] == max(bafden$y[bafden$x > .75]))]
#         BAFk.full <- (BAFk.full / (highedge-lowedge)) - (lowedge/highedge)
#       }
#       return(BAFk.full)
#     }
#     bafcdf.nna$value[bafcdf.nna$chrs == k] <- scaledBAF.k
#     return(bafcdf.nna$value[bafcdf.nna$chrs == k])
#   }
#   ASCATobj$Tumor_BAF_Unscaled <- ASCATobj$Tumor_BAF
#   ASCATobj$Tumor_BAF[,1][!is.na(ASCATobj$Tumor_BAF[,1])] <- scaledBAF
#   return(ASCATobj)
# }


## Rescale BAF
BAF.Rescale <- function(data = NULL, bafbin.size = 1E+07, toclustname = "BAF", out.dir = getwd(), return.data = FALSE, write.data = TRUE) {
  
  if (is.null(data)) stop(tmsg("data is NULL"), call. = FALSE)
  if (!is.null(data$data$Tumor_BAF_segmented)) stop(tmsg("BAF is already segmented : BAF.Rescale should be used on unsegmented data."), call. = FALSE)
  if (!dir.exists(out.dir)) stop(tmsg("out.dir does not exist !"), call. = FALSE)
  
  tmsg("Rescaling BAF ...")
  Tumor_BAF_noNA = data$data$Tumor_BAF[!is.na(data$data$Tumor_BAF[, 1]), 1]
  names(Tumor_BAF_noNA) = rownames(data$data$Tumor_BAF)[!is.na(data$data$Tumor_BAF[, 1])]
  Tumor_LogR_noNA = data$data$Tumor_LogR[names(Tumor_BAF_noNA), 1]
  names(Tumor_LogR_noNA) <- names(Tumor_BAF_noNA)
  bafcdf <- data.frame(idx = seq_len(nrow(data$data$Tumor_BAF)), data$data$SNPpos, value = data$data$Tumor_BAF[,1], stringsAsFactors = FALSE)
  tbsm <- bsm <- bafcdf$value
  tbsm[which(tbsm < 0)] <- -tbsm[which(tbsm < 0)]
  tbsm[which(tbsm > 1)] <- 2 - tbsm[which(tbsm > 1)]
  bsm <- ifelse(bsm < .5, bsm, 1 - bsm)
  if (toclustname == "BAF") bafcdf$value <- tbsm else if (toclustname == "mBAF") bafcdf$value <- bsm else stop(tmsg("Unknown toclustname value !"), call. = FALSE)
  bafcdf.nna <- bafcdf[!is.na(bafcdf$value),]
  
  `%do%` <- foreach::"%do%"
  scaledBAF <- foreach::foreach(k = unique(data$data$SNPpos$chrs), .combine = "c") %do% {
    WESk <- bafcdf.nna[bafcdf.nna$chrs == k,]
    if(nrow(WESk) == 0) return(NULL)
    krange <- range(WESk$pos, na.rm = TRUE)
    binmax <- ceiling(diff(krange) / bafbin.size)
    
    ## Scaling
    scaledBAF.k <- foreach::foreach(bk = seq_len(binmax), .combine = "c") %do% {
      bleft <- krange[1] + (bafbin.size * (bk - 1))
      bright <- krange[1] + (bafbin.size * bk)
      inbk.full <- WESk$pos >= bleft & WESk$pos < bright
      
      if(!any(inbk.full)) return()
      BAFk.full <- WESk$value[inbk.full]
      
      if (length(BAFk.full) >= 100) {
        bafden <- density(BAFk.full, adjust = .5)
        lowedge <- bafden$x[bafden$x < .25][which(bafden$y[bafden$x < .25] == max(bafden$y[bafden$x < .25]))]
        highedge <- bafden$x[bafden$x > .75][which(bafden$y[bafden$x > .75] == max(bafden$y[bafden$x > .75]))]
        BAFk.full <- (BAFk.full / (highedge-lowedge)) - (lowedge/highedge)
      }
      return(BAFk.full)
    }
    bafcdf.nna$value[bafcdf.nna$chrs == k] <- scaledBAF.k
    return(bafcdf.nna$value[bafcdf.nna$chrs == k])
  }
  data$data$Tumor_BAF_Unscaled <- data$data$Tumor_BAF
  data$data$Tumor_BAF[,1][!is.na(data$data$Tumor_BAF[,1])] <- scaledBAF
  data$meta$eacon$BAF.rescale <- "TRUE"
  
  if (return.data) return(data$data)
  
  ## Saving segmentation object
  if (write.data) {
    tmsg("Saving data ...")
    saveRDS(data, paste0(out.dir, "/", data$meta$basic$samplename, "_", data$meta$basic$type, "_", data$meta$basic$genome, "_processed.RDS"), compress = "xz")
  }
  
  tmsg("Done.")
}

## Rescale BAF

BAF.Rescale.ff <- function(RDSfile = NULL, ...) {
  if (is.null(RDS.file)) stop(tmsg("A RDS file is needed !"), call. = FALSE)
  if (!file.exists(RDS.file)) stop(tmsg(paste0("Could not find RDS file ", RDS.file, " !")), call. = FALSE)
  ## Data loading
  tmsg(paste0("Loading data from ", RDS.file, " ..."))
  my.data <- readRDS(RDS.file)
  Segment.ASCAT(data = my.data, out.dir = dirname(RDS.file), ...)
}

Segment.ASCAT.ff <- function(RDS.file = NULL, ...) {
  if (is.null(RDS.file)) stop(tmsg("A RDS file is needed !"), call. = FALSE)
  if (!file.exists(RDS.file)) stop(tmsg(paste0("Could not find RDS file ", RDS.file, " !")), call. = FALSE)
  ## Data loading
  tmsg(paste0("Loading data from ", RDS.file, " ..."))
  my.data <- readRDS(RDS.file)
  Segment.ASCAT(data = my.data, out.dir = dirname(RDS.file), ...)
}

EaCoN.Predict.Germline <- function(ASCATobj = NULL, bafbin.size = 1E+07, modelName = "E", toclustname = "BAF", mc.G = 2:4, nfactor = 4, BAF.filter = .9, BAF.cutter = 0, segmentLength = 5, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19") {

  # ## TEMP
  # setwd("/home/job/WORKSPACE/MP/CYTO/18R00201/")
  # ASCATobj = readRDS("/home/job/WORKSPACE/MP/CYTO/18R00201/20180220134039/18R00201.EaCoN.ASPCF.RDS")$data
  # bafbin.size = 1E+07
  # modelName = "E"
  # toclustname = "BAF"
  # mc.G = 2:4
  # nfactor = 4
  # segmentLength = 5
  # genome.pkg = "BSgenome.Hsapiens.UCSC.hg19"
  # BAF.filter <- .90
  # BAF.cutter <- 0
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # require(foreach)
  # require(mclust)

  ## Loading genome data
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  # data(list = genome, package = "chromosomes", envir = environment())
  message(tmsg(paste0("Loading ", genome.pkg, " ...")))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  # genome <- BSgenome::providerVersion(BSg.obj)
  genome <- metadata(BSg.obj)$genome
  cs <- chromobjector(BSg.obj)
  
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
  if (toclustname == "BAF") bafcdf$value <- tbsm else if (toclustname == "mBAF") bafcdf$value <- bsm else stop("Unknown toclustname value !", call. = FALSE)
  bafcdf.nna <- bafcdf[!is.na(bafcdf$value),]

  ## BAF.cutter modif
  baf.Xin <- bafcdf.nna$value >= BAF.cutter & bafcdf.nna$value <= (1-BAF.cutter)

  ## Computing genome-wide prior (if needed)
  set.seed(123456)
  message(tmsg("Computing whole genome BAF clustering prior ..."))
  priorX <- mclust::defaultPrior(bafcdf.nna$value[baf.Xin], G = mc.G, modelName = modelName)

  ## Computing banded BAF clusters
  message(tmsg("BAF clustering ..."))
  `%do%` <- foreach::"%do%"
  mcreslist <- foreach::foreach(k = unique(ASCATobj$SNPpos$chrs), .combine = "c") %do% {
    WESk <- bafcdf.nna[bafcdf.nna$chrs == k,]
    if(nrow(WESk) == 0) return(NULL)
    krange <- range(WESk$pos, na.rm = TRUE)
    binmax <- ceiling(diff(krange) / bafbin.size)

    ## Legacy code when prior was chr-based
    # if (!is.null(prior)) {
    #   if (prior == "K") {
        set.seed(123456)
    #     
        priorX <- mclust::defaultPrior(WESk$value, G = mc.G, modelName = modelName)
    #   }
    # }

    ## Clustering
    mcresk <- foreach::foreach(bk = seq_len(binmax)) %do% {

      bleft <- krange[1] + (bafbin.size * (bk - 1))
      bright <- krange[1] + (bafbin.size * bk)
      
      inbk.full <- WESk$pos >= bleft & WESk$pos < bright
      
      if(!any(inbk.full)) return()
      BAFk.full <- WESk$value[inbk.full]
      
      ## Scaling 
      if (length(BAFk.full) >= 100) {
        bafden <- density(BAFk.full, adjust = .5)
        lowedge <- bafden$x[bafden$x < .25][which(bafden$y[bafden$x < .25] == max(bafden$y[bafden$x < .25]))]
        highedge <- bafden$x[bafden$x > .75][which(bafden$y[bafden$x > .75] == max(bafden$y[bafden$x > .75]))]
        BAFk.full <- (BAFk.full / (highedge-lowedge)) - (lowedge/highedge)
      }
      
      ## Low-pass filtering
      Xlo <- BAFk.full < BAF.cutter
      Xhi <- BAFk.full > (1 - BAF.cutter)
      Xin <- BAFk.full >= BAF.cutter & BAFk.full <= (1 - BAF.cutter)
      retvec <- rep(NA, length(BAFk.full))
      if(!any(Xin)) {
        # retvec[Xlo] <- min(mc.G)
        # retvec[Xhi] <- max(mc.G)
        retvec[!Xin] <- min(mc.G)
        return(retvec)
      }
      
      ## Clustering BAF fragment
      BAFk <- BAFk.full[Xin]
      len.BAFk <- length(BAFk)

      if (length(unique(BAFk)) == 1) return(retvec)
      set.seed(123456)
      # if (is.null(prior)) {
        # mcBAFk <- try(suppressWarnings(mclust::Mclust(BAFk, G = mcMin:mcMax, modelNames = modelName, verbose = FALSE)))
        # mcBAFk <- try(suppressWarnings(mclust::Mclust(BAFk, G = mc.G, modelNames = modelName, verbose = FALSE)))
      # } else  {
        # mcBAFk <- try(suppressWarnings(mclust::Mclust(BAFk, G = mcMin:mcMax, modelNames = modelName, verbose = FALSE, prior = priorControl(scale = priorX$scale))))
        mcBAFk <- try(suppressWarnings(mclust::Mclust(BAFk, G = mc.G, modelNames = modelName, verbose = FALSE, prior = priorControl(scale = priorX$scale))), silent = TRUE)
      # }
      ## Handling case of failed mclust
      if (is.null(mcBAFk)) return(retvec)
      if (is.character(mcBAFk)) return(retvec)

      kclassif <- mcBAFk$classification

      ### Rescue method to restore symmetry in BAF mode
      if (max(kclassif, na.rm = TRUE) >= 3) {
        cINrange <- range(BAFk[kclassif > min(kclassif, na.rm = TRUE) & kclassif < max(kclassif, na.rm = TRUE)])
        inv.cINrange <- 1 - cINrange
        objrange <- range(c(cINrange, inv.cINrange))
        inobj <- BAFk >= objrange[1] & BAFk <= objrange[2]
        kclassif[inobj] <- 2
      }

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

  hoObj <- list(germlinegenotypes = matrix(NA, nrow = nrow(ASCATobj$Tumor_BAF), ncol = 1, dimnames = list(rownames(ASCATobj$SNPpos), ASCATobj$samples[1])), failedarrays = NULL)
  hoObj$germlinegenotypes[bafcdf.nna$idx,1] <- unl.rescued
  mcreslist.col <- unl.mcres
  mcreslist.col[is.na(mcreslist.col)] <- 4
  mcreslist.col <- c("black", "red", "green4", "purple", "blue", "pink")[mcreslist.col]
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
  # data(list = genome, package = "chromosomes", envir = environment())
  # if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) {
  #   bafcdf.nna$genopos <- bafcdf.nna$pos + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[paste0("chr", bafcdf.nna$chrs)])]
  # } else bafcdf.nna$genopos <- bafcdf.nna$pos + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[as.character(bafcdf.nna$chrs)])]
  bafcdf.nna$genopos <- bafcdf.nna$pos + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[as.character(bafcdf.nna$chrs)])]
  # message(str(bafcdf.nna$genopos))
  # message('chk6')
  message(tmsg("Plotting BAF ..."))
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

