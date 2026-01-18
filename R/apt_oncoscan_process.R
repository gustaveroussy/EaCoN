## Performs OS CEL pairs processing
OS.Process <- function(ATChannelCel = NULL, GCChannelCel = NULL, samplename = NULL, dual.norm = TRUE, l2r.level = "weighted", gc.renorm = TRUE, gc.rda = NULL, wave.renorm = TRUE, wave.rda = NULL, mingap = 1E+06, out.dir = getwd(), oschp.keep = FALSE, force.OS = NULL, apt.version = "2.4.0", apt.build = "na33.r4", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  ## TEMP
  # setwd("/home/job/WORKSPACE/EaCoN_tests/OS")
  # ATChannelCel = "/home/job/WORKSPACE/EaCoN_tests/OS/18H03657_A.CEL.bz2"
  # GCChannelCel = "/home/job/WORKSPACE/EaCoN_tests/OS/18H03657_C.CEL.bz2"
  # samplename = "FACETSTEST"
  # dual.norm = TRUE
  # out.dir = getwd()
  # temp.files.keep = TRUE
  # force.OS = NULL
  # wave.renorm = TRUE
  # wave.rda <- NULL
  # gc.renorm = TRUE
  # gc.rda <- NULL
  # l2r.level <- "normal"
  # apt.build = "na33.r2"
  # apt.version = "2.4.0"
  # mingap = 1E+06
  # genome.pkg = "BSgenome.Hsapiens.UCSC.hg19"
  # BAF.filter = .75
  # force = FALSE
  # write.data = TRUE
  # plot = TRUE
  # return.data = FALSE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/renorm_functions.R")

  
  
  ## Early checks
  if (is.null(ATChannelCel)) stop(tmsg("An ATChannel CEL file is required !"), call. = FALSE)
  if (is.null(GCChannelCel)) stop(tmsg("A GCChannel CEL file is required !"), call. = FALSE)
  if (is.null(samplename)) stop(tmsg("A samplename is required !"), call. = FALSE)
  if (gc.renorm) { if (!is.null(gc.rda)) { if (!file.exists(gc.rda)) stop(tmsg(paste0("Could not find gc.rda file ", gc.rda)), call. = FALSE) } }
  if (wave.renorm) { if (!is.null(wave.rda)) { if (!file.exists(wave.rda)) stop(tmsg(paste0("Could not find wave.rda file ", wave.rda)), call. = FALSE) } }
  if(!file.exists(ATChannelCel)) stop(tmsg(paste0("Could not find ATChannelCel file ", ATChannelCel, " !")), call. = FALSE)
  if(!file.exists(GCChannelCel)) stop(paste0("Could not find GCChannelCel file ", GCChannelCel, " !"), call. = FALSE)
  if (ATChannelCel == GCChannelCel) stop(tmsg("ATChannelCel and GCChannelCel files are identical !"), call. = FALSE)
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !')), call. = FALSE) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  l2r.lev.conv <- list("normal" = "Log2Ratio", "weighted" = "WeightedLog2Ratio")
  if (!(l2r.level %in% names(l2r.lev.conv))) stop(tmsg("Option 'l2r.level' should be 'normal' or 'weighted' !"), call. = FALSE)
  
  ## Handling compressed files
  CEL.OS <- compressed_handler(c(ATChannelCel, GCChannelCel))
  ATChannelCel <- CEL.OS[1]
  GCChannelCel <- CEL.OS[2]
  
  ## Secondary checks
  sup.array <- c("OncoScan", "OncoScan_CNV")
  arraytype.celA = affxparser::readCelHeader(filename = ATChannelCel)$chiptype
  arraytype.celC = affxparser::readCelHeader(filename = GCChannelCel)$chiptype
  if (!arraytype.celA %in% sup.array) stop(tmsg(paste0("Identified array type for ATChannelCel file '", arraytype.celA, "' is not supported by this function !")), call. = FALSE)
  if (!arraytype.celC %in% sup.array) stop(tmsg(paste0("Identified array type for GCChannelCel file '", arraytype.celC, "' is not supported by this function !")), call. = FALSE)
  if (arraytype.celA != arraytype.celC) stop(tmsg(paste0("ATChannelCel and GCChannelCel files are not of the same design : ", arraytype.celA, " != ", arraytype.celC, " !")), call. = FALSE)
  
  ## Checking APT version compatibility
  valid.apt.versions <- c("2.4.0")
  if (!(apt.version %in% valid.apt.versions)) tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !"))
  
  ## Checking build compatibility
  valid.builds <- c("na33.r2", "na33.r4", "na36.r1")
  if (!(tolower(apt.build) %in% valid.builds)) tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !"))
  
  ## Checking apt-copynumber-onco-ssa package loc
  apt.onco.pkg.name <- paste0("apt.oncoscan.", apt.version)
  if (!(apt.onco.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", apt.onco.pkg.name, " not found !")), call. = FALSE)
  suppressPackageStartupMessages(require(package = apt.onco.pkg.name, character.only = TRUE))
  
  ## Processing CELs to an OSCHP file
  if (dir.exists(samplename) && force) unlink(samplename, recursive = TRUE, force = FALSE)
  oscf <- apt.oncoscan.process(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel, samplename = samplename, dual.norm = dual.norm, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)
  
  ## Reading OSCHP
  my.oschp <- oschp.load(file = oscf)
  if (length(names(my.oschp$Meta)) == 3) names(my.oschp$Meta) <- c("ATChannelCel", "GCChannelCel", "analysis")
  sex.chr <- c("chrX", "chrY")
  
  ## Processing : meta (and checks)
  if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA
  
  ### Loading genome info
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  # genome2 <- BSgenome::providerVersion(BSg.obj)
  genome2 <- metadata(BSg.obj)$genome
  cs <- chromobjector(BSg.obj)
  
  ### Getting basic meta
  genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  if (genome != genome2) stop(tmsg(paste0("Genome build name given with BSgenome package '", genome.pkg, "', (", genome2, ") is different from the genome build specified by provided APT build version '", apt.build, "' (", genome, ") !")), call. = FALSE)
  arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
  manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
  species <- getmeta("affymetrix-algorithm-param-genome-species", my.oschp$Meta$analysis)
  
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = "NA")
  pgender <- gender.conv[[(getmeta("affymetrix-chipsummary-Y-gender-call", my.oschp$Meta$analysis))]]
  
  if (!(arraytype %in% sup.array)) stop(tmsg(paste0("Unsupported array : '", arraytype, "' !")), call. = FALSE)
  
  meta.b <- list(
    samplename = samplename,
    source = "microarray",
    source.file = list(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel),
    type = arraytype,
    manufacturer = manufacturer,
    species = species,
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = pgender
  )
  
  ## Extracting data : L2R
  ao.df <- data.frame(ProbeSetName = my.oschp$ProbeSets$CopyNumber$ProbeSetName, chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R.ori = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), L2R = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), BAF = as.vector(my.oschp$ProbeSets$AllelicData$BAF), AD = as.vector(my.oschp$ProbeSets$AllelicData$AllelicDifference), CallF = as.vector(my.oschp$Genotyping$Calls$ForcedCall), ASignal = my.oschp$Genotyping$Calls$ASignal, BSignal = my.oschp$Genotyping$Calls$BSignal, stringsAsFactors = FALSE)
  affy.chrom <- my.oschp$Chromosomes$Summary
  ak <- affy.chrom$Display
  names(ak) <- affy.chrom$Chromosome
  ao.df$chrA <- as.vector(ak[as.character(ao.df$chrs)])
  ao.df$chr <- paste0("chr", ao.df$chrA)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos, ao.df$ProbeSetName),]
  ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)),]

  ## Adding specific data for FACETS
  rcmat <- round(cbind(ao.df$BAF * (ao.df$ASignal + ao.df$BSignal), (1-ao.df$BAF)*(ao.df$ASignal + ao.df$BSignal)))
  ao.df$LOR <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
  ao.df$LORvar <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
  rm(rcmat)
  gc()
  
  ## L2R renormalizations
  smo <- round(nrow(ao.df) / 550)
  if(smo%%2 == 0) smo <- smo+1
  
  ### Wave  
  if (wave.renorm) {
    tmsg("Wave re-normalization ...")
    
    ren.res <- renorm.go(input.data = ao.df, renorm.rda = wave.rda, track.type = "Wave", smo = smo, arraytype = arraytype, genome = genome)
    
    ao.df <- ren.res$data
    fitted.l2r <- ren.res$renorm$l2r$l2r
    
    if(is.null(ren.res$renorm$pos)) {
      meta.b <- setmeta("wave.renorm", "None", meta.b)
      tmsg(" No positive fit.")
    } else {
      ## Tweaking sex chromosomes
      sex.idx <- ao.df$chr %in% sex.chr
      auto.ori.med <- median(ao.df$L2R[!sex.idx], na.rm = TRUE)
      auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
      if (any(sex.idx)) {
        for (k in sex.chr) {
          k.idx <- ao.df$chr == k
          if (any(k.idx)) {
            k.ori.diffmed <- median(ao.df$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
            fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
          }
        }
      }
      meta.b <- setmeta("wave.renorm", paste0(ren.res$mrenorm$pos, collapse = ","), meta.b)
    }
    rm(ren.res)
    ao.df[["L2R.Wave"]] <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
    ao.df$L2R <- ao.df[["L2R.Wave"]]
  } else {
    meta.b <- setmeta("wave.renorm", "FALSE", meta.b)
  }
  
  ### GC
  if (gc.renorm) {
    tmsg("GC renormalization ...")
    
    ren.res <- renorm.go(input.data = ao.df, renorm.rda = gc.rda, track.type = "GC", smo = smo, arraytype = arraytype, genome = genome)
    ao.df <- ren.res$data
    fitted.l2r <- ren.res$renorm$l2r$l2r
    
    if(is.null(ren.res$renorm$pos)) {
      meta.b <- setmeta("gc.renorm", "None", meta.b)
      message(tmsg(" No positive fit."))
    } else {
      ## Tweaking sex chromosomes
      sex.idx <- ao.df$chr %in% sex.chr
      auto.ori.med <- median(ao.df$L2R[!sex.idx], na.rm = TRUE)
      auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
      if (any(sex.idx)) {
        for (k in sex.chr) {
          k.idx <- ao.df$chr == k
          if (any(k.idx)) {
            k.ori.diffmed <- median(ao.df$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
            fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
          }
        }
      }
      meta.b <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, collapse = ","), meta.b)
    }
    rm(ren.res)
    ao.df[["L2R.GC"]] <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
    ao.df$L2R <- ao.df[["L2R.GC"]]
  } else {
    meta.b <- setmeta("gc.renorm", "FALSE", meta.b)
  }
  
  ## Rough median-centering of L2R
  ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  ao.df$AD[is.nan(ao.df$AD)] <- NA
  
  ## Identifying gaps and clustering chromosomal portions
  gaps <- which(diff(ao.df$pos) >= mingap)
  kends <- vapply(unique(ao.df$chrs), function(k) { max(which(ao.df$chrs == k)) }, 1L)
  kbreaks <- sort(unique(c(gaps, kends)))
  ao.df$chrgap <- rep(seq_along(kbreaks), times = c(kbreaks[1], diff(kbreaks)))
  
  ## Preparing germline
  germ <- ao.df$CallF
  
  
  if(as.numeric(version$major) >= 4) {
    germ[germ %in% as.raw(c(8,11))] <- as.raw(0)
    germ[germ != 0 ] <- as.raw(1)
  } else {
    germ[germ %in% c(8,11)] <- 0
    germ[germ != 0 ] <- 1
  }
  
  ## Building ASCAT-like object
  tmsg("Building normalized object ...")
  my.ch <- sapply(unique(ao.df$chrs), function(x) { which(ao.df$chrs == x) })
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR.ori = data.frame(sample = ao.df$L2R.ori, row.names = ao.df$ProbeSetName),
      Tumor_LogR = data.frame(sample = ao.df$L2R, row.names = ao.df$ProbeSetName),
      Tumor_BAF = data.frame(sample = ao.df$BAF, row.names = ao.df$ProbeSetName),
      Tumor_AD = data.frame(sample = ao.df$AD, row.names = ao.df$ProbeSetName),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      # SNPpos = data.frame(chrs = ao.df$chrA, pos = ao.df$pos, row.names = rownames(ao.df)),
      SNPpos = data.frame(chrs = ao.df$chr, pos = ao.df$pos, row.names = ao.df$ProbeSetName),
      ch = sapply(unique(ao.df$chr), function(x) { which(ao.df$chr == x) }),
      chr = sapply(unique(ao.df$chrgap), function(x) { which(ao.df$chrgap == x) }),
      chrs = unique(ao.df$chr),
      samples = samplename,
      gender = as.vector(meta.b$predicted.gender),
      sexchromosomes = sex.chr,
      failedarrays = NULL,
      additional = data.frame(RD.test = ao.df$ASignal, RD.ref = ao.df$BSignal, LOR = ao.df$AD, LORvar = ao.df$LORvar, stringsAsFactors = FALSE)
    ),
    meta = list(
      basic = meta.b,
      affy = my.oschp$Meta
    ),
    germline = list(
      germlinegenotypes = matrix(as.logical(germ), ncol = 1),
      failedarrays = NULL
    ),
    CEL = list(
      ATChannelCel = affxparser::readCel(filename = ATChannelCel),
      GCChannelCel = affxparser::readCel(filename = GCChannelCel)
    )
  )
  colnames(my.ascat.obj$germline$germlinegenotypes) <- colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  rownames(my.ascat.obj$germline$germlinegenotypes) <- ao.df$ProbeSetName
  my.ascat.obj$CEL$ATChannelCel$intensities <- as.integer(my.ascat.obj$CEL$ATChannelCel$intensities)
  my.ascat.obj$CEL$GCChannelCel$intensities <- as.integer(my.ascat.obj$CEL$GCChannelCel$intensities)
  gc()
  
  if (write.data) saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "bzip2")
  
  ## Rough plot
  # if (plot) {
  #   tmsg("Plotting ...")
  #   ao.df$genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  #   ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  #   ao.df$L2R.ori <- ao.df$L2R.ori - median(ao.df$L2R.ori, na.rm = TRUE)
  #   kend <- ao.df$genopos[vapply(unique(ao.df$chrN), function(k) { max(which(ao.df$chrN == k))}, 1)]
  #   l2r.notna <- which(!is.na(ao.df$L2R))
  #   l2r.rm <- runmed(ao.df$L2R[l2r.notna], smo)
  #   l2r.ori.rm <- runmed(ao.df$L2R.ori[l2r.notna], smo)
  #   png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
  #   par(mfrow = c(3,1))
  #   plot(ao.df$genopos, ao.df$L2R.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " raw L2R profile (median-centered) / ", round(sum(abs(diff(l2r.ori.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
  #   lines(ao.df$genopos[l2r.notna], l2r.ori.rm, col = 1)
  #   abline(v = kend, col = 4, lty = 3, lwd = 2)
  #   abline(h = 0, col = 2, lty = 2, lwd = 2)
  #   plot(ao.df$genopos, ao.df$L2R, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " normalized L2R profile (median-centered) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
  #   lines(ao.df$genopos[l2r.notna], l2r.rm, col = 1)
  #   abline(v = kend, col = 4, lty = 3, lwd = 2)
  #   abline(h = 0, col = 2, lty = 2, lwd = 2)
  #   plot(ao.df$genopos, ao.df$BAF, pch = ".", cex = 3, col = ao.df$CallF-5, xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
  #   abline(v = kend, col = 4, lty = 3, lwd = 2)
  #   abline(h = .5, col = 2, lty = 2, lwd = 2)
  #   dev.off()
  # }
  
  genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  rm(ao.df)
  gc()
  if (plot) {
    tmsg("Plotting ...")
    kend <- genopos[vapply(my.ascat.obj$data$ch, max, 1L)]
    l2r.notna <- which(!is.na(my.ascat.obj$data$Tumor_LogR[,1]))
    l2r.rm <- runmed(my.ascat.obj$data$Tumor_LogR[,1][l2r.notna], smo)
    l2r.ori.rm <- runmed(my.ascat.obj$data$Tumor_LogR.ori[,1][l2r.notna], smo)
    png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
    par(mfrow = c(3,1))
    plot(genopos, my.ascat.obj$data$Tumor_LogR.ori[,1], pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " raw L2R profile (median-centered) / ", round(sum(abs(diff(l2r.ori.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
    lines(genopos[l2r.notna], l2r.ori.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(genopos, my.ascat.obj$data$Tumor_LogR[,1], pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " L2R profile (", l2r.level, ", median-centered)) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
    lines(genopos[l2r.notna], l2r.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
    points(genopos[germ == 0], my.ascat.obj$data$Tumor_BAF[germ == 0,1], pch = ".", cex = 3, col = 2)
    # plot(ao.df$genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = ao.df$ForcedCall-5, xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = .5, col = 2, lty = 2, lwd = 2)
    dev.off()
  }
  
  ## Cleaning
  if(!oschp.keep) {
    tmsg("Removing temporary OSCHP file ...")
    file.remove(oscf)
  }
  
  message(tmsg("Done."))
  gc()
  if(return.data) return(my.ascat.obj)
}


## Same in batch mode
OS.Process.Batch <- function(pairs.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the pairs.file
  if (!file.exists(pairs.file)) stop("Could not find pairs.file !", call. = FALSE)
  message("Reading and checking pairs.file ...")
  mypairs <- read.table(file = pairs.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("ATChannelCel", "GCChannelCel", "SampleName")
  head.chk <- all(colnames(pairs.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in pairs.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(mypairs)))
    stop("Invalid header.", call. = FALSE)
  }
  
  at.chk <- file.exists(mypairs$ATChannelCel)
  gc.chk <- file.exists(mypairs$GCChannelCel)
  
  if (!all(at.chk) || !all(at.chk)) {
    message("Some CEL files from the pairs.file could not be found (wrong path ?) !")
    message("Missing AT CEL files :")
    message(paste0(mypairs$ATChannelCel[which(!at.chk)], collapse = ", "))
    message("Missing GC CEL files :")
    message(paste0(mypairs$GCChannelCel[which(!gc.chk)], collapse = ", "))
    stop("Missing CEL.", call. = FALSE)
  }
  sn.chk <- duplicated(mypairs$SampleName)
  if (any(sn.chk)) {
    message("pairs.file contains duplicated SampleNames !")
    message(mypairs$SampleName[which(duplicated(mypairs$SampleName))])
    stop("Duplicated SampleNames.", call. = FALSE)
  }
  if(any(mypairs$ATChannelCel == mypairs$GCChannelCel)) {
    message("Some ATChannelCel and GCChannelCel are the same for at least one sample !")
    stop("Identical CEL files for AT and GC.", call. = FALSE)
  }
  
  message(paste0("Found ", nrow(mypairs), " samples to process."))
  
  ## Adjusting cores/threads
  message("Adjusting number of cores if needed ...")
  if (is.null(nthread)) nthread <- parallel::detectCores(logical = TRUE) -1
  if (nrow(mypairs) < nthread) nthread <- nrow(mypairs)
  
  current.bitmapType <- getOption("bitmapType")
  message("Running APT ...")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  
  p <- 0
  osres <- foreach(p = seq_len(nrow(mypairs)), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    OS.Process(ATChannelCel = mypairs$ATChannelCel[p], GCChannelCel = mypairs$GCChannelCel[p], samplename = mypairs$SampleName[p], ...)
  }
  message("Stopping cluster ...")
  parallel::stopCluster(cl)
  
  message("Done.")
}

