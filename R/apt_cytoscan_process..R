## Performs CS CEL processing
CS.Process <- function(CEL = NULL, samplename = NULL, dual.norm = FALSE, normal.diploid = FALSE, l2r.level = "weighted", gc.renorm = TRUE, gc.rda = NULL, wave.renorm = TRUE, wave.rda = NULL, mingap = 1E+06, out.dir = getwd(), oschp.keep = FALSE, force.OS = NULL, apt.version = "2.4.0", apt.build = "na33.r4", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  setwd("/home/job/WORKSPACE/EaCoN_tests/CSHD")
  CEL = "M2568_K02.CEL.bz2"
  samplename = "M2568_K02"
  dual.norm = FALSE
  normal.diploid = FALSE
  l2r.level = "normal"
  wave.renorm = TRUE
  wave.rda <- NULL
  gc.renorm = TRUE
  gc.rda <- NULL
  BAF.filter <- .75
  # BAF.binsize = 1E+07
  mingap = 1E+06
  out.dir = getwd()
  oschp.keep = TRUE
  force.OS = NULL
  apt.version = "2.4.0"
  apt.build = "na33.r4"
  return.data = FALSE
  write.data = TRUE
  plot = TRUE
  force = FALSE
  genome.pkg = "BSgenome.Hsapiens.UCSC.hg19"
  require(foreach)
  source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  source("~/git_gustaveroussy/EaCoN/R/renorm_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/germline_functions.R")

  
  ## Early checks
  if (is.null(CEL)) stop(tmsg("A CEL file is required !"))
  if (is.null(samplename)) stop(tmsg("A samplename is required !"))
  if (!file.exists(CEL)) stop(tmsg(paste0("Could not find CEL file ", CEL, " !")))
  if (gc.renorm) { if (!is.null(gc.rda)) { if (!file.exists(gc.rda)) stop(tmsg(paste0("Could not find gc.rda file ", gc.rda))) } }
  if (wave.renorm) { if (!is.null(wave.rda)) { if (!file.exists(wave.rda)) stop(tmsg(paste0("Could not find wave.rda file ", wave.rda))) } }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"))
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !'))) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  l2r.lev.conv <- list("normal" = "Log2Ratio", "weighted" = "WeightedLog2Ratio")
  if (!(l2r.level %in% names(l2r.lev.conv))) stop(tmsg("Option 'l2r.level' should be 'normal' or 'weighted' !"))
  
  
  ## Handling compressed files
  CEL <- compressed_handler(CEL)
  
  ## Secondary checks
  sup.array <- c("CytoScanHD_Array", "CytoScan750k_Array")
  arraytype.cel = affxparser::readCelHeader(filename = CEL)$chiptype
  if (!arraytype.cel %in% sup.array) stop(tmsg(paste0("Identified array type '", arraytype.cel, "' is not supported by this function !")))
  
  ## Checking APT version compatibility
  valid.apt.versions <- c("2.4.0")
  if (!(apt.version %in% valid.apt.versions)) tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !"))
  
  ## Checking build compatibility
  valid.builds <- c("na33.r1", "na33.r2", "na33.r4", "na36.r1")
  if (!(tolower(apt.build) %in% valid.builds)) tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !"))
  
  ## Checking apt-copynumber-cyto-ssa package loc
  apt.cyto.pkg.name <- paste0("apt.cytoscan.", apt.version)
  if (!(apt.cyto.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", apt.cyto.pkg.name, " not found !")))
  suppressPackageStartupMessages(require(package = apt.cyto.pkg.name, character.only = TRUE))
  
  ## Processing CEL to an OSCHP file
  # if (dir.exists(samplename) && force) unlink(samplename, recursive = TRUE, force = FALSE)
  oscf <- apt.cytoscan.process(CEL = CEL, samplename = samplename, dual.norm = dual.norm, normal.diploid = normal.diploid, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)
  
  ## Reading OSCHP
  my.oschp <- oschp.load(file = oscf)
  sex.chr <- c("chrX", "chrY")
  
  ## Processing : meta (and checks)
  # meta.a1.df <- my.oschp[["Dset_IO_HDF5_Gdh"]][["_&keyvals"]][,1:2]
  # meta.a1.df <- meta.a1.df[!duplicated(meta.a1.df$key),]
  # meta.a1 <- meta.df2list(meta.a1.df)
  # meta.a1[["L2R.level"]] <- l2r.level
  # rm(meta.a1.df)
  if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA
  
  ### Loading genome info
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome2 <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)
  
  ### Getting genome build version
  genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  if (genome != genome2) stop(tmsg(paste0("Genome build name given with BSgenome package '", genome.pkg, "', (", genome2, ") is different from the genome build specified by provided APT build version '", apt.build, "' (", genome, ") !")))
  # data(list = genome, package = "chromosomes", envir = environment())
  arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
  manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
  species <- getmeta("affymetrix-algorithm-param-genome-species", my.oschp$Meta$analysis)
  
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = "NA")
  pgender <- gender.conv[[(getmeta("affymetrix-chipsummary-Y-gender-call", my.oschp$Meta$analysis))]]
  
  if (!(arraytype %in% sup.array)) stop(tmsg(paste0("Unsupported array : '", arraytype, "' !")))
  # meta.a2.df <- my.oschp[["Dset_IO_HDF5_Gdh"]][["Dset_IO_HDF5_Gdh:0"]][["Dset_IO_HDF5_Gdh:0:0"]][["_&keyvals"]][,1:2]
  # meta.a3.df <- my.oschp[["Dset_IO_HDF5_Gdh"]][["Dset_IO_HDF5_Gdh:0"]][["Dset_IO_HDF5_Gdh:0:0"]][["Dset_IO_HDF5_Gdh:0:0:0"]][["_&keyvals"]][,1:2]
  # meta.a2 <- meta.df2list(meta.a2.df)
  # rm(meta.a2.df)
  # meta.a3 <- meta.df2list(meta.a3.df)
  # rm(meta.a3.df)
  
  ## Reconstructing missing meta
  if (!"CEL1" %in% names(my.oschp$Meta)) {
    datheader.split <- unlist(strsplit(x =  affxparser::readCelHeader(filename = CEL)$datheader, split = "\\s+"))
    my.oschp$Meta$CEL1$acquisition <- list("affymetrix-scanner-id" = datheader.split[8], "affymetrix-scan-date" = paste0(datheader.split[6:7], collapse = " "))
    my.oschp$Meta$CEL1$array <- list("affymetrix-array-id" = NA, "affymetrix-array-barcode" = NA)
  }
  
  meta.b <- list(
    samplename = samplename,
    source = "microarray",
    source.file = CEL,
    type = arraytype,
    manufacturer = manufacturer,
    species = species,
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = pgender
  )
  
  ## Extracting data : L2R
  # ao.df <- data.frame(chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R.ori = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), L2R = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), BAF = NA, AD = NA, CallF = NA, stringsAsFactors = FALSE)
  # # ao.df <- if (l2r.level == "normal") {
  # #   data.frame(chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R = as.vector(my.oschp$ProbeSets$CopyNumber$Log2Ratio), BAF = NA, stringsAsFactors = FALSE)
  # # } else if (l2r.level == "weighted") {
  # #   data.frame(chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R = as.vector(my.oschp$ProbeSets$CopyNumber$WeightedLog2Ratio), BAF = NA, stringsAsFactors = FALSE)
  # # } else stop(tmsg("Unrecognized value for [l2r.level] !"))
  # rownames(ao.df) <- my.oschp$ProbeSets$CopyNumber$ProbeSetName
  # affy.chrom <- my.oschp$Chromosomes$Summary
  # ak <- affy.chrom$Display
  # names(ak) <- affy.chrom$Chromosome
  # ao.df$chrA <- as.vector(ak[as.character(ao.df$chrs)])
  # ao.df$chr <- paste0("chr", ao.df$chrA)
  # ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  # ao.df <- ao.df[order(ao.df$chrN, ao.df$pos, rownames(ao.df)),]
  
  ao.df <- dplyr::as.tbl(data.frame(my.oschp$ProbeSets$CopyNumber[,c(1:3)], L2R.ori = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]])))
  ao.df$L2R <- ao.df$L2R.ori
  affy.chrom <- my.oschp$Chromosomes$Summary
  ak <- affy.chrom$Display
  names(ak) <- affy.chrom$Chromosome
  ao.df$chrA <- as.vector(ak[as.character(ao.df$Chromosome)])
  ao.df$chr <- paste0("chr", ao.df$chrA)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  
  ## Normalizing SNPs
  tmsg("Normalizing SNP data (using rcnorm) ...")
  baf.df <- rcnorm::rcnorm.snp(myCEL = CEL, genome = genome, allSNPs = FALSE)
  baf.df$chr <- paste0("chr", baf.df$chrs)
  baf.df$chrN <- unlist(cs$chrom2chr[baf.df$chr])
  baf.df <- baf.df[order(baf.df$chrN, baf.df$pos),]
  baf.df <- baf.df[!is.na(baf.df$BAF),]
  gc()
  
  ao.df <- suppressWarnings(Reduce(function(t1, t2) dplyr::left_join(t1, t2, by = "ProbeSetName"), list(ao.df, dplyr::as.tbl(data.frame(ProbeSetName = rownames(baf.df), BAF = baf.df$BAF)), dplyr::as.tbl(my.oschp$ProbeSets$AllelicData[,c(1,4)]), dplyr::as.tbl(my.oschp$Genotyping$Calls[,c(2,5)]))))
  rm(baf.df)
  gc()
  
  ## Sorting by position
  ao.df <- dplyr::arrange(ao.df, chrN, Position, ProbeSetName)
  colnames(ao.df)[3] <- "pos"
  
  ## Filtering positions without L2R nor BAF
  ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)),]
  
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
  
  ## Identifying gaps and clustering chromosomal portions
  gaps <- which(diff(ao.df$pos) >= mingap)
  kends <- vapply(unique(ao.df$Chromosome), function(k) { max(which(ao.df$Chromosome == k)) }, 1L)
  kbreaks <- sort(unique(c(gaps, kends)))
  ao.df$chrgap <- rep(seq_along(kbreaks), times = c(kbreaks[1], diff(kbreaks)))
  
  ## Rescaling
  # if (BAF.rescale) {
  ao.df$BAF.unscaled <- ao.df$BAF
  for(g in unique(ao.df$chrgap)) {
    ing <- ao.df$chrgap == g
    lmed <- median(ao.df$BAF[ing][ao.df$ForcedCall[ing] == 6], na.rm = TRUE)
    umed <- median(ao.df$BAF[ing][ao.df$ForcedCall[ing] == 7], na.rm = TRUE)
    ao.df$BAF[ing] <- (lmed - ao.df$BAF[ing]) / (lmed - umed)
  }
  # }
  
  ## Preparing germline
  ao.df$germ <- ao.df$ForcedCall
  ao.df$germ[ao.df$germ %in% c(8,11)] <- 0
  ao.df$germ[ao.df$germ !=0 ] <- 1
  
  ## Filtering BAF
  # called <- which(ao.df$germ == 0 & !is.na(ao.df$BAF))
  # ao.df$mBAF <- BAF2mBAF(ao.df$BAF)
  # smoB <- round(length(called) / 3300)
  # if(smoB%%2 == 0) smoB <- smoB+1
  # mBAF.rm <- runmed(ao.df$mBAF[called], smoB)
  # mBAF.diff <- abs(ao.df$mBAF[called] - mBAF.rm)
  # Bfiltered <- mBAF.diff <= quantile(mBAF.diff, BAF.filter)
  # if (any(Bfiltered)) ao.df$germ[called][!Bfiltered] <- 1
  # rm(mBAF.rm, mBAF.diff, Bfiltered)
  # 
  ## Building ASCAT-like object
  tmsg("Building normalized object ...")
  # my.ch <- sapply(unique(ao.df$chr), function(x) { which(ao.df$chr == x) })
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR.ori = data.frame(sample = ao.df$L2R.ori, row.names = ao.df$ProbeSetName),
      Tumor_LogR = data.frame(sample = ao.df$L2R, row.names = ao.df$ProbeSetName),
      Tumor_BAF = data.frame(sample = ao.df$BAF, row.names = ao.df$ProbeSetName),
      Tumor_AD = data.frame(sample = ao.df$AllelicDifference, row.names = ao.df$ProbeSetName),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      SNPpos = data.frame(chrs = ao.df$chr, pos = ao.df$pos, row.names = ao.df$ProbeSetName),
      ch = sapply(unique(ao.df$chr), function(x) { which(ao.df$chr == x) }),
      chr = sapply(unique(ao.df$chrgap), function(x) { which(ao.df$chrgap == x) }),
      chrs = unique(ao.df$chr),
      samples = samplename,
      gender = as.vector(meta.b$predicted.gender),
      sexchromosomes = sex.chr,
      failedarrays = NULL
    ), 
    germline = list(
      germlinegenotypes = matrix(as.logical(ao.df$germ), ncol = 1),
      failedarrays = NULL
    )
  )
  colnames(my.ascat.obj$germline$germlinegenotypes) <- colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  rownames(my.ascat.obj$germline$germlinegenotypes) <- ao.df$ProbeSetName
  my.ascat.obj$data$Tumor_BAF.unscaled = data.frame(sample = ao.df$BAF.unscaled, row.names = ao.df$ProbeSetName)
  colnames(my.ascat.obj$data$Tumor_BAF.unscaled) <- samplename
  genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  rm(ao.df)
  gc()
  
  ## Adding meta
  my.ascat.obj$meta = list(
    basic = meta.b,
    affy = my.oschp$Meta
  )
  
  ## Adding CEL intensities
  my.ascat.obj$CEL = list(
    CEL1 = affxparser::readCel(filename = CEL)
  )
  my.ascat.obj$CEL$CEL1$intensities <- as.integer(my.ascat.obj$CEL$CEL1$intensities)
  
  
  
  # if(BAF.rescale) my.ascat.obj$data <- BAF.Rescale(data = my.ascat.obj, bafbin.size = BAF.binsize, toclustname = "BAF", return.data = TRUE, write.data = FALSE)
  # 
  if(write.data) saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "bzip2")

  ## PLOT
  if (plot) {
    tmsg("Plotting ...")
    kend <- genopos[vapply(my.ascat.obj$data$ch, max, 1L)]
    l2r.notna <- which(!is.na(my.ascat.obj$data$Tumor_LogR[,1]))
    l2r.rm <- runmed(my.ascat.obj$data$Tumor_LogR[,1][l2r.notna], smo)
    l2r.ori.rm <- runmed(my.ascat.obj$data$Tumor_LogR.ori[,1][l2r.notna], smo)
    png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
    par(mfrow = c(3,1))
    plot(genopos, my.ascat.obj$data$Tumor_LogR.ori[,1], pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " raw L2R profile (median-centered) / ", round(sum(abs(diff(l2r.ori.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
    lines(genopos[l2r.notna], l2r.ori.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(genopos, my.ascat.obj$data$Tumor_LogR[,1], pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " L2R profile (", l2r.level, ", median-centered)) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
    lines(genopos[l2r.notna], l2r.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = "grey80", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
    points(genopos[my.ascat.obj$germline$germlinegenotypes[,1] == 0], my.ascat.obj$data$Tumor_BAF[my.ascat.obj$germline$germlinegenotypes[,1] == 0,1], pch = ".", cex = 3, col = "grey25")
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
  
  tmsg("Done.")
  
  if(return.data) return(my.ascat.obj) else rm(my.ascat.obj)
}

CS.Process.Batch <- function(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the CEL.list.file
  if (is.null(CEL.list.file)) stop("A CEL.list.file is required !")
  if (!file.exists(CEL.list.file)) stop("Could not find CEL.list.file !")
  message("Reading and checking CEL.list.file ...")
  myCELs <- read.table(file = CEL.list.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("cel_files", "SampleName")
  head.chk <- all(colnames(CEL.list.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in CEL.list.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(myCELs)))
    stop("Invalid header.")
  }
  sn.chk <- duplicated(myCELs$SampleName)
  if (any(sn.chk)) {
    message("CEL.list.file contains duplicated SampleNames !")
    message(myCELs$SampleName[which(duplicated(myCELs$SampleName))])
    stop("Duplicated SampleNames.")
  }
  fecheck <- !vapply(myCELs$cel_files, file.exists, TRUE)
  fecheck.pos <- which(fecheck)
  if (length(fecheck.pos) > 0) stop(paste0("\n", "CEL file could not be found : ", myCELs$cel_files[fecheck.pos], collapse = ""))

  message(paste0("Found ", nrow(myCELs), " samples to process."))

  ## Adjusting cores/threads
  message("Adjusting number of threads if needed ...")
  avail.cores <- parallel::detectCores(logical = TRUE)
  if (is.null(nthread)) { nthread <- avail.cores -1; message(paste0("Reset nthread to ", nthread)) }
  if (nrow(myCELs) < nthread) { nthread <- nrow(myCELs); message(paste0("Reset nthread to ", nthread)) }
  if (avail.cores <= nthread) message(paste0(" WARNING : nthread set to ", nthread, " while available logical threads number is ", avail.cores, " !"))

  ## Building cluster
  current.bitmapType <- getOption("bitmapType")
  `%dopar%` <- foreach::"%dopar%"
  `%do%` <- foreach::"%do%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)

  p <- 0
  csres <- foreach::foreach(p = seq_len(nrow(myCELs)), .inorder = FALSE, .errorhandling = "pass") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    CS.Process(CEL = myCELs$cel_files[p], samplename = myCELs$SampleName[p], ...)
  }

  ## Stopping cluster
  message("Stopping cluster ...")
  parallel::stopCluster(cl)

  message("Done.")
  # if (return.data) return(csres)
}
