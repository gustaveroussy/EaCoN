## Performs CS CEL processing
SNP6.Process <- function(CEL = NULL, samplename = NULL, l2r.level = "normal", gc.renorm = TRUE, gc.rda = NULL, wave.renorm = TRUE, wave.rda = NULL, mingap = 1E+06, out.dir = getwd(), oschp.keep = FALSE, force.OS = NULL, apt.version = "1.20.0", apt.build = "na35.r1", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {

  # setwd("/home/job/WORKSPACE/EaCoN_tests/SNP6")
  # CEL <- "GSM820994.CEL.bz2"
  # samplename <- "BITES_TEST"
  # l2r.level <- "normal"
  # wave.renorm <- TRUE
  # wave.rda <- NULL
  # gc.renorm <- TRUE
  # gc.rda <- NULL
  # mingap <- 1E+06
  # out.dir <- getwd()
  # oschp.keep <- TRUE
  # force.OS <- NULL
  # apt.version <- "1.20.0"
  # apt.build <- "na35.r1"
  # genome.pkg <- "BSgenome.Hsapiens.UCSC.hg19"
  # return.data <- FALSE
  # write.data <- TRUE
  # plot <- TRUE
  # force <- FALSE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/renorm_functions.R")


  ## Early checks
  if (is.null(CEL)) stop(tmsg("A CEL file is required !"), call. = FALSE)
  if (is.null(samplename)) stop(tmsg("A samplename is required !"), call. = FALSE)
  if (!file.exists(CEL)) stop(tmsg(paste0("Could not find CEL file ", CEL, " !")), call. = FALSE)
  if (gc.renorm) { if (!is.null(gc.rda)) { if (!file.exists(gc.rda)) stop(tmsg(paste0("Could not find gc.rda file ", gc.rda)), call. = FALSE) } }
  if (wave.renorm) { if (!is.null(wave.rda)) { if (!file.exists(wave.rda)) stop(tmsg(paste0("Could not find wave.rda file ", wave.rda)), call. = FALSE) } }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"), call. = FALSE)
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !')), call. = FALSE) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  l2r.lev.conv <- list("normal" = "Log2Ratio", "weighted" = "SmoothSignal")
  if (!(l2r.level %in% names(l2r.lev.conv))) stop(tmsg("Option 'l2r.level' should be 'normal' or 'weighted' !"), call. = FALSE)
  
  ## Handling compressed files
  CEL <- compressed_handler(CEL)
  
  ## Secondary checks
  sup.array <- c("GenomeWideSNP_6")
  arraytype.cel = affxparser::readCelHeader(filename = CEL)$chiptype
  if (!arraytype.cel %in% sup.array) stop(tmsg(paste0("Identified array type '", arraytype.cel, "' is not supported by this function !")), call. = FALSE)

  ## Checking APT version compatibility
  valid.apt.versions <- c("1.20.0")
  if (!(apt.version %in% valid.apt.versions)) warning(tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !")))

  ## Checking build compatibility
  valid.builds <- c("na35.r1")
  if (!(tolower(apt.build) %in% valid.builds)) warning(tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !")))

  ## Checking apt-copynumber-cyto-ssa package loc
  apt.snp6.pkg.name <- paste0("apt.snp6.", apt.version)
  if (!(apt.snp6.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", apt.snp6.pkg.name, " not found !")), call. = FALSE)
  suppressPackageStartupMessages(require(package = apt.snp6.pkg.name, character.only = TRUE))

  ## Processing CEL to an OSCHP file
  oscf <- apt.snp6.process(CEL = CEL, samplename = samplename, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)

  ## Reading OSCHP
  my.oschp <- oschp.load(file = oscf)
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
  
  snp6.conv <- list("1" = "male", "2" = "female", "NA" = "NA", "0" = "NA")
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = "NA")
  pgender <- gender.conv[[snp6.conv[[as.character(as.numeric(getmeta("affymetrix-chipsummary-Gender", my.oschp$Meta$analysis)))]]]]

  if (!(arraytype %in% sup.array)) stop(paste0("Unsupported array : '", arraytype, "' !"), call. = FALSE)

  ## Reconstructing missing meta
  if (!"CEL1" %in% names(my.oschp$Meta)) {
    datheader.split <- unlist(strsplit(x =  affxparser::readCelHeader(filename = CEL)$datheader, split = "\\s+"))
    my.oschp$Meta$CEL1$acquisition <- list("affymetrix-scanner-id" = datheader.split[8], "affymetrix-scan-date" = paste0(datheader.split[6:7], collapse = " "))
    my.oschp$Meta$CEL1$array <- list("affymetrix-array-id" = NA, "affymetrix-array-barcode" = NA)
  }
  # meta.a2 <- list("affymetrix-scanner-id" = scanner.id, "affymetrix-scan-date" = scan.date)
  # meta.a2 <- list("affymetrix-scanner-id" = scanner.id, "affymetrix-scan-date" = scan.date)
  # meta.a3 <- list("affymetrix-array-id" = NA, "affymetrix-array-barcode" = NA)
  # meta.a3 <- list("affymetrix-array-id" = NA, "affymetrix-array-barcode" = NA)

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
  ao.df <- dplyr::as.tbl(data.frame(my.oschp$MultiData$CopyNumber[,c(2:4)], L2R.ori = as.vector(my.oschp$MultiData$CopyNumber[[l2r.lev.conv[[l2r.level]]]])))
  ao.df$Chromosome <- as.integer(ao.df$Chromosome) ### Patching the Chromosome column : on R4+, it is read as 'raw', we need ints
  affy.meta <- my.oschp$Meta
  affy.chrom <- my.oschp$MultiData[["CopyNumber_&keyvals"]][seq.int(3, nrow(my.oschp$MultiData[["CopyNumber_&keyvals"]]), 3),1:2]
  ao.df$L2R <- ao.df$L2R.ori
  ak <- affy.chrom$val
  names(ak) <- as.numeric(sub(":display", "", affy.chrom$key))
  ao.df$chrA <- as.vector(ak[as.character(ao.df$Chromosome)])
  ao.df$chr <- paste0("chr", ao.df$chrA)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  
  ## Normalizing SNPs
  tmsg("Normalizing SNP data (using rcnorm) ...")
  baf.df <- rcnorm::rcnorm.snp(myCEL = CEL, genome.pkg = genome.pkg, allSNPs = FALSE)
  baf.df$chr <- paste0("chr", baf.df$chrs)
  baf.df$chrN <- unlist(cs$chrom2chr[baf.df$chr])
  baf.df <- baf.df[order(baf.df$chrN, baf.df$pos),]
  baf.df <- baf.df[!is.na(baf.df$BAF),]
  gc()
  
  ao.df <- suppressWarnings(Reduce(function(t1, t2) dplyr::left_join(t1, t2, by = "ProbeSetName"), list(ao.df, dplyr::as.tbl(data.frame(ProbeSetName = rownames(baf.df), BAF.ori = baf.df$BAF, BAF = baf.df$BAF)), dplyr::as.tbl(my.oschp$MultiData$CopyNumber[,c(2,9)]))))
  rm(my.oschp, baf.df)
  gc()
  
  ## Hacking Type of the 'Allele Difference' column (from array to vector)
  ao.df[['Allele Difference']] <- as.vector(ao.df[['Allele Difference']])
  ao.df <- dplyr::arrange(ao.df, chrN, Position, ProbeSetName)
  colnames(ao.df)[3] <- "pos"
  
  ## Merging L2R and BAF data
  ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)),]
  
  
  ######################
  ######################
  ##### WAR ZONE ! #####
  ######################
  ######################
  
  ########################
  ## Segmentation-based ##
  ########################
  minseglen <- 50
  # nna <- !is.na(ao.df$BAF)
  nna <- !is.na(ao.df$BAF) & !is.na(ao.df[["Allele Difference"]])
  ao.df$mBAF <- BAF2mBAF(ao.df$BAF)
  # str(peltres <- changepoint::cpt.meanvar(ao.df$mBAF[nna], penalty = "MBIC", method = "PELT", minseglen = minseglen)@cpts)
  peltres <- changepoint::cpt.var(ao.df[["Allele Difference"]][nna], penalty = "BIC", method = "PELT", minseglen = minseglen)@cpts
  
  ## Refining with chr ends
  kends <- vapply(unique(ao.df$chrN[nna]), function(k) { max(which(ao.df$chrN[nna] == k))}, 1L)
  peltres <- sort(unique(c(peltres, kends)))
  prdiff <- diff(peltres)
  toremove <- vector()
  if (any(prdiff < 50)) {
    for (b in which(prdiff < 50)) {
      if (peltres[b+1] %in% kends) toremove <- c(toremove, b) else toremove <- c(toremove, b+1)
    }
    peltres <- peltres[-toremove]
  }
  
  ## Clustering BAF segments
  bs.end <- peltres
  bs.start <- c(1, peltres[-length(peltres)]+1)
  
  ao.df$cluster <- NA
  ao.df$uni <- FALSE
  suppressPackageStartupMessages(library(mclust))
  mc.G <- 4
  mc.mN <- "E"
  smeds <- hrates <- vector()
  for(i in seq_along(bs.start)) {
    mcresBIC <- mclust::mclustBIC(data = ao.df$BAF[nna][bs.start[i]:bs.end[i]], G = mc.G, modelNames = mc.mN, verbose = FALSE)
    mcres <- mclust::Mclust(data = ao.df$BAF[nna][bs.start[i]:bs.end[i]], G = mc.G, modelNames = mc.mN, verbose = FALSE, x = mcresBIC)$classification
    ao.df$cluster[nna][bs.start[i]:bs.end[i]] <- mcres
    if (length(unique(mcres)) == 2) ao.df$uni[nna][bs.start[i]:bs.end[i]] <- TRUE
    smeds <- c(smeds, median(ao.df$mBAF[nna][bs.start[i]:bs.end[i]][!mcres %in% c(1,4)], na.rm = TRUE))
    hrates <- c(hrates, length(which(!mcres %in% c(1,4))) / length(mcres))
  }
  
  ## Rescaling
  tmsg("Rescaling BAF ...")
  ao.df$BAF.unscaled <- ao.df$BAF
  for(i in seq_along(bs.start)) {
    lmed <- median(ao.df$BAF.unscaled[nna][bs.start[i]:bs.end[i]][ao.df$cluster[nna][bs.start[i]:bs.end[i]] == 1], na.rm = TRUE)
    umed <- median(ao.df$BAF.unscaled[nna][bs.start[i]:bs.end[i]][ao.df$cluster[nna][bs.start[i]:bs.end[i]] == 4], na.rm = TRUE)
    ao.df$BAF[nna][bs.start[i]:bs.end[i]] <- (lmed - ao.df$BAF[nna][bs.start[i]:bs.end[i]]) / (lmed - umed)
  }
  ao.df$mBAF <- BAF2mBAF(ao.df$BAF)
  
  ## Rescue
  ao.df$cluster2 <- ao.df$cluster
  ao.df$cluster2[ao.df$cluster2 == 4] <- 1
  ao.df$cluster2[ao.df$cluster2 > 1] <- 2
  target.hrate <- .3
  for(i in seq_along(hrates)) {
    if (hrates[i] < .2 & smeds[i] < .45) {
      ao.df$cluster2[nna][bs.start[i]:bs.end[i]] <- 2
      ao.df$uni[nna][bs.start[i]:bs.end[i]] <- TRUE
    }
  }
  
  ## L2R renormalizations
  tmsg("Normalizing L2R data ...")
  
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
  
  ## Building ASCAT-like object
  tmsg("Building normalized object ...")
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR.ori = data.frame(sample = ao.df$L2R.ori, row.names = ao.df$ProbeSetName),
      Tumor_LogR = data.frame(sample = ao.df$L2R, row.names = ao.df$ProbeSetName),
      Tumor_BAF = data.frame(sample = ao.df$BAF, row.names = ao.df$ProbeSetName),
      Tumor_AD = data.frame(sample = ao.df[["Allele Difference"]], row.names = ao.df$ProbeSetName),
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
      germlinegenotypes = matrix(as.logical(abs(ao.df$cluster2 - 2L)), ncol = 1),
      failedarrays = NULL
    )
  )
  colnames(my.ascat.obj$germline$germlinegenotypes) <- colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  my.ascat.obj$data$Tumor_BAF.unscaled = data.frame(sample = ao.df$BAF.unscaled, row.names = ao.df$ProbeSetName)
  colnames(my.ascat.obj$data$Tumor_BAF.unscaled) <- samplename
  my.ascat.obj$data$Tumor_BAF.unisomy = data.frame(sample = ao.df$uni, row.names = ao.df$ProbeSetName)
  colnames(my.ascat.obj$data$Tumor_BAF.unisomy) <- samplename
  rownames(my.ascat.obj$germline$germlinegenotypes) <- ao.df$ProbeSetName
  genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  rm(ao.df)
  gc()
  
  ## Adding meta
  my.ascat.obj$meta = list(
    basic = meta.b,
    affy = affy.meta
  )
  
  ## Adding CEL intensities
  my.ascat.obj$CEL = list(
    CEL1 = affxparser::readCel(filename = CEL)
  )
  my.ascat.obj$CEL$CEL1$intensities <- as.integer(my.ascat.obj$CEL$CEL1$intensities)
  
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
  gc()
  if(return.data) return(my.ascat.obj)
}

SNP6.Process.Batch <- function(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the CEL.list.file
  if (is.null(CEL.list.file)) stop("A CEL.list.file is required !", call. = FALSE)
  if (!file.exists(CEL.list.file)) stop("Could not find CEL.list.file !", call. = FALSE)
  message("Reading and checking CEL.list.file ...")
  myCELs <- read.table(file = CEL.list.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("cel_files", "SampleName")
  head.chk <- all(colnames(CEL.list.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in CEL.list.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(myCELs)))
    stop("Invalid header.", call. = FALSE)
  }
  sn.chk <- duplicated(myCELs$SampleName)
  if (any(sn.chk)) {
    message("CEL.list.file contains duplicated SampleNames !")
    message(myCELs$SampleName[which(duplicated(myCELs$SampleName))])
    stop("Duplicated SampleNames.", call. = FALSE)
  }
  fecheck <- !vapply(myCELs$cel_files, file.exists, TRUE)
  fecheck.pos <- which(fecheck)
  if (length(fecheck.pos) > 0) stop(paste0("\n", "CEL file could not be found : ", myCELs$cel_files[fecheck.pos], collapse = ""), call. = FALSE)

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
  s6res <- foreach::foreach(p = seq_len(nrow(myCELs)), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    SNP6.Process(CEL = myCELs$cel_files[p], samplename = myCELs$SampleName[p], return.data = FALSE, ...)
  }

  ## Stopping cluster
  message("Stopping cluster ...")
  parallel::stopCluster(cl)

  message("Done.")
}

