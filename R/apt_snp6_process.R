## Performs CS CEL processing
EaCoN.SNP6.Process <- function(CEL = NULL, samplename = NULL, l2r.level = "normal", wave.renorm = TRUE, wave.renorm.rda = NULL, gc.renorm = TRUE, gc.renorm.rda = NULL, out.dir = getwd(), oschp.keep = TRUE, force.OS = NULL, apt.version = "1.20.0", apt.build = "na35.r1", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE) {

  # setwd("/home/job/WORKSPACE/MP/SNP6")
  # CEL <- "/home/job/WORKSPACE/MP/SNP6/AMAZE_p_TCGASNP_b86_87_88_N_GenomeWideSNP_6_H02_735486.CEL"
  # samplename <- "AMAZE_p_TCGASNP_b86_87_88_N_GenomeWideSNP_6_H02_735486"
  # l2r.level <- "normal"
  # wave.renorm <- TRUE
  # wave.renorm.rda <- NULL
  # gc.renorm <- TRUE
  # gc.renorm.rda <- NULL
  # out.dir <- getwd()
  # oschp.keep <- TRUE
  # force.OS <- NULL
  # apt.version <- "1.20.0"
  # apt.build <- "na35.r1"
  # genome.pkg <- "BSgenome.Hsapiens.UCSC.hg19"
  # return.data <- FALSE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/renorm_functions.R")


  ## Early checks
  if (is.null(CEL)) stop(tmsg("A CEL file is required !"))
  if (is.null(samplename)) stop(tmsg("A samplename is required !"))
  if (!file.exists(CEL)) stop(tmsg(paste0("Could not find CEL file ", CEL, " !")))
  if (gc.renorm) { if (!is.null(gc.renorm.rda)) { if (!file.exists(gc.renorm.rda)) stop(tmsg(paste0("Could not find gc.renorm.rda file ", gc.renorm.rda))) } }
  if (wave.renorm) { if (!is.null(wave.renorm.rda)) { if (!file.exists(wave.renorm.rda)) stop(tmsg(paste0("Could not find wave.renorm.rda file ", wave.renorm.rda))) } }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"))
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }

  sup.array <- c("GenomeWideSNP_6")
  arraytype.cel = affxparser::readCelHeader(filename = CEL)$chiptype
  if (!arraytype.cel %in% sup.array) stop(tmsg(paste0("Identified array type '", arraytype.cel, "' is not supported by this function !")))

  ## Checking APT version compatibility
  valid.apt.versions <- c("1.20.0")
  if (!(apt.version %in% valid.apt.versions)) warning(tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !")))

  ## Checking build compatibility
  valid.builds <- c("na35.r1")
  if (!(tolower(apt.build) %in% valid.builds)) warning(tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !")))

  ## Checking apt-copynumber-cyto-ssa package loc
  apt.snp6.pkg.name <- paste0("apt.snp6.", apt.version)
  if (!(apt.snp6.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", apt.snp6.pkg.name, " not found !")))
  suppressPackageStartupMessages(require(package = apt.snp6.pkg.name, character.only = TRUE))

  ## Processing CEL to an OSCHP file
  oscf <- apt.snp6.process(CEL = CEL, samplename = samplename, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)

  ## Reading OSCHP
  print(tmsg("Loading OSCHP file ..."))
  my.oschp <- oschp.load(file = oscf)
  sex.chr <- c("chrX", "chrY")
  
  ### Loading genome info
  message(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome2 <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)
  
  ## Processing : meta (and checks)
  if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA

  ### Getting basic meta
  genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  if (genome != genome2) stop(tmsg(paste0("Genome build name given with BSgenome package '", genome.pkg, "', (", genome2, ") is different from the genome build specified by provided APT build version '", apt.build, "' (", genome, ") !")))
  arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
  manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
  species <- getmeta("affymetrix-algorithm-param-genome-species", my.oschp$Meta$analysis)
  
  snp6.conv <- list("1" = "male", "2" = "female", "NA" = "NA", "0" = "NA")
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = "NA")
  pgender <- gender.conv[[snp6.conv[[as.character(as.numeric(getmeta("affymetrix-chipsummary-Gender", my.oschp$Meta$analysis)))]]]]

  if (!(arraytype %in% sup.array)) stop(paste0("Unsupported array : '", arraytype, "' !"))

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
  if (l2r.level == "normal") {
    ao.df <- data.frame(chrs = as.vector(my.oschp$MultiData$CopyNumber$Chromosome), pos = as.vector(my.oschp$MultiData$CopyNumber$Position), L2R.ori = as.vector(my.oschp$MultiData$CopyNumber$Log2Ratio), L2R = as.vector(my.oschp$MultiData$CopyNumber$Log2Ratio), BAF = NA, stringsAsFactors = FALSE)
  } else if (l2r.level == "weighted") {
    ao.df <- data.frame(chrs = as.vector(my.oschp$MultiData$CopyNumber$Chromosome), pos = as.vector(my.oschp$MultiData$CopyNumber$Position), L2R.ori = as.vector(my.oschp$MultiData$CopyNumber$SmoothSignal), L2R = as.vector(my.oschp$MultiData$CopyNumber$SmoothSignal), BAF = NA, stringsAsFactors = FALSE)
  } else stop(tmsg("Unrecognized value for [l2r.level] !"))
  rownames(ao.df) <- my.oschp$MultiData$CopyNumber$ProbeSetName
  affy.chrom <- my.oschp$MultiData[["CopyNumber_&keyvals"]][seq.int(3, nrow(my.oschp$MultiData[["CopyNumber_&keyvals"]]), 3),1:2]
  ak <- affy.chrom$val
  names(ak) <- as.numeric(sub(":display", "", affy.chrom$key))
  ao.df$chrA <- as.vector(ak[as.character(ao.df$chrs)])
  ao.df$chr <- paste0("chr", ao.df$chrA)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos),]

  ## Normalizing SNPs
  print(tmsg("Normalizing SNP data (using rcnorm) ..."))
  baf.df <- rcnorm::rcnorm.snp(myCEL = CEL, genome = genome, allSNPs = FALSE)
  baf.df$chr <- paste0("chr", baf.df$chrs)
  baf.df$chrN <- unlist(cs$chrom2chr[baf.df$chr])
  baf.df <- baf.df[order(baf.df$chrN, baf.df$pos),]
  baf.df <- baf.df[!is.na(baf.df$BAF),]
  gc()

  # ## Merging L2R and BAF data
  ao.df <- ao.df[order(rownames(ao.df)),]
  baf.df <- baf.df[order(rownames(baf.df)),]
  bina <- which(rownames(baf.df) %in% rownames(ao.df))
  ainb <- which(rownames(ao.df) %in% rownames(baf.df))
  if(!all(rownames(ao.df)[ainb] == rownames(baf.df)[bina])) stop(tmsg("Could not synch L2R and (rcnorm-processed) BAF data !"))
  ao.df$BAF[ainb] <- baf.df$BAF[bina]
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos, rownames(ao.df)),]

  # ao.df <- ao.df[!duplicated(ao.df$pos),]
  ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)),]

  ## L2R renormalizations
  smo <- round(nrow(ao.df) / 550)
  if(smo%%2 == 0) smo <- smo+1
  
  ### Wave  
  if (wave.renorm) {
    message(tmsg("Wave re-normalization ..."))
    
    ren.res <- renorm.go(input.data = ao.df, renorm.rda = wave.renorm.rda, track.type = "Wave", smo = smo, arraytype = arraytype, genome = genome)
    
    ao.df <- ren.res$data
    fitted.l2r <- ren.res$renorm$l2r$l2r
    
    if(is.null(ren.res$renorm$pos)) {
      meta.b <- setmeta("wave.renorm", "None", meta.b)
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
    message(tmsg("GC renormalization ..."))
    
    ren.res <- renorm.go(input.data = ao.df, renorm.rda = gc.renorm.rda, track.type = "GC", smo = smo, arraytype = arraytype, genome = genome)
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
  
  ## Building ASCAT-like object
  print(tmsg("Building normalized object ..."))
  my.ch <- sapply(unique(ao.df$chrs), function(x) { which(ao.df$chrs == x) })
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR = data.frame(sample = ao.df$L2R, row.names = rownames(ao.df)),
      Tumor_BAF = data.frame(sample = ao.df$BAF, row.names = rownames(ao.df)),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      # SNPpos = data.frame(chrs = ao.df$chrA, pos = ao.df$pos, row.names = rownames(ao.df)),
      SNPpos = data.frame(chrs = ao.df$chr, pos = ao.df$pos, row.names = rownames(ao.df)),
      ch = my.ch,
      chr = my.ch,
      # chrs = unique(ao.df$chrA),
      chrs = unique(ao.df$chr),
      samples = samplename,
      gender = as.vector(meta.b$predicted.gender),
      sexchromosomes = sex.chr,
      failedarrays = NULL
    ),
    meta = list(
      basic = meta.b,
      affy = my.oschp$Meta
    ),
    CEL = list(
      CEL1 = affxparser::readCel(filename = CEL)
    )
  )
  colnames(my.ascat.obj$data$Tumor_LogR) <- samplename
  colnames(my.ascat.obj$data$Tumor_BAF) <- samplename

  saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "xz")

  ## Rough plot
  print(tmsg("Plotting ..."))
  ao.df$genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  ao.df$L2R.ori <- ao.df$L2R.ori - median(ao.df$L2R.ori, na.rm = TRUE)
  kend <- ao.df$genopos[vapply(unique(ao.df$chrN), function(k) { max(which(ao.df$chrN == k))}, 1)]
  l2r.notna <- which(!is.na(ao.df$L2R))
  l2r.rm <- runmed(ao.df$L2R[l2r.notna], smo)
  l2r.ori.rm <- runmed(ao.df$L2R.ori[l2r.notna], smo)
  png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
  # par(mfrow = c(2,1))
  par(mfrow = c(3,1))
  plot(ao.df$genopos, ao.df$L2R.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " raw L2R profile (median-centered) / ", round(sum(abs(diff(l2r.ori.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.ori.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$L2R, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " L2R profile (median-centered) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$BAF, pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = .5, col = 2, lty = 2, lwd = 2)
  dev.off()

  ## Cleaning
  if(!oschp.keep) {
    print(tmsg("Removing temporary OSCHP file ..."))
    file.remove(oscf)
  }

  print(tmsg("Done."))
  gc()
  if(return.data) return(my.ascat.obj)
}

EaCoN.SNP6.Process.Batch <- function(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
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
  s6res <- foreach::foreach(p = seq_len(nrow(myCELs)), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    EaCoN.SNP6.Process(CEL = myCELs$cel_files[p], samplename = myCELs$SampleName[p], return.data = FALSE, ...)
  }

  ## Stopping cluster
  message("Stopping cluster ...")
  parallel::stopCluster(cl)

  message("Done.")
}

