## Performs OS CEL pairs processing
EaCoN.OS.Process <- function(ATChannelCel = NULL, GCChannelCel = NULL, samplename = NULL, dual.norm = TRUE, l2r.level = "weighted", wave.renorm = TRUE, wave.renorm.rda = NULL, gc.renorm = TRUE, gc.renorm.rda = NULL, out.dir = getwd(), oschp.keep = TRUE, force.OS = NULL, apt.version = "2.4.0", apt.build = "na33.r2", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", return.data = FALSE) {

  # ## TEMP
  # setwd("/home/job/svn/genomics/CGH/R/00_PIPELINE/TEST_ZONE/OS")
  # ATChannelCel = "17P00005_A.CEL"
  # GCChannelCel = "17P00005_C.CEL"
  # samplename = "17P00005"
  # dual.norm = TRUE
  # out.dir = getwd()
  # temp.files.keep = FALSE
  # force.OS = NULL
  # wave.renorm = TRUE
  # wave.renorm.rda <- NULL
  # gc.renorm = TRUE
  # gc.renorm.rda <- NULL
  # l2r.level <- "normal"
  # temp.files.keep = FALSE
  # apt.build = "na33.r2"
  # apt.version = "2.4.0"
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/renorm_functions.R")

  ## Early checks
  if (is.null(ATChannelCel)) stop(tmsg("An ATChannel CEL file is required !"))
  if (is.null(GCChannelCel)) stop(tmsg("A GCChannel CEL file is required !"))
  if (is.null(samplename)) stop(tmsg("A samplename is required !"))
  if (gc.renorm) { if (!is.null(gc.renorm.rda)) { if (!file.exists(gc.renorm.rda)) stop(tmsg(paste0("Could not find gc.renorm.rda file ", gc.renorm.rda))) } }
  if (wave.renorm) { if (!is.null(wave.renorm.rda)) { if (!file.exists(wave.renorm.rda)) stop(tmsg(paste0("Could not find wave.renorm.rda file ", wave.renorm.rda))) } }
  if(!file.exists(ATChannelCel)) stop(tmsg(paste0("Could not find ATChannelCel file ", ATChannelCel, " !")))
  if(!file.exists(GCChannelCel)) stop(paste0("Could not find GCChannelCel file ", GCChannelCel, " !"))
  if (ATChannelCel == GCChannelCel) stop(tmsg("ATChannelCel and GCChannelCel files are identical !"))

  sup.array <- c("OncoScan", "OncoScan_CNV")
  arraytype.celA = affxparser::readCelHeader(filename = ATChannelCel)$chiptype
  arraytype.celC = affxparser::readCelHeader(filename = GCChannelCel)$chiptype
  if (!arraytype.celA %in% sup.array) stop(tmsg(paste0("Identified array type for ATChannelCel file '", arraytype.celA, "' is not supported by this function !")))
  if (!arraytype.celC %in% sup.array) stop(tmsg(paste0("Identified array type for GCChannelCel file '", arraytype.celC, "' is not supported by this function !")))
  if (arraytype.celA != arraytype.celC) stop(tmsg(paste0("ATChannelCel and GCChannelCel files are not of the same design : ", arraytype.celA, " != ", arraytype.celC, " !")))

  ## Checking APT version compatibility
  valid.apt.versions <- c("2.4.0")
  if (!(apt.version %in% valid.apt.versions)) message(tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !")))

  ## Checking build compatibility
  valid.builds <- c("na33.r2", "na33.r4", "na36.r1")
  if (!(tolower(apt.build) %in% valid.builds)) message(tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !")))

  ## Checking apt-copynumber-onco-ssa package loc
  apt.onco.pkg.name <- paste0("apt.oncoscan.", apt.version)
  if (!(apt.onco.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", apt.onco.pkg.name, " not found !")))
  suppressPackageStartupMessages(require(package = apt.onco.pkg.name, character.only = TRUE))

  ## Processing CELs to an OSCHP file
  oscf <- apt.oncoscan.process(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel, samplename = samplename, dual.norm = dual.norm, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)

  ## Reading OSCHP
  my.oschp <- oschp.load(file = oscf)
  if (length(names(my.oschp$Meta)) == 3) names(my.oschp$Meta) <- c("ATChannelCel", "GCChannelCel", "analysis")
  sex.chr <- c("chrX", "chrY")


  ## Processing : meta (and checks)
  if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA
  
  ### Loading genome info
  message(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome2 <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)

  ### Getting basic meta
  # genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  if (genome != genome2) stop(tmsg(paste0("Genome build name given with BSgenome package '", genome.pkg, "', (", genome2, ") is different from the genome build specified by provided APT build version '", apt.build, "' (", genome, ") !")))
  # data(list = genome, package = "chromosomes", envir = environment())
  arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
  manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
  species <- getmeta("affymetrix-algorithm-param-genome-species", my.oschp$Meta$analysis)
  
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = "NA")
  pgender <- gender.conv[[(getmeta("affymetrix-chipsummary-Y-gender-call", my.oschp$Meta$analysis))]]

  if (!(arraytype %in% sup.array)) stop(tmsg(paste0("Unsupported array : '", arraytype, "' !")))

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
  ao.df <- if (l2r.level == "normal") {
    data.frame(chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R = as.vector(my.oschp$ProbeSets$CopyNumber$Log2Ratio), BAF = as.vector(my.oschp$ProbeSets$AllelicData$BAF), stringsAsFactors = FALSE)
  } else if (l2r.level == "weighted") {
    data.frame(chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), L2R = as.vector(my.oschp$ProbeSets$CopyNumber$WeightedLog2Ratio), BAF = as.vector(my.oschp$ProbeSets$AllelicData$BAF), stringsAsFactors = FALSE)
  } else stop(tmsg("Unrecognized value for [l2r.level] !"))
  rownames(ao.df) <- my.oschp$ProbeSets$CopyNumber$ProbeSetName
  affy.chrom <- my.oschp$Chromosomes$Summary
  ak <- affy.chrom$Display
  names(ak) <- affy.chrom$Chromosome
  ao.df$chrA <- as.vector(ak[as.character(ao.df$chrs)])
  ao.df$chr <- paste0("chr", ao.df$chrA)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos, rownames(ao.df)),]
  ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)),]
  ao.df$L2R.ori <- ao.df$L2R

  ## L2R renormalizations
  # smo <- round(nrow(ao.df) / 550)
  # if(smo%%2 == 0) smo <- smo+1
  # if (renormalize) {
  #   message(tmsg("Re-normalization ..."))
  #   if (!is.null(renorm.rda)) {
  #     load(renorm.rda, envir = environment())
  #     gcd.arraytype <- GC.data$info$value[GC.data$info$key == "array_type"]
  #     gcd.genome <- GC.data$info$value[GC.data$info$key == "genome-version"]
  #     if ((gcd.arraytype != arraytype) | (gcd.genome != genome)) stop(tmsg(paste0("Provided renorm.rda is not compatible with currently analyzed sample : expected [", arraytype, ", ", genome, "], got [", gcd.arraytype, ", ", gcd.genome, "] !")))
  #   } else {
  #     GC.pkg.name <- "affy.CN.norm.data"
  #     if (!(GC.pkg.name %in% installed.packages())) stop(paste0("Package ", GC.pkg.name, " not found !"))
  #     GC.file <- system.file(paste0("data/", arraytype, ".", genome, ".GC.rda"), package = GC.pkg.name)
  #     if (GC.file == "") stop(tmsg(paste0("Could not find a GC data package for [", arraytype, ", ", genome, "] in package '", GC.pkg.name, "' ! Please build your own GC data pack with ", GC.pkg.name, "::affy.gc.compute() and submit it using the 'GC.rda' option.")))
  #     data(list = paste0(arraytype, ".", genome, ".GC"), package = GC.pkg.name, envir = environment())
  #   }
  #   gcdata <- GC.data$GC[GC.data$GC$ProbeSetName %in% rownames(ao.df),]
  #   ao.df <- ao.df[rownames(ao.df) %in% gcdata$ProbeSetName,]
  #   if (!all(unique(gcdata$ProbeSetName == rownames(ao.df)))) stop(tmsg("GC data and L2R data are not synched, or ordered differently !"))
  #   ndata <- data.frame(chr = paste0("chr", ao.df$chrs), start = ao.df$pos, end = ao.df$pos, name = rownames(ao.df), gcdata[,-c(1:4)], stringsAsFactors = FALSE)
  #   # ndata <- cbind(ndata, gcdata[,-c(1:4)])
  #   rm(gcdata, GC.data)
  #   my.rm.mad <- sum(abs(diff(as.numeric(runmed(ao.df$L2R[!is.na(ao.df$L2R)], smo)))))
  #   normloop.res <- l2r.fitloop(l2rObj = list(l2r=ao.df$L2R, rm.mad = my.rm.mad), tfd = ndata, smo = smo)
  #   ao.df$L2R <- normloop.res$l2r$l2r + median(ao.df$L2R, na.rm = TRUE)
  #   # if(is.null(normloop.res$pos)) meta.df <- rbind(meta.df, c("a2p-renorm", "None")) else meta.df <- rbind(meta.df, c("a2p-renorm", paste0(normloop.res$pos, collapse = ",")))
  #   if(is.null(normloop.res$pos)) meta.b <- setmeta("renorm", "None", meta.b) else meta.b <- setmeta("renorm", paste0(normloop.res$pos, collapse = ","), meta.b)
  # } else {
  #   meta.b <- setmeta("renorm", "FALSE", meta.b)
  # }

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
  message(tmsg("Building normalized object ..."))
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
      sexchromosomes = c("X", "Y"),
      failedarrays = NULL
    ),
    meta = list(
      basic = meta.b,
      affy = my.oschp$Meta
    ),
    CEL = list(
      ATChannelCel = affxparser::readCel(filename = ATChannelCel),
      GCChannelCel = affxparser::readCel(filename = GCChannelCel)
    )
  )
  colnames(my.ascat.obj$data$Tumor_LogR) <- samplename
  colnames(my.ascat.obj$data$Tumor_BAF) <- samplename

  saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "xz")

  ## Rough plot
  message(tmsg("Plotting ..."))
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
  plot(ao.df$genopos, ao.df$L2R, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " normalized L2R profile (median-centered) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$BAF, pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " ", arraytype, " BAF profile"), xlab = "Genomic position", ylab = "BAF")
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = .5, col = 2, lty = 2, lwd = 2)
  dev.off()

  ## Cleaning
  if(!oschp.keep) {
    message(tmsg("Removing temporary OSCHP file ..."))
    file.remove(oscf)
  }

  message(tmsg("Done."))
  gc()
  if(return.data) return(my.ascat.obj)
}

## Same in batch mode
EaCoN.OS.Process.Batch <- function(pairs.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the pairs.file
  if (!file.exists(pairs.file)) stop("Could not find pairs.file !")
  message("Reading and checking pairs.file ...")
  mypairs <- read.table(file = pairs.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("ATChannelCel", "GCChannelCel", "SampleName")
  head.chk <- all(colnames(pairs.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in pairs.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(mypairs)))
    stop("Invalid header.")
  }

  at.chk <- file.exists(mypairs$ATChannelCel)
  gc.chk <- file.exists(mypairs$GCChannelCel)

  if (!all(at.chk) || !all(at.chk)) {
    message("Some CEL files from the pairs.file could not be found (wrong path ?) !")
    message("Missing AT CEL files :")
    message(paste0(mypairs$ATChannelCel[which(!at.chk)], collapse = ", "))
    message("Missing GC CEL files :")
    message(paste0(mypairs$GCChannelCel[which(!gc.chk)], collapse = ", "))
    stop("Missing CEL.")
  }
  sn.chk <- duplicated(mypairs$SampleName)
  if (any(sn.chk)) {
    message("pairs.file contains duplicated SampleNames !")
    message(mypairs$SampleName[which(duplicated(mypairs$SampleName))])
    stop("Duplicated SampleNames.")
  }
  if(any(mypairs$ATChannelCel == mypairs$GCChannelCel)) {
    message("Some ATChannelCel and GCChannelCel are the same for at least one sample !")
    stop("Identical CEL files for AT and GC.")
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
    EaCoN.OS.Process(ATChannelCel = mypairs$ATChannelCel[p], GCChannelCel = mypairs$GCChannelCel[p], samplename = mypairs$SampleName[p], ...)
  }
  message("Stopping cluster ...")
  parallel::stopCluster(cl)

  message("Done.")
}

