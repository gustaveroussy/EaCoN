## Performs CS CEL processing
EaCoN.SNP6.Process <- function(CEL = NULL, samplename = NULL, l2r.level = "normal", renormalize = TRUE, renorm.rda = NULL, out.dir = getwd(), oschp.keep = TRUE, force.OS = NULL, apt.version = "1.20.0", apt.build = "na35.r1", return.data = FALSE) {

  # setwd("/home/job/svn/genomics/CGH/R/00_PIPELINE/TEST_ZONE/SNP6")
  # CEL <- "BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_C03_772214.CEL"
  # samplename <- "BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_C03_772214"
  # l2r.level <- "normal"
  # renormalize <- TRUE
  # renorm.rda <- NULL
  # out.dir <- getwd()
  # oschp.keep <- TRUE
  # force.OS <- NULL
  # apt.version <- "1.20.0"
  # apt.build <- "na35.r1"
  # return.data <- FALSE


  ## Early checks
  if (is.null(CEL)) stop(tmsg("A CEL file is required !"))
  if (is.null(samplename)) stop(tmsg("A samplename is required !"))
  if (!file.exists(CEL)) stop(tmsg(paste0("Could not find CEL file ", CEL, " !")))
  if (renormalize) { if (!is.null(renorm.rda)) { if (!file.exists(renorm.rda)) stop(tmsg(paste0("Could not find renorm.rda file ", renorm.rda))) } }

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
  require(package = apt.snp6.pkg.name, character.only = TRUE)

  ## Processing CEL to an OSCHP file
  oscf <- apt.snp6.process(CEL = CEL, samplename = samplename, out.dir = out.dir, temp.files.keep = FALSE, force.OS = force.OS, apt.build = apt.build)

  ## Reading OSCHP
  print(tmsg("Loading OSCHP file ..."))
  my.oschp <- oschp.load(file = oscf)

  ## Processing : meta (and checks)
  # meta.a1.df <- my.oschp[["Dset_IO_HDF5_Gdh"]][["_&keyvals"]][,1:2]
  # meta.a1.df <- meta.a1.df[!duplicated(meta.a1.df$key),]
  # meta.a1 <- meta.df2list(meta.a1.df)
  # meta.a1[["L2R.level"]] <- l2r.level
  # if (!("affymetrix-chipsummary-snp-qc" %in% names(meta.a1))) meta.a1[["affymetrix-chipsummary-snp-qc"]] <- NA
  if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA
  # rm(meta.a1.df)

  ### Getting basic meta
  # genome <- getmeta("affymetrix-algorithm-param-genome-version", meta.a1)
  genome <- getmeta("affymetrix-algorithm-param-genome-version", my.oschp$Meta$analysis)
  data(list = genome, package = "chromosomes", envir = environment())
  # arraytype <- getmeta("affymetrix-array-type", meta.a1)
  arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
  # manufacturer <- getmeta("program-company", meta.a1)
  manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
  # species <- getmeta("affymetrix-algorithm-param-genome-species", meta.a1)
  species <- getmeta("affymetrix-algorithm-param-genome-species", my.oschp$Meta$analysis)
  # pgender <- as.character(as.numeric(getmeta("affymetrix-chipsummary-Gender", meta.a1)))

  snp6.conv <- list("1" = "male", "2" = "female", "NA" = NA, "0" = NA)
  gender.conv <- list("female" = "XX", "male" = "XY", "NA" = NA)
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
    predicted.gender = pgender
  )

  ## Extracting data : L2R
  ao.df <- if (l2r.level == "normal") {
    ao.df <- data.frame(chrs = as.vector(my.oschp$MultiData$CopyNumber$Chromosome), pos = as.vector(my.oschp$MultiData$CopyNumber$Position), L2R = as.vector(my.oschp$MultiData$CopyNumber$Log2Ratio), BAF = NA, stringsAsFactors = FALSE)
  } else if (l2r.level == "weighted") {
    ao.df <- data.frame(chrs = as.vector(my.oschp$MultiData$CopyNumber$Chromosome), pos = as.vector(my.oschp$MultiData$CopyNumber$Position), L2R = as.vector(my.oschp$MultiData$CopyNumber$SmoothSignal), BAF = NA, stringsAsFactors = FALSE)
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
  if (renormalize) {
    print(tmsg("Re-normalization ..."))
    if (!is.null(renorm.rda)) {
      load(renorm.rda, envir = environment())
      gcd.arraytype <- GC.data$info$value[GC.data$info$key == "array_type"]
      gcd.genome <- GC.data$info$value[GC.data$info$key == "genome-version"]
      if ((gcd.arraytype != arraytype) | (gcd.genome != genome)) stop(tmsg(paste0("Provided renorm.rda ins not compatible with currently analyzed sample : expected [", arraytype, ", ", genome, "], got [", gcd.arraytype, ", ", gcd.genome, "] !")))
    } else {
      GC.pkg.name <- "affy.CN.norm.data"
      if (!(GC.pkg.name %in% installed.packages())) stop(tmsg(paste0("Package ", GC.pkg.name, " not found !")))
      GC.file <- system.file(paste0("data/", arraytype, ".", genome, ".GC.rda"), package = GC.pkg.name)
      if (GC.file == "") stop(tmsg(paste0("Could not find a GC data package for [", arraytype, ", ", genome, "] in package '", GC.pkg.name, "' ! Please build your own GC data pack with ", GC.pkg.name, "::affy.gc.compute() and submit it using the 'GC.rda' option.")))
      data(list = paste0(arraytype, ".", genome, ".GC"), package = GC.pkg.name, envir = environment())
    }
    gcdata <- GC.data$GC[GC.data$GC$ProbeSetName %in% rownames(ao.df),]
    ao.df <- ao.df[rownames(ao.df) %in% gcdata$ProbeSetName,]
    if (!all(unique(gcdata$ProbeSetName == rownames(ao.df)))) stop(tmsg("GC data and L2R data are not synched, or ordered differently !"))
    ndata <- data.frame(chr = paste0("chr", ao.df$chrs), start = ao.df$pos, end = ao.df$pos, name = rownames(ao.df), gcdata[,-c(1:4)], stringsAsFactors = FALSE)
    # ndata <- cbind(ndata, gcdata[,-c(1:4)])
    rm(gcdata, GC.data)
    my.rm.mad <- sum(abs(diff(as.numeric(runmed(ao.df$L2R[!is.na(ao.df$L2R)], smo)))))
    normloop.res <- l2r.fitloop(l2rObj = list(l2r=ao.df$L2R, rm.mad = my.rm.mad), tfd = ndata, smo = smo)
    ao.df$L2R <- normloop.res$l2r$l2r + median(ao.df$L2R, na.rm = TRUE)
    # if(is.null(normloop.res$pos)) meta.df <- rbind(meta.df, c("a2p-renorm", "None")) else meta.df <- rbind(meta.df, c("a2p-renorm", paste0(normloop.res$pos, collapse = ",")))
    if(is.null(normloop.res$pos)) meta.b <- setmeta("renorm", "None", meta.b) else meta.b <- setmeta("renorm", paste0(normloop.res$pos, collapse = ","), meta.b)
  } else {
    meta.b <- setmeta("renorm", "FALSE", meta.b)
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
      SNPpos = data.frame(chrs = ao.df$chrA, pos = ao.df$pos, row.names = rownames(ao.df)),
      ch = my.ch,
      chr = my.ch,
      chrs = unique(ao.df$chrA),
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
      CEL1 = affxparser::readCel(filename = CEL)
    )
  )
  colnames(my.ascat.obj$data$Tumor_LogR) <- samplename
  colnames(my.ascat.obj$data$Tumor_BAF) <- samplename

  saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "xz")

  ## Rough plot
  print(tmsg("Plotting ..."))
  ao.df$genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  # ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  kend <- ao.df$genopos[vapply(unique(ao.df$chrN), function(k) { max(which(ao.df$chrN == k))}, 1)]
  l2r.notna <- which(!is.na(ao.df$L2R))
  l2r.rm <- runmed(ao.df$L2R[l2r.notna], smo)
  png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
  par(mfrow = c(2,1))
  plot(ao.df$genopos, ao.df$L2R, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " ", arraytype, " L2R profile (median-centered)"), xlab = "Genomic position", ylab = "L2R")
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
  s6res <- foreach::foreach(p = seq_len(nrow(myCELs)), .inorder = FALSE, .errorhandling = "pass") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    EaCoN.SNP6.Process(CEL = myCELs$cel_files[p], samplename = myCELs$SampleName[p], return.data = FALSE, ...)
  }

  ## Stopping cluster
  message("Stopping cluster ...")
  parallel::stopCluster(cl)

  message("Done.")
}

