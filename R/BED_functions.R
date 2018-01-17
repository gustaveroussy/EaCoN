EaCoN.BedGC.fasta <- function(binned.bed.file = NULL, genome = "hg19", fasta = NULL, na.to0 = TRUE, nt.add = c(0,50,100, 250, 500, 1000, 2500, 5000), out.dir = getwd(), return.data = FALSE, nthread = 1) {

  print("WARNING ! On hg19, this function consumes 3 GB per thread (to add to your current R session) !")
  print("When possible, prefer the BedGC.chr() function that performs the same task using sequences as single-chromosome FASTA files. This consumes only 0.3 GB of RAM per thread.")
  if (is.null(binned.bed.file)) stop("A BED file is required !")
  if (!file.exists(binned.bed.file)) stop("Could not find the BED file !")
  if (is.null(fasta)) stop("A FASTA file is required !")
  if (!file.exists(fasta)) stop("Could not find FASTA file !")
  if (!all(is.numeric(nt.add))) stop("nt.add must be a numeric vector !")
  if (is.null(out.dir)) print("NOTE : Checked / cleaned bed will be written in the same directory as source.")
  if (!file.exists(out.dir)) stop("Could not find the output directory !")
  if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !")

  print("Loading genome information ...")
  data(list = genome, package = "chromosomes", envir = environment())

  print("Loading BED ...")
  bed.data <- read.table(file = binned.bed.file, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(bed.data) < 3) stop("BED file must contain at least 3 columns !")
  bed.data <- bed.data[,1:3]
  colnames(bed.data) <- c("chr", "start", "end")
  if (!all(unique(bed.data$chr) %in% cs$chromosomes$chrom)) stop(paste0("BED file contains chromosome(s) not defined in the ", genome, " genome !"))

  klen <- Biostrings::fasta.seqlengths(filepath = fasta)
  kdata <- Biostrings::readDNAStringSet(filepath = fasta, format = "fasta")

  print("Starting cluster ...")
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK")
  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::"%dopar%"
  a <- 0
  Agcpc <- foreach::foreach(a = nt.add, .combine = "cbind") %dopar% {
    print(paste0("Computing GC for frame +", a, " bp"))

    querystart <- bed.data$start - a
    queryend <- bed.data$end + a
    querystart[querystart < 1] <- 1
    kmore.idx <- which(queryend > klen[bed.data$chr])
    if (length(kmore.idx) > 0) queryend[kmore.idx] <- klen[names(kmore.idx)]

    kbs.idx <- vapply(bed.data$chr, function(k) { which(kdata@ranges@NAMES == k) }, 1L)

    gcpc <- vapply(1:nrow(bed.data), function(x) {
      aFreq <- Biostrings::alphabetFrequency(Biostrings::subseq(kdata[kbs.idx[x]], querystart[x], queryend[x]), baseOnly = TRUE)
      acgt.count <- sum(aFreq[colnames(aFreq) != "other"])
      gc.count <- sum(aFreq[colnames(aFreq) %in% c("G", "C")])
      gcpc.x <- gc.count / acgt.count
      return(gcpc.x)
    }, .1)
    if (na.to0) gcpc[is.na(gcpc)] <- 0
    return(gcpc)
  }

  print("Stopping cluster ...")
  parallel::stopCluster(cl)

  ## Building the output object
  colnames(Agcpc) <- paste0("GC", nt.add, "b")
  out.df <- cbind(bed.data, Agcpc)
  outname <- paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = ".gc", x = basename(binned.bed.file), ignore.case = TRUE))
  write.table.fast(x = out.df, file = outname)
  print("Done.")
  if(return.data) return(out.df) else return(NULL)
}

EaCoN.BedGC.fasta.chr <- function(binned.bed.file = NULL, genome = "hg19", fasta.dir = NULL, na.to0 = TRUE, nt.add = c(0,50,100, 250, 500, 1000, 2500, 5000), out.dir = getwd(), return.data = FALSE, nthread = 2) {

  if (is.null(binned.bed.file)) stop("A BED file is required !")
  if (!file.exists(binned.bed.file)) stop("Could not find the BED file !")
  if (is.null(fasta.dir)) stop("A directory containing indexed FASTA file(s) for each chromosome is required !")
  if (!file.exists(fasta.dir)) stop("Could not find the FASTA directory !")
  if (!file.info(fasta.dir)$isdir) stop("fasta.dir is not a directory !")
  if (!all(is.numeric(nt.add))) stop("nt.add must be a numeric vector !")
  if (is.null(out.dir)) stop("An output directory is required !")
  if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !")

  print("Loading genome information ...")
  data(list = genome, package = "chromosomes", envir = environment())

  print("Loading BED ...")
  bed.data <- read.table(file = binned.bed.file, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(bed.data) < 3) stop("BED file must contain at least 3 columns !")
  bed.data <- bed.data[,1:3]
  colnames(bed.data) <- c("chr", "start", "end")
  if (!all(unique(bed.data$chr) %in% cs$chromosomes$chrom)) stop(paste0("BED file contains chromosome(s) not defined in the ", genome, " genome !"))

  kcoords <- sapply(unique(bed.data$chr), function(k) { return(bed.data[bed.data$chr == k, ]) }, simplify = FALSE)

  print("Starting cluster ...")
  # require(foreach)
  # require(doParallel)
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK")
  doParallel::registerDoParallel(cl)

  a <- 0
  `%do%` <- foreach::"%do%"
  `%dopar%` <- foreach::"%dopar%"

  Agcpc <- foreach(a=nt.add, .combine = "cbind") %do% {
    print(paste0("Computing GC for frame +", a, " bp"))

    k <- 0
    gcpc <- foreach(k=kcoords, .combine = "c", .noexport = c("kcoords", "bed.data")) %dopar% {
      kchr <- unique(k$chr)
      kfafile <- paste0(fasta.dir, "/", kchr, ".fa")
      if (!file.exists(kfafile)) stop(paste0("Could ont find : ", kfafile))
      cat("Loading sequence for", kchr, "...\n")
      # require(Biostrings)
      klen <- Biostrings::fasta.seqlengths(filepath = kfafile)
      kdata <- Biostrings::readDNAStringSet(filepath = kfafile, format = "fasta")

      querystart <- k$start - a
      queryend <- k$end + a
      querystart[querystart < 1] <- 1
      queryend[queryend > klen] <- klen

      gcpc.k <- vapply(1:nrow(k), function(x) {
        aFreq <- Biostrings::alphabetFrequency(Biostrings::subseq(kdata, querystart[x], queryend[x]), baseOnly = TRUE)
        acgt.count <- sum(aFreq[colnames(aFreq) != "other"])
        gc.count <- sum(aFreq[colnames(aFreq) %in% c("G", "C")])
        gcpc.x <- gc.count / acgt.count
        return(gcpc.x)
      }, .1)
      return(gcpc.k)
    }
    if (na.to0) gcpc[is.na(gcpc)] <- 0
    return(gcpc)
  }
  print("Stopping cluster ...")
  parallel::stopCluster(cl)
  ## Building the output object
  # colnames(Agcpc) <- paste0("GC", nt.add, "b")
  colnames(Agcpc) <- paste0("GC", nt.add, "b")
  out.df <- cbind(bed.data, Agcpc)
  # outname <- paste0(out.dir,  "/", sub(pattern = "\\.bed$", replacement = ".gc", x = basename(binned.bed.file), ignore.case = TRUE))
  outname <- paste0(out.dir,  "/", sub(pattern = "\\.bed$", replacement = "_GC.RDS", x = basename(binned.bed.file), ignore.case = TRUE))
  saveRDS(out.df, file = outname, compress = "bzip2")
  if(return.data) return(out.df) else return(NULL)
  print("Done.")
  return(cbind(bed.data, Agcpc))
}

EaCoN.BedGC.R <- function(binned.bed.file = NULL, human.genome.build = "hg19", na.to0 = TRUE, nt.add = c(0,50,100, 250, 500, 1000, 2500, 5000), out.dir = getwd(), return.data = FALSE, nthread = 1) {

  print("WARNING ! On hg19, this function consumes 3 GB per thread (to add to your current R session) !")
  print("When possible, prefer the BedGC.chr() function that performs the same task using sequences as single-chromosome FASTA files. This consumes only 0.3 GB of RAM per thread.")
  if (is.null(binned.bed.file)) stop("A BED file is required !")
  if (is.null(human.genome.build)) stop("A genome name is required !")
  if (!file.exists(binned.bed.file)) stop("Could not find the BED file !")
  # if (is.null(fasta)) stop("A FASTA file is required !")
  # if (!file.exists(fasta)) stop("Could not find FASTA file !")
  if (!all(is.numeric(nt.add))) stop("nt.add must be a numeric vector !")
  if (is.null(out.dir)) print("NOTE : Checked / cleaned bed will be written in the same directory as source.")
  if (!file.exists(out.dir)) stop("Could not find the output directory !")
  if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !")

  print("Loading genome information ...")
  data(list = human.genome.build, package = "chromosomes", envir = environment())

  if (human.genome.build == "hg19") genome.package <- "BSgenome.Hsapiens.UCSC.hg19"
  if (human.genome.build == "hg38") genome.package <- "BSgenome.Hsapiens.UCSC.hg38"

  print("Loading BED ...")
  bed.data <- read.table(file = binned.bed.file, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(bed.data) < 3) stop("BED file must contain at least 3 columns !")
  bed.data <- bed.data[,1:3]
  colnames(bed.data) <- c("chr", "start", "end")
  if (!all(unique(bed.data$chr) %in% cs$chromosomes$chrom)) stop(paste0("BED file contains chromosome(s) not defined in the ", genome, " genome !"))

  require(package = genome.package, character.only = TRUE)
  organisms <- ls(paste0("package:",genome.package))
  orga.id <- organisms[length(organisms)]

  klen <- GenomeInfoDb::seqlengths(Hsapiens)

  print("Starting cluster ...")
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK")
  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::"%dopar%"
  a <- 0
  Agcpc <- foreach::foreach(a = nt.add, .combine = "cbind") %dopar% {
    print(paste0("Computing GC for frame +", a, " bp"))

    Kgcpc <- foreach::foreach(k = unique(bed.data$chr), .combine = "rbind") %dopar% {
      kbed <- bed.data[bed.data$chr == k,]
      kbed$start <- kbed$start - a
      kbed$end <- kbed$end + a
      kbed$start[kbed$start < 1] <- 1
      kmore.idx <- which(kbed$end > klen[k])
      if (length(kmore.idx) > 0) kbed$end[kmore.idx] <- klen[k]

      kseq <- Hsapiens[[k]]

      gcpc <- vapply(1:nrow(kbed), function(x) {
        aFreq <- Biostrings::alphabetFrequency(Biostrings::subseq(kseq, kbed$start[x], kbed$end[x]), baseOnly = TRUE)
        acgt.count <- sum(aFreq[names(aFreq) != "other"])
        gc.count <- sum(aFreq[names(aFreq) %in% c("G", "C")])
        gcpc.x <- gc.count / acgt.count
        return(gcpc.x)
      }, .1)
      if (na.to0) gcpc[is.na(gcpc)] <- 0
      return(gcpc)
    }
    return(Kgcpc)
  }

  print("Stopping cluster ...")
  parallel::stopCluster(cl)

  ## Building the output object
  colnames(Agcpc) <- paste0("GC", nt.add, "b")
  out.df <- cbind(bed.data, Agcpc)
  outname <- paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = ".gc", x = basename(binned.bed.file), ignore.case = TRUE))
  write.table.fast(x = out.df, file = outname)
  print("Done.")
  if(return.data) return(out.df) else return(NULL)
}

EaCoN.BedCheck <- function(bed.file = NULL, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", out.dir = getwd(), return.data = FALSE) {

  if (is.null(bed.file)) stop("A BED file is required !")
  if (!file.exists(bed.file)) stop("Could not find the BED file !")
  if (is.null(out.dir)) print("NOTE : Checked / cleaned bed will be written in the same directory as source.") else {
    if (!file.exists(out.dir)) stop("Could not find the output directory !")
    if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !")
  }
  print("Checking BED ...")
  # data(list = genome, package = "chromosomes", envir = environment())
  print(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)

  bed.data <- read.table(file = bed.file, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(bed.data) < 3) stop("BED file must contain at least 3 columns !")
  bed.data <- bed.data[,1:3]

  colnames(bed.data) <- c("chr", "start", "end")

  if (!all(unique(bed.data$chr) %in% BSgenome::seqnames(BSg.obj))) stop(paste0("BED file contains chromosome(s) not defined in the ", genome.pkg, " genome !"))

  myGR <- GenomicRanges::GRanges(bed.data)

  print(" Removing inproper coordinates (start >= end) ...")
  myGR <- myGR[!GenomicRanges::width(myGR) < 1,]

  print(" Removing replicates ...")
  myGR <- myGR[!GenomicRanges::duplicated(myGR),]

  print(" Looking for overlaps ...")
  myGR <- GenomicRanges::reduce(myGR)

  print(" Sorting ...")
  myGR <- GenomicRanges::sort(myGR)

  ## Converting GR back to df
  bed.data <- as.data.frame(myGR)
  bed.data$seqnames <- as.character(bed.data$seqnames)
  bed.data <- bed.data[,1:3]

  if (!return.data) {
    print("Writing clean bed ...")
    colnames(bed.data)[1] <- "#chr"
    outname <- paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = paste0("_", genome, "_merged_sorted.bed"), x = basename(bed.file), ignore.case = TRUE))
    write.table.fast(x = bed.data, file = outname)
  }

  colnames(bed.data)[1] <- "chr"

  print("Done.")
  if(return.data) return(bed.data)
}



