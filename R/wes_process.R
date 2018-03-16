## Generates a BINpack (bins and pre-computed GC tracks) from a capture bed
EaCoN.BINpack.Maker <- function(bed.file = NULL, bin.size = 50, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", extend.multi = c(0, 50, 100, 200, 400, 800, 1600, 3200, 6400), blocksize = 1E+04, nthread = 1, out.dir = getwd(), return.data = FALSE) {

  ## Checks
  if (is.null(bed.file)) stop("A BED file is required !")
  if (!file.exists(bed.file)) stop("Could not find the BED file !")
  if (is.null(out.dir)) message("NOTE : Checked / cleaned bed will be written in the same directory as source.") else {
    if (!file.exists(out.dir)) stop("Could not find the output directory !")
    if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !")
  }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"))
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }
  

  ### Loading genome
  message(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)


  ## Cleaning BED
  bed.clean <- EaCoN.BedCheck(bed.file = bed.file, genome = genome, out.dir = out.dir, return.data = TRUE)

  ## Binning
  message(paste0("Performing binning (", bin.size, ") ..."))
  bed.binned <- bedBinner(bed = bed.clean, bin.size = bin.size, nthread = nthread)

  bed.binned <- data.frame(ProbeSetName = paste0("B", seq_len(nrow(bed.binned))), bed.binned, stringsAsFactors = FALSE)

  message("Generating GC% tracks ...")
  wes.gc <- loc.nt.gcc.hs.multi(loc.df = bed.binned, genome.pkg = genome.pkg, extend.multi = extend.multi, blocksize = blocksize, nthread = nthread)

  GC.data <- list(GC = wes.gc, info = list(genome.pkg = genome.pkg, genome = genome, bin.size = bin.size), bed.clean = bed.clean, bed.binned = bed.binned)
  save("GC.data", file = paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = paste0("_", genome, "_b", bin.size, ".rda"), x = basename(bed.file), ignore.case = TRUE)), compress = "xz")

  message("Done.")
  if (return.data) return(GC.data)
}

## Performs the whole processing for WES data from BAMs and a BINpack
EaCoN.WES.Process <- function(testBAM = NULL, refBAM = NULL, BINpack = NULL, samplename = "SAMPLE", Q = 20, L2R.RD.min.Ref = 10, L2R.RD.min.Test = 10, BAF.RD.min = 25, out.dir = getwd(), return.data = FALSE) {
  ifret <- EaCoN.WES.Bin(refBAM = refBAM, testBAM = testBAM, BINpack = BINpack, samplename = samplename, Q = Q, out.dir = out.dir, return.data = TRUE)
  ifret2 <- EaCoN.WES.Normalize(data = ifret, BINpack = BINpack, L2R.RD.min.Ref = L2R.RD.min.Ref, L2R.RD.min.Test = L2R.RD.min.Test, BAF.RD.min = BAF.RD.min, out.dir = out.dir, return.data = return.data)
  message("EaCoN.WES.Process : done.")
  if (return.data) return(ifret2)
}

## Performs the binning of BAMs using a BINpack with Rsamtools
EaCoN.WES.Bin <- function(testBAM = NULL, refBAM = NULL, BINpack = NULL, samplename = "SAMPLE", Q = 20, nsubthread = 1, cluster.type = "PSOCK", out.dir = getwd(), return.data  = FALSE) {
  
  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/MATCHR/TEST")
  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/EaCoN_v0.2.8")
  # refBAM <- "/mnt/data_cigogne/job/PUBLI_EaCoN/MATCHR/BAMS/MR009.normal.MARKdup.bam"
  # refBAM <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/BAMS/TCGA-A7-A0CE-11A-21W-A100-09_IlluminaGA-DNASeq_exome.bam"
  # testBAM <- "/mnt/data_cigogne/job/PUBLI_EaCoN/MATCHR/BAMS/MR009.MARKdup.bam"
  # testBAM <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/BAMS/TCGA-A7-A0CE-01A-11W-A019-09_IlluminaGA-DNASeq_exome.bam"
  # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/MATCHR/RESOURCES/SureSelect_ClinicalResearchExome.padded_hg19_b50.rda"
  # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_hs37d5_b50.rda"
  # samplename <- "MR09TEST_TUMORBOOST"
  # samplename <- "TCGAtest"
  # Q <- 20
  # out.dir = getwd()
  # nsubthread = 3
  # cluster.type = "PSOCK"
  # source("/home/job/svn/genomics/CGH/R/00_PIPELINE/MAIN_MODULES/EaCoN/R/mini_functions.R")
  # require(foreach)
  # require(magrittr)
  # `%do%` <- foreach::"%do%"
  # `%dopar%` <- foreach::"%dopar%"

  ## CHECKS (files/parameters)
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"))
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"))
  if (is.null(refBAM)) stop(tmsg("A reference BAM file is required !"))
  if (!file.exists(refBAM)) stop(tmsg("Could not find the refBAM file !"))
  if (is.null(testBAM)) stop(tmsg("A test BAM file is required !"))
  if (!file.exists(testBAM)) stop(tmsg("Could not find the testBAM file !"))
  if (!is.numeric(Q)) stop(tmsg("Q must be numeric !"))
  if (is.null(out.dir)) stop(tmsg("An output directory is required !"))
  if (!file.exists(out.dir)) stop(tmsg("Could not find the output directory !"))
  if (!file.info(out.dir)$isdir) stop(tmsg("out.dir is not a directory"))
  if (Q < 0) stop(tmsg("Q should be positive !"))
  
  
  # data(list = BSgenome::providerVersion(BSg.obj), package = "chromosomes", envir = environment())
  
  # suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  # genome <- unique(genome(Hsapiens))
  # valid.genomes <- get.valid.genomes()
  # if (!genome %in% c(names(valid.genomes), unlist(valid.genomes))) message(tmsg("Genome is not in the valid list ..."))
  
  
  ## Loading binpack
  message(tmsg("Loading BINpack ..."))
  load(BINpack)
  
  ## CHECKS (genome)
  genome.pkg <- GC.data$info$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }
  
  ### Loading genome
  message(tmsg(paste0("Loading ", genome.pkg, " ...")))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  
  ## Files controls
  message(tmsg("Checking BINpack and BAMs compatibility ..."))
  ## Inspecting BAM headers
  refBAM.h <- Rsamtools::scanBamHeader(refBAM)
  testBAM.h <- Rsamtools::scanBamHeader(testBAM)
  bed.data <- GC.data$bed.binned
  if (!all(names(refBAM.h[[1]]$targets) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and Test BAM are not compatible (different chr names) !"))
  if (!all(bed.data$chr %in% names(refBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and BED are not compatible (different chr names) !"))
  if (!all(unique(bed.data$chr) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Test BAM and BED are not compatible (different chr names) !"))
  if (!all(unique(bed.data$chr) %in% BSgenome::seqnames(BSg.obj))) stop(tmsg("BED and BSgenome are not compatible (different chr names) !"))
  
  # data(list = genome, package = "chromosomes", envir = environment())
  
  ## Identifying the platform
  testBAM.h.unl <- unlist(testBAM.h)
  manuf.hentry <- grep(pattern = "^PL:", testBAM.h.unl)[1]
  manufacturer <- if(!is.na(manuf.hentry)) sub(pattern = "^PL:", replacement = "", testBAM.h.unl[[manuf.hentry]]) else "NA"
  
  meta.b <- list(
    samplename = samplename,
    source = "WES",
    source.file = list(refBAM = refBAM, testBAM = testBAM, BINpack = BINpack),
    type = "WES",
    manufacturer = manufacturer,
    species = GenomeInfoDb::organism(BSg.obj),
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = "NA"
  )
  
  meta.w <- list(
    testBAM.header = paste0(testBAM.h, collapse = " "),
    refBAM.header = paste0(refBAM.h, collapse = " "),
    samtools.Q = Q,
    bin.size = GC.data$info$bin.size
  )
  
  ### Indexing BAM if needed
  if (!file.exists(paste0(testBAM, ".bai"))) Rsamtools::indexBam(testBAM)
  if (!file.exists(paste0(refBAM, ".bai"))) Rsamtools::indexBam(refBAM)
  ### Preparing BAM loading
  param.FLAG <- Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE, isDuplicate=FALSE)
  param.PILEUP.test <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = TRUE)
  # param.PILEUP.ref <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = FALSE)
  param.PILEUP.ref <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = TRUE)
  openBAM.test <- Rsamtools::BamFile(testBAM)
  openBAM.ref <- Rsamtools::BamFile(refBAM)
  ### Loading BED
  # bed.data <- maelstrom::read.table.fast(bed.file)
  
  # colnames(bed.data)[1:3] <- c("chr", 'start', 'end')
  # bzcon <- bzfile(paste0(samplename, "_TEST.bz2"), open = "w")
  # foreach::foreach(k = rev(unique(bed.data$chr))[1:4], .combine = "rbind") %do% {
  # `%do%` <- foreach::"%do%"
  cl <- parallel::makeCluster(spec = nsubthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  k <- 0
  WESdata <- foreach::foreach(k = unique(bed.data$chr), .inorder = TRUE) %dopar% {
    message(tmsg(paste0("Sequence : ", k)))
    ### Loading reference nucleotides
    bed.data.k <- bed.data[bed.data$chr == k,]
    bed.gr <- GenomicRanges::makeGRangesFromDataFrame(bed.data.k, seqnames.field = "chr")
    ### Making pileup
    message(tmsg(" TEST : Getting pileup ..."))
    param.BAM <- Rsamtools::ScanBamParam(which = bed.gr, flag = param.FLAG)
    pres <- dplyr::as.tbl(Rsamtools::pileup(openBAM.test, scanBamParam = param.BAM, pileupParam = param.PILEUP.test))
    pres$nucleotide <- as.character(pres$nucleotide)
    gc()
    ### Building genomic sequence of reference block
    message(tmsg(" TEST : Getting reference genome sequence ..."))
    refblock <- pres[!duplicated(pres$pos), c(1,2)]
    pres <- dplyr::group_by(pres, pos)
    refblock$tot_count <- dplyr::summarize(pres, tot_count = sum(count))$tot_count
    refblock$nucleotide <- unlist(strsplit(as.character(BSgenome::getSeq(BSg.obj, names = GenomicRanges::makeGRangesFromDataFrame(refblock, start.field = "pos", end.field = "pos"))), split = ""))
    ### Merging blocks
    message(tmsg(" TEST : Computing raw BAF ..."))
    refblock <- dplyr::group_by(refblock, pos, nucleotide)
    pres <- dplyr::group_by(pres, pos, nucleotide)
    merged <- dplyr::left_join(refblock, pres, by = c("seqnames", "pos", "nucleotide"))
    rm(pres, refblock)
    gc()
    merged$BAF <- 1 - (merged$count / merged$tot_count)
    merged$nucleotide <- merged$count <- merged$which_label <- NULL
    merged$BAF[merged$BAF == 0] <- NA
    colnames(merged) <- c("chr", "pos", "depth", "BAF")
    
    ### Making pileup
    message(tmsg(" REF : Getting pileup ..."))
    presR <- dplyr::as.tbl(Rsamtools::pileup(openBAM.ref, scanBamParam = param.BAM, pileupParam = param.PILEUP.ref))
    presR$nucleotide <- as.character(presR$nucleotide)
    gc()
    ### Building genomic sequence of reference block
    message(tmsg(" REF : Getting reference genome sequence ..."))
    refblockR <- presR[!duplicated(presR$pos), c(1,2)]
    presR <- dplyr::group_by(presR, pos)
    refblockR$tot_count <- dplyr::summarize(presR, tot_count = sum(count))$tot_count
    refblockR$nucleotide <- unlist(strsplit(as.character(BSgenome::getSeq(BSg.obj, names = GenomicRanges::makeGRangesFromDataFrame(refblockR, start.field = "pos", end.field = "pos"))), split = ""))
    ### Merging blocks
    message(tmsg(" REF : Computing raw BAF ..."))
    refblockR <- dplyr::group_by(refblockR, pos, nucleotide)
    presR <- dplyr::group_by(presR, pos, nucleotide)
    mergedR <- dplyr::left_join(refblockR, presR, by = c("seqnames", "pos", "nucleotide"))
    rm(presR, refblockR)
    gc()
    mergedR$BAF <- 1 - (mergedR$count / mergedR$tot_count)
    mergedR$nucleotide <- mergedR$count <- mergedR$which_label <- NULL
    mergedR$BAF[mergedR$BAF == 0] <- NA
    colnames(mergedR) <- c("chr", "pos", "depth", "BAF")
    
    bed.based <- data.frame(chr = k, pos = unlist(seq.int2(from = bed.data.k$start, to = bed.data.k$end, by = 1)), bin = rep(x = seq_along(bed.data.k$ProbeSetName), times = bed.data.k$end - bed.data.k$start +1), depth.test = 0, depth.ref = 0, BAF.test = NA, BAF.ref = NA)
    binT <- bed.based$pos %in% merged$pos
    Tinb <- merged$pos %in% bed.based$pos
    if (length(binT) > 0) {
      bed.based$depth.test[binT] <- merged$depth[Tinb]
      bed.based$BAF.test[binT] <- merged$BAF[Tinb]
    }
    binR <- bed.based$pos %in% mergedR$pos
    Rinb <- mergedR$pos %in% bed.based$pos
    if (length(binR) > 0) {
      bed.based$depth.ref[binR] <- mergedR$depth[Rinb]
      bed.based$BAF.ref[binR] <- mergedR$BAF[Rinb]
    }
    
    CN.table <- data.frame(dplyr::summarize(dplyr::group_by(dplyr::as.tbl(bed.based), bin), chr = unique(chr), start = min(pos), end = max(pos), RD.test.mean = mean(depth.test), RD.ref.mean = mean(depth.ref)))
    CN.table$start <- bed.data.k$start
    CN.table$end <- bed.data.k$end
    BAF.table <- data.frame(bed.based[!is.na(bed.based$BAF.test), c(3,1,2,4:7)])
    colnames(BAF.table)[c(4,5)] <- c("RD.test", "RD.ref")
    # bed.based$pos[which(diff(bed.based$pos) == 1)+1] <- NA
    
    # write.table(x = bed.based, file = bzcon, col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
    rm(merged, mergedR, bed.based)
    gc()
    return(list(CN = CN.table, BAF = BAF.table))
    
    # oddz <- (BAF.table$pos %% 2) != 0
    # BAF.table$BAF.test2 <- BAF.table$BAF.test
    # BAF.table$BAF.ref2 <- BAF.table$BAF.ref
    # BAF.table$BAF.test2[oddz] <- 1 - BAF.table$BAF.test2[oddz]
    # BAF.table$BAF.ref2[oddz] <- 1 - BAF.table$BAF.ref2[oddz]
    # BAF.table$geno.ref <- BAF.table$BAF.ref2
    # BAF.table$geno.ref[BAF.table$geno.ref < .25] <- 0
    # BAF.table$geno.ref[BAF.table$geno.ref > .75] <- 1
    # BAF.table$geno.ref[BAF.table$geno.ref >= .25 & BAF.table$geno.ref <= .75 ] <- .5
    # 
    # plot(BAF.table$BAF.test2, pch = ".", cex = 3)
    # aroma.cn
    
  }
  # message("Stopping cluster ...")
  parallel::stopCluster(cl)
  
  message("Done.")
  # close(bzcon)
  
  ## Untangling WESobj
  CN <- foreach::foreach(x = seq_len(length(WESdata)), .combine = "rbind") %do% {
    toret <- WESdata[[x]]$CN
    WESdata[[x]]$CN <- NULL
    return(toret)
  }
  BAF <- foreach::foreach(x = seq_len(length(WESdata)), .combine = "rbind") %do% {
    toret <- WESdata[[x]]$BAF
    WESdata[[x]]$BAF <- NULL
    return(toret)
  }
  rm(WESdata)
  
  message(tmsg("Saving counts data ..."))
  meta.w$RD.test.mean.summary <- as.vector(summary(CN$RD.test.mean))
  meta.w$RD.ref.mean.summary <- as.vector(summary(CN$RD.ref.mean))
  meta.w$BAF.RD.test.summary <- as.vector(summary(BAF$RD.test))
  names(meta.w$BAF.RD.test.summary) <- names(meta.w$RD.test.mean.summary) <- names(meta.w$RD.ref.mean.summary) <- c("min", "q25", "median", "mean", "q75", "max")
  WESobj <- list(RD = CN, SNP = BAF, meta = list(basic = meta.b, WES = meta.w))
  rm(list = c("CN", "BAF"))
  
  ## QC : Computing coverages
  message(tmsg("Computing coverages ..."))
  gw.rd <- sum(WESobj$RD$end - WESobj$RD$start +1)
  gw.snp <- nrow(WESobj$SNP)
  rd.cov <- data.frame(cuts = c(1, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200), stringsAsFactors = FALSE)
  rd.cov <- cbind(rd.cov, t(foreach::foreach(x = rd.cov$cuts, .combine = "cbind") %do% {
    test.rd.in <- WESobj$RD$RD.test.mean >= x
    ref.rd.in <- WESobj$RD$RD.ref.mean >= x
    test.snprd.in <- WESobj$SNP$RD.test >= x
    ref.snprd.in <- WESobj$SNP$RD.ref >= x
    test.cut.cov <- if(!any(test.rd.in)) NA else (sum(WESobj$RD$end[test.rd.in] - WESobj$RD$start[test.rd.in] +1)/gw.rd)
    ref.cut.cov <- if(!any(ref.rd.in)) NA else (sum(WESobj$RD$end[ref.rd.in] - WESobj$RD$start[ref.rd.in] +1)/gw.rd)
    test.snpcut.cov <- if(!any(test.snprd.in)) NA else (length(which(test.snprd.in))/gw.snp)
    ref.snpcut.cov <- if(!any(ref.snprd.in)) NA else (length(which(ref.snprd.in))/gw.snp)
    return(c(test.cut.cov, ref.cut.cov, test.snpcut.cov, ref.snpcut.cov))
  }))
  colnames(rd.cov) <- c("MinDepth", "TestBINCoverage", "RefBINCoverage", "TestBAFCoverage", "RefBAFCoverage")
  
  ## QC : Plotting coverages
  dir.create(paste0(out.dir, "/", samplename))
  png(paste0(out.dir, "/", samplename, "/", samplename, "_WES_", WESobj$meta$basic$genome, "_coverage.png"), 800, 640)
  plot(rd.cov$MinDepth, rd.cov$TestBAFCoverage, type = "b", col = 2, lty = 3, pch = 20, main = paste0(WESobj$meta$basic$samplename, "\nCoverage Plot"), xlab = "Minimum depth", ylab = "Coverage", ylim = c(0,1), xaxp = c(0,200,10))
  abline(v = rd.cov$MinDepth, lty = 2, col = "grey75")
  abline(h = seq(0,1,.1), lty = 2, col = "grey75")
  lines(rd.cov$MinDepth, rd.cov$RefBAFCoverage, type = "b", col = 1, lty = 3, pch = 20)
  lines(rd.cov$MinDepth, rd.cov$TestBINCoverage, type = "b", col = 2)
  lines(rd.cov$MinDepth, rd.cov$RefBINCoverage, type = "b", col = 1)
  abline(h = .5, lty = 2)
  legend("topright", legend = c("Test BAF", "Ref BAF", "Test BIN", "Ref BIN"), inset = .02, col = c(2,1,2,1), lty = c(3,3,1,1), pch = c(20,20,1,1))
  dev.off()
  
  ## Saving
  saveRDS(WESobj, file = paste0(out.dir, "/", samplename, "/", samplename, "_", genome, "_b", meta.w$bin.size, "_binned.RDS"), compress = "xz")
  if (return.data) return(WESobj)
}

## DEPRECATED. Performs the binning of BAMs using a BINpack with Rsamtools
EaCoN.WES.Bin.OLD <- function(testBAM = NULL, refBAM = NULL, BINpack = NULL, samplename = "SAMPLE", Q = 20, out.dir = getwd(), return.data  = FALSE) {

  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/EaCoN_v0.2.8")
  # refBAM <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/BAMS/TCGA-A7-A0CE-11A-21W-A100-09_IlluminaGA-DNASeq_exome.bam"
  # testBAM <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/BAMS/TCGA-A7-A0CE-01A-11W-A019-09_IlluminaGA-DNASeq_exome.bam"
  # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_GRCh37-lite_b50.rda"
  # samplename <- "BAFCHR"
  # genome.pkg <- "BSgenome.Hsapiens.TCGA.GRCh37.lite"
  # Q <- 20
  # out.dir = getwd()
  # source("/home/job/svn/genomics/CGH/R/00_PIPELINE/MAIN_MODULES/EaCoN/R/mini_functions.R")
  # require(foreach)
  # #
  ## CLEANED !
  # require(magrittr)
  `%do%` <- foreach::"%do%"

  ## CHECKS (files/parameters)
  # if (is.null(bed.file)) stop("A BED file is required !")
  # if (!file.exists(bed.file)) stop("Could not find the BED file !")
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"))
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"))
  # if (is.null(fasta)) stop("A FASTA file is required !")
  # if (!file.exists(fasta)) stop("Could not find FASTA file !")
  # if (!file.exists(fasta.fai)) stop("Could not find the index for FASTA file ! Please use FastaIndex().")
  if (is.null(refBAM)) stop(tmsg("A reference BAM file is required !"))
  if (!file.exists(refBAM)) stop(tmsg("Could not find the refBAM file !"))
  if (is.null(testBAM)) stop(tmsg("A test BAM file is required !"))
  if (!file.exists(testBAM)) stop(tmsg("Could not find the testBAM file !"))
  if (!is.numeric(Q)) stop(tmsg("Q must be numeric !"))
  if (is.null(out.dir)) stop(tmsg("An output directory is required !"))
  if (!file.exists(out.dir)) stop(tmsg("Could not find the output directory !"))
  if (!file.info(out.dir)$isdir) stop(tmsg("out.dir is not a directory"))
  if (Q < 0) stop(tmsg("Q should be positive !"))


  # data(list = BSgenome::providerVersion(BSg.obj), package = "chromosomes", envir = environment())

  # suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  # genome <- unique(genome(Hsapiens))
  # valid.genomes <- get.valid.genomes()
  # if (!genome %in% c(names(valid.genomes), unlist(valid.genomes))) message(tmsg("Genome is not in the valid list ..."))


  ## Loading binpack
  message(tmsg("Loading BINpack ..."))
  load(BINpack)

  ## CHECKS (genome)
  genome.pkg <- GC.data$info$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }

  ### Loading genome
  message(tmsg(paste0("Loading ", genome.pkg, " ...")))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)

  ## Files controls
  message(tmsg("Checking BINpack and BAMs compatibility ..."))
  ## Inspecting BAM headers
  refBAM.h <- Rsamtools::scanBamHeader(refBAM)
  testBAM.h <- Rsamtools::scanBamHeader(testBAM)
  bed.data <- GC.data$bed.binned
  if (!all(names(refBAM.h[[1]]$targets) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and Test BAM are not compatible (different chr names) !"))
  if (!all(bed.data$chr %in% names(refBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and BED are not compatible (different chr names) !"))
  if (!all(unique(bed.data$chr) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Test BAM and BED are not compatible (different chr names) !"))
  if (!all(unique(bed.data$chr) %in% BSgenome::seqnames(BSg.obj))) stop(tmsg("BED and BSgenome are not compatible (different chr names) !"))

  # data(list = genome, package = "chromosomes", envir = environment())

  ## Identifying the platform
  testBAM.h.unl <- unlist(testBAM.h)
  manuf.hentry <- grep(pattern = "^PL:", testBAM.h.unl)[1]
  manufacturer <- if(!is.na(manuf.hentry)) sub(pattern = "^PL:", replacement = "", testBAM.h.unl[[manuf.hentry]]) else "NA"

  meta.b <- list(
    samplename = samplename,
    source = "WES",
    source.file = list(refBAM = refBAM, testBAM = testBAM, BINpack = BINpack),
    type = "WES",
    manufacturer = manufacturer,
    species = GenomeInfoDb::organism(BSg.obj),
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = "XX"
  )

  meta.w <- list(
    testBAM.header = paste0(testBAM.h, collapse = " "),
    refBAM.header = paste0(refBAM.h, collapse = " "),
    samtools.Q = Q,
    bin.size = GC.data$info$bin.size
  )

  ### Indexing BAM if needed
  if (!file.exists(paste0(testBAM, ".bai"))) Rsamtools::indexBam(testBAM)
  if (!file.exists(paste0(refBAM, ".bai"))) Rsamtools::indexBam(refBAM)
  ### Preparing BAM loading
  param.FLAG <- Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE, isDuplicate=FALSE)
  param.PILEUP.test <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = TRUE)
  param.PILEUP.ref <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = FALSE)
  openBAM.test <- Rsamtools::BamFile(testBAM)
  openBAM.ref <- Rsamtools::BamFile(refBAM)
  ### Loading BED
  # bed.data <- maelstrom::read.table.fast(bed.file)

  # colnames(bed.data)[1:3] <- c("chr", 'start', 'end')
  # bzcon <- bzfile(paste0(samplename, "_TEST.bz2"), open = "w")
  # foreach::foreach(k = rev(unique(bed.data$chr))[1:4], .combine = "rbind") %do% {
  `%do%` <- foreach::"%do%"
  WESdata <- foreach::foreach(k = unique(bed.data$chr)) %do% {
    message(tmsg(paste0("Sequence : ", k)))
    ### Loading reference nucleotides
    bed.data.k <- bed.data[bed.data$chr == k,]
    bed.gr <- GenomicRanges::makeGRangesFromDataFrame(bed.data.k, seqnames.field = "chr")
    ### Making pileup
    message(tmsg(" TEST : Getting pileup ..."))
    param.BAM <- Rsamtools::ScanBamParam(which = bed.gr, flag = param.FLAG)
    pres <- dplyr::as.tbl(Rsamtools::pileup(openBAM.test, scanBamParam = param.BAM, pileupParam = param.PILEUP.test))
    pres$nucleotide <- as.character(pres$nucleotide)
    gc()
    ### Building reference block
    message(tmsg(" TEST : Getting reference genome sequence ..."))
    refblock <- pres[!duplicated(pres$pos), c(1,2)]
    pres <- dplyr::group_by(pres, pos)
    refblock$tot_count <- dplyr::summarize(pres, tot_count = sum(count))$tot_count
    refblock$nucleotide <- unlist(strsplit(as.character(BSgenome::getSeq(BSg.obj, names = GenomicRanges::makeGRangesFromDataFrame(refblock, start.field = "pos", end.field = "pos"))), split = ""))
    ### Merging blocks
    message(tmsg(" TEST : Computing raw BAF ..."))
    refblock <- dplyr::group_by(refblock, pos, nucleotide)
    pres <- dplyr::group_by(pres, pos, nucleotide)
    merged <- dplyr::left_join(refblock, pres, by = c("seqnames", "pos", "nucleotide"))
    rm(pres, refblock)
    gc()
    merged$BAF <- 1 - (merged$count / merged$tot_count)
    merged$nucleotide <- merged$count <- merged$which_label <- NULL
    merged$BAF[merged$BAF == 0] <- NA
    colnames(merged) <- c("chr", "pos", "depth", "BAF")

    message(tmsg(" REF : Getting pileup ..."))
    presR <- dplyr::as.tbl(Rsamtools::pileup(openBAM.ref, scanBamParam = param.BAM, pileupParam = param.PILEUP.ref))
    # presR$seqnames <- as.character(presR$seqnames)
    # presR$which_label <- as.character(presR$which_label)


    # bed.based <- data.frame(chr = k, pos = unlist(seq.int2(from = bed.data.k$start, to = bed.data.k$end, by = 1)), bin = rep(x = bed.data.k$ProbeSetName, times = bed.data.k$end - bed.data.k$start +1), depth.test = 0, depth.ref = 0, BAF = NA)
    bed.based <- data.frame(chr = k, pos = unlist(seq.int2(from = bed.data.k$start, to = bed.data.k$end, by = 1)), bin = rep(x = seq_along(bed.data.k$ProbeSetName), times = bed.data.k$end - bed.data.k$start +1), depth.test = 0, depth.ref = 0, BAF = NA)
    binT <- bed.based$pos %in% merged$pos
    Tinb <- merged$pos %in% bed.based$pos
    if (length(binT) > 0) {
      bed.based$depth.test[binT] <- merged$depth[Tinb]
      bed.based$BAF[binT] <- merged$BAF[Tinb]
    }
    binR <- bed.based$pos %in% presR$pos
    Rinb <- presR$pos %in% bed.based$pos
    if (length(binR) > 0) bed.based$depth.ref[binR] <- presR$count[Rinb]

    CN.table <- data.frame(dplyr::summarize(dplyr::group_by(dplyr::as.tbl(bed.based), bin), chr = unique(chr), start = min(pos), end = max(pos), RD.test.mean = mean(depth.test), RD.ref.mean = mean(depth.ref)))
    CN.table$start <- bed.data.k$start
    CN.table$end <- bed.data.k$end
    BAF.table <- data.frame(bed.based[!is.na(bed.based$BAF), c(3,1,2,4,6)])
    colnames(BAF.table)[4] <- "RD.test"
    # bed.based$pos[which(diff(bed.based$pos) == 1)+1] <- NA

    # write.table(x = bed.based, file = bzcon, col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
    rm(merged, bed.based, presR)
    return(list(CN = CN.table, BAF = BAF.table))
  }
  # close(bzcon)

  ## Untangling WESobj
  CN <- foreach::foreach(x = seq_len(length(WESdata)), .combine = "rbind") %do% {
    toret <- WESdata[[x]]$CN
    WESdata[[x]]$CN <- NULL
    return(toret)
  }
  BAF <- foreach::foreach(x = seq_len(length(WESdata)), .combine = "rbind") %do% {
    toret <- WESdata[[x]]$BAF
    WESdata[[x]]$BAF <- NULL
    return(toret)
  }
  rm(WESdata)

  ## Saving
  message(tmsg("Saving counts data ..."))
  meta.w$RD.test.mean.summary <- as.vector(summary(CN$RD.test.mean))
  meta.w$RD.ref.mean.summary <- as.vector(summary(CN$RD.ref.mean))
  meta.w$BAF.RD.test.summary <- as.vector(summary(BAF$RD.test))
  names(meta.w$BAF.RD.test.summary) <- names(meta.w$RD.test.mean.summary) <- names(meta.w$RD.ref.mean.summary) <- c("min", "q25", "median", "mean", "q75", "max")
  WESobj <- list(RD = CN, SNP = BAF, meta = list(basic = meta.b, WES = meta.w))
  dir.create(paste0(out.dir, "/", samplename))
  saveRDS(WESobj, file = paste0(out.dir, "/", samplename, "/", samplename, "_", genome, "_b", meta.w$bin.size, "_binned.RDS"), compress = "xz")
  if (return.data) return(WESobj)
}

## Performs the binning of BAMs using a BINpack, batch mode
EaCoN.WES.Bin.Batch <- function(BAM.list.file = NULL, BINpack = NULL, nthread = 1, cluster.type = "PSOCK", ...) {

  if (!file.exists(BAM.list.file)) stop("Could not find BAM.list.file !")
  message("Reading and checking BAM.list.file ...")
  myBAMs <- read.table(file = BAM.list.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("testBAM", "refBAM", "SampleName")
  head.chk <- all(colnames(BAM.list.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in BAM.list.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(myBAMs)))
    stop("Invalid header.")
  }

  tb.chk <- file.exists(myBAMs$testBAM)
  rb.chk <- file.exists(myBAMs$refBAM)

  if (!all(tb.chk) || !all(tb.chk)) {
    message("Some BAM files from the BAM.list.file could not be found (wrong path or filename ?) !")
    message("Missing testBAM file(s) :")
    message(myBAMs$testBAM[which(!tb.chk)])
    message("Missing refBAM file(s) :")
    message(myBAMs$refBAM[which(!rb.chk)])
    stop("Missing BAM file(s).")
  }
  sn.chk <- duplicated(myBAMs$SampleName)
  if (any(sn.chk)) {
    message("BAM.list.file contains duplicated samplenames !")
    message(myBAMs$SampleName[which(duplicated(myBAMs$SampleName))])
    stop("Duplicated samplenames")
  }
  if(any(myBAMs$testBAM == myBAMs$refBAM)) {
    message("Some testBAM and refBAM are identical for at least one sample !")
    stop("Identical BAM files for Test and Ref.")
  }

  ## Adjusting cores/threads
  message("Adjusting number of cores if needed ...")
  if (is.null(nthread)) nthread <- parallel::detectCores(logical = TRUE) -1
  if (nrow(myBAMs) < nthread) nthread <- nrow(myBAMs)

  message("Running EaCoN.WES.Bin() in batch mode ...")
  message(paste0("Found ", nrow(myBAMs), " samples to process ..."))
  current.bitmapType <- getOption("bitmapType")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  eacon.batchres <- foreach::foreach(r = seq_len(nrow(myBAMs)), .inorder = TRUE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    EaCoN.WES.Bin(testBAM = myBAMs$testBAM[r], refBAM = myBAMs$refBAM[r], BINpack = BINpack, samplename = myBAMs$SampleName[r], cluster.type = cluster.type, ...)
  }
  parallel::stopCluster(cl)
}

## Performs the normalization of WES L2R and BAF signals
EaCoN.WES.Normalize <- function(data = NULL, BINpack = NULL, L2R.RD.min.Ref = 20, L2R.RD.min.Test = 20, BAF.RD.min = 25,TumorBoost = TRUE, out.dir = getwd(), return.data = FALSE) {

  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/EaCoN_v0.2.8/TCGA-A7-A0CE-01A_vs_11A/")
  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/EaCoN_v0.2.8/TCGA-BH-A0DT-01A_vs_11A/")
  # setwd("/mnt/data_cigogne/job/OS2006/WES/OS2K6/EaCoN_0.2.8/00_TONORM/OS_046")
  # data <- readRDS("TCGA-BH-A0DT-01A_vs_11A_hs37d5_b50_binned.RDS")
  # data <- readRDS("OS_046_hg19_b50_binned.RDS")
  # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_hs37d5_b50.rda"
  # BINpack <- "/mnt/data_cigogne/job/OS2006/WES/RESOURCES/SureSelect_ClinicalResearchExome.padded_hg19_b50.rda"
  # L2R.RD.min.Ref = 20
  # L2R.RD.min.Ref = 10
  # L2R.RD.min.Test = 20
  # L2R.RD.min.Test = 10
  # BAF.RD.min = 25
  # BAF.RD.min = 10
  # TumorBoost = TRUE
  # out.dir = getwd()
  # return.data = FALSE
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("/home/job/git_gustaveroussy/EaCoN/R/fit_functions.R")
  # require(foreach)

  
  
  ## CHECKS
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"))
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"))

  ## Loading BINpack
  load(BINpack)
  genome.pkg <- GC.data$info$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }
  
  ### Loading genome
  message(tmsg(paste0("Loading ", genome.pkg, " ...")))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)
  
  ## DEPTH controls
  ### L2R.TEST
  LT <- length(which(data$RD$RD.test.mean < L2R.RD.min.Test)) / nrow(data$RD)
  LR <- length(which(data$RD$RD.ref.mean < L2R.RD.min.Ref)) / nrow(data$RD)
  LB <- length(which(data$RD$RD.test.mean < L2R.RD.min.Test | data$RD$RD.ref.mean < L2R.RD.min.Ref)) / nrow(data$RD)
  message(tmsg(paste0("Using L2R.RD.min.Test = ", L2R.RD.min.Test, ", ", round(LT*100, digits = 2), "% of bins are discarded.")))
  message(tmsg(paste0("Using L2R.RD.min.Ref = ", L2R.RD.min.Ref, ", ", round(LR*100, digits = 2), "% of bins are discarded.")))
  message(tmsg(paste0("Using both, ", round(LB*100, digits = 2), "% of bins are discarded.")))
  if (LB >= .5) {
    message(tmsg("More than 50% of bins would have been set to get imputed, which is not possible ! Please lower the L2R.RD.min.Test and/or L2R.RD.min.Ref parameters, when possible. Note that this will increase profile noise."))
    stop("More than 50% of bins would have been set to get imputed, which is not possible ! Please lower the L2R.RD.min.Test and/or L2R.RD.min.Ref parameters, when possible. Note that this will increase profile noise.")
  }

  meta.b <- data$meta$basic
  meta.w <- data$meta$WES
  samplename <- meta.b$samplename

  meta.w$L2R.RD.min.Ref <- L2R.RD.min.Ref
  meta.w$L2R.RD.min.Test <- L2R.RD.min.Test
  meta.w$TumorBoost <- as.character(TumorBoost)

  ## Normalizing L2R
  # data(list = data$meta$basic$genome, package = "chromosomes", envir = environment())
  rownames(data$RD) <- paste0(data$RD$chr, ":", data$RD$start, "-", data$RD$end)
  # data$RD$chrN <- unlist(cs$chrom2chr[data$RD$chr])
  data$RD$chrN <- unclass(data$RD$chr)
  smo <- round(nrow(data$RD)/200)
  data$RD$L2R <- log2((data$RD$RD.test.mean+1) / (data$RD$RD.ref.mean+1))

  ## Handling low-depth bins
  message(tmsg("Imputing flagged bins L2r (low-depth, GC%-outiler) ..."))
  ### low depth
  ldbins <- data$RD$RD.test.mean < L2R.RD.min.Test | data$RD$RD.ref.mean < L2R.RD.min.Ref
  message(tmsg(paste0(" Found ", length(which(ldbins)), " low-depth bin(s).")))
  
  ### GC outliers
  gcobins <- GC.data$GC[,5] < .2 | GC.data$GC[,5] > .8
  message(tmsg(paste0(" Found ", length(which(gcobins)), " GC%-outlier bin(s).")))
  ### Pooled
  fbins <- ldbins + gcobins > 0

  if (any(ldbins)) {
    message(tmsg(paste0(" Imputed ", length(which(fbins)), " bins L2R.")))
    l2r.tmp <- data$RD$L2R
    l2r.tmp[fbins] <- NA
    data$RD$L2Ri <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  } else data$RD$L2Ri <- data$RD$L2R

  ## Pseudo-centering on the median
  data$RD[,ncol(data$RD)] <- data$RD[,ncol(data$RD)] - median(data$RD[,ncol(data$RD)], na.rm = TRUE)

  ## Normalizing L2R using GC content
  message(tmsg("GC-normalization ..."))
  rownames(GC.data$GC) <- paste0(GC.data$GC$chr, ":", GC.data$GC$start, "-", GC.data$GC$end)
  if (any(rownames(GC.data$GC) != rownames(data$RD))) stop(tmsg("GC data and L2R data are not synched, or ordered differently !"))
  ndata <- data.frame(data$RD[,1:3], name = rownames(data$RD), GC.data$GC[,-c(1:4)], stringsAsFactors = FALSE)
  my.rm.mad <- sum(abs(diff(as.numeric(runmed(data$RD[!is.na(data$RD$L2Ri),ncol(data$RD)], smo)))))
  normloop.res <- l2r.fitloop(l2rObj = list(l2r=data$RD[,ncol(data$RD)], rm.mad = my.rm.mad), tfd = ndata, smo = smo)
  rm(ndata)
  data$RD$L2Rrn <- normloop.res$l2r$l2r + median(data$RD[,ncol(data$RD)], na.rm = TRUE)
  meta.b$renorm <- if(is.null(normloop.res$pos)) "None" else paste0(normloop.res$pos, collapse = ",")

  ## Preparing BAF
  message(tmsg("Normalizing BAF data ..."))
  mid.medrefrd <- median(data$SNP$RD.ref, na.rm = TRUE)/2
  if(BAF.RD.min > mid.medrefrd) {
    message(tmsg(paste0("WARNING : set value for BAF.RD.min (", BAF.RD.min, ") is lower than half the median depth of reference SNPs (", mid.medrefrd, "), so its value was adjusted to the latter.")))
    BAF.RD.min <- mid.medrefrd
  }
  meta.w$BAF.RD.min <- BAF.RD.min
  ldbaf <- data$SNP$RD.test < BAF.RD.min
  BT <- length(which(ldbaf)) / nrow(data$SNP)
  if (any(ldbaf)) {
    data$SNP <- data$SNP[!ldbaf,]
    message(tmsg(paste0(" Removed ", length(which(ldbaf)), " low-depth SNPs (", round(BT*100, digits = 2), " %).")))
  }
  data$SNP$BAF.test[data$SNP$BAF.test < 0] <- 0
  data$SNP$BAF.test[data$SNP$BAF.test > 1] <- 1
  data$SNP$BAF.ref[data$SNP$BAF.ref < 0] <- 0
  data$SNP$BAF.ref[data$SNP$BAF.ref > 1] <- 1
  odd.idx <- which(data$SNP$pos %% 2 == 1)
  data$SNP$BAF.test[odd.idx] <- -data$SNP$BAF.test[odd.idx] +1
  data$SNP$BAF.ref[odd.idx] <- -data$SNP$BAF.ref[odd.idx] +1
  
  ## TumorBoost
  if (TumorBoost) {
    message(tmsg("Applying TumorBoost BAF normalization ..."))
    data$SNP$BAF.test <- as.numeric(aroma.light::normalizeTumorBoost(data$SNP$BAF.test, data$SNP$BAF.ref, flavor = "v4", preserveScale = FALSE))
  }
  
  ## Building ASCAT-like object
  message(tmsg("Building normalized object ..."))
  SNPz <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = data$SNP$chr, start = data$SNP$pos, end = data$SNP$pos, stringsAsFactors = FALSE))
  L2Rz <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = data$RD$chr, start = data$RD$start, end = data$RD$end, stringsAsFactors = FALSE))
  SinL <- GenomicRanges::findOverlaps(SNPz, L2Rz)
  data$SNP$L2R <- NA
  data$SNP$L2R[SinL@from] <- data$RD[, ncol(data$RD)][SinL@to]

  l2r.ao.df <- data.frame(chr = data$RD$chr, pos = round((data$RD$start + data$RD$end) / 2), L2R.ori = data$RD$L2R, L2R = data$RD[,ncol(data$RD)], BAF = NA, stringsAsFactors = FALSE)
  rownames(l2r.ao.df) <- paste0(l2r.ao.df$chr, ":", l2r.ao.df$pos)
  baf.ao.df <- data.frame(chr = data$SNP$chr, pos = data$SNP$pos, L2R.ori = data$SNP$L2R, L2R = data$SNP$L2R, BAF = data$SNP$BAF.test, stringsAsFactors = FALSE)
  rownames(baf.ao.df) <- rownames(data$SNP)
  ao.df <- rbind(l2r.ao.df, baf.ao.df)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  # ao.df$chrN <- unclass(as.vector(ao.df$chr))
  ao.df$chrs <- ao.df$chr
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos),]

  my.ch <- sapply(unique(ao.df$chrs), function(x) { which(ao.df$chrs == x) })

  my.ascat.obj <- list(
    data = list(
      Tumor_LogR = data.frame(sample = ao.df$L2R),
      Tumor_BAF = data.frame(sample = ao.df$BAF),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      SNPpos = data.frame(chrs = ao.df$chrs, pos = ao.df$pos),
      ch = my.ch,
      chr = my.ch,
      chrs = levels(ao.df$chrs),
      samples = samplename,
      gender = "NA",
      sexchromosomes = c("X", "Y"),
      failedarrays = NULL
    ),
    meta = list(
      basic = meta.b,
      wes = meta.w
    )
  )
  colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename

  ## Saving normalized object
  message(tmsg("Saving normalized / prepared data ..."))
  # new.rds.name <- sub(pattern = "\\.RDS$", replacement = "_normalized.RDS", x = data.file)
  # saveRDS(my.ascat.obj, file = new.rds.name, compress = "bzip2")
  saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "_", data$meta$basic$genome, "_b", data$meta$WES$bin.size, "_processed.RDS"), compress = "xz")

  ## Rough plot
  message(tmsg("Plotting ..."))
  l2r.notna <- which(!is.na(ao.df$L2R))
  l2r.rm <- runmed(ao.df$L2R[l2r.notna], smo)
  l2r.mad <- median(abs(diff(ao.df$L2R[l2r.notna])))
  l2r.ssad <- sum(abs(diff(l2r.rm)))
  l2r.ori.rm <- runmed(ao.df$L2R.ori[l2r.notna], smo)
  l2r.ori.mad <- median(abs(diff(ao.df$L2R.ori[l2r.notna])))
  l2r.ori.ssad <- sum(abs(diff(l2r.ori.rm)))
  ao.df$genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  ao.df$L2R.ori <- ao.df$L2R.ori - median(ao.df$L2R.ori, na.rm = TRUE)
  kend <- ao.df$genopos[vapply(unique(ao.df$chrN), function(k) { max(which(ao.df$chrN == k))}, 1)]
  png(paste0(out.dir, "/", samplename, "_WES_", data$meta$basic$genome, "_rawplot.png"), 1600, 1050)
  # par(mfrow = c(2,1))
  par(mfrow = c(3,1))
  plot(ao.df$genopos, ao.df$L2R.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") raw L2R profile (median-centered)\nMAD = ", round(l2r.ori.mad, digits = 2), " ; SSAD = ", round(l2r.ori.ssad, digits = 2)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.ori.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$L2R, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") normalized L2R profile (median-centered)\nMAD = ", round(l2r.mad, digits = 2), " ; SSAD = ", round(l2r.ssad, digits = 2)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$BAF, pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ")", if(TumorBoost) " TumorBoost-normalized", " BAF profile"), xlab = "Genomic position", ylab = "BAF")
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = .5, col = 2, lty = 2, lwd = 2)
  dev.off()

  message(tmsg("Done."))
  if(return.data) return(my.ascat.obj)
}

## Performs the normalization of WES L2R and BAF signals
EaCoN.WES.Normalize.OLD <- function(data = NULL, BINpack = NULL, L2R.RD.min.Ref = 20, L2R.RD.min.Test = 20, BAF.RD.min = 25, TumorBoost = TRUE, out.dir = getwd(), return.data = FALSE) {
  
  # # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/EaCoN_v0.2.8/TCGA-A7-A0CE-01A_vs_11A/")
  # setwd("/mnt/data_cigogne/job/OS2006/WES/OS2K6/EaCoN_0.2.8/00_TONORM/OS_046")
  # # data <- readRDS("TCGA-A7-A0CE-01A_vs_11A_hs37d5_b50_binned.RDS")
  # data <- readRDS("OS_046_hg19_b50_binned.RDS")
  # # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_hs37d5_b50.rda"
  # BINpack <- "/mnt/data_cigogne/job/OS2006/WES/RESOURCES/SureSelect_ClinicalResearchExome.padded_hg19_b50.rda"
  # # L2R.RD.min.Ref = 20
  # L2R.RD.min.Ref = 10
  # # L2R.RD.min.Test = 20
  # L2R.RD.min.Test = 10
  # # BAF.RD.min = 25
  # BAF.RD.min = 10
  # TumorBoost = TRUE
  # out.dir = getwd()
  # return.data = FALSE
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("/home/job/git_gustaveroussy/EaCoN/R/fit_functions.R")
  # require(foreach)
  
  
  
  ## CHECKS
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"))
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"))
  
  ## Loading BINpack
  load(BINpack)
  genome.pkg <- GC.data$info$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }
  
  ### Loading genome
  message(tmsg(paste0("Loading ", genome.pkg, " ...")))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  
  
  ## DEPTH controls
  ### L2R.TEST
  LT <- length(which(data$RD$RD.test.mean < L2R.RD.min.Test)) / nrow(data$RD)
  LR <- length(which(data$RD$RD.ref.mean < L2R.RD.min.Ref)) / nrow(data$RD)
  LB <- length(which(data$RD$RD.test.mean < L2R.RD.min.Test | data$RD$RD.ref.mean < L2R.RD.min.Ref)) / nrow(data$RD)
  message(tmsg(paste0("Using L2R.RD.min.Test = ", L2R.RD.min.Test, ", ", round(LT*100, digits = 2), "% of bins are discarded.")))
  message(tmsg(paste0("Using L2R.RD.min.Ref = ", L2R.RD.min.Ref, ", ", round(LR*100, digits = 2), "% of bins are discarded.")))
  message(tmsg(paste0("Using both, ", round(LB*100, digits = 2), "% of bins are discarded.")))
  if (LB >= .5) stop(tmsg("More than 50% of bins would have been set to get imputed, which is not possible ! Please lower the L2R.RD.min.Test and/or L2R.RD.min.Ref parameters, when possible. Note that this will increase profile noise."))
  
  meta.b <- data$meta$basic
  meta.w <- data$meta$WES
  samplename <- meta.b$samplename
  
  meta.w$L2R.RD.min.Ref <- L2R.RD.min.Ref
  meta.w$L2R.RD.min.Test <- L2R.RD.min.Test
  meta.w$BAF.RD.min <- BAF.RD.min
  
  ## Normalizing L2R
  data(list = data$meta$basic$genome, package = "chromosomes", envir = environment())
  rownames(data$RD) <- paste0(data$RD$chr, ":", data$RD$start, "-", data$RD$end)
  data$RD$chrN <- unlist(cs$chrom2chr[data$RD$chr])
  smo <- round(nrow(data$RD)/200)
  data$RD$L2R <- log2((data$RD$RD.test.mean+1) / (data$RD$RD.ref.mean+1))
  
  ## Handling low-depth bins
  message(tmsg("Imputing flagged bins L2r (low-depth, GC%-outiler) ..."))
  ### low depth
  # ldbins <- data$RD$RD_B < L2R.RD.min & data$RD$RD_A < L2R.RD.min
  ldbins <- data$RD$RD.test.mean < L2R.RD.min.Test | data$RD$RD.ref.mean < L2R.RD.min.Ref
  message(tmsg(paste0(" Found ", length(which(ldbins)), " low-depth bin(s).")))
  
  ### GC outlier
  gcobins <- GC.data$GC[,5] < .2 | GC.data$GC[,5] > .8
  message(tmsg(paste0(" Found ", length(which(gcobins)), " GC%-outlier bin(s).")))
  ### Pooled
  fbins <- ldbins + gcobins > 0
  
  if (any(ldbins)) {
    message(tmsg(paste0(" Imputed ", length(which(fbins)), " bins L2R.")))
    l2r.tmp <- data$RD$L2R
    l2r.tmp[fbins] <- NA
    data$RD$L2Ri <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  } else data$RD$L2Ri <- data$RD$L2R
  
  ## Pseudo-centering on the median
  data$RD[,ncol(data$RD)] <- data$RD[,ncol(data$RD)] - median(data$RD[,ncol(data$RD)], na.rm = TRUE)
  
  ## Normalizing L2R using GC content
  message(tmsg("GC-normalization ..."))
  rownames(GC.data$GC) <- paste0(GC.data$GC$chr, ":", GC.data$GC$start, "-", GC.data$GC$end)
  if (any(rownames(GC.data$GC) != rownames(data$RD))) stop(tmsg("GC data and L2R data are not synched, or ordered differently !"))
  ndata <- data.frame(data$RD[,1:3], name = rownames(data$RD), GC.data$GC[,-c(1:4)], stringsAsFactors = FALSE)
  my.rm.mad <- sum(abs(diff(as.numeric(runmed(data$RD[!is.na(data$RD$L2Ri),ncol(data$RD)], smo)))))
  normloop.res <- l2r.fitloop(l2rObj = list(l2r=data$RD[,ncol(data$RD)], rm.mad = my.rm.mad), tfd = ndata, smo = smo)
  rm(ndata)
  data$RD$L2Rrn <- normloop.res$l2r$l2r + median(data$RD[,ncol(data$RD)], na.rm = TRUE)
  meta.b$renorm <- if(is.null(normloop.res$pos)) "None" else paste0(normloop.res$pos, collapse = ",")
  
  ## Preparing BAF
  ###
  message(tmsg("Normalizing BAF data ..."))
  ldbaf <- data$SNP$RD.test < BAF.RD.min
  BT <- length(which(ldbaf)) / nrow(data$SNP)
  if (any(ldbaf)) {
    data$SNP <- data$SNP[!ldbaf,]
    message(tmsg(paste0(" Removed ", length(which(ldbaf)), " low-depth SNPs (", round(BT*100, digits = 2), " %).")))
  }
  # data$SNP$BAF <- data$SNP$BAF
  data$SNP$BAF[data$SNP$BAF < 0] <- 0
  data$SNP$BAF[data$SNP$BAF > 1] <- 1
  odd.idx <- which(data$SNP$pos %% 2 == 1)
  data$SNP$BAF[odd.idx] <- -data$SNP$BAF[odd.idx] +1
  
  ## Building ASCAT-like object
  message(tmsg("Building normalized object ..."))
  SNPz <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = data$SNP$chr, start = data$SNP$pos, end = data$SNP$pos, stringsAsFactors = FALSE))
  L2Rz <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = data$RD$chr, start = data$RD$start, end = data$RD$end, stringsAsFactors = FALSE))
  SinL <- GenomicRanges::findOverlaps(SNPz, L2Rz)
  data$SNP$L2R <- NA
  data$SNP$L2R[SinL@from] <- data$RD[, ncol(data$RD)][SinL@to]
  
  l2r.ao.df <- data.frame(chr = data$RD$chr, pos = round((data$RD$start + data$RD$end) / 2), L2R.ori = data$RD$L2R, L2R = data$RD[,ncol(data$RD)], BAF = NA, stringsAsFactors = FALSE)
  rownames(l2r.ao.df) <- paste0(l2r.ao.df$chr, ":", l2r.ao.df$pos)
  baf.ao.df <- data.frame(chr = data$SNP$chr, pos = data$SNP$pos, L2R.ori = data$SNP$L2R, L2R = data$SNP$L2R, BAF = data$SNP$BAF, stringsAsFactors = FALSE)
  rownames(baf.ao.df) <- rownames(data$SNP)
  ao.df <- rbind(l2r.ao.df, baf.ao.df)
  ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
  ao.df$chrs <- sub(pattern = "chr", replacement = "", ao.df$chr)
  ao.df <- ao.df[order(ao.df$chrN, ao.df$pos),]
  
  my.ch <- sapply(unique(ao.df$chrs), function(x) { which(ao.df$chrs == x) })
  
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR = data.frame(sample = ao.df$L2R),
      Tumor_BAF = data.frame(sample = ao.df$BAF),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      SNPpos = data.frame(chrs = ao.df$chrs, pos = ao.df$pos),
      ch = my.ch,
      chr = my.ch,
      chrs = unique(ao.df$chrs),
      samples = samplename,
      gender = "NA",
      sexchromosomes = c("X", "Y"),
      failedarrays = NULL
    ),
    meta = list(
      basic = meta.b,
      wes = meta.w
    )
  )
  colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  
  ## Saving normalized object
  message(tmsg("Saving normalized / prepared data ..."))
  # new.rds.name <- sub(pattern = "\\.RDS$", replacement = "_normalized.RDS", x = data.file)
  # saveRDS(my.ascat.obj, file = new.rds.name, compress = "bzip2")
  saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "_", data$meta$basic$genome, "_b", data$meta$WES$bin.size, "_processed.RDS"), compress = "xz")
  
  ## Rough plot
  message(tmsg("Plotting ..."))
  l2r.notna <- which(!is.na(ao.df$L2R))
  l2r.rm <- runmed(ao.df$L2R[l2r.notna], smo)
  l2r.mad <- median(abs(diff(ao.df$L2R[l2r.notna])))
  l2r.ssad <- sum(abs(diff(l2r.rm)))
  l2r.ori.rm <- runmed(ao.df$L2R.ori[l2r.notna], smo)
  l2r.ori.mad <- median(abs(diff(ao.df$L2R.ori[l2r.notna])))
  l2r.ori.ssad <- sum(abs(diff(l2r.ori.rm)))
  ao.df$genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
  ao.df$L2R.ori <- ao.df$L2R.ori - median(ao.df$L2R.ori, na.rm = TRUE)
  kend <- ao.df$genopos[vapply(unique(ao.df$chrN), function(k) { max(which(ao.df$chrN == k))}, 1)]
  png(paste0(out.dir, "/", samplename, "_WES_", data$meta$basic$genome, "_rawplot.png"), 1600, 1050)
  # par(mfrow = c(2,1))
  par(mfrow = c(3,1))
  plot(ao.df$genopos, ao.df$L2R.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") raw L2R profile (median-centered)\nMAD = ", round(l2r.ori.mad, digits = 2), " ; SSAD = ", round(l2r.ori.ssad, digits = 2)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.ori.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$L2R, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") normalized L2R profile (median-centered)\nMAD = ", round(l2r.mad, digits = 2), " ; SSAD = ", round(l2r.ssad, digits = 2)), xlab = "Genomic position", ylab = "L2R")
  lines(ao.df$genopos[l2r.notna], l2r.rm, col = 1)
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = 0, col = 2, lty = 2, lwd = 2)
  plot(ao.df$genopos, ao.df$BAF, pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") BAF profile"), xlab = "Genomic position", ylab = "BAF")
  abline(v = kend, col = 4, lty = 3, lwd = 2)
  abline(h = .5, col = 2, lty = 2, lwd = 2)
  dev.off()
  
  message("Done.")
  if(return.data) return(my.ascat.obj)
}

## Runs EaCoN.WES.Normalize using a RDS filename
EaCoN.WES.Normalize.ff <- function(BIN.RDS.file = NULL, ...) {
  if (is.null(BIN.RDS.file)) stop(tmsg("An RDS file from EaCoN::EaCoN.WES.Bin is required !"))
  if (!file.exists(BIN.RDS.file)) stop(tmsg(paste0("Could not find ", BIN.RDS.file, " .")))
  message(tmsg("Loading binned WES data ..."))
  my.data <- readRDS(BIN.RDS.file)
  EaCoN.WES.Normalize(data = my.data, out.dir = dirname(BIN.RDS.file), ...)
}

EaCoN.WES.Normalize.ff.Batch <- function(BIN.RDS.files = list.files(path = getwd(), pattern = "_binned.RDS$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE), nthread = 1, cluster.type = "PSOCK", ...) {
  if (length(BIN.RDS.files) == 0) stop("No file found to process !")
  message("Running EaCoN.WES.Normalize.ff() in batch mode ...")
  message(paste0("Found ", length(BIN.RDS.files), " samples to process ..."))
  current.bitmapType <- getOption("bitmapType")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  eacon.batchres <- foreach::foreach(r = seq_along(BIN.RDS.files), .inorder = TRUE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    EaCoN.WES.Normalize.ff(BIN.RDS.file = BIN.RDS.files[r], ...)
  }
  parallel::stopCluster(cl)
}

bedBinner <- function(bed = NULL, bin.size = 50, nthread = 1) {
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK", outfile = "")
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::"%dopar%"
  bed.binned <- foreach::foreach(b = seq_len(nrow(bed)), .combine = "rbind") %dopar% {
    # print(b)
    ### Smaller exon
    exon.length <- (bed$end[b] - bed$start[b] + 1)
    if (exon.length <= bin.size) return(bed[b,])
    mod.rest <- exon.length %% bin.size
    mod.count <- (exon.length - mod.rest) / bin.size

    bin.starts <- bed$start[b] + ((seq_len(mod.count)-1) * bin.size)
    bin.ends <- bin.starts + bin.size - 1

    ## Non-Round count
    if (mod.rest > 0) {
      ## Enough for a new bin
      if (mod.rest >= (bin.size / 2)) {
        bin.starts <- c(bin.starts, bin.starts[mod.count]+bin.size)
        bin.ends <- c(bin.ends, bed$end[b])
        mod.count <- mod.count+1
      } else { ## Dispatch to inside bins
        if (mod.rest >= mod.count) {
          mod.rest2 <- mod.rest %% mod.count
          mod.count2 <- (mod.rest - mod.rest2) / mod.count

          bin.starts <- bed$start[b] + ((seq_len(mod.count)-1) * (bin.size + mod.count2))
          bin.ends <- bin.starts + (bin.size + mod.count2) - 1
          bin.ends[mod.count] <- bin.ends[mod.count] + mod.rest2
        } else {
          # bin.ends[mod.count] <- bin.ends[mod.count] + mod.rest
          bin.ends[mod.count] <- bed$end[b]
        }
      }
    }
    # if(length(bin.starts) != length(bin.starts))
    chrs = rep(bed$chr[b], mod.count)
    return(data.frame(chr = chrs, start = bin.starts, end = bin.ends, stringsAsFactors = FALSE))
  }
  parallel::stopCluster(cl)
  return(bed.binned)
}

## Compute letter composition of nucleotidic sequences from a (chr, start, end) dataframe, with possible extension.
loc.nt.count.hs <- function(loc.df = NULL, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", extend = 0, blocksize = 1E+04, nthread = 5) {
  if (is.null(loc.df)) stop("loc.df is required !")
  # valid.genomes <- list("hg18" = "BSgenome.Hsapiens.UCSC.hg18", "hg19" = "BSgenome.Hsapiens.UCSC.hg19", "hg38" = "BSgenome.Hsapiens.UCSC.hg38", "GRCh37-lite" = "BSgenome.Hsapiens.TCGA.GRCh37.lite")

  if (extend < 0) stop("extend should be >= 0")
  if (blocksize <= 0) stop("blocksize should be > 0")
  if (!all(is.character(loc.df$chr))) stop("chr should be character !")
  if (!all(is.numeric(loc.df$start))) stop("start should be numeric !")
  if (!all(is.numeric(loc.df$end))) stop("end should be numeric !")
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }

  print(paste0("Loading ", genome.pkg, " sequence ..."))
  require(genome.pkg, character.only = TRUE)
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  # requireNamespace(GenomicRanges, quietly = TRUE)
  genome <- BSgenome::providerVersion(BSg.obj)
  # data(list = genome, package = "chromosomes", envir = environment())
  cs <- chromobjector(BSg.obj)


  print("Removing replicated locations ...")
  idz <- paste0(loc.df$chr, ":", loc.df$start, "-", loc.df$end)
  loc.df <- loc.df[!duplicated(idz),]

  print("Removing non-canonical sequences ...")
  loc.df <- loc.df[loc.df$chr %in% seqnames(BSg.obj),]

  print("Ordering data ...")
  loc.df <- loc.df[order(unlist(cs$chrom2chr[loc.df$chr]), loc.df$start, loc.df$ProbeSetName),]


  myGR.ex <- GenomicRanges::makeGRangesFromDataFrame(loc.df, seqinfo = seqinfo(BSg.obj))
  if (extend > 0) myGR.ex <- GenomicRanges::trim(myGR.ex + extend)

  instep <- as.numeric(cut(seq_along(myGR.ex), seq.int(0, length(myGR.ex) + blocksize, blocksize)))

  print("Computing base composition ...")

  print("Starting cluster ...")
  if (length(unique(instep)) < nthread) nthread <- length(unique(instep))
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK", outfile = "")
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::"%dopar%"
  xcounts <- foreach::foreach(x = unique(instep), .combine = "rbind", .packages = c("Biostrings", "BSgenome"), .export = c("BSg.obj", "getSeq", "myGR.ex", "instep")) %dopar% {
    return(Biostrings::alphabetFrequency(BSgenome::getSeq(BSg.obj, myGR.ex[which(instep == x)]), baseOnly = TRUE))
  }
  print("Stopping cluster ...")
  parallel::stopCluster(cl)

  # rownames(xcounts) <- paste0(loc.df$chr, ":", loc.df$start, "-", loc.df$end)
  out.df <- cbind(loc.df, xcounts)

  return(out.df)
}

loc.nt.gcc.hs <- function(loc.counts = NULL) {
  gcc <- (loc.counts$C + loc.counts$G) / (loc.counts$A + loc.counts$C + loc.counts$G + loc.counts$T)
  return(data.frame(loc.counts, GC = gcc, stringsAsFactors = FALSE))
}

## Compute GC on a (chr, start, end) dataframe using multiple extend values
loc.nt.gcc.hs.multi <- function(loc.df = NULL, extend.multi = c(50, 100, 200, 400, 800, 1600, 3200, 6400), ...) {
  # require(foreach)
  requireNamespace(foreach, quietly = TRUE)
  gc.list <- foreach(nt.add = extend.multi) %do% {
    print(paste0("Computing GC +", nt.add, " ..."))
    adb.counts <- loc.nt.count.hs(loc.df = loc.df, extend = nt.add, ...)
    adb.gc <- loc.nt.gcc.hs(loc.counts = adb.counts)
    return(adb.gc)
  }
  base.df <- gc.list[[1]][,1:4]
  gc.df <- foreach(nt.add = 1:length(gc.list), .combine = "cbind") %do% { return(gc.list[[nt.add]][["GC"]]) }
  colnames(gc.df) <- paste0("GC", extend.multi)
  return(data.frame(base.df, gc.df, stringsAsFactors = FALSE))
}

genome.build.finder <- function(BAM.header = NULL, valid.genomes = NULL) {
  BAM.header <- unlist(BAM.header)
  query <- paste0("(", paste0(valid.genomes, collapse = "|"), ")")
  stvh.grep <- grep(query, unlist(BAM.header))
  if (length(stvh.grep) == 0) stop(tmsg("Could not automatically determine genome build ! Please specify it !"))

  stvh.regexec <- unique(vapply(stvh.grep, function(x) {
    rc.res <- regexec(query, BAM.header[x])[[1]]
    return(as.character(substr(BAM.header[x], start = rc.res[1], stop = rc.res[1]+attr(rc.res, "match.length")[1]-1)))
  }, "a"))

  ok.genome <- unique(stvh.regexec[stvh.regexec %in% valid.genomes])
  if (length(ok.genome) == 0) stop(tmsg(paste0("Identified a putative genome build (", ok.genome, "), but not a supported one !")))
  if (length(ok.genome) >= 2) stop(tmsg(paste0("Identified more than one putative genome build (", ok.genome, ") !")))
  return(ok.genome)
}

