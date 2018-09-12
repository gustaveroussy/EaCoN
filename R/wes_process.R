## Generates a BINpack (bins and pre-computed GC tracks) from a capture bed
BINpack.Maker <- function(bed.file = NULL, bin.size = 50, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", extend.multi = c(0, 50, 100, 200, 400, 800, 1600, 3200, 6400), blocksize = 1E+04, nthread = 1, out.dir = getwd(), return.data = FALSE) {

  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/test/")
  # bed.file <- "SureSelect_ClinicalResearchExome.padded_GRCh37-lite_merged_sorted.bed"
  # # bed.file = "test1.bed"
  # bin.size = 50
  # genome.pkg = "BSgenome.Hsapiens.1000genomes.hs37d5"
  # extend.multi = c(0, 100, 500, 1000, 2000, 4000)
  # blocksize = 1E+04
  # nthread = 5
  # out.dir = getwd()
  # return.data = FALSE
  # source("~/git_gustaveroussy/EaCoN/R/BED_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/wes_process.R")
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  
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
  
  bin.size <- as.integer(bin.size)
  extend.multi <- as.integer(extend.multi)
  
  ### Loading genome
  message(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  organism <- BSgenome::organism(BSg.obj)


  ## Cleaning BED
  bed.clean <- BedCheck(bed.file = bed.file, genome.pkg = genome.pkg, out.dir = out.dir, return.data = TRUE)

  ## Binning
  message(paste0("Performing binning (", bin.size, ") ..."))
  bed.binned <- bedBinner(bed = bed.clean, bin.size = bin.size, nthread = nthread)

  bed.binned <- data.frame(ProbeSetName = seq_len(nrow(bed.binned)), bed.binned, stringsAsFactors = TRUE)

  message("Generating GC% tracks ...")
  wes.gc <- loc.nt.gcc.hs.multi(loc.df = bed.binned, genome.pkg = genome.pkg, extend.multi = extend.multi, blocksize = blocksize, nthread = nthread)

  rm(bed.binned)
  
  
  wes.meta.key <- c("genome-species", "genome-version", "genome-package", "array_type", "track_type", "bin_size")
  wes.meta.value <- c(organism, genome, genome.pkg, "WES", "GC", bin.size)
  # renorm.data <- list(tracks = wes.gc, info = list("genome-version" = genome, "genome-package" = genome.pkg, bin.size = bin.size, track.type = "GC"), bed.clean = bed.clean, bed.binned = bed.binned)
  renorm.data <- list(tracks = wes.gc, info = data.frame(key = wes.meta.key, value = wes.meta.value, stringsAsFactors = FALSE), bed.clean = bed.clean)
  rm(wes.gc, wes.meta.value, wes.meta.key, bed.clean)
  
  save("renorm.data", file = paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = paste0("_", genome, "_b", bin.size, ".GC.rda"), x = basename(bed.file), ignore.case = TRUE)), compress = "xz")

  message("Done.")
  if (return.data) return(renorm.data)
}


WES.Bin <- function(testBAM = NULL, refBAM = NULL, BINpack = NULL, samplename = "SAMPLE", Q = 20, nsubthread = 1, cluster.type = "PSOCK", out.dir = getwd(), return.data  = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  # setwd("/home/job/WORKSPACE/EaCoN_tests/WES/AlexandreLefranc/Results")
  # # setwd("/home/job/WORKSPACE/EaCoN_tests")
  # testBAM <- "/home/job/WORKSPACE/EaCoN_tests/WES/AlexandreLefranc/DATA/Sample_PHEO_AG_HS_048_DNA.reord.sorted.dedup.recal.reheaded.bam"
  # refBAM <- "/home/job/WORKSPACE/EaCoN_tests/WES/AlexandreLefranc/DATA/Sample_PHEO_AG_HS_048G_DNA.reord.sorted.dedup.recal.reheaded.bam"
  # # # # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/MATCHR/RESOURCES/SureSelect_ClinicalResearchExome.padded_hg19_b50.rda"
  # # # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_hs37d5_b50.GC.rda"
  # BINpack <- "/home/job/WORKSPACE/EaCoN_tests/WES/AlexandreLefranc/Results/V4-UTRs.hg38.fragment_targets_minimal_sorted_longChr_hg38_b50.GC.rda"
  # samplename <- "ALEX1"
  # Q <- 20
  # out.dir = getwd()
  # nsubthread = 4
  # cluster.type = "PSOCK"
  # return.data = FALSE
  # write.data = TRUE
  # plot = TRUE
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("/home/job/git_gustaveroussy/EaCoN/R/BED_functions.R")
  # source("/home/job/git_gustaveroussy/EaCoN/R/wes_process.R")
  # suppressPackageStartupMessages(require(foreach))
  # require(magrittr)
  
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
  if (!return.data & !write.data) stop(tmsg("Data should be returned and/or written on disk, but not none !"))
  if (Q < 0) stop(tmsg("Q should be positive !"))
  
  ## WARNINGS
  if (return.data) tmsg("Data will be returned.")
  if (write.data) tmsg("Data will be written on disk.")
  if (!plot) tmsg("No plot will get drawn.")
  
  
  ## Loading binpack
  tmsg("Loading BINpack ...")
  load(BINpack)
  
  ## CHECKS (genome)
  # genome.pkg <- GC.data$info$genome.pkg
  genome.pkg <- renorm.data$info$value[renorm.data$info$key == "genome-package"]
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !'))) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  
  ## Files controls
  tmsg("Checking BINpack and BAMs compatibility ...")
  
  ## Inspecting BAM headers
  refBAM.h <- Rsamtools::scanBamHeader(refBAM)
  testBAM.h <- Rsamtools::scanBamHeader(testBAM)
  bed.data <- renorm.data$tracks[,1:4]
  renorm.data$tracks <- NULL
  gc()
  
  if (!all(names(refBAM.h[[1]]$targets) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and Test BAM are not compatible (different chr names) !"))
  if (!all(bed.data$chr %in% names(refBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and BED are not compatible (different chr names) !"))
  if (!all(unique(bed.data$chr) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Test BAM and BED are not compatible (different chr names) !"))
  if (!all(unique(bed.data$chr) %in% BSgenome::seqnames(BSg.obj))) stop(tmsg("BED and BSgenome are not compatible (different chr names) !"))
  
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
    # bin.size = GC.data$info$bin.size
    bin.size = as.numeric(renorm.data$info$value[renorm.data$info$key == "bin_size"])
  )
  rm(testBAM.h, refBAM.h)
  
  ### Common BAM flags
  param.FLAG <- Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE, isDuplicate=FALSE)
  param.PILEUP <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = TRUE)
  
  bedbam2pup <- function(BamFile = NULL, bed.data = NULL, scanBamFlag = NULL, pileupParam = NULL) {
    
    ### Making pileup
    tmsg("   Getting pileup ...")
    param.BAM <- Rsamtools::ScanBamParam(which = GenomicRanges::makeGRangesFromDataFrame(bed.data, seqnames.field = "chr"), flag = scanBamFlag)
    pres <- dplyr::as.tbl(Rsamtools::pileup(BamFile, scanBamParam = param.BAM, pileupParam = pileupParam))
    colnames(pres)[c(1,5)] <- c("chr", "bin")
    levels(pres$bin) <- bed.data$ProbeSetName
    pres$bin <- as.integer(as.character(pres$bin))
    
    # levels(pres$bin) <- rownames(bed.data)
    
    gc()
    
    if("=" %in% unique(pres$nucleotide)) tmsg("   Bam contains the '=' sign !")
    
    ### Building genomic sequence of reference block
    tmsg("   Getting reference genome sequence ...")
    refblock <- pres[!duplicated(pres$pos), c(1,2)]
    pres <- dplyr::group_by(pres, pos)
    refblock$tot_count <- dplyr::summarize(pres, tot_count = sum(count))$tot_count
    pres <- dplyr::ungroup(pres)
    refblock$nucleotide <- as.factor(BSgenome::getSeq(BSg.obj, names = GenomicRanges::makeGRangesFromDataFrame(refblock, start.field = "pos", end.field = "pos"), as.character = TRUE))
    
    ### Merging blocks
    tmsg("   Computing alt counts ...")
    # refblock <- dplyr::group_by(refblock, pos, nucleotide)
    # pres <- dplyr::group_by(pres, pos, nucleotide)
    merged <- suppressWarnings(dplyr::left_join(refblock, pres, by = c("chr", "pos", "nucleotide")))
    rm(pres, refblock)
    gc()
    # bed.based <- dplyr::as.tbl(data.frame(chr = unique(bed.data$chr), pos = as.integer(unlist(seq.int2(from = bed.data$start, to = bed.data$end, by = 1))), which_label = rep(paste0(bed.data$chr, ":", bed.data$start, "-", bed.data$end), times = (bed.data$end - bed.data$start +1))))
    # bed.based <- dplyr::as.tbl(data.frame(chr = unique(bed.data$chr), pos = as.integer(unlist(seq.int2(from = bed.data$start, to = bed.data$end, by = 1))), bin = rep(rownames(bed.data), times = (bed.data$end - bed.data$start +1))))
    bed.based <- dplyr::as.tbl(data.frame(chr = unique(bed.data$chr), pos = as.integer(unlist(seq.int2(from = bed.data$start, to = bed.data$end, by = 1))), bin = rep(bed.data$ProbeSetName, times = (bed.data$end - bed.data$start +1))))
    bed.joint <- suppressWarnings(dplyr::left_join(bed.based, merged, c("chr", "pos", "bin"))) ### yeah !
    rm(bed.based)
    gc()
    ### Cleaning and formatting
    bed.joint$nucleotide <- NULL
    bed.joint$chr <- as.factor(bed.joint$chr)
    # bed.joint$which_label <- as.factor(bed.joint$which_label)
    
    bed.joint$count <- bed.joint$tot_count - bed.joint$count
    colnames(bed.joint) <- c("chr", "pos", "bin", "tot_count", "alt_count")
    
    bed.joint$tot_count[is.na(bed.joint$tot_count)] <- 0L
    bed.joint$alt_count[is.na(bed.joint$alt_count)] <- 0L
    # bed.joint$nt[is.na(bed.joint$nt)] <- bed.based$nt[is.na(bed.joint$nt)]

    return(bed.joint)
  }
  
  pileup.go <- function(testBAM = NULL, refBAM = NULL, bed.data = NULL, scanBamFlag = NULL, pileupParam = NULL, nsubthread = 1, cluster.type = "PSOCK") {
    
    # BamFile <- testBAM
    # scanBamFlag <- param.FLAG
    # pileupParam <- param.PILEUP
    # nsubthread <- 1
    # cluster.type = "PSOCK"

    #### Indexing BAM if needed
    if (!file.exists(paste0(testBAM, ".bai"))) {
      tmsg("Indexing Test BAM ...")
      Rsamtools::indexBam(testBAM) 
    } else tmsg("Test BAM is already indexed.")
    if (!file.exists(paste0(refBAM, ".bai"))) {
      tmsg("Indexing Ref BAM ...")
      Rsamtools::indexBam(refBAM)
    } else tmsg("Ref BAM is already indexed.")
    
    #### Opening BAM connections
    # openBAM <- Rsamtools::BamFileList(c(testBAM, refBAM))
    
    #### Launching cluster
    if (length(unique(bed.data$chr)) < nsubthread) nsubthread <- length(unique(bed.data$chr))
    cl <- parallel::makeCluster(spec = nsubthread, type = cluster.type, outfile = "")
    doParallel::registerDoParallel(cl)
    k <- 0
    BAMcounts <- foreach::foreach(k = unique(as.character(bed.data$chr)), .inorder = TRUE, .export = c("tmsg", "BSg.obj", "bedbam2pup", "seq.int2")) %dopar% {
      
      tmsg(paste0(" Sequence : ", k))

      ### Computing counts for both BAMs
      bed.data.k <- bed.data[bed.data$chr == k,]
      
      ## Opening BAM connections
      openBAM <- Rsamtools::BamFileList(c(testBAM, refBAM))
      
      tmsg("  Computing TEST counts ...")
      testPUP <- bedbam2pup(BamFile = openBAM[[1]], bed.data = bed.data.k, scanBamFlag = scanBamFlag, pileupParam = pileupParam)
      tmsg("  Computing REF counts ...")
      refPUP <- bedbam2pup(BamFile = openBAM[[2]], bed.data = bed.data.k, scanBamFlag = scanBamFlag, pileupParam = pileupParam)
      
      ### Merge counts
      tmsg("  Merging TEST and REF ...")
      mPUP <- dplyr::inner_join(testPUP, refPUP, c("chr", "pos", "bin"))
      rm(testPUP, refPUP)
      gc()
      mPUP <- dplyr::group_by(mPUP, bin)
      
      ### Adding missing nucleotides
      # mPUP$nt <- as.factor(unlist(strsplit(as.character(BSgenome::getSeq(BSg.obj, names = GenomicRanges::makeGRangesFromDataFrame(mPUP, seqnames.field = "chr", start.field = "pos", end.field = "pos"))), split = "")))

      ### Binning
      tmsg("  Binning ...")
      # CN.table <- dplyr::summarize(mPUP, chr = unique(chr), start = as.integer(min(pos)), end = as.integer(max(pos)), tot_count.test = as.integer(round(mean(tot_count.x))), tot_count.ref = as.integer(round(mean(tot_count.y))), GCPC = as.integer(round(seqinr::GC(as.character(nt))*100)))
      CN.table <- dplyr::summarize(mPUP, chr = unique(chr), start = as.integer(min(pos)), end = as.integer(max(pos)), tot_count.test = as.integer(round(mean(tot_count.x))), tot_count.ref = as.integer(round(mean(tot_count.y))))
      mPUP <- dplyr::ungroup(mPUP)
      SNP.table <- mPUP[mPUP$alt_count.x > 0 | mPUP$alt_count.y > 0,]
      
      
      ### Cleaning
      rm(mPUP)
      gc()
      
      CN.table <- dplyr::arrange(CN.table, chr, start, end)
      # CN.table <- dplyr::select(CN.table, chr, start, end, bin, tot_count.test, tot_count.ref, GCPC)
      CN.table <- dplyr::select(CN.table, chr, start, end, bin, tot_count.test, tot_count.ref)
      # SNP.table$nt <- NULL
      colnames(SNP.table) <- c("chr", "pos", "bin", "tot_count.test", "alt_count.test", "tot_count.ref", "alt_count.ref")
      
      gc()
      return(list(CN = CN.table, SNP = SNP.table))
    }
    stopCluster(cl)
    return(BAMcounts)
  }
  
  COUNTS.all <- pileup.go(testBAM = testBAM, refBAM = refBAM, bed.data = bed.data, scanBamFlag = param.FLAG, pileupParam = param.PILEUP, nsubthread = nsubthread, cluster.type = cluster.type)
  
  CN.all <- foreach(k = seq_along(COUNTS.all), .combine = "rbind") %do% {
    k.tmp <- COUNTS.all[[k]]$CN
    COUNTS.all[[k]]$CN <- NULL
    gc()
    return(k.tmp)
  }
  
  SNP.all <- foreach(k = seq_along(COUNTS.all), .combine = "rbind") %do% {
    k.tmp <- COUNTS.all[[k]]$SNP
    COUNTS.all[[k]]$SNP <- NULL
    gc()
    return(k.tmp)
  }
  
  rm(COUNTS.all)
  gc()
  
  ## Building WESobj
  
  CN.all$bin <- as.integer(CN.all$bin)
  SNP.all$bin <- as.integer(SNP.all$bin)
  
  ### summaries (recoded function as R internal summary uses too much RAM !)
  my.summary <- function(myv = NULL) {
    vsum <- c(min(myv, na.rm = TRUE), quantile(myv, .25, na.rm = TRUE), median(myv, na.rm = TRUE), mean(myv, na.rm = TRUE), quantile(myv, .75, na.rm = TRUE), max(myv, na.rm = TRUE))
    names(vsum) <- c("min", "q25", "median", "mean", "q75", "max")
    return(vsum)
  }
  meta.w$BIN.tot.count.test.mean.summary <- my.summary(CN.all$tot_count.test[!is.na(CN.all$tot_count.test)])
  meta.w$BIN.tot.count.ref.mean.summary <- my.summary(CN.all$tot_count.ref[!is.na(CN.all$tot_count.ref)])
  meta.w$SNP.tot.count.test.summary <- my.summary(SNP.all$tot_count.test[!is.na(SNP.all$tot_count.test)])
  meta.w$SNP.tot.count.ref.summary <- my.summary(SNP.all$tot_count.ref[!is.na(SNP.all$tot_count.ref)])
  gc()
  
  WESobj <- list(RD = CN.all, SNP = SNP.all, meta = list(basic = meta.b, WES = meta.w))
  rm(CN.all, SNP.all)
  gc()
  

  ## QC : Computing coverages
  tmsg("Computing coverages ...")
  gw.rd <- sum(WESobj$RD$end - WESobj$RD$start +1)
  gw.snp <- nrow(WESobj$SNP)
  rd.cov <- data.frame(cuts = c(1, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200), stringsAsFactors = FALSE)
  rd.cov <- cbind(rd.cov, t(foreach::foreach(x = rd.cov$cuts, .combine = "cbind") %do% {
    test.rd.in <- WESobj$RD$tot_count.test >= x
    ref.rd.in <- WESobj$RD$tot_count.ref >= x
    test.snprd.in <- WESobj$SNP$tot_count.test >= x
    ref.snprd.in <- WESobj$SNP$tot_count.ref >= x
    test.cut.cov <- if(!any(test.rd.in)) NA else (sum(WESobj$RD$end[test.rd.in] - WESobj$RD$start[test.rd.in] +1)/gw.rd)
    ref.cut.cov <- if(!any(ref.rd.in)) NA else (sum(WESobj$RD$end[ref.rd.in] - WESobj$RD$start[ref.rd.in] +1)/gw.rd)
    test.snpcut.cov <- if(!any(test.snprd.in)) NA else (length(which(test.snprd.in))/gw.snp)
    ref.snpcut.cov <- if(!any(ref.snprd.in)) NA else (length(which(ref.snprd.in))/gw.snp)
    return(c(test.cut.cov, ref.cut.cov, test.snpcut.cov, ref.snpcut.cov))
  }))
  colnames(rd.cov) <- c("MinDepth", "TestBINCoverage", "RefBINCoverage", "TestBAFCoverage", "RefBAFCoverage")
  rm(test.snprd.in, ref.snprd.in, test.rd.in, ref.rd.in)
  
  if (write.data || plot) dir.create(paste0(out.dir, "/", samplename))
  if (write.data) write.table(rd.cov, file = paste0(out.dir, "/", samplename,  "/", samplename, '_WES_', genome, "_b", meta.w$bin.size, "_coverage.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  ## QC : Plotting coverages
  if (plot) {
    png(paste0(out.dir, "/", samplename, "/", samplename, "_WES_", genome, "_b", meta.w$bin.size, "_coverage.png"), 800, 640)
    plot(rd.cov$MinDepth, rd.cov$TestBAFCoverage, type = "b", col = 2, lty = 3, pch = 20, main = paste0(WESobj$meta$basic$samplename, "\nCoverage Plot"), xlab = "Minimum depth", ylab = "Coverage", ylim = c(0,1), xaxp = c(0,200,10))
    abline(v = rd.cov$MinDepth, lty = 2, col = "grey75")
    abline(h = seq(0,1,.1), lty = 2, col = "grey75")
    lines(rd.cov$MinDepth, rd.cov$RefBAFCoverage, type = "b", col = 1, lty = 3, pch = 20)
    lines(rd.cov$MinDepth, rd.cov$TestBINCoverage, type = "b", col = 2)
    lines(rd.cov$MinDepth, rd.cov$RefBINCoverage, type = "b", col = 1)
    abline(h = .5, lty = 2)
    legend("topright", legend = c("Test BAF", "Ref BAF", "Test CN", "Ref CN"), inset = .02, col = c(2,1,2,1), lty = c(3,3,1,1), pch = c(20,20,1,1))
    dev.off()
  }
  rm(rd.cov)
  
  ## Saving
  if (write.data) {
    tmsg("Saving counts data ...")
    saveRDS(WESobj, file = paste0(out.dir, "/", samplename, "/", samplename, "_", genome, "_b", meta.w$bin.size, "_binned.RDS"), compress = "bzip2")
  }
  if (return.data) return(WESobj)
}

## Performs the binning of BAMs using a BINpack, batch mode
WES.Bin.Batch <- function(BAM.list.file = NULL, BINpack = NULL, nthread = 1, cluster.type = "PSOCK", ...) {

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
  eacon.batchres <- foreach::foreach(r = seq_len(nrow(myBAMs)), .inorder = TRUE, .errorhandling = "stop", .export = c("EaCoN.set.bitmapType", "WES.Bin", "tmsg")) %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    WES.Bin(testBAM = myBAMs$testBAM[r], refBAM = myBAMs$refBAM[r], BINpack = BINpack, samplename = myBAMs$SampleName[r], cluster.type = cluster.type, ...)
  }
  parallel::stopCluster(cl)
}

## Performs the normalization of WES L2R and BAF signals
WES.Normalize <- function(data = NULL, BINpack = NULL, gc.renorm = TRUE, wave.renorm = FALSE, wave.rda = NULL, RD.tot.min = 20, RD.alt.min = 3, BAF.hetmin = .33, sex.chr = c("chrX", "chrY"), TumorBoost = FALSE, out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE) {

  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta2/WES/TCGA-A7-A0CE-01A_vs_10A")
  # setwd("/home/job/WORKSPACE/EaCoN_tests/WES/AlexandreLefranc")
  # # data <- readRDS("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta2/WES/TCGA-A7-A0CE-01A_vs_10A/TCGA-A7-A0CE-01A_vs_10A_hs37d5_b_binned.RDS")
  # # # # # # data <- readRDS("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_DEV/TCGA_A0CE_01.10/TCGA_A0CE_01.10_hs37d5_b50_binned.RDS")
  # # # # # # # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_hs37d5_b50.rda"
  # # # # # # # BINpack <- "/mnt/data_cigogne/job/OS2006/WES/RESOURCES/SureSelect_ClinicalResearchExome.padded_hg19_b50.rda"
  # # # data <- readRDS("TCGA-AC-A2BK-01A_vs_11A_hs37d5_b_binned.RDS")
  # data <- readRDS("Sample_PHEO_AG_HS_048_DNA_hg38_b50_binned.RDS")
  # # BINpack <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_GRCh37-lite_merged_sorted_hs37d5_b50.GC.rda"
  # BINpack <- "V4-UTRs.hg38.fragment_targets_minimal_sorted_longChr_hg38_b50.GC.rda"
  # gc.renorm <- TRUE
  # wave.renorm <- FALSE
  # # wave.rda <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/RESOURCES/SureSelect_ClinicalResearchExome.padded_GRCh37-lite_merged_sorted_hs37d5_b50.Wave.rda"
  # wave.rds <- NULL
  # RD.tot.min = 20
  # RD.alt.min = 3
  # TumorBoost = FALSE
  # # sex.chr <- c("X", "Y")
  # sex.chr <- c("chrX", "chrY")
  # out.dir = getwd()
  # return.data = FALSE
  # write.data = TRUE
  # plot = TRUE
  # BAF.hetmin <- .33
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("/home/job/git_gustaveroussy/EaCoN/R/renorm_functions.R")
  # require(foreach)

  ### AJOUTER UN CONTROLE DU RDS (pour que ce ne soit pas un _processed.RDS donné en entrée!)
  
  
  ## CHECKS
  if (!is.list(data)) stop(tmsg("data should be a list !"))
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"))
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"))
  if (wave.renorm) { if (!is.null(wave.rda)) { if (!file.exists(wave.rda)) stop(tmsg(paste0("Could not find wave.rda file ", wave.rda))) } }
  if (RD.tot.min < 0) stop(tmsg("RD.tot.min must be >= 0 !"))
  if (RD.alt.min <= 0) stop(tmsg("RD.alt.min must be > 0 !"))

  ## TAGS
  data$meta$WES$TumorBoost <- as.character(TumorBoost)
  data$meta$WES$RD.tot.min <- RD.tot.min
  data$meta$WES$GC.renorm <- as.character(gc.renorm)
  data$meta$WES$Wave.renorm <- as.character(wave.renorm)
  samplename <- data$meta$basic$samplename

  ## Loading BINpack
  load(BINpack)
  gc()

  genome.pkg <- renorm.data$info$value[renorm.data$info$key == "genome-package"]
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")))
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")))
    }
  }

  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  # requireNamespace(genome.pkg, quietly = TRUE)
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)

  ## BAF HANDLING
  tmsg("Processing BAF ...")
  
  BAF.adder <- function(Bdata = NULL, Bvalues = NULL, newname = NULL, type = "test") {
    Bdata[[paste0("BAF.", type)]] <- Bdata[[newname]] <- Bvalues
    Bdata[[paste0("mBAF.", type)]] <- BAF2mBAF(Bdata[[paste0("BAF.", type)]])
    return(Bdata)
  }

  ### BAF : Filtering (low depth)
  rd.ori <- nrow(data$SNP)
  #### 1) Filtering for low ref or test depth
  RDlow <- data$SNP$tot_count.test < RD.tot.min | data$SNP$tot_count.ref < RD.tot.min
  if (length(which(RDlow)) == nrow(data$SNP)) stop("All SNP positions were discarded for their low read count ! You may consider lowering the BAF.tot.min value.")
  tmsg(paste0("Removed ", length(which(RDlow)), " (", round(length(which(RDlow)) / rd.ori * 100, digits = 2), "%) SNP positions with low depth (<", RD.tot.min, ")"))
  data$SNP <- data$SNP[!RDlow,]
  #### 2) Filtering for low alt test depth
  BRDlow <- data$SNP$alt_count.test < RD.alt.min
  if (length(which(BRDlow)) == nrow(data$SNP)) stop("All SNP positions were discarded for their low alternative allele count ! You may consider lowering the RD.alt.min value.")
  tmsg(paste0("Removed ", length(which(BRDlow)), " (", round(length(which(BRDlow)) / rd.ori * 100, digits = 2), "%) SNP positions with low alt RD (<", RD.alt.min, ")"))
  data$SNP <- data$SNP[!BRDlow,]
  gc()

  ### Computing BAF
  data$SNP$BAF.test.ori <- data$SNP$alt_count.test / data$SNP$tot_count.test
  data$SNP$BAF.ref.ori <- data$SNP$alt_count.ref / data$SNP$tot_count.ref
  odd.idx <- which(data$SNP$pos %% 2 == 1)
  data$SNP$BAF.test.ori[odd.idx] <- -data$SNP$BAF.test.ori[odd.idx] +1L
  data$SNP$BAF.ref.ori[odd.idx] <- -data$SNP$BAF.ref.ori[odd.idx] +1L

  data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = data$SNP$BAF.test.ori, newname = "BAF.test.ori")
  data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = data$SNP$BAF.ref.ori, newname = "BAF.ref.ori", type = "ref")
  
  ###### Computing LOR (and variance)
  rcmat <- round(cbind(data$SNP$BAF.test.ori*data$SNP$tot_count.test, (1-data$SNP$BAF.test.ori)*data$SNP$tot_count.test))
  data$SNP$LOR <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
  data$SNP$LORvar <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
  rm(rcmat)
  gc()
  
  ### BAF : TumorBoost
  if (TumorBoost) {
    message(tmsg("Applying TumorBoost BAF normalization ..."))
    # data$SNP$BAF.test <- data$SNP$BAF.test.TB <- as.numeric(aroma.light::normalizeTumorBoost(data$SNP$BAF.test, data$SNP$BAF.ref, flavor = "v4", preserveScale = FALSE))
    BTB <- as.numeric(aroma.light::normalizeTumorBoost(data$SNP$BAF.test, data$SNP$BAF.ref, flavor = "v4", preserveScale = FALSE))
    data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = BTB, newname = "BAF.test.TB")
  }

  ### Getting heterozygous probes from Ref
  Ref.hetero <- data$SNP$mBAF.ref >= BAF.hetmin
  if (!any(Ref.hetero)) stop(tmsg("All SNP positions were tagged as homozygous in Ref : there may be a problem with your reference BAM ploidy !"))

  ### Keeping hetero positions
  data$SNP <- data$SNP[Ref.hetero,]
  ### Removing additional values per bin
  data$SNP <- data$SNP[!duplicated(data$SNP$bin),]
  gc()
  
  ### BAF filtering
  # smoB <- round(nrow(data$SNP) / 3300)
  # if(smoB%%2 == 0) smoB <- smoB+1
  # mBAF.rm <- runmed(data$SNP$mBAF.test, smoB)
  # mBAF.diff <- abs(data$SNP$mBAF.test - mBAF.rm)
  # Bfiltered <- mBAF.diff < quantile(mBAF.diff, BAF.filter)
  # data$SNP <- data$SNP[Bfiltered,]

  ### Adding LOR
  data$SNP$LOR.test <- log(data$SNP$alt_count.test / (data$SNP$tot_count.test - data$SNP$alt_count.test))

  ## L2R
  tmsg("Processing RD bins ...")

  #### Computing L2R
  data$RD$L2R <- data$RD$L2R.ori <- log2((data$RD$tot_count.test+1) / (data$RD$tot_count.ref+1))

  ### L2R : Filtering
  #### 1) low depth
  rd.ori <- nrow(data$RD)
  RDlow <- data$RD$tot_count.test < RD.tot.min | data$RD$tot_count.ref < RD.tot.min
  data$meta$WES$Imputed.lowdepth.bins <- length(which(RDlow))
  if (length(which(RDlow)) == nrow(data$RD)) stop("All RD bins were flagged for their low read count ! You may consider lowering the BAF.tot.min value.")
  tmsg(paste0("Flagged ", length(which(RDlow)), " (", round(length(which(RDlow)) / rd.ori * 100, digits = 2), "%) RD bins with low depth (<", RD.tot.min, ")"))

  #### 2) GC% outliers
  GCOL <- renorm.data$tracks[,5] < 200 | renorm.data$tracks[,5] > 800
  data$meta$WES$Imputed.GCoutlier.bins <- length(which(GCOL))
  if (length(which(GCOL)) == nrow(data$RD)) stop("All RD bins were flagged as GC% outliers  ! There may be something wrong with your reference genome and/or capture BED.")
  tmsg(paste0("Flagged ", length(which(GCOL)), " (", round(length(which(GCOL)) / rd.ori * 100, digits = 2), "%) RD bins as GC% outliers."))

  ### Pooling and imputing
  FLAGS <- RDlow + GCOL > 0

  if (any(FLAGS)) {
    tmsg(paste0(" Imputed ", length(which(FLAGS)), " (", round(length(which(FLAGS))/rd.ori*100, digits = 2), "%) L2R bins."))
    l2r.tmp <- data$RD$L2R
    l2r.tmp[FLAGS] <- NA
    data$RD$L2R <- data$RD$L2R.imp <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  } else data$RD$L2R.imp <- data$RD$L2R
  gc()

  ## L2R : Normalization
  smo <- round(nrow(data$RD) / 550)
  if(smo%%2 == 0) smo <- smo+1

  ### Wave
  if (wave.renorm) {
    tmsg("Wave normalization ...")

    l2r2norm <- data.frame(ProbeSetName = data$RD$bin, chr = as.character(data$RD$chr), pos = data$RD$start, L2R = data$RD$L2R)
    # rownames(l2r2norm) <- seq_len(nrow(l2r2norm))
    ren.res <- renorm.go(input.data = l2r2norm, renorm.rda = wave.rda, track.type = "Wave", smo = smo, arraytype = data$meta$basic$type, genome = genome)

    fitted.l2r <- ren.res$renorm$l2r$l2r

    GCF <- is.na(fitted.l2r)
    if (any(GCF)) {
      l2r.tmp <- fitted.l2r
      l2r.tmp[GCF] <- NA
      fitted.l2r <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
    }

    if(is.null(ren.res$renorm$pos)) {
      # meta.b <- setmeta("gc.renorm", "None", meta.b)
      data$meta$WES <- setmeta("wave.renorm", "None", data$meta$WES)
      tmsg(" No positive fit.")
    } else {
      ## Tweaking sex chromosomes
      sex.idx <- data$RD$chr %in% sex.chr
      auto.ori.med <- median(data$RD$L2R[!sex.idx], na.rm = TRUE)
      auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
      if (any(sex.idx)) {
        for (k in sex.chr) {
          k.idx <- data$RD$chr == k
          if (any(k.idx)) {
            k.ori.diffmed <- median(data$RD$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
            fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
          }
        }
      }
      # meta.b <- setmeta("wave.renorm", paste0(ren.res$renorm$pos, collapse = ","), meta.b)
      data$meta$WES <- setmeta("wave.renorm", paste0(ren.res$renorm$pos, collapse = ","), data$meta$WES)
    }
    rm(ren.res)

    data$RD$L2R.WAVE <- data$RD$L2R <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
  } else {
    # meta.b <- setmeta("wave.renorm", "FALSE", meta.b)
    data$meta$WES <- setmeta("wave.renorm", "FALSE", data$meta$WES)
  }

  #### GC%
  # message("GC% normalization ...")
  # data$RD$L2R.GC <- data$RD$L2R <- limma::loessFit(x = data$RD$GCPC, y = data$RD$L2R)$residuals
  # if (any(is.na(data$RD$L2R))) {
  #   l2r.tmp <- data$RD$L2R
  #   l2r.tmp[is.na(data$RD$L2R)] <- NA
  #   data$RD$L2R <- data$RD$L2R.GC <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  # }



  if (gc.renorm) {
    tmsg("GC% normalization ...")

    l2r2norm <- data.frame(ProbeSetName = data$RD$bin, chr = as.character(data$RD$chr), pos = data$RD$start, L2R = data$RD$L2R)
    # rownames(l2r2norm) <- seq_len(nrow(l2r2norm))
    ren.res <- renorm.go(input.data = l2r2norm, renorm.rda = BINpack, track.type = "GC", smo = smo, arraytype = data$meta$basic$type, genome = genome)

    fitted.l2r <- ren.res$renorm$l2r$l2r

    GCF <- is.na(fitted.l2r)
    if (any(GCF)) {
      l2r.tmp <- fitted.l2r
      l2r.tmp[GCF] <- NA
      fitted.l2r <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
    }

    if(is.null(ren.res$renorm$pos)) {
      # meta.b <- setmeta("gc.renorm", "None", meta.b)
      data$meta$eacon <- setmeta("gc.renorm", "None", data$meta$eacon)
      tmsg(" No positive fit.")
    } else {
      ## Tweaking sex chromosomes
      sex.idx <- data$RD$chr %in% sex.chr
      auto.ori.med <- median(data$RD$L2R.ori[!sex.idx], na.rm = TRUE)
      auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
      if (any(sex.idx)) {
        for (k in sex.chr) {
          k.idx <- data$RD$chr == k
          if (any(k.idx)) {
            # k.ori.diffmed <- median(data$RD$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.ori.diffmed <- median(data$RD$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
            k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
            fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
          }
        }
      }
      # meta.b <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, collapse = ","), meta.b)
      data$meta$eacon <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, collapse = ","), data$meta$eacon)
    }
    rm(ren.res)

    data$RD$L2R.GC <- data$RD$L2R <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
  } else {
    # meta.b <- setmeta("gc.renorm", "FALSE", meta.b)
    data$meta$eacon <- setmeta("gc.renorm", "FALSE", data$meta$eacon)
  }

  ## Merging
  data$RD$BAF <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "BAF.test"))], by = c("chr", "bin"))$BAF.test
  # data$RD$LOR <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LOR.test"))], by = c("chr", "bin"))$LOR.test
  data$RD$LOR <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LOR"))], by = c("chr", "bin"))$LOR
  data$RD$LORvar <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LORvar"))], by = c("chr", "bin"))$LORvar
  data$RD$RD.test <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "tot_count.test"))], by = c("chr", "bin"))$tot_count.test
  data$RD$RD.ref <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "tot_count.ref"))], by = c("chr", "bin"))$tot_count.ref
  
  ## Building ASCAT object
  tmsg("Building normalized object ...")

  my.ch <- sapply(unique(data$RD$chr), function(x) { which(data$RD$chr == x) })
  
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR.ori = data.frame(sample = data$RD$L2R.ori, row.names = data$RD$bin),
      Tumor_LogR = data.frame(sample = data$RD$L2R, row.names = data$RD$bin),
      Tumor_BAF = data.frame(sample = data$RD$BAF, row.names = data$RD$bin),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      SNPpos = data.frame(chrs = data$RD$chr, pos = round((data$RD$start + data$RD$end)/2)),
      ch = my.ch,
      chr = my.ch,
      chrs = levels(data$RD$chr),
      samples = samplename,
      gender = "NA",
      sexchromosomes = sex.chr,
      failedarrays = NULL,
      additional = data$RD[,colnames(data$RD) %in% c("RD.test", "RD.ref", "LOR", "LORvar")]
    ),
    meta = data$meta,
    germline = list(germlinegenotypes = matrix(is.na(data$RD$BAF), ncol = 1, dimnames = list(data$RD$bin, samplename)), failedarrays = NULL)
  )
  # colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- colnames(my.ascat.obj$data$Tumor_LOR) <- samplename
  colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  # rm(my.ch, data)
  gc()

  # plot(ares$Tumor_LogR[,1], pch = ".", cex = 3, xaxs = "i", ylim = c(-2,2))
  # points(ares$Tumor_LogR_segmented, pch = ".", cex = 3, col = 2)
  # plot(ares$Tumor_BAF[!is.na(ares$Tumor_BAF),1], pch = ".", xaxs = "i", cex = 3)
  # points(ares$Tumor_BAF_segmented[[1]], pch = ".", cex = 3, col = 2)
  # points(1 - ares$Tumor_BAF_segmented[[1]], pch = ".", cex = 3, col = 2)
  #
  ## Saving data
  if (write.data) {
    tmsg("Saving normalized data ...")
    saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "_", data$meta$basic$genome, "_b", data$meta$WES$bin.size, "_processed.RDS"), compress = "bzip2")
  }

  ## Plot
  tmsg("Plotting ...")
  if (plot) {
    l2r <- my.ascat.obj$data$Tumor_LogR[,1]
    l2r.rm <- runmed(l2r, smo)
    l2r.dif <- diff(l2r)
    l2r.mad <- median(abs(l2r.dif[l2r.dif != 0]))
    l2r.rm.dif <- diff(l2r.rm)
    l2r.ssad <- sum(abs(l2r.rm.dif[l2r.rm.dif != 0]))
    
    l2r.ori <- my.ascat.obj$data$Tumor_LogR.ori[,1]
    l2r.ori.rm <- runmed(l2r.ori, smo)
    l2r.ori.dif <- diff(l2r.ori)
    l2r.ori.mad <- median(abs(l2r.ori.dif[l2r.ori.dif != 0]))
    l2r.ori.rm.dif <- diff(l2r.ori.rm)
    l2r.ori.ssad <- sum(abs(l2r.ori.rm.dif[l2r.ori.rm.dif != 0]))
    
    l2r.genopos <- my.ascat.obj$data$SNPpos$pos + cs$chromosomes$chr.length.toadd[my.ascat.obj$data$SNPpos$chrs]
    
    l2r <- l2r - median(l2r, na.rm = TRUE)
    l2r.ori <- l2r.ori - median(l2r.ori, na.rm = TRUE)
    kend <- l2r.genopos[vapply(unique(my.ascat.obj$data$SNPpos$chr), function(k) { max(which(my.ascat.obj$data$SNPpos$chrs == k))}, 1)]
    
    png(paste0(out.dir, "/", samplename, "_WES_", data$meta$basic$genome, "_rawplot.png"), 1600, 1050)
    par(mfrow = c(3,1))
    plot(l2r.genopos, l2r.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") raw L2R profile (median-centered)\nMAD = ", round(l2r.ori.mad, digits = 2), " ; SSAD = ", round(l2r.ori.ssad, digits = 2)), xlab = "Genomic position", ylab = "L2R")
    lines(l2r.genopos, l2r.ori.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(l2r.genopos, l2r, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ") normalized L2R profile (median-centered)\nMAD = ", round(l2r.mad, digits = 2), " ; SSAD = ", round(l2r.ssad, digits = 2)), xlab = "Genomic position", ylab = "L2R")
    lines(l2r.genopos, l2r.rm, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(l2r.genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = "grey75", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " WES (", data$meta$basic$manufacturer, ")", if(TumorBoost) " TumorBoost-normalized", " BAF profile"), xlab = "Genomic position", ylab = "BAF")
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = .5, col = 2, lty = 2, lwd = 2)
    dev.off()
  }

  tmsg("Done.")
  if(return.data) return(my.ascat.obj)
}

## Runs WES.Normalize using a RDS filename
WES.Normalize.ff <- function(BIN.RDS.file = NULL, ...) {
  
  ## CHECKS
  if (is.null(BIN.RDS.file)) stop(tmsg("An RDS file from EaCoN::EaCoN.WES.Bin is required !"))
  if (!file.exists(BIN.RDS.file)) stop(tmsg(paste0("Could not find ", BIN.RDS.file, " .")))
  # if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"))
  # if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"))
  
  tmsg("Loading binned WES data ...")
  my.data <- readRDS(BIN.RDS.file)
  WES.Normalize(data = my.data, out.dir = dirname(BIN.RDS.file), ...)
}

## Runs WES.Normalize.ff, batch mode
WES.Normalize.ff.Batch <- function(BIN.RDS.files = list.files(path = getwd(), pattern = "_binned.RDS$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE), nthread = 1, cluster.type = "PSOCK", ...) {
  if (length(BIN.RDS.files) == 0) stop("No file found to process !")
  message("Running EaCoN.WES.Normalize.ff() in batch mode ...")
  message(paste0("Found ", length(BIN.RDS.files), " samples to process ..."))
  current.bitmapType <- getOption("bitmapType")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  eacon.batchres <- foreach::foreach(r = seq_along(BIN.RDS.files), .inorder = TRUE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    WES.Normalize.ff(BIN.RDS.file = BIN.RDS.files[r], ...)
  }
  parallel::stopCluster(cl)
}

bedBinner <- function(bed = NULL, bin.size = 50, nthread = 1) {
  
  bin.size <- as.integer(bin.size)
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK", outfile = "")
  suppressPackageStartupMessages(require(foreach))
  doParallel::registerDoParallel(cl)
  k <- 0
  bed.binned <- foreach::foreach(k = unique(bed$chr), .combine = "rbind", .export = "tmsg") %dopar% {
    tmsg(k)
    
    bedk <- bed[bed$chr == k,]
    
    b <- 0
    bbk <- foreach::foreach(b = seq_len(nrow(bedk)), .combine = "rbind", .export = "bin.size") %do% {
      
      ### Smaller exon
      exon.length <- (bedk$end[b] - bedk$start[b] + 1L)
      if (exon.length <= bin.size) return(bedk[b,])
      mod.rest <- exon.length %% bin.size
      mod.count <- as.integer((exon.length - mod.rest) / bin.size)
      
      bin.starts <- bedk$start[b] + ((seq_len(mod.count)-1L) * bin.size)
      bin.ends <- bin.starts + bin.size - 1L
      
      ## Non-Round count
      if (mod.rest > 0L) {
        ## Enough for a new bin
        if (mod.rest >= (bin.size / 2L)) {
          bin.starts <- c(bin.starts, bin.starts[mod.count]+bin.size)
          bin.ends <- c(bin.ends, bedk$end[b])
          mod.count <- mod.count+1L
        } else { ## Dispatch to inside bins
          if (mod.rest >= mod.count) {
            mod.rest2 <- mod.rest %% mod.count
            mod.count2 <- as.integer((mod.rest - mod.rest2) / mod.count)
            
            bin.starts <- bedk$start[b] + ((seq_len(mod.count)-1L) * (bin.size + mod.count2))
            bin.ends <- bin.starts + (bin.size + mod.count2) - 1L
            bin.ends[mod.count] <- bin.ends[mod.count] + mod.rest2
          } else {
            # bin.ends[mod.count] <- bin.ends[mod.count] + mod.rest
            bin.ends[mod.count] <- bedk$end[b]
          }
        }
      }
      # if(length(bin.starts) != length(bin.starts))
      chrs = rep(bedk$chr[b], mod.count)
      return(data.frame(chr = chrs, start = bin.starts, end = bin.ends, stringsAsFactors = FALSE))
    }
    
    return(bbk)
  }
  parallel::stopCluster(cl)
  return(bed.binned)
}

## Compute letter composition of nucleotidic sequences from a (chr, start, end) dataframe, with possible extension.
loc.nt.count.hs <- function(loc.df = NULL, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", extend = 0, blocksize = 1E+04, nthread = 5) {
  if (is.null(loc.df)) stop("loc.df is required !")
  
  if (extend < 0) stop("extend should be >= 0")
  if (blocksize <= 0) stop("blocksize should be > 0")
  if (!all(is.character(loc.df$chr) | is.factor(loc.df$chr))) stop("chr should be character !")
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
  # requireNamespace(genome.pkg, quietly = TRUE)
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  # require(genome.pkg, character.only = TRUE)
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)


  print("Removing replicated locations ...")
  idz <- paste0(loc.df$chr, ":", loc.df$start, "-", loc.df$end)
  loc.df <- loc.df[!duplicated(idz),]

  print("Removing non-canonical sequences ...")
  loc.df <- loc.df[loc.df$chr %in% seqnames(BSg.obj),]

  print("Ordering data ...")
  loc.df <- loc.df[order(unlist(cs$chrom2chr[loc.df$chr]), loc.df$start, loc.df$ProbeSetName),]


  
  myGR.ex <- suppressPackageStartupMessages(GenomicRanges::makeGRangesFromDataFrame(loc.df, seqinfo = seqinfo(BSg.obj)))
  if (extend > 0) myGR.ex <- GenomicRanges::trim(myGR.ex + extend)

  instep <- as.numeric(cut(seq_along(myGR.ex), seq.int(0, length(myGR.ex) + blocksize, blocksize)))

  print("Computing base composition ...")

  print("Starting cluster ...")
  
  if (length(unique(instep)) < nthread) nthread <- length(unique(instep))
  cl <- parallel::makeCluster(spec = nthread, type = "PSOCK", outfile = "")
  doParallel::registerDoParallel(cl)
  requireNamespace("foreach", quietly = TRUE)
  `%dopar%` <- foreach::"%dopar%"
  xcounts <- foreach::foreach(x = unique(instep), .combine = "rbind", .packages = c("Biostrings", "BSgenome"), .export = "getSeq") %dopar% {
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
  requireNamespace("foreach", quietly = TRUE)
  `%do%` <- foreach::"%do%"
  gc.list <- foreach::foreach(nt.add = extend.multi) %do% {
    print(paste0("Computing GC +", nt.add, " ..."))
    adb.counts <- loc.nt.count.hs(loc.df = loc.df, extend = nt.add, ...)
    adb.gc <- loc.nt.gcc.hs(loc.counts = adb.counts)
    return(adb.gc)
  }
  base.df <- gc.list[[1]][,1:4]
  gc.df <- foreach::foreach(nt.add = seq_along(gc.list), .combine = "cbind") %do% { return(as.integer(round(gc.list[[nt.add]][["GC"]] * 1000))) }
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

