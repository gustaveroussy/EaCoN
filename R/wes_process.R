## Generates a BINpack (bins and pre-computed GC tracks) from a capture bed
BINpack.Maker.old <- function(bed.file = NULL, bin.size = 50, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", extend.multi = c(0, 50, 100, 200, 400, 800, 1600, 3200, 6400), blocksize = 1E+04, nthread = 1, out.dir = getwd(), return.data = FALSE) {

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
  # source("~/git_gustaveroussy/EaCoN/R/hts_process.R")
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  
    ## Checks
  if (is.null(bed.file)) stop("A BED file is required !", call. = FALSE)
  if (!file.exists(bed.file)) stop("Could not find the BED file !", call. = FALSE)
  if (is.null(out.dir)) message("NOTE : Checked / cleaned bed will be written in the same directory as source.") else {
    if (!file.exists(out.dir)) stop("Could not find the output directory !", call. = FALSE)
    if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !", call. = FALSE)
  }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"), call. = FALSE)
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  
  bin.size <- as.integer(bin.size)
  extend.multi <- as.integer(extend.multi)
  
  ### Loading genome
  message(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  organism <- BSgenome::organism(BSg.obj)


  ## Cleaning BED
  bed.clean <- BedCheck(bed.file = bed.file, genome.pkg = genome.pkg, out.dir = out.dir, return.data = TRUE)

  ## Binning
  message(paste0("Performing binning (", bin.size, ") ..."))
  bed.binned <- bedBinner(bed = bed.clean, bin.size = bin.size, nthread = nthread)

  bed.binned <- data.frame(ProbeSetName = seq_len(nrow(bed.binned)), bed.binned, stringsAsFactors = TRUE)

  message("Generating GC% tracks ...")
  hts.gc <- loc.nt.gcc.hs.multi(loc.df = bed.binned, genome.pkg = genome.pkg, extend.multi = extend.multi, blocksize = blocksize, nthread = nthread)

  rm(bed.binned)
  
  hts.meta.key <- c("genome-species", "genome-version", "genome-package", "array_type", "track_type", "bin_size")
  hts.meta.value <- c(organism, genome, genome.pkg, "HTS", "GC", bin.size)
  # renorm.data <- list(tracks = hts.gc, info = list("genome-version" = genome, "genome-package" = genome.pkg, bin.size = bin.size, track.type = "GC"), bed.clean = bed.clean, bed.binned = bed.binned)
  renorm.data <- list(tracks = hts.gc, info = data.frame(key = hts.meta.key, value = hts.meta.value, stringsAsFactors = FALSE), bed.clean = bed.clean)
  rm(hts.gc, hts.meta.value, hts.meta.key, bed.clean)
  
  save("renorm.data", file = paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = paste0("_", genome, "_b", bin.size, ".GC.rda"), x = basename(bed.file), ignore.case = TRUE)), compress = "xz")

  message("Done.")
  if (return.data) return(renorm.data)
}

## Generates a BINpack (bins and pre-computed GC tracks) from a capture bed
## NEW : add distance to centromere track (and rank of it) as new tracks if a UCSC "cytoband/cytobandideo" table is given !
## NEW : add mappability track ! (from UCSC bigwig)
BINpack.Maker <- function(bed.file = NULL, bin.size = 50, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", gc.extend = c(0, 50, 100, 200, 400, 800, 1600, 3200, 6400), cytoband.file = NULL, mappability.bw = NULL, reptime.file = NULL, blacklist.file = NULL, blocksize = 1E+04, nthread = 1, out.dir = NULL, return.data = FALSE) {
  
  ## Checks
  if (is.null(bed.file)) stop("A BED file is required !", call. = FALSE)
  if (!file.exists(bed.file)) stop("Could not find the BED file !", call. = FALSE)
  if (is.null(out.dir)) {
    tmsg("NOTE : Checked / cleaned bed will be written in the same directory as source.")
    out.dir <- dirname(bed.file)
  } else {
    if (!file.exists(out.dir)) stop("Could not find the output directory !", call. = FALSE)
    if (!file.info(out.dir)$isdir) stop("out.dir is not a directory !", call. = FALSE)
  }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"), call. = FALSE)
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  ## GC (for GC% tracks)
  if (!is.null(gc.extend)) {
    if (any(!is.numeric(gc.extend))) stop("'gc.extend' should be a vector of positive ints !")
    if (any(gc.extend < 0)) stop("'gc.extend' should be a vector of positive ints !")
  }
  ## Cytoband (for dist-to-centromere tracks)
  if (!is.null(cytoband.file)) {
    if (!is.character(cytoband.file)) stop("'cytoband.file' (UCSC cytoband table file) is not a character string !")
    if (!file.exists(cytoband.file)) stop("Provided 'cytoband.file' (UCSC cytoband table file) does not exist !")
  }
  ## Mappability (from UCSC b36 bigwig)
  if (!is.null(mappability.bw)) {
    if (!is.character(mappability.bw)) stop("'mappability.bw' (UCSC mappability bigwig file) is not a character string !")
    if (!file.exists(mappability.bw)) stop("Provided 'mappability.bw' (UCSC mappability bigwig file) does not exist !")
  }
  ## Replication time tracks (from ASCAT Zenodo)
  if (!is.null(reptime.file)) {
    if (!is.character(reptime.file)) stop("'reptime.file' (ASCAT replication time file) is not a character string !")
    if (!file.exists(reptime.file)) stop("Provided 'reptime.file' (ASCAT replication time file) does not exist !")
  }
  ## Blacklist (from UCSC encode_blacklist_v2 table)
  if (!is.null(blacklist.file)) {
    if (!is.character(blacklist.file)) stop("'blacklist.file' (UCSC encode blacklist table file) is not a character string !")
    if (!file.exists(blacklist.file)) stop("Provided 'blacklist.file' (UCSC encode blacklist table file) does not exist !")
  }
  
  ## Force integers
  bin.size <- as.integer(bin.size)
  if(!is.null(gc.extend)) gc.extend <- round(gc.extend)
  
  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  organism <- BSgenome::organism(BSg.obj)
  
  ## Cleaning BED
  bed.clean <- BedCheck(bed.file = bed.file, genome.pkg = genome.pkg, out.dir = out.dir, return.data = TRUE)
  
  ## Remove blacklisted regions ?
  if (!is.null(blacklist.file)) {
    tmsg('Discarding blacklisted regions ...')
    bl_df <- data.table::fread(file = blacklist.file, sep = '\t', data.table = FALSE)
    colnames(bl_df)[1:3] <- c('chr', 'start', 'end')
    bl_df <- bl_df[bl_df$chr %in% bed.clean$chr,]
    # bl_gr <- GenomicRanges::makeGRangesFromDataFrame(df = bl_df)
    # bc_gr <- GenomicRanges::makeGRangesFromDataFrame(df = bed.clean)
    # sdf <- GenomicRanges::setdiff(x = bc_gr, y = bl_gr)
    bed.clean <- as.data.frame(GenomicRanges::setdiff(x = GenomicRanges::makeGRangesFromDataFrame(df = bed.clean), y = GenomicRanges::makeGRangesFromDataFrame(df = bl_df)))[,1:3]
    colnames(bed.clean) <- c('chr', 'start', 'end')
  }
  
  ## Binning
  tmsg(paste0("Performing binning (", bin.size, ") ..."))
  bed.binned <- bedBinner(bed = bed.clean, bin.size = bin.size, nthread = nthread)
  
  renorm.tracks <- list(pos = data.frame(ProbeSetName = seq_len(nrow(bed.binned)), bed.binned, stringsAsFactors = TRUE))
  rm(bed.binned)
  
  ## Add width
  tmsg("Adding regions width ...")
  renorm.tracks[['width']] <- data.frame(width = as.integer(renorm.tracks$pos$end - renorm.tracks$pos$start + 1))
  
  ## Add GC tracks ?
  if (!is.null(gc.extend)) {
    tmsg("Generating GC% tracks ...")
    renorm.tracks[['gc']] <- loc.nt.gcc.hs.multi(loc.df = renorm.tracks$pos, genome.pkg = genome.pkg, extend.multi = gc.extend, blocksize = blocksize, nthread = nthread)
  }
  
  ## Add distance to centroid (from UCSC cyrobandideo table) ?
  if (!is.null(cytoband.file)) {
    tmsg("... adding distance to centromere track.")
    
    ## Get cytobands
    cs <- read_ucsc_cytobandideo(file = cytoband.file, genomebuild = genome, species = organism)
    
    ## Compute centrotrack (and corresponding rank to avoid non-linear length effects)
    centro.tracks <- data.frame(mid = round((renorm.tracks$pos$start + renorm.tracks$pos$end)/2))
    centro.tracks$centrodist.rank <- centro.tracks$centrodist <- NA
    for (k in levels(renorm.tracks$pos$chr)) {
      kidx <- renorm.tracks$pos$chr == k
      centro.tracks$centrodist[kidx] <- abs(cs$chromosomes$centromere[cs$chromosomes$chrom == k] - centro.tracks$mid[kidx])
      centro.tracks$centrodist.rank[kidx] <- rank(centro.tracks$centrodist[kidx])
    }
    ## Add log10
    centro.tracks$centrodist.l10 <- log10(centro.tracks$centrodist)
    ## Reduce size
    centro.tracks$centrodist <- as.integer(centro.tracks$centrodist)
    centro.tracks$centrodist.rank <- as.integer(centro.tracks$centrodist.rank)
    centro.tracks$centrodist.l10 <- as.integer(round(centro.tracks$centrodist.l10 * 1E+06))
    
    ## Clean
    centro.tracks$mid <- NULL
    rm(cs, read_ucsc_cytobandideo)
    renorm.tracks[['centro']] <- centro.tracks
    rm(centro.tracks)
  }
  
  ## Add mappability (from UCSC bigwig b36) ?
  if (!is.null(mappability.bw)) {
    tmsg("... adding mappability track.")
    
    ## Import the mappaibility bw
    mapbw <- rtracklayer::import(con = mappability.bw, as = "GRanges")
    
    ## Convert to Granges
    rd_gr <- GenomicRanges::makeGRangesFromDataFrame(df = renorm.tracks$pos, keep.extra.columns = FALSE, ignore.strand = TRUE, start.field = 'start', end.field = 'end', seqnames.field = 'chr')
    
    ## Overlap
    fo <- tibble::as_tibble(GenomicRanges::findOverlaps(mapbw, rd_gr))
    fo$score <- mapbw$score[fo$queryHits]
    rm(rd_gr, mapbw)
    
    ## Merge
    fo <- dplyr::group_by(.data = fo, subjectHits)
    fom <- dplyr::summarize(.data = fo, score = median(score))
    rm(fo)
    
    ## Insert
    map.tracks <- data.frame(map = rep(NA, nrow(renorm.tracks$pos)))
    map.tracks$map[fom$subjectHits] <- fom$score
    rm(fom)
    
    ## Impute NAs ?
    # if (any(is.na(map.tracks$map))) map.tracks$map <- zoo::na.spline(map.tracks$map)
    if (any(is.na(map.tracks$map))) map.tracks$map <- stats::approxfun(seq_along(map.tracks$map), map.tracks$map, rule = 2)(seq_along(map.tracks$map))
    
    ## Reduce size
    map.tracks$map <- as.integer(round(map.tracks$map * 1000))
    
    renorm.tracks[['map']] <- map.tracks
    rm(map.tracks)
  }
  
  ## Add replication time tracks (from ASCAT) ?
  if (!is.null(reptime.file)) {
    tmsg("... adding time-replication tracks.")
    
    rtime <- data.table::fread(file = reptime.file, data.table = FALSE)
    colnames(rtime)[1] <- 'ID'
    
    rtime$Chrom <- paste0('chr', rtime$Chr)
    
    pos_gr <- GenomicRanges::makeGRangesFromDataFrame(df = renorm.tracks$pos[,c('chr', 'start', 'end')], keep.extra.columns = FALSE, ignore.strand = TRUE, start.field = 'start', end.field = 'end', seqnames.field = 'chr')
    rt_gr <- GenomicRanges::makeGRangesFromDataFrame(df = rtime[,c('Chrom', 'Position')], keep.extra.columns = FALSE, ignore.strand = TRUE, start.field = 'Position', end.field = 'Position', seqnames.field = 'Chrom')
  
    fo <- tibble::as_tibble(GenomicRanges::findOverlaps(rt_gr, pos_gr))
    rm(pos_gr, rt_gr)
    
    ## Synching each track
    tracknames <- colnames(rtime)[-c(1:3, ncol(rtime))]
    for (tn in tracknames) {
      message(tn)
      fotn <- fo
      fotn$tomerge <- rtime[[tn]][fotn$queryHits]
      
      ## Merge
      fotn <- dplyr::group_by(.data = fotn, subjectHits)
      fom <- dplyr::summarize(.data = fotn, tomerge = median(tomerge))
      
      ## Insert
      renorm.tracks[['reptime']][[tn]] <- rep(NA, nrow(renorm.tracks$pos))
      renorm.tracks[['reptime']][[tn]][fom$subjectHits] <- fom$tomerge
      track_range <- range(renorm.tracks[['reptime']][[tn]], na.rm = TRUE)
      
      ## Impute with splines
      renorm.tracks[['reptime']][[tn]] <- zoo::na.spline(renorm.tracks[['reptime']][[tn]])
      ## Prune outliers
      renorm.tracks[['reptime']][[tn]][renorm.tracks[['reptime']][[tn]] < min(track_range)] <- min(track_range)
      renorm.tracks[['reptime']][[tn]][renorm.tracks[['reptime']][[tn]] > max(track_range)] <- max(track_range)
      rm(fotn, fom)
    }
    renorm.tracks$reptime <- as.data.frame(renorm.tracks$reptime)
    ## Convert to int
    for (rt in colnames(renorm.tracks$reptime)) renorm.tracks$reptime[[rt]] <- as.integer(round(renorm.tracks$reptime[[rt]] * 1000))
    ## Clean
    rm(rtime)
  }
  
  tmsg("Aggregating results.")
  
  hts.meta.key <- c("genome-species", "genome-version", "genome-package", "array_type", "bin_size")
  hts.meta.value <- c(organism, genome, genome.pkg, "HTS", bin.size)
  renorm.data <- list(tracks = renorm.tracks, info = data.frame(key = hts.meta.key, value = hts.meta.value, stringsAsFactors = FALSE), bed.clean = bed.clean)
  rm(renorm.tracks, hts.meta.value, hts.meta.key, bed.clean)
  
  tmsg("Saving object.")
  
  save("renorm.data", file = paste0(out.dir, "/", sub(pattern = "\\.bed$", replacement = paste0("_", genome, "_b", bin.size, ".tracks.rda"), x = basename(bed.file), ignore.case = TRUE)), compress = "xz")
  
  tmsg("Done.")
  if (return.data) return(renorm.data)
}

## Bin WES data from a BAM and BINpack
WES.Bin_bak <- function(BAM = NULL, BINpack = NULL, samplename = "SAMPLE", Q = 20, blocksize = 1000, nsubthread = 1, cluster.type = "PSOCK", out.dir = getwd(), return.data  = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  ## CHECKS (files/parameters)
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"), call. = FALSE)
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"), call. = FALSE)
  if (is.null(BAM)) stop(tmsg("A BAM file is required !"), call. = FALSE)
  if (!file.exists(BAM)) stop(tmsg("Could not find the BAM file !"), call. = FALSE)
  if (!is.numeric(Q)) stop(tmsg("Q must be numeric !"), call. = FALSE)
  if (is.null(out.dir)) stop(tmsg("An output directory is required !"), call. = FALSE)
  if (!file.exists(out.dir)) stop(tmsg("Could not find the output directory !"), call. = FALSE)
  if (!file.info(out.dir)$isdir) stop(tmsg("out.dir is not a directory"), call. = FALSE)
  if (!return.data & !write.data) stop(tmsg("Data should be returned and/or written on disk, but not none !"), call. = FALSE)
  if (Q < 0) stop(tmsg("Q should be positive !"), call. = FALSE)
  
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
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !')), call. = FALSE) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(library(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  ## Files controls
  tmsg("Checking BINpack and BAMs compatibility ...")
  
  ## Inspecting BAM header
  BAM.h <- Rsamtools::scanBamHeader(BAM)
  bed.data <- renorm.data$tracks$pos
  renorm.data$tracks <- NULL
  gc()
  
  # if (!all(names(refBAM.h[[1]]$targets) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and Test BAM are not compatible (different chr names) !"), call. = FALSE)
  if (!all(unique(bed.data$chr) %in% names(BAM.h[[1]]$targets))) stop(tmsg("BAM and BED are not compatible (different chr names) !"), call. = FALSE)
  if (!all(unique(bed.data$chr) %in% BSgenome::seqnames(BSg.obj))) stop(tmsg("BED and BSgenome are not compatible (different chr names) !"), call. = FALSE)
  
  ## Identifying the platform
  BAM.h.unl <- unlist(BAM.h)
  manuf.hentry <- grep(pattern = "^PL:", BAM.h.unl)[1]
  manufacturer <- if(!is.na(manuf.hentry)) sub(pattern = "^PL:", replacement = "", BAM.h.unl[[manuf.hentry]]) else "NA"
  
  meta.b <- list(
    samplename = samplename,
    source = "WES",
    source.file = list(BAM = BAM, BINpack = BINpack),
    type = "WES",
    manufacturer = manufacturer,
    species = GenomeInfoDb::organism(BSg.obj),
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = "NA"
  )
  
  meta.w <- list(
    BAM.header = paste0(BAM.h, collapse = " "),
    samtools.Q = Q,
    # bin.size = GC.data$info$bin.size
    bin.size = as.numeric(renorm.data$info$value[renorm.data$info$key == "bin_size"])
  )
  rm(BAM.h)
  
  ### Common BAM flags
  param.FLAG <- Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
  param.PILEUP <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = TRUE, include_deletions = FALSE, include_insertions = FALSE)
  
  COUNTS.all <- pileup.go(BAM = BAM, genome.pkg = genome.pkg, bed.data = bed.data, scanBamFlag = param.FLAG, pileupParam = param.PILEUP, nsubthread = nsubthread, cluster.type = cluster.type)
  ## Merge CNs
  tmsg('Merging CN data ...')
  CN.all <- do.call(what = rbind, args = lapply(X = COUNTS.all, function(x) { x$CN } ))
  
  ## Merge SNPs
  tmsg('Merging SNP data ...')
  SNP.all <- do.call(what = rbind, args = lapply(X = COUNTS.all, function(x) { x$SNP } ))
  
  rm(COUNTS.all)
  gc()
  
  ## Building HTSobj
  tmsg('Building the output object ...')
  CN.all$bin <- as.integer(CN.all$bin)
  SNP.all$bin <- as.integer(SNP.all$bin)
  
  ### summaries (recoded function as R internal summary uses too much RAM !)
  my.summary <- function(myv = NULL) {
    vsum <- c(min(myv, na.rm = TRUE), quantile(myv, .25, na.rm = TRUE), median(myv, na.rm = TRUE), mean(myv, na.rm = TRUE), quantile(myv, .75, na.rm = TRUE), max(myv, na.rm = TRUE))
    names(vsum) <- c("min", "q25", "median", "mean", "q75", "max")
    return(vsum)
  }
  meta.w$BIN.tot.count.mean.summary <- my.summary(CN.all$tot_count[!is.na(CN.all$tot_count)])
  meta.w$SNP.tot.count.summary <- my.summary(SNP.all$tot_count[!is.na(SNP.all$tot_count)])
  gc()
  
  ## Cleaning uncovered chr levels
  CN.all$chr <- droplevels(CN.all$chr)
  SNP.all$chr <- droplevels(SNP.all$chr)
  
  HTSobj <- list(RD = CN.all, SNP = SNP.all, meta = list(basic = meta.b, WES = meta.w))
  rm(CN.all, SNP.all)
  gc()
  
  ## QC : Computing coverages
  tmsg("Computing coverages ...")
  gw.rd <- sum(HTSobj$RD$end - HTSobj$RD$start +1)
  gw.snp <- nrow(HTSobj$SNP)
  cov_cuts <- c(1, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200)
  rd.cov <- data.frame(cuts = cov_cuts, do.call(rbind, lapply(cov_cuts, function(x) {
    rd.in <- HTSobj$RD$tot_count >= x
    snprd.in <- HTSobj$SNP$tot_count >= x
    cut.cov <- if(!any(rd.in)) NA else (sum(HTSobj$RD$end[rd.in] - HTSobj$RD$start[rd.in] +1)/gw.rd)
    snpcut.cov <- if(!any(snprd.in)) NA else (length(which(snprd.in))/gw.snp)
    return(c(cut.cov, snpcut.cov))
  })))
  dimnames(rd.cov) <- list(cov_cuts, c("MinDepth", "BINCoverage", "BAFCoverage"))
  
  
  if (write.data || plot) dir.create(paste0(out.dir, "/", samplename))
  if (write.data) write.table(rd.cov, file = paste0(out.dir, "/", samplename,  "/", samplename, '_WES_', genome, "_b", meta.w$bin.size, "_coverage.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  ## QC : Plotting coverages
  if (plot) {
    tmsg('Plotting coverage ...')
    ### Coverage plot
    png(paste0(out.dir, "/", samplename, "/", samplename, "_WES_", genome, "_b", meta.w$bin.size, "_coverage.png"), 1200, 900)
    plot(rd.cov$MinDepth, rd.cov$BAFCoverage, type = "b", col = 2, lty = 3, pch = 20, main = paste0(HTSobj$meta$basic$samplename, "\nCoverage Plot"), xlab = "Minimum depth", ylab = "Coverage", ylim = c(0,1), xaxp = c(0,200,10))
    abline(v = rd.cov$MinDepth, lty = 2, col = "grey75")
    abline(h = seq(0,1,.1), lty = 2, col = "grey75")
    lines(rd.cov$MinDepth, rd.cov$BINCoverage, type = "b", col = 2)
    legend("topright", legend = c("SNP", "BIN"), inset = .02, col = c(2,2), lty = c(3,1), pch = c(20,1))
    dev.off()
    
    ### RD plots
    #### Create genomic positions
    chr_max <- vapply(X = levels(HTSobj$RD$chr), function(x) { max(HTSobj$RD$end[HTSobj$RD$chr == x])}, 1.0)
    chr2add <- c(0, cumsum(chr_max)[-length(chr_max)])
    names(chr2add) <- levels(HTSobj$RD$chr)
    tmsg('Plotting raw depth profile ...')
    png(paste0(out.dir, "/", HTSobj$meta$basic$samplename, "/", HTSobj$meta$basic$samplename, "_WES_", HTSobj$meta$basic$genome, "_b", HTSobj$meta$WES$bin.size, "_rawdepth.png"), width = 1600, height = 1050)
    par(mfrow = c(3,1))
    ## BIN
    l10 <- log10(HTSobj$RD$tot_count +1)
    l10.med <- median(l10, na.rm = TRUE)
    x.range <- c(0, max(chr2add)+chr_max[length(chr_max)])
    rm(chr_max)
    y.range <- c(l10.med -1.5, l10.med +1.5)
    plot(HTSobj$RD$start + chr2add[HTSobj$RD$chr], l10,  pch = ".", cex = 2, col = "grey66", xaxs = "i", yaxs = "i", xlab = "Position", ylab = "log10(BIN.RD+1)", main = paste0(HTSobj$meta$basic$samplename, " BIN depth (raw)"), xlim = x.range, ylim = y.range)
    abline(h = l10.med, lty = 2, col = "black")
    abline(v = chr2add, lty = 2, col = 4)
    lines(HTSobj$RD$start + chr2add[HTSobj$RD$chr], suppressWarnings(runmed(l10, 9999)), col = 2, lwd = 3)
    ## SNP
    l10 <- log10(HTSobj$SNP$tot_count +1)
    l10.med <- median(l10, na.rm = TRUE)
    plot.range <- c(l10.med -2, l10.med +2)
    plot(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], l10,  pch = ".", cex = 2, col = "grey66", xaxs = "i", yaxs = "i", xlab = "Position", ylab = "log10(SNP.RD+1)", main = paste0(HTSobj$meta$basic$samplename, " SNP depth (raw))"), xlim = x.range, ylim = y.range)
    abline(h = l10.med, lty = 2, col = "black")
    lines(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], suppressWarnings(runmed(l10, 9999)), col = 2, lwd = 3)
    abline(v = chr2add, lty = 2, col = 4)
    ## BAF
    BAF <- HTSobj$SNP$alt_count / HTSobj$SNP$tot_count
    mBAF <- BAF < .5
    BAF[mBAF] <- 1 - BAF[mBAF]
    plot(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], BAF,  pch = ".", cex = 2, col = "grey66", xaxs = "i", yaxs = "i", xlab = "Position", ylab = "BAF (raw)", main = paste0(HTSobj$meta$basic$samplename, " SNP BAF (raw))"), xlim = x.range, ylim = c(0,1))
    points(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], 1 - BAF,  pch = ".", cex = 2, col = "grey66")
    abline(h = c(0,.5,1), lty = 3, col = "black")
    abline(v = chr2add, lty = 2, col = 4)
    dev.off()
    rm(BAF, l10)
  }
  rm(rd.cov)
  
  ## Saving
  if (write.data) {
    tmsg('Writing the output object ...')
    tmsg("Saving counts data ...")
    saveRDS(HTSobj, file = paste0(out.dir, "/", samplename, "/", samplename, "_", genome, "_b", meta.w$bin.size, "_binned.RDS"), compress = "xz")
  }
  if (return.data) return(HTSobj)
}

## Bin WES data from a BAM and BINpack
HTS.Bin <- function(BAM = NULL, BINpack = NULL, samplename = "SAMPLE", Q = 20, blocksize = 1000, nsubthread = 1, cluster.type = "PSOCK", out.dir = getwd(), return.data  = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  ## CHECKS (files/parameters)
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"), call. = FALSE)
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"), call. = FALSE)
  if (is.null(BAM)) stop(tmsg("A BAM file is required !"), call. = FALSE)
  if (!file.exists(BAM)) stop(tmsg("Could not find the BAM file !"), call. = FALSE)
  if (!is.numeric(Q)) stop(tmsg("Q must be numeric !"), call. = FALSE)
  if (is.null(out.dir)) stop(tmsg("An output directory is required !"), call. = FALSE)
  if (!file.exists(out.dir)) stop(tmsg("Could not find the output directory !"), call. = FALSE)
  if (!file.info(out.dir)$isdir) stop(tmsg("out.dir is not a directory"), call. = FALSE)
  if (!return.data & !write.data) stop(tmsg("Data should be returned and/or written on disk, but not none !"), call. = FALSE)
  if (Q < 0) stop(tmsg("Q should be positive !"), call. = FALSE)
  
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
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !')), call. = FALSE) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(library(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  ## Files controls
  tmsg("Checking BINpack and BAMs compatibility ...")
  
  ## Inspecting BAM header
  BAM.h <- Rsamtools::scanBamHeader(BAM)
  bed.data <- renorm.data$tracks$pos
  renorm.data$tracks <- NULL
  gc()
  
  # if (!all(names(refBAM.h[[1]]$targets) %in% names(testBAM.h[[1]]$targets))) stop(tmsg("Reference BAM and Test BAM are not compatible (different chr names) !"), call. = FALSE)
  if (!all(unique(bed.data$chr) %in% names(BAM.h[[1]]$targets))) stop(tmsg("BAM and BED are not compatible (different chr names) !"), call. = FALSE)
  if (!all(unique(bed.data$chr) %in% BSgenome::seqnames(BSg.obj))) stop(tmsg("BED and BSgenome are not compatible (different chr names) !"), call. = FALSE)
  
  ## Identifying the platform
  BAM.h.unl <- unlist(BAM.h)
  manuf.hentry <- grep(pattern = "^PL:", BAM.h.unl)[1]
  manufacturer <- if(!is.na(manuf.hentry)) sub(pattern = "^PL:", replacement = "", BAM.h.unl[[manuf.hentry]]) else "NA"
  
  meta.b <- list(
    samplename = samplename,
    source = "HTS",
    source.file = list(BAM = BAM, BINpack = BINpack),
    type = "HTS",
    manufacturer = manufacturer,
    species = GenomeInfoDb::organism(BSg.obj),
    genome = genome,
    genome.pkg = genome.pkg,
    predicted.gender = "NA"
  )
  
  meta.w <- list(
    BAM.header = paste0(BAM.h, collapse = " "),
    samtools.Q = Q,
    # bin.size = GC.data$info$bin.size
    bin.size = as.numeric(renorm.data$info$value[renorm.data$info$key == "bin_size"])
  )
  rm(BAM.h)
  
  ### Common BAM flags
  param.FLAG <- Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isNotPassingQualityControls = FALSE, isDuplicate = FALSE)
  param.PILEUP <- Rsamtools::PileupParam(distinguish_strands = FALSE, max_depth = 5E+04, min_base_quality = Q, min_nucleotide_depth = 0, distinguish_nucleotides = TRUE, include_deletions = FALSE, include_insertions = FALSE)
  
  COUNTS.all <- pileup.go(BAM = BAM, genome.pkg = genome.pkg, bed.data = bed.data, scanBamFlag = param.FLAG, pileupParam = param.PILEUP, nsubthread = nsubthread, cluster.type = cluster.type)
  
  ## Retrieve order
  bin_order <- do.call(what = c, args = lapply(X = seq_along(COUNTS.all), function(x) { COUNTS.all[[x]]$id } ))
  
  ## Merge CNs
  tmsg('Merging CN data ...')
  # CN.all <- do.call(what = rbind, args = lapply(X = COUNTS.all, function(x) { x$CN } ))
  CN.all <- do.call(what = rbind, args = lapply(X = bin_order, function(x) { COUNTS.all[[x]]$CN } ))
  
  ## Merge SNPs
  tmsg('Merging SNP data ...')
  # SNP.all <- do.call(what = rbind, args = lapply(X = COUNTS.all, function(x) { x$SNP } ))
  SNP.all <- do.call(what = rbind, args = lapply(X = bin_order, function(x) { COUNTS.all[[x]]$SNP } ))
  
  rm(COUNTS.all)
  gc()
  
  ## Building HTS obj
  tmsg('Building the output object ...')
  CN.all$bin <- as.integer(CN.all$bin)
  SNP.all$bin <- as.integer(SNP.all$bin)
  
  meta.w$BIN.tot.count.mean.summary <- my.summary(CN.all$tot_count[!is.na(CN.all$tot_count)])
  meta.w$SNP.tot.count.summary <- my.summary(SNP.all$tot_count[!is.na(SNP.all$tot_count)])
  gc()
  
  ## Cleaning uncovered chr levels
  CN.all$chr <- droplevels(CN.all$chr)
  SNP.all$chr <- droplevels(SNP.all$chr)
  
  HTSobj <- list(RD = CN.all, SNP = SNP.all, meta = list(basic = meta.b, HTS = meta.w))
  rm(CN.all, SNP.all)
  gc()
  
  ## QC : Computing coverages
  tmsg("Computing coverages ...")
  gw.rd <- sum(HTSobj$RD$end - HTSobj$RD$start +1)
  gw.snp <- nrow(HTSobj$SNP)
  cov_cuts <- c(1, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200)
  rd.cov <- data.frame(cuts = cov_cuts, do.call(rbind, lapply(cov_cuts, function(x) {
    rd.in <- HTSobj$RD$tot_count >= x
    snprd.in <- HTSobj$SNP$tot_count >= x
    cut.cov <- if(!any(rd.in)) NA else (sum(HTSobj$RD$end[rd.in] - HTSobj$RD$start[rd.in] +1)/gw.rd)
    snpcut.cov <- if(!any(snprd.in)) NA else (length(which(snprd.in))/gw.snp)
    return(c(cut.cov, snpcut.cov))
  })))
  dimnames(rd.cov) <- list(cov_cuts, c("MinDepth", "BINCoverage", "BAFCoverage"))
  
  
  if (write.data || plot) dir.create(paste0(out.dir, "/", samplename))
  if (write.data) write.table(rd.cov, file = paste0(out.dir, "/", samplename,  "/", samplename, '_HTS_', genome, "_b", meta.w$bin.size, "_coverage.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  ## QC : Plotting coverages
  if (plot) {
    tmsg('Plotting coverage ...')
    ### Coverage plot
    png(paste0(out.dir, "/", samplename, "/", samplename, "_HTS_", genome, "_b", meta.w$bin.size, "_coverage.png"), 1200, 900)
    plot(rd.cov$MinDepth, rd.cov$BAFCoverage, type = "b", col = 2, lty = 3, pch = 20, main = paste0(HTSobj$meta$basic$samplename, "\nCoverage Plot"), xlab = "Minimum depth", ylab = "Coverage", ylim = c(0,1), xaxp = c(0,200,10))
    abline(v = rd.cov$MinDepth, lty = 2, col = "grey75")
    abline(h = seq(0,1,.1), lty = 2, col = "grey75")
    lines(rd.cov$MinDepth, rd.cov$BINCoverage, type = "b", col = 2)
    legend("topright", legend = c("SNP", "BIN"), inset = .02, col = c(2,2), lty = c(3,1), pch = c(20,1))
    dev.off()
    
    ### RD plots
    #### Create genomic positions
    chr_max <- vapply(X = levels(HTSobj$RD$chr), function(x) { max(HTSobj$RD$end[HTSobj$RD$chr == x])}, 1.0)
    chr2add <- c(0, cumsum(chr_max)[-length(chr_max)])
    names(chr2add) <- levels(HTSobj$RD$chr)
    tmsg('Plotting raw depth profile ...')
    png(paste0(out.dir, "/", HTSobj$meta$basic$samplename, "/", HTSobj$meta$basic$samplename, "_HTS_", HTSobj$meta$basic$genome, "_b", HTSobj$meta$HTS$bin.size, "_rawdepth.png"), width = 1600, height = 1050)
    par(mfrow = c(3,1))
    ## BIN
    l10 <- log10(HTSobj$RD$tot_count +1)
    l10.med <- median(l10, na.rm = TRUE)
    x.range <- c(0, max(chr2add)+chr_max[length(chr_max)])
    rm(chr_max)
    y.range <- c(l10.med -1.5, l10.med +1.5)
    plot(HTSobj$RD$start + chr2add[HTSobj$RD$chr], l10,  pch = ".", cex = 2, col = "grey66", xaxs = "i", yaxs = "i", xlab = "Position", ylab = "log10(BIN.RD+1)", main = paste0(HTSobj$meta$basic$samplename, " BIN depth (raw)"), xlim = x.range, ylim = y.range)
    abline(h = l10.med, lty = 2, col = "black")
    abline(v = chr2add, lty = 2, col = 4)
    lines(HTSobj$RD$start + chr2add[HTSobj$RD$chr], suppressWarnings(runmed(l10, 9999)), col = 2, lwd = 3)
    ## SNP
    l10 <- log10(HTSobj$SNP$tot_count +1)
    l10.med <- median(l10, na.rm = TRUE)
    plot.range <- c(l10.med -2, l10.med +2)
    plot(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], l10,  pch = ".", cex = 2, col = "grey66", xaxs = "i", yaxs = "i", xlab = "Position", ylab = "log10(SNP.RD+1)", main = paste0(HTSobj$meta$basic$samplename, " SNP depth (raw))"), xlim = x.range, ylim = y.range)
    abline(h = l10.med, lty = 2, col = "black")
    lines(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], suppressWarnings(runmed(l10, 9999)), col = 2, lwd = 3)
    abline(v = chr2add, lty = 2, col = 4)
    ## BAF
    BAF <- HTSobj$SNP$alt_count / HTSobj$SNP$tot_count
    mBAF <- BAF < .5
    BAF[mBAF] <- 1 - BAF[mBAF]
    plot(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], BAF,  pch = ".", cex = 2, col = "grey66", xaxs = "i", yaxs = "i", xlab = "Position", ylab = "BAF (raw)", main = paste0(HTSobj$meta$basic$samplename, " SNP BAF (raw)"), xlim = x.range, ylim = c(0,1))
    points(HTSobj$SNP$pos + chr2add[HTSobj$SNP$chr], 1 - BAF,  pch = ".", cex = 2, col = "grey66")
    abline(h = c(0,.5,1), lty = 3, col = "black")
    abline(v = chr2add, lty = 2, col = 4)
    dev.off()
    rm(BAF, l10)
  }
  rm(rd.cov)
  
  ## Saving
  if (write.data) {
    tmsg('Writing the output object ...')
    tmsg("Saving counts data ...")
    saveRDS(HTSobj, file = paste0(out.dir, "/", samplename, "/", samplename, "_", genome, "_b", meta.w$bin.size, "_binned.RDS"), compress = "xz")
  }
  if (return.data) return(HTSobj)
}

## Test
# system.time(WES.Bin_new(
#   BAM = "~/WORKSPACE/DEVELOPMENT/EACON/WES_data/1371AJM0006_hg19.bam"
#   ,
#   BINpack = "~/WORKSPACE/DEVELOPMENT/EACON/WES_data/Sureselect_v5.padded.order.collapse_hs37d5_b50.GC.rda"
#   ,
#   samplename = "1371AJM0006"
#   ,
#   Q = 20
#   ,
#   out.dir = dirname(path = BAM)
#   ,
#   nsubthread = 4
#   ,
#   cluster.type = "PSOCK"
#   ,
#   return.data = FALSE
#   ,
#   write.data = TRUE
#   ,
#   plot = TRUE
# ))

## Performs the binning of BAMs using a BINpack, batch mode
HTS.Bin.Batch <- function(BAM.list.file = NULL, BINpack = NULL, nthread = 1, cluster.type = "PSOCK", ...) {

  if (!file.exists(BAM.list.file)) stop("Could not find BAM.list.file !", call. = FALSE)
  message("Reading and checking BAM.list.file ...")
  myBAMs <- read.table(file = BAM.list.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("BAM", "SampleName")
  head.chk <- all(colnames(BAM.list.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in BAM.list.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(myBAMs)))
    stop("Invalid header.", call. = FALSE)
  }

  fbam.chk <- file.exists(myBAMs$BAM)
  
  if (!all(fbam.chk)) {
    message("Some BAM file(s) from the BAM.list.file could not be found (wrong path or filename ?) !")
    message("Missing BAM file(s) :")
    message(myBAMs$BAM[which(!fbam.chk)])
    stop("Missing BAM file(s).", call. = FALSE)
  }
  sn.chk <- duplicated(myBAMs$SampleName)
  if (any(sn.chk)) {
    message("BAM.list.file contains duplicated samplename(s) !")
    message(myBAMs$SampleName[which(sn.chk)])
    stop("Duplicated samplename(s).", call. = FALSE)
  }
  bamn.chk <- duplicated(myBAMs$BAM)
  if (any(bamn.chk)) {
    message("BAM.list.file contains duplicated BAM file(s) !")
    message(myBAMs$BAM[which(bamn.chk)])
    stop("Duplicated BAM(s).", call. = FALSE)
  }

  ## Adjusting cores/threads
  message("Adjusting number of cores if needed ...")
  if (is.null(nthread)) nthread <- parallel::detectCores(logical = TRUE) -1
  if (nrow(myBAMs) < nthread) nthread <- nrow(myBAMs)

  message("Running HTS.Bin() in batch mode ...")
  message(paste0("Found ", nrow(myBAMs), " samples to process ..."))
  current.bitmapType <- getOption("bitmapType")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  eacon.batchres <- foreach::foreach(r = seq_len(nrow(myBAMs)), .inorder = TRUE, .errorhandling = "stop", .export = c("EaCoN.set.bitmapType", "HTS.Bin", "tmsg")) %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    HTS.Bin(BAM = myBAMs$BAM[r], BINpack = BINpack, samplename = myBAMs$SampleName[r], cluster.type = cluster.type, ...)
  }
  parallel::stopCluster(cl)
}

## Pairs two binned HTS data into a two-tracks profile (by ex, T+N)
HTS.Pair <- function(test_binned = NULL, ref_binned = NULL, force_samplename = NULL, out.dir = getwd(), return.data = FALSE, write.data = TRUE) {
  ## Add missing bins
  
  ## RD
  tmsg('Pairing RD ...')
  test_binned$RD <- dplyr::full_join(x = test_binned$RD, y = ref_binned$RD, by = c('chr', 'start', 'end', 'bin'))
  colnames(test_binned$RD) <- c('chr', 'start', 'end', 'bin', 'tot_count.test', 'tot_count.ref')
  for (cn in c('tot_count.test', 'tot_count.ref')) {
    if (any(is.na(test_binned$RD[[cn]]))) test_binned$RD[[cn]][is.na(test_binned$RD[[cn]])] <- 0
  }
  ref_binned$RD <- NULL
  ## SNP
  tmsg('Pairing SNP ...')
  test_binned$SNP <- dplyr::full_join(x = test_binned$SNP, y = ref_binned$SNP, by = c('chr', 'pos', 'bin'))
  test_binned$SNP <- test_binned$SNP[order(test_binned$SNP$bin, test_binned$SNP$pos),]
  colnames(test_binned$SNP) <- c('chr', 'pos', 'bin', 'tot_count.test', 'alt_count.test', 'tot_count.ref', 'alt_count.ref')
  ref_binned$SNP <- NULL
  ## META
  tmsg('Pairing metadata ...')
  ### BASIC
  test_binned$meta$basic$samplename <- if (!is.null(force_samplename)) force_samplename else paste0(test_binned$meta$basic$samplename, '_vs_', ref_binned$meta$basic$samplename)
  names(test_binned$meta$basic$source.file)[names(test_binned$meta$basic$source.file) == 'BAM'] <- 'testBAM'
  test_binned$meta$basic$source.file$refBAM <- ref_binned$meta$basic$source.file[['BAM']]
  names(test_binned$meta$HTS)[names(test_binned$meta$HTS) == 'BAM.header'] <- 'testBAM.header'
  ref_binned$meta$basic <- NULL
  ### HTS
  test_binned$meta$HTS$refBAM.header <- ref_binned$meta$HTS$BAM.header
  ref_binned$meta$HTS <- NULL
  test_binned$meta$HTS$BIN.tot.count.test.mean.summary <- my.summary(test_binned$RD$tot_count.test[!is.na(test_binned$RD$tot_count.test)])
  test_binned$meta$HTS$BIN.tot.count.ref.mean.summary <- my.summary(test_binned$RD$tot_count.ref[!is.na(test_binned$RD$tot_count.ref)])
  test_binned$meta$HTS$SNP.tot.count.test.summary <- my.summary(test_binned$SNP$tot_count.test[!is.na(test_binned$SNP$tot_count.test)])
  test_binned$meta$HTS$SNP.tot.count.ref.summary <- my.summary(test_binned$SNP$tot_count.ref[!is.na(test_binned$SNP$tot_count.ref)])
  test_binned$meta$HTS$BIN.tot.count.mean.summary <- test_binned$meta$HTS$SNP.tot.count.summary <- NULL
  rm(ref_binned)
  if(write.data) {
    tmsg('Saving paired object ...')
    o2 <- paste(c(out.dir, test_binned$meta$basic$samplename), collapse = '/')
    dir.create(path = o2, recursive = TRUE)
    saveRDS(object = test_binned, file = paste0(o2, '/', paste(c(test_binned$meta$basic$samplename, test_binned$meta$basic$source, test_binned$meta$basic$genome, paste0('b', test_binned$meta$HTS$bin.size), 'paired.RDS'), collapse = '_')), compress = 'xz')
  }
  if(return.data) return(test_binned)
}

## Wrapper to HTS.Pair for files
HTS.Pair.ff <- function(test_binned.RDS = NULL, ref_binned.RDS = NULL, ...) {
  if (any(!file.exists(test_binned.RDS, ref_binned.RDS))) stop("'test_binned.RDS' and/or 'ref_binned.RDS' do(es) not exist !")
  HTS.Pair(test_binned = readRDS(test_binned.RDS), ref_binned = readRDS(ref_binned.RDS), ...)
}

### WIP
# target <- readRDS('/home/job/WORKSPACE/ElsaLAB/CNA/EACON_TESTS/RESULTS/P-FE8113-TMN-150X_vs_P-FE8113-REF-100X/P-FE8113-TMN-150X_vs_P-FE8113-REF-100X_hg38_b50_binned.RDS')
# ref_binned.RDS <- '/home/job/WORKSPACE/ElsaLAB/CNA/EACON_TESTS/DATA/TEST_EaCoN_v0.4/P-FE8113-TMN-150X/P-FE8113-TMN-150X_hg38_b50_binned.RDS'
# test_binned.RDS <- '/home/job/WORKSPACE/ElsaLAB/CNA/EACON_TESTS/DATA/TEST_EaCoN_v0.4/P-FE8113-REF-100X/P-FE8113-REF-100X_hg38_b50_binned.RDS'


## Performs the normalization of HTS L2R and BAF signals
HTS.Normalize_bak <- function(data = NULL, BINpack = NULL, use.tracks = c('gc'), method = 'loess', RD.tot.min = 20, RD.alt.min = 3, BAF.hetmin = .33, sex.chr = c("chrX", "chrY"), TumorBoost = FALSE, out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE) {

  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta2/WES/TCGA-A7-A0CE-01A_vs_10A")
  # data <- readRDS("Sample_PHEO_AG_HS_048_DNA_hg38_b50_binned.RDS")
  # BINpack <- "V4-UTRs.hg38.fragment_targets_minimal_sorted_longChr_hg38_b50.GC.rda"
  # use.tracks <- "gc"
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
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"), call. = FALSE)
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"), call. = FALSE)
  if (is.null(use.tracks)) stop('No track type to use (from BINPack) provided : normalization is not possible !')
  if (RD.tot.min < 0) stop(tmsg("RD.tot.min must be >= 0 !"), call. = FALSE)
  if (RD.alt.min <= 0) stop(tmsg("RD.alt.min must be > 0 !"), call. = FALSE)

  ## TAGS
  data$meta$HTS$TumorBoost <- as.character(TumorBoost)
  data$meta$HTS$RD.tot.min <- RD.tot.min
  data$meta$HTS$renorm.use.tracks <- use.tracks
  samplename <- data$meta$basic$samplename

  ## Loading BINpack
  load(BINpack)
  gc()
  
  ## Checks from BINpack
  avail.tracks <- names(renorm.data$tracks)
  avail.tracks <- avail.tracks[!avail.tracks == 'pos']
  
  if (!all(use.tracks %in% avail.tracks)) stop(paste0('Available normalization track type(s) is/are : ["', paste(avail.tracks, collapse = '", "'), '"], while requested track type(s) is/are : ["', paste(use.tracks, collapse = '", "'), '"].'))

  genome.pkg <- renorm.data$info$value[renorm.data$info$key == "genome-package"]
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }

  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  # requireNamespace(genome.pkg, quietly = TRUE)
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  cs <- chromobjector(BSg.obj)

  ## BAF HANDLING
  tmsg("Processing BAF ...")

  ### BAF : Filtering (low depth)
  rd.ori <- nrow(data$SNP)
  
  #### 1) Filtering for low ref or test depth
  RDlow <- data$SNP$tot_count.test < RD.tot.min | data$SNP$tot_count.ref < RD.tot.min
  if (length(which(RDlow)) == nrow(data$SNP)) stop(tmsg("All SNP positions were discarded for their low read count ! You may consider lowering the BAF.tot.min value."), call. = FALSE)
  tmsg(paste0("Removed ", length(which(RDlow)), " (", round(length(which(RDlow)) / rd.ori * 100, digits = 2), "%) SNP positions with low depth (<", RD.tot.min, ")"))
  data$SNP <- data$SNP[!RDlow,]
  
  #### 2) Filtering for low alt test depth
  BRDlow <- data$SNP$alt_count.test < RD.alt.min
  if (length(which(BRDlow)) == nrow(data$SNP)) stop(tmsg("All SNP positions were discarded for their low alternative allele count ! You may consider lowering the RD.alt.min value."), call. = FALSE)
  tmsg(paste0("Removed ", length(which(BRDlow)), " (", round(length(which(BRDlow)) / rd.ori * 100, digits = 2), "%) SNP positions with low alt RD (<", RD.alt.min, ")"))
  data$SNP <- data$SNP[!BRDlow,]
  gc()

  ### Computing BAF
  data$SNP$BAF.test.ori <- data$SNP$alt_count.test / data$SNP$tot_count.test
  data$SNP$BAF.ref.ori <- data$SNP$alt_count.ref / data$SNP$tot_count.ref
  odd.idx <- which(data$SNP$pos %% 2 == 1)
  data$SNP$BAF.test.ori[odd.idx] <- -data$SNP$BAF.test.ori[odd.idx] +1L
  data$SNP$BAF.ref.ori[odd.idx] <- -data$SNP$BAF.ref.ori[odd.idx] +1L

  ### Test
  data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = data$SNP$BAF.test.ori, newname = "BAF.test.ori")
  ### Ref
  data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = data$SNP$BAF.ref.ori, newname = "BAF.ref.ori", type = "ref")
  
  ###### Computing LOR (and variance)
  rcmat <- round(cbind(data$SNP$BAF.test.ori*data$SNP$tot_count.test, (1-data$SNP$BAF.test.ori)*data$SNP$tot_count.test))
  data$SNP$LOR <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
  data$SNP$LORvar <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
  rm(rcmat)
  gc()
  
  ### BAF : TumorBoost
  if (TumorBoost) {
    tmsg('Applying TumorBoost BAF normalization ...')
    # data$SNP$BAF.test <- data$SNP$BAF.test.TB <- as.numeric(aroma.light::normalizeTumorBoost(data$SNP$BAF.test, data$SNP$BAF.ref, flavor = "v4", preserveScale = FALSE))
    BTB <- as.numeric(aroma.light::normalizeTumorBoost(data$SNP$BAF.test, data$SNP$BAF.ref, flavor = "v4", preserveScale = FALSE))
    data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = BTB, newname = "BAF.test.TB")
  }

  ### Getting heterozygous probes from Ref
  Ref.hetero <- data$SNP$mBAF.ref >= BAF.hetmin
  if (!any(Ref.hetero)) stop(tmsg("All SNP positions were tagged as homozygous in Ref : there may be a problem with your reference BAM ploidy !"), call. = FALSE)

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

  #### Getting renorm.data bins in data 
  rind <- renorm.data$tracks$pos$ProbeSetName %in% data$RD$bin
  ### L2R : Filtering
  #### 1) low depth
  rd.ori <- nrow(data$RD)
  RDlow <- data$RD$tot_count.test < RD.tot.min | data$RD$tot_count.ref < RD.tot.min
  data$meta$HTS$Imputed.lowdepth.bins <- length(which(RDlow))
  if (length(which(RDlow)) == nrow(data$RD)) stop(tmsg("All RD bins were flagged for their low read count ! You may consider lowering the BAF.tot.min value."), call. = FALSE)
  tmsg(paste0("Flagged ", length(which(RDlow)), " (", round(length(which(RDlow)) / rd.ori * 100, digits = 2), "%) RD bins with low depth (<", RD.tot.min, ")"))
  FLAGS <- RDlow > 0
  

  #### 2) GC% outliers
  if ('gc' %in% tolower(use.tracks)) {
    # GCOL <- renorm.data$tracks$gc[,1] < 200 | renorm.data$tracks$gc[,1] > 800
    GCOL <- renorm.data$tracks$gc[rind,1] < 200 | renorm.data$tracks$gc[rind,1] > 800
    data$meta$HTS$Imputed.GCoutlier.bins <- length(which(GCOL))
    if (length(which(GCOL)) == nrow(data$RD)) stop(tmsg("All RD bins were flagged as GC% outliers  ! There may be something wrong with your reference genome and/or capture BED."), call. = FALSE)
    tmsg(paste0("Flagged ", length(which(GCOL)), " (", round(length(which(GCOL)) / rd.ori * 100, digits = 2), "%) RD bins as GC% outliers."))
    FLAGS <- RDlow + GCOL > 0
  }

  ### Imputing ?
  if (any(FLAGS)) {
    tmsg(paste0(" Imputed ", length(which(FLAGS)), " (", round(length(which(FLAGS))/rd.ori*100, digits = 2), "%) L2R bins."))
    l2r.tmp <- data$RD$L2R
    l2r.tmp[FLAGS] <- NA
    data$RD$L2R <- data$RD$L2R.imp <- stats::approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  } else data$RD$L2R.imp <- data$RD$L2R
  gc()

  ## L2R : Normalization
  smo <- round(nrow(data$RD) / 550)
  if(smo%%2 == 0) smo <- smo+1

  # ### Wave
  # if (wave.renorm) {
  #   tmsg("Wave normalization ...")
  # 
  #   l2r2norm <- data.frame(ProbeSetName = data$RD$bin, chr = as.character(data$RD$chr), pos = data$RD$start, L2R = data$RD$L2R)
  #   # rownames(l2r2norm) <- seq_len(nrow(l2r2norm))
  #   ren.res <- renorm.go(input.data = l2r2norm, renorm.rda = wave.rda, track.type = "Wave", smo = smo, arraytype = data$meta$basic$type, genome = genome)
  # 
  #   fitted.l2r <- ren.res$renorm$l2r$l2r
  # 
  #   GCF <- is.na(fitted.l2r)
  #   if (any(GCF)) {
  #     l2r.tmp <- fitted.l2r
  #     l2r.tmp[GCF] <- NA
  #     fitted.l2r <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  #   }
  # 
  #   if(is.null(ren.res$renorm$pos)) {
  #     # meta.b <- setmeta("gc.renorm", "None", meta.b)
  #     data$meta$WES <- setmeta("wave.renorm", "None", data$meta$WES)
  #     tmsg(" No positive fit.")
  #   } else {
  #     ## Tweaking sex chromosomes
  #     sex.idx <- data$RD$chr %in% sex.chr
  #     auto.ori.med <- median(data$RD$L2R[!sex.idx], na.rm = TRUE)
  #     auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
  #     if (any(sex.idx)) {
  #       for (k in sex.chr) {
  #         k.idx <- data$RD$chr == k
  #         if (any(k.idx)) {
  #           k.ori.diffmed <- median(data$RD$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
  #           k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
  #           fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
  #         }
  #       }
  #     }
  #     # meta.b <- setmeta("wave.renorm", paste0(ren.res$renorm$pos, collapse = ","), meta.b)
  #     data$meta$WES <- setmeta("wave.renorm", paste0(ren.res$renorm$pos, collapse = ","), data$meta$WES)
  #   }
  #   rm(ren.res)
  # 
  #   data$RD$L2R.WAVE <- data$RD$L2R <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
  # } else {
  #   # meta.b <- setmeta("wave.renorm", "FALSE", meta.b)
  #   data$meta$WES <- setmeta("wave.renorm", "FALSE", data$meta$WES)
  # }

  #### GC%
  # message("GC% normalization ...")
  # data$RD$L2R.GC <- data$RD$L2R <- limma::loessFit(x = data$RD$GCPC, y = data$RD$L2R)$residuals
  # if (any(is.na(data$RD$L2R))) {
  #   l2r.tmp <- data$RD$L2R
  #   l2r.tmp[is.na(data$RD$L2R)] <- NA
  #   data$RD$L2R <- data$RD$L2R.GC <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  # }

  # if (gc.renorm) {
  #   tmsg("Multi-tracks normalization ...")
  # 
  #   l2r2norm <- data.frame(ProbeSetName = data$RD$bin, chr = as.character(data$RD$chr), pos = data$RD$start, L2R = data$RD$L2R)
  #   # rownames(l2r2norm) <- seq_len(nrow(l2r2norm))
  #   ren.res <- renorm.go(input.data = l2r2norm, renorm.rda = BINpack, track.type = "GC", smo = smo, arraytype = data$meta$basic$type, genome = genome)
  # 
  #   fitted.l2r <- ren.res$renorm$l2r$l2r
  # 
  #   GCF <- is.na(fitted.l2r)
  #   if (any(GCF)) {
  #     l2r.tmp <- fitted.l2r
  #     l2r.tmp[GCF] <- NA
  #     fitted.l2r <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  #   }
  # 
  #   if(is.null(ren.res$renorm$pos)) {
  #     # meta.b <- setmeta("gc.renorm", "None", meta.b)
  #     data$meta$eacon <- setmeta("gc.renorm", "None", data$meta$eacon)
  #     tmsg(" No positive fit.")
  #   } else {
  #     ## Tweaking sex chromosomes
  #     sex.idx <- data$RD$chr %in% sex.chr
  #     auto.ori.med <- median(data$RD$L2R.ori[!sex.idx], na.rm = TRUE)
  #     auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
  #     if (any(sex.idx)) {
  #       for (k in sex.chr) {
  #         k.idx <- data$RD$chr == k
  #         if (any(k.idx)) {
  #           # k.ori.diffmed <- median(data$RD$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
  #           k.ori.diffmed <- median(data$RD$L2R.ori[k.idx], na.rm = TRUE) - auto.ori.med
  #           k.rn.diffmed <- median(fitted.l2r[k.idx], na.rm = TRUE) - auto.rn.med
  #           fitted.l2r[k.idx] <- fitted.l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
  #         }
  #       }
  #     }
  #     # meta.b <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, collapse = ","), meta.b)
  #     data$meta$eacon <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, collapse = ","), data$meta$eacon)
  #   }
  #   rm(ren.res)
  # 
  #   data$RD$L2R.GC <- data$RD$L2R <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
  # } else {
  #   # meta.b <- setmeta("gc.renorm", "FALSE", meta.b)
  #   data$meta$eacon <- setmeta("gc.renorm", "FALSE", data$meta$eacon)
  # }
  
  ## Multi-track normalization
  for (tt in use.tracks) {
  
    tmsg(paste0('Multi-tracks (', tt, ') normalization ...'))
    
    l2r2norm <- data.frame(ProbeSetName = data$RD$bin, chr = as.character(data$RD$chr), pos = data$RD$start, L2R = data$RD$L2R)
    ren.res <- renorm.go(input.data = l2r2norm, renorm.data = renorm.data, track.type = tt, method = method, smo = smo, arraytype = data$meta$basic$type, genome = genome, sex.chr = sex.chr)
    
    fitted.l2r <- ren.res$renorm$l2r
    
    # GCF <- is.na(fitted.l2r)
    # if (any(GCF)) {
    #   l2r.tmp <- fitted.l2r
    #   l2r.tmp[GCF] <- NA
    #   fitted.l2r <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
    # }
    
    if(is.null(ren.res$renorm$pos)) {
      # meta.b <- setmeta("gc.renorm", "None", meta.b)
      data$meta$eacon <- setmeta(paste0(tt, ".renorm"), "None", data$meta$eacon)
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
      data$meta$eacon <- setmeta(paste0(tt, ".renorm"), paste0(ren.res$renorm$pos, collapse = ","), data$meta$eacon)
    }
    rm(ren.res)
    data$RD$L2R <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
  }
  # else {
  #   # meta.b <- setmeta("gc.renorm", "FALSE", meta.b)
  #   data$meta$eacon <- setmeta("gc.renorm", "FALSE", data$meta$eacon)
  # }

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
    saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "_", data$meta$basic$genome, "_b", data$meta$HTS$bin.size, "_processed.RDS"), compress = "xz")
  }

  ## Plot
  tmsg("Plotting ...")
  if (plot) {
    l2r <- my.ascat.obj$data$Tumor_LogR[,1]
    l2r.rm <- suppressWarnings(runmed(l2r, smo))
    # l2r.dif <- diff(l2r)
    l2r.mad <- MAD.scorer(l2r)
    # l2r.rm.dif <- diff(l2r.rm)
    l2r.ssad <- SSAD.scorer(x = l2r, smo = smo)
    
    l2r.ori <- my.ascat.obj$data$Tumor_LogR.ori[,1]
    l2r.ori.rm <- suppressWarnings(runmed(l2r.ori, smo))
    l2r.ori.rm <- l2r.ori.rm - median(l2r.ori.rm, na.rm = TRUE)
    # l2r.ori.dif <- diff(l2r.ori)
    l2r.ori.mad <- MAD.scorer(x = l2r.ori)
    # l2r.ori.mad <- median(abs(l2r.ori.dif[l2r.ori.dif != 0]))
    # l2r.ori.rm.dif <- diff(l2r.ori.rm)
    l2r.ori.ssad <- SSAD.scorer(x = l2r.ori, smo = smo)
    
    l2r.genopos <- my.ascat.obj$data$SNPpos$pos + cs$chromosomes$chr.length.toadd[my.ascat.obj$data$SNPpos$chrs]
    
    l2r <- l2r - median(l2r, na.rm = TRUE)
    l2r.ori <- l2r.ori - median(l2r.ori, na.rm = TRUE)
    kend <- l2r.genopos[vapply(unique(my.ascat.obj$data$SNPpos$chr), function(k) { max(which(my.ascat.obj$data$SNPpos$chrs == k))}, 1)]
    
    png(paste0(out.dir, "/", samplename, "_HTS_", data$meta$basic$genome, "_rawplot.png"), 1600, 1050)
    par(mfrow = c(3,1), mar = c(5, 5, 4, 2) + 0.1)
    plot(l2r.genopos, l2r.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " HTS (", data$meta$basic$manufacturer, ") raw L2R profile (median-centered)\nMAD = ", round(l2r.ori.mad, digits = 3), " ; SSAD = ", round(l2r.ori.ssad, digits = 3)), xlab = "Genomic position", ylab = "L2R", cex.axis = 2, cex.lab = 2, cex.main = 2)
    points(x = l2r.genopos, y = l2r.ori.rm, pch = ".", cex = 5, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(l2r.genopos, l2r, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " HTS (", data$meta$basic$manufacturer, ") normalized L2R profile (median-centered)\nMAD = ", round(l2r.mad, digits = 3), " ; SSAD = ", round(l2r.ssad, digits = 3)), xlab = "Genomic position", ylab = "L2R", cex.axis = 2, cex.lab = 2, cex.main = 2)
    points(x = l2r.genopos, y = l2r.rm, pch = ".", cex = 5, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(l2r.genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " HTS (", data$meta$basic$manufacturer, ")", if(TumorBoost) " TumorBoost-normalized", " BAF profile"), xlab = "Genomic position", ylab = "BAF", cex.axis = 2, cex.lab = 2, cex.main = 2)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = .5, col = 2, lty = 2, lwd = 2)
    dev.off()
  }

  tmsg("Done.")
  if(return.data) return(my.ascat.obj)
}

## Add BAF tracks to the main object
BAF.adder <- function(Bdata = NULL, Bvalues = NULL, newname = NULL, type = "test") {
    Bdata[[paste0("BAF.", type)]] <- Bdata[[newname]] <- Bvalues
    Bdata[[paste0("mBAF.", type)]] <- BAF2mBAF(Bdata[[paste0("BAF.", type)]])
    return(Bdata)
  }

## Tweak gonosomes l2r profile that may be affected by regression (normalization) by re-establishing their median pre-norm levels
sex_tweaker <- function(l2r.ori = NULL, l2r = NULL, chr = NULL, sex.chr = c('chrX', 'chrY')) {
  sex.idx <- chr %in% sex.chr
  auto.ori.med <- median(l2r.ori[!sex.idx], na.rm = TRUE)
  auto.rn.med <- median(l2r[!sex.idx], na.rm = TRUE)
  if (any(sex.idx)) {
    for (k in sex.chr) {
      k.idx <- chr == k
      if (any(k.idx)) {
        k.ori.diffmed <- median(l2r.ori[k.idx], na.rm = TRUE) - auto.ori.med
        k.rn.diffmed <- median(l2r[k.idx], na.rm = TRUE) - auto.rn.med
        l2r[k.idx] <- l2r[k.idx] - k.rn.diffmed + k.ori.diffmed
      }
    }
  }
  return(l2r)
}

## Performs the normalization of HTS L2R and BAF signals
HTS.Normalize <- function(data = NULL, BINpack = NULL, use.tracks = c('gc'), method = 'spline', RD.tot.min = 20, RD.alt.min = 3, BAF.hetmin = .33, sex.chr = c("chrX", "chrY"), TumorBoost = FALSE, out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE) {
 
  ### AJOUTER UN CONTROLE DU RDS (pour que ce ne soit pas un _processed.RDS donné en entrée!)
  
  
  ## CHECKS
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  if (is.null(BINpack)) stop(tmsg("A BINpack file is required !"), call. = FALSE)
  if (!file.exists(BINpack)) stop(tmsg("Could not find the BINpack file !"), call. = FALSE)
  if (is.null(use.tracks)) stop('No track type to use (from BINPack) provided : normalization is not possible !')
  if (RD.tot.min < 0) stop(tmsg("RD.tot.min must be >= 0 !"), call. = FALSE)
  if (RD.alt.min <= 0) stop(tmsg("RD.alt.min must be > 0 !"), call. = FALSE)
  
  ## Fix bad order of bins (dunno why, maybe the out of order mode of foreach at the pairing step)
  data$RD <- data$RD[order(data$RD$bin),]
  
  ## TAGS
  data$meta$HTS$TumorBoost <- as.character(TumorBoost)
  data$meta$HTS$RD.tot.min <- RD.tot.min
  data$meta$HTS$renorm.use.tracks <- use.tracks
  samplename <- data$meta$basic$samplename
  
  ## Loading BINpack
  load(BINpack)
  gc()
  
  ## Checks from BINpack
  avail.tracks <- names(renorm.data$tracks)
  avail.tracks <- avail.tracks[!avail.tracks == 'pos']
  
  if (!all(use.tracks %in% avail.tracks)) stop(paste0('Available normalization track type(s) is/are : ["', paste(avail.tracks, collapse = '", "'), '"], while requested track type(s) is/are : ["', paste(use.tracks, collapse = '", "'), '"].'))
  
  genome.pkg <- renorm.data$info$value[renorm.data$info$key == "genome-package"]
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  
  ### Loading genome
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  # requireNamespace(genome.pkg, quietly = TRUE)
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  cs <- chromobjector(BSg.obj)
  
  ## BAF HANDLING
  tmsg("Processing BAF ...")
  
  # BAF.adder <- function(Bdata = NULL, Bvalues = NULL, newname = NULL, type = "test") {
  #   Bdata[[paste0("BAF.", type)]] <- Bdata[[newname]] <- Bvalues
  #   Bdata[[paste0("mBAF.", type)]] <- BAF2mBAF(Bdata[[paste0("BAF.", type)]])
  #   return(Bdata)
  # }
  
  ### BAF : Filtering (low depth)
  rd.ori <- nrow(data$SNP)
  
  #### 1) Filtering for low ref or test depth
  RDlow <- data$SNP$tot_count.test < RD.tot.min | data$SNP$tot_count.ref < RD.tot.min
  if (length(which(RDlow)) == nrow(data$SNP)) stop(tmsg("All SNP positions were discarded for their low read count ! You may consider lowering the BAF.tot.min value."), call. = FALSE)
  tmsg(paste0("Removed ", length(which(RDlow)), " (", round(length(which(RDlow)) / rd.ori * 100, digits = 2), "%) SNP positions with low depth (<", RD.tot.min, ")"))
  data$SNP <- data$SNP[!RDlow,]
  
  #### 2) Filtering for low alt test depth
  BRDlow <- data$SNP$alt_count.test < RD.alt.min
  if (length(which(BRDlow)) == nrow(data$SNP)) stop(tmsg("All SNP positions were discarded for their low alternative allele count ! You may consider lowering the RD.alt.min value."), call. = FALSE)
  tmsg(paste0("Removed ", length(which(BRDlow)), " (", round(length(which(BRDlow)) / rd.ori * 100, digits = 2), "%) SNP positions with low alt RD (<", RD.alt.min, ")"))
  data$SNP <- data$SNP[!BRDlow,]
  gc()
  
  ### Computing BAF
  data$SNP$BAF.test.ori <- data$SNP$alt_count.test / data$SNP$tot_count.test
  data$SNP$BAF.ref.ori <- data$SNP$alt_count.ref / data$SNP$tot_count.ref
  odd.idx <- which(data$SNP$pos %% 2 == 1)
  data$SNP$BAF.test.ori[odd.idx] <- -data$SNP$BAF.test.ori[odd.idx] +1L
  data$SNP$BAF.ref.ori[odd.idx] <- -data$SNP$BAF.ref.ori[odd.idx] +1L
  
  ### Test
  data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = data$SNP$BAF.test.ori, newname = "BAF.test.ori")
  ### Ref
  data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = data$SNP$BAF.ref.ori, newname = "BAF.ref.ori", type = "ref")
  
  ###### Computing LOR (and variance)
  rcmat <- round(cbind(data$SNP$BAF.test.ori*data$SNP$tot_count.test, (1-data$SNP$BAF.test.ori)*data$SNP$tot_count.test))
  data$SNP$LOR <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
  data$SNP$LORvar <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
  rm(rcmat)
  gc()
  
  ### BAF : TumorBoost
  if (TumorBoost) {
    ### Apply TumorBoost
    tmsg('Applying TumorBoost BAF normalization ...')
    BTB <- as.numeric(aroma.light::normalizeTumorBoost(data$SNP$BAF.test, data$SNP$BAF.ref, flavor = "v4", preserveScale = FALSE))
    ### Fix for weird shift sometimes
    BTB <- BTB - median(BTB, na.rm = TRUE) + .5
    data$SNP <- BAF.adder(Bdata = data$SNP, Bvalues = BTB, newname = "BAF.test.TB")
  }
  
  ### Getting heterozygous probes from Ref
  Ref.hetero <- data$SNP$mBAF.ref >= BAF.hetmin
  if (!any(Ref.hetero)) stop(tmsg("All SNP positions were tagged as homozygous in Ref : there may be a problem with your reference BAM ploidy !"), call. = FALSE)
  
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
  
  #### Getting renorm.data bins in data 
  
  RinI <- renorm.data$tracks$pos$ProbeSetName %in% data$RD$bin
  # rind <- renorm.data$tracks$pos$ProbeSetName %in% data$RD$bin
  ### L2R : Filtering
  #### 1) low depth
  rd.ori <- nrow(data$RD)
  RDlow <- data$RD$tot_count.test < RD.tot.min | data$RD$tot_count.ref < RD.tot.min
  data$meta$HTS$Imputed.lowdepth.bins <- length(which(RDlow))
  if (length(which(RDlow)) == nrow(data$RD)) stop(tmsg("All RD bins were flagged for their low read count ! You may consider lowering the BAF.tot.min value."), call. = FALSE)
  tmsg(paste0("Flagged ", length(which(RDlow)), " (", round(length(which(RDlow)) / rd.ori * 100, digits = 2), "%) RD bins with low depth (<", RD.tot.min, ")"))
  FLAGS <- RDlow > 0
  
  
  #### 2) GC% outliers
  if ('gc' %in% tolower(use.tracks)) {
    GCOL <- renorm.data$tracks$gc[RinI,1] < 200 | renorm.data$tracks$gc[RinI,1] > 800
    data$meta$HTS$Imputed.GCoutlier.bins <- length(which(GCOL))
    if (length(which(GCOL)) == nrow(data$RD)) stop(tmsg("All RD bins were flagged as GC% outliers  ! There may be something wrong with your reference genome and/or capture BED."), call. = FALSE)
    tmsg(paste0("Flagged ", length(which(GCOL)), " (", round(length(which(GCOL)) / rd.ori * 100, digits = 2), "%) RD bins as GC% outliers."))
    FLAGS <- RDlow + GCOL > 0
  }
  
  ### Imputing ?
  if (any(FLAGS)) {
    tmsg(paste0(" Imputed ", length(which(FLAGS)), " (", round(length(which(FLAGS))/rd.ori*100, digits = 2), "%) L2R bins."))
    l2r.tmp <- data$RD$L2R
    l2r.tmp[FLAGS] <- NA
    data$RD$L2R <- stats::approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  }
  gc()
  
  ## L2R : Normalization
  smo <- round(nrow(data$RD) / 550)
  if(smo%%2 == 0) smo <- smo+1
  
  ## Sync and NA checks 
  # RinI <- renorm.data$tracks$pos$ProbeSetName %in% data$RD$bin
  if (!all(unique(renorm.data$tracks$pos$ProbeSetName[RinI] == data$RD$bin))) stop(tmsg(paste0(track.type, " data and L2R data are not synched, or ordered differently !")), call. = FALSE)
  # nona <- !is.na(data$RD$L2R)
  
  ## Initialization
  init.marmmm <- MARMMM.scorer(x = data$RD$L2R.ori, smo = smo)
  tmsg(paste0('Global init, no filtering (', init.marmmm, ')'))
  
  ## Multi-track normalization
  tflag <- FALSE
  while(length(use.tracks) > 0 & !tflag) {
    ### Select tracks order
    tmsg(paste0('Selecting optimal track type in ["', paste(use.tracks, collapse = '", "'), '"] ...'))
    l2r.tmp <- data$RD$L2R
    init.marmmm <- MARMMM.scorer(x = l2r.tmp, smo = smo)
    tscorelist <- sapply(use.tracks, function(ut) {
      if (method == "loess") tempfit <- l2r.loess(l2r = l2r.tmp, tf = renorm.data$tracks[[ut]][RinI,1])
      if (method == "spline") tempfit <- l2r.spline(l2r = l2r.tmp, tf = renorm.data$tracks[[ut]][RinI,1])
      if (method == "pcs") tempfit <- l2r.pcs(l2r = l2r.tmp, tf = renorm.data$tracks[[ut]][RinI,1])
      ## Tweaking sex chromosomes
      tempfit <- sex_tweaker(l2r.ori = data$RD$L2R.ori, l2r = tempfit, chr = data$RD$chr, sex.chr = sex.chr)
      ## Center
      tempfit <- tempfit - median(tempfit, na.rm = TRUE)
      ## Compute new MARMMM
      z.marmmm <- MARMMM.scorer(x = tempfit, smo = smo)
      return(list(track.type = ut, ntracks = ncol(renorm.data$tracks[[ut]]), marmmm = z.marmmm, l2r = tempfit))
    }, simplify = FALSE)
    tscore <- setNames(object = vapply(tscorelist, function(x) { x$marmmm }, .1), nm = names(tscorelist))
    opti_tt <- which.min(tscore)
    tt <- names(opti_tt)
    
    # tt
    # tmsg(str(tscore))
    
    # plot(data$RD$start[data$RD$chr == 'chr16'], tscorelist[[names(opti_tt)]]$l2r[data$RD$chr == 'chr16'], pch = '.', xaxs = 'i', main = names(opti_tt))
    
    
    if (tscore[opti_tt] > init.marmmm) {
      ## No regression better than original L2R profile
      tmsg(' Regression did not provide a better profile on any assessed track type !')
      tflag <- TRUE
    } else {
      ## Some regression to test !
      if (tscorelist[[tt]]$ntracks == 1) {
        ## Track type with single track : regression already done !
        tmsg(paste0(' Single track ("', tt, '") found as optimal !'))
        tmsg(paste0('  Positive fit with ', colnames(renorm.data$tracks[[tt]]), ' (', tscorelist[[opti_tt]]$marmmm, ')'))
        ## Remove track type used
        use.tracks <- use.tracks[!use.tracks == tt]
        data$RD$L2R <- tscorelist[[tt]]$l2r
        rm(tscorelist, opti_tt, tt)
      } else {
        ## Some recursive regression(s) to perform !
        rm(tscorelist, tscore)
        tmsg(paste0(' Multi-tracks (', tt, ') normalization ...'))
        l2r2norm <- data.frame(ProbeSetName = data$RD$bin, chr = as.character(data$RD$chr), pos = data$RD$start, L2R = data$RD$L2R)
        ren.res <- renorm.go(input.data = l2r2norm, renorm.data = renorm.data, track.type = tt, method = method, smo = smo, arraytype = data$meta$basic$type, genome = genome, sex.chr = sex.chr)
        
        # plot(data$RD$start[data$RD$chr == 'chr16'], ren.res$data$L2R[data$RD$chr == 'chr16'], pch = '.', xaxs = 'i', main = ren.res$renorm$pos)
        
        fitted.l2r <- ren.res$renorm$l2r
        if(is.null(ren.res$renorm$pos)) {
          data$meta$eacon <- setmeta(paste0(tt, ".renorm"), "None", data$meta$eacon)
          tmsg("  No positive fit.")
        } else {
          ## Tweaking sex chromosomes
          fitted.l2r <- sex_tweaker(l2r.ori = data$RD$L2R.ori, l2r = fitted.l2r, chr = data$RD$chr, sex.chr = sex.chr)
          data$meta$eacon <- setmeta(paste0(tt, ".renorm"), paste0(ren.res$renorm$pos, collapse = ","), data$meta$eacon)
        }
        
        # plot(data$RD$start[data$RD$chr == 'chr16'], fitted.l2r[data$RD$chr == 'chr16'], pch = '.', xaxs = 'i', main = ren.res$renorm$pos)
        
        data$RD$L2R <- fitted.l2r - median(fitted.l2r, na.rm = TRUE)
        use.tracks <- use.tracks[!use.tracks == tt]
        
        # plot(data$RD$start[data$RD$chr == 'chr16'], data$RD$L2R[data$RD$chr == 'chr16'], pch = '.', xaxs = 'i', main = ren.res$renorm$pos)
        rm(ren.res)
        
      }
    }
  }

  ## Merging
  # data$RD$BAF <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "BAF.test"))], by = c("chr", "bin"))$BAF.test
  data$RD$BAF <- dplyr::left_join(data$RD[, c(1:4)], data$SNP[, c(1,3,which(colnames(data$SNP) == "BAF.test"))], by = c("chr", "bin"))$BAF.test
  # data$RD$LOR <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LOR"))], by = c("chr", "bin"))$LOR
  data$RD$LOR <- dplyr::left_join(data$RD[, c(1:4)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LOR"))], by = c("chr", "bin"))$LOR
  # data$RD$LORvar <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LORvar"))], by = c("chr", "bin"))$LORvar
  data$RD$LORvar <- dplyr::left_join(data$RD[, c(1:4)], data$SNP[, c(1,3,which(colnames(data$SNP) == "LORvar"))], by = c("chr", "bin"))$LORvar
  # data$RD$RD.test <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "tot_count.test"))], by = c("chr", "bin"))$tot_count.test
  data$RD$RD.test <- dplyr::left_join(data$RD[, c(1:4)], data$SNP[, c(1,3,which(colnames(data$SNP) == "tot_count.test"))], by = c("chr", "bin"))$tot_count.test
  # data$RD$RD.ref <- dplyr::left_join(data$RD[, c(1:4,9)], data$SNP[, c(1,3,which(colnames(data$SNP) == "tot_count.ref"))], by = c("chr", "bin"))$tot_count.ref
  data$RD$RD.ref <- dplyr::left_join(data$RD[, c(1:4)], data$SNP[, c(1,3,which(colnames(data$SNP) == "tot_count.ref"))], by = c("chr", "bin"))$tot_count.ref
  
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
  colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  # rm(my.ch, data)
  gc()
  
  # plot(ares$Tumor_LogR[,1], pch = ".", cex = 3, xaxs = "i", ylim = c(-2,2))
  # points(ares$Tumor_LogR_segmented, pch = ".", cex = 3, col = 2)
  # plot(ares$Tumor_BAF[!is.na(ares$Tumor_BAF),1], pch = ".", xaxs = "i", cex = 3)
  # points(ares$Tumor_BAF_segmented[[1]], pch = ".", cex = 3, col = 2)
  # points(1 - ares$Tumor_BAF_segmented[[1]], pch = ".", cex = 3, col = 2)
  #
  
  ## Plot
  tmsg("Plotting ...")
  if (plot) {
    l2r <- my.ascat.obj$data$Tumor_LogR[,1]
    l2r.rm <- suppressWarnings(runmed(l2r, smo))
    # l2r.dif <- diff(l2r)
    l2r.mad <- MAD.scorer(l2r)
    # l2r.rm.dif <- diff(l2r.rm)
    l2r.marmmm <- MARMMM.scorer(x = l2r, smo = smo)
    
    l2r.ori <- my.ascat.obj$data$Tumor_LogR.ori[,1]
    l2r.ori.rm <- suppressWarnings(runmed(l2r.ori, smo))
    l2r.ori.rm <- l2r.ori.rm - median(l2r.ori.rm, na.rm = TRUE)
    # l2r.ori.dif <- diff(l2r.ori)
    l2r.ori.mad <- MAD.scorer(x = l2r.ori)
    # l2r.ori.mad <- median(abs(l2r.ori.dif[l2r.ori.dif != 0]))
    # l2r.ori.rm.dif <- diff(l2r.ori.rm)
    l2r.ori.marmmm <- MARMMM.scorer(x = l2r.ori, smo = smo)
    
    l2r.genopos <- my.ascat.obj$data$SNPpos$pos + cs$chromosomes$chr.length.toadd[my.ascat.obj$data$SNPpos$chrs]
    
    l2r <- l2r - median(l2r, na.rm = TRUE)
    l2r.ori <- l2r.ori - median(l2r.ori, na.rm = TRUE)
    kend <- l2r.genopos[vapply(unique(my.ascat.obj$data$SNPpos$chr), function(k) { max(which(my.ascat.obj$data$SNPpos$chrs == k))}, 1)]
    
    png(paste0(out.dir, "/", samplename, "_HTS_", data$meta$basic$genome, "_rawplot.png"), 1600, 1050)
    par(mfrow = c(3,1), mar = c(5, 5, 4, 2) + 0.1)
    plot(l2r.genopos, l2r.ori, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " HTS (", data$meta$basic$manufacturer, ") raw L2R profile (median-centered)\nMAD = ", round(l2r.ori.mad, digits = 3), " ; MARMMM = ", round(l2r.ori.marmmm, digits = 5)), xlab = "Genomic position", ylab = "L2R", cex.axis = 2, cex.lab = 2, cex.main = 2)
    points(x = l2r.genopos, y = l2r.ori.rm, pch = ".", cex = 5, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(l2r.genopos, l2r, pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(-2,2), main = paste0(samplename, " HTS (", data$meta$basic$manufacturer, ") normalized L2R profile (median-centered)\nMAD = ", round(l2r.mad, digits = 3), " ; MARMMM = ", round(l2r.marmmm, digits = 5)), xlab = "Genomic position", ylab = "L2R", cex.axis = 2, cex.lab = 2, cex.main = 2)
    points(x = l2r.genopos, y = l2r.rm, pch = ".", cex = 5, col = 1)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = 0, col = 2, lty = 2, lwd = 2)
    plot(l2r.genopos, my.ascat.obj$data$Tumor_BAF[,1], pch = ".", cex = 3, col = "grey70", xaxs = "i", yaxs = "i", ylim = c(0,1), main = paste0(samplename, " HTS (", data$meta$basic$manufacturer, ")", if(TumorBoost) " TumorBoost-normalized", " BAF profile"), xlab = "Genomic position", ylab = "BAF", cex.axis = 2, cex.lab = 2, cex.main = 2)
    abline(v = kend, col = 4, lty = 3, lwd = 2)
    abline(h = .5, col = 2, lty = 2, lwd = 2)
    dev.off()
  }
  
  ## Saving data
  if (write.data) {
    tmsg("Saving normalized data ...")
    saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "_", data$meta$basic$genome, "_b", data$meta$HTS$bin.size, "_processed.RDS"), compress = "xz")
  }
  
  tmsg("Done.")
  if(return.data) return(my.ascat.obj)
}


## Runs HTS.Normalize using a RDS filename
HTS.Normalize.ff <- function(BIN.RDS.file = NULL, ...) {
  
  ## CHECKS
  if (is.null(BIN.RDS.file)) stop(tmsg("An RDS file from EaCoN::EaCoN.WES.Bin is required !"), call. = FALSE)
  if (!file.exists(BIN.RDS.file)) stop(tmsg(paste0("Could not find ", BIN.RDS.file, " .")), call. = FALSE)

  tmsg("Loading binned HTS data ...")
  my.data <- readRDS(BIN.RDS.file)
  message("Running HTS.Normalize() ...")
  HTS.Normalize(data = my.data, out.dir = dirname(BIN.RDS.file), ...)
}

## Runs HTS.Normalize.ff, batch mode
HTS.Normalize.ff.Batch <- function(BIN.RDS.files = list.files(path = getwd(), pattern = "_binned.RDS$", all.files = FALSE, full.names = TRUE, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE), nthread = 1, cluster.type = "PSOCK", ...) {
  if (length(BIN.RDS.files) == 0) stop("No file found to process !", call. = FALSE)
  message("Running HTS.Normalize.ff() in batch mode ...")
  message(paste0("Found ", length(BIN.RDS.files), " samples to process ..."))
  current.bitmapType <- getOption("bitmapType")
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  eacon.batchres <- foreach::foreach(r = seq_along(BIN.RDS.files), .inorder = TRUE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    HTS.Normalize.ff(BIN.RDS.file = BIN.RDS.files[r], ...)
  }
  parallel::stopCluster(cl)
}

bedBinner <- function(bed = NULL, bin.size = 50, nthread = 1, cluster.type = 'PSOCK') {
  
  bin.size <- as.integer(bin.size)
  cl <- if(nthread <= 1) BiocParallel::SerialParam() else parallel::makeCluster(spec = nthread, type = cluster.type, outfile = '')
  doParallel::registerDoParallel(cl)
  suppressPackageStartupMessages(require(foreach))
  `%mydo%` <- if (nthread <= 1) foreach::"%do%" else foreach::"%dopar%"
  `%do%` <- foreach::"%do%"
  k <- 0
  bed.binned <- foreach::foreach(k = unique(bed$chr), .combine = "rbind", .export = "tmsg") %mydo% {
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
  if(nthread <= 1) BiocParallel::bpstop() else parallel::stopCluster(cl)
  # parallel::stopCluster(cl)
  return(bed.binned)
}

## Get mpileup-like data from a BAM and positional (bed-like) data
bedbam2pup <- function(BamFile = NULL, bed.data = NULL, scanBamFlag = NULL, pileupParam = NULL) {
  
  ### Making pileup
  # tmsg("   Getting pileup ...")
  param.BAM <- Rsamtools::ScanBamParam(
    which = GenomicRanges::makeGRangesFromDataFrame(
      bed.data
      , seqnames.field = "chr")
    , flag = scanBamFlag)
  pres <- tibble::as_tibble(Rsamtools::pileup(BamFile, scanBamParam = param.BAM, pileupParam = pileupParam))
  if (nrow(pres) == 0) {
    ## No matched read ?
    out_tbl <- tibble::tibble(chr = factor(), pos = integer(), bin = integer(), tot_count = integer(), alt_count = integer())
    # out_tbl <- tibble::tibble(chr = bed.data$chr, pos = bed.data$start, bin = bed.data$ProbeSetName, tot_count = 0, alt_count = 0)
    return(out_tbl)
  } else {
    ## There are counts
    colnames(pres)[c(1,5)] <- c("chr", "bin")
    levels(pres$bin) <- bed.data$ProbeSetName
    pres$bin <- as.integer(as.character(pres$bin))
    
    if("=" %in% unique(pres$nucleotide)) tmsg("   Bam contains the '=' sign !")
    
    ### Building genomic sequence of reference block
    # tmsg("   Getting reference genome sequence ...")
    refblock <- pres[!duplicated(pres$pos), c(1,2)]
    pres <- dplyr::group_by(pres, pos)
    refblock$tot_count <- dplyr::summarize(pres, tot_count = sum(count))$tot_count
    pres <- dplyr::ungroup(pres)
    refblock$nucleotide <- as.factor(suppressWarnings(BSgenome::getSeq(BSg.obj, names = GenomicRanges::makeGRangesFromDataFrame(refblock, start.field = "pos", end.field = "pos"), as.character = TRUE)))
    
    ### Merging blocks
    # tmsg("   Computing alt counts ...")
    merged <- suppressWarnings(dplyr::left_join(refblock, pres, by = c("chr", "pos", "nucleotide")))
    rm(pres, refblock)
    gc()
    bed.based <- tibble::as_tibble(data.frame(chr = unique(bed.data$chr), pos = as.integer(unlist(EaCoN:::seq.int2(from = bed.data$start, to = bed.data$end, by = 1))), bin = rep(bed.data$ProbeSetName, times = (bed.data$end - bed.data$start +1))))
    bed.joint <- suppressWarnings(dplyr::left_join(bed.based, merged, c("chr", "pos", "bin"))) ### yeah !
    rm(bed.based)
    gc()
    ### Cleaning and formatting
    bed.joint$nucleotide <- NULL
    bed.joint$chr <- as.factor(bed.joint$chr)
    
    bed.joint$count <- bed.joint$tot_count - bed.joint$count
    colnames(bed.joint) <- c("chr", "pos", "bin", "tot_count", "alt_count")
    
    bed.joint$tot_count[is.na(bed.joint$tot_count)] <- 0L
    bed.joint$alt_count[is.na(bed.joint$alt_count)] <- 0L
    
    return(bed.joint)
  }
}

## Retrieve mpileup-like data from a BAM connection, in data blocks of limited size
pileup.go <- function(BAM = NULL, genome.pkg = NULL, bed.data = NULL, scanBamFlag = NULL, pileupParam = NULL, nsubthread = 1, cluster.type = "PSOCK", blocksize = 1000) {
  
  # BamFile <- testBAM
  # scanBamFlag <- param.FLAG
  # pileupParam <- param.PILEUP
  # nsubthread <- 1
  # cluster.type = "PSOCK"
  
  ## Loading genome
  message(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  
  #### Indexing BAM if needed
  if (!(file.exists(paste0(BAM, ".bai")) || file.exists(sub(pattern = "\\.bam$", replacement = ".bai", x = BAM, ignore.case = TRUE)))) {
    tmsg("Indexing BAM ...")
    Rsamtools::indexBam(BAM) 
  } else tmsg("BAM is already indexed.")
  
  #### Launching cluster
  if (length(unique(bed.data$chr)) < nsubthread) nsubthread <- length(unique(bed.data$chr))
  cl <- if(nsubthread <= 1) BiocParallel::SerialParam() else parallel::makeCluster(spec = nsubthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  `%mydo%` <- if (nsubthread <= 1) foreach::"%do%" else foreach::"%dopar%"
  
  ## Creating step-vector for block processing (limiting data to blocksize)
  bd.tbl <- table(bed.data$chr)
  instep <- lapply(names(bd.tbl), function(k) as.numeric(cut(seq_len(nrow(bed.data[bed.data$chr == k,])), seq.int(0, nrow(bed.data[bed.data$chr == k,]) + blocksize, blocksize))))
  for (k in 2:length(instep)) instep[[k]] <- instep[[k]] + max(instep[[k-1]])
  instep <- unlist(instep)
  
  s <- 0
  BAMcounts <- foreach::foreach(s = sort(unique(instep)), .inorder = FALSE, .export = c("tmsg", "BSg.obj", "bedbam2pup", "seq.int2")) %mydo% {
  # BAMcounts <- foreach::foreach(s = c(1,max(instep)), .inorder = TRUE, .export = c("tmsg", "BSg.obj", "bedbam2pup", "seq.int2")) %mydo% {
    
    # BAMcounts <- lapply(sort(unique(instep))[1:2], function(s) {
    # tmsg(paste0(' Block : ', s, ' / ', length(unique(instep))))
    
    ### Computing counts for BAM
    bed.data.s <- bed.data[instep == s,]
    
    ## Opening BAM connections
    openBAM <- Rsamtools::BamFile(BAM)
    
    mPUP <- bedbam2pup(BamFile = openBAM, bed.data = bed.data.s, scanBamFlag = scanBamFlag, pileupParam = pileupParam)
    mPUP <- dplyr::group_by(mPUP, bin)
    
    ### Binning
    # tmsg("  Binning ...")
    CN.table <- suppressWarnings(dplyr::summarize(mPUP, chr = unique(chr), start = min(pos), end = max(pos), tot_count = round(mean(tot_count))))
    mPUP <- dplyr::ungroup(mPUP)
    SNP.table <- mPUP[mPUP$alt_count > 0,]
    
    ### Cleaning
    rm(mPUP)
    
    CN.table <- dplyr::arrange(CN.table, chr, start, end)
    CN.table <- dplyr::select(CN.table, chr, start, end, bin, tot_count)
    colnames(SNP.table) <- c("chr", "pos", "bin", "tot_count", "alt_count")
    
    return(list(id = s, CN = CN.table, SNP = SNP.table))
  }
  if(nsubthread <= 1) BiocParallel::bpstop() else parallel::stopCluster(cl)
  return(BAMcounts)
}

## Compute letter composition of nucleotidic sequences from a bed-like (chr, start, end) dataframe, with possible extension.
loc.nt.count.hs <- function(loc.df = NULL, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", extend = 0, blocksize = 1E+04, nthread = 5, cluster.type = 'PSOCK') {
  if (is.null(loc.df)) stop("loc.df is required !", call. = FALSE)
  
  if (extend < 0) stop("extend should be >= 0", call. = FALSE)
  if (blocksize <= 0) stop("blocksize should be > 0", call. = FALSE)
  if (!all(is.character(loc.df$chr) | is.factor(loc.df$chr))) stop("chr should be character !", call. = FALSE)
  if (!all(is.numeric(loc.df$start))) stop("start should be numeric !", call. = FALSE)
  if (!all(is.numeric(loc.df$end))) stop("end should be numeric !", call. = FALSE)

  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }

  message(paste0("Loading ", genome.pkg, " sequence ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- S4Vectors::metadata(BSg.obj)$genome
  cs <- chromobjector(BSg.obj)

  message("Removing replicated locations ...")
  idz <- paste0(loc.df$chr, ":", loc.df$start, "-", loc.df$end)
  loc.df <- loc.df[!duplicated(idz),]

  message("Removing non-canonical sequences ...")
  loc.df <- loc.df[loc.df$chr %in% seqnames(BSg.obj),]

  message("Ordering data ...")
  loc.df <- loc.df[order(unlist(cs$chrom2chr[loc.df$chr]), loc.df$start, loc.df$ProbeSetName),]
  
  myGR.ex <- suppressPackageStartupMessages(GenomicRanges::makeGRangesFromDataFrame(loc.df, seqinfo = seqinfo(BSg.obj)))
  if (extend > 0) myGR.ex <- GenomicRanges::trim(myGR.ex + extend)

  instep <- as.numeric(cut(seq_along(myGR.ex), seq.int(0, length(myGR.ex) + blocksize, blocksize)))

  message("Computing base composition ...")

  message("Starting cluster ...")
  
  if (length(unique(instep)) < nthread) nthread <- length(unique(instep))
  
  cl <- if(nthread <= 1) BiocParallel::SerialParam() else parallel::makeCluster(spec = nthread, type = cluster.type, outfile = '')
  
  doParallel::registerDoParallel(cl)
  
  requireNamespace("foreach", quietly = TRUE)
  
  `%mydo%` <- if (nthread <= 1) foreach::"%do%" else foreach::"%dopar%"
  
  xcounts <- foreach::foreach(x = unique(instep), .combine = "rbind", .packages = c("Biostrings", "BSgenome"), .export = "getSeq") %mydo% {
    return(Biostrings::alphabetFrequency(BSgenome::getSeq(BSg.obj, myGR.ex[which(instep == x)]), baseOnly = TRUE))
  }
  message("Stopping cluster ...")
  if(nthread <= 1) BiocParallel::bpstop() else parallel::stopCluster(cl)
  
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
    message(paste0("Computing GC +", nt.add, " ..."))
    adb.counts <- loc.nt.count.hs(loc.df = loc.df, extend = nt.add, ...)
    adb.gc <- loc.nt.gcc.hs(loc.counts = adb.counts)
    return(adb.gc)
  }
  base.df <- gc.list[[1]][,1:4]
  gc.df <- foreach::foreach(nt.add = seq_along(gc.list), .combine = "cbind") %do% {
    my.gc <- round(gc.list[[nt.add]][,"GC", drop = FALSE] * 1000)
    my.gc$GC <- as.integer(my.gc$GC)
    colnames(my.gc) <- paste0('GC', extend.multi[nt.add])
    return(my.gc)
  }
  # return(data.frame(base.df, gc.df, stringsAsFactors = FALSE))
  return(gc.df)
}

genome.build.finder <- function(BAM.header = NULL, valid.genomes = NULL) {
  BAM.header <- unlist(BAM.header)
  query <- paste0("(", paste0(valid.genomes, collapse = "|"), ")")
  stvh.grep <- grep(query, unlist(BAM.header))
  if (length(stvh.grep) == 0) stop(tmsg("Could not automatically determine genome build ! Please specify it !"), call. = FALSE)

  stvh.regexec <- unique(vapply(stvh.grep, function(x) {
    rc.res <- regexec(query, BAM.header[x])[[1]]
    return(as.character(substr(BAM.header[x], start = rc.res[1], stop = rc.res[1]+attr(rc.res, "match.length")[1]-1)))
  }, "a"))

  ok.genome <- unique(stvh.regexec[stvh.regexec %in% valid.genomes])
  if (length(ok.genome) == 0) stop(tmsg(paste0("Identified a putative genome build (", ok.genome, "), but not a supported one !")), call. = FALSE)
  if (length(ok.genome) >= 2) stop(tmsg(paste0("Identified more than one putative genome build (", ok.genome, ") !")), call. = FALSE)
  return(ok.genome)
}

## Read a cytoband(ideo) table download from UCSC Tables
read_ucsc_cytobandideo <- function(file = NULL, genomebuild = NULL, species = NULL) {
  
  ## Import cytoband data
  cytob <- utils::read.table(file = file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE, check.names = FALSE)
  colnames(cytob) <- c("chrom", "start", "end", "cytoband", "giestain")
  
  ## Filtering out invalid lines (UCSC recently added incomplete lines which seem to refer to chrN_unknown sequence sets, which is quite totally illogical. BJ 20140227)
  cytob <- cytob[!is.na(cytob$start),]
  
  ## Filtering out non-canonical chromosomes
  cytob <- cytob[grep(pattern = "^chr([0-9]+|X|Y|M)$", x = cytob$chrom),]
  
  ## Adding alternative versions of chr names
  cytob$chrN <- cytob$chrA <- sub(pattern = "^chr", replacement = "", x = cytob$chrom)
  non.auto <- c("X", "Y", "M")
  num.chr.max <- max(as.numeric(cytob$chrN[!cytob$chrN %in% non.auto]))
  for (cx in 1:length(non.auto)) cytob$chrN[cytob$chrN == non.auto[cx]] <- num.chr.max + cx
  cytob$chrN <- as.integer(cytob$chrN)
  
  ## Ordering by chromosomal number
  cytob <- cytob[order(cytob$chrN, cytob$start, cytob$end),]
  
  ## Additional data for chromosomes
  chrom <- data.frame(chrom = unique(cytob$chrom), chrA = unique(cytob$chrA), chrN = unique(cytob$chrN), stringsAsFactors = FALSE)
  chrom$chr.length <- vapply(unique(cytob$chrN), function(x) { max((cytob[which(cytob$chrN == x),])$end) }, 1)
  lchrtoadd <- c(0, chrom$chr.length[1:length(chrom$chr.length)-1])
  chrom$chr.length.sum <- cumsum(chrom$chr.length)
  chrom$chr.length.toadd <- c(0, chrom$chr.length.sum[-c(length(chrom$chr.length.sum))])
  chrom$mid.chr <- round(chrom$chr.length / 2)
  chrom$mid.chr.geno <- chrom$mid.chr + chrom$chr.length.toadd[chrom$chrN]
  chrom$centromere <- cytob$end[vapply(unique(cytob$chrN), function(k) {
    centro.pos <- which(cytob$chrN == k & cytob$giestain == "acen" & strtrim(cytob$cytoband, 1) == "p")
    return(ifelse(length(centro.pos) == 0, NA, centro.pos))
  }, 1L)]
  names(chrom$centromere) <- unique(cytob$chrN)
  glen <- sum(chrom$chr.length)
  
  ## Adding genomic coordinates for cytob
  cytob$start.geno <- cytob$start + chrom$chr.length.toadd[cytob$chrN]
  cytob$end.geno <- cytob$end + chrom$chr.length.toadd[cytob$chrN]
  
  ## Adding genomic coordinates for centromeres too
  chrom$centromere.geno <- chrom$centromere + chrom$chr.length.toadd[unique(cytob$chrN)]
  
  ## Adding converters
  chrom2chr <- as.list(unique(cytob$chrN))
  names(chrom2chr) = unique(cytob$chrom)
  chr2chrom <- as.list(unique(cytob$chrom))
  names(chr2chrom) = unique(cytob$chrN)
  
  ## Final form
  cs <- list(species = species,
             genomebuild = genomebuild,
             cytobands = cytob,
             chromosomes = chrom,
             chrom2chr = chrom2chr,
             chr2chrom = chr2chrom,
             genome.length=glen)
  return(cs)
}


## Get a quick description of a BINpack content (especially meta data). Useful to get the available types of
# HTS.BINpack.descriptor <- function(BINpack = NULL) {
#   
# }


