## ASCAT L2R + BAF segmentation
Segment.ASCAT <- function(data = NULL, mingap = 5E+06, smooth.k = NULL, BAF.filter = .75, homoCut = .05, penalty = 50, recenter = "l2r.centeredpeak", calling.method = "mad", nrf = .5, SER.pen = 40, out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  # setwd("/home/job/Documents/ROSCOFF/Roscoff_2018/TP_CNV/WES/REDUX/A18R.11.17.18")
  # data <- readRDS("A18R.11.17.18_hs37d5_b50_processed.RDS")
  # mingap = 5E+06
  # smooth.k = 5
  # BAF.filter = .75
  # homoCut = .05
  # penalty = 50
  # recenter = "l2r.centeredpeak"
  # calling.method = "mad"
  # nrf = 1
  # SER.pen = 5
  # out.dir = getwd()
  # return.data = FALSE
  # write.data = TRUE
  # plot = TRUE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/plot_functions.R")

  `%do%` <- foreach::"%do%"
  
  calling.method <- tolower(calling.method)
  
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  if (!dir.exists(out.dir)) stop(tmsg(paste0("Output directory [", out.dir, "] does not exist !")), call. = FALSE)
  if (!(calling.method %in% c("mad", "density"))) stop(tmsg("calling.method should be 'MAD' or 'density' !"), call. = FALSE)
  if (calling.method == "mad" & is.null(nrf)) stop(tmsg("If calling.method is set to 'MAD', nrf is required !"), call. = FALSE)
  if (!is.null(SER.pen)) if (!is.character(SER.pen)) if (SER.pen <= 0) stop(tmsg("SER.pen should be NULL, a character or an integer/float > 0 !"), call. = FALSE)
  
  
  ## Extract samplename
  samplename <- data$meta$basic$samplename
  tmsg(paste0("Sample : ", samplename))
  
  ## Extract genome version and load corresponding data
  genome <- data$meta$basic$genome
  genome.pkg <- data$meta$basic$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  cs <- chromobjector(BSg.obj)
  
  ## Init graphical parameters
  oridir <- getwd()
  
  ## Create output dir
  if (write.data | plot) {
    odir <- paste0(out.dir, "/ASCAT/L2R")
    # if (dir.exists(odir)) stop(tmsg(paste0("A [", odir, "] dir already exists !"))) else dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(odir)) {
      if (force) {
        unlink(odir, recursive = TRUE, force = FALSE)
        dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
      } else stop(tmsg(paste0("A [", odir, "] dir already exists !")), call. = FALSE)
    } else dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
  } else odir <- out.dir
  
  setwd(odir)
  
  data$meta$eacon <- c(data$meta$eacon, list(
    "segmenter" = "ASCAT",
    "mingap" = mingap,
    "BAF.filter" = BAF.filter,
    "BAF.segments.homo.limit" = homoCut,
    "winsorize.k" = if(is.null(smooth.k)) "NA" else smooth.k,
    "ASCAT.penalty" = penalty,
    "calling.method" = calling.method,
    "calling.nrf" = if(is.null(nrf)) "NA" else nrf,
    "small.events.rescue.PELT.penalty" = if(is.null(SER.pen)) "NA" else SER.pen
  ))
  
  ## BAF scaling
  # if (BAF.rescale) {
  #   tmsg("BAF scaling ...")
  #   data$data <- EaCoN.BAF.Scaler(ASCATobj = data$data, bafbin.size = bafbin.size, toclustname = "BAF")
  # }

  ## Germline prediction
  # tmsg("Generating germline ...")
  # # data$germline <- EaCoN.Predict.Germline(ASCATobj = data$data, prior = prior, bafbin.size = bafbin.size, modelName = "E", toclustname = "BAF", mcMin = mc.range[1], mcMax = mc.range[2], nfactor = 4, BAF.filter = BAF.filter, BAF.cutter = BAF.cutter, segmentLength = segmentLength, genome.pkg = genome.pkg)
  # data$germline <- EaCoN.Predict.Germline(ASCATobj = data$data, bafbin.size = bafbin.size, modelName = "E", toclustname = "BAF", mc.G = mc.G, nfactor = nfactor, BAF.filter = BAF.filter, BAF.cutter = BAF.cutter, segmentLength = segmentLength, genome.pkg = genome.pkg)
  # 
  
  
  ## Winsorization
  if(!is.null(smooth.k)) {
   tmsg("Smoothing L2R outliers ...")
    cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
    l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
    cndf <- cndf[l2r.nona,]
    cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", method = "mad", k = smooth.k, tau = 1, verbose = FALSE)
    data$data$Tumor_LogR[l2r.nona,1] <- cndf.wins[, 3, drop = FALSE]
    rm(list = c("cndf", "cndf.wins", "l2r.nona"))
  }
  
  ## BAF filtering
  tmsg("Filtering BAF...")
  if ("Tumor_BAF.unisomy" %in% names(data$data)) {
    called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & !is.na(data$data$Tumor_BAF[,1]) & !data$data$Tumor_BAF.unisomy[,1])
  } else {
    called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & !is.na(data$data$Tumor_BAF[,1]))
  }
  mBAF <- BAF2mBAF(data$data$Tumor_BAF[,1])
  smoB <- round(length(called) / 3300)
  if(smoB%%2 == 0) smoB <- smoB+1
  mBAF.rm <- runmed(mBAF[called], smoB)
  mBAF.diff <- abs(mBAF[called] - mBAF.rm)
  Bfiltered <- mBAF.diff <= quantile(mBAF.diff, BAF.filter)
  if (any(Bfiltered)) data$germline$germlinegenotypes[called][!Bfiltered] <- TRUE
  rm(called, mBAF, smoB, mBAF.rm, mBAF.diff, Bfiltered)
  
  ## Edging BAF
  # bnna <- !is.na(data$data$Tumor_BAF)
  # data$data$Tumor_BAF[bnna,1][data$data$Tumor_BAF[bnna,1] > 1] <- 1
  # data$data$Tumor_BAF[bnna,1][data$data$Tumor_BAF[bnna,1] < 0] <- 0
  
  ## Rorschard plot
  if (plot) {
    png(paste0(samplename, ".Rorschach.png"), width = 980, height = 980)
    EaCoN.Rorschard.plot(data = data)
    dev.off()
  }
  
  ## Computing gaps
  if (!is.null(mingap)) {
    # `%do%` <- foreach::"%do%"
    data$data$chr <- foreach::foreach(k = data$data$ch, .combine = "c") %do% {
      gapz <- which(diff(data$data$SNPpos$pos[k]) >= mingap)
      return(unname(split(k, findInterval(k, k[gapz+1]))))
    }
  }
  
  ## Limiting BAF range
  # data$data$Tumor_BAF[[1]][data$data$Tumor_BAF[[1]] > 1] <- 1
  # data$data$Tumor_BAF[[1]][data$data$Tumor_BAF[[1]] < 0] <- 0
  
  ## ASPCF segmentation
  tmsg("ASPCF segmentation ...")
  aspcf.res <- ASCAT::ascat.aspcf(ASCATobj = data$data, ascat.gg = data$germline, penalty = penalty)
  aspcf.res$germline <- data$germline
  data$data <- aspcf.res
  rm(aspcf.res)
  
  ## Removing big, unused text files
  pcftxt <- list.files(path = ".", pattern = "\\.PCFed\\.txt", full.names = TRUE, recursive = FALSE, ignore.case = TRUE, include.dirs = FALSE)
  if (length(pcftxt) > 0) unlink(pcftxt)
  
  ## Adjusting out-of-range BAF segments
  data$data$Tumor_BAF_segmented[[1]][(data$data$Tumor_BAF_segmented[[1]][,1] < 0),1] <- 0
  
  ## Re-centering
  smo <- round(length(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1]) / 550)
  if(smo%%2 == 0) smo <- smo+1
  if (recenter %in% c("l2r.mainpeak", "l2r.centeredpeak", "l2r.median") | is.numeric(recenter)) {
    tmsg("Recentering ...")
    if (is.numeric(recenter)) {
      tmsg(paste0(" ... using direct value ", recenter, "."))
      shifter <- recenter
    } else if (recenter == "l2r.median") {
      tmsg(" ... using median.")
      shifter <- stats::median(data$data$Tumor_LogR[,1], na.rm = TRUE)
    } else {
      lrrm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1], smo)
      my.den <- stats::density(lrrm, bw = .015)
      if (recenter == "l2r.mainpeak") {
        tmsg(" ... using the main peak.")
        shifter <- my.den$x[which.max(my.den$y)]
      } else if (recenter == "l2r.centeredpeak") {
        tmsg(" ... using the most centered of populated peaks.")
        denf <- data.frame(x = my.den$x, y = my.den$y, sign = c(1, sign(diff(my.den$y))))
        denf$sign2 <- c(diff(denf$sign), 0)
        rrr <- rle(denf$sign2)
        repr <- data.frame(values = rrr$values, start = c(0, (cumsum(rrr$lengths[-c(length(rrr$lengths))]) + 1)), end = cumsum(rrr$lengths), stringsAsFactors = F)
        npk <- which(repr$values == -2)
        
        if (length(npk) == 1) {
          fx <- my.den$x[repr$start[npk[1]]]
        } else {
          if (1 %in% npk) npk <- npk[-1]
          if (nrow(repr) %in% npk) npk <- npk[nrow(repr)]
          parea <- sapply(npk, function(r) { sum(my.den$y[repr$start[r - 1]:repr$end[r + 1]])/sum(my.den$y) })
          parea.rel <- parea/max(parea)
          npk <- npk[parea > 0.1]
          peaklist <- repr$start[npk]
          shifter <- denf$x[peaklist[which.min(sapply(peaklist, function(p) { abs(sum(denf$y[1:(p - 1)]) - sum(denf$y[(p + 1):nrow(denf)]))/sum(denf$y) }))]]
        }
      }
    }
    data$data$Tumor_LogR[, 1] <- data$data$Tumor_LogR[,1] - shifter
    data$data$Tumor_LogR_segmented <- data$data$Tumor_LogR_segmented - shifter
    data$meta$eacon[["recenter-value"]] <- shifter
  } else if (is.null(recenter)) {
    tmsg("No recentering.")
  } else stop(tmsg("Invalid recentering method called !"), call. = FALSE)
  
  ## Winsorization (for aesthetics)
  tmsg("Smoothing L2R (for plots)...")
  cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
  l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
  cndf <- cndf[l2r.nona,]
  cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", method = "mad", k = 5, tau = 1, verbose = FALSE)
  data$data$Tumor_LogR_wins <- data$data$Tumor_LogR
  data$data$Tumor_LogR_wins[l2r.nona,] <- cndf.wins[, 3, drop = FALSE]
  colnames(data$data$Tumor_LogR_wins) <- samplename
  rm(list = c("cndf", "cndf.wins", "l2r.nona"))
  
  ## PELT rescue
  if (!is.null(SER.pen)) {
    tmsg("Rescuing small events ...")
    seg.maxwidth <- 5e+06
    seg.maxn <- 500
    mydf <- data.frame(data$data$SNPpos, l2r = as.vector(data$data$Tumor_LogR[,1]), idx.ori = 1:nrow(data$data$SNPpos))
    mydf$chrs <- as.character(mydf$chrs)
    mydf <- mydf[!is.na(mydf$l2r), ]
    chrends <- cumsum(rle(as.character(data$data$SNPpos$chrs[!is.na(data$data$Tumor_LogR[,1])]))$lengths)
    if (is.character(SER.pen)) {
      seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
    } else if (is.numeric(SER.pen)) {
      if (SER.pen < 1) {
        seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = 'Asymptotic', pen.value = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
      } else {
        SER.pen <- sum(abs(diff(runmed(mydf$l2r[!is.na(mydf$l2r)], smo)))) / SER.pen
        seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = 'Manual', pen.value = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
      }
    } else stop(tmsg("SER.pen should be a character or a numeric !"), call. = FALSE)
    if (is.character(seg.end)) {
      tmsg(" PELT segmentation failed with this combination of SER.pen and segmentLength options !")
      data$meta$eacon[["small.events.rescue.PELT.penalty"]] <- "ERROR"
    } else {
      seg.end <- sort(unique(c(seg.end, chrends)))
      seg.start <- c(1, seg.end[-length(seg.end)] + 1)
      seg.med <- vapply(1:length(seg.end), function(x) { return(median(mydf$l2r[seg.start[x]:seg.end[x]], na.rm = TRUE)) }, 0.1)
      seg.width <- mydf$pos[seg.end] - mydf$pos[seg.start] + 1
      rescued <- which(seg.width < seg.maxwidth)
      tmsg(paste0(" Found ", length(rescued), "."))
      if (length(rescued) > seg.maxn) tmsg("WARNING : Many small events found, profile may be noisy ! Consider using 'smooth.k', or for WES data, strengthen low depth filtering !")
      data$meta$eacon[["PELT-nseg"]] <- length(rescued)
      # `%do%` <- foreach::"%do%"
      foreach::foreach(re = rescued, .combine = "c") %do% {
        interv <- mydf$idx.ori[seg.start[re]]:mydf$idx.ori[seg.end[re]]
        data$data$Tumor_LogR_segmented[interv] <- median(data$data$Tumor_LogR[interv, 1], na.rm = TRUE)
        return()
      }
    }
  }
  
  ## BAF objects
  bafpos <- data$data$SNPpos[rownames(data$data$Tumor_LogR_segmented) %in% rownames(data$data$Tumor_BAF_segmented[[1]]),]
  bafpos$ProbeSet <- rownames(data$data$Tumor_BAF_segmented[[1]])
  bafpos$BAF <- data$data$Tumor_BAF_segmented[[1]][,1]
  bafpos$chrs <- as.character(bafpos$chrs)
  
  bafkend <- vapply(unique(bafpos$chrs), function(x) { max(which(bafpos$chrs == x))}, 1)
  bafbreaks <- cumsum(rle(bafpos$BAF)$lengths)
  baf.seg <- data.frame(end.idx = sort(unique(c(bafkend, bafbreaks))), stringsAsFactors = FALSE)
  baf.seg$start.idx <- c(1, baf.seg$end.idx[-nrow(baf.seg)]+1)
  baf.seg$length.idx <- baf.seg$end.idx - baf.seg$start.idx +1
  baf.seg$start.probeset <- bafpos$ProbeSet[baf.seg$start.idx]
  baf.seg$end.probeset <- bafpos$ProbeSet[baf.seg$end.idx]
  baf.seg$chrA <- bafpos$chrs[baf.seg$start.idx]
  # if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) {
  # baf.seg$Chr <- unlist(cs$chrom2chr[paste0("chr", baf.seg$chrA)])
  # } else baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Start <- bafpos$pos[baf.seg$start.idx]
  baf.seg$End <- bafpos$pos[baf.seg$end.idx]
  baf.seg$Width <- baf.seg$End - baf.seg$Start + 1
  baf.seg$Value <- bafpos$BAF[baf.seg$start.idx]
  
  ## BAF calling
  baf.homocut <- homoCut
  ### Calling
  baf.seg$Status <- "Unbalanced"
  baf.homo <- sort(which(baf.seg$Value <= baf.homocut))
  if (length(baf.homo) > 0) baf.seg$Status[baf.homo] <- "Homo"
  baf.hetero <- sort(which(baf.seg$Value == .5))
  if (length(baf.hetero) > 0) baf.seg$Status[baf.hetero] <- "Hetero"
  
  baf.seg.out <- baf.seg[,c(6,7:12,4:5,1:3)]
  colnames(baf.seg.out) <- c("Chrom", "Chr", "Start", "End", "Width", "BAF.Value", "BAF.Status", "Start.FeatureName", "End.FeatureName", "Start.FeatureIdx", "End.FeatureIdx", "Features")
  
  if (write.data) write.table(baf.seg.out, paste0(samplename, ".SegmentedBAF.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  ## L2R calling
  tmsg("Calling L2R ...")
  l2r.segments <- foreach::foreach(k = 1:length(data$data$ch), .combine = "rbind") %do% {
    segrle <- rle(data$data$Tumor_LogR_segmented[data$data$ch[[k]],1])
    segdf <- data.frame(Probes = segrle$lengths, Value = segrle$values, stringsAsFactors = FALSE)
    segdf$Chr <- unlist(cs$chrom2chr[data$data$chrs[k]])
    segdf$Chrom <- data$data$chrs[[k]]
    probecs <- cumsum(segdf$Probes)
    segdf$End.idx <- as.integer(cumsum(segdf$Probes) + data$data$ch[[k]][1] - 1)
    segdf$Start.idx <-as.integer(c(data$data$ch[[k]][1], segdf$End.idx[-nrow(segdf)] + 1))
    segdf$Start <- as.integer(data$data$SNPpos$pos[segdf$Start.idx])
    segdf$End <- as.integer(data$data$SNPpos$pos[segdf$End.idx])
    segdf <- segdf[, c(3, 4, 7, 8, 2, 1, 6, 5)]
  }
  
  my.mad <- get.mad(data$data$Tumor_LogR[,1])
  
  if (calling.method == "mad") {
    g.cut <- my.mad * nrf
    l.cut <- -g.cut
  }
  l2r.rm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1], smo)
  if (calling.method == "density") {
    lr.den <- stats::density(l2r.rm, bw = .015)
    atf <- lr.den$y > max(lr.den$y)/20
    newx <- lr.den$x[atf]
    newy <- lr.den$y[atf]
    sigdiff <- sign(diff(newy))
    negz <- newx < 0
    l.idx <- length(which(negz)) - rle(rev(sign(diff(newy[negz]))))$lengths[1] +1
    g.idx <- length(which(negz)) + rle(sign(diff(newy[!negz])))$lengths[1]
    g.cut <- newx[g.idx]
    l.cut <- newx[l.idx]
    if (abs(g.cut) > abs(l.cut)) l.cut <- -g.cut else g.cut <- -l.cut
  }
  data$meta$eacon[["L2R-segments-gain-cutoff"]] <- g.cut
  data$meta$eacon[["L2R-segments-loss-cutoff"]] <- l.cut
  
  
  ## Profile metrics
  data$meta$eacon$MAD <- my.mad
  data$meta$eacon$SSAD <- sum(abs(diff(runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]),1], smo))))

  ## Generating CBS
  gain.idx <- which(l2r.segments$Value > g.cut)
  loss.idx <- which(l2r.segments$Value < l.cut)
  normal.idx <- which(l2r.segments$Value >= l.cut & l2r.segments$Value <= g.cut)
  tmsg("Writing CBS files ...")
  data$cbs$nocut <- data.frame(Samplename = samplename, l2r.segments[,c(1, 3, 4, 6, 5)], stringsAsFactors = FALSE)
  colnames(data$cbs$nocut) <- c(samplename, "Chr", "Start", "End", "Probes", "Log2Ratio")
  my.cbs.cut <- data$cbs$nocut
  my.cbs.cut$Log2Ratio[my.cbs.cut$Log2Ratio > l.cut & my.cbs.cut$Log2Ratio < g.cut] <- 0
  data$cbs$cut <- my.cbs.cut
  rm(my.cbs.cut)
  if (write.data) {
    write.table(data$cbs$nocut, paste0(samplename, ".NoCut.cbs"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data$cbs$cut, paste0(samplename, ".Cut.cbs"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  ## REPLOTS
  if (plot) {
    tmsg("Plotting ...")
    l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, loss = loss.idx, normal = normal.idx), cutval = c(l.cut, g.cut))
    seg.col <- list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.loss = "darkred", normal = "black")
    l2r.chr <- unname(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]))
    l2r.value <- data.frame(Chr = l2r.chr,
                            Start = as.integer(data$data$SNPpos$pos),
                            End = as.integer(data$data$SNPpos$pos),
                            Value = data$data$Tumor_LogR_wins[,1],
                            # Value = data$data$Tumor_LogR[,1],
                            stringsAsFactors = FALSE)
    baf.value <- data.frame(Chr = l2r.chr,
                            Start = as.integer(data$data$SNPpos$pos),
                            End = as.integer(data$data$SNPpos$pos),
                            Value = data$data$Tumor_BAF[,1],
                            stringsAsFactors = FALSE)
    
    png(paste0(samplename, ".SEG.ASCAT.png"), width = 1850, height = 980)
    par(mar = c(1, 6, 3, 1), mfrow = c(2, 1))
    EaCoN.l2rplot.geno(l2r = l2r.value,
                       seg = l2r.seg.obj,
                       seg.col = seg.col,
                       seg.type = "block",
                       seg.normal = TRUE,
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " L2R"),
                       ylim = c(-1.5,1.5))
    
    EaCoN.bafplot.geno(baf = baf.value,
                       seg = baf.seg,
                       seg.type = "both",
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " BAF"))
    dev.off()
  }
  
  ## Saving segmentation object
  if (write.data) {
    tmsg("Saving data ...")
    saveRDS(data, paste0(samplename, ".SEG.ASCAT.RDS"), compress = "bzip2")
  }
  
  setwd(oridir)
  tmsg("Done.")
  if (return.data) return(data)
}

## FACETS L2R + BAF segmentation
Segment.FACETS <- function(data = NULL, smooth.k = NULL, BAF.filter = .75, homoCut = .05, penalty = 100, recenter = "l2r.centeredpeak", calling.method = "mad", nrf = .5, SER.pen = 2, out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta2/WES/TCGA-A7-A0CE-01A_vs_10A")
  # data = readRDS("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta2/WES/TCGA-A7-A0CE-01A_vs_10A/TCGA-A7-A0CE-01A_vs_10A_hs37d5_b_processed.RDS")
  # setwd("/home/job/WORKSPACE/EaCoN_tests/OS/FACETSTEST")
  # data <- readRDS("/home/job/WORKSPACE/EaCoN_tests/OS/FACETSTEST/FACETSTEST_OncoScan_CNV_hg19_processed.RDS")
  # smooth.k = 5
  # BAF.filter = .75
  # homoCut = .05
  # penalty = 150
  # recenter = "l2r.centeredpeak"
  # calling.method = "mad"
  # nrf = 0.5
  # SER.pen = 40
  # out.dir = getwd()
  # return.data = FALSE
  # write.data = TRUE
  # plot = TRUE
  # force = FALSE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/plot_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/FACETS_hack.R")

  
  calling.method <- tolower(calling.method)
  
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  if (!dir.exists(out.dir)) stop(tmsg(paste0("Output directory [", out.dir, "] does not exist !")), call. = FALSE)
  if (!(calling.method %in% c("mad", "density"))) stop(tmsg("calling.method should be 'MAD' or 'density' !"), call. = FALSE)
  if (calling.method == "mad" & is.null(nrf)) stop(tmsg("If calling.method is set to 'MAD', nrf is required !"), call. = FALSE)
  if (!is.null(SER.pen)) if (!is.character(SER.pen)) if (SER.pen <= 0) stop(tmsg("SER.pen should be NULL, a character or an integer/float > 0 !"), call. = FALSE)
  
  # valid.types <- c("WES", "WGS", "OncoScan_CNV")
  if(!("additional" %in% names(data$data))) stop(tmsg("FACETS is not available for this data source !"), call. = FALSE)
  
  ## Extract samplename
  samplename <- data$meta$basic$samplename
  tmsg(paste0("Sample : ", samplename))
  
  ## Extract genome version and load corresponding data
  genome <- data$meta$basic$genome
  genome.pkg <- data$meta$basic$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  cs <- chromobjector(BSg.obj)
  
  ## Init graphical parameters
  oridir <- getwd()
  
  ## Create output dir
  if (write.data | plot) {
    odir <- paste0(out.dir, "/FACETS/L2R")
    if (dir.exists(odir)) {
      if (force) {
        unlink(odir, recursive = TRUE, force = FALSE)
        dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
      } else stop(tmsg(paste0("A [", odir, "] dir already exists !")), call. = FALSE)
    } else dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
  } else odir <- out.dir
  
  setwd(odir)
  
  data$meta$eacon <- c(data$meta$eacon, list(
    "segmenter" = "FACETS",
    "BAF.filter" = BAF.filter,
    "BAF.segments.homo.limit" = homoCut,
    "winsorize.k" = if(is.null(smooth.k)) "NA" else smooth.k,
    "FACETS.penalty" = penalty,
    "calling.method" = calling.method,
    "calling.nrf" = if(is.null(nrf)) "NA" else nrf,
    "small.events.rescue.PELT.penalty" = if(is.null(SER.pen)) "NA" else SER.pen
  ))
  
  ## Winsorization
   if(!is.null(smooth.k)) {
    tmsg("Smoothing L2R outliers ...")
    cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
    l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
    cndf <- cndf[l2r.nona,]
    cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", method = "mad", k = smooth.k, tau = 1, verbose = FALSE)
    data$data$Tumor_LogR[l2r.nona,1] <- cndf.wins[, 3, drop = FALSE]
    rm(list = c("cndf", "cndf.wins", "l2r.nona"))
  }
  
  ## BAF filtering
  tmsg("Filtering BAF...")
  if ("Tumor_BAF.unisomy" %in% names(data$data)) {
    called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & !is.na(data$data$Tumor_BAF[,1]) & !data$data$Tumor_BAF.unisomy[,1])
  } else {
    called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & !is.na(data$data$Tumor_BAF[,1]))
  }
  mBAF <- BAF2mBAF(data$data$Tumor_BAF[,1])
  smoB <- round(length(called) / 3300)
  if(smoB%%2 == 0) smoB <- smoB+1
  mBAF.rm <- runmed(mBAF[called], smoB)
  mBAF.diff <- abs(mBAF[called] - mBAF.rm)
  Bfiltered <- mBAF.diff <= quantile(mBAF.diff, BAF.filter)
  if (any(Bfiltered)) data$germline$germlinegenotypes[called][!Bfiltered] <- TRUE
  rm(called, mBAF, smoB, mBAF.rm, mBAF.diff, Bfiltered)
  
  ## Edging BAF
  # bnna <- !is.na(data$data$Tumor_BAF)
  # data$data$Tumor_BAF[bnna,1][data$data$Tumor_BAF[bnna,1] > 1] <- 1
  # data$data$Tumor_BAF[bnna,1][data$data$Tumor_BAF[bnna,1] < 0] <- 0
  
  ## Rorschard plot
  if (plot) {
    png(paste0(samplename, ".Rorschach.png"), width = 980, height = 980)
    EaCoN.Rorschard.plot(data = data)
    dev.off()
  }
  
  ## Computing gaps
  # if (!is.null(mingap)) {
  #   data$data$chr <- foreach(k = data$data$ch, .combine = "c") %do% {
  #     gapz <- which(diff(data$data$SNPpos$pos[k]) >= mingap)
  #     return(unname(split(k, findInterval(k, k[gapz+1]))))
  #   }
  # }
  # 
  ## Limiting BAF range
  # data$data$Tumor_BAF[[1]][data$data$Tumor_BAF[[1]] > 1] <- 1
  # data$data$Tumor_BAF[[1]][data$data$Tumor_BAF[[1]] < 0] <- 0
  
  ## FACETS segmentation
  tmsg("FACETS segmentation ...")
  
  snp.nbhd <- 250
  hetscale <- TRUE
  deltaCN <- 0
  
  ### Data conversion
  # F.dmat <- data.frame(chrom = as.integer(data$data$SNPpos$chrs), maploc = as.numeric(data$data$SNPpos$pos), rCountT = data$data$additional$RD.test, rCountN = data$data$additional$RD.ref, vafT = data$data$Tumor_BAF[,1], vafN = NA, het = abs(as.numeric(data$germline$germlinegenotypes)-1), keep = NA, gcpct = NA, gcbias = NA, cnlr = data$data$Tumor_LogR[,1], valor = data$data$additional$LOR, lorvar = data$data$additional$LORvar)
  F.dmat <- data.frame(chrom = unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]), maploc = as.numeric(data$data$SNPpos$pos), rCountT = data$data$additional$RD.test, rCountN = data$data$additional$RD.ref, vafT = data$data$Tumor_BAF[,1], vafN = NA, het = abs(as.numeric(data$germline$germlinegenotypes)-1), keep = NA, gcpct = NA, gcbias = NA, cnlr = data$data$Tumor_LogR[,1], valor = data$data$additional$LOR, lorvar = data$data$additional$LORvar)
  
  ### Selection of SNP position to smooth their density
  F.dmat$keep <- if ("germline" %in% names(data)) as.numeric(!unname(data$germline$germlinegenotypes[,1])) else F.dmat$keep <- facets:::scanSnp(maploc = F.dmat$maploc, het = F.dmat$het, nbhd = snp.nbhd)
  facets.keep <- F.dmat$keep == 1 & !is.na(F.dmat$vafT) & !is.na(F.dmat$lorvar) & !is.na(F.dmat$cnlr)
  #### Hacking the selection ton include chromosomes edges
  facets.keep[c(vapply(data$data$chr, function(x) {head(x, 1)}, 1L), vapply(data$data$chr, function(x) {tail(x, 1)}, 1L))] <- TRUE
  #### Applying the selectionn
  F.dmat <- F.dmat[facets.keep,]
  
  ## OS hack
  if (any(is.na(F.dmat$lorvar))) F.dmat$lorvar[is.na(F.dmat$lorvar)] <- 0
  if (any(is.na(F.dmat$vafT))) F.dmat$vafT[is.na(F.dmat$vafT)] <- 0.5
  if (any(is.na(F.dmat$cnlr))) F.dmat$cnlr[is.na(F.dmat$cnlr)] <- 0
  
  ### Presegmentation (clustering positions)
  nX <- which(data$data$chrs == data$data$sexchromosomes[1])
  F.tmp <- facets:::segsnps(mat = F.dmat, cval = penalty, hetscale = hetscale, delta = deltaCN)
  rm(F.dmat)
  
  ### Segmentation
  F.proc.out <- procSample.hacked(x = c(gbuild = data$meta$basic$genome, nX = nX, F.tmp), cval = penalty, chromlevels = data$data$chrs)
  rm(F.tmp)
  
  ### Inserting results in data
  ori.idx <- which(facets.keep)
  seg.end.idx <- ori.idx[cumsum(F.proc.out$out$num.mark)]
  seg.start.idx <- c(1, seg.end.idx[-length(seg.end.idx)]+1)
  #### L2R
  # data$data$Tumor_LogR_segmented <- matrix(rep(F.proc.out$out$cnlr.median, diff(c(0,seg.end.idx))) - F.proc.out$dipLogR, ncol = 1)
  data$data$Tumor_LogR_segmented <- matrix(rep(F.proc.out$out$cnlr.median, diff(c(0,seg.end.idx))), ncol = 1)
  rownames(data$data$Tumor_LogR_segmented) <- rownames(data$data$SNPpos)
  #### BAF
  tmbaf <- BAF2mBAF(data$data$Tumor_BAF[,1])
  baf.tbl <- dplyr::as.tbl(data.frame(seg = F.proc.out$jointseg$seg, BAF = tmbaf[facets.keep]))
  baf.tbl <- dplyr::group_by(baf.tbl, seg)
  data$data$Tumor_BAF_segmented <- list(matrix(rep(dplyr::summarize(baf.tbl, baf.med = median(BAF, na.rm = TRUE))$baf.med, diff(c(0,seg.end.idx))), ncol = 1))
  rownames(data$data$Tumor_BAF_segmented[[1]]) <- rownames(data$data$SNPpos)
  #### Hacking BAF near 0.5 with homoCut
  pc50.idx <- data$data$Tumor_BAF_segmented[[1]][,1] > (.5 - homoCut) & data$data$Tumor_BAF_segmented[[1]][,1] < (.5 + homoCut)
  data$data$Tumor_BAF_segmented[[1]][pc50.idx,1] <- .5
  #### Additional
  data$data$additional <- F.proc.out
  
  ## Adjusting out-of-range BAF segments
  data$data$Tumor_BAF_segmented[[1]][(data$data$Tumor_BAF_segmented[[1]][,1] < 0),1] <- 0
  
  ## Re-centering
  smo <- round(length(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1]) / 550)
  if(smo%%2 == 0) smo <- smo+1
  if (recenter %in% c("l2r.mainpeak", "l2r.centeredpeak", "l2r.median") | is.numeric(recenter)) {
    tmsg("Recentering ...")
    if (is.numeric(recenter)) {
      tmsg(paste0(" ... using direct value ", recenter, "."))
      shifter <- recenter
    } else if (recenter == "l2r.median") {
      tmsg(" ... using median.")
      shifter <- stats::median(data$data$Tumor_LogR[,1], na.rm = TRUE)
    } else {
      lrrm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1], smo)
      my.den <- stats::density(lrrm, bw = .015)
      if (recenter == "l2r.mainpeak") {
        tmsg(" ... using the main peak.")
        shifter <- my.den$x[which.max(my.den$y)]
      } else if (recenter == "l2r.centeredpeak") {
        tmsg(" ... using the most centered of populated peaks.")
        denf <- data.frame(x = my.den$x, y = my.den$y, sign = c(1, sign(diff(my.den$y))))
        denf$sign2 <- c(diff(denf$sign), 0)
        rrr <- rle(denf$sign2)
        repr <- data.frame(values = rrr$values, start = c(0, (cumsum(rrr$lengths[-c(length(rrr$lengths))]) + 1)), end = cumsum(rrr$lengths), stringsAsFactors = F)
        npk <- which(repr$values == -2)
        
        if (length(npk) == 1) {
          fx <- my.den$x[repr$start[npk[1]]]
        } else {
          if (1 %in% npk) npk <- npk[-1]
          if (nrow(repr) %in% npk) npk <- npk[nrow(repr)]
          parea <- sapply(npk, function(r) { sum(my.den$y[repr$start[r - 1]:repr$end[r + 1]])/sum(my.den$y) })
          parea.rel <- parea/max(parea)
          npk <- npk[parea > 0.1]
          peaklist <- repr$start[npk]
          shifter <- denf$x[peaklist[which.min(sapply(peaklist, function(p) { abs(sum(denf$y[1:(p - 1)]) - sum(denf$y[(p + 1):nrow(denf)]))/sum(denf$y) }))]]
        }
      }
    }
    data$data$Tumor_LogR[, 1] <- data$data$Tumor_LogR[,1] - shifter
    data$data$Tumor_LogR_segmented <- data$data$Tumor_LogR_segmented - shifter
    data$meta$eacon[["recenter-value"]] <- shifter
  } else if (is.null(recenter)) {
    tmsg("No recentering.")
  } else stop(tmsg("Invalid recentering method called !"), call. = FALSE)
  
  ## Winsorization (for aesthetics)
  tmsg("Smoothing L2R (for plots)...")
  cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
  l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
  cndf <- cndf[l2r.nona,]
  cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", method = "mad", k = 5, tau = 1, verbose = FALSE)
  data$data$Tumor_LogR_wins <- data$data$Tumor_LogR
  data$data$Tumor_LogR_wins[l2r.nona,] <- cndf.wins[, 3, drop = FALSE]
  colnames(data$data$Tumor_LogR_wins) <- samplename
  rm(list = c("cndf", "cndf.wins", "l2r.nona"))
  
  
  ## PELT rescue
  if (!is.null(SER.pen)) {
    tmsg("Rescuing small events ...")
    seg.maxwidth <- 5e+06
    seg.maxn <- 500
    mydf <- data.frame(data$data$SNPpos, l2r = as.vector(data$data$Tumor_LogR[,1]), idx.ori = 1:nrow(data$data$SNPpos))
    mydf$chrs <- as.character(mydf$chrs)
    mydf <- mydf[!is.na(mydf$l2r), ]
    chrends <- cumsum(rle(as.character(data$data$SNPpos$chrs[!is.na(data$data$Tumor_LogR[,1])]))$lengths)
    if (is.character(SER.pen)) {
      seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
    } else if (is.numeric(SER.pen)) {
      if (SER.pen < 1) {
        seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = 'Asymptotic', pen.value = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
      } else {
        SER.pen <- sum(abs(diff(runmed(mydf$l2r[!is.na(mydf$l2r)], smo)))) / SER.pen
        seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = 'Manual', pen.value = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
      }
    } else stop(tmsg("SER.pen should be a character or a numeric !"), call. = FALSE)
    if (is.character(seg.end)) {
      tmsg(" PELT segmentation failed with this combination of SER.pen and segmentLength options !")
      data$meta$eacon[["small.events.rescue.PELT.penalty"]] <- "ERROR"
    } else {
      seg.end <- sort(unique(c(seg.end, chrends)))
      seg.start <- c(1, seg.end[-length(seg.end)] + 1)
      seg.med <- vapply(1:length(seg.end), function(x) { return(median(mydf$l2r[seg.start[x]:seg.end[x]], na.rm = TRUE)) }, 0.1)
      seg.width <- mydf$pos[seg.end] - mydf$pos[seg.start] + 1
      rescued <- which(seg.width < seg.maxwidth)
      tmsg(paste0(" Found ", length(rescued), "."))
      if (length(rescued) > seg.maxn) tmsg("WARNING : Many small events found, profile may be noisy ! Consider using 'smooth.k', or for WES data, strengthen low depth filtering !")
      data$meta$eacon[["PELT-nseg"]] <- length(rescued)
      `%do%` <- foreach::"%do%"
      foreach::foreach(re = rescued, .combine = "c") %do% {
        interv <- mydf$idx.ori[seg.start[re]]:mydf$idx.ori[seg.end[re]]
        data$data$Tumor_LogR_segmented[interv] <- median(data$data$Tumor_LogR[interv, 1], na.rm = TRUE)
        return()
      }
    }
  }
  
  ## BAF objects
  bafpos <- data$data$SNPpos[rownames(data$data$Tumor_LogR_segmented) %in% rownames(data$data$Tumor_BAF_segmented[[1]]),]
  bafpos$ProbeSet <- rownames(data$data$Tumor_BAF_segmented[[1]])
  bafpos$BAF <- data$data$Tumor_BAF_segmented[[1]][,1]
  bafpos$chrs <- as.character(bafpos$chrs)
  
  bafkend <- vapply(unique(bafpos$chrs), function(x) { max(which(bafpos$chrs == x))}, 1)
  bafbreaks <- cumsum(rle(bafpos$BAF)$lengths)
  baf.seg <- data.frame(end.idx = sort(unique(c(bafkend, bafbreaks))), stringsAsFactors = FALSE)
  baf.seg$start.idx <- c(1, baf.seg$end.idx[-nrow(baf.seg)]+1)
  baf.seg$length.idx <- baf.seg$end.idx - baf.seg$start.idx +1
  baf.seg$start.probeset <- bafpos$ProbeSet[baf.seg$start.idx]
  baf.seg$end.probeset <- bafpos$ProbeSet[baf.seg$end.idx]
  baf.seg$chrA <- bafpos$chrs[baf.seg$start.idx]
  # if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) {
  # baf.seg$Chr <- unlist(cs$chrom2chr[paste0("chr", baf.seg$chrA)])
  # } else baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Start <- bafpos$pos[baf.seg$start.idx]
  baf.seg$End <- bafpos$pos[baf.seg$end.idx]
  baf.seg$Width <- baf.seg$End - baf.seg$Start + 1
  baf.seg$Value <- bafpos$BAF[baf.seg$start.idx]
  
  ## BAF calling
  baf.homocut <- homoCut
  ### Calling
  baf.seg$Status <- "Unbalanced"
  baf.homo <- sort(which(baf.seg$Value <= baf.homocut))
  if (length(baf.homo) > 0) baf.seg$Status[baf.homo] <- "Homo"
  baf.hetero <- sort(which(baf.seg$Value == .5))
  if (length(baf.hetero) > 0) baf.seg$Status[baf.hetero] <- "Hetero"
  
  baf.seg.out <- baf.seg[,c(6,7:12,4:5,1:3)]
  colnames(baf.seg.out) <- c("Chrom", "Chr", "Start", "End", "Width", "BAF.Value", "BAF.Status", "Start.FeatureName", "End.FeatureName", "Start.FeatureIdx", "End.FeatureIdx", "Features")
  
  if (write.data) write.table(baf.seg.out, paste0(samplename, ".SegmentedBAF.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  ## L2R calling
  tmsg("Calling L2R ...")
  l2r.segments <- foreach::foreach(k = 1:length(data$data$ch), .combine = "rbind") %do% {
    segrle <- rle(data$data$Tumor_LogR_segmented[data$data$ch[[k]],1])
    segdf <- data.frame(Probes = segrle$lengths, Value = segrle$values, stringsAsFactors = FALSE)
    # segdf$Chr <- k
    segdf$Chr <- unlist(cs$chrom2chr[data$data$chrs[k]])
    segdf$Chrom <- data$data$chrs[[k]]
    probecs <- cumsum(segdf$Probes)
    segdf$End.idx <- cumsum(segdf$Probes) + data$data$ch[[k]][1] - 1
    segdf$Start.idx <- c(data$data$ch[[k]][1], segdf$End.idx[-nrow(segdf)] + 1)
    segdf$Start <- data$data$SNPpos$pos[segdf$Start.idx]
    segdf$End <- data$data$SNPpos$pos[segdf$End.idx]
    segdf <- segdf[, c(3, 4, 7, 8, 2, 1, 6, 5)]
  }
  
  my.mad <- get.mad(data$data$Tumor_LogR[,1])
  
  if (calling.method == "mad") {
    g.cut <- my.mad * nrf
    l.cut <- -g.cut
  }
  l2r.rm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1], smo)
  if (calling.method == "density") {
    lr.den <- stats::density(l2r.rm, bw = .015)
    atf <- lr.den$y > max(lr.den$y)/20
    newx <- lr.den$x[atf]
    newy <- lr.den$y[atf]
    sigdiff <- sign(diff(newy))
    negz <- newx < 0
    l.idx <- length(which(negz)) - rle(rev(sign(diff(newy[negz]))))$lengths[1] +1
    g.idx <- length(which(negz)) + rle(sign(diff(newy[!negz])))$lengths[1]
    g.cut <- newx[g.idx]
    l.cut <- newx[l.idx]
    if (abs(g.cut) > abs(l.cut)) l.cut <- -g.cut else g.cut <- -l.cut
  }
  data$meta$eacon[["L2R-segments-gain-cutoff"]] <- g.cut
  data$meta$eacon[["L2R-segments-loss-cutoff"]] <- l.cut
  
  
  ## Profile metrics
  data$meta$eacon$MAD <- my.mad
  data$meta$eacon$SSAD <- sum(abs(diff(runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]),1], smo))))
  
  ## Generating CBS
  gain.idx <- which(l2r.segments$Value > g.cut)
  loss.idx <- which(l2r.segments$Value < l.cut)
  normal.idx <- which(l2r.segments$Value >= l.cut & l2r.segments$Value <= g.cut)
  tmsg("Writing CBS files ...")
  data$cbs$nocut <- data.frame(Samplename = samplename, l2r.segments[,c(1, 3, 4, 6, 5)], stringsAsFactors = FALSE)
  colnames(data$cbs$nocut) <- c(samplename, "Chr", "Start", "End", "Probes", "Log2Ratio")
  my.cbs.cut <- data$cbs$nocut
  my.cbs.cut$Log2Ratio[my.cbs.cut$Log2Ratio > l.cut & my.cbs.cut$Log2Ratio < g.cut] <- 0
  data$cbs$cut <- my.cbs.cut
  rm(my.cbs.cut)
  if (write.data) {
    write.table(data$cbs$nocut, paste0(samplename, ".NoCut.cbs"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data$cbs$cut, paste0(samplename, ".Cut.cbs"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  ## REPLOTS
  if (plot) {
    tmsg("Plotting ...")
    # l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, loss = loss.idx, normal = normal.idx), cutval = pcut)
    l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, loss = loss.idx, normal = normal.idx), cutval = c(l.cut, g.cut))
    seg.col <- list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.loss = "darkred", normal = "black")
    # l2r.chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", as.character(data$data$SNPpos$chrs))]) else unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)])
    l2r.chr <- unname(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]))
    l2r.value <- data.frame(Chr = l2r.chr,
                            Start = data$data$SNPpos$pos,
                            End = data$data$SNPpos$pos,
                            Value = data$data$Tumor_LogR_wins[,1],
                            # Value = data$data$Tumor_LogR[,1],
                            stringsAsFactors = FALSE)
    # baf.chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", as.character(data$data$SNPpos$chrs))]) else unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)])
    baf.value <- data.frame(Chr = l2r.chr,
                            Start = data$data$SNPpos$pos,
                            End = data$data$SNPpos$pos,
                            Value = data$data$Tumor_BAF[,1],
                            stringsAsFactors = FALSE)
    
    png(paste0(samplename, ".SEG.FACETS.png"), width = 1850, height = 980)
    par(mar = c(1, 6, 3, 1), mfrow = c(2, 1))
    EaCoN.l2rplot.geno(l2r = l2r.value,
                       seg = l2r.seg.obj,
                       seg.col = seg.col,
                       seg.type = "block",
                       seg.normal = TRUE,
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " L2R"),
                       ylim = c(-1.5,1.5))
    
    EaCoN.bafplot.geno(baf = baf.value,
                       seg = baf.seg,
                       seg.type = "both",
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " BAF"))
    dev.off()
  }
  
  ## Saving segmentation object
  if (write.data) {
    tmsg("Saving data ...")
    saveRDS(data, paste0(samplename, ".SEG.FACETS.RDS"), compress = "bzip2")
  }
  
  setwd(oridir)
  tmsg("Done.")
  if (return.data) return(data)
}

## SEQUENZA L2R + BAF segmentation
Segment.SEQUENZA <- function(data = NULL, smooth.k = NULL, BAF.filter = .75, homoCut = .05, penalty = 50, recenter = "l2r.centeredpeak", calling.method = "mad", nrf = .5, SER.pen = 40, out.dir = getwd(), return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) {
  
  # setwd("/home/job/Documents/ROSCOFF/Roscoff_2018/TP_CNV/WES/REDUX/A18R.11.17.18")
  # data = readRDS("A18R.11.17.18_hs37d5_b50_processed.RDS")
  # smooth.k = 5
  # BAF.filter = .75
  # homoCut = .05
  # penalty = 50
  # recenter = "l2r.centeredpeak"
  # calling.method = "mad"
  # nrf = 1
  # SER.pen = 6
  # out.dir = getwd()
  # return.data = FALSE
  # write.data = TRUE
  # plot = TRUE
  # force = FALSE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/plot_functions.R")

  calling.method <- tolower(calling.method)
  
  `%do%` <- foreach::"%do%"
  
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  if (!dir.exists(out.dir)) stop(tmsg(paste0("Output directory [", out.dir, "] does not exist !")), call. = FALSE)
  if (!(calling.method %in% c("mad", "density"))) stop(tmsg("calling.method should be 'MAD' or 'density' !"), call. = FALSE)
  if (calling.method == "mad" & is.null(nrf)) stop(tmsg("If calling.method is set to 'MAD', nrf is required !"), call. = FALSE)
  if (!is.null(SER.pen)) if (!is.character(SER.pen)) if (SER.pen <= 0) stop(tmsg("SER.pen should be NULL, a character or an integer/float > 0 !"), call. = FALSE)
  
  ## Extract samplename
  samplename <- data$meta$basic$samplename
  tmsg(paste0("Sample : ", samplename))
  
  ## Extract genome version and load corresponding data
  genome <- data$meta$basic$genome
  genome.pkg <- data$meta$basic$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  cs <- chromobjector(BSg.obj)
  
  ## Init graphical parameters
  oridir <- getwd()
  
  ## Create output dir
  if (write.data | plot) {
    odir <- paste0(out.dir, "/SEQUENZA/L2R")
    if (dir.exists(odir)) {
      if (force) {
        unlink(odir, recursive = TRUE, force = FALSE)
        dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
      } else stop(tmsg(paste0("A [", odir, "] dir already exists !")), call. = FALSE)
    } else dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
  } else odir <- out.dir
  
  setwd(odir)
  
  data$meta$eacon <- c(data$meta$eacon, list(
    "segmenter" = "SEQUENZA",
    "BAF.filter" = BAF.filter,
    "BAF.segments.homo.limit" = homoCut,
    "winsorize.k" = if(is.null(smooth.k)) "NA" else smooth.k,
    "SEQUENZA.penalty" = penalty,
    "calling.method" = calling.method,
    "calling.nrf" = if(is.null(nrf)) "NA" else nrf,
    "small.events.rescue.PELT.penalty" = if(is.null(SER.pen)) "NA" else SER.pen
  ))
  
  ## Winsorization
  if(!is.null(smooth.k)) {
    tmsg("Smoothing L2R outliers ...")
    cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
    l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
    cndf <- cndf[l2r.nona,]
    cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", method = "mad", k = smooth.k, tau = 1, verbose = FALSE)
    data$data$Tumor_LogR[l2r.nona,1] <- cndf.wins[, 3, drop = FALSE]
    rm(list = c("cndf", "cndf.wins", "l2r.nona"))
  }
  
  ## BAF filtering
  tmsg("Filtering BAF...")
  if ("Tumor_BAF.unisomy" %in% names(data$data)) {
    called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & !is.na(data$data$Tumor_BAF[,1]) & !data$data$Tumor_BAF.unisomy[,1])
  } else {
    called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & !is.na(data$data$Tumor_BAF[,1]))
  }
  mBAF <- BAF2mBAF(data$data$Tumor_BAF[,1])
  smoB <- round(length(called) / 3300)
  if(smoB%%2 == 0) smoB <- smoB+1
  mBAF.rm <- runmed(mBAF[called], smoB)
  mBAF.diff <- abs(mBAF[called] - mBAF.rm)
  Bfiltered <- mBAF.diff <= quantile(mBAF.diff, BAF.filter)
  if (any(Bfiltered)) data$germline$germlinegenotypes[called][!Bfiltered] <- TRUE
  rm(called, mBAF, smoB, mBAF.rm, mBAF.diff, Bfiltered)
  
  ## Edging BAF
  # bnna <- !is.na(data$data$Tumor_BAF)
  # data$data$Tumor_BAF[bnna,1][data$data$Tumor_BAF[bnna,1] > 1] <- 1
  # data$data$Tumor_BAF[bnna,1][data$data$Tumor_BAF[bnna,1] < 0] <- 0
  
  ## Rorschard plot
  if (plot) {
    png(paste0(samplename, ".Rorschach.png"), width = 980, height = 980)
    EaCoN.Rorschard.plot(data = data)
    dev.off()
  }
  
  ## SEQUENZA segmentation
  tmsg("SEQUENZA segmentation ...")
  
  ### Data conversion
  tpos <- as.integer(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]))
  l2r.df <- data.frame(chrs = tpos, pos = as.integer(data$data$SNPpos$pos), s1 = data$data$Tumor_LogR[,1], stringsAsFactors = FALSE)
  tbaf <- data$data$Tumor_BAF[,1]
  tbaf[data$germline$germlinegenotypes] <- 0
  baf.df <- data.frame(chrs = tpos, pos = as.integer(data$data$SNPpos$pos), s1 = tbaf, stringsAsFactors = FALSE)
  colnames(l2r.df)[3] <- colnames(baf.df)[3] <- samplename
  
  ### Blanking bins uncovered by BAF
  isna <- is.na(baf.df[[samplename]])
  baf.df[[samplename]][isna] <- 1
  baf.df[[samplename]][isna][baf.df$pos[isna] %% 2 == 0] <- 0
  
  ### Imputing missing L2R values
  if(any(is.na(l2r.df[,3]))) {
    l2r.tmp <- l2r.df[,3]
    l2r.df[,3] <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  }
  
  # ### Data conversion
  # tpos <- as.integer(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]))
  # l2r.df <- data.frame(chrs = tpos, pos = as.integer(data$data$SNPpos$pos), s1 = data$data$Tumor_LogR[,1], stringsAsFactors = FALSE)
  # 
  # ### Imputing missing L2R values
  # if(any(is.na(l2r.df[,3]))) {
  #   l2r.tmp <- l2r.df[,3]
  #   l2r.df[,3] <- approxfun(seq_along(l2r.tmp), l2r.tmp, rule = 2)(seq_along(l2r.tmp))
  # }
  # 
  # tbaf <- data$data$Tumor_BAF[,1]
  # 
  # k.edges <- as.vector(vapply(unique(l2r.df$chrs), function(k) { c(min(which(l2r.df$chrs == k)), max(which(l2r.df$chrs == k))) }, c(1,1)))
  # p.keep <- data$germline$germlinegenotypes
  # p.keep[is.na(p.keep)] <- TRUE
  # p.keep[is.na(tbaf)] <- TRUE
  # p.keep[k.edges] <- FALSE
  # 
  # baf.df <- data.frame(chrs = tpos, pos = as.integer(data$data$SNPpos$pos), s1 = tbaf, stringsAsFactors = FALSE)
  # colnames(l2r.df)[3] <- colnames(baf.df)[3] <- samplename
  # 
  # # isna <- is.na(baf.df[[samplename]])
  # p.idx <- seq_len(length(data$data$SNPpos$pos))[!p.keep]
  # l2r.df <- l2r.df[!p.keep,]
  # baf.df <- baf.df[!p.keep,]
  # 
  # ### Imputing missing BAF values
  # if(any(is.na(baf.df[,3]))) {
  #   baf.tmp <- baf.df[,3]
  #   baf.df[,3] <- approxfun(seq_along(baf.tmp), baf.tmp, rule = 2)(seq_along(baf.tmp))
  # }
  
  ### Segmentation
  sqz.pcf <- copynumber::aspcf(logR = l2r.df, BAF = baf.df, kmin = 10, gamma = penalty, baf.thres = c(0.05,0.95), arms = rep("p", nrow(l2r.df)), verbose = FALSE)

  ### Hacking n.probes
  # seg.end.idx <- cumsum(sqz.pcf[["n.probes"]])
  # seg.start.idx <- c(1, seg.end.idx[-length(seg.end.idx)]+1)
  # glob.start.idx <- p.idx[seg.start.idx]
  # glob.end.idx <- p.idx[seg.end.idx]
  # new.n.probes <- p.idx[seg.end.idx] - p.idx[seg.start.idx] + 1
  
  ### Inserting results in data
  
  data$data$Tumor_LogR_segmented <- matrix(rep(sqz.pcf$logR.mean, sqz.pcf$n.probes), ncol = 1)
  rownames(data$data$Tumor_LogR_segmented) <- rownames(data$data$SNPpos)
  data$data$Tumor_BAF_segmented <- list(matrix(rep(sqz.pcf$BAF.mean, sqz.pcf$n.probes), ncol = 1))
  rownames(data$data$Tumor_BAF_segmented[[1]]) <- rownames(data$data$SNPpos)
  
  
  #### Additional
  # sqz.pcf$chrom <- data$data$chrs[sqz.pcf$chrom]
  # sqz.pcf$chrom <- unlist(cs$chr2chrom[sqz.pcf$chrom])
  data$data$additional <- sqz.pcf
  
  
  ## Adjusting out-of-range BAF segments
  data$data$Tumor_BAF_segmented[[1]][(data$data$Tumor_BAF_segmented[[1]][,1] < 0),1] <- 0
  
  ## Re-centering
  smo <- round(length(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1]) / 550)
  if(smo%%2 == 0) smo <- smo+1
  if (recenter %in% c("l2r.mainpeak", "l2r.centeredpeak", "l2r.median") | is.numeric(recenter)) {
    tmsg("Recentering ...")
    if (is.numeric(recenter)) {
      tmsg(paste0(" ... using direct value ", recenter, "."))
      shifter <- recenter
    } else if (recenter == "l2r.median") {
      tmsg(" ... using median.")
      shifter <- stats::median(data$data$Tumor_LogR[,1], na.rm = TRUE)
    } else {
      lrrm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1], smo)
      my.den <- stats::density(lrrm, bw = .015)
      if (recenter == "l2r.mainpeak") {
        tmsg(" ... using the main peak.")
        shifter <- my.den$x[which.max(my.den$y)]
      } else if (recenter == "l2r.centeredpeak") {
        tmsg(" ... using the most centered of populated peaks.")
        denf <- data.frame(x = my.den$x, y = my.den$y, sign = c(1, sign(diff(my.den$y))))
        denf$sign2 <- c(diff(denf$sign), 0)
        rrr <- rle(denf$sign2)
        repr <- data.frame(values = rrr$values, start = c(0, (cumsum(rrr$lengths[-c(length(rrr$lengths))]) + 1)), end = cumsum(rrr$lengths), stringsAsFactors = F)
        npk <- which(repr$values == -2)
        
        if (length(npk) == 1) {
          fx <- my.den$x[repr$start[npk[1]]]
        } else {
          if (1 %in% npk) npk <- npk[-1]
          if (nrow(repr) %in% npk) npk <- npk[nrow(repr)]
          parea <- sapply(npk, function(r) { sum(my.den$y[repr$start[r - 1]:repr$end[r + 1]])/sum(my.den$y) })
          parea.rel <- parea/max(parea)
          npk <- npk[parea > 0.1]
          peaklist <- repr$start[npk]
          shifter <- denf$x[peaklist[which.min(sapply(peaklist, function(p) { abs(sum(denf$y[1:(p - 1)]) - sum(denf$y[(p + 1):nrow(denf)]))/sum(denf$y) }))]]
        }
      }
    }
    data$data$Tumor_LogR[, 1] <- data$data$Tumor_LogR[,1] - shifter
    data$data$Tumor_LogR_segmented <- data$data$Tumor_LogR_segmented - shifter
    data$meta$eacon[["recenter-value"]] <- shifter
  } else if (is.null(recenter)) {
    tmsg("No recentering.")
  } else stop(tmsg("Invalid recentering method called !"), call. = FALSE)
  
  ## Winsorization (for aesthetics)
  tmsg("Smoothing L2R (for plots)...")
  cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
  l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
  cndf <- cndf[l2r.nona,]
  cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", method = "mad", k = 5, tau = 1, verbose = FALSE)
  data$data$Tumor_LogR_wins <- data$data$Tumor_LogR
  data$data$Tumor_LogR_wins[l2r.nona,] <- cndf.wins[, 3, drop = FALSE]
  colnames(data$data$Tumor_LogR_wins) <- samplename
  rm(list = c("cndf", "cndf.wins", "l2r.nona"))
  
  ## PELT rescue
  if (!is.null(SER.pen)) {
    tmsg("Rescuing small events ...")
    seg.maxwidth <- 5e+06
    seg.maxn <- 500
    mydf <- data.frame(data$data$SNPpos, l2r = as.vector(data$data$Tumor_LogR[,1]), idx.ori = 1:nrow(data$data$SNPpos))
    mydf$chrs <- as.character(mydf$chrs)
    mydf <- mydf[!is.na(mydf$l2r), ]
    chrends <- cumsum(rle(as.character(data$data$SNPpos$chrs[!is.na(data$data$Tumor_LogR[,1])]))$lengths)
    if (is.character(SER.pen)) {
      seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
    } else if (is.numeric(SER.pen)) {
      if (SER.pen < 1) {
        seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = 'Asymptotic', pen.value = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
      } else {
        SER.pen <- sum(abs(diff(runmed(mydf$l2r[!is.na(mydf$l2r)], smo)))) / SER.pen
        seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, penalty = 'Manual', pen.value = SER.pen, method = "PELT", param.estimates = FALSE, minseglen = 5)@cpts))
      }
    } else stop(tmsg("SER.pen should be a character or a numeric !"), call. = FALSE)
    if (is.character(seg.end)) {
      tmsg(" PELT segmentation failed with this combination of SER.pen and segmentLength options !")
      data$meta$eacon[["small.events.rescue.PELT.penalty"]] <- "ERROR"
    } else {
      seg.end <- sort(unique(c(seg.end, chrends)))
      seg.start <- c(1, seg.end[-length(seg.end)] + 1)
      seg.med <- vapply(1:length(seg.end), function(x) { return(median(mydf$l2r[seg.start[x]:seg.end[x]], na.rm = TRUE)) }, 0.1)
      seg.width <- mydf$pos[seg.end] - mydf$pos[seg.start] + 1
      rescued <- which(seg.width < seg.maxwidth)
      tmsg(paste0(" Found ", length(rescued), "."))
      if (length(rescued) > seg.maxn) tmsg("WARNING : Many small events found, profile may be noisy ! Consider using 'smooth.k', or for WES data, strengthen low depth filtering !")
      data$meta$eacon[["PELT-nseg"]] <- length(rescued)
      foreach::foreach(re = rescued, .combine = "c") %do% {
        interv <- mydf$idx.ori[seg.start[re]]:mydf$idx.ori[seg.end[re]]
        data$data$Tumor_LogR_segmented[interv] <- median(data$data$Tumor_LogR[interv, 1], na.rm = TRUE)
        return()
      }
    }
  }
  
  ## BAF objects
  bafpos <- data$data$SNPpos[rownames(data$data$Tumor_LogR_segmented) %in% rownames(data$data$Tumor_BAF_segmented[[1]]),]
  bafpos$ProbeSet <- rownames(data$data$Tumor_BAF_segmented[[1]])
  bafpos$BAF <- data$data$Tumor_BAF_segmented[[1]][,1]
  bafpos$chrs <- as.character(bafpos$chrs)
  
  bafkend <- vapply(unique(bafpos$chrs), function(x) { max(which(bafpos$chrs == x))}, 1)
  bafbreaks <- cumsum(rle(bafpos$BAF)$lengths)
  baf.seg <- data.frame(end.idx = sort(unique(c(bafkend, bafbreaks))), stringsAsFactors = FALSE)
  baf.seg$start.idx <- c(1, baf.seg$end.idx[-nrow(baf.seg)]+1)
  baf.seg$length.idx <- baf.seg$end.idx - baf.seg$start.idx +1
  baf.seg$start.probeset <- bafpos$ProbeSet[baf.seg$start.idx]
  baf.seg$end.probeset <- bafpos$ProbeSet[baf.seg$end.idx]
  baf.seg$chrA <- bafpos$chrs[baf.seg$start.idx]
  # if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) {
  # baf.seg$Chr <- unlist(cs$chrom2chr[paste0("chr", baf.seg$chrA)])
  # } else baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Start <- bafpos$pos[baf.seg$start.idx]
  baf.seg$End <- bafpos$pos[baf.seg$end.idx]
  baf.seg$Width <- baf.seg$End - baf.seg$Start + 1
  baf.seg$Value <- bafpos$BAF[baf.seg$start.idx]
  
  ## BAF calling
  baf.homocut <- homoCut
  ### Calling
  baf.seg$Status <- "Unbalanced"
  baf.homo <- sort(which(baf.seg$Value <= baf.homocut))
  if (length(baf.homo) > 0) baf.seg$Status[baf.homo] <- "Homo"
  baf.hetero <- sort(which(baf.seg$Value == .5))
  if (length(baf.hetero) > 0) baf.seg$Status[baf.hetero] <- "Hetero"
  
  baf.seg.out <- baf.seg[,c(6,7:12,4:5,1:3)]
  colnames(baf.seg.out) <- c("Chrom", "Chr", "Start", "End", "Width", "BAF.Value", "BAF.Status", "Start.FeatureName", "End.FeatureName", "Start.FeatureIdx", "End.FeatureIdx", "Features")
  
  if (write.data) write.table(baf.seg.out, paste0(samplename, ".SegmentedBAF.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  ## L2R calling
  tmsg("Calling L2R ...")
  l2r.segments <- foreach::foreach(k = 1:length(data$data$ch), .combine = "rbind") %do% {
    segrle <- rle(data$data$Tumor_LogR_segmented[data$data$ch[[k]],1])
    segdf <- data.frame(Probes = segrle$lengths, Value = segrle$values, stringsAsFactors = FALSE)
    # segdf$Chr <- k
    segdf$Chr <- unlist(cs$chrom2chr[data$data$chrs[k]])
    segdf$Chrom <- data$data$chrs[[k]]
    probecs <- cumsum(segdf$Probes)
    segdf$End.idx <- cumsum(segdf$Probes) + data$data$ch[[k]][1] - 1
    segdf$Start.idx <- c(data$data$ch[[k]][1], segdf$End.idx[-nrow(segdf)] + 1)
    segdf$Start <- data$data$SNPpos$pos[segdf$Start.idx]
    segdf$End <- data$data$SNPpos$pos[segdf$End.idx]
    segdf <- segdf[, c(3, 4, 7, 8, 2, 1, 6, 5)]
  }
  
  my.mad <- get.mad(data$data$Tumor_LogR[,1])
  
  if (calling.method == "mad") {
    g.cut <- my.mad * nrf
    l.cut <- -g.cut
  }
  l2r.rm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]), 1], smo)
  if (calling.method == "density") {
    lr.den <- stats::density(l2r.rm, bw = .015)
    atf <- lr.den$y > max(lr.den$y)/20
    newx <- lr.den$x[atf]
    newy <- lr.den$y[atf]
    sigdiff <- sign(diff(newy))
    negz <- newx < 0
    l.idx <- length(which(negz)) - rle(rev(sign(diff(newy[negz]))))$lengths[1] +1
    g.idx <- length(which(negz)) + rle(sign(diff(newy[!negz])))$lengths[1]
    g.cut <- newx[g.idx]
    l.cut <- newx[l.idx]
    if (abs(g.cut) > abs(l.cut)) l.cut <- -g.cut else g.cut <- -l.cut
  }
  data$meta$eacon[["L2R-segments-gain-cutoff"]] <- g.cut
  data$meta$eacon[["L2R-segments-loss-cutoff"]] <- l.cut
  
  
  ## Profile metrics
  data$meta$eacon$MAD <- my.mad
  data$meta$eacon$SSAD <- sum(abs(diff(runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[,1]),1], smo))))
  
  ## Generating CBS
  gain.idx <- which(l2r.segments$Value > g.cut)
  loss.idx <- which(l2r.segments$Value < l.cut)
  normal.idx <- which(l2r.segments$Value >= l.cut & l2r.segments$Value <= g.cut)
  tmsg("Writing CBS files ...")
  data$cbs$nocut <- data.frame(Samplename = samplename, l2r.segments[,c(1, 3, 4, 6, 5)], stringsAsFactors = FALSE)
  colnames(data$cbs$nocut) <- c(samplename, "Chr", "Start", "End", "Probes", "Log2Ratio")
  my.cbs.cut <- data$cbs$nocut
  my.cbs.cut$Log2Ratio[my.cbs.cut$Log2Ratio > l.cut & my.cbs.cut$Log2Ratio < g.cut] <- 0
  data$cbs$cut <- my.cbs.cut
  rm(my.cbs.cut)
  if (write.data) {
    write.table(data$cbs$nocut, paste0(samplename, ".NoCut.cbs"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data$cbs$cut, paste0(samplename, ".Cut.cbs"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  ## REPLOTS
  if (plot) {
    tmsg("Plotting ...")
    # l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, loss = loss.idx, normal = normal.idx), cutval = pcut)
    l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, loss = loss.idx, normal = normal.idx), cutval = c(l.cut, g.cut))
    seg.col <- list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.loss = "darkred", normal = "black")
    # l2r.chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", as.character(data$data$SNPpos$chrs))]) else unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)])
    l2r.chr <- unname(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]))
    l2r.value <- data.frame(Chr = l2r.chr,
                            Start = data$data$SNPpos$pos,
                            End = data$data$SNPpos$pos,
                            Value = data$data$Tumor_LogR_wins[,1],
                            # Value = data$data$Tumor_LogR[,1],
                            stringsAsFactors = FALSE)
    # baf.chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", as.character(data$data$SNPpos$chrs))]) else unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)])
    baf.value <- data.frame(Chr = l2r.chr,
                            Start = data$data$SNPpos$pos,
                            End = data$data$SNPpos$pos,
                            Value = data$data$Tumor_BAF[,1],
                            stringsAsFactors = FALSE)
    
    png(paste0(samplename, ".SEG.SEQUENZA.png"), width = 1850, height = 980)
    par(mar = c(1, 6, 3, 1), mfrow = c(2, 1))
    EaCoN.l2rplot.geno(l2r = l2r.value,
                       seg = l2r.seg.obj,
                       seg.col = seg.col,
                       seg.type = "block",
                       seg.normal = TRUE,
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " L2R"),
                       ylim = c(-1.5,1.5))
    
    EaCoN.bafplot.geno(baf = baf.value,
                       seg = baf.seg,
                       seg.type = "both",
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " BAF"))
    dev.off()
  }
  
  ## Saving segmentation object
  if (write.data) {
    tmsg("Saving data ...")
    saveRDS(data, paste0(samplename, ".SEG.SEQUENZA.RDS"), compress = "bzip2")
  }
  
  setwd(oridir)
  tmsg("Done.")
  if (return.data) return(data)
}


## Run segmentation, from a file
Segment.ff <- function(RDS.file = NULL, segmenter = "ASCAT", ...) {
  if (is.null(RDS.file)) stop(tmsg("A RDS file is needed !"), call. = FALSE)
  valid.segmenters <- c("ASCAT", "FACETS", "SEQUENZA")
  if (!(toupper(segmenter) %in% valid.segmenters)) stop(tmsg(paste0("Segmenter should be one of : ", paste0(valid.segmenters, collapse = ", "))), call. = FALSE)
  if (!file.exists(RDS.file)) stop(tmsg(paste0("Could not find RDS file ", RDS.file, " !")), call. = FALSE)

  ## Data loading
  tmsg(paste0("Loading data from ", RDS.file, " ..."))
  my.data <- readRDS(RDS.file)
  # Segment.ASCAT(data = my.data, out.dir = dirname(RDS.file), ...)
  do.call(paste0("Segment.", toupper(segmenter)), list(data = my.data, out.dir = dirname(RDS.file), ...))
}

## Run Segment.ff() in batch mode
Segment.ff.Batch <- function (RDS.files = list.files(path = getwd(), pattern = "_processed.RDS$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE), segmenter = "ASCAT", nthread = 1, cluster.type = "PSOCK", ...) {
  if (length(RDS.files) == 0) stop("A list of RDS files is required !", call. = FALSE)
  message("Running EaCoN.Segment.ff() in batch mode ...")
  message(paste0(" Found ", length(RDS.files), " files to process."))
  current.bitmapType <- getOption("bitmapType")
  if (length(RDS.files) < nthread) nthread <- length(RDS.files)
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  eacon.batchres <- foreach::foreach(r = seq_along(RDS.files), .inorder = TRUE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    Segment.ff(RDS.file = RDS.files[r], segmenter = segmenter, ...)
  }
  parallel::stopCluster(cl)
}

## ASCAT Total and Allele-Specific Copy Number
ASCN.ASCAT <- function(data = NULL, gammaRange = c(.35,.95), nsubthread = 1, cluster.type = "PSOCK", out.dir = getwd(), force = FALSE, ...) {
  
  # setwd("/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta1/WES/TCGA-A7-A0CE-01A_vs_10A/ASCAT/L2R")
  # data <- readRDS("TCGA-A7-A0CE-01A_vs_10A.EaCoN.ASPCF.RDS")
  # gammaRange <- c(0.35,.95)
  # nsubthread <- 1
  # cluster.type = "PSOCK"
  # out.dir <- "/mnt/data_cigogne/job/PUBLI_EaCoN/TCGA/ANALYSES/EaCoN_0.3.0_beta1/WES/TCGA-A7-A0CE-01A_vs_10A"
  # force <- TRUE
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # require(foreach)

  
  ## CHECKS
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  odir <- paste0(out.dir, "/ASCAT/ASCN")
  if (any(is.null(c(data$data$Tumor_LogR_segmented, data$data$Tumor_BAF_segmented)))) stop(tmsg("No segmentation data found in the provided RDS file !"), call. = FALSE)
  if (dir.exists(odir)) {
    if (force) {
      unlink(odir, recursive = TRUE, force = FALSE)
    } else stop(tmsg(paste0("Output directory [", out.dir, "] already exists !")), call. = FALSE)
  }
  
  samplename <- data$meta$basic$samplename
  ## Extract genome version and load corresponding data
  genome <- data$meta$basic$genome
  genome.pkg <- data$meta$basic$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  cs <- chromobjector(BSg.obj)
  
  tmsg("ASCN modeling (using ASCAT) ...")
  gammavec <- if(length(gammaRange) > 1 ) seq(gammaRange[1], gammaRange[2], 0.05) else gammaRange
  oridirx <- getwd()
  
  cls <- parallel::makeCluster(spec = nsubthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cls)
  gamma <- 0
  `%dopar%` <- foreach::"%dopar%"
  fit.val <- as.data.frame(foreach::foreach(gamma = gammavec, .combine = "rbind", .inorder = TRUE) %dopar% {
    tmsg(paste0(" gamma = ", gamma))
    odirg <- paste0(odir, "/gamma", sprintf("%.2f", gamma))
    dir.create(path = odirg, recursive = TRUE, showWarnings = FALSE)
    setwd(odirg)
    my.ascat.seg.ascn <- suppressWarnings(ASCAT::ascat.runAscat(ASCATobj = data$data, gamma = gamma, ...))
    
    if (is.null(my.ascat.seg.ascn$nA)) {
      tmsg("  ASCAT could not find an optimal ploidy / cellularity from the data.")
      setwd(oridirx)
      unlink(odirg, recursive = TRUE)
      return(rep(NA, 8))
    }
    else {
      
      unlink(paste0(samplename, ".ASCATprofile.png"))
      unlink(paste0(samplename, ".sunrise.png"))
      
      ## Correcting bugged ploidy in ASCAT
      # my.tcn <- my.ascat.seg.ascn$segments$nMajor + my.ascat.seg.ascn$segments$nMinor
      # median(my.tcn)
      # wt.tcn <- my.tcn * (my.ascat.seg.ascn$segments$endpos - my.ascat.seg.ascn$segments$startpos + 1)
      # sum(wt.tcn) / sum(my.ascat.seg.ascn$segments$endpos - my.ascat.seg.ascn$segments$startpos + 1)
      
      tcn.tbl.ung <- dplyr::as.tbl(cbind(my.ascat.seg.ascn$segments, nTotal = my.ascat.seg.ascn$segments$nMajor + my.ascat.seg.ascn$segments$nMinor, width = my.ascat.seg.ascn$segments$endpos - my.ascat.seg.ascn$segments$startpos + 1))
      tcn.tbl <- dplyr::group_by(tcn.tbl.ung, nTotal)
      tcn.tbl.prop <- dplyr::summarise(tcn.tbl, tot_width = sum(width))
      ascat.ploidy <- my.ascat.seg.ascn$ploidy
      median.ploidy <- limma::weighted.median(tcn.tbl.prop$nTotal, tcn.tbl.prop$tot_width)
      baseline.ploidy <- tcn.tbl.prop$nTotal[which.max(tcn.tbl.prop$tot_width)]
      if (baseline.ploidy == 0) baseline.ploidy <- tcn.tbl.prop$nTotal[order(tcn.tbl.prop$tot_width, decreasing = TRUE)[2]]
      weighted.ploidy <- sum(tcn.tbl.prop$nTotal * (tcn.tbl.prop$tot_width / sum(tcn.tbl.prop$tot_width)))
      
      my.ascat.seg.ascn$ploidy <- list(ascat = as.numeric(my.ascat.seg.ascn$ploidy), median = median.ploidy, most.width = baseline.ploidy, width.weighted = weighted.ploidy)
      
      ## Saving ASCN object
      saveRDS(my.ascat.seg.ascn, paste0(samplename, ".ASCN.ASCAT.RDS"), compress = "bzip2")

      ## Generating TCN-CBS
      outfile <- paste0(samplename, ".gamma", gamma, ".cn")
      outdf <- my.ascat.seg.ascn$segments
      outdf$chrom <- outdf$chr
      outdf$chr <- unlist(cs$chrom2chr[outdf$chrom])
      # outdf$chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", outdf$chr)]) else unlist(cs$chrom2chr[outdf$chr])
      outdf$Width <- outdf$endpos - outdf$startpos + 1
      outdf$TCN <- outdf$nMajor + outdf$nMinor
      # outdf$Ploidy.ASCAT <- my.ascat.seg.ascn$ploidy$ascat
      # outdf$Ploidy.Median <- my.ascat.seg.ascn$ploidy$median
      # outdf$Ploidy.MostWidth <- my.ascat.seg.ascn$ploidy$most.width
      # outdf$Ploidy.WidthWeighted <- my.ascat.seg.ascn$ploidy$width.weighted
      # outdf$GoF <- my.ascat.seg.ascn$goodnessOfFit
      # outdf <- outdf[,c(1:4,7,8,5,6,9:12)]
      outdf <- outdf[,c(1,2,7,3,4,8,9,5,6)]
      colnames(outdf)[1:5] <- c(samplename, "Chr", "Chrom", "Start", "End")
      write.table.fast(x = outdf, file = outfile)
      
      ## Generating cellularity + ploidy + metrics file
      outfile <- paste0(samplename, ".gamma", gamma, "_model.txt")
      modeldf <- data.frame(key = c("Sample", "Gamma", "Goodness.of.Fit", "Psi", "Ploidy.ASCAT", "Ploidy.Median", "Ploidy.Most.Width", "Ploidy.Width.weighted", "Cellularity"), value = c(samplename, gamma, unname(my.ascat.seg.ascn$goodnessOfFit), unname(my.ascat.seg.ascn$psi), my.ascat.seg.ascn$ploidy$ascat, my.ascat.seg.ascn$ploidy$median, my.ascat.seg.ascn$ploidy$most.width, my.ascat.seg.ascn$ploidy$width.weighted, unname(my.ascat.seg.ascn$aberrantcellfraction)), stringsAsFactors = FALSE)
      write.table.fast(x = modeldf, file = outfile, header = FALSE)
      
      ## Reploting
      ylim <- 6
      ### TCN
      segments.genostart <- cs$chromosomes$chr.length.toadd[outdf$Chr] + my.ascat.seg.ascn$segments$startpos
      segments.genoend <- cs$chromosomes$chr.length.toadd[outdf$Chr] + my.ascat.seg.ascn$segments$endpos
      ink <- cs$chromosomes$chrN %in% outdf$Chr
      png(paste0(samplename, ".ASCN.ASCAT.png"), width = 1850, height = 980)
      par(mar = c(1, 4, 5, 1), mfrow = c(2, 1))
      # plot(0, 0, type = "n", xlim = c(0,max(segments.genoend)), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " Allele-Specific Copy Number (Gamma = ", gamma, ")\n Cellularity = ", my.ascat.seg.ascn$aberrantcellfraction, " // Ploidy = (A:", round(my.ascat.seg.ascn$ploidy$ascat, digits = 2), ", M:", round(my.ascat.seg.ascn$ploidy$median, digits = 2), ", MW:", round(my.ascat.seg.ascn$ploidy$most.width, digits = 2), ", WW:", round(my.ascat.seg.ascn$ploidy$width.weighted, digits = 2), ") // Goodness of fit = ", round(my.ascat.seg.ascn$goodnessOfFit, digits = 2), "% // Psi = ", my.ascat.seg.ascn$psi, " // nSeg = ", nrow(my.ascat.seg.ascn$segments), " // Predicted gender = ", data$data$gender), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
      plot(0, 0, type = "n", xlim = c(0,cs$genome.length), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " Allele-Specific Copy Number (Gamma = ", gamma, ")\n Cellularity = ", my.ascat.seg.ascn$aberrantcellfraction, " // Ploidy = (A:", round(my.ascat.seg.ascn$ploidy$ascat, digits = 2), ", M:", round(my.ascat.seg.ascn$ploidy$median, digits = 2), ", MW:", round(my.ascat.seg.ascn$ploidy$most.width, digits = 2), ", WW:", round(my.ascat.seg.ascn$ploidy$width.weighted, digits = 2), ") // Goodness of fit = ", round(my.ascat.seg.ascn$goodnessOfFit, digits = 2), "% // Psi = ", my.ascat.seg.ascn$psi, " // nSeg = ", nrow(my.ascat.seg.ascn$segments), " // Predicted gender = ", data$data$gender), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
      abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
      segments(segments.genostart, my.ascat.seg.ascn$segments$nMajor + 0.05, segments.genoend, my.ascat.seg.ascn$segments$nMajor + 0.05, pch = ".", col = 2, lwd = 6)
      segments(segments.genostart, my.ascat.seg.ascn$segments$nMinor - 0.05, segments.genoend, my.ascat.seg.ascn$segments$nMinor - 0.05, pch = ".", col = 3, lwd = 6)
      up.nMajor <- which(my.ascat.seg.ascn$segments$nMajor > ylim)
      up.nMinor <- which(my.ascat.seg.ascn$segments$nMinor > ylim)
      dn.nMajor <- which(my.ascat.seg.ascn$segments$nMajor == 0)
      dn.nMinor <- which(my.ascat.seg.ascn$segments$nMinor == 0)
      if (length(up.nMajor) > 0)
        segments(segments.genostart[up.nMajor],
                 (ylim + 0.2) + 0.05, segments.genoend[up.nMajor],
                 (ylim + 0.2) + 0.05, pch = ".", col = "orange",
                 lwd = 6)
      if (length(up.nMinor) > 0)
        segments(segments.genostart[up.nMinor],
                 (ylim + 0.2) - 0.05, segments.genoend[up.nMinor],
                 (ylim + 0.2) - 0.05, pch = ".", col = "lightgreen",
                 lwd = 6)
      if (length(dn.nMajor) > 0)
        segments(segments.genostart[dn.nMajor],
                 0.05, segments.genoend[dn.nMajor], 0.05,
                 pch = ".", col = "darkred", lwd = 6)
      if (length(dn.nMinor) > 0)
        segments(segments.genostart[dn.nMinor],
                 -0.05, segments.genoend[dn.nMinor], -0.05,
                 pch = ".", col = "darkgreen", lwd = 6)
      abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
      abline(h = 0:ylim, col = "grey75", lty = 3)
      
      try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
      
      # graphics::plot(0, 0, type = "n", xlim = c(0, max(segments.genoend)), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
      graphics::plot(0, 0, type = "n", xlim = c(0, cs$genome.length), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
      abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
      abline(h = my.ascat.seg.ascn$ploidy$median, col = "red", lty = 2)
      abline(h = my.ascat.seg.ascn$ploidy$most.width, col = "cyan", lty = 2)
      nTotal <- my.ascat.seg.ascn$segments$nMajor + my.ascat.seg.ascn$segments$nMinor
      up.nTotal <- which(nTotal > ylim)
      dn.nTotal <- which(nTotal == 0)
      segments(segments.genostart, nTotal, segments.genoend, nTotal, pch = ".", col = 4, lwd = 6)
      if (length(up.nTotal) > 0) segments(segments.genostart[up.nTotal], (ylim + 0.2) + 0.05, segments.genoend[up.nTotal], (ylim + 0.2) + 0.05, pch = ".", col = "cyan", lwd = 6)
      if (length(dn.nTotal) > 0) segments(segments.genostart[dn.nTotal], 0, segments.genoend[dn.nTotal], 0, pch = ".", col = "midnightblue", lwd = 6)
      abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
      abline(h = 0:ylim, col = "grey75", lty = 3)
      try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
      dev.off()
      
      ### RAW CN
      segments.posN <- unlist(cs$chrom2chr[my.ascat.seg.ascn$segments_raw$chr])
      # segments.posN[segments.posN == "X"] <- 23
      # segments.posN[segments.posN == "Y"] <- 24
      # segments.posN <- as.numeric(segments.posN)
      segments.genostart <- cs$chromosomes$chr.length.toadd[segments.posN] + my.ascat.seg.ascn$segments_raw$startpos
      segments.genoend <- cs$chromosomes$chr.length.toadd[segments.posN] + my.ascat.seg.ascn$segments_raw$endpos
      png(paste0(samplename, ".rawprofile.png"), width = 1850, height = 980)
      par(mar = c(1, 4, 5, 1), mfrow = c(2, 1))
      # plot(0, 0, type = "n", xlim = c(0,max(segments.genoend)), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " RAW Allele-Specific Copy Number (Gamma = ", gamma, ")"), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
      plot(0, 0, type = "n", xlim = c(0,cs$genome.length), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " RAW Allele-Specific Copy Number (Gamma = ", gamma, ")"), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
      abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
      segments(segments.genostart, my.ascat.seg.ascn$segments_raw$nAraw + 0.05, segments.genoend, my.ascat.seg.ascn$segments_raw$nAraw + 0.05, pch = ".", col = 2, lwd = 6)
      segments(segments.genostart, my.ascat.seg.ascn$segments_raw$nBraw - 0.05, segments.genoend, my.ascat.seg.ascn$segments_raw$nBraw - 0.05, pch = ".", col = 3, lwd = 6)
      up.nMajor <- which(my.ascat.seg.ascn$segments_raw$nAraw > ylim)
      up.nMinor <- which(my.ascat.seg.ascn$segments_raw$nBraw > ylim)
      dn.nMajor <- which(my.ascat.seg.ascn$segments_raw$nAraw == 0)
      dn.nMinor <- which(my.ascat.seg.ascn$segments_raw$nBraw == 0)
      if (length(up.nMajor) > 0)
        segments(segments.genostart[up.nMajor],
                 (ylim + 0.2) + 0.05, segments.genoend[up.nMajor],
                 (ylim + 0.2) + 0.05, pch = ".", col = "orange",
                 lwd = 6)
      if (length(up.nMinor) > 0)
        segments(segments.genostart[up.nMinor],
                 (ylim + 0.2) - 0.05, segments.genoend[up.nMinor],
                 (ylim + 0.2) - 0.05, pch = ".", col = "lightgreen",
                 lwd = 6)
      if (length(dn.nMajor) > 0)
        segments(segments.genostart[dn.nMajor],
                 0.05, segments.genoend[dn.nMajor], 0.05,
                 pch = ".", col = "darkred", lwd = 6)
      if (length(dn.nMinor) > 0)
        segments(segments.genostart[dn.nMinor],
                 -0.05, segments.genoend[dn.nMinor], -0.05,
                 pch = ".", col = "darkgreen", lwd = 6)
      # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
      abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 3, lwd = 2)
      abline(h = 0:ylim, col = "grey75", lty = 3)
      try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
      # graphics::plot(0, 0, type = "n", xlim = c(0, max(segments.genoend)), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " RAW Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
      graphics::plot(0, 0, type = "n", xlim = c(0, cs$genome.length), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " RAW Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
      abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
      abline(h = my.ascat.seg.ascn$ploidy$median, col = 2, lty = 2)
      abline(h = my.ascat.seg.ascn$ploidy$most.width, col = "purple", lty = 2)
      # nTotal <- my.ascat.seg.ascn$segments_raw$nMajor + my.ascat.seg.ascn$segments_raw$nMinor
      nTotal <- my.ascat.seg.ascn$segments_raw$nAraw + my.ascat.seg.ascn$segments_raw$nBraw
      up.nTotal <- which(nTotal > ylim)
      dn.nTotal <- which(nTotal == 0)
      segments(segments.genostart, nTotal, segments.genoend, nTotal, pch = ".", col = 4, lwd = 6)
      if (length(up.nTotal) > 0) segments(segments.genostart[up.nTotal], (ylim + 0.2) + 0.05, segments.genoend[up.nTotal], (ylim + 0.2) + 0.05, pch = ".", col = "cyan", lwd = 6)
      if (length(dn.nTotal) > 0) segments(segments.genostart[dn.nTotal], 0, segments.genoend[dn.nTotal], 0, pch = ".", col = "midnightblue", lwd = 6)
      # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
      abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 3, lwd = 2)
      abline(h = 0:ylim, col = "grey75", lty = 3)
      try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
      dev.off()
      
      ### TCNvsL2R
      cnpTotal <- my.ascat.seg.ascn$nA + my.ascat.seg.ascn$nB
      my.xlim <- c(min(cnpTotal, na.rm = TRUE) - 0.5, max(cnpTotal, na.rm = TRUE) + 0.5)
      png(paste0(samplename, ".TCNvsL2R.png"), width = 980, height = 980)
      par(mar = c(4, 4, 5, 1))
      graphics::plot(cnpTotal, data$data$Tumor_LogR_segmented,
                     main = paste0(samplename, "\nCoherence of estimated total copy number (TCN) versus post-ASPCF segmented log2(ratio) (L2R)"),
                     xlab = "TCN", ylab = "L2R", xaxs = "i",
                     xlim = my.xlim, ylim = c(-1.5, 1.5), type = "n")
      for (k in sort(unique(cnpTotal))) {
        my.yval <- data$data$Tumor_LogR_segmented[cnpTotal == k]
        my.min <- min(my.yval, na.rm = TRUE)
        my.max <- max(my.yval, na.rm = TRUE)
        rect(xleft = my.xlim[1], ybottom = my.min, xright = my.xlim[2], ytop = my.max, border = adjustcolor(k + 1, alpha.f = 0.25), col = adjustcolor(k + 1, alpha.f = 0.25))
        segments(k, my.min, k, my.max, col = k + 1)
      }
      points(cnpTotal, data$data$Tumor_LogR_segmented, pch = ".", cex = 5, col = cnpTotal + 1)
      abline(h = 0, lty = 2)
      dev.off()
      png(paste0(samplename, ".Rorschach.clown.png"), width = 980, height = 980)
      EaCoN.Rorschard.plot(data = data, cnpTotal = cnpTotal)
      # par(mar = c(1, 1, 1, 1), mfrow = c(5, 5))
      # for (k in 1:length(data$data$ch)) {
      #   graphics::plot(data$data$Tumor_BAF[[1]][data$germline$germlinegenotypes],
      #                  data$data$Tumor_LogR[[1]][data$germline$germlinegenotypes],
      #                  pch = ".", cex = 2, xlim = c(0, 1), ylim = c(-2,2), col = "grey95")
      #   points(data$data$Tumor_BAF[[1]][!data$germline$germlinegenotypes], data$data$Tumor_LogR[[1]][!data$germline$germlinegenotypes], pch = ".", cex = 2, col = "grey50")
      #   kin <- data$data$ch[[k]][!(data$data$ch[[k]] %in% which(data$germline$germlinegenotypes))]
      #   points(data$data$Tumor_BAF[[1]][kin], data$data$Tumor_LogR[[1]][kin], pch = ".", cex = 4, col = cnpTotal[kin] + 1)
      #   try(text(x = 0.5, y = 2, labels = data$data$chrs[k], pos = 1, cex = 2))
      # }
      dev.off()

      tmsg(paste0("    ", round(my.ascat.seg.ascn$goodnessOfFit, digits = 3), " / ", my.ascat.seg.ascn$psi))
      setwd(oridirx)
      return(unname(c(gamma, unlist(my.ascat.seg.ascn$ploidy, use.names = FALSE), my.ascat.seg.ascn$aberrantcellfraction, my.ascat.seg.ascn$goodnessOfFit, my.ascat.seg.ascn$psi)))
    }
  }, stringsAsFactors = FALSE)
  parallel::stopCluster(cls)
  rownames(fit.val) <- gammavec
  # colnames(fit.val) <- c("gamma", "ploidy", "rounded.ploidy", "aberrant.cell.fraction", "GoF", "psi")
  colnames(fit.val) <- c("gamma", "ploidy.ASCAT", "ploidy.Median", "ploidy.Most.width", "ploidy.Width.weighted", "aberrant.cell.fraction", "GoF", "psi")
  if (any(!is.na(fit.val$gamma))) {
    fit.val[,1] <- gammavec
    gammaOpt.idx <- which.max(fit.val$GoF)
    gammaOpt <- fit.val$gamma[gammaOpt.idx]
    write.table(fit.val, file = paste0(odir, "/", samplename, ".gammaEval.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    png(paste0(odir, "/", samplename, ".gammaEval.png"), width = 1850, height = 980)
    par(mfrow = c(3, 1), mar = c(2, 4, 3, 1), cex = 1)
    graphics::plot(fit.val$gamma, fit.val$GoF, xlab = "Gamma", ylab = "Goodness of fit", main = "Goodness of fit curve", type = "b", pch = 20)
    points(fit.val$gamma[gammaOpt.idx], fit.val$GoF[gammaOpt.idx], pch = 20, col = 2)
    abline(v = fit.val$gamma[gammaOpt.idx], lty = 2, col = 2)
    abline(h = fit.val$GoF[gammaOpt.idx], lty = 2, col = 2)
    graphics::plot(fit.val$gamma, fit.val$psi, xlab = "Gamma", ylab = "Psi", main = "Psi curve", type = "b", pch = 20)
    points(fit.val$gamma[gammaOpt.idx], fit.val$psi[gammaOpt.idx], pch = 20, col = 2)
    abline(v = fit.val$gamma[gammaOpt.idx], lty = 2, col = 2)
    ploidy.mat <- as.matrix(fit.val[,2:5])
    plo.ymax <- max(ploidy.mat, na.rm = TRUE)
    graphics::plot(fit.val$gamma, fit.val$ploidy.ASCAT, xlab = "Gamma", ylab = "Ploidy", main = "Ploidy : ASCAT (A=black), Median (M=red), Most width (MW=cyan), Width-weighted (WW=yellow)", type = "b", pch = 20, col = "black", ylim = c(0,plo.ymax), lwd = 2)
    abline(h = seq.int(2, plo.ymax, 2), lty = 2, col = "grey50")
    lines(fit.val$gamma, fit.val$ploidy.Median, type = "b", pch = 20, col = "red", lwd = 2)
    lines(fit.val$gamma, fit.val$ploidy.Most.width, type = "b", pch = 20, col = "cyan", lwd = 2)
    lines(fit.val$gamma, fit.val$ploidy.Width.weighted, type = "b", pch = 20, col = "yellow", lwd = 2)
    abline(v = fit.val$gamma[gammaOpt.idx], lty = 2, col = 2)
    dev.off()
    
    # try(file.rename(from = paste0(out.dir, "/gamma", sprintf("%.2f", gammaOpt)), to = paste0(out.dir, "/gamma", sprintf("%.2f", gammaOpt), "_optimal")))
  } else {
    tmsg("WARNING : ASCN failed for all evaluated gamma values !")
  }
}

## FACETS Total and Allele-Specific Copy Number
ASCN.FACETS <- function(data = NULL, out.dir = getwd(), force = FALSE, ...) {
  
  # setwd("/home/job/Documents/ROSCOFF/Roscoff_2018/TP_CNV/WES/REDUX/A18R.11.17.18/FACETS/L2R")
  # data <- readRDS("A18R.11.17.18.SEG.FACETS.RDS")
  # out.dir <- "/home/job/Documents/ROSCOFF/Roscoff_2018/TP_CNV/WES/REDUX/A18R.11.17.18"
  # force <- TRUE
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # require(foreach)

  
  ## CHECKS
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  odir <- paste0(out.dir, "/FACETS/ASCN")
  if (any(is.null(c(data$data$Tumor_LogR_segmented, data$data$Tumor_BAF_segmented)))) stop(tmsg("No segmentation data found in the provided RDS file !"), call. = FALSE)
  if (dir.exists(odir)) {
    if (force) {
      unlink(odir, recursive = TRUE, force = FALSE)
    } else stop(tmsg(paste0("Output directory [", out.dir, "] already exists !")), call. = FALSE)
  }
  
  samplename <- data$meta$basic$samplename
  ## Extract genome version and load corresponding data
  genome <- data$meta$basic$genome
  genome.pkg <- data$meta$basic$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  cs <- chromobjector(BSg.obj)
  
  tmsg("ASCN modeling (using FACETS) ...")
  ascn.res <- facets::emcncf(x = data$data$additional, ...)
  ascn.res$cncf$lcn.em[is.na(ascn.res$cncf$lcn.em)] <- 0
  
  ## Handling purity
  purity <- ascn.res$purity
  if (!"loglik" %in% names(ascn.res)) ascn.res$loglik <- NA
  
  ## Handling ploidy
  tcn.tbl.ung <- dplyr::as.tbl(cbind(ascn.res$cncf, width = ascn.res$cncf$end - ascn.res$cncf$start + 1))
  tcn.tbl <- dplyr::group_by(tcn.tbl.ung, tcn.em)
  tcn.tbl.prop <- dplyr::summarise(tcn.tbl, tot_width = sum(width))
  
  facets.ploidy <- ascn.res$ploidy
  median.ploidy <- limma::weighted.median(tcn.tbl.prop$tcn.em, tcn.tbl.prop$tot_width)
  baseline.ploidy <- tcn.tbl.prop$tcn.em[which.max(tcn.tbl.prop$tot_width)]
  if (baseline.ploidy == 0) baseline.ploidy <- tcn.tbl.prop$tcn.em[order(tcn.tbl.prop$tot_width, decreasing = TRUE)[2]]
  weighted.ploidy <- sum(tcn.tbl.prop$tcn.em * (tcn.tbl.prop$tot_width / sum(tcn.tbl.prop$tot_width)))
  
  ascn.res$ploidy <- list(facets = as.numeric(facets.ploidy), median = median.ploidy, most.width = baseline.ploidy, width.weighted = weighted.ploidy)
  
  ## Saving ASCN object
  oridir <- getwd()
  dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
  setwd(odir)
  saveRDS(ascn.res, paste0(samplename, ".ASCN.FACETS.RDS"), compress = "bzip2")
  
  ## Generating TCN-CBS
  outfile <- paste0(samplename, ".cn")
  outdf <- ascn.res$cncf
  outdf$chr <- outdf$chrom
  outdf$chrom <- unlist(cs$chr2chrom[outdf$chr])
  outdf$Width <- outdf$end - outdf$start + 1
  outdf$Sample <- samplename
  outdf$nMajor <- outdf$tcn.em - outdf$lcn.em
  # outdf <- outdf[,c(1,2,7,3,4,8,9,5,6)]
  outdf <- outdf[,c(17,15,1,10,11,16,13,18,14)]
  colnames(outdf) <- c(samplename, "Chr", "Chrom", "Start", "End", "Width", "TCN", "nMajor", "nMinor")
  write.table.fast(x = outdf, file = outfile)
  
  ## Generating cellularity + ploidy + metrics file
  outfile <- paste0(samplename, "_model.txt")
  modeldf <- data.frame(
    key =   c("Sample", "LogLik", "dipLogR", "emFlags", "Ploidy.FACETS", "Ploidy.Median", "Ploidy.Most.Width", "Ploidy.Width.weighted", "Cellularity"),
    value = c(samplename, ascn.res$loglik, ascn.res$dipLogR, if(is.null(ascn.res$emflags)) "NA" else ascn.res$emflags, ascn.res$ploidy$facets, ascn.res$ploidy$median, ascn.res$ploidy$most.width, ascn.res$ploidy$width.weighted, ascn.res$purity), stringsAsFactors = FALSE)
  write.table.fast(x = modeldf, file = outfile, header = FALSE)
  
  ## Plots
  ylim <- 6
  ### TCN
  segments.genostart <- cs$chromosomes$chr.length.toadd[ascn.res$cncf$chrom] + ascn.res$cncf$start
  segments.genoend <- cs$chromosomes$chr.length.toadd[ascn.res$cncf$chrom] + ascn.res$cncf$end
  ink <- cs$chromosomes$chrN %in% ascn.res$cncf$chrom
  png(paste0(samplename, ".ASCN.FACETS.png"), width = 1850, height = 980)
  par(mar = c(1, 4, 5, 1), mfrow = c(2, 1))
  # plot(0, 0, type = "n", xlim = c(0,max(segments.genoend)), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " Allele-Specific Copy Number (Gamma = ", gamma, ")\n Cellularity = ", my.ascat.seg.ascn$aberrantcellfraction, " // Ploidy = (A:", round(my.ascat.seg.ascn$ploidy$ascat, digits = 2), ", M:", round(my.ascat.seg.ascn$ploidy$median, digits = 2), ", MW:", round(my.ascat.seg.ascn$ploidy$most.width, digits = 2), ", WW:", round(my.ascat.seg.ascn$ploidy$width.weighted, digits = 2), ") // Goodness of fit = ", round(my.ascat.seg.ascn$goodnessOfFit, digits = 2), "% // Psi = ", my.ascat.seg.ascn$psi, " // nSeg = ", nrow(my.ascat.seg.ascn$segments), " // Predicted gender = ", data$data$gender), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
  plot(0, 0, type = "n", xlim = c(0, cs$genome.length), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " Allele-Specific Copy Number\n Cellularity = ", round(ascn.res$purity, digits = 2), " // Ploidy = (F:", round(ascn.res$ploidy$facets, digits = 2), ", M:", round(ascn.res$ploidy$median, digits = 2), ", MW:", round(ascn.res$ploidy$most.width, digits = 2), ", WW:", round(ascn.res$ploidy$width.weighted, digits = 2), ") // nSeg = ", nrow(ascn.res$cncf), " // Predicted gender = ", data$data$gender), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
  abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
  segments(segments.genostart, ascn.res$cncf$tcn.em - ascn.res$cncf$lcn.em + 0.05, segments.genoend, ascn.res$cncf$tcn.em - ascn.res$cncf$lcn.em + 0.05, pch = ".", col = 2, lwd = 6)
  segments(segments.genostart, ascn.res$cncf$lcn.em - 0.05, segments.genoend, ascn.res$cncf$lcn.em - 0.05, pch = ".", col = 3, lwd = 6)
  up.nMajor <- which((ascn.res$cncf$tcn.em - ascn.res$cncf$lcn.em) > ylim)
  up.nMinor <- which(ascn.res$cncf$lcn.em > ylim)
  dn.nMajor <- which((ascn.res$cncf$tcn.em - ascn.res$cncf$lcn.em) == 0)
  dn.nMinor <- which(ascn.res$cncf$lcn.em == 0)
  if (length(up.nMajor) > 0)
    segments(segments.genostart[up.nMajor],
             (ylim + 0.2) + 0.05, segments.genoend[up.nMajor],
             (ylim + 0.2) + 0.05, pch = ".", col = "orange",
             lwd = 6)
  if (length(up.nMinor) > 0)
    segments(segments.genostart[up.nMinor],
             (ylim + 0.2) - 0.05, segments.genoend[up.nMinor],
             (ylim + 0.2) - 0.05, pch = ".", col = "lightgreen",
             lwd = 6)
  if (length(dn.nMajor) > 0)
    segments(segments.genostart[dn.nMajor],
             0.05, segments.genoend[dn.nMajor], 0.05,
             pch = ".", col = "darkred", lwd = 6)
  if (length(dn.nMinor) > 0)
    segments(segments.genostart[dn.nMinor],
             -0.05, segments.genoend[dn.nMinor], -0.05,
             pch = ".", col = "darkgreen", lwd = 6)
  # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 3, lwd = 2)
  abline(h = 0:ylim, col = "grey75", lty = 3)
  
  try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
  
  # graphics::plot(0, 0, type = "n", xlim = c(0, max(segments.genoend)), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
  graphics::plot(0, 0, type = "n", xlim = c(0, cs$genome.length), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
  abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
  abline(h = ascn.res$ploidy$median, col = "red", lty = 2)
  abline(h = ascn.res$ploidy$most.width, col = "cyan", lty = 2)
  up.nTotal <- which(ascn.res$cncf$tcn.em > ylim)
  dn.nTotal <- which(ascn.res$cncf$tcn.em == 0)
  segments(segments.genostart, ascn.res$cncf$tcn.em, segments.genoend, ascn.res$cncf$tcn.em, pch = ".", col = 4, lwd = 6)
  if (length(up.nTotal) > 0) segments(segments.genostart[up.nTotal], (ylim + 0.2) + 0.05, segments.genoend[up.nTotal], (ylim + 0.2) + 0.05, pch = ".", col = "cyan", lwd = 6)
  if (length(dn.nTotal) > 0) segments(segments.genostart[dn.nTotal], 0, segments.genoend[dn.nTotal], 0, pch = ".", col = "midnightblue", lwd = 6)
  # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 3, lwd = 2)
  abline(h = 0:ylim, col = "grey75", lty = 3)
  try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
  dev.off()
  
  ### TCNvsL2R
  my.xlim <- c(min(ascn.res$cncf$tcn.em, na.rm = TRUE) - 0.5, max(ascn.res$cncf$tcn.em, na.rm = TRUE) + 0.5)
  png(paste0(samplename, ".TCNvsL2R.png"), width = 980, height = 980)
  par(mar = c(4, 4, 5, 1))
  graphics::plot(ascn.res$cncf$tcn.em, ascn.res$cncf$cnlr.median,
                 main = paste0(samplename, "\nCoherence of estimated total copy number (TCN) versus post-ASPCF segmented log2(ratio) (L2R)"),
                 xlab = "TCN", ylab = "L2R", xaxs = "i",
                 xlim = my.xlim, ylim = c(-1.5, 1.5), type = "n")
  for (k in sort(unique(ascn.res$cncf$tcn.em))) {
    my.yval <- ascn.res$cncf$cnlr.median[ascn.res$cncf$tcn.em == k]
    my.min <- min(my.yval, na.rm = TRUE)
    my.max <- max(my.yval, na.rm = TRUE)
    rect(xleft = my.xlim[1], ybottom = my.min, xright = my.xlim[2], ytop = my.max, border = adjustcolor(k + 1, alpha.f = 0.25), col = adjustcolor(k + 1, alpha.f = 0.25))
    segments(k, my.min, k, my.max, col = k + 1)
  }
  points(ascn.res$cncf$tcn.em, ascn.res$cncf$cnlr.median, pch = ".", cex = 5, col = ascn.res$cncf$tcn.em + 1)
  abline(h = 0, lty = 2)
  dev.off()
  
  ### Fail
  # ## Rorschach clown
  # png(paste0(samplename, ".Rorschach.clown.png"), width = 980, height = 980)
  # par(mar = c(1, 1, 1, 1), mfrow = c(5, 5))
  # for (k in 1:length(data$data$ch)) {
  #   graphics::plot(data$data$Tumor_BAF[[1]][data$germline$germlinegenotypes],
  #                  data$data$Tumor_LogR[[1]][data$germline$germlinegenotypes],
  #                  pch = ".", cex = 2, xlim = c(0, 1), ylim = c(-2,2), col = "grey95")
  #   points(data$data$Tumor_BAF[[1]][!data$germline$germlinegenotypes], data$data$Tumor_LogR[[1]][!data$germline$germlinegenotypes], pch = ".", cex = 2, col = "grey50")
  #   kin <- data$data$ch[[k]][!(data$data$ch[[k]] %in% which(data$germline$germlinegenotypes))]
  #   points(data$data$Tumor_BAF[[1]][kin], data$data$Tumor_LogR[[1]][kin], pch = ".", cex = 4, col = ascn.res$cncf$tcn.em[kin] + 1)
  #   try(text(x = 0.5, y = 2, labels = data$data$chrs[k], pos = 1, cex = 2))
  # }
  # dev.off()
  setwd(oridir)
}

## Total and Allele-Specific Copy Number (sequenza)
### *** WARNING : Does not seem very efficient at predicting absolute copy number, unfortunately
### *** and especially for homozygous regions which are totaly wrong...
ASCN.SEQUENZA <- function(data = NULL, max.ploidy = 4, ploidy.step = .1, seg.min.size = 1E+06, out.dir = getwd(), force = FALSE) {
  
  # setwd("/home/job/Documents/ROSCOFF/Roscoff_2018/TP_CNV/WES/REDUX/A18R.11.17.18")
  # data <- readRDS("SEQUENZA/L2R/A18R.11.17.18.SEG.SEQUENZA.RDS")
  # max.ploidy = 4
  # ploidy.step = .1
  # seg.min.size = 1E+06
  # out.dir = getwd()
  # force = TRUE
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  
  ## CHECKS
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  odir <- paste0(out.dir, "/SEQUENZA/ASCN")
  if (any(is.null(c(data$data$Tumor_LogR_segmented, data$data$Tumor_BAF_segmented)))) stop(tmsg("No segmentation data found in the provided RDS file !"), call. = FALSE)
  if (dir.exists(odir)) {
    if (force) {
      unlink(odir, recursive = TRUE, force = FALSE)
    } else stop(tmsg(paste0("Output directory [", out.dir, "] already exists !")), call. = FALSE)
  }
  
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  if (any(is.null(c(data$data$Tumor_LogR_segmented, data$data$Tumor_BAF_segmented)))) stop(tmsg("No segmentation data found in the provided RDS file !"), call. = FALSE)
  
  samplename <- data$meta$basic$samplename
  
  ## Extract genome version and load corresponding data
  genome <- data$meta$basic$genome
  genome.pkg <- data$meta$basic$genome.pkg
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  tmsg(paste0("Loading ", genome.pkg, " ..."))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  cs <- chromobjector(BSg.obj)
  
  tmsg("ASCN modeling (using sequenza) ...")
  
  ## Loading segments table
  seg.df <- data$data$additional
  
  ## Fixing segments of length 1
  for (k in unique(seg.df$chrom)) {
    kidx <- seg.df$chrom == k
    kstart <- which(kidx)[1]
    # len1.idx <- which(kidx & seg.df$Probes == 1)
    len1.idx <- which(kidx & seg.df$n.probes < 3)
    if (length(len1.idx) > 0) {
      if (len1.idx[1] == kstart) {
        seg.df$start.pos[kstart+1] <- seg.df$start.pos[kstart]
        seg.df$n.probes[kstart+1] <- seg.df$n.probes[kstart+1] + seg.df$n.probes[kstart]
        seg.df <- seg.df[-kstart,]
        len1.idx <- len1.idx[-1]
      }
      if (length(len1.idx) > 0) {
        seg.df$end.pos[len1.idx-1] <- seg.df$end.pos[len1.idx]
        seg.df$n.probes[len1.idx-1] <- seg.df$n.probes[len1.idx-1] + seg.df$n.probes[len1.idx]
        seg.df <- seg.df[-len1.idx,]
      }
    }
  }
  
  ## Computing width
  seg.df$width <- seg.df$end.pos - seg.df$start.pos + 1
  seg.df$weight <- 0 + round(seg.df$width / 1E+06, 0)
  
  ## Ploidy model
  seg.min.size <- 1E+06
  max.ploidy <- 4
  ploidy.step <- .1
  sms.idx <- seg.df$width > seg.min.size
  CP <- sequenza::baf.model.fit(Bf = seg.df$BAF.mean[sms.idx], depth.ratio = 2^seg.df$logR.mean[sms.idx], weight.ratio = seg.df$weight[sms.idx], weight.Bf = seg.df$weight[sms.idx], avg.depth.ratio = 1, cellularity = seq(0.1,1,0.01), ploidy = seq(0,max.ploidy,ploidy.step))
  confint <- sequenza::get.ci(CP)
  ploidy   <- confint$max.ploidy
  cellularity <- confint$max.cellularity
  rm(CP)
  
  ## Modeling ASCN
  seq.ascn <- sequenza::baf.bayes(Bf = seg.df$BAF.mean, depth.ratio = 2^seg.df$logR.mean, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = 1)
  seg.df$chr <- seg.df$chrom
  seg.df$chr <- unlist(cs$chrom2chr[seg.df$chr])
  ascn.res <- list(data = cbind(seg.df, seq.ascn))
  rm(seq.ascn)
  
  
  ## Computing other ploidy measures
  ascn.res$purity <- cellularity
  seg.tbl <- dplyr::as.tbl(ascn.res$data)
  seg.tbl <- dplyr::group_by(seg.tbl, CNt)
  seg.tbl.prop <- dplyr::summarise(seg.tbl, tot_width = sum(width))
  median.ploidy <- limma::weighted.median(seg.tbl.prop$CNt, seg.tbl.prop$tot_width)
  # baseline.ploidy <- seg.tbl.prop$CNt[which.max(seg.tbl.prop$tot_width)]
  baseline.ploidy <- seg.tbl.prop$CNt[which.max(seg.tbl.prop$tot_width)]
  if (baseline.ploidy == 0) baseline.ploidy <- seg.tbl.prop$CNt[order(seg.tbl.prop$tot_width, decreasing = TRUE)[2]]
  weighted.ploidy <- sum(seg.tbl.prop$CNt * (seg.tbl.prop$tot_width / sum(seg.tbl.prop$tot_width)))
  ascn.res$ploidy <- list(sequenza = as.numeric(ploidy), median = median.ploidy, most.width = baseline.ploidy, width.weighted = weighted.ploidy)
  
  
  ## Saving ASCN object
  oridir <- getwd()
  dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
  setwd(odir)
  saveRDS(ascn.res, paste0(samplename, ".ASCN.SEQUENZA.RDS"), compress = "bzip2")
  
  ## Generating TCN-CBS
  outfile <- paste0(samplename, ".cn")
  # outdf <- ascn.res$data[,c(1,2,15,4,5,9,11:13)]
  outdf <- ascn.res$data[,c(1,11,2,4,5,9,12:14)]
  colnames(outdf) <- c(samplename, "Chr", "Chrom", "Start", "End", "Width", "TCN", "nMajor", "nMinor")
  write.table.fast(x = outdf, file = outfile)
  
  ## Generating cellularity + ploidy + metrics file
  outfile <- paste0(samplename, "_model.txt")
  modeldf <- data.frame(
    key =   c("Sample", "Ploidy.SEQUENZA", "Ploidy.Median", "Ploidy.Most.Width", "Ploidy.Width.weighted", "Cellularity"),
    value = c(samplename, ascn.res$ploidy$sequenza, ascn.res$ploidy$median, ascn.res$ploidy$most.width, ascn.res$ploidy$width.weighted, ascn.res$purity), stringsAsFactors = FALSE)
  write.table.fast(x = modeldf, file = outfile, header = FALSE)
  
  ## Plots
  ylim <- 6
  ### TCN
  segments.genostart <- cs$chromosomes$chr.length.toadd[ascn.res$data$chr] + ascn.res$data$start
  segments.genoend <- cs$chromosomes$chr.length.toadd[ascn.res$data$chr] + ascn.res$data$end
  ink <- cs$chromosomes$chrN %in% ascn.res$data$chr
  png(paste0(samplename, ".ASCN.SEQUENZA.png"), width = 1850, height = 980)
  par(mar = c(1, 4, 5, 1), mfrow = c(2, 1))
  plot(0, 0, type = "n", xlim = c(0, cs$genome.length), xaxs = "i", ylim = c(0,(ylim + 0.1)), main = paste0(samplename, " Allele-Specific Copy Number\n Cellularity = ", round(ascn.res$purity, digits = 2), " // Ploidy = (S:", round(ascn.res$ploidy$sequenza, digits = 2), ", M:", round(ascn.res$ploidy$median, digits = 2), ", MW:", round(ascn.res$ploidy$most.width, digits = 2), ", WW:", round(ascn.res$ploidy$width.weighted, digits = 2), ") // nSeg = ", nrow(ascn.res$cncf), " // Predicted gender = ", data$data$gender), xlab = "Genomic position", ylab = "ASCN", xaxt = "n", cex.main = 2)
  abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
  segments(segments.genostart, ascn.res$data$A  + 0.05, segments.genoend, ascn.res$data$A + 0.05, pch = ".", col = 2, lwd = 6)
  segments(segments.genostart, ascn.res$data$B - 0.05, segments.genoend, ascn.res$data$B - 0.05, pch = ".", col = 3, lwd = 6)
  up.nMajor <- which(ascn.res$data$A > ylim)
  up.nMinor <- which(ascn.res$data$B > ylim)
  dn.nMajor <- which(ascn.res$data$A == 0)
  dn.nMinor <- which(ascn.res$data$B == 0)
  if (length(up.nMajor) > 0)
    segments(segments.genostart[up.nMajor],
             (ylim + 0.2) + 0.05, segments.genoend[up.nMajor],
             (ylim + 0.2) + 0.05, pch = ".", col = "orange",
             lwd = 6)
  if (length(up.nMinor) > 0)
    segments(segments.genostart[up.nMinor],
             (ylim + 0.2) - 0.05, segments.genoend[up.nMinor],
             (ylim + 0.2) - 0.05, pch = ".", col = "lightgreen",
             lwd = 6)
  if (length(dn.nMajor) > 0)
    segments(segments.genostart[dn.nMajor],
             0.05, segments.genoend[dn.nMajor], 0.05,
             pch = ".", col = "darkred", lwd = 6)
  if (length(dn.nMinor) > 0)
    segments(segments.genostart[dn.nMinor],
             -0.05, segments.genoend[dn.nMinor], -0.05,
             pch = ".", col = "darkgreen", lwd = 6)
  # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 3, lwd = 2)
  abline(h = 0:ylim, col = "grey75", lty = 3)
  
  try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
  
  graphics::plot(0, 0, type = "n", xlim = c(0, cs$genome.length), xaxs = "i", ylim = c(0, (ylim + 0.1)), main = paste0(samplename, " Total Copy Number"), xlab = "Genomic position", ylab = "TCN", xaxt = "n", cex.main = 2)
  abline(v = c(segments.genostart, segments.genoend), col = "grey90", lwd = 1)
  abline(h = ascn.res$ploidy$median, col = "red", lty = 2)
  abline(h = ascn.res$ploidy$most.width, col = "cyan", lty = 2)
  up.nTotal <- which(ascn.res$data$CNt > ylim)
  dn.nTotal <- which(ascn.res$data$CNt == 0)
  segments(segments.genostart, ascn.res$data$CNt, segments.genoend, ascn.res$data$CNt, pch = ".", col = 4, lwd = 6)
  if (length(up.nTotal) > 0) segments(segments.genostart[up.nTotal], (ylim + 0.2) + 0.05, segments.genoend[up.nTotal], (ylim + 0.2) + 0.05, pch = ".", col = "cyan", lwd = 6)
  if (length(dn.nTotal) > 0) segments(segments.genostart[dn.nTotal], 0, segments.genoend[dn.nTotal], 0, pch = ".", col = "midnightblue", lwd = 6)
  # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 3, lwd = 2)
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 3, lwd = 2)
  abline(h = 0:ylim, col = "grey75", lty = 3)
  try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, labels = cs$chromosomes$chrom[ink], pos = 1, cex = 1))
  dev.off()
  
  ### TCNvsL2R
  my.xlim <- c(min(ascn.res$data$CNt, na.rm = TRUE) - 0.5, max(ascn.res$data$CNt, na.rm = TRUE) + 0.5)
  png(paste0(samplename, ".TCNvsL2R.png"), width = 980, height = 980)
  par(mar = c(4, 4, 5, 1))
  graphics::plot(ascn.res$data$CNt, ascn.res$data$logR.mean,
                 main = paste0(samplename, "\nCoherence of estimated total copy number (TCN) versus segmented log2(ratio) (L2R)"),
                 xlab = "TCN", ylab = "L2R", xaxs = "i",
                 xlim = my.xlim, ylim = c(-1.5, 1.5), type = "n")
  for (k in sort(unique(ascn.res$data$CNt))) {
    my.yval <- ascn.res$data$logR.mean[ascn.res$data$CNt == k]
    my.min <- min(my.yval, na.rm = TRUE)
    my.max <- max(my.yval, na.rm = TRUE)
    rect(xleft = my.xlim[1], ybottom = my.min, xright = my.xlim[2], ytop = my.max, border = adjustcolor(k + 1, alpha.f = 0.25), col = adjustcolor(k + 1, alpha.f = 0.25))
    segments(k, my.min, k, my.max, col = k + 1)
  }
  points(ascn.res$data$CNt, ascn.res$data$logR.mean, pch = ".", cex = 5, col = ascn.res$data$CNt + 1)
  abline(h = 0, lty = 2)
  dev.off()
  
  setwd(oridir)
}


## Run ASCN segmentation, from a file
ASCN.ff <- function(RDS.file = NULL, ...) {
  if (is.null(RDS.file)) stop(tmsg("A RDS file is needed !"), call. = FALSE)
  if (!file.exists(RDS.file)) stop(tmsg(paste0("Could not find RDS file ", RDS.file, " !")), call. = FALSE)

    ## Data loading
  tmsg(paste0("Loading data from ", RDS.file, " ..."))
  my.data <- readRDS(RDS.file)
  segmenter <- toupper(my.data$meta$eacon$segmenter)
  tmsg(paste0("Provided RDS.file was segmented with ", my.data$meta$eacon$segmenter, "."))
  
  ## Calling
  out.dir <- sub(pattern = paste0(toupper(segmenter), "/L2R"), replacement = "", x = dirname(RDS.file))
  if (out.dir == "") out.dir <- "."
  do.call(paste0("ASCN.", toupper(segmenter)), list(data = my.data, out.dir = out.dir, ...))
}

## Run ASCN.ff() in batch mode
ASCN.ff.Batch <- function(RDS.files = list.files(path = getwd(), pattern = "\\.SEG\\..*\\.RDS$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE), nthread = 1, cluster.type = "PSOCK", ...) {
  if (length(RDS.files) == 0) stop("A list of RDS files is required !", call. = FALSE)
  message("Running EaCoN.ASCN.ff() in batch mode ...")
  message(paste0(" Found ", length(RDS.files), " files to process."))
  current.bitmapType <- getOption("bitmapType")
  if (length(RDS.files) < nthread) nthread <- length(RDS.files)
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  r <- ""
  eacon.batchres <- foreach::foreach(r = seq_along(RDS.files), .inorder = TRUE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    ASCN.ff(RDS.file = RDS.files[r], ...)
  }
  parallel::stopCluster(cl)
}

## Generate the HTML report
Annotate <- function(data = NULL, refGene.table = NULL, targets.table = NULL, report = TRUE, author.name = "", out.dir = getwd(), solo = FALSE, ldb = "/mnt/data_cigogne/bioinfo/") {
  
  # setwd("/home/job/Documents/ROSCOFF/Roscoff_2018/TP_CNV/WES/REDUX/A18R.11.17.18/ASCAT/L2R")
  # data <- readRDS("A18R.11.17.18.SEG.ASCAT.RDS")
  # targets.table <- NULL
  # out.dir <- getwd()
  # refGene.table = NULL
  # solo = FALSE
  # report = TRUE
  # ldb = "/mnt/data_cigogne/bioinfo/"
  # source("/home/job/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("/home/job/git_gustaveroussy/EaCoN/R/plot_functions.R")
  # require(foreach)
  
  oridir <- getwd()
  
  `%do%` <- foreach::"%do%"
  
  if (!is.list(data)) stop(tmsg("data should be a list !"), call. = FALSE)
  
  valid.genomes <- get.valid.genomes()
  my.ascat.seg <- data$data
  
  samplename <- data$meta$basic$samplename
  genome <- data$meta$basic$genome
  manufacturer <- data$meta$basic$manufacturer
  source <- data$meta$basic$source
  
  tmsg("Loading genome data ...")
  ## No refGene provided ...
  if (is.null(refGene.table)) {
    ## ... but some in bank : loading !
    if (genome %in% names(valid.genomes)) {
      self.pkg.name <- "EaCoN"
      data(list = paste0("refGene_cl_", genome), package = self.pkg.name, envir = environment())
      data(list = genome, package = self.pkg.name, envir = environment())
    } else {
      ## ... and none in bank : no genes !
      tmsg("WARNING : No refGene table provided, and none banked for current genome : Gene information won't exist in the report!")
    }
  } else if (!file.exists(refGene.table)) {
    ## A refGene was provided, but could not be found : stopping !
    stop(tmsg(paste0("Could not open file ", refGene.table, " !")), call. = FALSE)
  } else {
    ## A refGene was provided and file exists : loading !
    rg.df <- read.table.fast(file = refGene.table)
    rg.df <- rg.df[grep(pattern = "^chr([0-9]+|X|Y)$", x = rg.df$chrom),]
    gen.df <- foreach::foreach(k = sort(unique(rg.df$chrom)), .combine = "rbind") %do% {
      rg.k <- rg.df[rg.df$chrom == k, ]
      rg.gr <- GenomicRanges::GRanges(seqnames = rg.k$name2, IRanges::IRanges(rg.k$txStart, rg.k$txEnd), strand = rg.k$strand)
      rdx.k <- as.data.frame(GenomicRanges::reduce(rg.gr, ignore.strand = FALSE))
      rdx.df <- data.frame(symbol = as.character(rdx.k$seqnames), chrom = k, chrN = as.numeric(cs$chrom2chr[k]), start = rdx.k$start, end = rdx.k$end, width = rdx.k$width, strand = as.character(rdx.k$strand), stringsAsFactors = FALSE)
      return(rdx.df)
    }
    gen.df <- gen.df[order(gen.df$chrN, gen.df$start, gen.df$end),]
  }
  
  if (!("cbs" %in% names(data))) stop(tmsg("CBS slot not found in RDS object ! Are you sure it is a valid one ?"), call. = FALSE)
  cbs.df <- foreach::foreach(seg = 1:nrow(data$cbs$cut), .combine = "rbind") %do% {
    ingenz <- if(exists("gen.df")) gen.df$symbol[gen.df$chrN == data$cbs$cut$Chr[seg] & gen.df$start <= data$cbs$cut$End[seg] & gen.df$end >= data$cbs$cut$Start[seg]] else ingenz <- "NA"
    ingenz.len <- if(exists("gen.df")) length(ingenz) else "NA"
    return(data.frame(data$cbs$cut[seg, ], Genes = ingenz.len, Symbol = paste0(ingenz, collapse = ","), stringsAsFactors = FALSE))
  }
  
  tmsg("Building L2R segmentation table ...")
  l2r.segments <- foreach::foreach(k = 1:length(my.ascat.seg$ch), .combine = "rbind") %do% {
    segrle <- rle(my.ascat.seg$Tumor_LogR_segmented[my.ascat.seg$ch[[k]],1])
    segdf <- data.frame(Probes = segrle$lengths, Value = segrle$values, stringsAsFactors = FALSE)
    # segdf$Chr <- k
    segdf$Chr <- unlist(cs$chrom2chr[data$data$chrs[k]])
    segdf$Chrom <- my.ascat.seg$chrs[[k]]
    probecs <- cumsum(segdf$Probes)
    segdf$End.idx <- as.integer(cumsum(segdf$Probes) + my.ascat.seg$ch[[k]][1] - 1)
    segdf$Start.idx <- as.integer(c(my.ascat.seg$ch[[k]][1], segdf$End.idx[-nrow(segdf)] + 1))
    segdf$Start <- as.integer(my.ascat.seg$SNPpos$pos[segdf$Start.idx])
    segdf$End <- as.integer(my.ascat.seg$SNPpos$pos[segdf$End.idx])
    segdf <- segdf[, c(3, 4, 7, 8, 2, 1, 6, 5)]
  }
  
  tmsg("Building BAF segmentation table ...")
  bafpos <- my.ascat.seg$SNPpos[rownames(my.ascat.seg$Tumor_LogR_segmented) %in% rownames(my.ascat.seg$Tumor_BAF_segmented[[1]]),]
  bafpos$ProbeSet <- rownames(my.ascat.seg$Tumor_BAF_segmented[[1]])
  bafpos$BAF <- my.ascat.seg$Tumor_BAF_segmented[[1]][,1]
  bafpos$chrs <- as.character(bafpos$chrs)
  bafkend <- vapply(unique(bafpos$chrs), function(x) { max(which(bafpos$chrs == x))}, 1)
  bafbreaks <- cumsum(rle(bafpos$BAF)$lengths)
  baf.seg <- data.frame(end.idx = sort(unique(c(bafkend, bafbreaks))), stringsAsFactors = FALSE)
  baf.seg$start.idx <- c(1, baf.seg$end.idx[-nrow(baf.seg)]+1)
  baf.seg$length.idx <- baf.seg$end.idx - baf.seg$start.idx +1
  baf.seg$start.probeset <- bafpos$ProbeSet[baf.seg$start.idx]
  baf.seg$end.probeset <- bafpos$ProbeSet[baf.seg$end.idx]
  baf.seg$chrA <- bafpos$chrs[baf.seg$start.idx]
  # baf.seg$Chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", baf.seg$chrA)]) else unlist(cs$chrom2chr[baf.seg$chrA])
  baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
  # baf.seg$Chr <- unlist(cs$chrom2chr[paste0("chr", baf.seg$chrA)])
  baf.seg$Start <- bafpos$pos[baf.seg$start.idx]
  baf.seg$End <- bafpos$pos[baf.seg$end.idx]
  baf.seg$Width <- baf.seg$End - baf.seg$Start + 1
  baf.seg$Value <- bafpos$BAF[baf.seg$start.idx]
  
  ## BAF calling
  # baf.homocut <- as.numeric(data$meta$eacon[["BAF-segments-homo-limit"]])
  baf.homocut <- as.numeric(data$meta$eacon[["BAF.segments.homo.limit"]])
  
  ### Calling
  baf.seg$Status <- "Unbalanced"
  baf.homo <- sort(which(baf.seg$Value <= baf.homocut))
  if (length(baf.homo) > 0) baf.seg$Status[baf.homo] <- "Homo"
  baf.hetero <- sort(which(baf.seg$Value == .5))
  if (length(baf.hetero) > 0) baf.seg$Status[baf.hetero] <- "Hetero"
  
  ## Writing ACBS (cut)
  setwd(out.dir)
  if (exists("gen.df")) {
    tmsg("Building Targets table ...")
    write.table(cbs.df, file = paste0(samplename, ".Cut.acbs"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  cbsNC.df <- data$cbs$nocut
  cbsNC.df <- cbind(cbsNC.df, cbs.df[, c(7, 8)])
  cbsNC.df$Status <- "Normal"
  cbsNC.df$Status[cbs.df$Log2Ratio > 0] <- "Gain"
  cbsNC.df$Status[cbs.df$Log2Ratio < 0] <- "Loss"
  
  ## Loading targets
  if (is.null(targets.table)) {
    if (data$meta$basic$species %in% c("Homo sapiens")) {
      self.pkg.name <- "EaCoN"
      data("targetlist_475", package = self.pkg.name, envir = environment())
    }
  } else load(targets.table)
  if(exists("targetlist")) {
    gen.targ <- gen.df[gen.df$symbol %in% targetlist, ]
    targ.regz <- foreach::foreach(g = 1:nrow(gen.targ), .combine = "rbind") %do% {
      ginreg.idx <- which(cbsNC.df$Chr == gen.targ$chrN[g] & cbsNC.df$Start <= gen.targ$end[g] & cbsNC.df$End >= gen.targ$start[g])
      ginreg <- foreach::foreach(gg = ginreg.idx, .combine = "rbind") %do% {
        match.start <- max(cbsNC.df$Start[gg], gen.targ$start[g])
        match.end <- min(cbsNC.df$End[gg], gen.targ$end[g])
        return(cbind(gen.targ[g, c(1, 2, 4:7)], match.start = match.start, match.end = match.end, cbsNC.df[gg,c(6, 9)], l2rwidth = cbsNC.df$End[gg] - cbsNC.df$Start[gg] + 1))
      }
      return(ginreg)
    }
    targ.regz <- foreach::foreach(g = 1:nrow(targ.regz), .combine = "rbind") %do% {
      # ginreg.idx <- which(paste0("chr", baf.seg$chrA) == targ.regz$chrom[g] & baf.seg$Start <= targ.regz$match.end[g] & baf.seg$End >= targ.regz$match.start[g])
      ginreg.idx <- which(baf.seg$chrA == targ.regz$chrom[g] & baf.seg$Start <= targ.regz$match.end[g] & baf.seg$End >= targ.regz$match.start[g])
      ginreg <- foreach::foreach(gg = ginreg.idx, .combine = "rbind") %do% {
        match.start <- max(baf.seg$Start[gg], targ.regz$match.start[g])
        match.end <- min(baf.seg$End[gg], targ.regz$match.end[g])
        return(cbind(targ.regz[g,c(1:6)], match.start = match.start, match.end = match.end, match.width = match.end - match.start +1, targ.regz[g,c(10,9,11)], baf.seg[gg,c(12,11)], bafwidth = baf.seg$End[gg] - baf.seg$Start[gg] + 1))
      }
      return(ginreg)
    }
    targ.regz$Cytoband <- vapply(1:nrow(targ.regz), function(x) {
      scb <- cs$cytobands$chrom == targ.regz$chrom[x] & cs$cytobands$start <= targ.regz$start[x] & cs$cytobands$end >= targ.regz$start[x]
      return(paste0(cs$cytobands$chrA[scb], cs$cytobands$cytoband[scb]))
    }, "a")
    
    targ.regz <- targ.regz[order(as.numeric(unlist(cs$chrom2chr[targ.regz$chrom])), targ.regz$match.start, targ.regz$match.end), c(1:6,16,7:15)]
    colnames(targ.regz) <- c("Target Symbol", "Chr", "Gene Start", "Gene End", "Gene Width", "Gene Strand", "Gene Cytoband", "Match Start", "Match End", "Match Width", "L2R Status", "L2R Value", "L2R Segment Width", "BAF Status", "BAF Value", "BAF Segment Width")
    # write.table(targ.regz, file = paste0(out.dir, "/", samplename, ".TargetGenes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(targ.regz, file = paste0(samplename, ".TargetGenes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  tmsg("Building Truncated table ...")
  # gen.trunk.idx <- which(gen.df$chrN == cbsNC.df$Chr & gen.df$start < cbsNC.df$End & gen.df & gen.df$end > cbsNC.df$Start)
  gen.trunk.idx.list <- sapply(1:nrow(gen.df), function(g) {
    gidx <- which(cbsNC.df$Chr == gen.df$chrN[g] & ((cbsNC.df$Start > gen.df$start[g] & cbsNC.df$Start < gen.df$end[g]) | (cbsNC.df$End > gen.df$start[g] & cbsNC.df$End < gen.df$end[g])))
    return(length(gidx))
    # if(length(gen.trunk.idx) > 0) return(gidx) else return(NULL)
  })
  gen.trunk.idx <- which(gen.trunk.idx.list > 1)
  if (length(gen.trunk.idx) > 0) {
    gen.trunk <- gen.df[gen.trunk.idx,]
    trunc.regz <- foreach::foreach(g = 1:nrow(gen.trunk), .combine = "rbind") %do% {
      ginreg.idx <- which(cbsNC.df$Chr == gen.trunk$chrN[g] & cbsNC.df$Start <= gen.trunk$end[g] & cbsNC.df$End >= gen.trunk$start[g])
      ginreg <- foreach::foreach(gg = ginreg.idx, .combine = "rbind") %do% {
        match.start <- max(cbsNC.df$Start[gg], gen.trunk$start[g])
        match.end <- min(cbsNC.df$End[gg], gen.trunk$end[g])
        return(cbind(gen.trunk[g, c(1, 2, 4:7)], match.start = match.start, match.end = match.end, cbsNC.df[gg,c(6, 9)], l2rwidth = cbsNC.df$End[gg] - cbsNC.df$Start[gg] + 1))
      }
      return(ginreg)
    }
    trunc.regz <- foreach::foreach(g = 1:nrow(trunc.regz), .combine = "rbind") %do% {
      # ginreg.idx <- which(paste0("chr", baf.seg$chrA) == trunc.regz$chrom[g] & baf.seg$Start <= trunc.regz$match.end[g] & baf.seg$End >= trunc.regz$match.start[g])
      # ginreg.idx <- which(paste0("chr", baf.seg$chrA) == trunc.regz$chrom[g] & ((baf.seg$Start <= trunc.regz$match.start[g] & baf.seg$End >= trunc.regz$match.start[g]) | (baf.seg$Start <= trunc.regz$match.end[g] & baf.seg$End >= trunc.regz$match.end[g])))
      # ginreg.idx <- which(paste0("chr", baf.seg$chrA) == trunc.regz$chrom[g] & baf.seg$Start <= trunc.regz$match.end[g] & baf.seg$End >= trunc.regz$match.start[g])
      ginreg.idx <- which(baf.seg$chrA == trunc.regz$chrom[g] & baf.seg$Start <= trunc.regz$match.end[g] & baf.seg$End >= trunc.regz$match.start[g])
      if (length(ginreg.idx > 0)) {
        ginreg <- foreach::foreach(gg = ginreg.idx, .combine = "rbind") %do% {
          match.start <- max(baf.seg$Start[gg], trunc.regz$match.start[g])
          match.end <- min(baf.seg$End[gg], trunc.regz$match.end[g])
          return(cbind(trunc.regz[g,c(1:6)], match.start = match.start, match.end = match.end, match.width = match.end - match.start +1, trunc.regz[g,c(10,9,11)], baf.seg[gg,c(12,11)], bafwidth = baf.seg$End[gg] - baf.seg$Start[gg] + 1))
        }
      } else {
        ginreg <- data.frame(trunc.regz[g,c(1:8)], match.width = trunc.regz$match.end[g] - trunc.regz$match.start[g] + 1, trunc.regz[g,c(10,9,11)], Status = NA, Value = NA, bafwidth = NA, stringsAsFactors = FALSE, check.names = FALSE)
      }
      return(ginreg)
    }
    trunc.regz$Cytoband <- vapply(1:nrow(trunc.regz), function(x) {
      scb <- cs$cytobands$chrom == trunc.regz$chrom[x] & cs$cytobands$start <= trunc.regz$start[x] & cs$cytobands$end >= trunc.regz$start[x]
      return(paste0(cs$cytobands$chrA[scb], cs$cytobands$cytoband[scb]))
    }, "a")
    trunc.regz <- trunc.regz[order(as.numeric(unlist(cs$chrom2chr[trunc.regz$chrom])), trunc.regz$match.start, trunc.regz$match.end), c(1:6,16,7:15)]
    colnames(trunc.regz) <- c("Gene Symbol", "Chr", "Gene Start", "Gene End", "Gene Width", "Gene Strand", "Gene Cytoband", "Match Start", "Match End", "Match Width", "L2R Status", "L2R Value", "L2R Segment Width", "BAF Status", "BAF Value", "BAF Segment Width")
    # write.table(trunc.regz, file = paste0(out.dir, "/", samplename, ".TruncatedGenes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(trunc.regz, file = paste0(samplename, ".TruncatedGenes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  if (report) {
    tmsg("Preparing HTML report ...")
    
    ## Locating RMD
    self.pkg.name <- "EaCoN"
    rmd.path <- system.file("extdata", "html_report.Rmd", package = self.pkg.name)
    
    if ( (source == "microarray") & (manufacturer == "Affymetrix") ) {
      
      ## Getting meta data
      tmsg("Getting META keys ...")
      meta.tags <- data.frame(
        key = c("affymetrix-algorithm-param-option-set-analysis-name",
                "affymetrix-array-type",
                "affymetrix-array-id",
                "affymetrix-array-barcode",
                "affymetrix-scanner-id",
                "affymetrix-scan-date",
                "affymetrix-chipsummary-MAPD",
                "affymetrix-chipsummary-iqr",
                "affymetrix-chipsummary-snp-qc",
                "affymetrix-chipsummary-waviness-sd",
                "predicted.gender"
        ),
        entry = c("analysis",
                  "analysis",
                  "array",
                  "array",
                  "acquisition",
                  "acquisition",
                  "analysis",
                  "analysis",
                  "analysis",
                  "analysis",
                  "basic"
        ),
        sig = c("Sample Name",
                "Array Type",
                "Array ID",
                "Array Barcode",
                "Scanner ID",
                "Scan Date",
                "MAPD",
                "IQR",
                "SNPQC",
                "WavinessSd",
                "Predicted Gender"
        ),
        stringsAsFactors = FALSE
      )
      
      chip.list <- names(data$meta$affy)[names(data$meta$affy) != "analysis"]
      meta.tags$val <- foreach(t = seq_len(nrow(meta.tags)), .combine = "c") %do% {
        if (meta.tags$entry[t] == "analysis") {
          ext.meta <- getmeta(key = meta.tags$key[t], meta = data$meta$affy$analysis)
        } else if (meta.tags$entry[t] %in% c("eacon", "basic")) {
          ext.meta <- getmeta(key = meta.tags$key[t], meta = data$meta[[meta.tags$entry[t]]])
        } else {
          ext.meta <- paste0(sapply(chip.list, function(p) { getmeta(key = meta.tags$key[t], meta = data$meta$affy[[p]][[meta.tags$entry[t]]])}), collapse = " ")
        }
        return(ext.meta)
      }
      
      arraytype <- meta.tags$val[meta.tags$key == "affymetrix-array-type"]
      meta.tags$val[meta.tags$key == "affymetrix-scan-date"] <- gsub(replacement = "/", x = gsub(replacement = "_", x = gsub(replacement = "", x = meta.tags$val[meta.tags$key == "affymetrix-scan-date"], pattern = "Z", ignore.case = FALSE), pattern = "T", ignore.case = FALSE), pattern = "-", ignore.case = FALSE)
      meta.tags$val[meta.tags$key == "affymetrix-array-id"] <- gsub(replacement = "_", x = meta.tags$val[meta.tags$key == "affymetrix-array-id"], pattern = "-", ignore.case = FALSE)
      
      ## Loading CEL files for intensity data/plot
      tmsg("Plotting CEL file(s) intensity map ...")
      # intplotf <- paste0(out.dir, "/", samplename, ".INT.png")
      intplotf <- paste0(samplename, ".INT.png")
      narray <- length(names(data$CEL))
      plot.grid <- num2mat(narray)
      png(intplotf, width = 600*plot.grid[1], height = 600*plot.grid[2])
      par(mar=c(1,1,1,1),xaxs='i',yaxs='i',las=1, mfrow=rev(plot.grid))
      co <- grDevices::colorRampPalette(c('black','white'))
      for(p in seq_len(narray)) {
        celobj <- celstruc(celdat = data$CEL[[p]])
        graphics::image(log2(celobj$intensities[,ncol(celobj$intensities):1]),col=co(20),zlim=c(5,14),axes=F,xlab='',ylab='',useRaster=T)
        box(lwd=3)
      }
      dev.off()
      intplotf <- tools::file_path_as_absolute(intplotf)
      # meta.tags <- rbind(meta.tags, data.frame(key = c("a2p-cel-median-intensity", "a2p-cel-outliers", "a2p-cel-outliers-pc"), sig = c("Median Raw Intensity", "# Outliers", "Outliers (%)"), val = c(celobj$int.median, celobj$n.outliers, paste0(round(celobj$pc.outliers / 100, digits = 3), " %")), stringsAsFactors = FALSE))
      
      ## Building array & QC tables
      tmsg("Building Array & QC tables ...")
      array.df <- data.frame(t(meta.tags$val[1:6]), stringsAsFactors = FALSE)
      colnames(array.df) <- meta.tags$sig[1:6]
      qc.df <- data.frame(t(meta.tags$val[7:11]), stringsAsFactors = FALSE)
      colnames(qc.df) <- meta.tags$sig[7:11]
      qc.df[["Median Raw Intensity"]] <- paste0(vapply(chip.list, function(p) { return(median(data$CEL[[p]]$intensities, na.rm = TRUE)) }, 1), collapse = " ; ")
      qc.df[["Outlier Probes"]] <- paste0(vapply(chip.list, function(p) { return(length(data$CEL[[p]]$outliers)) }, 1), collapse = " ; ")
      
      qc.df$mapd.flag <- if(is.na(as.numeric(qc.df$MAPD))) "#e0e0e0" else if(as.numeric(qc.df$MAPD) < .2) "#00cc00" else if(as.numeric(qc.df$MAPD) < .25) "#ff9900" else "#f7543b"
      qc.df$snpqc.flag <- if(is.na(as.numeric(qc.df$SNPQC))) "#e0e0e0" else if(as.numeric(qc.df$SNPQC) >= 26) "#00cc00" else if(as.numeric(qc.df$SNPQC) >= 15) "#ff9900" else "#f7543b"
      qc.df$waviness.flag <- if(is.na(as.numeric(qc.df$WavinessSd))) "#e0e0e0" else if(as.numeric(qc.df$WavinessSd) < .12) "#00cc00" else "#f7543b"
    }
    ## PLOTS
    ### Genomic L2R plot
    tmsg("Plotting L2R (geno) ...")
    mad.diff <- abs(diff(my.ascat.seg$Tumor_LogR[,1]))
    my.mad <- stats::median(mad.diff[mad.diff>0], na.rm = TRUE)
    
    # g.cut <- getmeta(key = "L2R-segments-gain-cutoff", meta = data$meta$eacon)
    g.cut <- data$meta$eacon[["L2R-segments-gain-cutoff"]]
    # l.cut <- getmeta(key = "L2R-segments-loss-cutoff", meta = data$meta$eacon)
    l.cut <- data$meta$eacon[["L2R-segments-loss-cutoff"]]
    
    gain.idx <- which(l2r.segments$Value > g.cut)
    loss.idx <- which(l2r.segments$Value < l.cut)
    normal.idx <- which(l2r.segments$Value >= l.cut & l2r.segments$Value <= g.cut)
    l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, loss = loss.idx, normal = normal.idx), cutval = c(l.cut, g.cut))
    seg.col <- list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.loss = "darkred", normal = "black")
    # l2r.chr <- if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", as.character(my.ascat.seg$SNPpos$chrs))]) else unlist(cs$chrom2chr[as.character(my.ascat.seg$SNPpos$chrs)])
    l2r.chr <- unlist(cs$chrom2chr[as.character(my.ascat.seg$SNPpos$chrs)])
    l2r.value <- data.frame(Chr = l2r.chr,
                            Start = my.ascat.seg$SNPpos$pos,
                            End = my.ascat.seg$SNPpos$pos,
                            Value = my.ascat.seg$Tumor_LogR_wins[,1],
                            stringsAsFactors = FALSE)
    genome.pkg <- data$meta$basic$genome.pkg
    
    # png(paste0(out.dir, "/", samplename, ".L2R.G.png"), width = 1850, height = 980)
    png(paste0(samplename, ".L2R.G.png"), width = 1850, height = 980)
    par(mar = c(1, 6, 3, 1))
    EaCoN.l2rplot.geno(l2r = l2r.value,
                       seg = l2r.seg.obj,
                       seg.col = seg.col,
                       seg.type = "block",
                       seg.normal = TRUE,
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " L2R"),
                       ylim = c(-1.5,1.5))
    dev.off()
    
    ### Genomic BAF plot
    tmsg("Plotting BAF ...")
    # baf.chr <-if(length(grep(pattern = "chr", x = names(cs$chrom2chr), ignore.case = TRUE)) > 0) unlist(cs$chrom2chr[paste0("chr", as.character(my.ascat.seg$SNPpos$chrs))]) else unlist(cs$chrom2chr[as.character(my.ascat.seg$SNPpos$chrs)])
    baf.chr <- unlist(cs$chrom2chr[as.character(my.ascat.seg$SNPpos$chrs)])
    baf.value <- data.frame(Chr = baf.chr,
                            Start = my.ascat.seg$SNPpos$pos,
                            End = my.ascat.seg$SNPpos$pos,
                            Value = my.ascat.seg$Tumor_BAF[,1],
                            stringsAsFactors = FALSE)
    
    # png(paste0(out.dir, "/", samplename, ".BAF.png"), width = 1850, height = 980)
    png(paste0(samplename, ".BAF.png"), width = 1850, height = 980)
    par(mar = c(1, 6, 3, 1))
    EaCoN.bafplot.geno(baf = baf.value,
                       seg = baf.seg,
                       seg.type = "both",
                       genome.pkg = genome.pkg,
                       title = paste0(samplename, " BAF")
    )
    dev.off()
    
    ### L2R Karyotypic plot
    tmsg("Plotting L2R (karyotype) ...")
    # karyoplotf <- paste0(out.dir, "/", samplename, ".L2R.K.png")
    karyoplotf <- paste0(samplename, ".L2R.K.png")
    png(karyoplotf, width = 1850, height = 980)
    EaCoN.l2rplot.karyo(l2r = l2r.value, seg = l2r.seg.obj, seg.col = seg.col, seg.type = "block", seg.normal = TRUE, ylim = c(-1.5,1.5), genome.pkg = genome.pkg)
    dev.off()
    karyoplotf <- tools::file_path_as_absolute(karyoplotf)
    
    ### Chromosomal plots
    tmsg("Plotting L2R (chromosomes) ...")
    # kpdir <- paste0(out.dir, "/chromosomes/")
    kpdir <- "chromosomes/"
    dir.create(kpdir)
    krN <- unique(l2r.seg.obj$pos$Chr)
    krA <- unlist(cs$chr2chrom[krN])
    for (kr in 1:length(krN)) {
      
      png(paste0(kpdir, "/", krA[kr], ".png"), width=1350, height = 750)
      EaCoN.l2rplot.chromo(chr = krN[kr],
                           l2r = l2r.value,
                           l2r.seg = l2r.seg.obj,
                           baf = baf.value,
                           baf.seg = baf.seg,
                           l2r.seg.type = "block",
                           baf.seg.type = "both",
                           seg.normal = TRUE,
                           genome.pkg = genome.pkg
      )
      dev.off()
    }
    # kplotlist <- fpaav(list.files(path = "chromosomes", pattern = "\\.png", full.names = TRUE, ignore.case = FALSE, include.dirs = FALSE, recursive = FALSE))
    kplotlist <- fpaav(vapply(krA, function(k) { paste0("chromosomes/", k, ".png") }, "a"))
    
    ## Instability metrics
    tmsg("Building instability metrics table ...")
    gi.df <- data.frame(Status = c("Gain", "Loss", "Aberrant (Gain + Loss)", "Normal", "TOTAL"), stringsAsFactors = FALSE)
    gi.idx <- list(gi.gain.idx = gain.idx,
                   gi.loss.idx = loss.idx,
                   gi.aberrant.idx = sort(c(gain.idx, loss.idx)),
                   gi.normal.idx = normal.idx,
                   gi.all.idx = 1:nrow(cbsNC.df))
    gi.df$Segments <- vapply(gi.idx, length, 1)
    gi.df$SumWidth <- vapply(names(gi.idx), function(x) { if(length(gi.idx[[x]]) > 0) return(sum(cbsNC.df$End[gi.idx[[x]]] - cbsNC.df$Start[gi.idx[[x]]] +1)) else return(0) }, 1)
    gi.df$Frac <- gi.df$SumWidth / cs$genome.length
    gi.df$MedianWidth <- vapply(names(gi.idx), function(x) { if(length(gi.idx[[x]]) > 0) return(median(cbsNC.df$End[gi.idx[[x]]] - cbsNC.df$Start[gi.idx[[x]]] +1)) else return(0)}, 1)
    gi.df$MedianL2R <- round(vapply(names(gi.idx), function(x) { if(length(gi.idx[[x]]) > 0) return(median(cbsNC.df$Log2Ratio[gi.idx[[x]]])) else return(NA)}, 1), digits = 2)
    
    gi.df$SumWidth <- format(gi.df$SumWidth, big.mark = ",")
    gi.df$Frac <- paste0(format(gi.df$Frac * 100, digits = 2), " %")
    gi.df$MedianWidth <- format(ceiling(gi.df$MedianWidth), big.mark = ",")
    colnames(gi.df) <- c("CNA Status", "# Segments", "Total Width (nt)", "Genome Fraction", "Median Width (nt)", "Median log2(ratio)")
    # write.table(gi.df, file = paste0(out.dir, "/", samplename, ".Instab.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(gi.df, file = paste0(samplename, ".Instab.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    ## Segments table
    ### L2R
    tmsg("Polishing L2R segments table ...")
    minbold <- 1
    segtab.df <- l2r.segments
    segtab.df$Ratio <- 2^segtab.df$Value
    # segtab.df$ChromFull <- paste0("chr", segtab.df$Chrom)
    # segtab.df$ChromFull <- segtab.df$Chrom
    segtab.df$Width <- segtab.df$End - segtab.df$Start +1
    segtab.df$Start.band <- vapply(1:nrow(segtab.df), function(x) {
      scb <- cs$cytobands$chrN == segtab.df$Chr[x] & cs$cytobands$start <= segtab.df$Start[x] & cs$cytobands$end >= segtab.df$Start[x]
      return(paste0(cs$cytobands$chrA[scb], cs$cytobands$cytoband[scb]))
    }, "a")
    segtab.df$End.band <- vapply(1:nrow(segtab.df), function(x) {
      ecb <- cs$cytobands$chrN == segtab.df$Chr[x] & cs$cytobands$start <= segtab.df$End[x] & cs$cytobands$end >= segtab.df$End[x]
      return(paste0(cs$cytobands$chrA[ecb], cs$cytobands$cytoband[ecb]))
    }, "a")
    segtab.df$Status <- "Normal"
    segtab.df$Status[gain.idx] <- "Gain"
    segtab.df$Status[loss.idx] <- "Loss"
    segtab.df$Symbol <- cbs.df$Symbol
    segtab.df$NrSymbol <- vapply(1:nrow(segtab.df), function(x) { length(unlist(strsplit(x = segtab.df$Symbol[x], split = ","))) }, 1)
    segtab.df$l2rbold <- 0
    segtab.df$l2rbold[abs(segtab.df$Value) > minbold] <- 1
    # segtab.df <- segtab.df[,c(10,3,4,11:14,5,9,6,16,15,17)]
    
    # segtab.df <- segtab.df[,c(1,3,4,10:13,5,9,6,15,14,16)]
    segtab.df <- segtab.df[,c(2,3,4,10:13,5,9,6,15,14,16)]
    
    colnames(segtab.df) <- c("Chr", "Start", "End", "Width", "Start Cytoband", "End Cytoband", "Status", "L2R", "Ratio", "# Markers", "# Genes", "symbols", "l2rbold")
    # segtab.df$Start <- format(segtab.df$Start, big.mark = ",")
    # segtab.df$End <- format(segtab.df$End, big.mark = ",")
    # segtab.df$Width <- width.unit.conv(segtab.df$Width)
    segtab.df$L2R <- round(segtab.df$L2R, digits = 3)
    segtab.df$Ratio <- round(segtab.df$Ratio, digits = 3)
    
    l2r.segstat <- c("Gain", "Loss", "Normal")
    l2r.segcol <- c(seg.col$gain, seg.col$loss, seg.col$normal)
    
    ### BAF
    tmsg("Polishing BAF segmentation table ...")
    segbaf.df <- baf.seg
    # segbaf.df$Chrom <- paste0("chr", segbaf.df$chrA)
    segbaf.df$Chrom <- segbaf.df$chrA
    segbaf.df$Start.band <- vapply(1:nrow(segbaf.df), function(x) {
      scb <- cs$cytobands$chrN == segbaf.df$Chr[x] & cs$cytobands$start <= segbaf.df$Start[x] & cs$cytobands$end >= segbaf.df$Start[x]
      return(paste0(cs$cytobands$chrA[scb], cs$cytobands$cytoband[scb]))
    }, "a")
    segbaf.df$End.band <- vapply(1:nrow(segbaf.df), function(x) {
      ecb <- cs$cytobands$chrN == segbaf.df$Chr[x] & cs$cytobands$start <= segbaf.df$End[x] & cs$cytobands$end >= segbaf.df$End[x]
      return(paste0(cs$cytobands$chrA[ecb], cs$cytobands$cytoband[ecb]))
    }, "a")
    segbaf.df <- foreach::foreach(seg = 1:nrow(segbaf.df), .combine = "rbind") %do% {
      ingenz <- gen.df$symbol[gen.df$chrN == segbaf.df$Chr[seg] & gen.df$start <= segbaf.df$End[seg] & gen.df$end >= segbaf.df$Start[seg]]
      return(data.frame(segbaf.df[seg, ], Genes = length(ingenz), Symbol = paste0(ingenz, collapse = ","), stringsAsFactors = FALSE))
    }
    segbaf.df <- segbaf.df[,c(13,8:10,14,15,12,11,3,16,17)]
    colnames(segbaf.df) <- c("Chr", "Start", "End", "Width", "Start Cytoband", "End Cytoband", "Status", "BAF", "# Markers", "# Genes", "symbols")
    # segbaf.df$Start <- format(segbaf.df$Start, big.mark = ",")
    # segbaf.df$End <- format(segbaf.df$End, big.mark = ",")
    # segbaf.df$Width <- format(segbaf.df$Width, big.mark = ",")
    # segbaf.df$Width <- width.unit.conv(segbaf.df$Width)
    segbaf.df$BAF <- round(segbaf.df$BAF, digits = 3)
    
    baf.segstat <- c("Hetero", "Homo", "Unbalanced")
    baf.segcol <- c(paste0("rgb(", paste0(col2rgb("black"), collapse = ","), ")"), paste0("rgb(", paste0(col2rgb("cadetblue4"), collapse = ","), ")"), paste0("rgb(", paste0(col2rgb("coral1"), collapse = ","), ")"))
    
    
    ## Targets table
    tmsg("Polishing targets table ...")
    # targ.regz[["Gene Start"]] <- format(targ.regz[["Gene Start"]], big.mark = ",")
    # targ.regz[["Gene End"]] <- format(targ.regz[["Gene End"]], big.mark = ",")
    # targ.regz[["Gene Width"]] <- width.unit.conv(targ.regz[["Gene Width"]])
    # targ.regz[["Match Start"]] <- format(targ.regz[["Match Start"]], big.mark = ",")
    # targ.regz[["Match End"]] <- format(targ.regz[["Match End"]], big.mark = ",")
    # targ.regz[["Match Width"]] <- width.unit.conv(targ.regz[["Match Width"]])
    # targ.regz[["L2R Segment Width"]] <- width.unit.conv(targ.regz[["L2R Segment Width"]])
    # targ.regz[["BAF Segment Width"]] <- width.unit.conv(targ.regz[["BAF Segment Width"]])
    targ.regz[["L2R Value"]] <- round(targ.regz[["L2R Value"]], digits = 3)
    targ.regz[["BAF Value"]] <- round(targ.regz[["BAF Value"]], digits = 3)
    
    targ.regz$l2rbold <- 0
    targ.regz$l2rbold[abs(targ.regz[["L2R Value"]]) > minbold] <- 1
    
    ## Truncated table
    if (length(gen.trunk.idx) > 0) {
      
      # trunc.regz[["Gene Start"]] <- format(trunc.regz[["Gene Start"]], big.mark = ",")
      # trunc.regz[["Gene End"]] <- format(trunc.regz[["Gene End"]], big.mark = ",")
      # trunc.regz[["Gene Width"]] <- width.unit.conv(trunc.regz[["Gene Width"]])
      # trunc.regz[["Match Start"]] <- format(trunc.regz[["Match Start"]], big.mark = ",")
      # trunc.regz[["Match End"]] <- format(trunc.regz[["Match End"]], big.mark = ",")
      # trunc.regz[["Match Width"]] <- width.unit.conv(trunc.regz[["Match Width"]])
      # trunc.regz[["L2R Segment Width"]] <- width.unit.conv(trunc.regz[["L2R Segment Width"]])
      # trunc.regz[["BAF Segment Width"]] <- width.unit.conv(trunc.regz[["BAF Segment Width"]])
      trunc.regz[["L2R Value"]] <- round(trunc.regz[["L2R Value"]], digits = 3)
      trunc.regz[["BAF Value"]] <- round(trunc.regz[["BAF Value"]], digits = 3)
      
      trunc.regz$l2rbold <- 0
      trunc.regz$l2rbold[abs(trunc.regz[["L2R Value"]]) > minbold] <- 1
    }
    
    ## Deactivate SOLO if no abnormal segment available
    if (solo & (length(which(segtab.df$Status != "Normal")) == 0)) {
      solo <- FALSE
      tmsg("No aberrant segment identified, deactivating SOLO ...")
    }
    
    ## Render
    message(tmsg("Rendering report ..."))
    # htmlout <- paste0(out.dir, "/", samplename, ".REPORT.html")
    # gplotlist <- fpaav(paste0(samplename, ".", c("ASPCF", "L2R.G", "BAF", "PredictedGermline"), ".png"))
    gplotlist <- fpaav(paste0(samplename, ".", c(paste0("SEG.", data$meta$eacon$segmenter), "L2R.G", "BAF"), ".png"))
    htmlout <- paste0(samplename, ".REPORT.html")
    show.flag <- if ((data$meta$basic$source == "microarray") & (data$meta$basic$manufacturer == "Affymetrix")) TRUE else FALSE
    rmarkdown::render(input = rmd.path, output_format = "html_document", output_file = htmlout, output_dir = getwd(), intermediates_dir = getwd())
    
    ## Hacking HTML
    message(tmsg("Hacking HTML ..."))
    htmltmp <- readLines(htmlout)
    ### Page width
    ltcidx <- grep(pattern = ".main-container \\{", x = htmltmp)
    for (x in ltcidx) htmltmp[x+1] <- "  max-width: 2000px;"
    ### Page margin
    mtocidx1 <- grep(pattern = "<div class=\"col-xs-12 col-sm-4 col-md-3\">", x = htmltmp)
    if(length(mtocidx1) > 0) htmltmp[mtocidx1] <- "<div class=\"col-xs-12 col-sm-3 col-md-2\">"
    mtocidx2 <- grep(pattern = "<div class=\"toc-content col-xs-12 col-sm-8 col-md-9\">", x = htmltmp)
    if(length(mtocidx2) > 0) htmltmp[mtocidx2] <- "<div class=\"toc-content col-xs-12 col-sm-9 col-md-10\">"
    writeLines(htmltmp, htmlout)
    tmsg("Done.")
    
  }
  
  if (solo) {
    tmsg("Making SOLO ...")
    # tmpcbsname <- paste0(out.dir, "/", samplename, "_tmp.cbs")
    tmpcbsname <- paste0(samplename, "_solo.cbs")
    write.table(data$cbs$cut, tmpcbsname, sep="\t", quote = FALSE, row.names = FALSE)
    # EaCoN.Annotate.solo(out.dir = out.dir, samplename = samplename, genome = genome, ldb = ldb)
    Annotate.solo(cbs.file = tmpcbsname, genome = genome, ldb = ldb)
    file.remove(tmpcbsname)
  }
  setwd(oridir)
  tmsg("Done.")
}

## Run EaCoN.Annotate(), from a file
Annotate.ff <- function (RDS.file = NULL, ...) {
  if (is.null(RDS.file)) stop(tmsg("A RDS file is needed !"), call. = FALSE)
  if (!file.exists(RDS.file)) stop(tmsg(paste0("Could not find RDS file ", RDS.file, " !")), call. = FALSE)
  ## Data loading
  tmsg(paste0("Loading data from ", RDS.file, " ..."))
  my.data <- readRDS(RDS.file)
  Annotate(data = my.data, out.dir = dirname(RDS.file), ...)
}

## Run EaCoN.Annotate() in batch mode
Annotate.ff.Batch <- function(RDS.files = list.files(path = getwd(), pattern = "\\.SEG\\..*\\.RDS$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE), nthread = 1, cluster.type = "PSOCK", ...) {

  if (length(RDS.files) == 0) stop("A list of RDS files is required !", call. = FALSE)
  message("Running Annotate.ff() in batch mode ...")
  message(paste0(" Found ", length(RDS.files), " files to process."))
  current.bitmapType <- getOption("bitmapType")
  if (length(RDS.files) < nthread) nthread <- length(RDS.files)
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  s <- ""
  targ.all <- foreach::foreach(s = 1:length(RDS.files), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    EaCoN.set.bitmapType(type = current.bitmapType)
    Annotate.ff(RDS.file = RDS.files[s], ...)
  }
  parallel::stopCluster(cl)
}

Annotate.solo <- function(cbs.file = NULL, genome = NULL, ldb = "/mnt/data_cigogne/bioinfo/") {
  if (is.null(genome)) stop(tmsg("A genome is required !"), call. = FALSE)
  if (!file.exists(cbs.file)) stop(tmsg(paste0("Could not find CBS file [", cbs.file, "] !")), call. = FALSE)
  cbs.df <- read.table.fast(cbs.file)
  self.pkg.name <- "EaCoN"
  data(list = genome, package = self.pkg.name, envir = environment())
  ## Converting CBS to SOLO
  message(tmsg(" Converting to SOLO ..."))
  solo.df <- data.frame(Loc = paste0(unlist(cs$chr2chrom[cbs.df$Chr]), ":", cbs.df$Start, "-", cbs.df$End), Probes = cbs.df$Probes, Status = sign(cbs.df$Log2Ratio), L2R = cbs.df$Log2Ratio, Ratio = 2^cbs.df$Log2Ratio, stringsAsFactors = FALSE)
  solo.df <- solo.df[!solo.df$Status == 0,]
  solo.df$Status[solo.df$Status == 1] <- "G"
  solo.df$Status[solo.df$Status == -1] <- "L"
  solofile <- sub(pattern = paste0(tools::file_ext(cbs.file), "$"), replacement = "solo", x = cbs.file, ignore.case = TRUE)
  write.table(solo.df, file = solofile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  gb <- as.numeric(sub(pattern = "[a-z]+", replacement = "", x = genome, ignore.case = TRUE))
  sp <- as.character(sub(pattern = "[0-9]+", replacement = "", x = genome, ignore.case = TRUE))
  if(sp == "hg") sp <- "hs"
  suppressMessages(try(system(paste0("grd --sp ", sp, " --gb ", gb, " --ldb ", ldb, ' -o ', dirname(cbs.file), " -m solo ", solofile), wait = TRUE)))
}


