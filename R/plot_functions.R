EaCoN.l2rplot.geno <- function(l2r = NULL, seg = NULL, seg.col = list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.red = "darkred"),
                               seg.type = "block", seg.normal = TRUE, ylim = c(-1.5,1.5), genome.pkg = NULL, title = "L2RPLOT") {

  # message(tmsg(paste0("Loading ", genome.pkg, " ...")))
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)
  
  # data(list = genome, package = "chromosomes", envir = environment())
  l2r <- l2r[!is.na(l2r$Value),]
  l2r$Start.geno <- l2r$Start + cs$chromosomes$chr.length.toadd[l2r$Chr]
  l2r$End.geno <- l2r$End + cs$chromosomes$chr.length.toadd[l2r$Chr]

  seg$pos$Start.geno <- seg$pos$Start + cs$chromosomes$chr.length.toadd[seg$pos$Chr]
  seg$pos$End.geno <- seg$pos$End + cs$chromosomes$chr.length.toadd[seg$pos$Chr]
  graphics::plot(l2r$Start.geno, l2r$Value,
                 pch = ".", cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,cs$genome.length), ylim = ylim,
                 main = title, cex.main = 2, ylab = "Log2(ratio)",
                 cex.lab = 2, col = "grey80", xaxt = "n")
  ink <- cs$chromosomes$chrN %in% l2r$Chr
  yrange = abs(diff(ylim))
  m.pos <- c(ylim[2] - (sign(ylim[2]) * yrange * .05), ylim[1] + (sign(ylim[2]) * yrange * .05))
  m.mod <- -(cs$chromosomes$chrN[ink] %% 2) +2
  try(text(x = cs$chromosomes$mid.chr.geno[ink], y = m.pos[m.mod], labels = cs$chromosomes$chrom[ink], cex = 1))
  abline(h = 0, col = 1, lwd = 2, lty = 3)
  smo <- round(nrow(l2r)/200)
  if (smo%%2 == 0) smo <- smo + 1
  l2r.rm <- stats::runmed(l2r$Value, smo)

  for (k in unique(l2r$Chr)) {
    kplot <- l2r$Chr == k
    lines(l2r$Start.geno[kplot], l2r.rm[kplot], col="grey25", lwd=1)
  }
  
  if(seg.type %in% c("block", "both")) {
    if (seg.normal & length(seg$idx[["normal"]]) > 0) rect(seg$pos$Start.geno[seg$idx[["normal"]]], 0, seg$pos$End.geno[seg$idx[["normal"]]], seg$pos$Value[seg$idx[["normal"]]], lwd = 1, col=adjustcolor(seg.col[["normal"]], alpha.f = .3), border = seg.col[["normal"]])
    for (stype in c("gain", "loss")) {
      if (length(seg$idx[[stype]]) > 0) rect(seg$pos$Start.geno[seg$idx[[stype]]], 0, seg$pos$End.geno[seg$idx[[stype]]], seg$pos$Value[seg$idx[[stype]]], lwd = 1, col=adjustcolor(seg.col[[stype]], alpha.f = .5), border = seg.col[[stype]])
    }
  }
  
  if(seg.type %in% c("line", "both")) {
    seg$idx$outscale.gain <- which(seg$pos$Value >= ylim[2])
    seg$idx$outscale.loss <- which(seg$pos$Value <= ylim[1])
    if (seg.normal & length(seg$idx[["normal"]]) > 0) segments(seg$pos$Start.geno[seg$idx[["normal"]]], seg$pos$Value[seg$idx[["normal"]]], seg$pos$End.geno[seg$idx[["normal"]]], seg$pos$Value[seg$idx[["normal"]]], col = seg.col[["normal"]], lwd = 5)
    for (stype in names(seg$idx)[names(seg$idx) != "normal"]) {
      if (length(seg$idx[[stype]]) > 0) segments(seg$pos$Start.geno[seg$idx[[stype]]], seg$pos$Value[seg$idx[[stype]]], seg$pos$End.geno[seg$idx[[stype]]], seg$pos$Value[seg$idx[[stype]]], col = seg.col[[stype]], lwd = 5)
    }
  }
  abline(h = seg$cutval, lty = 2, col = c(seg.col[["loss"]], seg.col[["gain"]]))
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 2, lwd = 2)
}

EaCoN.bafplot.geno <- function(baf = NULL, seg = NULL, seg.col = list(Hetero = "black", Homo = "cadetblue4", Unbalanced = "coral1"),
                               seg.type = "line", ylim = c(-.01,1.01), genome.pkg = NULL, title = "BAFPLOT") {

  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  cs <- chromobjector(BSg.obj)
  
  # data(list = genome, package = "chromosomes", envir = environment())
  baf$Start.geno <- baf$Start + cs$chromosomes$chr.length.toadd[baf$Chr]
  baf$End.geno <- baf$End + cs$chromosomes$chr.length.toadd[baf$Chr]

  seg$Start.geno <- seg$Start + cs$chromosomes$chr.length.toadd[seg$Chr]
  seg$End.geno <- seg$End + cs$chromosomes$chr.length.toadd[seg$Chr]
  graphics::plot(baf$Start.geno, baf$Value,
                 pch = ".", cex = 2, xaxs = "i", yaxs = "i", xlim = c(0,cs$genome.length), ylim = ylim,
                 main = title, cex.main = 2, ylab = "B-Allele Frequency",
                 cex.lab = 2, col = "grey80", xaxt = "n")
  ink <- cs$chromosomes$chrN %in% baf$Chr
  yrange = abs(diff(ylim))
  m.pos <- c(ylim[2] - (sign(ylim[2]) * yrange * .075), ylim[1] + (sign(ylim[2]) * yrange * .075))
  m.mod <- -(cs$chromosomes$chrN[ink] %% 2) +2
  try(text(x = cs$chromosomes$mid.chr.geno[ink], y = m.pos[m.mod], labels = cs$chromosomes$chrom[ink], cex = 1))
  abline(h = 0.5, col = 1, lwd = 2, lty = 3)

  if(seg.type %in% c("block", "both")) {
    rect(seg$Start.geno, 0.5, seg$End.geno, seg$Value, lwd = 1, col=adjustcolor(unlist(seg.col[seg$Status]), alpha.f = .2), border = NA)
    rect(seg$Start.geno, 0.5, seg$End.geno, -seg$Value + 1, lwd = 1, col=adjustcolor(unlist(seg.col[seg$Status]), alpha.f = .2), border = NA)
  }
  if(seg.type %in% c("line", "both")) {
    segments(seg$Start.geno, seg$Value, seg$End.geno, seg$Value, col = unlist(seg.col[seg$Status]), lwd = 7)
    segments(seg$Start.geno, 1 - seg$Value, seg$End.geno, 1 - seg$Value, col = unlist(seg.col[seg$Status]), lwd = 7)
  }
  # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, lty = 2, lwd = 2)
  abline(v = cs$chromosomes$chr.length.sum, col = 1, lty = 2, lwd = 2)
}

EaCoN.l2rplot.karyo <- function(l2r = NULL, seg = NULL, seg.col = list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.red = "darkred"),
                                seg.type = "block", seg.normal = TRUE, ylim = c(-1.5,1.5), genome.pkg = NULL) {

  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  # cs <- chromobjector(BSg.obj)
  
  self.pkg.name <- "EaCoN"
  data(list = genome, package = self.pkg.name, envir = environment())
  l2r <- l2r[!is.na(l2r$Value),]

  par(mgp=c(0,0,0), omi=c(0,0.2,0.1,0), xaxt="n", yaxt="n", bty="n")
  nk <- length(cs$chromosomes$chrom[!cs$chromosomes$chrom %in% c("chrM", "M")])
  bigk1 <- which.max(cs$chromosomes$chr.length[1:ceiling(nk/2)])
  bigk2 <- ceiling(nk/2) + which.max(cs$chromosomes$chr.length[ceiling(nk/2)+1:nk])
  zone <- matrix(seq(1,ceiling(nk/2)*4,1), nrow=2, ncol=nk, byrow=T)
  layW = rep(c(0.01, .031667), 12)
  layH = c(0.6,0.4)
  chridx <- 1:nk
  layout(zone, widths=layW, heights = layH)
  smo <- round(nrow(l2r)/200)
  if(smo %% 2 == 0) smo <- smo + 1
  l2r.rm <- stats::runmed(l2r$Value, smo)

  for (k in chridx) {

    kx <- cs$chromosomes$chrA[k]
    kylim <- ifelse(k<=ceiling(nk/2), -cs$chromosomes$chr.length[bigk1], -cs$chromosomes$chr.length[bigk2])
    par(mar=c(1,1,1,0))
    graphics::plot(0, 0, xlim=c(0,1), ylim=c(kylim, 0), type="n", xlab="", ylab="")
    cytochr <- which(cs$cytobands$chrN == k)
    rect(cs$cytobands$x1[cytochr], -cs$cytobands$start[cytochr], cs$cytobands$x2[cytochr], -cs$cytobands$end[cytochr], col=cs$cytobands$giestaincol[cytochr])

    kinseg <- which(seg$pos$Chr == k)
    sli <- intersect(seg$idx$loss, kinseg)
    sgi <- intersect(seg$idx$gain, kinseg)
    sni <- intersect(seg$idx$normal, kinseg)

    if (length(sli) > 0) segments(0, -seg$pos$Start[sli], 0, -seg$pos$End[sli], col=seg.col$loss, lwd=4)
    if (length(sgi) > 0) segments(1, -seg$pos$Start[sgi], 1, -seg$pos$End[sgi], col=seg.col$gain, lwd=4)

    kinprob <- which(l2r$Chr == k)
    par(mar=c(1,0,1,1))
    graphics::plot(0, 0, xlim=ylim, ylim=c(kylim, 0), type="n", xlab="")
    points(l2r$Value[kinprob], -l2r$Start[kinprob], cex=0.2, col="grey80")
    lines(l2r.rm[kinprob], -l2r$Start[kinprob], col="grey25")
    if (length(sli) > 0) rect(0, -seg$pos$Start[sli], seg$pos$Value[sli], -seg$pos$End[sli], col=adjustcolor(seg.col$loss, alpha.f = .5), border = adjustcolor(seg.col$loss, alpha.f = .5))
    if (length(sgi) > 0) rect(0, -seg$pos$Start[sgi], seg$pos$Value[sgi], -seg$pos$End[sgi], col=adjustcolor(seg.col$gain, alpha.f = .5), border = adjustcolor(seg.col$gain, alpha.f = .5))
    if (seg.normal) if (length(sni) > 0) rect(0, -seg$pos$Start[sni], seg$pos$Value[sni], -seg$pos$End[sni], col=adjustcolor(seg.col$normal, alpha.f = .30), border = adjustcolor(seg.col$normal, alpha.f = .5))
    abline(v=0, col="grey50", lty=2)
    try(text(-.8, -1e+06, kx, pos=3, cex=2))
  }
}

EaCoN.l2rplot.chromo <- function(chr = NULL, l2r = NULL, l2r.seg = NULL, baf = NULL, baf.seg = NULL,
                                 l2r.seg.col = list(gain = "blue", outscale.gain = "midnightblue", loss = "red", outscale.red = "darkred"),
                                 l2r.seg.type = "block", baf.seg.col = list(Hetero = "black", Homo = "cadetblue4", Unbalanced = "coral1"),
                                 baf.seg.type = "both", seg.normal = TRUE, genome.pkg = NULL, l2r.ylim = c(-1.5,1.5), baf.ylim = c(-.01,1.01)) {

  if (is.null(chr)) return()
  
  suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
  BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
  genome <- BSgenome::providerVersion(BSg.obj)
  # cs <- chromobjector(BSg.obj)
  
  self.pkg.name <- "EaCoN"
  data(list = genome, package = self.pkg.name, envir = environment())

  l2r <- l2r[!is.na(l2r$Value),]
  kLprob.idx <- which(l2r$Chr == chr)
  kLseg.idx <- which(l2r.seg$pos$Chr == chr)
  kBprob.idx <- which(baf$Chr == chr)
  kBseg.idx <- which(baf.seg$Chr == chr)

  bigkn <- max(cs$chromosomes$chr.length)

  zone <- matrix(1:3, ncol=1, byrow=T)
  layH = c(.47,0.06,.47)
  par(mgp=c(0,0,0), omi=c(0.1,0.2,0.1,0.1), mar = c(0,1,1,1), xaxt="n", xaxs = "i", yaxs = "i", bty = "n")
  layout(zone, heights = layH)

  graphics::plot(l2r$Start[kLprob.idx], l2r$Value[kLprob.idx], pch = ".", cex = 4, col = "grey80", ylim = l2r.ylim, xlim=c(0, bigkn), xlab = "", ylab = "", cex.axis = 2)
  abline(h = 0, col = 1, lty = 3)
  smo <- round(nrow(l2r)/200)
  if(smo %% 2 == 0) smo <- smo + 1
  l2r.rm <- stats::runmed(l2r$Value, smo)
  lines(l2r$Start[kLprob.idx], l2r.rm[kLprob.idx], col = "grey25")
  if(l2r.seg.type %in% c("block", "both")) {
    if (seg.normal) {
      normidx <- intersect(l2r.seg$idx[["normal"]],kLseg.idx)
      if(length(normidx) > 0) rect(l2r.seg$pos$Start[normidx], 0, l2r.seg$pos$End[normidx], l2r.seg$pos$Value[normidx], lwd = 1, col=adjustcolor(l2r.seg.col[["normal"]], alpha.f = .3), border = l2r.seg.col[["normal"]])
    }
    for (stype in c("gain", "loss")) {
      stidx <- intersect(l2r.seg$idx[[stype]],kLseg.idx)
      if (length(stidx) > 0) rect(l2r.seg$pos$Start[stidx], 0, l2r.seg$pos$End[stidx], l2r.seg$pos$Value[stidx], lwd = 1, col=adjustcolor(l2r.seg.col[[stype]], alpha.f = .3), border = l2r.seg.col[[stype]])
    }
  }
  if(l2r.seg.type %in% c("line", "both")) {
    if (seg.normal) {
      if (length(normidx) > 0) segments(l2r.seg$pos$Start[normidx], l2r.seg$pos$Value[normidx], l2r.seg$pos$End[normidx], l2r.seg$pos$Value[normidx], col=l2r.seg.col[["normal"]], lwd = 5)
    }
    for (stype in names(seg$idx)[names(seg$idx) != "normal"]) {
      stidx <- intersect(l2r.seg$idx[[stype]],kLseg.idx)
      if (length(stidx) > 0) segments(l2r.seg$pos$Start[stidx], l2r.seg$pos$Value[stidx], l2r.seg$pos$End[stidx], l2r.seg$pos$Value[stidx], col=l2r.seg.col[[stype]], lwd = 5)
    }
  }
  abline(h = l2r.seg$cutval, lty = 2, col = c(l2r.seg.col[["loss"]], l2r.seg.col[["gain"]]))

  cytochr <- which(cs$cytobands$chrN == chr)
  graphics::plot(0, 0, ylim=c(.15,.85), xlim=c(0, bigkn), type="n", xlab="", ylab="", yaxt = "n")
  rect(cs$cytobands$start[cytochr], cs$cytobands$x1[cytochr], cs$cytobands$end[cytochr], cs$cytobands$x2[cytochr], col=cs$cytobands$giestaincol[cytochr])

  graphics::plot(baf$Start[kBprob.idx], baf$Value[kBprob.idx], pch = ".", cex = 4, col = "grey80", ylim = baf.ylim, xlim=c(0, bigkn), xlab = "", ylab = "", cex.axis = 2)
  if (length(kBseg.idx) > 0) {
    if(baf.seg.type %in% c("block", "both")) {
      rect(baf.seg$Start[kBseg.idx], 0.5, baf.seg$End[kBseg.idx], baf.seg$Value[kBseg.idx], lwd = 1, col=adjustcolor(unlist(baf.seg.col[baf.seg$Status[kBseg.idx]]), alpha.f = .2), border = NA)
      rect(baf.seg$Start[kBseg.idx], 0.5, baf.seg$End[kBseg.idx], -baf.seg$Value[kBseg.idx] + 1, lwd = 1, col=adjustcolor(unlist(baf.seg.col[baf.seg$Status[kBseg.idx]]), alpha.f = .2), border = NA)
    }
    if(baf.seg.type %in% c("line", "both")) {
      segments(baf.seg$Start[kBseg.idx], baf.seg$Value[kBseg.idx], baf.seg$End[kBseg.idx], baf.seg$Value[kBseg.idx], lwd = 7, col = unlist(baf.seg.col[baf.seg$Status[kBseg.idx]]))
      segments(baf.seg$Start[kBseg.idx], 1 - baf.seg$Value[kBseg.idx], baf.seg$End[kBseg.idx], 1 - baf.seg$Value[kBseg.idx], lwd = 7, col = unlist(baf.seg.col[baf.seg$Status[kBseg.idx]]))
    }
  }
}

## Rorchard plot (BAF vs L2R)
EaCoN.Rorschard.plot <- function(data = NULL, cnpTotal = NULL) {
  k.sqrt <- ceiling(sqrt(length(data$data$chrs)))
  par(mar = c(1, 1, 1, 1), mfrow = c(k.sqrt, k.sqrt))
  for (k in 1:length(data$data$ch)) {
    graphics::plot(data$data$Tumor_BAF[[1]][data$germline$germlinegenotypes],
                   data$data$Tumor_LogR[[1]][data$germline$germlinegenotypes],
                   pch = ".", cex = 2, xlim = c(0, 1), ylim = c(-2,2), col = "grey95")
    points(data$data$Tumor_BAF[[1]][!data$germline$germlinegenotypes],
           data$data$Tumor_LogR[[1]][!data$germline$germlinegenotypes],
           pch = ".", cex = 2, col = "grey50")
    kin <- data$data$ch[[k]][!(data$data$ch[[k]] %in% which(data$germline$germlinegenotypes))]
    r.col <- if(is.null(cnpTotal)) 3 else cnpTotal[kin] + 1
    points(data$data$Tumor_BAF[[1]][kin], data$data$Tumor_LogR[[1]][kin], pch = ".", cex = 4, col = r.col)
    try(text(x = 0.5, y = 2, labels = data$data$chrs[k], pos = 1, cex = 2))
  }
}
