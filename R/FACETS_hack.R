## Hacked from facets package (0.5.14) to handle any genome which X chromosome can be obtained from the data$data$sexchromosomes.
procSample.hacked <- function(x, cval=150, min.nhet=15, dipLogR=NULL, chromlevels=c(1:22, "X", "Y")) {
  # ensure availability of seg.tree
  if (is.null(x$seg.tree)) stop("seg.tree is not available")
  # get the numeric value of chromosome X
  nX <- x$nX
  # make sure that original cval is smaller than current one
  cval.fit <- attr(x$seg.tree, "cval")
  if (cval.fit > cval) stop("original fit used cval = ", cval.fit)
  # jointseg etc
  jseg <- x$jointseg
  jseg <- jseg[is.finite(jseg$cnlr),]
  # chromosomes with data and their counts
  chrs <- x$chromlevels
  nchr <- length(chrs)
  # get the segment summary for the fit in seg.tree
  nsegs <- 0
  for (i in 1:nchr) {
    seg.widths <- diff(c(0, sort(unique(x$seg.tree[[i]][,3]))))
    jseg$seg[jseg$chrom==chrs[i]] <- nsegs + rep(1:length(seg.widths), seg.widths)
    nsegs <- nsegs + length(seg.widths)
  }
  focalout <- facets:::jointsegsummary(jseg)
  # cnlr.median to the left and right
  cnlr.med.l <- c(0, focalout$cnlr.median[-nsegs])
  cnlr.med.r <- c(focalout$cnlr.median[-1], 0)
  # mad of cnlr noise
  cnlr.mad <- mad(jseg$cnlr - rep(focalout$cnlr.median, focalout$num.mark))
  # segments that show focal changes have big jump in cnlr.median
  focalout$focal <- 1*(focalout$cnlr.median > pmax(cnlr.med.l, cnlr.med.r)+3*cnlr.mad) + 1*(focalout$cnlr.median < pmin(cnlr.med.l, cnlr.med.r)-3*cnlr.mad)
  # get the segments for the specified cval 
  nsegs <- 0
  for (i in 1:nchr) {
    seg.widths <- diff(facets:::prune.cpt.tree(x$seg.tree[[i]], cval))
    jseg$seg[jseg$chrom==chrs[i]] <- nsegs + rep(1:length(seg.widths), seg.widths)
    nsegs <- nsegs + length(seg.widths)
  }
  # adding the focal change segments - need a jump at the beginning and end
  jseg$seg0 <- jseg$seg # detected segments
  # jump at the beginning (twice the height)
  jseg$seg <- jseg$seg + rep(cumsum(2*focalout$focal), focalout$num.mark)
  # drop back for the focal segment to get the steps right
  jseg$seg <- jseg$seg - rep(focalout$focal, focalout$num.mark)
  # focal segment could already be in; so change seg indicator
  jseg$seg <- cumsum(c(1, 1*(diff(jseg$seg) > 0)))
  # segment summaries
  out <- facets:::jointsegsummary(jseg)
  # cluster the segments
  out <- facets:::clustersegs(out, jseg, min.nhet)
  # put in the clustered values for snps
  jseg$segclust[is.finite(jseg$cnlr)] <- rep(out$segclust, out$num.mark)
  # find dipLogR and fit cncf
  if (is.null(dipLogR)) {
    oo <- facets:::findDiploidLogR(out, jseg$cnlr)
  } else {
    oo <- list()
    oo$out0 <- "empty"
    oo$dipLogR <- dipLogR
  }
  out <- facets:::fitcncf(out, oo$dipLogR, nX)
  c(list(jointseg=jseg, out=out, nX=nX, chromlevels=chromlevels), oo[-1])
}