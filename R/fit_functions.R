## L2R tf fit loop function
l2r.fitloop <- function(l2rObj, tfd, smo = 400, method = "loess") {

  ### FITLOOP
  minitf <- tfd[,-c(1:4), drop = FALSE]
  tfheads <- colnames(minitf)
  b <- ncol(minitf)+1
  posfit <- c()

  message(tmsg(paste0("Init (", l2rObj$rm.mad, ")")))

  while ( (b != 1) & (ncol(minitf) != 0) ) {

    biggy <- list()
    biggy <- append(biggy, list(l2rObj))
    rmtest <- l2rObj$rm.mad
    for (z in 1:length(minitf)) {
      if (method == "loess") tempfit <- l2r.fit(biggy[[1]]$l2r, minitf[,z], smo)
      if (method == "pcs") tempfit <- l2r.pcs(biggy[[1]]$l2r, minitf[,z], smo)
      biggy <- append(biggy, list(tempfit))
      rmtest <- c(rmtest, tempfit$rm.mad)
      # message(paste0(z, " / ", tempfit$rm.mad))
    }
    b <- which.min(rmtest)
    if (b > 1) {
      message(tmsg(paste0(" Positive fit with ", tfheads[b-1], " (", min(rmtest), ")")))
      posfit <- c(posfit, tfheads[b-1])
      l2rObj <- biggy[[b]]
      minitf <- as.data.frame(minitf[,-c(b-1)])
      tfheads <- tfheads[-c(b-1)]
    }
  }
  return(list(l2r = l2rObj, pos = posfit))
}

## FONCTION LOESSFIT
l2r.fit <- function(l2r, tf, smo) {
  l2fN <- limma::loessFit(l2r, tf)
  l2N <- l2r-l2fN$fitted
  Nrmspread <- sum(abs(diff(as.numeric(runmed(l2N[!is.na(l2N)], smo)))))
  return(list(l2r=l2N, rm.mad=Nrmspread))
}

## FONCTION PERCENTILE SCALING
l2r.pcs <- function(l2r, tf, smo) {
  l2N <- GCnorm.pcs(l2r, tf)
  Nrmspread <- sum(abs(diff(as.numeric(runmed(l2N[!is.na(l2N)], smo)))))
  return(list(l2r=l2N, rm.mad=Nrmspread))
}

GCnorm.pcs <- function(measures=NULL, gc=NULL) {
  # message("Performing GC PC-scaling ...")
  gcpc.idx <- sapply(sort(unique(gc)), function(x) { return(which(gc == x)) })
  names(gcpc.idx) <- sort(unique(gc))
  gcpc.value <- sapply(1:length(gcpc.idx), function(x) { return(measures[gcpc.idx[[x]]]) })
  gcpc.med <- sapply(gcpc.value, median)
  measures.new <- rep(NA, length(gc))
  for (x in 1:length(gcpc.med))  measures.new[gcpc.idx[[x]]] <- measures[gcpc.idx[[x]]] - gcpc.med[x]
  return(measures.new)
}


