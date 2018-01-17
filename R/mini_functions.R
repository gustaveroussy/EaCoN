celstruc <- function(celdat = NULL) {
  celcols <- celdat$header$cols
  celrows <- celdat$header$rows
  intensities <- matrix(celdat$intensities, nrow=celcols, ncol=celrows)
  int.median <- median(celdat$intensities, na.rm = TRUE)
  n.outliers <- length(celdat$outliers)
  pc.outliers <- n.outliers / prod(dim(intensities))
  rm(celdat)
  return(list(intensities = intensities, int.median = int.median, n.outliers = n.outliers, pc.outliers = pc.outliers))
}

width.unit.conv <- function(coord = NULL, digits = 3) {
  coord.num <- as.numeric(coord)
  coord.na <- is.na(coord.num)
  coord.new <- NA
  if (length(coord.na) > 0) coord.new[coord.na] <-coord[coord.na]
  coord <- coord.num[!coord.na]
  coord.sign <- sign(coord)
  coord <- coord * coord.sign
  unit.power <- seq.int(0,18,3)
  unit.list <- list("1" = "b", "2" = "Kb", "3" = "Mb", "4" = "Gb", "Tb", "Pb")
  coord.log <- log10(coord)
  coord.idx <- vapply(coord.log, function(x) { if (is.na(x)) NA else max(which((x - unit.power) >= 0)) }, 1)
  coord.unit <- unit.list[coord.idx]
  coord.div <- unit.power[coord.idx]
  coord.new.temp <- round(coord*coord.sign / 10^coord.div, digits = digits)
  coord.new.temp <- sapply(1:length(coord.new.temp), function(x) { if (is.na(coord.new.temp[x])) NA else paste0(coord.new.temp[x], " ", coord.unit[x])})
  coord.new[!coord.na] <- coord.new.temp
  return(coord.new)
}

get.mad <- function(val = NULL) {
  ad <- abs(diff(val))
  ad <- ad[ad != 0]
  return(stats::median(ad, na.rm = TRUE))
}

get.valid.genomes <- function() {
  return(list("hg18" = "BSgenome.Hsapiens.UCSC.hg18", "hg19" = "BSgenome.Hsapiens.UCSC.hg19", "hg38" = "BSgenome.Hsapiens.UCSC.hg38", "GRCh37-lite" = "BSgenome.Hsapiens.TCGA.GRCh37.lite"))
}

num2mat <- function(num=NULL) {
  x = ceiling(num**.5)
  y = ceiling(num/x)
  return(c(x, y))
}

getmeta <- function(key = NULL, meta = NULL) {
  if (key %in% names(meta)) {
    val <- meta[[key]]
    if(is.character(val)) val <- sub(replacement = "   ", x = val, pattern = " // ")
    return(val)
  } else return(NA)
}

setmeta <- function(key = NULL, val = NULL, meta = NULL) {
  for (x in 1:length(key)) meta[[key[x]]] <- val[x]
  return(meta)
}

meta.df2list <- function(meta.df = NULL) {
  return(sapply(seq_len(nrow(meta.df)), function(x) {
    l <- meta.df[x,]
    return(setNames(l[2], l[1]))
  }))
}

list.depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, list.depth)), 0L)

oschp.load <- function(file = NULL) {
  if (is.null(file)) stop(tmsg("Please provide an OSCHP file !"))
  if (!file.exists(file)) stop(tmsg("Provided OSCHP file does not exist !"))
  h5.data <- rhdf5::h5read(file = file, name = "/")
  h5.mlist <- h5.data$Dset_IO_HDF5_Gdh
  if(length(h5.mlist) > 1) {
    `%do%` <- foreach::"%do%"
    h5.meta <- foreach::foreach (a = 1:(length(h5.mlist)-1)) %do% {
      h5.meta.c <- foreach (l = 1:(list.depth(h5.mlist[[a]])-1)) %do% {
        tmp.list <- h5.mlist[[a]][["_&keyvals"]]
        h5.mlist[[a]] <- h5.mlist[[a]][[1]]
        return(meta.df2list(tmp.list))
      }
      names(h5.meta.c) <- foreach (l = h5.meta.c, .combine = "c") %do% {
        return(rev(unlist(strsplit(x = l[["data_source"]], split = "-")))[1])
      }
      return(h5.meta.c)
    }
    names(h5.meta) <- paste0("CEL", 1:(length(h5.mlist)-1))
  } else h5.meta <- list()
  h5.meta$analysis = meta.df2list(h5.mlist[["_&keyvals"]])
  h5.data$Meta <- h5.meta
  h5.data$Dset_IO_HDF5_Gdh <- NULL
  return(h5.data)
}

fpaav <- Vectorize(tools::file_path_as_absolute)

tmsg <- function(text = NULL) { return(paste0(text, " [", Sys.info()[['nodename']], ":", Sys.getpid(), "]")) }

## Vectorization of seq.default()
seq.int2 <- Vectorize(seq.default, SIMPLIFY = FALSE)

## Fast file writing using iotools::write.csv.raw
write.table.fast <- function(x, file = NULL, header = TRUE, sep = "\t", fileEncoding="", row.names = FALSE, ...) {
  if (header) write.table(x = x[NULL,], file = file, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = fileEncoding)
  if(!row.names) rownames(x) <- NULL
  trychk <- try(iotools::write.csv.raw(x = x, file = file, sep = sep, col.names=FALSE, fileEncoding=fileEncoding, append = TRUE, ...))
  if (!is.null(trychk)) {
    print("Fast write failed, using canonical write.table ...")
    write.table(x = x, file = file, sep = sep, row.names = row.names, quote = FALSE)
  }
  gc()
}

## Fast file reader using data.table::fread
read.table.fast <- function(file = NULL, header = TRUE, sep= "\t", row.names = FALSE, ...) {
  if (row.names) {
    if (header) h.df <- read.table(file = file, sep = sep, header = header, nrows = 1, check.names = FALSE)
    data.df <- data.table::fread(input = file, sep = sep, header = FALSE, skip = 1, data.table = FALSE, ...)
    rownames(data.df) <- data.df[,1]
    data.df[,1] <- NULL
    if (header) colnames(data.df) <- colnames(h.df)
  } else {
    data.df <- data.table::fread(input = file, sep = sep, header = header, data.table = FALSE, ...)
  }
  return(data.df)
}

## Function to load data from a HDF5 file, using rhdf5. Returns a list of tables.
hdf5.load <- function (file = NULL) {
  if (is.null(file)) stop("Please provide a HDF5 file !")
  if (!file.exists(file)) stop("Provided HDF5 file does not exist !")
  return(rhdf5::h5read(file = file, name = "/"))
}

## Function to change global option "bitmapType" for PNG plotting on stations without X installed or launched
EaCoN.set.bitmapType <- function(type = "cairo") {
  options(bitmapType = type)
}
