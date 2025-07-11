.onLoad <- function(libname, pkgname){
   eacon.dir <- system.file(package = 'EaCoN')
   cache.dir <- paste0(eacon.dir, '/extdata/cache')
   options(EACON_CACHE=cache.dir)
}