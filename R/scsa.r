#' Sperm Chromatin Structure Assay (SCSA) FACs Analyzer
#'
#' @param x square matrix of raw FACs data from a sample
#' @param ... additional arguments to be passed. 
#' @return Returns a scsa object
#' @references Evenson DP, Larson KL, Jost LK (2002) Sperm Chromatin Structure Assay: Its Clinical Use for Detecting Sperm DNA Fragmentation in Male Infertility and Comparisons With Other Techniques. Journal of Andrology, 23, 25–43.
#' @seealso scsaDat, scsaStats, findHDS, findAlphat
#' @export
#' @examples
#' data(zf_facs)
#'
#' \dontrun{
#' scsa(zf_facs)
#' # Runs user through the SCSA pipeline.
#' }
#' dat <- scsa(zf_facs,id="sample1", thresh.alphat=0.25,thresh.HDS=500, max.green=1023)
#' # Returns SCSA object with manual configuration
#'
#' print(dat)
#' # Returns sample ID, DFI and HDS
#'
#' stats.scsa(dat)
#' # Returns SCSA statistics 
#' 
#' plot(dat)
#' # Returns a plot of FACs data for the sample.

scsa <- function(x, ...) UseMethod("scsa", x)

#' @export
scsa.default <- function(x, id=deparse(substitute(x)), ...){
 obj <- scsaDat(x, id, ...)
 obj$stats <- scsaStats(obj, ...)
 class(obj) <- "scsa"
 return(obj)
}

#' @export
print.scsa <- function(x, ...){
 cat("// SCSA object ///\n") 
 cat(paste("ID:",x$id,"\n"))
 cat(paste("DFI:",round(x$stats$DFI,digits=3),'\n'))
 cat(paste("HDS:",round(x$stats$HDS,digits=3),'\n'))
}

#' @export
plot.scsa <- function(x, ...){
 scsaPlot(x, ...)
}

#' Statistics for SCSA objects
#' @param x an object of class \code{"scsa"}
#' @param ... additional arguments to be passed
#' @return Returns a data frame of summary statistics. See scsaStats for variable descriptions.
#' @seealso scsaDat, scsaStats, findHDS, findAlphat
#' @export
#' @examples
#' data(zf_facs)
#' dat <- scsa(zf_facs,thresh.alphat=0.25,thresh.HDS=500)
#' stats.scsa(dat)
stats.scsa <- function(x, ...){
 obj <- data.frame(x$stats)
 return(obj)
}
