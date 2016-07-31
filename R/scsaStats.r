#' SCSA statistics calculator
#'
#' Calcuates SCSA statistics for a sample.
#' @param x list from scsaDat output
#' @param id sample id
#' @param thresh.alphat Alphat threshold value. If not specified, a plot will be generated for the user to select the threshold.
#' @param thresh.HDS High DNA Staining threshold value. If not specified, a plot will be generated for the user to select the threshold.
#' @param ... additional arguments to be passed
#' @details 
#' DNA Fragmentation index (DFI) = proportion of cells above the Alphat threshold
#'
#' High DNA Stainability (HDS) = proportion of cells above the HDS threshold
#'
#' N = total number of cells analyzed
#'
#' X.alphat = Weighted mean DFI of all cells analyzed
#'
#' sd.alphat = Standard deviation of the DFI of all cells analyzed
#'
#' n.comp = number of cells above the Alphat threshold
#'
#' X.comp = weighted mean DFI of cells above the Alphat threshold
#'
#' sd.comp = standard deviation of DFI above the Alphat threshold
#'
#' thresh.alphat = user-selected Alphat threshold
#'
#' thresh.HDS = user-selected HDS threshold
#'
#' X.red = mean red fluorescence of all cells
#'
#' X.green = mean green fluorescence of all cells
#' @return Returns a list of SCSA statistics
#' @references Evenson DP, Larson KL, Jost LK (2002) Sperm Chromatin Structure Assay: Its Clinical Use for Detecting Sperm DNA Fragmentation in Male Infertility and Comparisons With Other Techniques. Journal of Andrology, 23, 25–43.
#'
#' @seealso scsa, scsaDat, findHDS, findAlphat
#' @import stats utils
#' @export
#' @examples
#' data(zf_facs)
#' x <- scsaDat(zf_facs)
#' scsaStats(x,id="Sample 1",thresh.alphat=0.27,thresh.HDS=527)
#' # Values for SCSA thresholds can be manually entered to generate the statistics.
scsaStats <- function(x, id=x$id, thresh.alphat, thresh.HDS, ...){
 dat <- x$data
 alphat <- dat$alphat
 green <- rep(dat$green,dat$freq)
 red <- rep(dat$red,dat$freq)
 channels <- x$channels

# If thresh.alpha or thresh.HDS are not specified, prompt user to provide them.
 if(missing(thresh.alphat)){
  thresh.alphat <- findAlphat(dat)
  user <- menu(c("Accept ","Do over "))
  if(user==2){
   while(user==2){
   thresh <- findAlphat(dat)
   user <- menu(c("Accept ","Do over "))
			}
		}	
	}

 if(missing(thresh.HDS)){
  thresh.HDS <- findHDS(dat, channels)
  user <- menu(c("Accept ","Do over "))
  if(user==2){
   while(user==2){
   thresh.HDS <- findHDS(dat, channels)
   user <- menu(c("Accept ","Do over "))
			}
		}	
	}
# Generate stats table
 dat.comp <- subset(dat, alphat > thresh.alphat)
 dat2.comp <- rep(dat.comp$alphat,dat.comp$freq)
 dat.hds <- subset(dat, green > thresh.HDS)
list(
 id = id,
 DFI = sum(dat.comp$freq/sum(dat$freq)),
 HDS = sum(dat.hds$freq/sum(dat$freq)),
 X.alphat = weighted.mean(dat$alphat,w=dat$freq),
 sd.alphat = sd(rep(dat$alphat,dat$freq)),
 N = sum(dat$freq),
 X.comp = weighted.mean(dat.comp$alphat,w=dat.comp$freq),
 sd.comp = sd(dat2.comp),
 n.comp = sum(dat.comp$freq),
 thresh.alphat = thresh.alphat,
 thresh.HDS = thresh.HDS,
 X.green=mean(dat$green),
 X.red=mean(dat$red)
 )
}
