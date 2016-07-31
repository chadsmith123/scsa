#' Alphat calculator
#' Plots a histogram of FACs input and allows the user to select the Alphat threshold.
#'
#' @param x data frame from scsaDat
#' @param channels number of channels used in FACs acquisition
#' @return Returns Alphat threshold value selected by the user.
#' @references Evenson DP, Larson KL, Jost LK (2002) Sperm Chromatin Structure Assay: Its Clinical Use for Detecting Sperm DNA Fragmentation in Male Infertility and Comparisons With Other Techniques. Journal of Andrology, 23, 25–43.
#' @seealso scsa, scsaDat, findHDS
#' @import grDevices graphics
#' @export
#' @examples
#' \dontrun{
#' data(zf_facs)
#' x <- scsaDat(dat)
#' findAlphat(x)
#' }
findAlphat <- function(x,channels){
 cat("Select Alphat threshold\n")
 par(mfrow=c(1,2))
 plot(freq ~ alphat,type="h",
	data=aggregate(freq ~ alphat,data=x,FUN=sum),
	xlab="Frequency",ylab="Alpha t",main="Select Alphat threshold",
	xlim=c(0,.5))
  scsa.i <- identify(x=x$alphat,y=x$freq,n=1,labels="")
  thresh.alphat <- x[scsa.i,]$alphat
  legend("topleft", paste("Alphat Threshold:",round(thresh.alphat,digits=3)), bty="n")
  abline(v=thresh.alphat) 
  plot(total ~ alphat,data=x, type="n",
	ylab="Total Fluorescence",xlim=c(0,1),xlab="Alpha t")
  points(total ~ alphat,data=x,pch=20,cex=0.001)
  abline(v=thresh.alphat) 
return(thresh.alphat)
}
