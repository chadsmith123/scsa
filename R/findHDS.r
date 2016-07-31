#' High DNA stainability (HDS) calculator
#'
#' Plots a histogram of FACs input and allows the user to select the HDS threshold.
#' @param x data frame from scsaDat
#' @param channels number of channels used in FACs acquisition
#' @param ... additional arguments to be passed
#' @return Returns the HDS threshold selected by the user.
#' @references Evenson DP, Larson KL, Jost LK (2002) Sperm Chromatin Structure Assay: Its Clinical Use for Detecting Sperm DNA Fragmentation in Male Infertility and Comparisons With Other Techniques. Journal of Andrology, 23, 25–43.
#' @seealso scsa, scsaDat, findAlphat 
#' @import grDevices graphics
#' @export
#' @examples
#' \dontrun{
#' data(zf_facs)
#' x <- scsaDat(zf_facs)
#' findHDS(x)
#' }
findHDS <- function(x, channels, ...){
  cat("Select HDS threshold\n")
  par(mfrow=c(1,2))
  palette(gplots::rich.colors(n=25))
  plot(freq ~ green,type="h",
	data=aggregate(freq ~ green,data=x,FUN=sum),
	xlim=c(0,channels),xlab="Green",ylab="Frequency",main="Select HDS threshold")
  hds.i <- identify(x=x$green,y=x$freq,n=1,labels=round(x$green,digits=2))
  thresh.HDS <- x[hds.i,]$green
  abline(v=thresh.HDS)
  legend("topleft", paste("HDS Threshold:",round(thresh.HDS,digits=3)), bty="n")
  plot(x$red ~ x$green, type="n",
	xlim=c(0,channels),ylim=c(0,channels),	
	pch=20,ylab="Green",xlab="Red")
  points(x$green ~ x$red,pch=20,cex=.0001)
  points(x$green ~ x$red,pch=20,cex=.0001,col=x$freq)
  abline(h=thresh.HDS)
return(thresh.HDS)
	}
