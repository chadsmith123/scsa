#' SCSA plot
#'
#' @param x list from scsaDat function
#' @param title title of the plot
#' @param ... additional arguments to be passed
#' @return Returns summary plots SCSA data
#' @references Evenson DP, Larson KL, Jost LK (2002) Sperm Chromatin Structure Assay: Its Clinical Use for Detecting Sperm DNA Fragmentation in Male Infertility and Comparisons With Other Techniques. Journal of Andrology, 23, 25–43.
#'
#' @seealso scsa, scsaDat, findHDS, findAlphat
#' @import grDevices graphics stats
#' @export
#' @examples
#' \dontrun{
#' scsaPlot(x)
#' }
scsaPlot <- function(x,title=x$id, ...){
  dat <- x$data
  alphat <- x$alphat
  freq <- x$freq
  stats <- x$stats
  thresh.alphat <- x$stats$thresh.alphat
  thresh.HDS <- x$stats$thresh.HDS

  par(mar=c(4.5,4.5,1.5,2),oma=c(0,0,3,0))
  #par(mar=c(4.5,4,1.5,2),oma=c(0,0,3,0))
  par(cex.axis=1.5,cex.lab=1.5)
  layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow = TRUE),widths=c(3,3,1))
  palette(gplots::rich.colors(n=25))

# Scatterplot of green ~ red   
  plot(green ~ red, type="n",data=dat,
	xlim=c(0,x$channels),ylim=c(0,x$channels),	
	pch=20,ylab="Green",xlab="Red")
  with(dat, points(green ~ red,pch=20,col=freq, ...))
  points(green ~ red,pch=20,col="red",data=subset(dat,alphat >= thresh.alphat), ...)
  legend("topleft", bty="n",paste(
	" DFI:",round(x$stats$DFI,digits=2),"\n",
	"HDS:",round(x$stats$HDS,digits=2))
	) 
  abline(h=x$stats$thresh.HDS)

# Frequency histogram of Alphat
  plot(freq ~ alphat,type="h",
	data=aggregate(freq ~ alphat,data=dat,FUN=sum),
	col="#000040FF",ylab="Frequency",xlab="Alpha t",xlim=c(0,1))
  abline(v=thresh.alphat)
  legend("topright",bty="n",paste("Alphat Threshold",round(thresh.alphat,digits=3)))

# Frequency histogram of HDS
  plot(freq ~ green,type="h",
	data=aggregate(freq ~ green,data=dat,FUN=sum),
	col="green4",ylab="Frequency",xlab="Green",xlim=c(0,x$channels), ...)
  abline(v=x$stats$thresh.HDS)
  legend("topright",bty="n",paste("HDS Threshold",thresh.HDS))
  mtext(title,outer=T,cex=1.5)

# Scatterplot of Total Fluorescence ~ Alphat
  plot( total ~ alphat,data=dat, type="n",
	ylab="Total Fluorescence",xlab="Alpha t",xlim=c(0,1))
  points(total ~ alphat,data=dat,pch=20,col=freq, ...)
  points(total ~ alphat,col="red",pch=16,data=subset(dat,alphat >= thresh.alphat), ...)
  abline(v=thresh.alphat)
}
