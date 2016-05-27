# Updated 15.05.12
# SCSA analysis file
require(gplots)
 scsa.alphat <- function(red,green){
	red/(green+red)
	}
scsa <- function(x, ...) UseMethod("scsa")
#scsa.default <- function(x,id,channels,max.green,thresh.alphat,thresh.HDS, ...){
scsa.default <- function(x, ...){
# DAT <- scsaDat(x,channels,max.green, ...) 
# DAT$stats <- scsaStats(DAT,thresh.alphat,thresh.HDS, ...)
 DAT <- scsaDat(x,id=deparse(substitute(x)),channels=1024) 
# DAT$stats <- scsaStats(DAT,thresh.alphat=NULL,thresh.HDS=NULL)
 DAT$stats <- scsaStats(DAT, ...)
 class(DAT) <- "scsa"
 DAT
}
print.scsa <- function(x, ...){
 cat(paste("ID:",x$id,"\n"))
 cat(paste("DFI:",round(x$stats$DFI,digits=3),'\n'))
 cat(paste("HDS:",round(x$stats$HDS,digits=3),'\n'))
}
summary.scsa <- function(x){
	data.frame(x$stats)
}

findAlphat <- function(x,channels){
 par(mfrow=c(1,2))
 plot(freq ~ alphat,type="h",
	data=aggregate(freq ~ alphat,data=x,FUN=sum),
	xlab="Frequency",ylab="Alpha t",
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

findHDS <- function(x,channels){
  par(mfrow=c(1,2))
  palette(rich.colors(n=25))
  plot(freq ~ green,type="h",
	data=aggregate(freq ~ green,data=x,FUN=sum),
	xlim=c(0,channels),xlab="Green",ylab="Frequency")
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

scsaDat <- function(x,id,max.green,channels){
# x= matrix of flow cytometry values
# skip.row=number of rows to skip in flow cytometry file
# max.green = maximum of green values to retain
if(missing(id)){id=deparse(substitute(x))}
# Merge imported data from FACS with SCSA parameters
 x <- as.matrix(x)
 rows <- length(x[1,])
 cols <- length(x[,1])
 if(rows != cols){
	stop(paste(rows,"and",cols,"columns. These values are not equal. Input the number of channels used to collect the flow cytometry data by setting the channels variable."))
 	}
 freq <- as.numeric(paste(x))

# Make a vector of all combinations of red and green values and calculate alphat for each
 red <- rep(1:cols,each=cols)
 color <- data.frame(red,green=seq(1,cols,1))
 color$alphat <- scsa.alphat(color$red,color$green)
 color$total <- 0.5*(color$red + color$green)

# Merge flow data and color values
 dat.tmp <- cbind(color,freq)
 dat <- subset(dat.tmp,freq != 0) # remove all cells in matrix with zero
 if(missing(max.green)){max.green=cols}
 dat <- subset(dat,green <= max.green)
 out <- list("id"=id,"data"=dat,"channels"=rows,"max.green"=max.green)
return(out)
}
scsaStats <- function(x,id,thresh.alphat,thresh.HDS){
 if(missing(id)){id <- x$id}
 dat <- x$data
# Calculate statistics. Plots saved to directory png/. Click on plot to set COMP & HGRN threshold.
 green <- rep(dat$green,dat$freq)
 red <- rep(dat$red,dat$freq)
 if(missing(thresh.alphat)){
  thresh.alphat <- findAlphat(dat,channels=x$channels)
  user <- menu(c("Accept ","Do over "))
  if(user==2){
   while(user==2){
   thresh <- findAlphat(dat)
   user <- menu(c("Accept ","Do over "))
			}
		}	
	}
 if(missing(thresh.HDS)){
  thresh.HDS <- findHDS(dat,channels=x$channels)
  user <- menu(c("Accept ","Do over "))
  if(user==2){
   while(user==2){
   thresh.HDS <- findHDS(dat,thresh.HDS)
   user <- menu(c("Accept ","Do over "))
			}
		}	
	}
 dat.comp <- subset(dat, alphat > thresh.alphat)
 dat2.comp <- rep(dat.comp$alphat,dat.comp$freq)
 dat.hds <- subset(dat, green > thresh.HDS)
list(
 id = id,
 DFI = sum(dat.comp$freq/sum(dat$freq)),
 HDS = sum(dat.hds$freq/sum(dat$freq)),
 N = sum(dat$freq),
 X.comp = weighted.mean(dat.comp$alphat,w=dat.comp$freq),
 sd.comp = sd(dat2.comp),
 n.comp = sum(dat.comp$freq),
 X.alphat = weighted.mean(dat$alphat,w=dat$freq),
 sd.alphat = sd(rep(dat$alphat,dat$freq)),
 thresh.alphat = thresh.alphat,
 thresh.HDS = thresh.HDS,
 X.green=mean(dat$green),
 X.red=mean(dat$red)
 )
}
plot.scsa <- function(x,title=x$id, ...){
# Make plots of green vs red, alphat histogram, total fluor. vs alphat
  par(mar=c(4.5,4,1.5,2),oma=c(0,0,3,0),cex=0.5)
  layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow = TRUE),width=c(3,3,1))
  palette(rich.colors(n=25))
  plot(green ~ red, type="n",data=x$data,
	xlim=c(0,x$channels),ylim=c(0,x$channels),	
	pch=20,ylab="Green",xlab="Red")
  points(green ~ red,pch=20,col=freq,data=x$data, ...)
  points(green ~ red,pch=20,col="red",data=subset(x$data,alphat >= x$stats$thresh.alphat), ...)
  legend("topleft", bty="n",paste(
	" DFI:",round(x$stats$DFI,digits=2),"\n",
	"HDS:",round(x$stats$HDS,digits=2))
	) 
  abline(h=x$stats$thresh.HDS)
  plot(freq ~ alphat,type="h",
	data=aggregate(freq ~ alphat,data=x$data,FUN=sum),
	col="#000040FF",ylab="Frequency",xlab="Alpha t",xlim=c(0,1))
  abline(v=x$stats$thresh.alphat)
  legend("topright",bty="n",paste("Alphat Threshold",round(x$stats$thresh.alphat,digits=3)))
  plot( total ~ alphat,data=x$data, type="n",
	ylab="Total Fluorescence",xlab="Alpha t",xlim=c(0,1))
  points(total ~ alphat,data=x$data,pch=20,col=freq, ...)
  points(total ~ alphat,col="red",pch=16,data=subset(x$data,alphat >= x$stats$thresh.alphat), ...)
  abline(v=x$stats$thresh.alphat)
  plot(freq ~ green,type="h",
	data=aggregate(freq ~ green,data=x$data,FUN=sum),
	col="green4",ylab="Frequency",xlab="Green",xlim=c(0,x$channels), ...)
  abline(v=x$stats$thresh.HDS)
  legend("topright",bty="n",paste("HDS Threshold",x$stats$thresh.HDS))
  mtext(title,outer=T,cex=1.5)
}
# Updated 12.08.13
# Calculates peak green. Uncomment below if you want to manually specify it from a plot.  
scsa.pkgrn <- function(x,title=as.character(substitute(x)),pk.grn,plot=F){
  par(mfrow=c(1,1))
  palette(rich.colors(n=25))
  x.agg <- aggregate(freq ~ green,data=x,FUN=sum)	# Adds the frequencies of each channel
  pk.grn.o <- x[which(x$freq==max(x$freq)),]$green	# Finds the channels with the maximum frequency to locate the peaks
  pk.grn <- mean(pk.grn.o)				# Calculates the mean of the channels to estimate the peak, as max(x$freq) can have the same value for more than one channel
  median.grn <- median(x.agg$green)			# Alternately, the median
  X.green <- mean(x.agg$green)				# Just the mean green as calculate in scsa.stats
if(plot=="T"){
  plot(freq ~ green,type="h",
	data=x.agg,
	col="grey30",main=title,xlim=c(0,1024))
# Uncomment next two lines to set peak green  manually from a plot
  #pkgrn.i <- identify(x=x$green,y=x$freq,n=1,labels=round(x$green,digits=2)) # Next two lines allow manual setting of pkgrn
  #pkgrn <- x[pkgrn.i,]$green
  abline(v=pk.grn,lty=1)
  abline(v=median.grn,lty=2)
  abline(v=X.green,lty=3)
  text(x=800,y=30,pos=4,cex=0.75,label="Solid\tPeak")
  text(x=800,y=20,pos=4,cex=0.75,label="Dashed\tMedian")
  text(x=800,y=10,pos=4,cex=0.75,label="Dotted\tMean")
  savePlot(filename=paste("pkgrn/",title,".png",sep=""),type="png")
  	}
out <- data.frame(title,pk.grn,median.grn)
return(out)
	}
