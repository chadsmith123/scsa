# Updated 15.05.12
# SCSA analysis file
alphat <- function(red,green){
        red/(green+red)
        }

## Calculate alphat, green, red and total flourescence for a 1024 x 1024 matrix
# Make vector of all possible alphat values in 1024 x 1024 resolution
seq1024 <- seq(1,1024,1)
alphat.1024 <- NULL

# Make a vector of all combinations of red and green values
red <-NULL
for(i in 1:1024){red <- append(red,rep(i,1024))}
color <- data.frame(red,green=seq1024)
color$alphat <- alphat(color$red,color$green)
color$total <- 0.5*(color$red + color$green)

scsa.dat <- function(x){
# Merge imported data from FACS with SCSA parameters
 x1 <- as.matrix(x)
 freq <- as.numeric(paste(x1))
 dat.tmp <- cbind(color,freq)
 dat.tmp <- subset(dat.tmp,freq != 0)
 dat <- subset(dat.tmp,green < 1023)  
 return(dat)
        }

plot.alphat <- function(x,id,thresh){
  par(mfrow=c(1,2))
  #palette(rich.colors(n=25))
  #plot(x$red ~ x$green, type="n",
#	xlim=c(0,1024),ylim=c(0,1024),	
#	pch=20,col=x$freq,main=id,,ylab="FITC",xlab="TX Red")
 # points(x$green ~ x$red,pch=20,cex=.0001,col=x$freq)
  plot(freq ~ alphat,type="h",
	data=aggregate(freq ~ alphat,data=x,FUN=sum),
	col="grey30",main=id,xlim=c(0,.5))
	#col="#000040FF",main=id,xlim=c(0,.5))
  scsa.i <- identify(x=x$alphat,y=x$freq,n=1,labels=round(x$alphat,digits=2))
  thresh <- x[scsa.i,]$alphat
  abline(v=thresh) 
  plot(total ~ alphat,data=x, type="n",
	ylab="Total Fluorescence",xlim=c(0,1))
  #points(total ~ alphat,data=x,pch=20,col="#000040FF",cex=0.001)
  points(total ~ alphat,data=x,pch=20,col="grey30",cex=0.001)
  abline(v=thresh) 
return(thresh)
}

plot.hgrn <- function(x,id,thresh.hgrn){
  par(mfrow=c(1,2))
  palette(rich.colors(n=25))
  plot(freq ~ green,type="h",
	data=aggregate(freq ~ green,data=x,FUN=sum),
	#col="green4",main=id,xlim=c(0,1024))
	col="grey30",main=id,xlim=c(0,1024))
	#col="#000040FF",main=id,xlim=c(0,1024))
  hgrn.i <- identify(x=x$green,y=x$freq,n=1,labels=round(x$green,digits=2))
  thresh.hgrn <- x[hgrn.i,]$green
  abline(v=thresh.hgrn)
  plot(x$red ~ x$green, type="n",
	xlim=c(0,1024),ylim=c(0,1024),	
	pch=20,col=x$freq,main=id,,ylab="FITC",xlab="TX Red")
  points(x$green ~ x$red,pch=20,cex=.0001,col=x$freq)
  abline(h=thresh.hgrn)
return(thresh.hgrn)
	}

scsa.stats <- function(x,id=as.character(substitute(x)),thresh=1,thresh.hgrn=1,save.plot=1){
# Calculate statistics. Plots saved to directory png/. Click on plot to set COMP & HGRN threshold.
 green <- rep(x$green,x$freq)
 red <- rep(x$red,x$freq)
 if(thresh == 1){
  thresh <- plot.alphat(x,id,thresh)
  user <- menu(c("Accept ","Do over "))
  if(user==2){
   while(user==2){
   thresh <- plot.alphat(x,id,thresh)
   user <- menu(c("Accept ","Do over "))
			}
		}	
	}
 if(thresh.hgrn == 1){
  thresh.hgrn <- plot.hgrn(x,id,thresh.hgrn)
  user <- menu(c("Accept ","Do over "))
  if(user==2){
   while(user==2){
   thresh.hgrn <- plot.hgrn(x,id,thresh.hgrn)
   user <- menu(c("Accept ","Do over "))
			}
		}	
	}
 x.comp <<- subset(x, alphat > thresh)
 x2.comp <<- rep(x.comp$alphat,x.comp$freq)
 p.comp <- sum(x.comp$freq/sum(x$freq))*100
 X.comp <- weighted.mean(x.comp$alphat,w=x.comp$freq)*100
 n.comp <- sum(x.comp$freq)
 sd.comp <- sd(x2.comp)*100

 x.hgrn <- subset(x, green > thresh.hgrn)
 hgrn <- sum(x.hgrn$freq/sum(x$freq))*100

 X.alphat <- weighted.mean(x$alphat,w=x$freq)*100
 sd.alphat <- sd(rep(x$alphat,x$freq)*100)
 X.green <- weighted.mean(x$green,w=x$freq)
 X.red <- weighted.mean(x$red,w=x$freq)
 N <- sum(x$freq)
 out <- data.frame(
	id,
	thresh,
	p.comp,
	X.comp,
	sd.comp,
	n.comp,
	X.alphat,
	sd.alphat,
	thresh.hgrn,
	hgrn,
	X.green,
	X.red,
	N
	)
# Make plots of green vs red, alphat histogram, total fluor. vs alphat
  #layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE),width=c(3,3,1))
  par(mar=c(4.5,4,1.5,2))
  layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow = TRUE),width=c(3,3,1))
  palette(rich.colors(n=25))
  plot(x$green ~ x$red, type="n",
	xlim=c(0,1024),ylim=c(0,1024),	
	pch=20,col=x$freq,main=id,,ylab="FITC",xlab="TX Red")
  points(x$green ~ x$red,pch=20,cex=.0001,col=x$freq)
  #text(x=750,y=1000,
  text(x=0,y=1000,
	labels=paste(sep="",
	"X.alphat\t\t",round(X.alphat,digits=2),"%","\n",
	"COMP.alphat\t",round(p.comp,digits=2),"%","\n",
	"HGRN\t\t",round(hgrn,digits=2),"%"),adj=c(0,0))
  abline(h=thresh.hgrn)
  plot(freq ~ alphat,type="h",
	data=aggregate(freq ~ alphat,data=x,FUN=sum),
	col="#000040FF",ylab="Frequency",xlab="AlphaT",xlim=c(0,1))
  abline(v=thresh)
  plot(total ~ alphat,data=x, type="n",
	ylab="Total Fluorescence",xlab="AlphaT",xlim=c(0,1))
  points(total ~ alphat,data=x,pch=20,col="#000040FF",cex=0.001)
  abline(v=thresh)
  plot(freq ~ green,type="h",
	data=aggregate(freq ~ green,data=x,FUN=sum),
	col="#000040FF",ylab="Frequency",xlab="Green",xlim=c(0,1024))
  abline(v=thresh.hgrn)
  par(mfrow=c(1,1))
 if(save.plot==1){
  savePlot(filename=paste("png/",id,".png",sep=""),type="png")
		}
 return(out)
}

#scsa.dat.old <- function(x){
# x1 <- as.matrix(x)
# out=NULL
# for(i in 1:1024){tmp=c(x1[i,]);out=append(out,tmp)}
#        dat.tmp <<- cbind(alphat=alphat.1024,freq=as.numeric(out))
#        dat2.tmp <- as.data.frame(dat.tmp)
#        dat <<- subset(dat2.tmp,freq != 0)
#	return(dat)
#        #print("Output in 'dat'")
#        }

## Old way of getting alphat values
# increment the denominator +1 and calculate alphat 
#for(i in 1:1024){tmp <- alphat(i,seq1024); alphat.1024  <<- append(alphat.1024,tmp)}
## calculate total flourescence
#total.1024 <- NULL
#for(i in 1:1024){tmp <- 0.5*(i+seq1024); total.1024  <<- append(total.1024,tmp)}

## Create data frame of SCSA frequency distribution and alphat
#scsa.dat <- function(x){
# x1 <- as.matrix(x)
# freq <- as.numeric(paste(x1))
## dat.tmp <<- cbind(alphat=alphat.1024,freq)
# dat.tmp <<- cbind(total=total.1024,alphat=alphat.1024,freq)
# dat2.tmp <- as.data.frame(dat.tmp)
# dat <<- subset(dat2.tmp,freq != 0)
# return(dat)
#        }

# Updated 12.08.13
# Calculates peak green. Uncomment below if you want to manually specify it from a plot.  
scsa.pkgrn <- function(x,id=as.character(substitute(x)),pk.grn,plot=F){
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
	col="grey30",main=id,xlim=c(0,1024))
# Uncomment next two lines to set peak green  manually from a plot
  #pkgrn.i <- identify(x=x$green,y=x$freq,n=1,labels=round(x$green,digits=2)) # Next two lines allow manual setting of pkgrn
  #pkgrn <- x[pkgrn.i,]$green
  abline(v=pk.grn,lty=1)
  abline(v=median.grn,lty=2)
  abline(v=X.green,lty=3)
  text(x=800,y=30,pos=4,cex=0.75,label="Solid\tPeak")
  text(x=800,y=20,pos=4,cex=0.75,label="Dashed\tMedian")
  text(x=800,y=10,pos=4,cex=0.75,label="Dotted\tMean")
  savePlot(filename=paste("pkgrn/",id,".png",sep=""),type="png")
  	}
out <- data.frame(id,pk.grn,median.grn)
return(out)
	}
