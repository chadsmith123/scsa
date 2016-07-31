#' FACs data importer
#'
#' Imports raw data from a FACs sample. 
#' @param x matrix of cell counts with the "green" values in rows and "red" values in columns
#' @param id sample ID name
#' @param max.green maximum value of green to include in the analysis 
#' @param channels number of channels used in FACs acquisition
#' @param ... additional arguments to be passed
#' @return Returns a list of the frequency of cells in each green/red channel and calculates alphat
#' @references Evenson DP, Larson KL, Jost LK (2002) Sperm Chromatin Structure Assay: Its Clinical Use for Detecting Sperm DNA Fragmentation in Male Infertility and Comparisons With Other Techniques. Journal of Andrology, 23, 25–43.
#'
#' @seealso scsa, scsaDat, findHDS, findAlphat
#' @export
#' @examples
#' data(zf_facs)
#' scsaDat <- scsaDat(x=zf_facs,id="Sample 1",max.green=1023,channels=1024)

scsaDat <- function(x,id,max.green=dim(x)[1],channels=dim(x)[2], ...){
# Function to calculate alphat
scsa.alphat <- function(red,green){
		red/(green+red)
	}
if(missing(id)){id=deparse(substitute(x))}
# Merge imported data from FACS with SCSA parameters
 x <- as.matrix(x)
 rows <- length(x[1,])
 cols <- length(x[,1])
 if(rows != cols){
	stop(paste("The number of rows and columns are not equal. Check that all channels are present in your data file.")) 
 	}
 freq <- as.numeric(paste(x))

# Make a vector of all combinations of red and green values and calculate alphat for each
 red <- rep(1:cols,each=cols)
 green <- seq(1,cols,1)
 color <- data.frame(red,green=green)
 color$alphat <- scsa.alphat(color$red,color$green)
 color$total <- 0.5*(color$red + color$green)

# Merge flow data and color values
 dat.tmp <- cbind(color,freq)
 dat <- subset(dat.tmp,freq != 0) # remove all cells in matrix with zero
# if(missing(max.green)){max.green=cols}
 dat <- subset(dat,green <= max.green)
 out <- list("id"=id,"data"=dat,"channels"=channels,"max.green"=max.green)
return(out)
}
