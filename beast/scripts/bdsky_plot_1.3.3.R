# This script plots the epidemiological parameters of a BDSKY analysis with their HPD intervals. 
# Written by Denise Kühnert (denise.kuehnert@gmail.com)
#
# start R from Terminal and run script with commandline: source(bdsky_plot.R) 
#
# Input: 
#	- a text file listing the log files that should be analyzed (see example file loglist.txt)
#	- burnin (the percentage of samples that should be ignored from the log file)
#	- date of the most recent sample (for plotting from past to present)
#
# Assumptions:
# - The given log file contains parameters named "reproductiveNumber" or "reproductiveNumberx" with index x,
#	"becomeUninfectiousRate" / "becomeUninfectiousRatex" and 
#	"samplingProportion" / "samplingProportionx"
# - The number of reproductiveNumber parameters determines the intervalnumber
# - The parameters becomeUninfectiousRate and samplingProportion are assumed to be either constant or to have the same dimension as reproductiveNumber.


library(s20x)
library(boa)
library(Hmisc)
library(miscTools)
# input : a matrix M and a column ascii name 
# output : the numeric index of the column 
colnameindex = function(M , colname0) 
{ 
colsnames = names(M[1,]); 
theindex = which(colsnames==colname0); 
return(theindex); 
}

# /* get user input */
cat("Please enter the name of the file that stores the path to all log files to be analyzed: ")
x <- readLines(file("stdin"),1)

cat("Please enter the percentage you would like to use as burnin (e.g. 10 for 10%): ")
burninpercent <- as.integer(readLines(file("stdin"),1))

cat("Please enter the sampling date of the most recent sample: ")
recent <- as.integer(readLines(file("stdin"),1))

cat("Please enter the grid size (can be bigger than the intervalnumber to get a smoother plot): ")
gridSize <- as.integer(readLines(file("stdin"),1))

# /* read log list */ 
loglist = read.table(x, as.is=TRUE, header=FALSE)

closeAllConnections()


# todo: ask if delta / s constant

for(i in 1:length(loglist[,1])){

	# /* read and assign file from log list */
	assign(paste("log", i, sep=''), read.table(loglist[i,], header=T))
	attach(get(paste("log", i, sep='')))
	
	R0_names = names(get(paste("log", i, sep='')))[which(regexpr("reproductiveNumber.", names(get(paste("log", i, sep=''))))>0)]
	delta_names = names(get(paste("log", i, sep='')))[which(regexpr("becomeUninfectiousRate", names(get(paste("log", i, sep=''))))>0)]
	sampling_names = names(get(paste("log", i, sep='')))[which(regexpr("samplingProportion", names(get(paste("log", i, sep=''))))>0)]
	rhosampling=0
	if (length(sampling_names)==0){
    	sampling_names = names(get(paste("log", i, sep='')))[which(regexpr("rho", names(get(paste("log", i, sep=''))))>0)]
    	rhosampling=1
	}
	temp <- get(paste("log", i, sep=''))
	treeheights = get(paste("log", i, sep=''))$Tree.height # [,match("Tree.height", tolower(names(get(paste("log", i, sep='')))))]
	origins = get(paste("log", i, sep=''))$origin
	width=median(origins)
	
	nsamples= length(get(R0_names[1]))
	burnin = round(burninpercent*nsamples/100)

	intervalNumber = length(R0_names)
	if (intervalNumber > gridSize) {gridSize = intervalNumber}

	medians = matrix(data=NA, nrow= 1, ncol = gridSize)
	medians_G = matrix(data=NA, nrow= 1, ncol = gridSize)
	medians_H = matrix(data=NA, nrow= 1, ncol = gridSize)
	
	hpd_F = matrix(data=NA, nrow= 2, ncol = gridSize)
	hpd_G = matrix(data=NA, nrow= 2, ncol = gridSize)
	hpd_H = matrix(data=NA, nrow= 2, ncol = gridSize)
	
	F = matrix(data=NA, nrow=nsamples-burnin, ncol = gridSize) #reproductiveNumber
	G = matrix(data=NA, nrow=nsamples-burnin, ncol = gridSize) #becomeuninfectiousRate
	H = matrix(data=NA, nrow=nsamples-burnin, ncol = gridSize) #samplingProportion
	
	step = width/(gridSize-1)
	F_times = seq(recent-width, recent, step)
	
	for(k in 1:(nsamples-burnin)){
	
		time = origins[k+burnin]
			
		for (l in 1:length(F_times)){
			currentWidth = time / intervalNumber
			index = ceiling(intervalNumber - (recent - F_times[l])/currentWidth )
			
			F[k,l] = get(R0_names[max(index,1)])[k+burnin]
			if (length(delta_names)==length(R0_names))	
							G[k,l] = get(delta_names[max(index,1)])[k+burnin]
					else	G[k,l] = get(delta_names[1])[k+burnin]
			if (length(sampling_names)==length(R0_names))	
							H[k,l] = get(sampling_names[max(index,1)])[k+burnin]
					else	H[k,l] = get(sampling_names[1])[k+burnin]
				
		}
	}

	for(j in 1:gridSize){
		if (length(which(F[,j]!="NA")) > (nsamples / 10)) {
			medians[1,j] = median(F[,j],na.rm=T)
			medians_G[1,j] = median(G[,j],na.rm=T)
			medians_H[1,j] = median(H[,j],na.rm=T)
			hpd_F[,j] = boa.hpd(F[which(F[,j]!="NA"),j], 0.05)[1:2]
			hpd_G[,j] = boa.hpd(G[which(G[,j]!="NA"),j], 0.05)[1:2]
			hpd_H[,j] = boa.hpd(H[which(H[,j]!="NA"),j], 0.05)[1:2]
		}
	}
	
	#layout20x(3,1)
	pdf(paste0('plot', i, '.pdf'))

# /* plot reproductiveNumber */
	plot(1, ylab='Reproductive number', xlim=c(recent-width, recent), ylim=c(0,max(hpd_F[2,],na.rm=T)*1.1), xlab="Year", col='white', main = '')
	minor.tick(nx=5, ny=2, tick.ratio=.2)
	polygon(c(F_times,rev(F_times)), c(hpd_F[2,], rev(hpd_F[1,])),col = "grey90", border = NA)
	lines(c(F_times), c(medians[1,]), type='l')
	abline(1,0,col='grey')

	if(1==1){
	# /* plot become non-infectious rate */ 
		plot(1, ylab=expression(delta), xlim=c(recent-width, recent), ylim=c(0,max(hpd_G[2,],na.rm=T)*1.1), xlab="Year", col='white', main = '')
		minor.tick(nx=5, ny=2, tick.ratio=.2)
		polygon(c(F_times,rev(F_times)), c(hpd_G[2,], rev(hpd_G[1,])),col = "grey90", border = NA)
		lines(F_times, medians_G[1,], type='l')

	if(rhosampling==0){
	# /* plot samplingProportion */
		plot(0, ylab=expression(s), xlim=c(recent-width, recent), ylim=c(0,max(hpd_H[2,],na.rm=T)*1.1),  xlab="Year", col='white', main = '')
		polygon(c(F_times,rev(F_times)), c(hpd_H[2,], rev(hpd_H[1,])),col = "grey90", border = NA)
		lines(F_times, medians_H[1,], type='l')
		}
	}
	dev.off()
}
