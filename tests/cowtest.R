


######################

#Check calculated masses using: 
#http://education.expasy.org/student_projects/isotopident/htdocs/

######################

library("q2e")
library("bioarch")


######################
#For comparison: 
library("OrgMassSpecR")


oIsoDist <- function(seq){

	m <-ConvertPeptide(seq)
	d <-IsotopicDistribution(m)
	return(d)

}




#load the q23 library
require(q2e)


#Specify the `margins' of the sample 
#moff = 1




#Deamidation changes Qs to Es
#Oxidation adds 16 to the mass


#############################################
# Let's pick an example spot
spot <- 6
tseq <- "IGQPGAVGPAGIR"
tseq <- "GSAGPPGATGFPGAAGR"



#tseq<- "GIPGEFGLPGPAGAR"

#tseq <- "IGEPGAVGPAGIR"

#tseq <- "GVPGPPGAVGPAGK"
#         GVPGPPGAVGPAGK
#tseq = "GIPGEFGLPGPAGAR"
#tseq = "GSAGPPGATGFPGAAGR"
 
#a.a.s used: 
#tseq = "A GIPQR  V"
#tseq = "GIPGEFGLPGPAGAR"
#tseq = "AFG P RST"
###############################################
#Now load the line data: 

#load the test data
if(!exists("testdata"))
testdata <- bioarch_loadBrukerXML("/home/sjh/Desktop/sjh/bioarch/SF/SFData/20140307_SF_UPenn164-181/")


###############################################
#This should be a function: 


message(sprintf("Generating peaks from sequence %s",tseq))

message(sprintf("Working on spot %d",spot))



###############################################
q2e_plotalign <- function(tseq,spot,myby = 0.005,oxidation = 0, deamidation = FALSE){

	message("deamidation is ",deamidation)

	#deamidation
	if(deamidation){
		tseq2 <- gsub("Q","E",tseq)
		tseq <- tseq2
	}

	mylagmax <- 1/myby
	cd1 <- q2e_isodists(tseq)
	
	#handle the oxidation
	if(oxidation > 0)
		cd1$mass <- cd1$mass + (oxidation * 16)
	
	
	#set up the plot area: (4 rows, 1 column)
	par(mar=c(0.9,2.3,0.9,.3), mfrow = c(4,1), oma=c(5,0,2,0))

	#Let's do some setting up:
	myxlim = c(lbl = floor(min(cd1$mass))-1,ubl = ceiling(max(cd1$mass))+1)

	#grab the data for the range we are intested in 
	xm <-       spot@mass[myxlim[1] <= spot@mass & spot@mass < myxlim[2] ]
	yi <-  spot@intensity[myxlim[1] <= spot@mass & spot@mass < myxlim[2] ]

	###############
	#First plot: Raw data	
	message("Calling q2e_plotseqpeaks")
	q2e_plotseqpeaks(cd1,myxlim)
	lines(xm,yi/max(yi),type="l",col= "red")		
	title("Raw data (red), isodists (grey with red dots)", line = "-2")


	###############
	#Second plot: Interpolation
	#create an interpolation (isodists is accurate to 2 decimal places)
	xout = seq(from = myxlim[1], to = myxlim[2], by = myby)

	#sample:
	yii <- approx(x=xm,y=yi,xout=xout, method="linear", rule = 2)
	yii$y = yii$y/max(yii$y)
	plot(yii,type="l",col="red")	
	
	#theoretical: 
	# Map the peaks to an array of x values: 
	yri <- approx(x=cd1$mass,y=cd1$prob,xout=xout, method="linear", rule = 2)
	#set non-peak values to zero
	yri$y[] <-0
	#set the nearest peaks within this:
	for(i in 1:length(cd1$prob)){
		#find the closest x
		idx <- which.min(abs(yri$x-cd1$mass[i]))
		#set the peak value:
		yri$y[idx] <- cd1$prob[i]
	}
	#plot it
	lines(yri)
	title(sprintf("Raw data resampled to resolution %0.4f Da",myby), line = "-2")
	
	###############
	#Third plot: Cross correlation	
	ccd <- ccf(yri$y,yii$y,ylim=c(-0.5,1.0),plot=T,axes=F, lag.max = mylagmax)
		#message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
		#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
	cor = ccd$acf[,,1] 
	lag = ccd$lag[,,1] 
	res = data.frame(cor,lag) 
	res_max = res[which.max(res$cor),] 

	axis(1, at = c(-mybyplotlim,0,mybyplotlim), labels = c(-0.5,0,0.5))
	axis(2)

	td <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f",res_max$cor,res_max$lag*myby)
	message(td)
	title(td, line = "-2")
	
}
###############################################

cd1 <- q2e_isodists(tseq)



myby = 0.005
#TODO: Calculate this from the 'myby' value
mylagmax <- 1/myby




#set up the plot area: (4 rows, 1 column)
par(mar=c(0.9,2.3,0.9,.3), mfrow = c(4,1), oma=c(5,0,2,0))

#Let's do some setting up:
myxlim = c(lbl = floor(min(cd1$mass))-1,ubl = ceiling(max(cd1$mass))+1)

#grab the data for the range we are intested in 
xm <-       testdata[[spot]]@mass[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]
yi <-  testdata[[spot]]@intensity[myxlim[1] <= testdata[[spot]]@mass & testdata[[spot]]@mass < myxlim[2] ]


message("Calling q2e_plotseqpeaks")
q2e_plotseqpeaks(cd1,myxlim)
lines(xm,yi/max(yi),type="l",col= "red")		
title("Raw data (red), isodists (grey with red dots)", line = "-2")

#create an interpolation (isodists is accurate to 2 decimal places)
xout = seq(from = myxlim[1], to = myxlim[2], by = myby)


yii <- approx(x=xm,y=yi,xout=xout, method="linear", rule = 2)
yii$y = yii$y/max(yii$y)
plot(yii,type="l",col="red")


#This doesn't work - the interpolation is all wrong - gives convex hull of peaks, but useful for getting structure:
yri <- approx(x=cd1$mass,y=cd1$prob,xout=xout, method="linear", rule = 2)
#set yvals to zero
yri$y[] <-0
#go through each peak
for(i in 1:length(cd1$prob)){
	idx <- which.min(abs(yri$x-cd1$mass[i]))
	yri$y[idx] <- cd1$prob[i]
}
lines(yri)

title(sprintf("Data resampled to resolution %0.4f Da",myby), line = "-2")

ccd <- ccf(yri$y,yii$y,ylim=c(-0.5,1.0),plot=T,axes=F, lag.max = mylagmax)
	#message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
	#plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )
cor = ccd$acf[,,1] 
lag = ccd$lag[,,1] 
res = data.frame(cor,lag) 
res_max = res[which.max(res$cor),] 

td <- sprintf("Cross-Correlation\n max c=%.2f, at lag=%0.3f",res_max$cor,res_max$lag*myby)
message(td)
title(td, line = "-2")

#TODO: calculate this from the "myby" value
mybyplotlim = 500 #5000 * myby

axis(1, at = c(-mybyplotlim,0,mybyplotlim), labels = c(-0.5,0,0.5))
axis(2)
#text(0,-0.4,td)



#Plot the calculated peaks as probabilities 
q2e_plotseqpeaks(cd1,myxlim)

#PLOT OCOW HERE....

ocow = oIsoDist(tseq)
xx = ocow[,1] + 1 #add 1 to account for the charge
yy = ocow[,3]/100
points(xx, yy, xlim = myxlim)
segments(x0=xx, y0=yy, y1=0, col=8)
points(xx, yy, pch=21, col=3, bg=3)


lines(xm,yi/max(yi),type="l",col= "gray80")
lines(xm+(res_max$lag*myby),yi/max(yi),type="l",col= "green")
title("Shifted spectrum (in green), original (grey)",line = "-2")



message(sprintf("Working on spot %d, sequence %s",spot,tseq))


