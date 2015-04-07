#' Mimic the original q2e program (for testing the individual functions)
#' 
#' @param fitpeaks The number of peaks
#' @param galim The limit value for the genetic algorithm
#' @param snrlim The limit value for the signal to noise ratio
#' @param firstmass The smallest mass value
#' @param lastmass The largest mass value
#' @param csv flag to state whether to produce an ouput in csv format
#' @keywords q2e parameters
#' @export
#' @examples
#' q2eparams()



q2eparams <- function(fitpeaks=4,galim=0.02,snrlim=3.0,firstmass=800.0,lastmass=4000.0,csv=0){

	#TODO: Type checking on set values:

	params <- list(fitpeaks=fitpeaks,galim=galim,snrlim=snrlim,firstmass=firstmass,lastmass=lastmass,csv=csv)
	
	return(params)



}
