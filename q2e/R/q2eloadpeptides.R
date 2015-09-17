#' Create a structure that holds a list of files to be processed
#' 
#' @param listfile name of a file holding the list of files
#' @keywords q2e peptidelist
#' @export
#' @examples
#' q2eloadpeptides()



q2eloadpeptides <- function(listfile,verbose=F){

	#TODO: Type checking on set values:

	x <- read.table(listfile,fill=T)

	ndatafiles = as.integer(as.character(x[1,1]))

	if(ndatafiles != (nrow(x)-1)){
		message("ERROR:\n\tnumber of files specified is ",ndatafiles," but number of files listed is ", (nrow(x)-1) )
	}

	#TODO: Coerce this data into a form suitable for the C input.


	return(x)

}
