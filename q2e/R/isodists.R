#' Test the 'R_iso_main' q2e c function
#' @param pepfn peptide file name
#' @keywords test
#' @useDynLib q2e
#' @export
#' @examples
#' isodists()

isodists <- function(pepfn){

	failedflag <- 0

	#load pepfn
	data <- read.table(pepfn,fill=T,header=F,sep= " ")

	npep <- data$V1[1]

	#do some error checking
	if (npep +1 != nrow(data)){
		message("ERROR: First number in file doesn't match the number of peptides in the rest of the file")
		message("exiting isodists")
		return
	}

    message(sprintf("Calculating isodists for %d peptides",npep))

	#create the array for the result (npeptides * 5, 5)
	#result<- matrix(nrow = 5*npep,ncol=2)

	#todo:change the c code to match!

	result <- .C("R_iso_main", infile=as.character(pepfn),mass=as.double(1:(5*npep)),prob=as.double(1:(5*npep)),errflag=as.integer(failedflag))


	message("exiting isodists")
	#result <- c(mass,prob)

	return (result)

}
