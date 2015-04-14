#' Test the 'Rmain' q2e c function
#' @param dfn data file name
#' @param pepfn peptide file name
#' @param paramfn parameter file name
#' @keywords test
#' @useDynLib q2e
#' @export
#' @examples
#' testq2e()



testq2e <- function(dfn,pepfn,paramfn){

.C("Rmain", as.character(dfn), as.character(pepfn), as.character(paramfn))

}
