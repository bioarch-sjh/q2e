



#' run q2e in the cpp style
#' arg1 is the dataFile
#' arg2 is the peptideFile
#' arg3 is the parameterFile
#' @export
runq2e <- function(fileListFile){
  args = c("arg0","arg1","arg2","arg3")

  #REMEMBER! R arrays are 1-indexed, C arrays are 0-indexed!
  args[2] <- fileListFile
  args[3] <- system.file("extdata", "peptideList", package = "q2e")
  args[4] <- system.file("extdata", "parameterFile", package = "q2e")

  rq2e(args)

}
