



#' run q2e in the cpp style
#' fileList a vector of filenames for the MALDI data
#' arg2 is the peptideFile
#' arg3 is the parameterFile
#' @export
runq2e <- function(fileList,outFile = "rq2e", nReplicates = 1, parameterFile = NA, peptideList = NA,verbose = F){


  args = c("arg0","rq2e","arg2","arg3")

  if(is.character(outFile) && length(outFile) == 1 ){
    args[2] = outFile
  }
  else{
    warning("bad type for outFile, using 'rq2e' instead")
  }

  #REMEMBER! R arrays are 1-indexed, C arrays are 0-indexed!

  if(is.na(peptideList))
    args[3] <- system.file("extdata", "peptideList", package = "q2e")
  else
    args[3] <- peptideList

  if(is.na(parameterFile))
    args[4] <- system.file("extdata", "parameterFile", package = "q2e")
  else
    args[4] <- parameterFile


  if(verbose){
    message("args is:")
    print(args)

    message("fileList is:")
    print(fileList)
  }

  return (rq2e(args,fileList,nReplicates))

}
