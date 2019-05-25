



#' run q2e in the cpp style
#' arg1 is the dataFile
#' arg2 is the peptideFile
#' arg3 is the parameterFile
#' @export
runq2e <- function(fileList,nReplicates = 1, parameterFile = NA, peptideList = NA,verbose = F){
  args = c("arg0","arg1","arg2","arg3")

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
