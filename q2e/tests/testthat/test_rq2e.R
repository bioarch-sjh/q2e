context("Basic operation of rq2e")


removeFile <- function(fn){
  if(file.exists(fn))
    file.remove(fn)

}


# remove the datafiles that the tests produce
removeOutputs <- function(fnroot){

  removeFile(sprintf("%sBetas.csv",fnroot))
  removeFile(sprintf("%sResults.txt",fnroot))

}



test_that("bad paramFile path is detected",{



  removeOutputs("rq2e")
  datafile <- system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e")

  expect_equal(runq2e(datafile,parameterFile = "badparamfile"),2)

  removeOutputs("rq2e")

})





test_that("A single datafile is parsed",{


  datafile <- system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e")

  expect_equal(runq2e(datafile),0)

  removeOutputs("rq2e")
})



test_that("Julie Wilson's example datafiles are parsed with 3 replicates",{


  datafiles <- c(system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1052a2_65_1_plastic_A1_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1052a3_65_1_plastic_A16_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b1_65_1_plastic_M20_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b2_65_1_plastic_A2_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b3_65_1_plastic_A17_04132012Kinetics1rr.txt", package = "q2e"))

  expect_equal(runq2e(datafiles,nReplicates = 3),0)

  removeOutputs("rq2e")
})



test_that("Julie Wilson's example datafiles are parsed",{
  datafiles <- c(system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1052a2_65_1_plastic_A1_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1052a3_65_1_plastic_A16_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b1_65_1_plastic_M20_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b2_65_1_plastic_A2_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b3_65_1_plastic_A17_04132012Kinetics1rr.txt", package = "q2e"))

  expect_equal(runq2e(datafiles),0)

  removeOutputs("rq2e")
})




