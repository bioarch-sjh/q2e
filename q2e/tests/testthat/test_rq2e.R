context("Basic operation of rq2e")




test_that("bad paramFile path is detected",{



  datafile <- system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e")

  expect_equal(runq2e(datafile,parameterFile = "badparamfile"),2)

})





test_that("A single datafile is parsed",{


  datafile <- system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e")

  expect_equal(runq2e(datafile),0)



})



test_that("Julie Wilson's example datafiles are parsed",{


  datafiles <- c(system.file("extdata", "1052a1_65_1_plastic_M19_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1052a2_65_1_plastic_A1_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1052a3_65_1_plastic_A16_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b1_65_1_plastic_M20_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b2_65_1_plastic_A2_04132012Kinetics1rr.txt", package = "q2e"),
                 system.file("extdata", "1053b3_65_1_plastic_A17_04132012Kinetics1rr.txt", package = "q2e"))

  expect_equal(runq2e(datafiles),0)



})




