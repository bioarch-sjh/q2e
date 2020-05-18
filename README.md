# q2e

q2e is an R package that uses a genetic algorithm to determine the extent of deamidation that gives the best fit to a measured distribution of masses for glutamine.

# Installing

You'll need the `devtools` package to install q2e so either do this:

 `> require("devtools")`
 
or this:
 
 `> install.packages("devtools")`
 
 `> library(devtools)`
 
Then you can use the `install_github` function to install q2e. As this package is in development, it's best to use the `force` option. Also, the documentation for the package is saved in vignettes, so you probably want to set `build_vignettes=T` so that these are installed:
 
 `> devtools::install_github("https://github.com/bioarch-sjh/q2e",subdir = "q2e",force=T,build_vignettes=T)`
 
Of course, you'll only need to do that once for your R installation. After that, whenever you want to use q2e, simply enter: 
 
 `> library("q2e")`


## Documentation

The documentation (still in development) for q2e is available in vignettes for the package. To see these guides, enter:

`> browseVignettes("q2e")`

This will open a page in your browser that will list the available documentation

## C Code

For reference, Julie Wilson's original C code is in the directory 'jwc'. A modified version of this code is in the q2e 'src' subdirectory, but the modifocations simply allow the Rcpp compiler to access the functions. 


## Academic Referencing

*This software is free, but if you use it in an academic paper, please cite the following:*

Julie Wilson, Nienke L. van Doorn
and Matthew J. Collins (2012). Assessing the extent of bone degradation using glutamine deamidation in collagen.  *Anal. Chem.* 84(21), 9041-9048.

## folders

- *crinkle* contains the original isodists C code developed by Julie Wilson
- *raw* contains the original q2e C code developed by Julie Wilson
