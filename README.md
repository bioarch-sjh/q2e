# q2e

q2e is an R package that uses a genetic algorithm to determine the extent of deamidation that gives the best fit to a measured distribution of masses for glutamine.

## Installing from the R command-line

You'll need the `devtools` package to install q2e so either do this:

 `> require("devtools")`
 
or this:
 
 `> install.packages("devtools")`
 
 `> library(devtools)`
 
Then you can use the `install_github` function to install q2e:
 
 `> install_github("franticspider/q2e")`
 
Of course, you'll only need to do that once for your R installation. After that, whenever you want to use q2e, simply enter: 
 
 `> library("q2e")`


## C Code

Julie Wilson's original C code is in the directory 'jwc'. A modified version of this code is in the q2e 'src' subdirectory

###Compiling, debugging and running

To compile the code, go to the 'src' subdirectory and use the command:

 `R CMD SHLIB r_iso.c isodists.c getline.c`

Then you'll need to load the compiled library in R:

 `> dyn.load("r_iso.so")`
 
Now you can call the C directly: 

 `.C("R_iso_main", infile=as.character(pepfn),mass=as.double(1:(5*npep)),prob=as.double(1:(5*npep)),errflag=as.integer(failedflag))`


##Obsolete functions

- *q2e_isodists*: function has been renamed to *q2e_tpeaks*, as the function returns theoretical peak data from a peptide sequence - nothing to do with isodists


## Academic Referencing

*This software is free, but if you use it in an academic paper, please cite the following:*

Julie Wilson, Nienke L. van Doorn
and Matthew J. Collins (2012). Assessing the extent of bone degradation using glutamine deamidation in collagen.  *Anal. Chem.* 84(21), 9041-9048.

## folders

- *crinkle* contains the original C code developed by Julie Wilson

