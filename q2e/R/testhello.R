#' Test the 'hello' c calling function
#' 
#' @param n The number of times to say hello
#' @keywords test
#' @useDynLib q2e
#' @export
#' @examples
#' hello2()


hello2 <- function(n) {
.C("hello", as.integer(n))
}


#' Test string passing in the 'hello person' c function 
#' 
#' @param name The name of the person to say hello to
#' @keywords test
#' @useDynLib q2e
#' @export
#' @examples
#' hello3("Simon")


hello3 <- function(name) {
.C("helloname", as.character(name))
}
