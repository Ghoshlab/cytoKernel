# R CMD check does not like the variable used
# by foreach::foreach()
#
#' @importFrom utils globalVariables
   if (getRversion() >= "4.1")
   globalVariables("h")
   globalVariables("l")
