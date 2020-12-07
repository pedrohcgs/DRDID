#' @title Summary
#'
#' @description Summary of a drdid object
#'
#' @param object A drdid object
#' @param ... Other params (required as generic function, but not used)
#'
#' @export
#' @noRd
# Define new summary function
summary.drdid <- function(object, ...){
  drdid.obj <- object
  print(drdid.obj)

}
