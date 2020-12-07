#' @title print.matrix1
#'
#' @description Helper function to print a matrix; used by the print methods
#'
#' @param m Some matrix
#'
#' @noRd
#' @importFrom utils write.table
#'

print.matrix.drdid <- function(m){
  utils::write.table(format(m, justify= "centre", digits=2, nsmall=2),
                     row.names=F, col.names=F, quote=F, sep=" ")
}
