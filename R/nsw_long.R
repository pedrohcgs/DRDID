#' @title National Supported Work Demonstration dataset, in long format
#'
#' @description \code{nsw_long} is the same dataset as \code{nsw} but in a long format.
#'
#'
#'
#' @format A data frame in "long" format with 38408 observations on the following and
#' 15 variables:
#' \describe{
#'   \item{id}{unique identifier for each cross-sectional unit (worker).}
#'   \item{year}{year. 1975 is the pre-treatment and 1978 is the post-treatment}
#'   \item{treated}{an indicator variable for treatment status. Missing if not part of the
#'   NSW experimental sample.}
#'   \item{age}{age in years.}
#'   \item{educ}{years of schooling.}
#'   \item{black}{indicator variable for blacks.}
#'   \item{married}{indicator variable for martial status.}
#'   \item{nodegree}{indicator variable for high school diploma.}
#'   \item{dwincl}{indicator variable for inclusion in Dehejia and Wahba sample.
#'   Missing if not part of the experimental sample}
#'   \item{re74}{real earnings in 1974 (pre-treatment).}
#'   \item{hisp}{indicator variable for Hispanics.}
#'   \item{early_ra}{indicator variable for inclusion in the early random assignment
#'   sample in Smith and Todd (2005). Missing if not part of the experimental sample}
#'   \item{sample}{1 if NSW (experimental sample), 2 if CPS comparison group, 3 if PSID
#'    comparison group.}
#'   \item{re}{real earnings (outcome of interest).}
#'   \item{experimental}{1 if in experimental sample, 0 otherwise.}
#'
#' }
#' @source {https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/23407/DYEWLO&version=1.0.}
#' @references
#'   \cite{Diamond, Alexis, and Sekhon, Jasjeet S. (2013),
#'   'Genetic Matching for Estimating Causal Effects: A General Multivariate Matching Method for Achieving Balance
#'    in Observational Studies' Review of Economics and Statistics, vol. 95 , pp. 932-945, \url{https://doi.org/10.1162/REST_a_00318}}
#'
#'
#'   \cite{Smith, Jeffrey, and Todd, Petra (2005),
#'    Does matching overcome LaLonde's critique of nonexperimental estimators?' Journal of Econometrics,
#'    vol. 125, pp. 305-353, \url{https://doi.org/10.1016/j.jeconom.2004.04.011}}
#'
"nsw_long"
