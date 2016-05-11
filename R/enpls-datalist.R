#' Methylalkanes Retention Index Dataset
#'
#' Methylalkanes retention index dataset from Liang et, al.
#'
#' This dataset contains 207 methylalkanes' chromatographic retention index (y) 
#' which have been modeled by 21 molecular descriptors (x).
#'
#' Molecular descriptor types:
#' \itemize{
#' \item Chi path, cluster and path/cluster indices
#' \item Kappa shape indices
#' \item E-state indices
#' \item Molecular electricity distance vector index
#' }
#'
#' @docType data
#' @name alkanes
#' @usage data(alkanes)
#'
#' @format
#' A list with 2 components:
#' \itemize{
#' \item x - data frame with 207 rows (samples) and 21 columns (predictors)
#' \item y - numeric vector of length 207 (response)
#' }
#'
#' @references
#' Yizeng Liang, Dalin Yuan, Qingsong Xu, and Olav Martin Kvalheim. 
#' "Modeling based on subspace orthogonal projections for 
#' QSAR and QSPR research." 
#' Journal of Chemometrics 22, no. 1 (2008): 23--35.
#'
#' @examples
#' data(alkanes)
#' str(alkanes)
NULL


#' The aqueous solubility data
#'
#' The aqueous solubility data reported by Hou et al.
#'
#' This dataset contains 549 aqueous solubility samples (y) which have been modeled
#' by 166 molecular descriptors (x) for training set, 274 aqueous solubility samples 
#' (y) which have been modeled by 166 molecular descriptors (x) for test1 set and 446 
#' aqueous solubility samples (y) which have been modeled by 166 molecular descriptors 
#' (x) for test2 set.
#'
#' Molecular descriptor types:
#' \itemize{
#' \item Constitution
#' \item Topology
#' \item Connectivity
#' \item Kappa
#' \item EState
#' \item Autocorrelation-moran
#' \item Autocorrelation-geary
#' \item Autocorrelation-broto
#' \item Molecular properties
#' \item Charge
#' \item Moe-Type descriptors
#' }
#'
#' @docType data
#' @name logS
#' @usage data(logS)
#'
#' @format
#' A list with 6 components:
#' \itemize{
#' \item x - data frame with 549 rows (samples) and 166 columns (predictors)
#' \item y - numeric vector of length 549 (response)
#' \item x.test1 - data frame with 274 rows (samples) and 166 columns (predictors)
#' \item y.test1 - numeric vector of length 274 (response)
#' \item x.test2 - data frame with 446 rows (samples) and 166 columns (predictors)
#' \item y.test2 - numeric vector of length 446 (response)
#' }
#'
#' @references
#' Hou, T. J.; Xia, K.; Zhang, W.; Xu, X. J. 
#' "ADME Evaluation in Drug Discovery. 4. Prediction of Aqueous Solubility 
#' Based on Atom Contribution Approach" 
#' Journal of chemical information and computer sciences 44, no. 1 (2004): 266-275.
#' 
#' @examples
#' data(logS)
#' str(logS)
NULL
