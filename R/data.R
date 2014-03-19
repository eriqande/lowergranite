#' A variety of data that is useful for testing functions in the lowergranite package
#'
#' This is a list that holds a variety of quantities that are useful to have
#' to carry out tests of the functions in the package.  Mostly they were
#' created from the steelhead data.
#' 
#'
#' 
#' @docType data
#' @name lg_test_data
#' @usage data(lg_test_data)
#' @format This is a list with the following components:
#' \describe{
#'  \item{RU.list}{An example reporting unit list.}
#'  \item{Originnames}{An example of an Originnames character vector.}
#'  \item{big.gsi.out}{Example gsi_sim standard output (with --multi-fix-mix output) in a character vector.}
#'  \item{multi_fix_df_compare}{Data frame: correct output of \code{\link{extract_multi_fix_sim_to_df}} when run with the
#'  test data above. Allows comparison in tests. This data frame was created with 
#'  \code{expect_that(extract_multi_fix_sim_to_df(lg_test_data$gsi.big.out, lg_test_data$RU.list, lg_test_data$Originnames)}}
#'  \item{known_stock_result15}{output from run_boot_gsi_analysis.R that should look like some test results.}
#'  \item{gsi_stock_result15}{output from run_boot_gsi_analysis.R that should look like some test results.}
#' }
#' @keywords datasets
NULL
