

# This file just has a few wrapper functions that call run_boot_analysis
# to do the various analyses that need to be done.



#' Do the steelhead stock analysis in the current working directory
#'
#' Note that you have to be sure that stock_group_start_col is appropriate for the
#' input file.  The input file can be modified to do stock-sex or stock-age if desired.
#' @inheritParams run_boot_gsi_analysis
#' @param STOCK.DATA.XLSX The basename of the file with the data in it for running the simulation
#' The function assumes that the file is in the "data_files" directory of the package. 
#' @export
do_steelhead <- function(STOCK.DATA.XLSX, stock_group_start_col, nsim, B, DO_GSI_ON_PROP) {
  dat.dir = system.file("data_files", package="lowergranite", mustWork=T)
  run_boot_gsi_analysis(
    W = NULL,
    stock_group_start_col = stock_group_start_col,
    DAT.DIR = system.file("data_files", package="lowergranite", mustWork=T),
    WORK.DIR = getwd(),
    STOCK.DATA.XLSX = file.path(dat.dir, STOCK.DATA.XLSX),
    collaps = c(1,1,1,1,1,2,2,2,3,4,5,6,7,8,9,10,10,10,10,11,11,11,11,11,11,11,11),
    DO_GSI_ON_PROP  = DO_GSI_ON_PROP,
    GSISIM = gsi_simBinaryPath(),
    GSI_SEEDS = c(NA, NA),
    BLFILE  = file.path(dat.dir, "sthd_base_v3_187.txt"),
    RUFILE  = file.path(dat.dir, "sthd_base_v3_rg.txt"),
    alph = 0.10,
    B = B,
    nsim = nsim,
    console_messages_to="",
    reset_booty_seed = 0
  )
}



#' Do the chinook stock analysis in the current working directory
#'
#' Note that you have to be sure that stock_group_start_col is appropriate for the
#' input file.  The input file can be modified to do stock-sex or stock-age if desired.
#' @inheritParams run_boot_gsi_analysis
#' @param STOCK.DATA.XLSX The basename of the file with the data in it for running the simulation
#' The function assumes that the file is in the "data_files" directory of the package.
#' @export
do_chinook <- function(STOCK.DATA.XLSX, stock_group_start_col, nsim, B, DO_GSI_ON_PROP) {
  dat.dir = system.file("data_files", package="lowergranite", mustWork=T)
  run_boot_gsi_analysis(
    W = NULL,
    stock_group_start_col = stock_group_start_col,
    DAT.DIR = system.file("data_files", package="lowergranite", mustWork=T),
    WORK.DIR = getwd(),
    STOCK.DATA.XLSX = file.path(dat.dir, STOCK.DATA.XLSX),
    collaps = c(1,1,2,3,4,5,6,7,8,9,10,11,11,11),
    DO_GSI_ON_PROP  = DO_GSI_ON_PROP,
    GSISIM = gsi_simBinaryPath(),
    GSI_SEEDS = c(NA, NA),
    BLFILE  = file.path(dat.dir, "chnk_base_v3_180.txt"),
    RUFILE  = file.path(dat.dir, "chnk_base_v3_repunits.txt"),
    alph = 0.10,
    B = B,
    nsim = nsim,
    console_messages_to="",
    reset_booty_seed = 0
  )
}



