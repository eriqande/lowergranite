

# This file just has a few wrapper functions that call run_boot_analysis
# to do the various analyses that need to be done.



#' Do the steelhead stock analysis in the current working directory
#'
#' @param nsim number of simulations
#' @param B number of bootstrap replicates
#' @param DoGSI True means do the GSI simulation
do_steelhead_stock <- function(nsim, B, DoGSI=FALSE) {
  dat.dir = system.file("data_files", package="lowergranite", mustWork=T)
  run_boot_gsi_analysis(
    W = NULL,
    stock_group_start_col = 14,
    DAT.DIR = system.file("data_files", package="lowergranite", mustWork=T),
    WORK.DIR = getwd(),
    STOCK.DATA.XLSX = file.path(dat.dir, "SH11SIMPOPstock.xlsx"),
    collaps = c(1,1,1,1,1,2,2,2,3,4,5,6,7,8,9,10,10,10,10,11,11,11,11,11,11,11,11),
    DO_GSI_ON_PROP  = DoGSI,
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
#' @param nsim number of simulations
#' @param B number of bootstrap replicates
#' @param DoGSI True means do the GSI simulation
do_chinook_stock <- function(nsim, B, DoGSI=FALSE) {
  dat.dir = system.file("data_files", package="lowergranite", mustWork=T)
  run_boot_gsi_analysis(
    W = NULL,
    stock_group_start_col = 14,
    DAT.DIR = system.file("data_files", package="lowergranite", mustWork=T),
    WORK.DIR = getwd(),
    STOCK.DATA.XLSX = file.path(dat.dir, "CH11SIMPOP.xlsx"),
    collaps = c(1,1,2,3,4,5,6,7,8,9,10,11,11,11),
    DO_GSI_ON_PROP  = DoGSI,
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