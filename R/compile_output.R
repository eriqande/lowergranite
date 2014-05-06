# a file for functions to extract output from runs done on Kanaloa



#' Read all the results from all the runs with output in Dirs into a named list
#' 
#' @param Dirs  a character vector of directories to slurp results from. Ideally 
#' you will just give the basenae of them all as these will become the names
#' of the output list.
#' @export
#' @examples
#' \dontrun{
#' setwd("~/prj/LowerGranite/runs/lowergranite_big_run_CincoDeMayo_2014/")
#' big_run <- slurp_all_results(dir(pattern = "^[CS].*"))
#' }
slurp_all_results <- function(Dirs) {
  
  get_results <- function(dd) {
    ff <- file.path(dd, "results.rda")  # this is the file to get
    suppressWarnings(rm(results))  # make sure this is removed from the current environment
    load(ff)
    if(!exists("results")) {
      warning(paste("No results object found in: ", ff))
      return(NULL)
    } else {
      load(ff)
      return(results)
    }
  }
  
  ret <- lapply(Dirs, function(x) get_results(x))
  names(ret) <- Dirs
  ret
}