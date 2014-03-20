
#' given a file path x (or vector of them) return the path with quotes around it
#' 
#' The system path on Windows has spaces in it which causes problems with 
#' \code{\link{pipe}} and other such system calls.  This just puts those
#' paths inside quotation marks
#' @param x  a file path (or a vector of paths)
#' @export
quoteProtect <- function(x) {
  paste("\"", x, "\"", sep="")
}

gsi_simBinaryPath <- function() {
  if(.Platform$OS.type=="windows") {
    gspath <- file.path(system.file("bin", package="lowergranite"), "gsi_sim.exe")
  } 
  else {
    if(.Platform$OS.type=="unix") {
      if(Sys.info()["sysname"]=="Darwin") {
        gspath <- file.path(system.file("bin", package="lowergranite"), "gsi_sim")
      } 
      else {
        stop(paste("This appears to be a non-Mac Unix architecture.  Not supported by gpiper currently. Sys.info()[\"sysname\"] =", Sys.info()["sysname"]))
      }
    }
  }
  if(!file.exists(gspath)) stop(paste("gsi_sim executable should be installed at", gspath,"but does not seem to be there"))
  return(gspath)
}



#' reads a gsi_sim file and returns a vector of population names
#' @param path to gsi_sim baseline file
gsisim2PopsList <- function(bl.file) {
	bf <- readLines(bl.file)
	bf <- bf[grep("^POP[[:blank:]]", bf)]  # pull out just the POP lines
	sapply(strsplit(bf, "[[:blank:]]+"), function(x) x[2])  # now get and return the population names
}	


#' reads a gsi_sim reporting units file and returns and list of reporting units
#' 
#' They are ordered as they are ordered in the rep-units file.  Names of the components
#' are the names of the reporting units and the constituent populations are the 
#' contained character vector, in the order they appear in the reporting units file
#' @export
gsisimRepUnits2List <- function(ru.file="/Users/eriq/Documents/xp_dev_svn_checkouts/gsi_sim/snpset/2010_SNPset_GSI_TOOLS/Baseline/snpset_ReportingUnits.txt") {
	ru <- readLines(ru.file)
	ru <- ru[grepl("[[:alnum:]]", ru)] # drop any lines that don't have letters or digits on them
	rsp <- grep("^REPUNIT", ru) # indexes at which the REPUNIT lines are (rsp: rep-start-points)
	names(rsp) <- sapply(strsplit(ru[rsp], "[[:blank:]]+"), function(x) x[2])  # names of the rep units starting at those lines
	beg<-rsp+1  # beginning subscript
	final<-c(rsp[-1],length(ru)+1)-1 # ending subscripts for each reporting unit
	names(final)<-names(beg)
	be <- cbind(beg, final)  # be: beginning and endings
	
	# now get the populations in each one (and strip any white space there)
	ru.list <- lapply(1:nrow(be), function(x) gsub("[[:blank:]]", "", ru[be[x,1]:be[x,2]]))
	names(ru.list) <- names(rsp)
	
	ru.list  # return this
}

