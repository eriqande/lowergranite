
#' given a file path x (or vector of them) return the path with quotes around it
#' 
#' The system path on window has spaces in it which causes problems with 
#' \code{\link{pipe}} and other such system calls.  This just puts those
#' paths inside quotation marks
#' @param x  a file path (or a vector of paths)
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



# reads a gsi_sim file and returns a vector of population names
gsisim2PopsList <- function(bl.file="/Users/eriq/Documents/xp_dev_svn_checkouts/gsi_sim/snpset/2010_SNPset_GSI_TOOLS/Baseline/snpset_Baseline.txt") {
	bf <- readLines(bl.file)
	bf <- bf[grep("^POP[[:blank:]]", bf)]  # pull out just the POP lines
	sapply(strsplit(bf, "[[:blank:]]+"), function(x) x[2])  # now get and return the population names
}	


# reads a gsi_sim reporting units file and returns and list of reporting units.
# They are ordered as they are ordered in the rep-units file.  Names of the components
# are the names of the reporting units and the constituent populations are the 
# contained character vector, in the order they appear in the reporting units file
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


# here is a function to read in our California SNPset baseline file in gsi-sim format and
# turn the contents into a data frame that can be easily handled in R

# bl.file. is the path to the baseline file.
snpset2dataframe <- function(bl.file="/Users/eriq/Documents/xp_dev_svn_checkouts/gsi_sim/snpset/2010_SNPset_GSI_TOOLS/Baseline/snpset_Baseline.txt") {
	bf <- readLines(bl.file)


	# now process that baseline file into a useful format
	bfs <- strsplit(bf, split="[[:space:]]")  # split the lines on white space and make a list of the result
	loci<-unlist(bfs[sapply(bfs, function(x) length(x)==1)])  # get the names of the loci
	dat <- data.frame(t(simplify2array(bfs[sapply(bfs, function(x) length(x) > 100)])), stringsAsFactors=F)  # get the actual data part and save as a data frame
	dat<-cbind(t(simplify2array(strsplit(dat[,1],"[-:]")))[,c(1,3)], dat)  # add columns for reporting unit and population
	dat<-cbind(paste(dat[,1], dat[,2], sep="--"), dat)
	# now we put names on those columns, adding a .1 onto each locus name for the second allele:
	names(dat) <- make.unique(c("RepPop", "RepUnit", "Pop", "ID", rep(loci, each=2)))
	dat$RepPop <- factor(dat$RepPop, levels=unique(dat$RepPop))  # order the populations the  way we want to in the factor


	list(dat=dat, loci=loci)   # return the data frame
}

# make a data frame out of GENEPOP stuff that I processed slightly already
beacham.genepop2dataframe <- function(bl.file="/Users/eriq/Documents/work/prj/BaselinePaper/msat-snp-compare/data/from-terry-beacham/terry_data.txt") {
	bf <- readLines(bl.file)


	# now process that baseline file into a useful format
	bfs <- strsplit(bf, split="\t")  # split the lines on white space and make a list of the result
	bfs <- bfs[!sapply(bfs, function(x) x[1]=="POP")]  # get rid of the POP lines
	loc.logical <- sapply(bfs, function(x) x[2]=="") # these are the locus lines as a logical vector
	loci <- sapply(bfs[loc.logical], function(x) x[1])  # get the vector of locus names
	dat <- data.frame(t(simplify2array(bfs[!loc.logical])))[-4]  # get it as a data frame and drop an empty column (what was a comma)
	# this gets a matrix of numeric data values
	sep.locs <- data.frame(do.call(cbind, lapply(dat[-(1:3)], function(x) cbind(as.numeric(substr(x, 1, 3)), as.numeric(substr(x, 4,6))))))
	
	ff <- cbind(dat[1:3],sep.locs)
	names(ff) <- c("Pop", "Year", "Idx", paste(rep(loci, each=2), c("", ".1"), sep=""))
	
	list(dat=ff, loci=loci)
}




# write a genotype data frame F to GSI_SIM format.
# the gsi_sim part must start at a field named "ID" and proceed to the right with the genotypes
# from there.  ff is the file.  Note that it needs to have a field called RepPop
Frame2gsisim <- function(FR, ff="", POP.SPEC="As.Mixture") {
	IDx <- which(names(FR)=="ID")
	G <- FR[-(1:(IDx-1))]  # drop the extraneous columns for now
	
	write((dim(G) - c(0,1)) / c(1,2), file=ff, append=F) # write the number of inds and loci
	write(names(G[seq(2,ncol(G),by=2)]), file=ff, append=T) # write the locus names
	if(POP.SPEC=="As.Mixture") {
		write("POP Mixture", file=ff, append=T)
		write.table(G, append=T, quote=F, row.names=F, col.names=F, file=ff)
	}
	else {  # here we do things a little more slowly than if we used some other slick function, but we
	        # will get the order right
	  for(pop in unique(FR$RepPop)) {
	  	write(paste("POP", pop), file=ff, append=T)
	  	write.table(G[FR$RepPop==pop,], append=T, quote=F, row.names=F, col.names=F, file=ff)
	  }
		
	}
	
}
