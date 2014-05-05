



#' return a population drawn at random from those in a reporting unit
#' 
#' pass this a vector of reporting group names and a list of the 
#' populations in each of those, and this creates a vector of populations
#' sampled uniformly from within each rg
#' @param ru.list  list whose names are reporting units and components are vectors of population names
#' @param sampled_rgs a vector of reporting units (denoting the stock of each fish for example)
#' @export
get_pop_for_rg <- function(ru.list, sampled_rgs) {
  res <- rep(NA, length(sampled_rgs))
  tab <- table(sampled_rgs) # count up the occurrences
  tab <- tab[tab>0] # only keep the non-zero ones
  pops_list <- lapply(names(tab), function(x) {
    if(is.null(ru.list[[x]])) stop(paste("Big problems my friend.", x, "is not a valid reporting unit."))
    sample(ru.list[[x]], tab[x], replace=T)
    })
  names(pops_list) <- names(tab)
  for(i in names(tab)) {
    res[sampled_rgs==i] <- pops_list[[i]]
  }
  res
}



#' modify the Sim.List so that Prop are gsi assignments
#' 
#' @export
gsi_ize_the_Sim.List <- function(Sim.List, ru.list, GSISIM, BLFILE, Originnames, BL.pops) {
  
  # simulate a population of origin for each fish.  If doing stock-sex of stock-age with the ".."
  # column naming convention, then x$Prop$Groop will have stock-sex or stock-age, but we want
  # just plain stock for this section. Ultimately we will want to tweeze it off like this:
  # sub("\\.\\..*", "", x$Prop$Groop)
  # but I am going to commit my other stuff first.
  prop_with_pops <- lapply(Sim.List, function(x) cbind(x$Prop, gsi_pop = get_pop_for_rg(ru.list, sub("\\.\\..*", "", x$Prop$Groop))))  
  
  # now split each component of that on Stratum, and name the result rows of each data from 1:nrow
  # and then order it by the population of origin that will be chosen for each fish, where the ordering
  # of those pops is as they appear in the baseline file, which is BL.pops.  We have to order them like
  # this so that they will line up with the output from gsi_sim.  We use the rownames to re-order them
  # back into the order they were in in Sim.List.
  pwp_split_on_stratum <- lapply(prop_with_pops, function(x) {
                             y<-split(x, x$Stratum); lapply(y, function(z) {rownames(z)<-1:nrow(z); z[order(factor(z$gsi_pop, levels=BL.pops)), ]})
                           })
  
  # now we need to count up the number of fish from each population in there, and we want to include 0's for
  # populations that have no 0's and we want to put those in the order in which they appear in the baseline file
  # and we want to name them X_Y where X is the simulation number and Y is the stratum
  
  # first, get a list of matrices with stratum as the rowname and pop as the column
  tmp <- lapply(pwp_split_on_stratum, function(x) t(sapply(x, function(y) table(factor(y$gsi_pop, levels=BL.pops)))))
  
  # then add the simulation number to each rowname so it is like "X_Y",
  tmp2 <- lapply(1:length(tmp), function(x) {z<-tmp[[x]]; rownames(z) <- paste(x, rownames(z), sep="_"); z})
  
  # then rbind those
  gsi_sim_nums <- do.call(rbind, tmp2)
  
  # and now we can write gsi_sim multi-fix-mix commands from those
  unlink("gsi_multi_fix_calls.txt")  # start by removing the comman file if it is there already
  for(i in 1:nrow(gsi_sim_nums)) {
    cat("--multi-fix-mix", rownames(gsi_sim_nums)[i], as.integer(gsi_sim_nums[i, ]), "\n", file="gsi_multi_fix_calls.txt", append=T)
  }
  
  # now compile up the text for the arguments for the gsi_sim call:
  gsi.args <- c("-b", quoteProtect(BLFILE), "--command-file gsi_multi_fix_calls.txt") 
  
  # and then make a system call to gsi_sim and capture the standard output
  # note that if we wanted to condense on the fly in unix we could pipe it through this:
  # awk '/^MULTI_FIX_MIX_MIXFISH_NAMES_HEADER:/ && got_it==0 {print; got_it=1; fields=NF} /^MULTI_FIX_MIX_MIXED_FISH_INDIVS:/ {for(i=1;i<=fields;i++) printf("%s ", $i); printf("\n")}'
  # but for the PC folks, we will just do all this in R.  It is a little slower but not too bad
  suppressWarnings(gsi.big.out <- system2(GSISIM, args=gsi.args, stdout=T))  # suppressing warnings about chopping long lines. See notes.


  # now extract the gsi simulation results:
  gsi_output_df <- extract_multi_fix_sim_to_df(gsi.big.out, ru.list, Originnames)
  
  # Now we put all that info into a list structure to be transferred back
  # to the pwp_split_on_stratum and ultimately to the Sim.List
  gsi_output_list <- lapply(split(gsi_output_df, gsi_output_df$SimRep.GSI), function(x) split(x, x$Stratum.GSI))
  
  merged_list <- lapply(1:length(pwp_split_on_stratum), 
    function(x) {
      lapply(1:length(pwp_split_on_stratum[[x]]), 
        function(y) {
          z <- cbind(pwp_split_on_stratum[[x]][[y]], gsi_output_list[[x]][[y]])  # pwp must go first to preserve its rownames
          
          # here we deal with attaching the observed sex or age of the fish to the GSI assignment
          # if the stocknames have a ".." in them brokenNames is an array with 2 rows, the second
          # of which is the Age (or sex, etc)-determined category of the fish.  We will paste
          # that only the GSI-inferred group and put that into Groop.GSI. 
          brokenNames <- simplify2array(strsplit(as.character(z$Groop), split = "\\.\\.")) 
          
          if(is.null(dim(brokenNames))) {  # if there is no ".." in the names, just copy across Groop.GSI.As.String and make it a factor using the Originnames
            z$Groop.GSI <- factor(z$Groop.GSI.As.String, levels=Originnames)
          } else if(nrow(brokenNames)==2)  { # in this case we did have the ".." naming convention and we need to deal with the observed Ages/Sexes/Etc.
              z$Groop.GSI <- factor(paste(z$Groop.GSI.As.String, brokenNames[2, ], sep=".."), levels=Originnames)
          } else {
            stop("Wrong number of \"..\"\'s in stock names.  Should be at most one \"..\" per name.")
          }

          
          
          z <- z[as.character(1:nrow(z)), ]  # put it back into its original order
          z$Groop.Truth <- z$Groop  # add a column to preserve the truth
          z$Groop <- z$Groop.GSI   # set Groop to what you got with GSI, to be used in the downstream estimates
          z
        }
      )}) 
  
  merged_list
}


#' read in gsi_sim ouput and extract the multi-fix-mix lines into a data frame
#' 
#' takes a character vector of gsi_sim output from readLines, extracts the multi-fix-mix 
#' output lines and gets the header that
#' it needs to make a huge data frame of them all, then, using the popnames
#' and the ru.list, finds the MLE reporting group for each individual.
#' @param gsi.big.out  A character vector of gsi_sim's standard output.  
#' @export
#' 
extract_multi_fix_sim_to_df <- function(gsi.big.out, ru.list, Originnames) {
  # grab the names header, split it and remove there is our header 
  if(length(grep("^MULTI_FIX_MIX_MIXFISH_NAMES_HEADER:", gsi.big.out))==0) stop("Apparently didn't find any ^MULTI_FIX_MIX_MIXFISH_NAMES_HEADER: in the gsi_sim output file.")
  header <- strsplit(gsi.big.out[min(grep("^MULTI_FIX_MIX_MIXFISH_NAMES_HEADER:", gsi.big.out))], split="  *")[[1]]
  
  
  num_fields <- length(header)
  
  # then get the lines of output
  if(length(grep("^MULTI_FIX_MIX_MIXED_FISH_INDIVS:", gsi.big.out))==0) stop("Didn't find any ^MULTI_FIX_MIX_MIXED_FISH_INDIVS: lines in gsi.big.out")
  outlines <- t(sapply(strsplit(gsi.big.out[grep("^MULTI_FIX_MIX_MIXED_FISH_INDIVS:", gsi.big.out)], split="  *"), function(x) x[1:num_fields]))
  
  # now, get the MultMixFish column to use as a factor for splitting it up later
  MMF.id <- outlines[,2]
  FromPop <- outlines[,6] # and get the pop the fish was simulated from
  # now drop those and other not so useful fields from the header and the output lines
  drop <- c(1, 2, 3, 4, 5, 6)
  header <- header[-drop]
  outlines <- outlines[, -drop]
  
  # check at this point that we have the right number of columns
  if(ncol(outlines) != length(unlist(ru.list))) stop(paste("Wrong number of columns,",ncol(outlines), "in outlines. should be", length(unlist(ru.list))))
  
  # also check to make sure the names match
  if(length(union(setdiff(unlist(ru.list), header), setdiff(header, unlist(ru.list)) ) ) > 0) stop("Pop name mismatch between header and ru.list")

  # now, make outlines a numeric matrix and put colnames on it
  class(outlines) <- "numeric"
  colnames(outlines) <- header
  
  # here is the RG each fish is assigned to.
  # note that I make it a factor with explicit levels to be sure it is commensurate with the results in
  # Sim.List
  rg_by_gsi <- factor(names(ru.list)[apply(sapply(ru.list, function(x) rowSums(cbind(outlines[,x]))), 1, which.max)], levels=Originnames)
  
  # here we see about making the same thing, but we just store them as strings, as we will have
  # to muck with them later for stock-sex and stock-age...
  rg_strings_by_gsi <- names(ru.list)[apply(sapply(ru.list, function(x) rowSums(cbind(outlines[,x]))), 1, which.max)]
  
  # get the Simulation replicate and the Stratum for each 
  sr_and_strat <- matrix(as.numeric(unlist(strsplit(MMF.id, "_"))), byrow=T, ncol=2)
  colnames(sr_and_strat) <- c("SimRep.GSI", "Stratum.GSI")
  
  # and now send it back as a data frame
  data.frame(sr_and_strat, Groop.GSI=rg_by_gsi, Pop.Simmed.From=FromPop, Groop.GSI.As.String=rg_strings_by_gsi)
}



#' create fake "super-informative-data" for a set of baseline pops
#' 
#' Output is written as a gsi_sim file.
#' 
#' @param BLFILE  a baseline file in gsi_sim format that has all the population names
#' that you want in the order that you want them
#' @param outfile  path to where you want the output written
#' @export
make_super_informo_data <- function(BLFILE, outfile="super_informo_baseline.txt", nperpop=100) {
  blr <- readLines(BLFILE) 
  pops <- matrix(unlist(strsplit(x=blr[grep("^POP", blr)], split="\\s+", perl=T)), nrow=2)[2,]
  names(pops) <- pops
  npops <- length(pops)
  
  # now start writing.  there will be nperpop in each pop and one locus for each pop is 
  # fixed for an alternate allele
  cat(c(100*npops, npops, "\n"), file=outfile)
  cat(paste("Locus", 1:npops, sep="_"), sep="\n", file=outfile, append=T)
  
  # lapply over a bunch of tables
  ret <- lapply(1:npops, function(z) {
      x <- names(pops)[z]
      y<-pops[[z]]
      cat(c("POP", x, "\n"), file=outfile, append=T)
      mat <- matrix(1, ncol=2*npops, nrow=nperpop)
      mat[,c(2*z-1, 2*z)] <- 2  # here is their special fixed locus
      rownames(mat) <- paste(x,1:nperpop, sep="_")
      write.table(mat, quote=F, row.names=T, col.names=F, sep=" ", file=outfile, append=T)
  })
  
}



#' Test that all STOCK..CHARACTERISTIC are present.
#' 
#' When using the ".." notation to define subroups of stocks (for example stock-by-sex or
#' stock-by-age) it is required that any subgroup designations be applied to all stocks.
#' This function just checks that this is the case.  If not, it throws an error.  Of course, if none of the
#' stock group designations have a ".." in them, then there is no problem.
#' @param Onames  the names of the different origin groups.  This will be Originnames, typically.
#' @export
check_subgroup_level_completeness <- function(Onames) {
  if(any(grepl("\\.\\.", Onames))) {
    no_dots <- Onames[!grepl("\\.\\.", Onames)]  # check if there are any missing the two dots
    if(length(no_dots) > 0) stop(paste("Missing double dots in stock-name:", no_dots, "\n  "))
    
    df <- as.data.frame(matrix(unlist(strsplit(Onames, "\\.\\.")), ncol=2, byrow=T))
    df$V1 <- factor(df$V1, levels=unique(df$V1)) # put them in the right order
    df$V2 <- factor(df$V2, levels=unique(df$V2)) # put them in the right order
    
    name_counts <- as.data.frame(table(df$V1, df$V2))
    missing <- name_counts[name_counts$Freq==0, ]
    if(nrow(missing) > 0) {
      miss_str <- paste(missing$Var1, missing$Var2, sep="..", collapse="  ")
      stop("Missing certain stock-by-subgroup levels: ", miss_str)
    }
  }
}