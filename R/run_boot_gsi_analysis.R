# Main function for the package is here.


#' Run the bootstrap analysis, with gsi module if specified in parameters
#' 
#' This is basically Kirk's script that eric has wrapped up into a function, and made 
#' some significant changes to.  Notably, the interface to gsi_sim has been rewritten
#' to make it much more time efficient.  The simulated "true" data sets are all simulated
#' first and stored in a list, and then, if doing the gsi_sim part, all of that is passed
#' to gsi_sim in one fell swoop and the Groop of each fish is replaced by its "gsi-inferred"
#' Groop.
#' 
#' Note that the default values here are set up to do an analysis of the steelhead data,
#' and, by default, to not use gsi_sim assigments.
#' 
#' @note Running gsi_sim this produced warnings that look like this: 
#' \code{Warning message:}
#' \code{In system2(GSISIM, args = gsi.args, stdout = T) :} 
#' \code{line 98 may be truncated in call to system(, intern = TRUE)}.  
#' This was \emph{not} a problem.
#' The cause of this is that the output file from gsi_sim has a line it
#' that shows that the
#' command line looked like after all the --multi-fix-mix commands 
#' were stuck onto it using
#' gsi_sim's --command-file option.  That line is so large that system() 
#' chops it into two
#' (or many, if you are doing large nsim).  This line is 
#' far away from any of the important
#' lines we will grep out, however, so it is not of any consequence.
#' I wrapped it in a \code{\link{suppressWarnings}} to not bark too much.
#' 
#' Also note that this function still writes stuff out to an xlsx file of a name
#' that is not specified in the parameters yet.  It will probably be better 
#' to just return it as a data frame or matrix anyway.
#' @return returns the \code{sumrys} object from Kirk's script.
#' @param W a data frame that holds all the stock data in it used to drive the simulation.  This is the data
#' frame that we typically get by reading in the .xlsx file named in STOCK.DATA.XLSX.  However, we can also
#' just pass a data frame in directly.  By default W is NULL, in which case the stock data will be obtained
#' by reading in the STOCK.DATA.XLSX file.  If W is not NULL, then W will be used instead
#' of STOCK.DATA.XLSX.  If both W and STOCK.DATA.XLSX are NULL, that is an error.
#' @param stock_group_start_col The column at which the stock groups start in the xlsx file (or, equivalently,
#' in the data frame W).  That data frame, (whether it was passed in as W or read from an xlsx file) must have column
#' names that correspond to the groupings of stocks that you want to do the bootstrapping for.  (For example,
#' "UPSALM", "MFSALM", "SFSALM", ... ).  The column at which those names start must be passed to this function as
#' \code{stock_group_start_col}. You must specify a value for this.  There is no default.  \emph{Note:} there must be no other columns in the data frame after the stock group
#' columns.
#' @param DAT.DIR  The directory where all the data files are.  Defaults to the directory
#' "data_files" in the installed package
#' @param WORK.DIR  The working directory to do this in.  Default = current working directory.
#' Note that gsi_sim will also be run in this directory.
#' @param STOCK.DATA.XLSX the path of the file that has the stock data in it used to drive the simulations. This can be NULL, in which
#' case parameter \code{W} must be specified. 
#' @param drop.these.groups  A character vector of the names of the stocks or stock-by-age or stock-by-sex groups that you want
#' to drop from the analysis (typically because they are at such low numbers that there are bootstrap reps when none of
#' them occur).  Must follow the convention of column naming in the file.  For example \code{c("UPSALM..BY04", "MFSALM..BY08")}.
#' If this is non-null then only columns not matching any of the entries in this character vector will be retained.
#' @param collaps A vector of numbers in 1,...,N telling which weeks should be lumped together into "statistical weeks" from Kirk's code. It looks 
#' like this gets used in a lot of the bootstrapping functions, but is not a formal parameter of the bootstrapping functions.
#' @param DO_GSI_ON_PROP if set to TRUE then gsi_sim is used to create assignments that replace the assignments in the variable Prop.  If FALSE
#' then the true origins are used.
#' @param GSISIM path to the gsi_sim executable. 
#' @param GSI_SEEDS vector of two positive integers that will be written to the 
#' gsi_sim_seeds file to make reproducible results. If any elements of the vector
#' are NA, then gsi_sim_seeds is not modified.
#' @param BLFILE path to the gsi_sim baseline file.
#' @param RUFILE path to the gsi_sim reporting group file. 
#' @param alph One minus the size of the desired confidence interval to be calculated
#' @param B number of bootstrap replicates. The default is 5 --- much lower than it should be 
#' (should be more like 500)
#' because it takes a long time and this is better for testing.
#' @param nsim number of simulation replicates to do. The default is 2 --- much lower than it should be 
#' (should be more like 500)
#' because it takes a long time and this is better for testing
#' @param console_messages_to  path to a file you want the console messages written to.  Note that it will
#' always append these to a file.  Default is "" which means send it to the console.
#' @param reset_booty_seed for some reason, calling gsi_sim seems to get the random number generator out
#' of state in a way that cannot be restored by saving .Random.seed and then setting that value back to 
#' itself.  It is odd and vexing.  Anyway, in order to test that comparable results are obtained with
#' the super-informo gsi data, we have this.  It should be an integer.  If >0 then it will be passed
#' to \code{\link{set.seed()}} after the gsi code has been run (before entering the bootstrap loop.)
#' For an example of its use, see the test files.  
#' @export
#' @examples 
#' # Do a very short run with known stock of origin:
#' set.seed(5)
#' known_stock_result1 <- run_boot_gsi_analysis(stock_group_start_col = 14, nsim = 10, B = 50, DO_GSI_ON_PROP = F)
#' 
#' # Do a short run using the gsi assignments
#' gsi_result1 <- run_boot_gsi_analysis(stock_group_start_col = 14, nsim = 5, B = 10, DO_GSI_ON_PROP = T)
#' 
run_boot_gsi_analysis <- function(
  W=NULL,
  stock_group_start_col,
	DAT.DIR = system.file("data_files", package="lowergranite", mustWork=T),
  WORK.DIR = getwd(),
  STOCK.DATA.XLSX = file.path(DAT.DIR, "SH11SIMPOP_StockSex.xlsx"),
  drop.these.groups = NULL,
	collaps = c(1,1,1,1,1,2,2,2,3,4,5,6,7,8,9,10,10,10,10,11,11,11,11,11,11,11,11),
	DO_GSI_ON_PROP  = FALSE,
  GSISIM = gsi_simBinaryPath(),
  GSI_SEEDS = c(NA, NA),
  BLFILE  = file.path(DAT.DIR, "sthd_base_v3_187.txt"),
  RUFILE  = file.path(DAT.DIR, "sthd_base_v3_rg.txt"),
  alph = 0.10,
  B = 5,
  nsim = 2,
	console_messages_to="",
	reset_booty_seed = 0, 
  GroupMin = 0,
  Run = "2011 Steelhead Stock"
) {


  
  if(is.null(STOCK.DATA.XLSX) && is.null(W)) stop("At least one of STOCK.DATA.XLSX or W must be non-NULL")
  
  # Start by setting the working directory.
  setwd(WORK.DIR)
  
  
  ############ KIRK'S SCRIPT SET UP TO RUN WITH GSI ASSIGNMENTS  ##################
  ##
  ## Modified by MA on 2-6-14 using code provided by Eric Anderson on 8/7/13
  ##
  ## THIS SECTION HAS VALUES THAT NEED TO BE SET FOR DIFFERENT DATA SETS, ETC
  ## NOTE THAT ADDITIONAL MACHINATIONS ARE NEEDED FOR SETTING STOCKS FOR CHINOOK vs STHD

  #############  PREAMBLE ADDITIONS TO GATHER GSI BASELINE FILES AND REPORTING UNITS, ETC ############ 
  DO_GSI <- DO_GSI_ON_PROP   # this will do GSI stuff that must be done for either GSI method
  if(DO_GSI==TRUE)  {
  
    # get the populations in the order they appear in the baseline file:
    BL.pops <- gsisim2PopsList(BLFILE)
  
    # and make a vector of each pop's index in the baseline file:
    BL.pop.idx <- 1:length(BL.pops)
    names(BL.pop.idx) <- BL.pops
  
    # get the reporting units as a named list of character vectors of population names
    RU.list <- gsisimRepUnits2List(RUFILE)
  
    # check to make sure all pops are present in baseline file:
    tidx <- BL.pop.idx[unlist(RU.list)]  # temp index vector
    if(any(is.na(tidx))) {
      stop("Did not find pops in baseline file for pops in rep-unit file: ", paste(unlist(RU.list)[is.na(tidx)], collapse=" "))
    }
    # deal with gsi_sim seeds if need be;
    if(!any(is.na(GSI_SEEDS))) {
      gseeds <- as.integer(GSI_SEEDS)
      if(any(gseeds<=0)) stop("gsi_sim seeds must be positive")
      if(length(gseeds)!=2) stop("gsi_sim seeds must be an integer vector of length 2")
      if(!any(is.na(GSI_SEEDS))) {
        cat(gseeds, file="gsi_sim_seeds")
      }
    }
  
  }  # close if(DO_GSI==TRUE)
  ####################################################################################################

  

  ######################################################################
  # Define function to find weighted wild by week
   WBW <- function(trp){

  # Calculates and returns wild estimates by statistical week
     Wildweek <-  numeric(nrow(windata))
     for ( h in as.vector(windata[,1])) {
       Wildweek[h] <- trp[h]*windata[h,2]
     }
     return( Wildweek)
   }

  ################# Bootstrap confidence interval function #######################
  # Define the parametric bootstrap confidence interval function
  # Inputs are the estimated proportions of wild by statistical week
  # and estimated proportions for each brood year by collapsed statistical week
  BootLGRrun <- function(Trp,Pr) {

  # Set up storage for bootstrap results

    theta.b <- matrix(numeric(p*B),ncol=p)
    CI <- matrix(numeric(p*3),ncol=3)

  ################################################################################
  # bootstrap loop  -- this is a parametric weighted, weighted/collapsed version
    for (b in 1:B) {


      tstar <- rbinom(length(nT), size = nT, prob = Trp) # These are number of wild fish each week
      nHstar <- round( tstar*W$PWhandled ) # This is the number of bootstrap wild fish handled each week
      tstar <- tstar/nT # Convert binomial count into a proportion

      nHcollaps <- mApply(nHstar,collaps,sum)
      pstar <- matrix(numeric(nCollaps*nGrps),ncol=nGrps)
      for ( h in 1:nCollaps )
        pstar[h,] <- rmultinom(1,nHcollaps[h],Pr[h,])/nHcollaps[h]   # These are the proportions by sex/age year for this stratum

  # Calculate wild by week for the bootstrap data
      WiByWe <- WBW(tstar)  # This is wild by week for all statistical weeks
      WeWi <- mApply(WiByWe,collaps,sum)  # This is the wild by week collapsed for origin calculation
      OW <- WeWi %*% pstar  # Bootstrap numbers of wild fish by origin
      theta.b[b,] <- c(sum(WiByWe),OW)
    } # End of bootstrap loop

    theta.grps <- theta.b[,-1]
      
    # Find one-at-a-time confidence intervals for each statistic
    CI <- matrix(numeric(p*2),ncol=2)
    for  (j in 1:p) {
      CI[j,] <- quantile(theta.b[,j],c(alph/2,1-alph/2))
    }
    CI <- round(CI)
    TotalWildCI <- CI[1,]
    OneCI <- CI[2:p,]
  
  # Apr 26 '14
  # Before finding joint confidence intervals omit group columns with few members
  # based on the average of the bootstrap estimates < GroupMin
  # First set all joint CIs to NA 
    MBCI <- matrix(NA,nrow = nGrps,ncol=2) 
    theta.grps <- theta.b[,-1]                   # Delete total wild column
    mn <- apply(theta.grps,2,mean)
    keap <- (1:ncol(theta.grps))[mn>GroupMin]
    if ( length(keap) > 1) {                     # Do joint CI only if 2 or more groups remain
      theta.grps <- theta.grps[,keap]
  
  # Find simultaneous rectangular confidence intervals via Mandel/Betensky 2008
      temp <- SCSrank(theta.grps,conf.level=1-alph)
      temp <- round(temp$conf.int)
      MBCI[keap,] <- temp
    }
    return( list(TotalWildCI,OneCI,MBCI) )
  
  }
  #################### End of bootstrap CI function ##############################


  #################### MAIN  #####################################################

  # Define basic run parameters

  cat("\nStart time: ",date(),"\n", file=console_messages_to, append=T)
  cat("\nThis is a run of ", Run, "\n", file=console_messages_to, append=T)
  cat("\nParameter input file is ",STOCK.DATA.XLSX,"\n", file=console_messages_to, append=T)
  cat("\nParametric bootstrap: B = ",B,"Number Simulations = ",nsim,"\n", file=console_messages_to, append=T)
  cat("\nAlpha =",alph,"\n", file=console_messages_to, append=T)
  cat("\nMinimum group size for inclusion in joint CIs is ",GroupMin,"\n", file=console_messages_to, append=T)
  
  # Read in weekly counts and define the true POPULATION
  if(is.null(W)) {  # get this as the data frame already passed in, or, if not, read it from the xlsx file
    W <- read.xlsx(STOCK.DATA.XLSX,1)
  }
  
  # drop any of the populations in drop.these.groups
  if(!is.null(drop.these.groups)) {
    W <- W[!(names(W) %in% drop.these.groups)]
  }
  
  # here we pull the originnames out of the file
  Originnames <- names(W)[stock_group_start_col:ncol(W)]

  W$SimPop <- round(W$SimPop)
  TrueWild <- sum(W$PopWild)
  TrueOrigin <- as.matrix(W[Originnames]) 
  nGrps <- length(Originnames)
  nw <- length(unique(W$Stratum))
  Probab <- matrix(numeric(nw*nGrps),ncol=nGrps)
  for (i in 1:nw)
  Probab[i,] <- TrueOrigin[i,] 
  Truth <- Probab*W$PopWild
  TrueGroups <- mApply(t(Truth),1:nGrps,sum)
  names(TrueGroups) <- Originnames
  windata <- data.frame(W$Stratum,W$SimPop)
  windata

  p <- nGrps + 1 
  cols <- p*3 + 2*(p-1)
  simstf <- matrix(numeric(nsim*cols),ncol= cols)
  cn <- c("Total","TL","TU")
  for (ii in 1:nGrps)  cn <- c(cn,Originnames[ii],paste(Originnames[ii],"L",sep=""),paste(Originnames[ii],"U",sep=""))
  for (ii in 1:nGrps)  cn <- c(cn,paste(Originnames[ii],"Ls",sep=""),paste(Originnames[ii],"Us",sep=""))
  colnames(simstf) <- cn

  
  # eric builds up this whole list first, so that he can change stock of origin
  # of all the fish in one fell swoop in gsi_sim
	Sim.List <- lapply(1:nsim, function(x) {
        #  generate random trap and female proportions-- by week
        nT <- numeric(nw)
        nW <- numeric(nw)
        nH <- numeric(nw)  
        for (ii in 1:nw) {
          nT[ii] <- rbinom(1,W$SimPop[ii],W$PTrp[ii])  # Sample 4 to 12 percent of the run.  This is the number of fish trapped.
          tStrat <- rep(ii,nT[ii]) # Define the trap strata column
          nW[ii] <- rbinom(1,nT[ii],W$Pwild[ii])   # This is the number of trapped fish that are wild.
          grp <- c(rep("WILD",nW[ii]),rep("H",nT[ii]-nW[ii]))  # Define the column of WILD and H(atchery) data
          nH[ii] <- round(nW[ii]*W$PWhandled[ii]) # This is the number of trapped fish sexed or aged or sampled for genetic analysis
          Stratum <- rep(ii,nH[ii])  # Define the stratum column for the proportion data
          Groups <- rmultinom(1,nH[ii],Probab[ii,])
          Groop <- "   " # This is just a programming device to avoid an if statement.
          for (jj in 1:nGrps) {
             if(Groups[jj] >0) Groop <- c(Groop,rep(Originnames[jj],Groups[jj]))}
             Groop <- Groop[2:length(Groop)] # Now delete the placeholder entry
          if (ii == 1) {
             Trap <- data.frame(tStrat,grp)
             Prop <- data.frame(Stratum,Groop) }
          else { 
            Trap <- rbind(Trap,data.frame(tStrat,grp))
            Prop <-rbind(Prop,data.frame(Stratum,Groop))  
          }
        }
        list(Prop=Prop, Trap=Trap, nT=nT)  # return this into a list  
      }
    )  # end of the lapply function
  
  
  
  if(DO_GSI==TRUE)  {
    # this is just calling the function here, but not incorporating the changes in the real Sim.List
    Sim.List.Orig <- Sim.List  # just keeping this around for good measure
    # get gsi assignments for each simulated fish
    gsi.list <-  gsi_ize_the_Sim.List(Sim.List, RU.list, GSISIM, BLFILE, Originnames, BL.pops)
    # and then copy those back to the Sim.List structure
    Modified <- lapply(1:length(gsi.list), function(x) {y<-Sim.List[[x]]; y$Prop<-do.call(rbind, args=gsi.list[[x]]); y})
    Sim.List <- Modified  # this is where we put the gsi assignments back in place of the "truth"
  }

   if(reset_booty_seed > 0) {
      set.seed(reset_booty_seed) 
   }
  
  
    # Convert the random trap and proportion data sets to wild and sex/origin proportions
    # here is what needs to get passed below this point:
    # Trap, Prop 
    # nGrps, 
    
  for (ss in 1:nsim) {  # ECA put this simulation loop here which pulls stuff out from 
                       # the results from the last lapply.
  
  	Trap <- Sim.List[[ss]]$Trap   # here is where we pull out the results from the above lapply
  	Prop <- Sim.List[[ss]]$Prop
  	nT   <- Sim.List[[ss]]$nT
  	
    ntrap <- nrow(Trap)
    nprop <- nrow(Prop)
    Wtable <- table(Trap)
    
    
    Wprop <- Wtable/apply(Wtable,1,sum)  # This is the proportion hatchery and wild
    
    # Eric replaced Kirk's kluge here by forcing empty factors to get counted too
    Ptable <- table(factor(Prop$Stratum, levels=1:length(unique(W$Stratum))), factor(Prop$Groop, levels=Originnames))  
  
  # Although the proportion data are generated for each statistical week, the data are collapsed to insure adequate sample sizes

    # here we can test to make sure that collaps makes sense:
    if(length(collaps) != length(unique(W$Stratum))) stop("collaps should have length equal to the number of strata")
    missy <- setdiff(range(collaps)[1]:range(collaps)[2], collaps)
    if(length(missy)>0) stop(paste("collapse is not continuous.  Seems to be missing", paste(missy, collapse=", ") ))
    
  
  
    collapsed <- mApply(Ptable[,1],collaps,sum)
    for ( jj in 2:nGrps ) collapsed <- cbind(collapsed,mApply(Ptable[,jj],collaps,sum))
    Pprop <- collapsed/apply(collapsed,1,sum)
    nCollaps <- nrow(Pprop)

  # Calculate wild by week, then by origin
    WildByWeek <- WBW(Wprop[,2])
    WeeklyWild <- mApply(WildByWeek,collaps,sum)
    GroupWild <- WeeklyWild %*% Pprop
  # Put together estimate of total wild and wild by group
    ests <- c(sum(WildByWeek),GroupWild)

    # Call bootstrap routine to find one-at-a-time and simultaneous confidence intervals
    CIs <- BootLGRrun( Wprop[,2],Pprop )
  
    # Now put simulation results into rows of an array for later analysis
    TotalWild <- c(ests[1],CIs[[1]])
    simstf[ss,1:3] <- TotalWild
    single <- cbind(ests[2:p],CIs[[2]])
    joint <- CIs[[3]][1:p-1,]
    for ( jj in 2:p ) simstf[ss,(1+(jj-1)*3):(1+(jj-1)*3+2)] <- single[jj-1,1:3]
    for ( jj in 2:p ) simstf[ss,(3*p+(jj-2)*2+1):(3*p+(jj-2)*2+2)] <- joint[jj-1,]
  }   # end of simulation loop

  # Save simulation results in nsim rows
  simstf <- as.data.frame(simstf) # Use this the first time
 
  #  Summarize the results
  nS <- nrow(simstf)

  # Total wild
  TotHat <- simstf$Total
  biasTot <- mean(TotHat) - TrueWild
  pctbiasTot <- 100*biasTot/TrueWild
  varnTot<- var(TotHat)
  mseTot <- varnTot + biasTot^2
  biasTot <- round(biasTot,digits=3)
  pctbiasTot <- round(pctbiasTot,digits=3)
  rmseTot <- round(sqrt(mseTot),digits=3)
  seTot <- round(sqrt(varnTot),digits=3)
  varnTot <- round(varnTot,digits=3)
  mseTot <- round(mseTot,digits=3)
  cover <- sum( simstf$TL < TrueWild & simstf$TU > TrueWild )
  coverTot <- round(cover/nS,digits=3)
  EwidthTot <- round(mean(simstf$TU - simstf$TL),digits=3)

  # Group summary
  # Apr 26 '14 Add a summary row for number of simulation iterations used in joint coverage
  sumrys <- matrix(numeric(13*nGrps),ncol=nGrps)

  srn <- c("Truth        ","Bias         ","PercentBias   ","Variance    ","MeanSquareError  ","RootMSE    ","StandardError  ","CIcover    ","E(width)")
  srn <- c(srn, "PctHalfCI/Truth", "No. iter. in joint CIs", "E(sim_width)", "SimWidth/Width")
  colnames(sumrys) <- Originnames
  rownames(sumrys) <- srn
  # Apr 26 '14 Get saved joint conf intervals from simstf
  JointCIs <- as.matrix(simstf[,(3*p+1):(3*p+2*nGrps)])
  coverit <- matrix(NA,nrow=nS,ncol=nGrps)           # coverit tracks coverage row-by-row and group-by-group
  for (by in 1:nGrps){
    sumrys[1,by] <- TrueGroups[by]
    GroupsHat <- simstf[,4+(by-1)*3]
    sumrys[2,by] <- mean(GroupsHat) - TrueGroups[by] # Bias
    sumrys[3,by] <- 100*sumrys[2,by]/TrueGroups[by]  # Percent bias
    sumrys[4,by] <- var(GroupsHat)                   # Variance  
    sumrys[5,by] <- sumrys[4,by] + sumrys[2,by]^2    # Mean square error
    sumrys[6,by] <- sqrt(sumrys[5,by])               # Root mean square error
    sumrys[7,by] <- sqrt(sumrys[4,by])               # Standard error
    sumrys[8,by] <- sum( simstf[,5+(by-1)*3] < TrueGroups[by] & simstf[,6+(by-1)*3]> TrueGroups[by] ) 
    sumrys[8,by] <- sumrys[8,by]/nS                  # One-at-a-time coverage
    sumrys[9,by] <- mean( simstf[,6+(by-1)*3] - simstf[,5+(by-1)*3] ) # Expected width
    sumrys[10,by] <- 100*sumrys[9,by]/2/TrueGroups[by] # Half E(CI width) as percent of truth
    # Apr 26 '14 New calculation of joint confidence using all non-NA data in MBCI
    #  Get expected joint widths, ratios  and coverage for group 'by' based on joint CIs that are not (NA,NA)
    Ewidth <- NA  # Set expected joint width to NA initially
    nrvalid <- 0     # Set the number of simulation iterations used in calculating this by's joint coverage to 0
    OneJointCI <- JointCIs[,(2*(by-1)+1):(2*(by-1)+2)] # Look at joint CI in all rows of this group 
    validrows <- (1:nS)[!is.na(OneJointCI[,1])]        # Retain rows that are not (NA,NA)
    if ( length(validrows)>0 ) { # Check to see if there are any valid joint CIs
      validCIs <-OneJointCI[validrows,]
#      cat("\ndim validCIs ",dim(validCIs),"\n", file=console_messages_to, append=T)
      nrvalid <- nrow(validCIs)
      Ewidth <-mean(validCIs[,2] - validCIs[,1])
      coverit[validrows,by] <- validCIs[,1] < TrueGroups[by] & validCIs[,2] > TrueGroups[by]
# coverit is 1 or 0 if true value is in or out of joint CI for this grp; otherwise, it is NA
    }
    sumrys[11,by] <- nrvalid                            # Number of non(NA,NA) pairs in the current group
    sumrys[12,by] <- Ewidth                         # E(joint CI width)
    sumrys[13,by] <- sumrys[12,by]/sumrys[9,by]      # Ratio of E(joint CI width) to E(CI width)
  }   # End of nGrps loop
  
  jointcover <- rep(1,nS)
  for ( s in 1:nS ) {                          # Now go line by line through coverit
    veccover <- coverit[s,!is.na(coverit[s,])] # grap the non-NA individual coverages for all groups
    if ( length(veccover) > 0 ){               # in this row
      for ( j in 1:length(veccover) ) jointcover[s] <- jointcover[s] * veccover[j]
    }                                     # jointcover is 1 if ALL joint CIs in this row contain truth; o.w. it is 0
    else  jointcover[s] <- NA }           # It is NA if there are no valid entries in this row.
  coverALL <- sum(jointcover,na.rm=TRUE)
  coverALL <- coverALL/nS
  # End of new code Apr 26 '14

  sumrys = round(sumrys,digits=3)

  cat("\nNumber of simulations ",nS,"\n", file=console_messages_to, append=T)
  cat("\nPROPERTY                Total Wild", file=console_messages_to, append=T)
  cat("\nPopulation Truth       ",TrueWild, file=console_messages_to, append=T)
  cat("\nBias                   ",biasTot, file=console_messages_to, append=T)
  cat("\nPercent Bias           ",pctbiasTot, file=console_messages_to, append=T)
  cat("\nVariance               ",varnTot, file=console_messages_to, append=T)
  cat("\nMean Square Error      ",mseTot, file=console_messages_to, append=T)
  cat("\nRoot Mean Square Error ",rmseTot, file=console_messages_to, append=T)
  cat("\nStandard Error         ",seTot, file=console_messages_to, append=T)
  cat("\nCI coverage            ",coverTot, file=console_messages_to, append=T)
  cat("\nE(widthTotalWild)      ",EwidthTot, file=console_messages_to, append=T)
  cat("\nJoint CI coverage     ",coverALL,"\n", file=console_messages_to, append=T)
  # print(sumrys)  # eric commented this out because the function now returns this
  cat("\nEnd time: ",date(),"\n", file=console_messages_to, append=T)

  write.csv(sumrys, file = "summary.csv",row.names=TRUE)


  ######################## End of MAIN ##########################################
  return(sumrys)

}

