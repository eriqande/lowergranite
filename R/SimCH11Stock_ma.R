#' Run the chinook script
#' 
#' This is kirk's script for chinook that eric wrapped up in a function
#' @param WORK.DIR  The working directory to do this in.  Default = current working directory
#' @param STOCK.DATA.XLSX the file that has the stock data in it used to drive the simulations
#' @param DO_GSI if set to TRUE then it uses GSI assignments. Otherwise it uses the true origins
run_chinook_script <- function(
  WORK.DIR = getwd(),  
  STOCK.DATA.XLSX = "S:/Eagle Fish Genetics Lab/Ackerman/LGR/Kirk/ms/stock_sims/inputs/CH11SIMPOP.xlsx",
  DO_GSI = TRUE,
  BLFILE = "C:/GenSoftware/gsi_sim_pc/analysis/chnk_base_v3_180.txt",
  RUFILE = "C:/GenSoftware/gsi_sim_pc/analysis/chnk_base_v3_repunits.txt", 
  GSISIM = "C:/GenSoftware/gsi_sim_pc/gsi_sim.exe",  # path to the gsi_sim executable
  GSISIM.FOLDER = "C:/GenSoftware/gsi_sim_pc/analysis/"          # path to the gsi_sim analysis folder
) {

  ########## DONE WITH SETTING VALUES  ###########################################################
  ########## ALLOW FOR SETTING DIFFERENT VALUES ON THE COMMAND LINE #######
  comm.args <- commandArgs(TRUE) # this grabs the string of command args (if there is one)
  comm.args  # just print them here for debugging
  for(i in comm.args) {eval(parse(text=i))}  # this parses the command line arguments and sets values

 
  #############  PREAMBLE ADDITIONS TO GATHER GSI BASELINE FILES AND REPORTING UNITS, ETC ############ 
  if(DO_GSI==TRUE)  {
    # source some files of necessary code:
    source(file.GSI.SIM.R.STUFF)
    source(file.LOWER.GRANITE.GSI.FUNCS)
  
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
  
  }  # close if(DO_GSI==TRUE)
  ####################################################################################################

  # Get data. Start by setting the default working directory.
  setwd(WORK.DIR) 

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
      tstar <- numeric(nw)
      nHstar <- numeric(nw)
      for ( h in as.vector(windata[,1])) {
        tstar[h] <- rbinom(1,nT[h],Trp[h]) # These are number of wild fish for week h
        nHstar[h] <- round( tstar[h]*W$PWhandled[h] ) # This is the number of bootstrap wild fish handled for week h
        tstar[h] <- tstar[h]/nT[h] # Convert binomial count into a proportion
      }
      nHcollaps <- mApply(nHstar,collaps,sum)
      pstar <- matrix(numeric(nCollaps*nGrps),ncol=nGrps)
      for ( h in 1:nCollaps )
        pstar[h,] <- rmultinom(1,nHcollaps[h],Pr[h,])/nHcollaps[h]   # These are the proportions by sex/origin year for this stratum

  # Calculate wild by week for the bootstrap data
      WiByWe <- WBW(tstar)  # This is wild by week for all statistical weeks
      WeWi <- mApply(WiByWe,collaps,sum)  # This is the wild by week collapsed for origin calculation
      GW <- WeWi %*% pstar  # Bootstrap numbers of wild fish by sex/origin
      theta.b[b,] <- c(sum(WiByWe),GW)
       } # End of bootstrap loop

  # Find confidence intervals for each statistic
    for  (j in 1:p) {
      CIj <- quantile(theta.b[,j],c(alph/2,1-alph/2))
      CI[j,] <- c(ests[j],CIj)
    }
    CI <- round(CI)

  # Find simultaneous confidence intervals for each statistic
  # But first save the one-at-a-time intervals
    WOATCI <- CI

      ctr <- 0
      chk <- 0
      while(chk < (1-alph)){
        ctr <- ctr + 1
        chk <- 0
        for(i in 1:B){
          chkr <- 1
          for(j in 2:p){
            if( (theta.b[i,j] < CI[j,2]) | (theta.b[i,j] > CI[j,3]) ) chkr <- chkr * 0
          } # for over parameters
          chk <- chk + chkr/B
        } # for over bootstrap samples
        for (j in 2:p){
          CI[j,2] <- CI[j,2]*(1-0.001)
          CI[j,3] <- CI[j,3]*(1+0.001)
        }
      }  # while chk

  #    CI <- round(CI)

    return( rbind(WOATCI,CI) )
  }
  #################### End of bootstrap CI function ##############################


  #################### MAIN  #####################################################

  # Define basic run parameters
  Run <- "2011 Chinook Origin -- weighted, weighted/collapsed"

  alph <- 0.10
  B    <- 500
  nsim <- 500

  cat("\nStart time: ",date(),"\n")
  cat("\nThis is a run of ", Run, "\n")
  cat("\nParametric bootstrap: B = ",B,"Alpha= ",alph,"\n")

  # Read in weekly counts and define the true POPULATION

  W <- read.xlsx(STOCK.DATA.XLSX,1)
  #W <- read.xlsx("CH11SIMPOP.xlsx",1)
  W
  W$SimPop <- round(W$SimPop)
  TrueWild <- sum(W$PopWild)
  #TrueSex <- cbind(W$Pfemale, 1 - W$Pfemale)
  #SexNames <- c("Female","Male")
  #TrueAge <- cbind(W$BY04,W$BY05,W$BY06,W$BY07,W$BY08)
  #BYnames <- c("BY04","BY05","BY06","BY07","BY08")
  #nA <- length(BYnames)
  TrueOrigin <- cbind(W$CHMBLN,W$FALL,W$HELLSC,W$MFSALM,W$SFSALM,W$TUCANO,W$UPSALM)
  Originnames <- c("CHMBLN","FALL","HELLSC","MFSALM","SFSALM","TUCANO","UPSALM")
  nO <- length(Originnames)
  nGrps <- nO
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

  # Simulation loop
  for (ss in 1:nsim){

  #  generate random trap and female proportions-- by week
    nT <- numeric(nw)
    nW <- numeric(nw)
    nH <- numeric(nw)
    for (ii in 1:nw){
      nT[ii] <- rbinom(1,W$SimPop[ii],W$PTrp[ii])  # Sample 4 to 12 percent of the run.  This is the number of fish trapped.
      tStrat <- rep(ii,nT[ii]) # Define the trap strata column
      nW[ii] <- rbinom(1,nT[ii],W$Pwild[ii])   # This is the number of trapped fish that are wild.
      grp <- c(rep("WILD",nW[ii]),rep("H",nT[ii]-nW[ii]))  # Define the column of WILD and H(atchery) data
      nH[ii] <- round(nW[ii]*W$PWhandled[ii]) # This is the number of trapped fish sexed or aged or sampled for genetic analysis
      Stratum <- rep(ii,nH[ii])  # Define the stratum column for the proportion data
      Groups <- rmultinom(1,nH[ii],Probab[ii,])
      Groop <- "   " # This is just a programming device to avoid an if statement
      for (jj in 1:nGrps) {
         if(Groups[jj] >0) Groop <- c(Groop,rep(Originnames[jj],Groups[jj]))}
         Groop <- Groop[2:length(Groop)] # Now delete the placeholder entry
      if (ii == 1) {
         Trap <- data.frame(tStrat,grp)
         Prop <- data.frame(Stratum,Groop) }
      else { Trap <- rbind(Trap,data.frame(tStrat,grp))
           Prop <-rbind(Prop,data.frame(Stratum,Groop))  }
    }
    # Convert the random trap and proportion data sets to wild and origin proportions
    ntrap <- nrow(Trap)
    nprop <- nrow(Prop)
    Wtable <- table(Trap)
    Wprop <- Wtable/apply(Wtable,1,sum)  # This is the proportion hatchery and wild
  #  Ptable <- table(Prop)
  # Kluge the Ptable because not all groups occur in every statistical week
    Ptable <- matrix(numeric(nw*nGrps),ncol=nGrps)
    for ( h in 1:nw){
      justh <- Prop[Prop[,1]==h,2]
      for ( hh in 1:nGrps ) Ptable[h,hh] <- sum(justh==Originnames[hh])
    }
  
  ########################  MORE GSI STUFF.  If DO.GSI==T THIS RECONSTRUCTS Ptable ###################
    if(DO_GSI==TRUE) {
      Ptable.tmp <- Ptable  # make a new variable so as to not mess with the original Ptable
      colnames(Ptable.tmp) <- Originnames
    
      GSItable <- matrix(NA, nrow=nrow(Ptable.tmp), ncol=ncol(Ptable.tmp))  # a place to store the results
      colnames(GSItable) <- Originnames
    
      # compile the gsi_sim commands that tell how many fish from each population should be simulated
      SeedToRestore <- .Random.seed  # this is a bit of a hack to make sure that the same true origins of fish are simulated whether you use GSI or the true identity of the fish, given the same initial seed
      for(pti in 1:nrow(Ptable.tmp)) {
        fixed.pi.comms <- gsi.sim.Write.Fixed.Pi(RU.list, BL.pops, Ptable.tmp[pti,])
      
        # here is the call we want to make with GSI_SIM
        setwd(GSISIM.FOLDER)
        GSI.CALL <- paste(GSISIM, "-b", BLFILE, "-x 1 1", fixed.pi.comms)
      
        # then make the call and get the results
        GSItable[pti, ] <- gsi.simCallAndCondense(GSI.CALL,  BL.pops,  RU.list,  Originnames)
      
        Ptable <- GSItable # here we hand the GSI results back as if they were the true counts to be operated on....
      
      }
      .Random.seed <- SeedToRestore  # restores the seed to where it was (rolling back the multinom calls in gsi.sim.Write.Fixed.Pi)
    
    } # closes if(DO_GSI==TRUE)
  ####################################################################################################
  
  # Although the proportion data are generated for each statistical week, the data are collapsed to insure adequate sample sizes
    collaps <- c(1,1,2,3,4,5,6,7,8,9,10,11,11,11)
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
    simstf[ss,1:3] <- CIs[1,1:3]
    for ( jj in 2:p ) simstf[ss,(1+(jj-1)*3):(1+(jj-1)*3+2)] <- CIs[jj,1:3]
    for ( jj in 2:p ) simstf[ss,(3*p+(jj-2)*2+1):(3*p+(jj-2)*2+2)] <- CIs[p+jj,2:3]
  
  }   # end of simulation loop

  # Save simulation results to excel file
  simstf <- as.data.frame(simstf) # Use this the first time
  allstf <- as.data.frame(simstf) # Use this the first time 
  res <- write.xlsx(allstf,"CH11stock.xlsx",col.names=TRUE,row.names=FALSE,append=FALSE) #  Use this the first time

  #  Summarize the results
  #allstf <- read.xlsx("CH11OriginDec13.xlsx",1) # Get stored results
  nS <- nrow(allstf)

  # Total wild
  TotHat <- allstf$Total
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
  cover <- sum( allstf$TL < TrueWild & allstf$TU > TrueWild )
  coverTot <- round(cover/nS,digits=3)
  EwidthTot <- round(mean(allstf$TU - allstf$TL),digits=3)

  # Origin summary
  sumrys <- matrix(numeric(12*nGrps),ncol=nGrps)

  srn <- c("Truth        ","Bias         ","Percent bias   ","Variance    ","Mean Square Error  ","Root MSE    ","Standard Error  ","CI cover    ","E(width)")
  srn <- c(srn, "PctHalfCI/Truth", "E(sim_width)", "SimWidth/Width")
  colnames(sumrys) <- Originnames
  rownames(sumrys) <- srn
  coverit <- matrix(numeric(nS*nGrps),ncol=nGrps)
  for (by in 1:nGrps){
    sumrys[1,by] <- TrueGroups[by]
    GroupsHat <- allstf[,4+(by-1)*3]
    sumrys[2,by] <- mean(GroupsHat) - TrueGroups[by]
    sumrys[3,by] <- 100*sumrys[2,by]/TrueGroups[by]
    sumrys[4,by] <- var(GroupsHat)
    sumrys[5,by] <- sumrys[4,by] + sumrys[2,by]^2
    sumrys[6,by] <- sqrt(sumrys[5,by])
    sumrys[7,by] <- sqrt(sumrys[4,by])
    sumrys[8,by] <- sum( allstf[,5+(by-1)*3] < TrueGroups[by] & allstf[,6+(by-1)*3]> TrueGroups[by] )
    sumrys[8,by] <- sumrys[8,by]/nS
    sumrys[9,by] <- mean( allstf[,6+(by-1)*3] - allstf[,5+(by-1)*3] )
    sumrys[10,by] <- 100*sumrys[9,by]/2/TrueGroups[by]
    sumrys[11,by] <- mean( allstf[,26+(by-1)*2] - allstf[,25+(by-1)*2]) # This line only works for Origin
    sumrys[12,by] <- sumrys[11,by]/sumrys[9,by]    
    coverit[,by] <-  allstf[,3*(nGrps+1)+1+(by-1)*2] < TrueGroups[by] & allstf[,3*(nGrps+1)+2+(by-1)*2] > TrueGroups[by]
  }
  jointcover <- rep(1,nS)
  for (by in 1:nGrps)
    jointcover <- jointcover & coverit[,by]
  coverALL <- sum(jointcover)
  coverALL <- coverALL/nS

  sumrys = round(sumrys,digits=3)

  cat("\nNumber of simulations ",nS,"\n")
  cat("\nPROPERTY                Total Wild")
  cat("\nPopulation Truth       ",TrueWild)
  cat("\nBias                   ",biasTot)
  cat("\nPercent Bias           ",pctbiasTot)
  cat("\nVariance               ",varnTot)
  cat("\nMean Square Error      ",mseTot)
  cat("\nRoot Mean Square Error ",rmseTot)
  cat("\nStandard Error         ",seTot)
  cat("\nCI coverage            ",coverTot)
  cat("\nJoint CI coverage     ",coverALL,"\n")
  print(sumrys)
  cat("\nEnd time: ",date(),"\n")

  write.csv(sumrys, file = "summary.csv",row.names=TRUE)

  ######################## End of MAIN ##########################################
}

