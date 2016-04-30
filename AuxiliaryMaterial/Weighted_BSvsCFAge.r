# Simulate data for Steelhead 2011 stock June 29, 2015 weighted, POOLED, parametric

# Get data. Start by setting the default working directory.
setwd("C:\\Kirk\\NES\\IDFG\\TimSim\\SimulationsMarchApril2016")
library(xlsx)    # Note.  You have to install the xlsx package.
library(MCPAN)   # This package contains the code for the Mandel-Betensky confidence intervals
######################################################################

ClosedForm <- function(Tr,Pr) {
#  Tr <- Trp
#  Pr <- Pprop
  
# Find asymptotic CI for total wild

  wildweek <- C_h*Tr[,2]
  TotalWild <- sum(wildweek)
  VarTerm <- (C_h - nT)*Trp[,2]*(1-Trp[,2])/C_h/(nT-1)
  VarSeries <- C_h*C_h*VarTerm
  VarW <- sum(VarSeries)
  seW <- sqrt(VarW)
  HalfWidth <- 1.645*seW
  Wild <- c(TotalWild, TotalWild-HalfWidth, TotalWild+HalfWidth)
  
# Now find asymptotic CIs for genetic groups
#  ngg <- ncol(Pr)
  GenGrpWild <- wildweek %*% Pr

  for( gg in 1:nGrps ) {
    VarPiHat <- (wildweek - nH)*Pr[,gg]*(1-Pr[,gg])/wildweek/(nH-1)
    Goodman <-  VarPiHat*Tr[,2]*Tr[,2] + VarTerm*Pr[,gg]*Pr[,gg] - VarPiHat*VarTerm
    VarSeries <- C_h*C_h*Goodman
    VarNg <- sum(VarSeries)
    seNg <- sqrt(VarNg)
    HalfWidth <- 1.645*seNg
    Wild <- c(Wild,GenGrpWild[gg],GenGrpWild[gg]-HalfWidth,GenGrpWild[gg]+HalfWidth)
  }
  
  return(Wild)
}
  
  


## Bootstrap function
AllBoot <- function(Tr,Pr) {

# Set up storage for bootstrap results
  theta.b <- matrix(numeric((1+nGrps)*(B+1)),ncol=1+nGrps)
  thetaNames <- c("Wild",WgrpNames)
  colnames(theta.b) <- thetaNames

# bootstrap loop
  for (b in 1:(B+1)) {
#b <- 1
    if ( b == 1 ) {
      tstar <- Tr
      pstar <- Pr
    }
    else {
# First generate a tstar (a bootstrap version of Trp)
      tstar <- matrix(numeric(2*length(Strata)),ncol=2)
      for ( h in Strata ) {
        Eq0 <- which(Tr[h,] == 0)
        if( length(Eq0) == 2 )  tstar[h,] == c(0,0)  # if all probabilities are 0 then this stratum has no fish
        if( length(Eq0) == 1 )  { tstar[h,] <- Tr[h,] # if there is a 0, then probabilities are (1,0) or (0,1)
        } else {  # if there are 2 nonzero entries then we can generate bootstrap fish
        tstar[h,] <- t(rmultinom(1,nT[h],Tr[h,])) # These are number of hatchery and wild fish for week h
        tstar[h,] <- tstar[h,]/nT[h] # Convert multinomial count into a proportion
        }  # end of else
      }  # end of for
# Then generate a pstar (a bootstrap version of Pprop)
        pstar <- matrix(numeric(length(Strata)*nGrps),ncol=nGrps)
      for ( h in Strata ) {
        Eq0 <- which(Pr[h,] == 0)
        if( length(Eq0) == nGrps ) pstar[h,] <- rep(0,nGrps) # if all probabilities are 0 then this stratum has no fish of type RTYPE
        else if( length(Eq0) == nGrps-1 ) { pstar[h,] <- Pr[h,] # There is only one group with fish. No randomness here.
        } else pstar[h,] <- rmultinom(1,nH[h],Pr[h,])/nH[h]   # These are the bootstrap proportions by Wgrp for this stratum
      }  # end of for
    } # end of else from b==1 if statement

# Whether b = 1 or b > 1

# Calculate RTYPE by week for the bootstrap data
    RBW <- c(C_h)*tstar  # This is numbers by week for all rearing types
#    NhatRear <- apply(RearByWe,2,sum)
    Cats <- RBW[,2] %*% pstar  # Bootstrap numbers of wild fish by group

    theta.b[b,]  <- c(sum(RBW[,2]),Cats)


  } # End of bootstrap loop

   theta.grps <- theta.b[,-1]
#  pp <- ncol(theta.grps)

# Find one-at-a-time confidence intervals for each statistic
  CI <- matrix(numeric(p*2),ncol=2)
  for  (j in 1:p) {
    CI[j,] <- quantile(theta.b[,j],c(alph/2,1-alph/2))
  }
  CI <- round(CI)
  TotalWildCI <- CI[1,]
  OneCI <- CI[2:p,]

# Find simultaneous rectangular confidence intervals via Mandel/Betensky 2008
  MBCI <- SCSrank(theta.grps,conf.level=1-alph)
  MBCI <- MBCI$conf.int
  MBCI <- round(MBCI)
 
  return( list(theta.b[1,],TotalWildCI,OneCI,MBCI) ) 
}
#################### End of bootstrap CI function ##############################


#################### MAIN  #####################################################

# Define basic run parameters
#Run <- "2011 Chinook Sex -- weighted, weighted/collapsed"
Run <- "2011 Chinook Age -- weighted, weighted/collapsed"
#Run <- "2011 Steelhead Stock -- weighted, weighted/collapsed parametric"

alph <- 0.1
B <- 500
nsim <- 500
collaps <- c(1,1,1,1,1,2,2,2,3,4,5,6,7,8,9,10,10,10,10,11,11,11,11,11,11,11,11)
Strata <- unique(collaps)
cat("\nStart time: ",date(),"\n")
cat("\nThis is a run of ", Run, "\n")
cat("\nParametric bootstrap: B = ",B,"Alpha = ", alph, " with Closed Form CIs\n")

# Read in weekly counts and define the true POPULATION

W <- read.xlsx("SH11SIMPOP.xlsx",1)
W
W$SimPop <- round(W$SimPop)
TrueWild <- sum(W$PopWild)
TrueWgrp <- cbind(W$BY04,W$BY05,W$BY06,W$BY07,W$BY08)
WgrpNames <- c("BY04","BY05","BY06","BY07","BY08")
nGrps <- length(WgrpNames)
nw <- length(unique(W$Stratum))
Probab <- matrix(numeric(nw*nGrps),ncol=nGrps)
for (i in 1:nw)
Probab[i,] <- TrueWgrp[i,]
Truth <- Probab*W$PopWild
TrueGroups <- mApply(t(Truth),1:nGrps,sum)
names(TrueGroups) <- WgrpNames

windata <- data.frame(W$Stratum,W$SimPop)
windata
C_h <- mApply(windata[,2],collaps,sum)

p <- nGrps + 1 
cols <- p*3 + 2*(p-1)
simstf <- matrix(numeric(nsim*cols),ncol= cols)
cn <- c("Total","TL","TU")
for (ii in 1:nGrps)  cn <- c(cn,WgrpNames[ii],paste(WgrpNames[ii],"L",sep=""),paste(WgrpNames[ii],"U",sep=""))
for (ii in 1:nGrps)  cn <- c(cn,paste(WgrpNames[ii],"Ls",sep=""),paste(WgrpNames[ii],"Us",sep=""))
colnames(simstf) <- cn

# Simulation loop
for (ss in 1:nsim){
#ss <- 1

#  generate random trap and group proportions-- by week
  nTrp <- numeric(nw)
  nW <- numeric(nw)
  nHand <- numeric(nw)
  for (ii in 1:nw){
    nTrp[ii] <- rbinom(1,W$SimPop[ii],W$PTrp[ii])  # Sample 4 to 12 percent of the run.  This is the number of fish trapped.
    tStrat <- rep(ii,nTrp[ii]) # Define the trap strata column
    nW[ii] <- rbinom(1,nTrp[ii],W$Pwild[ii])   # This is the number of trapped fish that are wild.
    grp <- c(rep("WILD",nW[ii]),rep("H",nTrp[ii]-nW[ii]))  # Define the column of WILD and H(atchery) data
    nHand[ii] <- round(nW[ii]*W$PWhandled[ii]) # This is the number of trapped fish sexed or aged or sampled for genetic analysis
    Stratum <- rep(ii,nHand[ii])  # Define the stratum column for the proportion data
    Groups <- rmultinom(1,nHand[ii],Probab[ii,])
    Groop <- "   " # This is just a programming device to avoid an if statement.
    for (jj in 1:nGrps) {
       if(Groups[jj] >0) Groop <- c(Groop,rep(WgrpNames[jj],Groups[jj]))}
       Groop <- Groop[2:length(Groop)] # Now delete the placeholder entry
    if (ii == 1) {
       Trap <- data.frame(tStrat,grp)
       Prop <- data.frame(Stratum,Groop) }
    else { Trap <- rbind(Trap,data.frame(tStrat,grp))
         Prop <-rbind(Prop,data.frame(Stratum,Groop))  }
  }

# Convert the random trap and proportion data sets to wild and gen stock proportions
# for PARAMETRIC bootstrap  

#  ntrap <- nrow(Trap)
  Wtable <- table(Trap) # Frequency of H and W by original week
#  Wprop <- Wtable/apply(Wtable,1,sum)  # This is the proportion hatchery and wild by original week
  collapsed <- mApply(Wtable[,1],collaps,sum)
  collapsed <- cbind(collapsed,mApply(Wtable[,2],collaps,sum))
  nT <- apply(collapsed,1,sum)
  Trp <- collapsed/nT  # This is the proportion hatchery and wild by statistical week
  
  nprop <- nrow(Prop)
  Ptable <- table(factor(Prop$Stratum, levels=1:length(unique(W$Stratum))), factor(Prop$Groop, levels=WgrpNames))
# Although the genetic frequency data are generated for each statistical week, the data are collapsed to insure adequate sample sizes
  collapsed <- mApply(Ptable[,1],collaps,sum)
  for ( jj in 2:nGrps ) collapsed <- cbind(collapsed,mApply(Ptable[,jj],collaps,sum))
  nH <- apply(collapsed,1,sum) # Number of wild fish handled
#  nHan <- mApply(nHand,collaps,sum)  # Also number of wild fish handled by Strata
  Pprop <- collapsed/nH # Proportions of fish by genetic group by statistical week 

# Call ClosedForm routine to find aymptotic CIs

  CF <- ClosedForm(Trp,Pprop)
  if( ss == 1 ) CFstf <- CF  else CFstf <- rbind(CFstf,CF)

# Call bootstrap routine to find one-at-a-time and simultaneous confidence intervals

  CIs <- AllBoot( Trp,Pprop )
  

# Now put simulation results into rows of an array for later analysis
  TotalWild <- c(CIs[[1]][1],CIs[[2]])
  simstf[ss,1:3] <- TotalWild 
  single <- cbind(CIs[[1]][-1],CIs[[3]])
  joint <- CIs[[4]][1:p-1,]
  for ( jj in 2:p ) simstf[ss,(1+(jj-1)*3):(1+(jj-1)*3+2)] <- single[jj-1,1:3]
  for ( jj in 2:p ) simstf[ss,(3*p+(jj-2)*2+1):(3*p+(jj-2)*2+2)] <- joint[jj-1,]
  
}   # end of simulation loop

# Save simulation results to excel file
simstf <- as.data.frame(simstf) 
res <- write.xlsx(simstf,"SH11.xlsx",col.names=TRUE,row.names=FALSE,append=FALSE) #  Use this the first time

#  Summarize the results
#simstf <- read.xlsx("CH11SexDec13.xlsx",1) # Get stored results
  nS <- nrow(simstf)
  
  # Bootstrap Total wild
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
  P1Tot <- 100*EwidthTot/2/TrueWild
  
  # Closed Form Total wild
  CFTotHat <- CFstf[,1]
  CFbiasTot <- mean(CFTotHat) - TrueWild
  CFpctbiasTot <- 100*CFbiasTot/TrueWild
  CFvarnTot<- var(CFTotHat)
  CFmseTot <- CFvarnTot + CFbiasTot^2
  CFbiasTot <- round(CFbiasTot,digits=3)
  CFpctbiasTot <- round(CFpctbiasTot,digits=3)
  CFrmseTot <- round(sqrt(CFmseTot),digits=3)
  CFseTot <- round(sqrt(CFvarnTot),digits=3)
  CFvarnTot <- round(CFvarnTot,digits=3)
  CFmseTot <- round(CFmseTot,digits=3)
  CFcover <- sum( CFstf[,2] < TrueWild & CFstf[,3] > TrueWild )
  CFcoverTot <- round(CFcover/nS,digits=3)
  CFEwidthTot <- round(mean(CFstf[,3] - CFstf[,2]),digits=3)
  CFP1Tot <- 100*CFEwidthTot/2/TrueWild
  
  cat("\nNumber of simulations ",nS,"\n")
  cat("\nPROPERTY                Boot    ClosedForm")
  cat("\nPopulation Truth       ",TrueWild)
  cat("\nBias                   ",biasTot,CFbiasTot)
  cat("\nPercent Bias           ",pctbiasTot,CFpctbiasTot)
  cat("\nVariance               ",varnTot,CFvarnTot)
  cat("\nMean Square Error      ",mseTot,CFmseTot)
  cat("\nRoot Mean Square Error ",rmseTot,CFrmseTot)
  cat("\nStandard Error         ",seTot,CFseTot)
  cat("\nCI coverage            ",coverTot,CFcoverTot)
  cat("\nExpected width wild    ",EwidthTot,CFEwidthTot,"\n")
  cat("\nP1 total wild    ",P1Tot,CFP1Tot,"\n")
  
  # Bootstrap Group summary
  
  sumrys <- matrix(numeric(12*nGrps),ncol=nGrps)
  
  srn <- c("Truth        ","Bias         ","Percent bias   ","Variance    ","Mean Square Error  ","Root MSE    ","Standard Error  ","CI cover    ","E(width)")
  srn <- c(srn, "PctHalfCI/Truth", "E(sim_width)", "SimWidth/Width")
  colnames(sumrys) <- WgrpNames
  rownames(sumrys) <- srn
  coverit <- matrix(numeric(nS*nGrps),ncol=nGrps)
  for (by in 1:nGrps){
    sumrys[1,by] <- TrueGroups[by]
    GroupsHat <- simstf[,4+(by-1)*3]
    sumrys[2,by] <- mean(GroupsHat) - TrueGroups[by]
    sumrys[3,by] <- 100*sumrys[2,by]/TrueGroups[by]
    sumrys[4,by] <- var(GroupsHat)
    sumrys[5,by] <- sumrys[4,by] + sumrys[2,by]^2
    sumrys[6,by] <- sqrt(sumrys[5,by])
    sumrys[7,by] <- sqrt(sumrys[4,by])
    sumrys[8,by] <- sum( simstf[,5+(by-1)*3] < TrueGroups[by] & simstf[,6+(by-1)*3]> TrueGroups[by] )
    sumrys[8,by] <- sumrys[8,by]/nS
    sumrys[9,by] <- mean( simstf[,6+(by-1)*3] - simstf[,5+(by-1)*3] )
    sumrys[10,by] <- 100*sumrys[9,by]/2/TrueGroups[by]
    sumrys[11,by] <- mean( simstf[,3*(nGrps+1)+2+(by-1)*2] - simstf[,3*(nGrps+1)+1+(by-1)*2]) 
    sumrys[12,by] <- sumrys[11,by]/sumrys[9,by]    
    coverit[,by] <-  simstf[,3*(nGrps+1)+1+(by-1)*2] < TrueGroups[by] & simstf[,3*(nGrps+1)+2+(by-1)*2] > TrueGroups[by]
  }
  jointcover <- rep(1,nS)
  for (by in 1:nGrps)
    jointcover <- jointcover & coverit[,by]
  coverALL <- sum(jointcover)
  coverALL <- coverALL/nS
  cat("\nJoint CI coverage      ",coverALL,"\n")
  
  sumrys = round(sumrys,digits=3)
  
  cat("\nSimulation summary for bootstrap CIs\n")
  print(sumrys)
  
  saveBScoverage <- sumrys[8,]
  saveBSEW <- sumrys[9,]
  
  # Closed Form Group summary
  cat("\nC L O S E D  F O R M results\n")
  
  
  sumrys <- matrix(numeric(10*nGrps),ncol=nGrps)
  
  srn <- c("Truth        ","Bias         ","Percent bias   ","Variance    ","Mean Square Error  ","Root MSE    ","Standard Error  ","CI cover    ","E(width)")
  srn <- c(srn, "PctHalfCI/Truth")
  colnames(sumrys) <- WgrpNames
  rownames(sumrys) <- srn
  for (by in 1:nGrps){
    sumrys[1,by] <- TrueGroups[by]
    GroupsHat <- CFstf[,4+(by-1)*3]
    sumrys[2,by] <- mean(GroupsHat) - TrueGroups[by]
    sumrys[3,by] <- 100*sumrys[2,by]/TrueGroups[by]
    sumrys[4,by] <- var(GroupsHat)
    sumrys[5,by] <- sumrys[4,by] + sumrys[2,by]^2
    sumrys[6,by] <- sqrt(sumrys[5,by])
    sumrys[7,by] <- sqrt(sumrys[4,by])
    sumrys[8,by] <- sum( CFstf[,5+(by-1)*3] < TrueGroups[by] & CFstf[,6+(by-1)*3]> TrueGroups[by] )
    sumrys[8,by] <- sumrys[8,by]/nS
    sumrys[9,by] <- mean( CFstf[,6+(by-1)*3] - CFstf[,5+(by-1)*3] )
    sumrys[10,by] <- 100*sumrys[9,by]/2/TrueGroups[by]
  }
   sumrys = round(sumrys,digits=3)
  cat("\nSimulation summary for closed form CIs\n")
  print(sumrys)
  cat("\n")
  
  saveCFcoverage <- sumrys[8,]
  saveCFEW <- sumrys[9,]
  
  compear <-  rbind(round(TrueGroups),saveBScoverage,saveCFcoverage,saveBScoverage-saveCFcoverage,saveBSEW,saveCFEW,saveBSEW-saveCFEW)
  rownames(compear) <- c("Truth","BScover","CFcover","Diff","BSwidth","CFwidth","Diff")
  print(compear)
  
  cat("\nEnd time: ",date(),"\n")

######################## End of MAIN ##########################################