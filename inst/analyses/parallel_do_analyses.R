


# establish a place to put the output.  This is not machine independent!
out_path <- "/tmp/lowergranite_outputs"
dir.create(out_path)

shs <- steelhead_run_settings
chs <- chinook_run_settings

# put these together with a new column that will tell us
# which function to use.
shs$Type <- "steelhead"
chs$Type <- "chinook"

ps <- rbind(shs, chs)

# this is just here for testing
ps$nsim <- 3
ps$B <- 2




# since I can't do the run, I will at least write an R script that I can launch with
# GNU parallel...
stage_a_run <- function(w) {
  if(w$Type=="steelhead") {
    f <- "do_steelhead"
  } else if(w$Type=="chinook") {
    f <- "do_chinook"
  }
  newdir <- file.path(out_path, paste(w$File, w$start_col, w$nsim, w$B, w$DoGSI, sep="_"))
  dir.create(newdir)
  setwd(newdir)
  
  pream <- "\nlibrary(lowergranite)\n"
  call <-  sprintf("\nresults <- %s(\"%s\", %s, %s, %s, %s)\n", f, w$File, w$start_col, w$nsim, w$B, w$DoGSI)
  wrap_up <- "\nsave(results, file=\"results.rda\")\n "
  
  cat(pream, call, wrap_up, file = "rscript.R")
}


# so, do this.
lapply(1:nrow(ps), function(x) stage_a_run(ps[x, ]))


# then on the unix command line in out_path do:
#       for i in *; do echo "cd $i; R CMD BATCH rscript.R;" ; done > for_para.txt 
# then you can launch it all with gnu parallel like this:
#       parallel -P 2 < for_para.txt
#



##### SINCE MCLAPPLY DOESN'T SEEM TO REALLY WORK I CAN'T DO THIS AND INSTEAD JUST WRITE
#### A BUNCH OF SEPARATE R SCRIPTS ABOVE AND THEN I WILL RUN THEM WITH 
#### GNU PARALLEL.
# do_a_run <- function(f, w) {
#   newdir <- file.path(out_path, paste(w$File, w$start_col, w$nsim, w$B, w$DoGSI, sep="_"))
#   dir.create(newdir)
#   setwd(newdir)
#   results <- f(w$File, w$start_col, w$nsim, w$B, w$DoGSI)
#   save(results, file="results.rda")
# }
# 
# 
# # here do all the runs. Hey.  This does fine with lapply, but mclapply fails miserably and I 
# # do not know why.  Sucks!
# mclapply(1:nrow(ps), mc.set.seed = TRUE, mc.preschedule = FALSE, mc.cleanup = FALSE,  function(x) {
#   w <- ps[x, ]
#   if(w$Type=="steelhead") {
#     f <- do_steelhead
#   } else if(w$Type=="chinook") {
#     f <- do_chinook
#   }
#   else {
#     stop(paste("unknown type of analysis", w$Type))
#   }
#   do_a_run(f, w)
# }
# )
