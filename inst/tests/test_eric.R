## Here are some basic arguments for the tests that we will tweak little 
## by little as we go along.
curdir <- getwd()
test_dir <- tempfile()
dir.create(test_dir)
test_data_dir <- file.path(system.file("data_files", package="lowergranite", mustWork=T), "data_for_testing")
basic_args <- list(
  stock_group_start_col = 14,
  DAT.DIR = test_data_dir,
  WORK.DIR = test_dir,
  nsim = 2, 
  B = 5, 
  DO_GSI_ON_PROP = F,
  BLFILE = file.path(test_data_dir, "sthd_super_informo_baseline.txt"),
  RUFILE = file.path(test_data_dir, "sthd_base_v3_rg.txt"),
  STOCK.DATA.XLSX = file.path(test_data_dir, "SH11SIMPOPstock.xlsx"),
  console_messages_to = "TestConsoleOutput.txt",
  reset_booty_seed=325
)



context("Test some of eric's small functions")


test_that("extract_multi_fix_sim_to_df catches no multi-fix headers in the input",{
  expect_that(extract_multi_fix_sim_to_df(c("bip", "boing"), lg_test_data$RU.list, lg_test_data$Originnames),
              throws_error("Apparently did"))
  expect_that(extract_multi_fix_sim_to_df(c("MULTI_FIX_MIX_MIXFISH_NAMES_HEADER:", "bip", "boing"), lg_test_data$RU.list, lg_test_data$Originnames),
              throws_error("Didn\'t find any"))
})


test_that("extract_multi_fix_sim_to_df catches produces correct result on known data",{
  expect_that(extract_multi_fix_sim_to_df(lg_test_data$gsi.big.out, lg_test_data$RU.list, lg_test_data$Originnames),
              is_identical_to(lg_test_data$multi_fix_df_compare))
})



context("Test that passing a data frame into run_boot_gsi_analysis() gives the same results as passing in the identical xlsx file")
test_that("specifying data frame and xlsx give identical results using a simple no-gsi example", {
  set.seed(15)
  xlsx_args <- basic_args
  xlsx_result <- do.call(run_boot_gsi_analysis, args=xlsx_args)
  
  set.seed(15)
  df_args <- xlsx_args
  df_args$W <- read.xlsx(xlsx_args$STOCK.DATA.XLSX, 1)  # pass in a data frame
  df_args$STOCK.DATA.XLSX <- "NoRealPath"  # wipe out the xlsx path information and put something bogus in there
  df_result <- do.call(run_boot_gsi_analysis, args=df_args)
  
  expect_that(df_result, is_identical_to(xlsx_result))
})


context("Test for consistency of known results on short analyses")


test_that("The no-gsi case gives expected results on short analysis of test data", {
  set.seed(15)
  cons_arg <- basic_args
  cons_arg$reset_booty_seed <- 0
  known_stock_result15 <- do.call(run_boot_gsi_analysis, args=cons_arg)
  
  expect_that(known_stock_result15, is_identical_to(lg_test_data$known_stock_result15))
})


test_that("The with-gsi case gives expected results on short analysis of test data", {
  set.seed(15)
  cons_arg <- basic_args
  cons_arg$reset_booty_seed <- 0
  cons_arg$BLFILE <- file.path(test_data_dir, "sthd_base_v3_187.txt")
  cons_arg$DO_GSI_ON_PROP <- T
  cons_arg$GSI_SEEDS <- c(1234, 5678) 
  gsi_stock_result15 <- do.call(run_boot_gsi_analysis, args=cons_arg)
  
  expect_that(gsi_stock_result15, is_identical_to(lg_test_data$gsi_stock_result15))
})


test_that("GSI with different seeds gives different results", {
  set.seed(15)
  cons_arg <- basic_args
  cons_arg$reset_booty_seed <- 0
  cons_arg$BLFILE <- file.path(test_data_dir, "sthd_base_v3_187.txt")
  cons_arg$DO_GSI_ON_PROP <- T
  cons_arg$GSI_SEEDS <- c(1234, 5678) 
  gsi_stock_result_a <- do.call(run_boot_gsi_analysis, args=cons_arg)
  
  cons_arg$GSI_SEEDS <- c(765, 987) 
  gsi_stock_result_b <- do.call(run_boot_gsi_analysis, args=cons_arg)
  
  expect_that(identical(gsi_stock_result_a, gsi_stock_result_b), is_false())
})







context("Test eric's gsi-sim replacements")

test_that("gsi_sim replacements with super-informo-data give same results as known stock", {
  # we just do a very short run with no gsi and the same run with gsi using super informative data
  # which should return exactly the same result.
  rep_args <- basic_args
  set.seed(5)
  known_stock_result <- do.call(run_boot_gsi_analysis, args=rep_args)
  
  rep_args$DO_GSI_ON_PROP <- T  # set it up to make and use gsi_assignments
  set.seed(5)
  super_informo_gsi_result <- do.call(run_boot_gsi_analysis, args=rep_args)
                                                    
  expect_that(known_stock_result, is_identical_to(super_informo_gsi_result))
  
})





# at the end change back to the working directory
setwd(curdir)
cat(c("\n\nTesting was done in directory:", test_dir, "\n"))
