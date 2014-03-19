


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



context("Test eric's gsi-sim replacements")

test_that("gsi_sim replacements with super-informo-data give same results as known stock", {
  set.seed(5)
  known_stock_result <- run_boot_gsi_analysis(nsim=2, B=5, DO_GSI_ON_PROP = F)
  set.seed(5)
  super_informo_gsi_result <- run_boot_gsi_analysis(nsim=2, B=5, DO_GSI_ON_PROP = T,
                                                   BLFILE = file.path(DAT.DIR, "data_for_testing", "sthd_test_super_informo_baseline.txt")
                                                   )
  expect_that(known_stock_result[1:7,], is_identical_to(super_informo_gsi_result[1:7,]))
  
})