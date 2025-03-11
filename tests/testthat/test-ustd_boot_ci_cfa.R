library(testthat)

library(lavaan)

test_that("parameterEstimates_boot: CFA", {

# Example from https://lavaan.ugent.be/tutorial/cfa.html
data(HolzingerSwineford1939)
model <-
"
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
"

system.time(fit <- cfa(model,
                       data = HolzingerSwineford1939,
                       estimator = "MLR"))
suppressWarnings(fit <- store_boot(fit,
                                   R = 50,
                                   iseed = 1234,
                                   do_bootstrapping = TRUE))
ci_boot <- parameterEstimates_boot(fit,
                                   boot_pvalue_min_size = 50)
fit2 <- cfa(model,
            data = HolzingerSwineford1939,
            se = "none")
boot_ci_test <- suppressWarnings(bootstrapLavaan(fit2,
                                                 R = 50,
                                                 FUN = coef,
                                                 iseed = 1234))

expect_equal(
  attr(ci_boot, "boot_est_ustd"),
  boot_ci_test,
  ignore_attr = TRUE
)

})
