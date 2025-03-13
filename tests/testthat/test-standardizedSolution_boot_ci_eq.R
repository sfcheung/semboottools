library(testthat)

library(lavaan)

test_that("Models with equality constraints", {

# Adapted from https://lavaan.ugent.be/tutorial/syntax2.html

data(HolzingerSwineford1939)
model <-
"
visual  =~ x1 + v2*x2 + v2*x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
"

system.time(fit <- cfa(model = model,
                       data = HolzingerSwineford1939,
                       se = "boot",
                       bootstrap = 50,
                       iseed = 1234))
ci_boot <- standardizedSolution_boot(fit,
                                     boot_pvalue_min_size = 50)

fit2 <- cfa(model = model,
            data = HolzingerSwineford1939,
            se = "none")
boot_ci_test <- bootstrapLavaan(fit2, R = 50,
                                FUN = get_std,
                                iseed = 1234)
bl_est <- bootstrapLavaan(fit2, R = 50,
                                FUN = coef,
                                iseed = 1234)
fit_bl_est <- fit
fit_bl_est@boot$coef <- bl_est
ci_boot_bl_est <- standardizedSolution_boot(fit_bl_est,
                                            boot_pvalue_min_size = 50)

# Cannot compare with get_std results because, even with same seed,
# bootstrapLavaan and se="boot" does not result in the same set
# of bootstrap samples.

expect_equal(
  attr(ci_boot_bl_est, "boot_est_std"),
  boot_ci_test,
  ignore_attr = TRUE
)

})
