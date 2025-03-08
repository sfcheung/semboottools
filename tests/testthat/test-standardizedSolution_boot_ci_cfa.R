library(testthat)

library(lavaan)

test_that("CFA", {

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
                       se = "boot",
                       bootstrap = 100,
                       iseed = 1234,
                       warn = FALSE))

ci_boot <- standardizedSolution_boot(fit)

fit2 <- cfa(model,
            data = HolzingerSwineford1939,
            se = "none")
boot_ci_test <- suppressWarnings(bootstrapLavaan(fit2,
                                                 R = 100,
                                                 FUN = get_std,
                                                 iseed = 1234))

expect_equal(
  attr(ci_boot, "boot_est_std"),
  boot_ci_test,
  ignore_attr = TRUE
)

})
