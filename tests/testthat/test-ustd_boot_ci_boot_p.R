library(testthat)

library(lavaan)

test_that("parameterEstimates_boot: boot.p", {

# Example from https://lavaan.ugent.be/tutorial/mediation.html

set.seed(1234)
n <- 1000
x <- runif(n) - .5
m <- 0.20 * x + rnorm(n)
y <- 0.17 * m + rnorm(n)
dat <- data.frame(x, y, m)
mod <-
"
m ~ a*x
y ~ b*m + cp*x
ab := a*b
total := a*b + cp
"

suppressWarnings(system.time(fit <- sem(model = mod,
                                        data = dat,
                                        estimator = "MLR")))
suppressWarnings(fit <- store_boot(fit,
                                   R = 50,
                                   iseed = 1234,
                                   do_bootstrapping = TRUE))
suppressWarnings(ci_boot <- parameterEstimates_boot(fit,
                                     level = .90))
expect_true(is.null(ci_boot$boot.p))

ci_boot <- parameterEstimates_boot(fit,
                                   level = .90,
                                   boot_pvalue_min_size = 50)
ci_boot_p_chk <- est2p(fit@external$sbt_boot_def[, "ab"],
                       min_size = 50)
expect_equal(ci_boot$boot.p[7],
             ci_boot_p_chk)
})
