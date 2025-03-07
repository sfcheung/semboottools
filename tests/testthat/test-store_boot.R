library(testthat)
library(lavaan)

test_that("store_boot", {

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

# One bootstrap replication failed. Kept intentionally.
suppressWarnings(system.time(fit <- sem(model = mod,
                                        data = dat,
                                        se = "boot",
                                        bootstrap = 100,
                                        iseed = 4567)))
fit_noboot <- sem(model = mod,
                  data = dat)
fit <- store_boot(fit)

fit_noboot <- store_boot(fit_noboot,
                         R = 100,
                         iseed = 4567,
                         do_bootstrapping = TRUE)

expect_equal(lavInspect(fit, "boot"),
             fit_noboot@external$sbt_boot_ustd,
             ignore_attr = TRUE)
expect_equal(fit@external$sbt_boot_std,
             fit_noboot@external$sbt_boot_std)
expect_equal(fit@external$sbt_boot_def,
             fit_noboot@external$sbt_boot_def)

})
