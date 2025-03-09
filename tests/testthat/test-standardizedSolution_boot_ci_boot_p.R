library(testthat)

library(lavaan)

test_that("standardizedSolution_boot: boot.p", {

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
                                        se = "boot",
                                        bootstrap = 100,
                                        iseed = 1234)))

suppressWarnings(ci_boot <- standardizedSolution_boot(fit,
                                     level = .90))
expect_true(is.null(ci_boot$boot.p))

ci_boot <- standardizedSolution_boot(fit,
                                     level = .90,
                                     boot_pvalue_min_size = 99)

fit <- store_boot(fit)
ci_boot_p_chk <- est2p(fit@external$sbt_boot_std[, "ab"],
                       min_size = 99)
expect_equal(ci_boot$boot.p[7],
             ci_boot_p_chk)
})
