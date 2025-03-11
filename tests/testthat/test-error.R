library(testthat)

library(lavaan)

test_that("Error if no boot estimates", {

# Example from https://lavaan.ugent.be/tutorial/mediation.html

set.seed(1234)
n <- 100
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

suppressWarnings(system.time(fit_noboot <- sem(model = mod,
                                               data = dat)))
expect_error(standardizedSolution_boot(fit_noboot))
expect_error(parameterEstimates_boot(fit_noboot))
})
