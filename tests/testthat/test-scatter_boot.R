skip("To be run in an interactive mode")

library(testthat)
library(semboottools)

# Example from https://lavaan.ugent.be/tutorial/mediation.html

library(lavaan)
set.seed(1234)
n <- 100
X <- runif(n) - .5
M <- 0.20 * X + rnorm(n)
Y <- 0.17 * M + rnorm(n)
GP <- sample(c("GpA", "GpB"), replace = TRUE)
Data <- data.frame(X = 10 * X, Y = Y, M = 20 * M, GP = GP)
model <- ' # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
modelgp <- ' # direct effect
             Y ~ c(c1, c2)*X
           # mediator
             M ~ c(a1, a2)*X
             Y ~ c(b1, b2)*M
           # indirect effect (a*b)
             ab1 := a1*b1
             ab2 := a2*b2
           # total effect
             total1 := c1 + (a1*b1)
             total1 := c2 + (a2*b2)
         '
# One bootstrap replication failed. Kept intentionally.
suppressWarnings(system.time(fit <- sem(model,
                                        data = Data,
                                        se = "boot",
                                        bootstrap = 100,
                                        iseed = 1234)))
suppressWarnings(system.time(fitgp <- sem(modelgp,
                                          data = Data,
                                          se = "boot",
                                          group = "GP",
                                          bootstrap = 100,
                                          iseed = 1234)))

fit_boot <- store_boot(fit)
fitgp_boot <- store_boot(fitgp)

scatter_boot(fit_boot, c("ab", "a", "b"), standardized = TRUE)
scatter_boot(fit_boot, c("ab", "a", "b"), standardized = FALSE)
scatter_boot(fit_boot, c("c", "a", "b"), standardized = TRUE)
scatter_boot(fit_boot, c("c", "a", "b"), standardized = FALSE)
scatter_boot(fit_boot, c("c", "a", "b"), standardized = FALSE)

scatter_boot(fitgp_boot, c("ab1", "a1", "b1"), standardized = TRUE)
scatter_boot(fitgp_boot, c("ab1", "a1", "b1"), standardized = FALSE)
scatter_boot(fitgp_boot, c("ab1", "a1", "b1", "ab2", "a2", "b2"), standardized = TRUE)
scatter_boot(fitgp_boot, c("ab1", "a1", "b1", "ab2", "a2", "b2"), standardized = FALSE)

test_that("Expect errors", {
    expect_error(scatter_boot(fit_boot, "X~~X", standardized = TRUE))
    expect_error(scatter_boot(fit_boot, "X~~X", standardized = FALSE))
    expect_error(scatter_boot(fitgp_boot, "X~~X", standardized = TRUE))
    expect_error(scatter_boot(fitgp_boot, "X~~X", standardized = FALSE))
  })

# Support standardizedSolution_boot()

std <- standardizedSolution_boot(fit)
stdgp <- standardizedSolution_boot(fitgp)

# Examine interactively

scatter_boot(std, c("ab", "a", "b"))
scatter_boot(fit_boot, c("ab", "a", "b"), standardized = TRUE)

scatter_boot(std, c("c", "a", "b"))
scatter_boot(fit_boot, c("c", "a", "b"), standardized = TRUE)

scatter_boot(stdgp, c("M~~M", "M~~M.g2"))
scatter_boot(fitgp_boot, c("M~~M", "M~~M.g2"), standardized = TRUE)

