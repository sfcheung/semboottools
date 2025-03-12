skip_on_cran()

library(testthat)

test_that("std print", {

# Example from https://lavaan.ugent.be/tutorial/mediation.html

library(lavaan)
set.seed(1234)
n <- 1000
X <- runif(n) - .5
M <- 0.20*X + rnorm(n)
Y <- 0.17*M + rnorm(n)
gp <- sample(c("Group1", "Group2"),
             n,
             replace = TRUE)
Data <- data.frame(X = X, Y = Y, M = M, gp = gp)
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
model2 <- ' # direct effect
             Y ~ c(c1, c2)*X
           # mediator
             M ~ c(a1, a2)*X
             Y ~ c(b1, b2)*M
           # indirect effect (a*b)
             a1b1 := a1*b2
             a2b2 := a2*b2
           # total effect
             total1 := c1 + (a1*b1)
             total2 := c2 + (a2*b2)
         '
suppressWarnings(system.time(fit <- sem(model,
                                        data = Data,
                                        se = "boot",
                                        bootstrap = 200,
                                        iseed = 1234)))
suppressWarnings(system.time(fit2 <- sem(model2,
                                         data = Data,
                                         se = "boot",
                                         bootstrap = 200,
                                         group = "gp",
                                         iseed = 1234)))

suppressWarnings(system.time(fit_noboot <- sem(model,
                                               data = Data)))
suppressWarnings(fit_noboot <- store_boot(fit_noboot,
                                          R = 200,
                                          iseed = 1234))

ci_boot <- standardizedSolution_boot(fit, boot_pvalue_min_size = 199)
ci_boot2 <- standardizedSolution_boot(fit2, boot_pvalue_min_size = 200)
ci_boot_noboot <- standardizedSolution_boot(fit_noboot, boot_pvalue_min_size = 199)

print(ci_boot)
print(ci_boot, standardized_only = FALSE, nd = 5)
print(ci_boot, boot_ci_only = TRUE)
print(ci_boot_noboot, nd = 5)
print(ci_boot_noboot, boot_ci_only = TRUE)
print(ci_boot_noboot, standardized_only = FALSE)

print(ci_boot, output = "text")
print(ci_boot, output = "text", boot_ci_only = TRUE)
print(ci_boot_noboot, output = "text")
print(ci_boot_noboot, output = "text", boot_ci_only = TRUE)
print(ci_boot_noboot, output = "text", standardized_only = FALSE)

expect_output(print(ci_boot, nd = 5),
              " bSE ")
expect_output(print(ci_boot, output = "text"),
              "boot.se")
expect_output(print(ci_boot, output = "text", standardized_only = FALSE),
              "Estimate")

expect_output(print(ci_boot2, nd = 5),
              "0.00000")
expect_output(print(ci_boot2, output = "text"),
              "Defined Parameter")
expect_output(print(ci_boot2, output = "text", standardized_only = FALSE),
              "Estimate")

})
