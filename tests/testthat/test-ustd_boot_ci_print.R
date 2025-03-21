skip_on_cran()

library(testthat)

test_that("ustd print", {

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
                                        estimator = "MLR")))
suppressWarnings(system.time(fit2 <- sem(model2,
                                         data = Data,
                                         estimator = "MLR",
                                         group = "gp")))
suppressWarnings(fit <- store_boot(fit,
                                   R = 100,
                                   iseed = 1234,
                                   do_bootstrapping = TRUE))
suppressWarnings(fit2 <- store_boot(fit2,
                                   R = 100,
                                   iseed = 1234,
                                   do_bootstrapping = TRUE))
ci_boot <- parameterEstimates_boot(fit, boot_pvalue_min_size = 100)
ci_boot2 <- parameterEstimates_boot(fit2, boot_pvalue_min_size = 100)

print(ci_boot, nd = 4)
print(ci_boot, output = "text")
print(ci_boot, output = "table")

print(ci_boot2, nd = 2)
print(ci_boot2, output = "text")
print(ci_boot2, output = "table")

expect_output(print(ci_boot, nd = 5),
              "bSE")
expect_output(print(ci_boot, output = "text"),
              "Defined Parameter")
expect_output(print(ci_boot, output = "table"),
              "lhs")

expect_output(print(ci_boot2, nd = 2),
              "0.00 ")
expect_output(print(ci_boot2, output = "text"),
              "Defined Parameter")
expect_output(print(ci_boot2, output = "table"),
              "block")

})