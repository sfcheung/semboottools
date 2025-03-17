library(testthat)

library(lavaan)

test_that("parameterEstimates_boot: boot_org_ratio", {

set.seed(1234)
n <- 1000
x <- runif(n) - .5
m <- 0.20 * x + rexp(n)
y <- 0.17 * m + rexp(n)
dat <- data.frame(x, y, m)
mod <-
"
m ~ a*x
y ~ b*m + cp*x
ab := a*b
total := a*b + cp
"

fit <- sem(model = mod,
           data = dat,
           estimator = "MLR")
suppressWarnings(fit <- store_boot(fit,
                                   R = 100,
                                   iseed = 1234,
                                   do_bootstrapping = TRUE))
ci_boot2 <- parameterEstimates_boot(fit,
                                    level = .95,
                                    boot_pvalue_min_size = 50,
                                    boot_org_ratio = TRUE)
tmp <- as.data.frame(ci_boot2)
ratiolo <- abs(tmp$boot.ci.lower - tmp$est) /
           abs(tmp$ci.lower - tmp$est)
ratiohi <- abs(tmp$boot.ci.upper - tmp$est) /
           abs(tmp$ci.upper - tmp$est)
expect_identical(tmp$ratio.lower,
                 ratiolo)
expect_identical(tmp$ratio.upper,
                 ratiohi)
})
