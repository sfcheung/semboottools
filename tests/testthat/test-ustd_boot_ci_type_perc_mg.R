library(testthat)

library(lavaan)

test_that("parameterEstimates_boot", {

# Example from https://lavaan.ugent.be/tutorial/mediation.html

set.seed(1234)
n <- 1000
x <- runif(n) - .5
m <- 0.20 * x + rexp(n)
y <- 0.17 * m + rexp(n)
gp <- sample(c("Group1", "Group2"),
             n,
             replace = TRUE)
dat <- data.frame(X = x, Y = y, M = m, gp)
mod <- ' # direct effect
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
suppressWarnings(system.time(fit <- sem(model = mod,
                                        data = dat,
                                        estimator = "MLR",
                                        group = "gp")))
suppressWarnings(fit <- store_boot(fit,
                                   R = 100,
                                   iseed = 1234,
                                   do_bootstrapping = TRUE))
ci_boot <- parameterEstimates_boot(fit,
                                   level = .90,
                                   boot_pvalue_min_size = 50)
ci_boot

fit2 <- sem(model = mod,
            data = dat,
            se = "none",
            group = "gp")
suppressWarnings(boot_ci_test <- bootstrapLavaan(fit2,
                                                 R = 100,
                                                 FUN = coef,
                                                 iseed = 1234))

# For lavaan 0.9-13 or later
boot_ci_test_error_idx <- attr(boot_ci_test, "error.idx")
if (!is.null(boot_ci_test_error_idx)) {
  if (length(boot_ci_test_error_idx) > 0) {
    boot_ci_test <- boot_ci_test[-boot_ci_test_error_idx, ]
  }
}

chk1 <- boot_ci_perc(t0 = ci_boot$est[1],
                     t = boot_ci_test[, 1],
                     level = .90)
chk2 <- unlist(ci_boot[1, c("boot.ci.lower", "boot.ci.upper")])
expect_equal(chk1,
             chk2,
             ignore_attr = TRUE)
chk1 <- boot_ci_perc(t0 = ci_boot$est[20],
                     t = boot_ci_test[, "a2"] * boot_ci_test[, "b2"],
                     level = .90)
chk2 <- unlist(ci_boot[20, c("boot.ci.lower", "boot.ci.upper")])
expect_equal(chk1,
             chk2,
             ignore_attr = TRUE)

})
