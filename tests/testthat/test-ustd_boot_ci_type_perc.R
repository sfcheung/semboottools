skip("WIP")

library(testthat)

library(lavaan)

test_that("parameterEstimates_boot", {

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
                                        data = dat)))
fit <- store_boot(fit,
                  R = 100,
                  iseed = 1234,
                  do_bootstrapping = TRUE)
ci_boot <- parameterEstimates_boot(fit,
                                   level = .90,
                                   boot_pvalue_min_size = 50)
ci_boot

fit2 <- sem(model = mod,
            data = dat,
            se = "none")
suppressWarnings(boot_ci_test <- bootstrapLavaan(fit2,
                                                 R = 100,
                                                 FUN = get_std,
                                                 iseed = 1234))

# For lavaan 0.9-13 or later
boot_ci_test_error_idx <- attr(boot_ci_test, "error.idx")
if (!is.null(boot_ci_test_error_idx)) {
  if (length(boot_ci_test_error_idx) > 0) {
    boot_ci_test <- boot_ci_test[-boot_ci_test_error_idx, ]
  }
}

chk1 <- boot_ci_perc(t0 = ci_boot$est.std[1],
                     t = boot_ci_test[, 1],
                     level = .90)
chk2 <- unlist(ci_boot[1, c("boot.ci.lower", "boot.ci.upper")])
expect_equal(chk1,
             chk2,
             ignore_attr = TRUE)
chk1 <- boot_ci_perc(t0 = ci_boot$est.std[7],
                     t = boot_ci_test[, 7],
                     level = .90)
chk2 <- unlist(ci_boot[7, c("boot.ci.lower", "boot.ci.upper")])
expect_equal(chk1,
             chk2,
             ignore_attr = TRUE)

})
