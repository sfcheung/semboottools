library(testthat)

library(lavaan)
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

test_that("standardizedSolution_boot_ci", {

# Example from https://lavaan.ugent.be/tutorial/mediation.html

# One bootstrap replication failed. Kept intentionally.
suppressWarnings(system.time(fit <- sem(model = mod,
                                        data = dat,
                                        se = "boot",
                                        bootstrap = 100,
                                        iseed = 1234)))

ci_boot <- standardizedSolution_boot_ci(fit,
                                        level = .90)

get_std <- function(object) {
    tmp <- lavaan::standardizedSolution(object,
                                        se = FALSE)
    pnames <- lavaan::lav_partable_labels(tmp)
    out <- tmp$est.std
    names(out) <- pnames
    out
  }
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

test_that("Compare boot estimates directly", {

# Example from https://lavaan.ugent.be/tutorial/mediation.html

mod <-
"
m ~ a*x
y ~ b*m + cp*x
ab := a*b
total := a*b + cp

# Self-computed standardized coefficient
x ~~ v_x*x
m ~~ ev_m*m
y ~~ ev_y*y
sd_x := sqrt(v_x)
sd_m := sqrt((a^2)*(v_x) + ev_m)
sd_y := sqrt((b^2)*(sd_m^2) + ev_y)
a_std := a * sd_x / sd_m
b_std := b * sd_m / sd_y
ab_std := ab * sd_x / sd_y
"

# One bootstrap replication failed. Kept intentionally.
suppressWarnings(system.time(fit <- sem(model = mod,
                                        data = dat,
                                        se = "boot",
                                        bootstrap = 200,
                                        iseed = 1234)))
est <- parameterEstimates(fit, boot.ci.type = "bca.simple")
# print(est, nd = 5)
# print(standardizedSolution(fit), nd = 5)

suppressWarnings(ci_boot_bc <- standardizedSolution_boot_ci(fit,
                                                            boot_ci_type = "bc"))
suppressWarnings(ci_boot_bca_simple <- standardizedSolution_boot_ci(fit,
                                                                    boot_ci_type = "bca.simple"))

    expect_equal(
        ci_boot_bc[1, "boot.ci.lower"],
        est[12, "ci.lower"],
        ignore_attr = TRUE
      )
    expect_equal(
        ci_boot_bc[1, "boot.ci.upper"],
        est[12, "ci.upper"],
        ignore_attr = TRUE
      )
    expect_equal(
        ci_boot_bc[12, "boot.ci.lower"],
        est[12, "ci.lower"],
        ignore_attr = TRUE
      )
    expect_equal(
        ci_boot_bc[12, "boot.ci.upper"],
        est[12, "ci.upper"],
        ignore_attr = TRUE
      )
    expect_equal(
        ci_boot_bca_simple[12, "boot.ci.lower"],
        est[12, "ci.lower"],
        ignore_attr = TRUE
      )
    expect_equal(
        ci_boot_bca_simple[12, "boot.ci.upper"],
        est[12, "ci.upper"],
        ignore_attr = TRUE
      )

})

