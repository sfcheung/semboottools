library(testthat)

library(lavaan)

test_that("BC", {

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

suppressWarnings(system.time(fit <- sem(model = mod,
                                        data = dat,
                                        se = "boot",
                                        bootstrap = 200,
                                        iseed = 1234)))
est <- parameterEstimates(fit, boot.ci.type = "bca.simple")

suppressWarnings(ci_boot_bc <- standardizedSolution_boot(fit,
                                                         boot_ci_type = "bc"))
suppressWarnings(ci_boot_bca_simple <- standardizedSolution_boot(fit,
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

