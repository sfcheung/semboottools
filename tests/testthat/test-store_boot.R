skip("To be revised")

# skip_on_cran()
# skip_if(!interactive(),
#         message = "standardizedSolution_boot_ci not tested if not interactive")


# Input:
# - A lavaan object
# Output:
# - A lavaan object with bootstrap estimates stored.
store_boot <- function(object,
                       type = "std.all",
                       do_bootstrapping = FALSE,
                       R = 1000,
                       boot_type = "ordinary",
                       parallel = c("no", "multicore", "snow"),
                       ncpus = 1,
                       cl = NULL,
                       iseed = NULL,
                       bootstrapLavaan_args = list(),
                       ...) {
  if (!inherits(object, "lavaan")) {
    stop("The object must be a lavaan-class object.")
  }
  boot_est0 <- try(lavaan::lavTech(object, "boot"),
                   silent = TRUE)
  if (inherits(boot_est0, "try-error") && !do_bootstrapping) {
    stop("Bootstrapping estimates not found. Was se = 'boot' or 'bootstrap'?")
  }
  if (inherits(boot_est0, "try-error") && do_bootstrapping) {
    bootstrapLavaan_args1 <- utils::modifyList(bootstrapLavaan_args,
                                               list(object = object,
                                                    R = R,
                                                    type = boot_type,
                                                    parallel = parallel,
                                                    ncpus = ncpus,
                                                    cl = cl,
                                                    iseed = iseed,
                                                    keep.idx = TRUE))
    boot_ustd_out <- do.call(lavaan::bootstrapLavaan,
                             bootstrapLavaan_args1)
    object@external$sbt_boot_ustd <- boot_ustd_out
  }

  boot_std_out <- boot_est_std(object = object,
                               type = type,
                               ...)
  object@external$sbt_boot_std <- boot_std_out
  object@external$sbt_boot_std_type <- type
  boot_def_out <- boot_def(object = object)
  object@external$sbt_boot_def <- boot_def_out
  object
}

boot_lavaan <- function(object,
                        ...) {
  out <- lavaan::bootstrapLavaan(object = object,
                                 ...)
  out
}

library(testthat)
library(semboottools)

# Example from https://lavaan.ugent.be/tutorial/mediation.html

library(lavaan)
set.seed(1234)
n <- 1000
X <- runif(n) - .5
M <- 0.20*X + rnorm(n)
Y <- 0.17*M + rnorm(n)
Data <- data.frame(X = X, Y = Y, M = M)
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
# One bootstrap replication failed. Kept intentionally.
suppressWarnings(system.time(fit <- sem(model,
                       data = Data,
                       se = "boot",
                       bootstrap = 100,
                       iseed = 1234)))
fit_noboot <- sem(model,
                  data = Data)
fit <- store_boot(fit)
fit@external$sbt_boot_std
fit@external$sbt_boot_def

fit_noboot <- store_boot(fit_noboot, R = 100, iseed = 1234, do_bootstrapping = TRUE)
fit_noboot@external$sbt_boot_std
fit_noboot@external$sbt_boot_def

expect_equal(lavInspect(fit, "boot"),
             fit_noboot@external$sbt_boot_ustd,
             ignore_attr = TRUE)
expect_equal(fit@external$sbt_boot_std,
             fit_noboot@external$sbt_boot_std)
expect_equal(fit@external$sbt_boot_def,
             fit_noboot@external$sbt_boot_def)

get_std <- function(object) {
    lavaan::standardizedSolution(object)$est.std
  }
# fit2 <- update(fit, se = "none")
fit2 <- sem(model,
            data = Data,
            se = "none",
            bootstrap = 100)
set.seed(1234)
suppressWarnings(boot_ci_test <- bootstrapLavaan(fit2, R = 100,
                                FUN = get_std))

# For lavaan 0.9-13
boot_ci_test_error_idx <- attr(boot_ci_test, "error.idx")
if (!is.null(boot_ci_test_error_idx)) {
    if (length(boot_ci_test_error_idx) > 0) {
        boot_ci_test <- boot_ci_test[-boot_ci_test_error_idx, ]
      }
  }

test_that("Compare boot estimates directly", {
    expect_equal(
        attr(ci_boot, "boot_est_std"),
        boot_ci_test,
        ignore_attr = TRUE
      )
  })

# Test store_boot_est_std()

tmp <- store_boot_est_std(fit)
tmp_boot_est_std <- get_boot_est_std(tmp)
test_that("store_boot_est_std", {
    expect_equal(tmp_boot_est_std,
                 attr(ci_boot, "boot_est_std"),
                 ignore_attr = TRUE)
    expect_equal(tmp@external$shh_boot_est_std_type,
                 "std.all")
  })

# Test Bias-corrected CI

# Example from https://lavaan.ugent.be/tutorial/mediation.html

model <- ' # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
           # Self-computed standardized coefficient
             X ~~ v_x*X
             M ~~ ev_m*M
             Y ~~ ev_y*Y
             sd_x := sqrt(v_x)
             sd_m := sqrt((a^2)*(v_x) + ev_m)
             sd_y := sqrt((b^2)*(sd_m^2) + ev_y)
             a_std := a * sd_x / sd_m
             b_std := b * sd_m / sd_y
             ab_std := ab * sd_x / sd_y
         '
set.seed(1234)
# One bootstrap replication failed. Kept intentionally.
suppressWarnings(system.time(fit <- sem(model,
                       data = Data,
                       se = "boot",
                       bootstrap = 100)))
est <- parameterEstimates(fit, boot.ci.type = "bca.simple")
# print(est, nd = 5)
# print(standardizedSolution(fit), nd = 5)

ci_boot_bc <- standardizedSolution_boot_ci(fit, save_boot_est_std = TRUE, boot_ci_type = "bc")
ci_boot_bca_simple <- standardizedSolution_boot_ci(fit, save_boot_est_std = TRUE, boot_ci_type = "bca.simple")

test_that("Compare boot estimates directly", {
    expect_equal(
        ci_boot_bc[1, "boot.ci.lower"],
        est[10, "ci.lower"],
        ignore_attr = TRUE
      )
    expect_equal(
        ci_boot_bc[1, "boot.ci.upper"],
        est[10, "ci.upper"],
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

