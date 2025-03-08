#' @title Bootstrap CIs for Standardized
#' Solution
#'
#' @description Functions for forming
#' bootstrap confidence intervals
#' for the standardized solution.
#'
#' @details
#'
#' [standardizedSolution_boot_ci()]
#' receives a
#' [lavaan::lavaan-class] object fitted
#' with bootstrapping standard errors
#' requested and forms the confidence
#' intervals for the standardized
#' solution.
#'
#' It works by calling
#'  [lavaan::standardizedSolution()]
#' with the bootstrap estimates
#' of free parameters in each bootstrap sample
#' to compute the standardized estimates
#' in each sample.
#'
#' A more reliable way is to use
#' function like
#' [lavaan::bootstrapLavaan()].
#' Nevertheless, this simple function is
#' good enough for some simple scenarios,
#' and does not require repeating
#' the bootstrapping step.
#'
#' [store_boot_est_std()] computes the
#' standardized solution for each bootstrap
#' sample, stores them the
#' [lavaan::lavaan-class] object, and
#' returns it. These estimates can be used
#' by other functions, such as [plot_boot()],
#' to examine the
#' estimates, without the need
#' to repeat the computation.
#'
#' [get_boot_est_std()] retrieves
#' the bootstrap estimates of the
#' standardized solution stored by
#' [store_boot_est_std()].
#'
#' @return The output of
#' [lavaan::standardizedSolution()],
#' with bootstrap confidence intervals
#' appended to the right, with class
#' set to `std_solution_boot` (since
#' version 0.1.8.4). It has
#' a print method
#' ([print.std_solution_boot()]) that
#' can be used to print the standardized
#' solution in a format similar to
#' that of the printout of
#' the [summary()] of a [lavaan::lavaan-class] object.
#'
#' [store_boot_est_std()] returns
#' the fit object set to
#' `object`, with the bootstrap values
#' of standardized solution in the
#' bootstrap samples, as a matrix,
#' stored in the
#' slot `external` under the name
#' `shh_boot_est_std`.
#'
#' [get_boot_est_std()] returns a matrix
#' of the stored bootstrap estimates
#' of standardized solution. If none is
#' stored, `NULL` is returned.
#'
#' [store_boot_est_std()] is usually used
#' with diagnostic functions such
#' as [plot_boot()].
#'
#' @param object A 'lavaan'-class
#' object, fitted with 'se = "boot"'.
#'
#' @param level The level of confidence
#' of the confidence intervals. Default
#' is .95.
#'
#' @param type The type of standard
#' estimates. The same argument of
#' [lavaan::standardizedSolution()],
#' and support all values supported by
#' [lavaan::standardizedSolution()].
#' Default is `"std.all"`.
#'
#' @param boot_delta_ratio The ratio of
#' (a) the distance of the bootstrap
#' confidence limit from the point
#' estimate to (b) the distance of the
#' delta-method limit from the point
#' estimate. Default is `FALSE`.
#'
#' @param boot_ci_type The type of the
#' bootstrapping confidence intervals.
#' Support percentile confidence intervals
#' (`"perc"`, the default) and
#' bias-corrected confidence intervals
#' (`"bc"` or `"bca.simple"`).
#'
#' @param ... Other arguments to be
#' passed to
#' [lavaan::standardizedSolution()].
#'
#' @author Shu Fai Cheung
#' <https://orcid.org/0000-0002-9871-9448>.
#' Originally proposed in an issue at GitHub
#' <https://github.com/simsem/semTools/issues/101#issue-1021974657>,
#' inspired by a discussion at
#' the Google group for lavaan
#' <https://groups.google.com/g/lavaan/c/qQBXSz5cd0o/m/R8YT5HxNAgAJ>.
#' [boot::boot.ci()] is used to form the
#' percentile confidence intervals in
#' this version.
#'
#'
#' @seealso [lavaan::standardizedSolution()], [plot_boot()]
#'
#' @examples
#'
#' library(lavaan)
#' set.seed(5478374)
#' n <- 50
#' x <- runif(n) - .5
#' m <- .40 * x + rnorm(n, 0, sqrt(1 - .40))
#' y <- .30 * m + rnorm(n, 0, sqrt(1 - .30))
#' dat <- data.frame(x = x, y = y, m = m)
#' model <-
#' '
#' m ~ a*x
#' y ~ b*m
#' ab := a*b
#' '
#'
#' # Should set bootstrap to at least 2000 in real studies
#' fit <- sem(model, data = dat, fixed.x = FALSE,
#'            se = "boot",
#'            bootstrap = 100)
#' summary(fit)
#'
#' std <- standardizedSolution_boot_ci(fit)
#' std
#'
#' # Print in a friendly format with only standardized solution
#' print(std, output = "text")
#'
#' # Print in a friendly format with both unstandardized
#' # and standardized solution
#' print(std, output = "text", standardized_only = FALSE)
#'
#' # plot_boot() can be used to examine the bootstrap estimates
#' # of a parameter
#' plot_boot(std, param = "ab")
#'
# @name standardizedSolution_boot_ci
# NULL

# @rdname standardizedSolution_boot_ci
#' @noRd

standardizedSolution_boot_ci <- function(object,
                                         level = .95,
                                         type = "std.all",
                                         boot_delta_ratio = FALSE,
                                         boot_ci_type = c("perc", "bc", "bca.simple"),
                                         ...) {
  if (!inherits(object, "lavaan")) {
    stop("The object must be a lavaan-class object.")
  }
  # out_all <- boot_est_std(object = object,
  #                         type = type,
  #                         ...)
  out_all <- object@external$sbt_boot_std
  type0 <- object@external$sbt_boot_std_type
  if (is.null(out_all)) {
    # If customizing the call to store_boot()
    # is desired, should be done by calling
    # store_boot(), to avoid having too many
    # arguments for this function.
    object <- store_boot(object = object,
                         type = type,
                         do_bootstrapping = FALSE)
  } else {
    if (isFALSE(type0 != type)) {
      # Stored type and type in call do not match
      object <- store_boot(object = object,
                          type = type,
                          do_bootstrapping = FALSE)
    }
  }
  out <- lavaan::standardizedSolution(object,
                                      type = type,
                                      level = level,
                                      ...)
  est_org <- out$est.std

  # Adapted from boot
  boot_ci <- sapply(seq_along(est_org), function(x,
                                                  boot_type = boot_ci_type) {
                        if (isTRUE(all.equal(stats::var(out_all[, x], na.rm = TRUE), 0)) ||
                            all(is.na(out_all[, x]))) {
                            return(c(NA, NA))
                          }
                        boot_ci_internal(t0 = est_org[x],
                                          t = out_all[, x],
                                          level = level,
                                          boot_type = boot_type)
                      })
  boot_ci <- t(boot_ci)
  colnames(boot_ci) <- c("boot.ci.lower", "boot.ci.upper")
  boot_se <- apply(out_all, 2, stats::sd, na.rm = TRUE, simplify = TRUE)
  # Parameters with SE == NA are assumed to be fixed parameters
  boot_se[boot_se < .Machine$double.eps] <- NA
  out_final <- cbind(out, boot_ci, `boot.se` = boot_se)
  if (boot_delta_ratio) {
    tmp1 <- abs(out_final$boot.ci.lower - out_final$est.std) /
                              abs(out_final$ci.lower - out_final$est.std)
    tmp2 <- abs(out_final$boot.ci.upper - out_final$est.std) /
                              abs(out_final$ci.upper - out_final$est.std)
    tmp1[is.infinite(tmp1) | is.nan(tmp1)] <- NA
    tmp2[is.infinite(tmp2) | is.nan(tmp2)] <- NA
    out_final$ratio.lower <- tmp1
    out_final$ratio.upper <- tmp2
  }
  class(out_final) <- c("sbt_std_boot", class(out))
  fit_summary <- lavaan::summary(object)
  attr(out_final, "pe_attrib") <- attributes(fit_summary$pe)
  attr(out_final, "partable") <- lavaan::parameterTable(object)
  attr(out_final, "est") <- lavaan::parameterEstimates(object)
  attr(out_final, "level") <- level
  attr(out_final, "type") <- type
  attr(out_final, "call") <- match.call()
  out_final
}