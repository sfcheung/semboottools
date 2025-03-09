#' @title Bootstrap CIs for Standardized
#' Solution
#'
#' @description Functions for forming
#' bootstrap confidence intervals
#' for the standardized solution.
#'
#' @details
#'
#' [standardizedSolution_boot()]
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
#' Alternative, call [store_boot()] to
#' computes and store bootstrap estimates
#' of the standardized solution.
#' This function will then retrieve them,
#' even if `se` was not set to
#' `"boot"` or `"bootstrap"` when fitting
#' the model.
#'
#' @return The output of
#' [lavaan::standardizedSolution()],
#' with bootstrap confidence intervals
#' appended to the right, with class
#' set to `sbt_std_boot`. It has
#' a print method
#' ([print.sbt_std_boot()]) that
#' can be used to print the standardized
#' solution in a format similar to
#' that of the printout of
#' the [summary()] of a [lavaan::lavaan-class] object.
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
#' @param save_boot_est_std Whether the
#' bootstrap estimates of the
#' standardized solution are saved. If
#' saved, they will be stored in the
#' attribute `boot_est_std`. Default is
#' `TRUE`.
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
#' @seealso [lavaan::standardizedSolution()], [store_boot()]
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
#' std <- standardizedSolution_boot(fit)
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
#' # scatter_boot() can be used to examine the bootstrap estimates
#' # of two or more parameters
#' scatter_boot(std, params = c("ab", "a", "b"))
#'
#' @name standardizedSolution_boot
NULL

#' @rdname standardizedSolution_boot
#' @export

standardizedSolution_boot <- function(object,
                                         level = .95,
                                         type = "std.all",
                                         boot_delta_ratio = FALSE,
                                         boot_ci_type = c("perc", "bc", "bca.simple"),
                                         save_boot_est_std = TRUE,
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
    out_all <- object@external$sbt_boot_std
  } else {
    if (isFALSE(type0 != type)) {
      # Stored type and type in call do not match
      object <- store_boot(object = object,
                          type = type,
                          do_bootstrapping = FALSE)
      out_all <- object@external$sbt_boot_std
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

  if (save_boot_est_std) {
    colnames(out_all) <- std_names(object, ...)
    attr(out_final, "boot_est_std") <- out_all
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