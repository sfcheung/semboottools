#' @title Bootstrap CIs for Parameter
#' Estimates
#'
#' @description Functions for forming
#' bootstrap confidence intervals
#' for the parameter estimates.
#'
#' @details
#'
#' [parameterEstimates_boot()]
#' receives a
#' [lavaan::lavaan-class] object and
#' form bootstrap confidence intervals
#' for the parameter estimates.
#'
#' The function [store_boot()] should
#' be called first to
#' compute and store bootstrap estimates.
#' This function will then retrieve them.
#'
#' ## Bootstrap Confidence Intervals
#'
#' It supports percentile and
#' bias-corrected bootstrap confidence
#' intervals.
#'
#' ## Bootstrap Standard Errors
#'
#' The standard errors are the
#' standard deviation of the bootstrap
#' estimates.
#'
#' ## Bootstrap Asymmetric *p*-Values
#'
#' If percentile bootstrap confidence
#' interval is requested, asymmetric
#' bootstrap *p*-values are also
#' computed, using the method presented
#' in Asparouhov and Muthén (2021).
#'
#' @return The output of
#' [lavaan::parameterEstimates()],
#' with bootstrap confidence intervals
#' appended to the right, with class
#' set to `sbt_ustd_boot`. It has
#' a print method
#' ([print.sbt_ustd_boot()]) that
#' can be used to print the parameter
#' estimates in a format similar to
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
#' @param standardized The type of standardized
#' estimates. The same argument of
#' [lavaan::parameterEstimates()],
#' and support all values supported by
#' [lavaan::parameterEstimates()].
#'
#' @param boot_org_ratio The ratio of
#' (a) the distance of the bootstrap
#' confidence limit from the point
#' estimate to (b) the distance of the
#' original confidence limit in
#' `object` from the point
#' estimate. Default is `FALSE`.
#'
#' @param boot_ci_type The type of the
#' bootstrapping confidence intervals.
#' Support percentile confidence intervals
#' (`"perc"`, the default) and
#' bias-corrected confidence intervals
#' (`"bc"` or `"bca.simple"`).
#'
#' @param save_boot_est Whether the
#' bootstrap estimates of the
#' parameter estimates are saved. If
#' saved, the bootstrap estimates
#' of the free parameters will be stored
#' in the attribute `boot_est_ustd`,
#' while the bootstrap estimates of
#' user-defined parameters, if any,
#' will be stored in the attribute
#' `boot_def`. Default is
#' `TRUE`.
#'
#' @param boot_pvalue Whether asymmetric
#' bootstrap *p*-values are computed.
#' Default is `TRUE`.
#'
#' @param boot_pvalue_min_size Integer.
#' The asymmetric bootstrap *p*-values
#' will be computed only if the number
#' of valid bootstrap estimates is at
#' least this value. Otherwise, `NA`
#' will be returned. If the number of
#' valid bootstrap samples is less than
#' this value, then `boot_pvalue` will
#' be set to `FALSE`.
#'
#' @param ... Other arguments to be
#' passed to
#' [lavaan::parameterEstimates()].
#'
#' @author Shu Fai Cheung
#' <https://orcid.org/0000-0002-9871-9448>.
#'
#' @references
#' Asparouhov, A., & Muthén, B. (2021). Bootstrap p-value computation.
#' Retrieved from https://www.statmodel.com/download/FAQ-Bootstrap%20-%20Pvalue.pdf
#'
#' @seealso [lavaan::parameterEstimates()], [store_boot()]
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
#' fit <- sem(model, data = dat, fixed.x = FALSE)
#' summary(fit)
#' fit <- store_boot(fit,
#'                   do_bootstrapping = TRUE,
#'                   R = 100,
#'                   iseed = 1234)
#' est <- parameterEstimates_boot(fit)
#' est
#'
#' @name parameterEstimates_boot
#' NULL

#' @rdname parameterEstimates_boot
#' @export

parameterEstimates_boot <- function(object,
                                    level = .95,
                                    boot_org_ratio = FALSE,
                                    boot_ci_type = c("perc", "bc", "bca.simple"),
                                    save_boot_est = TRUE,
                                    boot_pvalue = TRUE,
                                    boot_pvalue_min_size = 1000,
                                    standardized = FALSE,
                                    ...) {
  type <- standardized
  if (is.logical(type)) type <- "std.all"
  boot_ci_type <- match.arg(boot_ci_type)
  if (!inherits(object, "lavaan")) {
    stop("The object must be a lavaan-class object.")
  }
  # out_all <- boot_est_std(object = object,
  #                         type = type,
  #                         ...)
  boot_est0 <- try(lavaan::lavTech(object, "boot"),
                   silent = TRUE)
  if (!inherits(boot_est0, "try-error")) {
    stop("This function should not be used if bootstrap SEs and CIs are already requested.")
  }
  out_all <- object@external$sbt_boot_ustd
  if (is.null(out_all)) {
    # If customizing the call to store_boot()
    # is desired, should be done by calling
    # store_boot(), to avoid having too many
    # arguments for this function.
    stop("Bootstrap estimates not stored. ",
         "Please call store_boot() first, with ",
         "do_bootstrapping set to 'TRUE'.")
    # object <- store_boot(object = object,
    #                      type = type,
    #                      do_bootstrapping = TRUE)
    # out_all <- object@external$sbt_boot_ustd
  }
  out <- lavaan::parameterEstimates(object,
                                    level = level,
                                    standardized = standardized,
                                    ...)
  ptable <- lavaan::parameterTable(object)
  out$id_est <- seq_len(nrow(out))
  if (!is.null(out$group)) {
    out0 <- merge(out,
                  ptable[, c("lhs", "op", "rhs", "group", "free")],
                  by = c("lhs", "op", "rhs", "group"))
  } else {
    out0 <- merge(out,
                  ptable[, c("lhs", "op", "rhs", "group", "free")],
                  by = c("lhs", "op", "rhs"),
                  all.x = TRUE,
                  all.y = FALSE,
                  sort = FALSE)
  }

  if (boot_ci_type != "perc") {
    # Asymmetric p-values computed only if percentile CIs are requested
    boot_pvalue <- FALSE
  }
  tmp <- nrow(out_all)
  if ((tmp < boot_pvalue_min_size) && boot_pvalue) {
    # Asymmetric p-values not computed because the number of bootstrap samples
    # is less than boot_pvalue_min_size
    tmp1 <- paste0("The number of bootstrap samples (%1$d) ",
                   "is less than 'boot_pvalue_min_size' (%2$d). ",
                   "Bootstrap p-values are not computed.")
    tmp2 <- sprintf(tmp1, tmp, boot_pvalue_min_size)
    warning(tmp2)
    boot_pvalue <- FALSE
  }

  # Adapted from boot
  boot_ci <- matrix(NA,
                    nrow = nrow(out),
                    ncol = 2)
  colnames(boot_ci) <- c("boot.ci.lower", "boot.ci.upper")
  boot_se <- rep(NA, nrow(out))
  boot_p <- rep(NA, nrow(out))
  # Form boot CI for free parameters
  for (i in seq_len(nrow(boot_ci))) {
    i_id <- out0[i, "free"]
    # Is x a fixed parameter?
    if (i_id < 1) next
    out_all_i <- out_all[, i_id]
    # Is x fixed in bootstrapping?
    tmp <- isTRUE(all.equal(stats::var(out_all_i, na.rm = TRUE), 0)) ||
           all(is.na(out_all_i))
    if (isTRUE(tmp)) next
    boot_ci_i <- boot_ci_internal(t0 = out0[i, "est"],
                                  t = out_all_i,
                                  level = level,
                                  boot_type = boot_ci_type)
    boot_ci[i, ] <- boot_ci_i
    boot_se[i] <- stats::sd(out_all_i,
                            na.rm = TRUE)
    if (boot_pvalue) {
      boot_p[i] <- est2p(x = out_all_i,
                         h0 = 0,
                         min_size = boot_pvalue_min_size,
                         warn = FALSE)
    }
  }
  # Form boot CI for user-defined parameters
  if (isTRUE(":=" %in% out0$op)) {
    out_all_def <- object@external$sbt_boot_def
    i0 <- which(out0$op %in% ":=")
    for (i in i0) {
      i_label <- out0$label[i]
      out_all_i <- out_all_def[, i_label]
      # Is x fixed in bootstrapping?
      tmp <- isTRUE(all.equal(stats::var(out_all_i, na.rm = TRUE), 0)) ||
            all(is.na(out_all_i))
      if (isTRUE(tmp)) next
      boot_ci_i <- boot_ci_internal(t0 = out0[i, "est"],
                                    t = out_all_i,
                                    level = level,
                                    boot_type = boot_ci_type)
      boot_ci[i, ] <- boot_ci_i
      boot_se[i] <- stats::sd(out_all_i,
                              na.rm = TRUE)
      if (boot_pvalue) {
        boot_p[i] <- est2p(x = out_all_i,
                          h0 = 0,
                          min_size = boot_pvalue_min_size,
                          warn = FALSE)
      }
    }
  } else {
    out_all_def <- NULL
  }
  out$id_est <- NULL
  if (boot_pvalue) {
    out_final <- cbind(out,
                       `boot.se` = boot_se,
                       `boot.p` = boot_p,
                       boot_ci)
  } else {
    out_final <- cbind(out,
                       `boot.se` = boot_se,
                       boot_ci)
  }
  if (boot_org_ratio) {
    tmp1 <- abs(out_final$boot.ci.lower - out_final$est.std) /
                              abs(out_final$ci.lower - out_final$est.std)
    tmp2 <- abs(out_final$boot.ci.upper - out_final$est.std) /
                              abs(out_final$ci.upper - out_final$est.std)
    tmp1[is.infinite(tmp1) | is.nan(tmp1)] <- NA
    tmp2[is.infinite(tmp2) | is.nan(tmp2)] <- NA
    out_final$ratio.lower <- tmp1
    out_final$ratio.upper <- tmp2
  }

  if (save_boot_est) {
    attr(out_final, "boot_est_ustd") <- out_all
    attr(out_final, "boot_def") <- out_all_def
  }

  class(out_final) <- c("sbt_ustd_boot", class(out))
  fit_summary <- lavaan::summary(object)
  attr(out_final, "pe_attrib") <- attributes(fit_summary$pe)
  attr(out_final, "partable") <- lavaan::parameterTable(object)
  attr(out_final, "level") <- level
  attr(out_final, "type") <- type
  attr(out_final, "boot_ci_type") <- boot_ci_type
  attr(out_final, "call") <- match.call()
  out_final
}