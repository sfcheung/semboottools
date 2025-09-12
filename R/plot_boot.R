#' @title Diagnostic Plots of Bootstrap Estimates in 'lavaan'
#'
#' @description Plots for examining the
#' distribution of bootstrap estimates
#' in a model fitted by `lavaan`.
#'
#' @details Rousselet, Pernet, and Wilcox (2021)
#' argued that when using bootstrapping,
#' it is necessary to examine the distribution
#' of bootstrap estimates. This can be
#' done when [boot::boot()] is used
#' because it has a `plot` method for
#' its output. This cannot be easily
#' done in model fitted by [lavaan::lavaan()],
#' such as [lavaan::sem()] and
#' [lavaan::cfa()].
#'
#' The function [hist_qq_boot()] is used for
#' plotting the distribution of bootstrap
#' estimates for a model fitted by
#' `lavaan` in a format similar to that
#' of the output of [boot::boot()], with
#' a histogram on the left and a normal
#' QQ-plot on the right.
#'
#' For free parameters in a model
#' (unstandardized), it can be called
#' directly on the output of `lavaan`
#' and retrieves the stored estimates.
#'
#' For estimates of user-defined parameters,
#' call [store_boot()] first to compute
#' and store the bootstrap estimates
#' first.
#'
#' For estimates in standardized solution,
#' for both free and user-defined
#' parameters, call [store_boot()]
#' first to compute and store the bootstrap
#' estimates in the standardized solution.
#'
#' It can also
#' plot bootstrap estimates in the output
#' of [standardizedSolution_boot()]
#' or [parameterEstimates_boot()].
#'
#' @references
#' Rousselet, G. A., Pernet, C. R., & Wilcox, R. R. (2021).
#' The percentile bootstrap: A primer with step-by-step
#' instructions in R.
#' *Advances in Methods and Practices in Psychological Science*,
#' *4*(1), 1--10. \doi{10.1177/2515245920911881}
#'
#' @return Return the original
#' [lavaan::lavaan-class] object
#' invisibly. Called for its side-effect
#' (plotting the graphs).
#'
#' @param object Either
#' a `lavaan-class`
#' object with bootstrap estimates
#' stored, or the output of
#' [standardizedSolution_boot()].
#' For standardized solution
#' and user-defined parameters, if
#' the object is a `lavaan-class``
#' object, the
#' estimates need to be stored by
#' [store_boot()].
#'
#' @param param String. The name of
#' the parameter to be plotted, which
#' should be the name as appeared in
#' a call to `coef()`.
#'
#' @param standardized Logical. Whether
#' the estimates from the standardized
#' solution are to be plotted. Default
#' is `NULL`. If `object` is a
#' `lavaan` object, then
#' this is a required parameter
#' and users need to explicitly set it
#' to `TRUE` or `FALSE`. If `object` is
#' the output of
#' [standardizedSolution_boot()],
#' then this argument is ignored (
#' forced to be `TRUE` internally).
#'
#' @param nclass The number of breaks.
#' This argument will be passed to
#' [hist()]. Default is `NULL`.
#'
#' @param hist_color String. The color of the
#' bars in the histogram. It will be
#' passed to [hist()] for the argument
#' `col`. Default is light blue (`scales::alpha("#5DADE2", 0.2)`).
#'
#' @param hist_linewidth The width of
#' the borders of the bars in the
#' histogram. Default is 1.5.
#'
#' @param hist_border_color String.
#' The color of the borders (outline) of the bars
#' in the histogram. It will be passed to [hist()]
#' for the argument `border`.
#' Default is a dark blue color (`"#1B4F72"`).

#' @param density_line_type String.
#' The type of the line of the density
#' curve in the histogram. It will be
#' passed to [lines()] for the argument
#' `lty`. Default is
#' `"solid"`.
#'
#' @param density_line_color String.
#' The color of the density curve in
#' the histogram. It will be
#' passed to [lines()] for the argument
#' `col`. Default is `"blue"`.
#'
#' @param density_line_linewidth The width
#' of the density curve in the histogram.
#' It will be
#' passed to [lines()] for the argument
#' `lwd`.
#' Default is 2.
#'
#' @param est_line_type String. The
#' type of the vertical line in the
#' histogram showing the point estimate
#' of the parameter. It will be
#' passed to [abline()] for the argument
#' `lty`. Default is
#' `"dashed"`,
#'
#' @param est_line_color String. The
#' color of the vertical line showing
#' the point estimate in the histogram.
#' It will be
#' passed to [abline()] for the argument
#' `col`.
#'
#' @param est_line_linewidth The width
#' of the vertical line showing the
#' point estimate in the histogram.
#' It will be
#' passed to [hist()] for the argument
#' `lwd`.  Default is 2.
#'
#' @param qq_dot_size The size of the
#' points in the normal QQ-plot.
#' It will be
#' passed to [qqnorm()] for the argument
#' `cex`. Default is 2.
#'
#' @param qq_dot_color String. The color
#' of the points in the normal QQ-plot.
#' It will be
#' passed to [qqnorm()] for the argument
#' `col`.
#'
#' @param qq_dot_pch Numeric. The shape
#' of the points in the normal QQ-plot.
#' It will be
#' passed to [qqnorm()] for the argument
#' `pch`. Default is 21.
#'
#' @param qq_dot_fill String.
#' The fill color of the points in the normal QQ-plot.
#' Only applicable when `qq_dot_pch` is set to a symbol
#' that allows fill color (e.g., `pch = 21`).
#' It will be passed to [qqnorm()] for the argument `bg`.
#' Default is a semi-transparent light blue
#' (`scales::alpha("#5DADE2", 0.2)`).

#' @param qq_line_linewidth The width
#' of the diagonal line to be drawn in
#' the normal QQ-plot.
#' It will be
#' passed to [qqline()] for the argument
#' `lwd`. Default is 2.1.
#'
#' @param qq_line_color String. The color
#' of the diagonal line to be drawn in
#' the normal QQ-plot.
#' It will be
#' passed to [qqline()] for the argument
#' `col`.
#'
#' @param qq_line_linetype The type of
#' the diagonal line to be drawn in the
#' normal QQ-plot. Default is `"solid"`.
#'
#' @author Shu Fai Cheung <https://orcid.org/0000-0002-9871-9448>
#'
#' @seealso [store_boot()]
#' and [standardizedSolution_boot()].
#'
#' @examples
#'
#' library(lavaan)
#'
#' set.seed(5478374)
#' n <- 50
#' x <- runif(n) - .5
#' m <- .40 * x + rnorm(n, 0, sqrt(1 - .40))
#' y <- .30 * m + rnorm(n, 0, sqrt(1 - .30))
#' dat <- data.frame(x = x, y = y, m = m)
#'
#' mod <-
#' "
#' m ~ a * x
#' y ~ b * m + x
#' ab := a * b
#' "
#' fit <- sem(mod,
#'            data = dat,
#'            se = "bootstrap",
#'            bootstrap = 50,
#'            iseed = 985714)
#'
#' # Can plot bootstrap estimates for
#' # free parameters directly
#' # Note that 'standardized' must be always be set to
#' # either TRUE or FALSE. No default value.
#' hist_qq_boot(fit, "a", standardized = FALSE)
#'
#' # For estimates of user-defined parameters,
#' # call store_boot() first.
#' fit <- store_boot(fit)
#' hist_qq_boot(fit, "ab", standardized = FALSE)
#'
#' # For estimates in standardized solution,
#' # call store_boot() first.
#' fit <- store_boot(fit)
#' hist_qq_boot(fit, "a", standardized = TRUE)
#' hist_qq_boot(fit, "ab", standardized = TRUE)
#'
#' # It can also plot the estimates stored
#' # in the output of standardizedSolution_boot().
#' std_boot <- standardizedSolution_boot(fit)
#' hist_qq_boot(std_boot, "ab")
#' hist_qq_boot(fit, "ab", standardized = TRUE)
#'
#'
#' @importFrom graphics abline hist lines par
#' @importFrom stats qqline qqnorm setNames
#' @export

hist_qq_boot <- function(object,
                         param,
                         standardized = NULL,
                         nclass = NULL,
                         hist_color = "#5DADE233",
                         hist_linewidth = 1.5,
                         hist_border_color = "#1B4F72",
                         density_line_type = "solid",
                         density_line_color = "#8B0000CC",
                         density_line_linewidth = 2,
                         est_line_color = "#154360",
                         est_line_type = "dashed",
                         est_line_linewidth = 2,
                         qq_dot_pch = 21,
                         qq_dot_color = "#1B4F72",
                         qq_dot_fill = "#5DADE233",
                         qq_dot_size = 1.3,
                         qq_line_color = "#8B0000CC",
                         qq_line_linewidth = 2.1,
                         qq_line_linetype = "solid"
) {
  if (is.null(standardized) &&
      !(inherits(object, "sbt_std_boot") ||
        inherits(object, "sbt_ustd_boot"))) {
    stop("'standardized' must be TRUE or FALSE.")
  }
  if (inherits(object, "sbt_std_boot") && is.null(standardized)) {
    standardized <- TRUE
  }
  if (inherits(object, "sbt_ustd_boot") && is.null(standardized)) {
    standardized <- FALSE
  }
  boot_out <- param_find_boot(object = object,
                              param = param,
                              standardized = standardized)
  if (any(is.na(boot_out))) {
    if (standardized) {
      stop("Bootstrap standardized estimates not found or not stored. ",
           "Please call 'store_boot_est_std()' first if ",
           "bootstrapping has been requested.")
    } else {
      stop("Bootstrap estimates not found or not stored.")
    }
  }
  t0 <- boot_out$t0
  t <- boot_out$t
  tmp <- range(t)
  tmp <- tmp[2] - tmp[1]
  if (tmp == 0) {
    stop("Identical estimates in all bootstrap samples.")
  }
  # From plot.boot()
  if (is.null(nclass)) {
    nclass <- min(max(ceiling(length(t) / 25), 10), 100)
  }

  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old))

  ## --- Apply desired settings inside function
  par(family = "serif",   # Times New Roman 系列
      cex.main = 1,
      font.main = 2,
      cex.lab = 0.9,
      cex.axis = 0.8,
      mgp = c(2, 0.6, 0),
      tcl = -0.3,
      lwd = 0.8)


  # From plot.boot()
  #  Calculate the breakpoints for the histogram so that one of them is
  #  exactly t0.
  rg <- range(t)
  rg[1] <- min(rg[1], t0)
  rg[2] <- max(rg[2], t0)
  rg <- rg + 0.05 * c(-1, 1) * diff(rg)
  lc <- diff(rg) / (nclass - 2)
  n1 <- ceiling((t0 - rg[1]) / lc)
  n2 <- ceiling((rg[2] - t0) / lc)
  bks <- t0 + (-n1:n2) * lc
  parold <- par(mfrow = c(1, 2))
  parold2 <- par(lwd = hist_linewidth)
  hist(t,
       probability = TRUE,
       breaks = bks,
       col = hist_color,
       border = hist_border_color,
       main = paste0("Histogram of ",
                     param),
       xlab = param,
       ylab = "Density")
  par(parold2)
  lines(stats::density(t),
        lwd = density_line_linewidth,
        col = density_line_color,
        lty = density_line_type)
  abline(v = t0,
         lwd = est_line_linewidth,
         col = est_line_color,
         lty = est_line_type)
  qqnorm(t,
         cex = qq_dot_size,
         col = qq_dot_color,
         pch = qq_dot_pch,
         bg = qq_dot_fill,
         main = paste0("Normal QQ-Plot of ",
                       param),
         xlab = "Quantiles of Standard Normal",
         ylab = param)
  qqline(t,
         lwd = qq_line_linewidth,
         col = qq_line_color,
         lty = qq_line_linetype
  )
  par(parold)
  invisible(object)
}


# Find and return the stored bootstrap estimates,if available.
#' @noRd

param_find_boot <- function(object,
                            param,
                            standardized) {
  boot_t0 <- NA
  boot_t <- NA
  is_sbt_std_boot <- inherits(object, "sbt_std_boot")
  if (is_sbt_std_boot) {
    standardized <- TRUE
  }
  is_sbt_ustd_boot <- inherits(object, "sbt_ustd_boot")
  if (is_sbt_ustd_boot) {
    standardized <- FALSE
  }
  if (standardized) {
    if (is_sbt_std_boot) {
      coef_names <- lavaan::lav_partable_labels(object)
      boot_i <- attr(object, "boot_est_std")
      if (is.null(boot_i)) {
        stop("Bootstrap estimates not found in the object.")
      }
      colnames(boot_i) <- coef_names
    } else {
      boot_i <- get_boot_est_std(object)
    }
    if (param %in% colnames(boot_i)) {
      if (is_sbt_std_boot) {
        i <- match(param, colnames(boot_i))
        boot_t0 <- object[i, "est.std"]
      } else {
        i <- match(param, std_names(object))
        boot_t0 <- setNames(lavaan::standardizedSolution(object,
                                                         se = FALSE)[i, "est.std"], param)
      }
      boot_t <- boot_i[, param, drop = TRUE]
      out <- list(t0 = boot_t0,
                  t = boot_t)
    }
  } else {
    if (is_sbt_ustd_boot) {
      coef_names <- lavaan::lav_partable_labels(object)
      boot_i <- attr(object, "boot_est_ustd")
      if (is.null(boot_i)) {
        stop("Bootstrap estimates not found in the object.")
      }
    } else {
      boot_i0 <- try(lavaan::lavInspect(object, "boot"), silent = TRUE)
      boot_i1 <- object@external$sbt_boot_ustd
      if (inherits(boot_i0, "try-error") && is.null(boot_i1)) {
        stop("Bootstrapping estimates not found.")
      }
      if (!inherits(boot_i0, "try-error")) {
        boot_i <- boot_i0
      } else {
        boot_i <- boot_i1
      }
    }
    error_idx <- attr(boot_i, "error.idx")
    if (length(error_idx) != 0) {
      boot_i <- boot_i[-error_idx, ]
    }
    if (param %in% colnames(boot_i)) {
      if (is_sbt_ustd_boot) {
        i <- match(param, coef_names)
        boot_t0 <- object[i, "est"]
      } else {
        boot_t0 <- lavaan::coef(object)[param]
      }
      boot_t <- boot_i[, param, drop = TRUE]
    } else {
      boot_i <- get_boot_def(object)
      if (param %in% colnames(boot_i)) {
        if (is_sbt_ustd_boot) {
          i <- match(param, object[, "label"])
          boot_t0 <- object[i, "est"]
        } else {
          boot_t0 <- lavaan::coef(object,
                                  type = "user")[param]
        }
        boot_t <- boot_i[, param, drop = TRUE]
      }
    }
  }
  out <- list(t0 = boot_t0,
              t = boot_t)
  return(out)
}


#' @details
#' The function [scatter_boot()] is
#' used to generate a scatterplot
#' matrix of the bootstrap estimates of
#' two or more parameters. The function
#' [psych::pairs.panels()] from the
#' package `psych` is used.
#'
#' Like [hist_qq_boot()], it can also
#' be used on the output
#' of [standardizedSolution_boot()]
#' or [parameterEstimates_boot()].
#'
#' @param params The vector of the names of
#' the parameters to be plotted, which
#' should be the names as appeared in
#' a call to `coef()`. The function
#' [scatter_boot()] requires two or more
#' parameters selected by this argument.
#'
#' @param main The title of the
#' scatterplot matrix. Default is
#' `"Bootstrap Estimates"`.
#'
#' @param ... Arguments to be passed to
#' [psych::pairs.panels()]. Please refer
#' to the help page of [psych::pairs.panels()]
#' for arguments to customize the
#' plot.
#'
#' @examples
#'
#' # Scatterplot matrix of bootstrap estimates for
#' # two or more free parameters
#' scatter_boot(fit, c("a", "b", "ab"), standardized = FALSE)
#'
#' # Can include user-defined parameters in
#' # scatterplot matrix, if their bootstrap
#' # estimates have been stored
#' scatter_boot(fit, c("ab", "a", "b"), standardized = FALSE)
#'
#' # scatter_boot also supports the
#' # standardized solution
#' scatter_boot(fit, c("a", "b", "ab"), standardized = TRUE)
#'
#' @rdname hist_qq_boot
#' @export

scatter_boot <- function(object,
                         params,
                         standardized = NULL,
                         main = "Bootstrap Estimates",
                         ...) {
  if (is.null(standardized) &&
      !(inherits(object, "sbt_std_boot") ||
        inherits(object, "sbt_ustd_boot"))) {
    stop("'standardized' must be TRUE or FALSE.")
  }
  if (inherits(object, "sbt_std_boot") && is.null(standardized)) {
    standardized <- TRUE
  }
  if (inherits(object, "sbt_ustd_boot") && is.null(standardized)) {
    standardized <- FALSE
  }
  if (length(params) < 2) {
    stop("Need to select two or more parameters.")
  }
  boot_out_list <- sapply(params,
                          param_find_boot,
                          object = object,
                          standardized = standardized,
                          simplify = FALSE,
                          USE.NAMES = TRUE)

  if (any(sapply(boot_out_list, is.na))) {
    if (standardized) {
      stop("Bootstrap standardized estimates not found or not stored. ",
           "Please call 'store_boot_est_std()' first if ",
           "bootstrapping has been requested.")
    } else {
      stop("Bootstrap estimates not found or not stored.")
    }
  }

  t0_tmp <- sapply(boot_out_list,
                   function(x) x$t0,
                   simplify = TRUE,
                   USE.NAMES = TRUE)
  names(t0_tmp) <- names(boot_out_list)
  t_tmp <- sapply(boot_out_list,
                  function(x) x$t,
                  simplify = TRUE,
                  USE.NAMES = TRUE)
  colnames(t_tmp) <- names(boot_out_list)
  boot_out <- list(t0 = t0_tmp,
                   t = t_tmp)

  # # Probably no need for this check
  # tmp <- range(t)
  # tmp <- tmp[2] - tmp[1]
  # if (tmp == 0) {
  #   stop("Identical estimates in all bootstrap samples.")
  # }

  psych::pairs.panels(boot_out$t,
                      main = main,
                      ...)

  invisible(object)
}
