#' @title Print a 'sbt_std_boot' Object
#'
#' @description Print method for a
#' 'sbt_std_boot' object, which
#' is the output of
#' [standardizedSolution_boot()].
#'
#' @details
#' The default format of the printout
#' is that of [lavaan::standardizedSolution()],
#' which is compact but not easy to
#' read. Users can request a format
#' similar to that of the printout
#' of the summary of a `lavaan` output
#' by setting `output` to `"text"`.
#'
#' For the `"text"` format, users can
#' also select whether
#' only the standardized solution is
#' printed (the default) or whether
#' the standardized solution is appended
#' to the right of the printout.
#'
#' @param x Object of the class
#' `sbt_std_boot`, the output of
#' [standardizedSolution_boot()].
#'
#' @param ... Optional arguments to be
#' passed to [print()] methods.
#'
#' @param nd The number of digits
#' after the decimal place. Default
#' is 3.
#'
#' @param output String. How the results
#' are printed. If set to `"table"`,
#' the results are printed in a table
#' format similar to that of
#' [lavaan::standardizedSolution()].
#' If set to `"text"`, the results will be
#' printed in a text format similar to
#' the printout of the output of
#' [summary()] of
#' a 'lavaan'-class object. Unlike
#' [lavaan::standardizedSolution()],
#' the default is `"text"`.
#'
#' @param standardized_only Logical.
#' If `TRUE`, the default, only the
#' results for the standardized solution
#' will be printed. If `FALSE`, then
#' the standardized solution is printed
#' alongside the unstandardized solution,
#' as in the printout of the output
#' of [summary()] of a 'lavaan'-class
#' object.
#'
#' @seealso [standardizedSolution_boot()]
#'
#' @examples
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
#'            bootstrap = 50)
#' std_out <- standardizedSolution_boot(fit)
#' std_out
#' print(std_out)
#' print(std_out, standardized_only = FALSE)
#'
#' @return
#'  `x` is returned invisibly. Called for its side effect.
#'
#' @author Shu Fai Cheung
#' <https://orcid.org/0000-0002-9871-9448>
#'
#' @export

print.sbt_std_boot <- function(x,
                               ...,
                               nd = 3,
                               output = c("text", "table"),
                               standardized_only = TRUE) {
    output <- match.arg(output)
    x_call <- attr(x, "call")
    if (output == "table") {
        NextMethod()
        return(invisible(x))
      }
    ptable <- attr(x, "partable")
    est0 <- attr(x, "est")
    est1 <- est0
    est1$id <- seq_len(nrow(est1))
    i0 <- colnames(x) %in% c("se", "z", "pvalue",
                             "ci.lower", "ci.upper")
    est1 <- merge(est1,
                  x[, !i0])
    i0 <- colnames(ptable) %in% c("est", "se",
                                  "user", "free",
                                  "ustart", "plabel",
                                  "start",
                                  "id")
    est1 <- merge(est1, ptable[, !i0])
    est1 <- est1[order(est1$id), ]
    est1$id <- NULL
    class(est1) <- class(est0)
    pe_attrib <- attr(x, "pe_attrib")
    tmp <- !(names(pe_attrib) %in% names(attributes(est1)))
    attributes(est1) <- c(attributes(est1),
                          pe_attrib[tmp])
    class(est1) <- c("lavaan.parameterEstimates", class(est1))
    if (!standardized_only) {
        tmp <- colnames(est1)
        tmp[tmp == "est.std"] <- "Standardized"
        tmp[tmp == "boot.ci.lower"] <- "ci.std.lower"
        tmp[tmp == "boot.ci.upper"] <- "ci.std.upper"
        tmp[tmp == "boot.se"] <- "Std.Err.std"
        tmp[tmp == "boot.p"] <- "pvalue.std"
        colnames(est1) <- tmp
        print(est1, ..., nd = nd)
        return(invisible(x))
      } else {
        level <- attr(x, "level")
        est2 <- est1
        est2$est <- est2$est.std
        est2$ci.lower <- est2$boot.ci.lower
        est2$ci.upper <- est2$boot.ci.upper
        est2$se <- est2$boot.se
        est2$boot.se <- NULL
        est2$z <- NULL
        # if (!is.null(est2$boot.p)) {
        #   est2$pvalue <- est2$boot.p
        # } else {
        #   est2$pvalue <- NULL
        # }
        est2$est.std <- NULL
        est2$boot.ci.lower <- NULL
        est2$boot.ci.upper <- NULL
        out <- utils::capture.output(print(est2, nd = nd))
        i <- grepl("Parameter Estimates:", out, fixed = TRUE)
        out[i] <- "Standardized Estimates Only"
        i <- grepl("  Standard errors  ", out, fixed = TRUE)
        j <- unlist(gregexpr("Bootstrap", out[i]))[1]
        tmp <- "  Confidence interval"
        st1 <- paste0(tmp,
                      paste0(rep(" ", j - nchar(tmp) - 1),
                             collapse = ""),
                      "Bootstrap")
        j <- nchar(out[i])
        tmp <- "  Confidence Level"
        tmp2 <- paste0(formatC(level * 100, digits = 1, format = "f"),
                       "%")
        st2 <- paste0(tmp,
                      paste0(rep(" ", j - nchar(tmp) - nchar(tmp2)),
                             collapse = ""),
                      tmp2)
        tmp <- "  Bootstrap CI Type"
        tmp2 <- switch(attr(x, "boot_ci_type"),
                       perc = "Percentile",
                       bc = "Bias-Corrected",
                       bca.simple = "Bias-Corrected")
        st2b <- paste0(tmp,
                       paste0(rep(" ", j - nchar(tmp) - nchar(tmp2)),
                              collapse = ""),
                       tmp2)
        if (!is.null(est2$boot.p)) {
          tmp <- "  Bootstrap P-Value"
          tmp2 <- "Asymmetric P-Value"
          st2c <- paste0(tmp,
                        paste0(rep(" ", j - nchar(tmp) - nchar(tmp2)),
                                collapse = ""),
                        tmp2)
        } else {
          st2c <- NULL
        }
        tmp <- "  Standardization Type"
        tmp2 <- attr(x, "type")
        st3 <- paste0(tmp,
                      paste0(rep(" ", j - nchar(tmp) - nchar(tmp2)),
                             collapse = ""),
                      tmp2)
        out <- c(out[seq_len(which(i))],
                 st1,
                 st2,
                 st2b,
                 st2c,
                 st3,
                 out[-seq_len(which(i))])
        out <- gsub("    Estimate  Std.Err",
                    "Standardized  Std.Err",
                    out)
        cat(out, sep = "\n")
        return(invisible(x))
      }
  }