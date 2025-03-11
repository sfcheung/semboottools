#' @title Compute and Store Bootstrap Estimates
#'
#' @description This function computes bootstrap estimates
#' of a fitted structural equation model
#' and stores the estimates for further
#' processing.
#'
#' @details
#'
#' The function [store_boot()]
#' receives a
#' [lavaan::lavaan-class] object, optionally
#' fitted with bootstrapping standard errors
#' requested, and compute and store
#' the bootstrap estimates of user-defined
#' parameters and estimates in the
#' standardized solution.
#'
#' If bootstrapping was not requested
#' when fitting the model (i.e., `se`
#' not set to `"boot"` or `"bootstrap"`),
#' then bootstrapping will be conducted
#' using [lavaan::bootstrapLavaan()] to
#' compute bootstrap estimates of free
#' parameters. Otherwise, the stored
#' bootstrap estimates will be used in
#' subsequent steps.
#'
#' For standardized solution bootstrap
#' estimates, it works by calling
#' [lavaan::standardizedSolution()]
#' with the bootstrap estimates
#' of free parameters in each bootstrap sample
#' to compute the standardized estimates
#' in each sample.
#'
#' For user-defined parameters, it works
#' by calling the function used to
#' compute user-defined parameters with
#' the bootstrap estimates of free
#' parameters in each bootstrap samples
#' to compute the user-defined parameters.
#'
#' The bootstrap estimates are then
#' stored in the `external` slot
#' of the fit object for further
#' processing.
#'
#' @return
#' The original `lavaan` object is
#' returned with the following objects
#' stored in the `external` slot:
#'
#' - `sbt_boot_std`: The matrix of
#' bootstrap estimates in the
#' standardized solution.
#'
#' - `sbt_boot_def`: The matrix of
#' bootstrap estimates of user-defined
#' parameters, if any.
#'
#' - `sbt_boot_ustd`: The matrix of
#' bootstrap estimates of free
#' parameters, if bootstrapping is
#' not requested when fitting the
#' model (i.e., `se` is not set to
#' `"boot"` or `"bootstrap"` when
#' fitting the model in `lavaan`).
#'
#' @param object A 'lavaan'-class
#' object, fitted with 'se = "boot"'.
#'
#' @param type The type of standard
#' estimates. The same argument of
#' [lavaan::standardizedSolution()],
#' and support all values supported by
#' [lavaan::standardizedSolution()].
#' Default is `"std.all"`.
#'
#' @param do_bootstrapping If `TRUE` and
#' bootstrapping was not requested when
#' fitting the model, bootstrapping
#' will be done using
#' [lavaan::bootstrapLavaan()]. Default
#' is `TRUE`.
#'
#' @param R If [lavaan::bootstrapLavaan()]
#' is called (see `do_bootstrapping`),
#' this is the number of bootstrap
#' samples, to be used by
#' [lavaan::bootstrapLavaan()].
#'
#' @param boot_type If [lavaan::bootstrapLavaan()]
#' is called (see `do_bootstrapping`),
#' this is type of bootstrapping,
#' to be passed to the argument `type`
#' of
#' [lavaan::bootstrapLavaan()].
#' Default is `"ordinary"`. See the
#' help page of [lavaan::bootstrapLavaan()]
#' for details.
#'
#' @param parallel If [lavaan::bootstrapLavaan()]
#' is called (see `do_bootstrapping`),
#' whether parallel processing will
#' be used.
#' to be passed to the argument of the
#' same name in [lavaan::bootstrapLavaan()].
#' Default is `"no"`. Can be
#' `"snow"` or `"multicore"`. See the
#' help page of [lavaan::bootstrapLavaan()]
#' for details.
#'
#' @param ncpus If [lavaan::bootstrapLavaan()]
#' is called (see `do_bootstrapping`),
#' and parallel processing is to be used,
#' this is the number of CPU cores to
#' use, to be passed to the argument of the
#' same name in [lavaan::bootstrapLavaan()].
#' Default is `parallel::detectCores(logical = FALSE) - 1`,
#' the number of physical cores minus 1,
#' different from the default of
#' [lavaan::bootstrapLavaan()] but identical
#' to the default of [lavaan::sem()] and
#' [lavaan::cfa()].
#'
#' @param iseed If [lavaan::bootstrapLavaan()]
#' is called (see `do_bootstrapping`),
#' this should be an integer used to
#' generate reproducible bootstrap
#' results, to be passed to the argument of the
#' same name in [lavaan::bootstrapLavaan()].
#' Default is `NULL` but it should nearly
#' always be set to an arbitrary integer.
#' See the
#' help page of [lavaan::bootstrapLavaan()]
#' for details.
#'
#' @param keep.idx Whether the indices
#' of cases selected in each bootstrap
#' sample is to be stored. To be passed
#' to the argument of the same name
#' in [lavaan::bootstrapLavaan()].
#' Default is `FALSE`.
#'
#' @param bootstrapLavaan_args A named
#' list of additional arguments to be
#' passed to [lavaan::bootstrapLavaan()].
#' Note that the other arguments in
#' [store_boot()] takes precedence,
#' overriding arguments of the same
#' names in this list, if any.
#'
#' @author Shu Fai Cheung
#' <https://orcid.org/0000-0002-9871-9448>.
#' Based on [semhelpinghands::standardizedSolution_boot_ci()]，
#' which was originally proposed in an issue at GitHub
#' <https://github.com/simsem/semTools/issues/101#issue-1021974657>,
#' inspired by a discussion at
#' the Google group for lavaan
#' <https://groups.google.com/g/lavaan/c/qQBXSz5cd0o/m/R8YT5HxNAgAJ>.
#' Unlike [semhelpinghands::standardizedSolution_boot_ci()],
#' this function only computes and stores
#' the bootstrap estimates.
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
#' fit <- store_boot(fit)
#'
#' @export

store_boot <- function(object,
                       type = "std.all",
                       do_bootstrapping = TRUE,
                       R = 1000,
                       boot_type = "ordinary",
                       parallel = c("no", "multicore", "snow"),
                       ncpus = parallel::detectCores(logical = FALSE) - 1,
                       iseed = NULL,
                       keep.idx = FALSE,
                       bootstrapLavaan_args = list()) {
  parallel <- match.arg(parallel)
  if (!inherits(object, "lavaan")) {
    stop("The object must be a lavaan-class object.")
  }
  boot_est0 <- try(lavaan::lavTech(object, "boot"),
                   silent = TRUE)
  if (inherits(boot_est0, "try-error") && !do_bootstrapping) {
    stop("Bootstrapping estimates not found ",
         "probably because se is not 'boot' or 'bootstrap'. ",
         "Call store_boot() with do_bootstrapping = TRUE.")
  }
  if (inherits(boot_est0, "try-error") && do_bootstrapping) {
    bootstrapLavaan_args1 <- utils::modifyList(bootstrapLavaan_args,
                                               list(object = object,
                                                    R = R,
                                                    type = boot_type,
                                                    parallel = parallel,
                                                    ncpus = ncpus,
                                                    iseed = iseed,
                                                    keep.idx = keep.idx))
    boot_ustd_out <- do.call(lavaan::bootstrapLavaan,
                             bootstrapLavaan_args1)
    object@external$sbt_boot_ustd <- boot_ustd_out
  }

  boot_std_out <- boot_est_std(object = object,
                               type = type)
  object@external$sbt_boot_std <- boot_std_out
  object@external$sbt_boot_std_type <- type
  boot_def_out <- boot_def(object = object)
  object@external$sbt_boot_def <- boot_def_out
  # call is stored as a backup plan
  object@external$sbt_args <- list(R = R,
                                   boot_type = boot_type,
                                   parallel = parallel,
                                   ncpus = ncpus,
                                   iseed = iseed,
                                   bootstrapLavaan_args = bootstrapLavaan_args,
                                   call = match.call())
  object
}