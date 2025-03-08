# Generate bootstrap estimates
# Called by:
# - store_boot()

boot_est_std <- function(object,
                         type,
                         ...) {
    # For lavaan 0.6-13
    # Remove bootstrap replications with error
    boot_est0 <- try(lavaan::lavTech(object, "boot"), silent = TRUE)
    boot_est1 <- object@external$sbt_boot_ustd
    if (inherits(boot_est0, "try-error") && is.null(boot_est1)) {
        stop("Bootstrapping estimates not found. Was se = 'boot' or 'bootstrap'?")
      }
    if (inherits(boot_est0, "try-error")) {
      boot_est0 <- boot_est1
    }
    boot_error_idx <- attr(boot_est0, "error.idx")
    if (!is.null(boot_error_idx)) {
        if (length(boot_error_idx) > 0) {
            boot_est0 <- boot_est0[-boot_error_idx, ]
          }
      }
    std_args <- list(...)
    ptable <- lavaan::parameterTable(object)
    p_free <- ptable$free > 0
    p_est  <- ptable$est
    boot_est <- split(boot_est0, row(boot_est0))
    out_all <- t(sapply(boot_est, std_i,
                        p_est = p_est,
                        p_free = p_free,
                        object = object,
                        type = type,
                        std_args = std_args))
    return(out_all)
  }

# Generate the function for bootstrapping.
# Called by:
# - boot_est_std

std_i <- function(est_i,
                  p_est,
                  p_free,
                  object,
                  std_args,
                  type) {
  p_est[p_free] <- est_i
  GLIST_i <- lavaan::lav_model_set_parameters(object@Model,
                                              est_i)@GLIST
  std_args1 <- utils::modifyList(std_args,
                                  list(object = object,
                                       type = type,
                                       est = p_est,
                                       GLIST = GLIST_i,
                                       se = FALSE,
                                       zstat = FALSE,
                                       pvalue = FALSE,
                                       ci = FALSE,
                                       output = "data.frame"))
  do.call(lavaan::standardizedSolution, std_args1)$est.std
}


# Compute bootstrap estimates of user-defined parameters
# Called by:
# - store_boot

boot_def <- function(object) {
  # For lavaan 0.6-13
  # Remove bootstrap replications with error
  if (!(":=" %in% lavaan::parameterTable(object)$op)) {
      return(NULL)
  }
  boot_est0 <- try(lavaan::lavTech(object,
                                   "boot"),
                   silent = TRUE)
  boot_est1 <- object@external$sbt_boot_ustd
  if (inherits(boot_est0, "try-error") && is.null(boot_est1)) {
      stop("Bootstrapping estimates not found. Was se = 'boot' or 'bootstrap'?")
  }
  if (inherits(boot_est0, "try-error")) {
    boot_est0 <- boot_est1
  }
  boot_error_idx <- attr(boot_est0,
                         "error.idx")
  if (!is.null(boot_error_idx)) {
    if (length(boot_error_idx) > 0) {
      boot_est0 <- boot_est0[-boot_error_idx, , drop = FALSE]
    }
  }
  boot_est <- split(boot_est0,
                    row(boot_est0))
  out_all <- lapply(boot_est,
                    object@Model@def.function)
  out_all <- do.call(rbind,
                     out_all)
  return(out_all)
}

# Called by:
# - plot_boot()

get_boot_est_std <- function(object) {
  return(object@external$sbt_boot_std)
}

# Called by:
# - plot_boot()

get_boot_def <- function(object) {
    return(object@external$sbt_boot_def)
  }

# Generate names for standardized solution
# Called by:
# - plot_boot()

std_names <- function(object, ...) {
  std <- lavaan::standardizedSolution(object, se = FALSE, ...)
  std$id <- seq_len(nrow(std))
  ptable <- lavaan::parameterTable(object)
  std1 <- merge(std, ptable,
                all.y = FALSE)
  std1 <- std1[order(std1$id), ]
  std1$lavlabel <- lavaan::lav_partable_labels(std1,
                      blocks = c("group", "level"),
                      group.equal = "",
                      group.partial = "",
                      type = "user")
  return(std1$lavlabel)
}

check_std_i <- function(object, type, std_args) {
  # Work-in-progress
  # Not used for now
  # Do one bootstrap with bootstrapLavaan(),
  #   with est and std
  # Put est as boot, and see if std_i can reproduce std
  fct <- function(fit, std_type, std_args) {
      args0 <- utils::modifyList(std_args,
                                  list(object = fit,
                                      type = std_type,
                                      se = FALSE,
                                      zstat = FALSE,
                                      pvalue = FALSE,
                                      ci = FALSE,
                                      output = "data.frame"))
      list(coef = lavaan::coef(fit),
            est.std = do.call(lavaan::standardizedSolution, args0)$est.std)
    }
  object_noboot <- lavaan::update(object, se = "none")
  out_test <- lavaan::bootstrapLavaan(object_noboot,
                                      R = 1,
                                      type = "ordinary",
                                      FUN = fct,
                                      warn = -1L,
                                      std_type = type,
                                      std_args = std_args)
  object_test <- object
  object_test@boot$coef <- out_test[[1]]
  ptable <- lavaan::parameterTable(object)
  boot_std_test <- std_i(est_i = out_test[[1]],
                          p_est = ptable$est,
                          p_free = ptable$free > 0,
                          object = object,
                          std_args = std_args,
                          type = type)
  if (!isTRUE(all.equal(boot_std_test, out_test[[2]]))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
