#' Internal function to form bootstrap CI
#' @noRd

# Internal function for different types
# of bootstrap CI.
# Only work for one statistic.
# Adapted from manymome

boot_ci_internal <- function(t0,
                             t,
                             level = .95,
                             boot_type = c("perc", "bc", "bca.simple")) {
  boot_type <- match.arg(boot_type)
  out <- switch(boot_type,
                perc = boot_ci_perc(t0 = t0, t = t, level = level),
                bc = boot_ci_bc(t0 = t0, t = t, level = level),
                bca.simple = boot_ci_bc(t0 = t0, t = t, level = level))
  out
}

# Adapted from manymome
#' @noRd

boot_ci_perc <- function(t0,
                         t,
                         level = .95) {
  boot_tmp <- list(t0 = t0,
                    t = matrix(t, ncol = 1),
                    R = length(t))
  ci <- boot::boot.ci(boot_tmp,
                      type = "perc",
                      conf = level)
  ci$perc[4:5]
}

# Adapted from manymome
#' @noRd

boot_ci_bc <- function(t0,
                       t,
                       level = .95) {
  # Compute the bias correction
  z0_hat <- stats::qnorm(mean(t < t0,
                              na.rm = TRUE))
  z_alpha_l <- stats::qnorm((1 - level) / 2)
  z_alpha_u <- -z_alpha_l
  # BC percentiles
  bc_a_l <- stats::pnorm(2 * z0_hat + z_alpha_l)
  bc_a_u <- stats::pnorm(2 * z0_hat + z_alpha_u)
  # Corresponding confidence level for percentile CIs
  bc_conf_l <- 1 - bc_a_l * 2
  bc_conf_u <- 1 - (1 - bc_a_u) * 2
  boot_ci_bc_l <- boot_ci_perc(t0 = t0,
                                t = t,
                                level = bc_conf_l)
  boot_ci_bc_u <- boot_ci_perc(t0 = t0,
                                t = t,
                                level = bc_conf_u)
  c(boot_ci_bc_l[1], boot_ci_bc_u[2])
}