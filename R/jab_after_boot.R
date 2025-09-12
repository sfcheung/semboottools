#' Jackknife-after-Bootstrap (JAB) influence for lavaan bootstrap results
#'
#' @description
#' Compute Jackknife-after-Bootstrap (JAB) influence values for a single
#' model parameter and diagnose observation-level influence by comparing
#' the full bootstrap distribution with leave-one-out (LOO) subdistributions
#' obtained via \code{boot.idx} (stored by \code{store_boot(keep.idx = TRUE)}).
#'
#' The ALL-summary (mean/SE/CI) is taken directly and consistently from
#' \code{standardizedSolution_boot()} (standardized case) or
#' \code{parameterEstimates_boot()} (unstandardized case). JAB influence is
#' defined as \eqn{I_j = n\{\bar\theta_{(-j)} - \bar\theta\}}, where
#' \eqn{\bar\theta} is the bootstrap mean over all replicates, and
#' \eqn{\bar\theta_{(-j)}} is the bootstrap mean over the replicates that
#' exclude observation \eqn{j}.
#'
#' @param fit A `lavaan` object for which `store_boot(keep.idx = TRUE)` has been called.
#'   Must contain `fit@external$sbt_boot_ustd` (with `boot.idx` attribute) and typically `sbt_boot_std`.
#' @param param Character(1). Target parameter. Accepts `"lhs op rhs"`, a `":="` label, or a free-parameter label.
#' @param standardized Logical. If `TRUE`, use standardized bootstrap; else unstandardized.
#' @param top_k Integer. Number of top cases (by `|JAB_value|`) to report.
#' @param ci_level Numeric (0,1). Percentile CI level for LOO subdistributions.
#' @param min_keep Integer. Minimum number of bootstrap replicates kept in each LOO subset;
#'   default `max(30, floor(0.2 * B))`.
#' @param plot Logical. If `TRUE`, draw a JAB diagnostic plot.
#' @param plot_engine Character. Plot engine for `plot = TRUE`. One of `"ggplot2"` or `"base"`.   # <-- 新增
#'   Default `"ggplot2"` if available, otherwise `"base"`.
#' @param ylab_override Optional character. Override the default y-axis label in the plot.        # <-- 新增
#' @param verbose Logical. If `TRUE`, print compact summaries.
#' @param font_family Character. Graphics font family (e.g., `"serif"`, `"sans"`, `"Times New Roman"`).
#'
#'
#' @returns A list with:
#' \describe{
#'   \item{\code{param}}{The target parameter string.}
#'   \item{\code{standardized}}{Logical flag as input.}
#'   \item{\code{full_summary}}{Data frame with ALL distribution summary:
#'     \code{scope="ALL"}, \code{param}, \code{mean}, \code{SE}, \code{CI.Lo}, \code{CI.Up}.}
#'   \item{\code{cases_summary}}{Data frame (top \code{top_k}) with columns:
#'     \code{case}, \code{JAB_value}, \code{mean}, \code{SE}, \code{CI.Lo}, \code{CI.Up}.}
#'   \item{\code{F}}{The occurrence matrix \eqn{B \times n}.}
#'   \item{\code{tstar}}{Full bootstrap vector for \code{param}.}
#' }
#'
#' @details
#' Parameter identification is robust: it accepts \code{"lhs op rhs"},
#' \code{":="} user-defined name, and \emph{label}. For unstandardized cases,
#' if \code{label} is given and the column is missing in
#' \code{attr(unstd_boot, "boot_est_ustd")}, the function will fall back to
#' the raw matrix \code{fit@external$sbt_boot_ustd} if a column with that label
#' exists (e.g., \code{"b"}), keeping SE/CI from the \code{parameterEstimates_boot()}
#' row if available.
#'
#' @seealso \code{\link{store_boot}}, \code{\link{standardizedSolution_boot}},
#'   \code{\link{parameterEstimates_boot}}
#' @importFrom rlang .data
#' @export
jab_after_boot <- function(
    fit,
    param,
    standardized = TRUE,
    top_k = 5L,
    ci_level = 0.95,
    min_keep = NULL,
    plot = FALSE,
    plot_engine = c("ggplot2","base"),
    ylab_override = NULL,
    verbose = TRUE,
    font_family = "sans"
){
  stopifnot(inherits(fit, "lavaan"))
  # --- need boot.idx from unstd bootstrap ---
  boot_ustd_raw <- fit@external$sbt_boot_ustd
  if (is.null(boot_ustd_raw)) {
    stop("fit@external$sbt_boot_ustd is NULL; call store_boot(keep.idx = TRUE) first.")
  }

  # Build B x n occurrence matrix F
  idx_list <- attr(boot_ustd_raw, "boot.idx")
  if (is.null(idx_list) || !is.list(idx_list) || length(idx_list) != 1L) {
    stop("sbt_boot_ustd lacks attribute 'boot.idx'; ensure keep.idx = TRUE in store_boot().")
  }
  idxmat <- idx_list[[1]]
  if (!is.matrix(idxmat)) stop("'boot.idx' is not a matrix.")
  B <- nrow(idxmat); n <- ncol(idxmat)
  Fmat <- matrix(0L, nrow = B, ncol = n)
  for (b in seq_len(B)) Fmat[b, ] <- tabulate(idxmat[b, ], nbins = n)
  if (is.null(min_keep)) min_keep <- max(30L, floor(0.2 * B))

  # ---- ALL summary & full tstar (consistent with std/unstd boot objects) ----
  if (isTRUE(standardized)) {
    std_boot <- standardizedSolution_boot(fit)
    all_out  <- .get_all_from_std_boot(std_boot, param, ci_level = ci_level)
    level_all <- default_if_null(attr(std_boot, "level"), ci_level)
  } else {
    unstd_boot <- parameterEstimates_boot(fit)
    all_out    <- .get_all_from_unstd_boot(unstd_boot, param, ci_level = ci_level,
                                           boot_ustd_raw = boot_ustd_raw)
    level_all <- default_if_null(attr(unstd_boot, "level"), ci_level)
  }
  tstar <- all_out$tstar
  full_summary <- all_out$summary
  full_summary$scope <- "ALL"  # normalized scope title
  mean_all <- full_summary$mean[1]

  # ---- JAB over leave-one-out subdistributions ----
  q_lo <- (1 - level_all) / 2; q_hi <- 1 - q_lo
  mean_excl <- se_excl <- ci_lo_excl <- ci_up_excl <- rep(NA_real_, n)
  I_vals <- rep(NA_real_, n)

  for (j in seq_len(n)) {
    keep <- which(Fmat[, j] == 0L)
    if (length(keep) < min_keep) next
    tj <- tstar[keep]
    mean_excl[j]  <- mean(tj, na.rm = TRUE)
    se_excl[j]    <- stats::sd(tj, na.rm = TRUE)
    ci_j          <- stats::quantile(tj, c(q_lo, q_hi), na.rm = TRUE)
    ci_lo_excl[j] <- as.numeric(ci_j[1])
    ci_up_excl[j] <- as.numeric(ci_j[2])
    I_vals[j]     <- n * (mean_excl[j] - mean_all) # JAB influence
  }

  ord <- order(abs(I_vals), decreasing = TRUE, na.last = NA)
  top_idx <- utils::head(ord, top_k)

  cases_summary <- data.frame(
    case      = top_idx,
    JAB_value = I_vals[top_idx],
    mean      = mean_excl[top_idx],
    SE        = se_excl[top_idx],
    CI.Lo     = ci_lo_excl[top_idx],
    CI.Up     = ci_up_excl[top_idx]
  )
  rownames(cases_summary) <- NULL

  if (isTRUE(verbose)) {
    cat("\n=== Full-sample bootstrap summary (ALL) ===\n")
    print(full_summary, row.names = FALSE)

    cat("\n=== Leave-one-out (LOO) summaries ===\n")
    print(cases_summary, row.names = FALSE)
  }


  plot_obj <- NULL
  plot_engine <- match.arg(plot_engine)

  if (isTRUE(plot)) {
    if (identical(plot_engine, "ggplot2")) {
      plot_obj <- .make_jab_plot_gg(
        mean_excl = mean_excl,
        mean_all  = mean_all,
        top_idx   = top_idx,
        param     = param,
        ylab_override = ylab_override,
        font_family   = font_family
      )
      if (is.null(plot_obj)) {
        warning("No valid leave-one-out means (all NA); skip plotting.")
      } else {
        print(plot_obj)
      }
    } else {
      # --------- 原 base 图的简化版（y 轴文案更短、图例移出绘图区右侧） ----------
      opar <- graphics::par(no.readonly = TRUE)
      on.exit({
        ro <- c("cin","cra","csi","cxy","din","page","pin")
        ok <- setdiff(names(opar), ro)
        graphics::par(opar[ok])
      }, add = TRUE)

      valid <- which(!is.na(mean_excl))
      if (length(valid) == 0L) {
        warning("No valid leave-one-out means (all NA); skip plotting.")
      } else {
        y_all <- mean_excl[valid]
        ylim <- range(c(y_all, mean_all), na.rm = TRUE)
        d <- diff(ylim)
        if (is.finite(d) && d == 0) {
          bump <- max(1e-6, abs(mean_all) * 1e-3)
          ylim <- ylim + c(-1, 1) * bump
        }

        # 让图例在图外：扩大右边距，并用 xpd=NA
        graphics::par(family = font_family, xpd = NA, mar = graphics::par("mar") + c(0,0,0,3))

        ylab_txt <- if (is.null(ylab_override)) sprintf("LOO bootstrap mean of %s", param) else ylab_override
        graphics::plot(valid, y_all, pch = 16, cex = 0.7,
                       xlab = "Observation index",
                       ylab = ylab_txt,
                       main = "Jackknife-after-Bootstrap diagnostic",
                       ylim = ylim)
        graphics::abline(h = mean_all, col = "red", lwd = 2)

        top_idx_ok <- top_idx[!is.na(top_idx) & !is.na(mean_excl[top_idx])]
        if (length(top_idx_ok)) {
          graphics::points(top_idx_ok, mean_excl[top_idx_ok], pch = 16, cex = 1.0, col = "blue")
          graphics::legend("right", inset = c(-0.12, 0), xpd = NA,
                           legend = c("Full-sample bootstrap mean", "Top influential cases"),
                           col = c("red", "blue"),
                           lwd = c(2, NA), pch = c(NA, 16), bty = "n", y.intersp = 1)
        } else {
          graphics::legend("right", inset = c(-0.12, 0), xpd = NA,
                           legend = "Full-sample bootstrap mean",
                           col = "red", lwd = 2, bty = "n", y.intersp = 1)
        }
      }
    }
  }


  invisible(list(
    param = param,
    standardized = standardized,
    full_summary = full_summary,
    cases_summary = cases_summary,
    F = Fmat,
    tstar = tstar,
    plot_obj = if (isTRUE(plot)) plot_obj else NULL
  ))
}

# ---------------- internal helpers (not exported) ----------------

#' @keywords internal
default_if_null <- function(x, default) if (is.null(x)) default else x

#' @keywords internal
.parse_lhs_op_rhs <- function(s){
  m <- regexec("^\\s*([^~:=]+)\\s*(=~|~~|~)\\s*([^~:=]+)\\s*$", s)
  r <- regmatches(s, m)[[1]]
  if (length(r) == 4L) list(lhs = trimws(r[2]), op = r[3], rhs = trimws(r[4])) else NULL
}

#' @keywords internal
.get_all_from_std_boot <- function(std_boot, param, ci_level = NULL) {
  stopifnot(is.data.frame(std_boot))
  boots <- attr(std_boot, "boot_est_std")
  if (is.null(boots)) stop("standardizedSolution_boot() output lacks attribute 'boot_est_std'.")
  if (is.null(ci_level)) ci_level <- default_if_null(attr(std_boot, "level"), 0.95)

  # 先尝试：param 本身是列名
  col_ok <- if (!is.null(colnames(boots)) && (param %in% colnames(boots))) param else NULL

  # 其次：param 是 "lhs op rhs"
  if (is.null(col_ok)) {
    p3 <- .parse_lhs_op_rhs(param)
    if (!is.null(p3)) {
      key <- paste0(p3$lhs, p3$op, p3$rhs)  # e.g., "Y~M", "Y~~Y", "X=~x1"
      if (key %in% colnames(boots)) col_ok <- key
    }
  }

  # 再次兜底：param 是 label（非 :=）
  if (is.null(col_ok)) {
    row_lab <- which(std_boot$label == param & std_boot$op != ":=")
    if (length(row_lab) == 1L) {
      lhs <- std_boot$lhs[row_lab]; op <- std_boot$op[row_lab]; rhs <- std_boot$rhs[row_lab]
      key <- paste0(lhs, op, rhs)
      if (key %in% colnames(boots)) col_ok <- key
    }
  }

  if (is.null(col_ok)) stop(sprintf("Parameter '%s' not found in 'boot_est_std' columns or by label.", param))

  tstar   <- as.numeric(boots[, col_ok])
  mean_all <- mean(tstar, na.rm = TRUE)
  se_all   <- stats::sd(tstar, na.rm = TRUE)

  # 在 std_boot 的行中定位（优先 lhs/op/rhs，其次 label）
  row_id <- integer(0)
  p3 <- .parse_lhs_op_rhs(col_ok)
  if (!is.null(p3)) {
    row_id <- which(std_boot$lhs == p3$lhs & std_boot$op == p3$op & std_boot$rhs == p3$rhs)
  }
  if (length(row_id) != 1L) row_id <- which(std_boot$label == col_ok)
  if (length(row_id) != 1L) {
    # 如果 col_ok 是由 label 映射来的，按 label 的那一行定位
    row_id <- which(std_boot$label == param & std_boot$op != ":=")
  }
  if (length(row_id) != 1L) stop(sprintf("Cannot locate '%s' row in standardizedSolution_boot().", param))

  ci_lo <- as.numeric(std_boot$boot.ci.lower[row_id])
  ci_up <- as.numeric(std_boot$boot.ci.upper[row_id])

  list(
    tstar = tstar,
    summary = data.frame(scope = "ALL(std_boot)", param = param,
                         mean = mean_all, SE = se_all, CI.Lo = ci_lo, CI.Up = ci_up)
  )
}

#' @keywords internal
.get_all_from_unstd_boot <- function(unstd_boot, param, ci_level = NULL, boot_ustd_raw = NULL) {
  stopifnot(is.data.frame(unstd_boot))
  boots_ustd <- attr(unstd_boot, "boot_est_ustd")  # free params
  boots_def  <- attr(unstd_boot, "boot_def")       # := params
  if (is.null(ci_level)) ci_level <- default_if_null(attr(unstd_boot, "level"), 0.95)
  q_lo <- (1 - ci_level) / 2; q_hi <- 1 - q_lo

  # --- Case A: free parameter "lhs op rhs" ---
  p3 <- .parse_lhs_op_rhs(param)
  if (!is.null(p3)) {
    if (is.null(boots_ustd)) stop("parameterEstimates_boot() output lacks 'boot_est_ustd'.")
    key <- paste0(p3$lhs, p3$op, p3$rhs)
    if (!(key %in% colnames(boots_ustd)))
      stop(sprintf("Column '%s' not found in 'boot_est_ustd'.", key))

    tstar <- as.numeric(boots_ustd[, key])
    mean_all <- mean(tstar, na.rm = TRUE)
    se_all   <- stats::sd(tstar, na.rm = TRUE)

    row_id <- which(unstd_boot$lhs == p3$lhs & unstd_boot$op == p3$op & unstd_boot$rhs == p3$rhs)
    if (length(row_id) != 1L)
      stop(sprintf("Cannot uniquely locate row for '%s %s %s' in parameterEstimates_boot.", p3$lhs, p3$op, p3$rhs))
    ci_lo <- as.numeric(unstd_boot$boot.ci.lower[row_id])
    ci_up <- as.numeric(unstd_boot$boot.ci.upper[row_id])

    return(list(
      tstar   = tstar,
      summary = data.frame(scope = "ALL(unstd_boot)", param = param,
                           mean = mean_all, SE = se_all, CI.Lo = ci_lo, CI.Up = ci_up)
    ))
  }

  # --- Case B: := user-defined (e.g., "ab") ---
  if (!is.null(boots_def) && !is.null(colnames(boots_def)) && (param %in% colnames(boots_def))) {
    tstar <- as.numeric(boots_def[, param])
    mean_all <- mean(tstar, na.rm = TRUE)
    se_all   <- stats::sd(tstar, na.rm = TRUE)

    row_id <- which(unstd_boot$label == param & unstd_boot$op == ":=")
    if (length(row_id) == 1L &&
        !is.na(unstd_boot$boot.ci.lower[row_id]) && !is.na(unstd_boot$boot.ci.upper[row_id])) {
      ci_lo <- as.numeric(unstd_boot$boot.ci.lower[row_id])
      ci_up <- as.numeric(unstd_boot$boot.ci.upper[row_id])
    } else {
      qs <- stats::quantile(tstar, probs = c(q_lo, q_hi), na.rm = TRUE)
      ci_lo <- as.numeric(qs[1]); ci_up <- as.numeric(qs[2])
    }

    return(list(
      tstar   = tstar,
      summary = data.frame(scope = "ALL(unstd_boot)", param = param,
                           mean = mean_all, SE = se_all, CI.Lo = ci_lo, CI.Up = ci_up)
    ))
  }

  # --- Case C: label (e.g., "b") ---
  row_lab <- which(unstd_boot$label == param & unstd_boot$op != ":=")
  if (length(row_lab) == 1L) {
    lhs <- unstd_boot$lhs[row_lab]; op <- unstd_boot$op[row_lab]; rhs <- unstd_boot$rhs[row_lab]
    key <- paste0(lhs, op, rhs)

    if (!is.null(boots_ustd) && (key %in% colnames(boots_ustd))) {
      tstar <- as.numeric(boots_ustd[, key])
    } else if (!is.null(boot_ustd_raw) && !is.null(colnames(boot_ustd_raw)) && (param %in% colnames(boot_ustd_raw))) {
      # fallback to raw matrix if it indeed has a column named by label (e.g., "b")
      tstar <- as.numeric(boot_ustd_raw[, param])
    } else {
      stop(sprintf("Neither 'boot_est_ustd' has column '%s' nor raw matrix has column '%s'.", key, param))
    }

    mean_all <- mean(tstar, na.rm = TRUE)
    se_all <- if (!is.na(unstd_boot$boot.se[row_lab])) as.numeric(unstd_boot$boot.se[row_lab]) else stats::sd(tstar, na.rm = TRUE)
    if (!is.na(unstd_boot$boot.ci.lower[row_lab]) && !is.na(unstd_boot$boot.ci.upper[row_lab])) {
      ci_lo <- as.numeric(unstd_boot$boot.ci.lower[row_lab])
      ci_up <- as.numeric(unstd_boot$boot.ci.upper[row_lab])
    } else {
      qs <- stats::quantile(tstar, probs = c(q_lo, q_hi), na.rm = TRUE)
      ci_lo <- as.numeric(qs[1]); ci_up <- as.numeric(qs[2])
    }

    return(list(
      tstar   = tstar,
      summary = data.frame(scope = "ALL(unstd_boot)", param = param,
                           mean = mean_all, SE = se_all, CI.Lo = ci_lo, CI.Up = ci_up)
    ))
  }

  stop(sprintf(
    "Could not resolve '%s': not a 'lhs op rhs' free parameter, not in ':=' names, and no matching label row.",
    param
  ))
}


#' @keywords internal
.make_jab_plot_gg <- function(
    mean_excl, mean_all, top_idx, param, ylab_override = NULL, font_family = "serif"
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_engine='ggplot2'. Install via install.packages('ggplot2').")
  }

  valid <- which(!is.na(mean_excl))
  if (length(valid) == 0L) return(NULL)

  df <- data.frame(
    idx   = valid,
    y     = mean_excl[valid],
    group = "Other cases",
    stringsAsFactors = FALSE
  )
  top_ok <- top_idx[!is.na(top_idx) & top_idx %in% valid]
  if (length(top_ok)) df$group[df$idx %in% top_ok] <- "Top influential cases"
  df$group <- factor(df$group, levels = c("Other cases", "Top influential cases"))

  # 均值线（进入线型图例）
  line_df <- data.frame(y = mean_all, lab = "Full-sample bootstrap mean")

  ylab_text <- if (is.null(ylab_override)) sprintf("LOO bootstrap mean of %s", param) else ylab_override

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$y)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$group, size = .data$group), shape = 16,
      alpha = 0.95, stroke = 0, show.legend = TRUE
    ) +
    ggplot2::geom_hline(
      data = line_df,
      ggplot2::aes(yintercept = .data$y, linetype = .data$lab),
      linewidth = 0.7, show.legend = TRUE, color = "black"
    ) +
    ggplot2::scale_color_manual(
      values = c("Other cases" = "#222222", "Top influential cases" = "#1f77b4"),
      name   = NULL
    ) +
    # 仅用 size 区分点大小，不出单独图例
    ggplot2::scale_size_manual(
      values = c("Other cases" = 1.6, "Top influential cases" = 2.6),
      guide = "none"
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Full-sample bootstrap mean" = "solid"),
      name   = NULL
    ) +
    ggplot2::labs(
      title = "Jackknife-after-Bootstrap diagnostic",
      x = "Observation index",
      y = ylab_text
    ) +
    ggplot2::guides(
      # 颜色图例只显示点
      color = ggplot2::guide_legend(
        order = 1, nrow = 1, byrow = TRUE,
        override.aes = list(linetype = "blank", shape = 16, size = 2.2, alpha = 1),
        keyheight = grid::unit(6, "mm"),
        keywidth  = grid::unit(10, "mm")
      ),
      # 线型图例只显示线
      linetype = ggplot2::guide_legend(
        order = 2, nrow = 1,
        override.aes = list(shape = NA, color = "black", size = 0.7),
        keyheight = grid::unit(6, "mm"),
        keywidth  = grid::unit(10, "mm")
      )
    ) +
    ggplot2::theme_minimal(base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust  = 0.5, face = "bold",
        margin = ggplot2::margin(b = 14)
      ),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),

      # 统一三个图例键尺寸并水平对齐
      legend.position      = "bottom",
      legend.box           = "horizontal",
      legend.margin        = ggplot2::margin(t = 4, b = 6),
      legend.key.height    = grid::unit(6, "mm"),
      legend.key.width     = grid::unit(10, "mm"),
      legend.box.just      = "center",
      legend.justification = "center",
      legend.spacing.x     = grid::unit(6, "pt"),

      plot.margin = ggplot2::margin(t = 10, r = 8, b = 14, l = 8)
    )

  p
}



