# gg_boot_plots.R
# ------------------------------------------------------------
# This file defines ggplot2-based diagnostics for bootstrap
# estimates stored/returned by semboottools.
# ------------------------------------------------------------

#' Diagnostic Plots of Bootstrap Estimates (ggplot2, Modular)
#'
#' @description
#' Produce diagnostic plots of bootstrap estimates using **ggplot2**.
#' Compared to [hist_qq_boot()] (base graphics), this function is *modular* and *modern*:
#' optional layers (histogram, mean line, CI, SD arrow, QQ), and the ability to
#' return only the data or the ggplot object for further customization.
#'
#' @details
#' - For free (unstandardized) parameters, bootstrap draws are extracted directly
#'   from the `lavaan` object;
#' - For user-defined parameters or standardized solutions, [store_boot()] or
#'   [standardizedSolution_boot()] must be called in advance;
#' - Supports objects of class `sbt_std_boot` (standardized) and `sbt_ustd_boot` (unstandardized).
#'
#' **Return mode**
#' - `output = "draw"`: draw the plot (default) and invisibly return the ggplot object;
#' - `output = "ggplot"`: return the ggplot object without drawing;
#' - `output = "data"`: return a list containing the data frames for histogram,
#'   density, and CI positions (class `"sbt_boot_plotdata"`).
#'
#' **Axis control**
#' - By default, the x-axis has a small margin beyond the observed range;
#'   set `trim_quantiles = c(.01, .99)` to crop extreme tails by quantiles.
#'
#' **Optional layers**
#' - `show_hist`, `show_mean`, `show_ci`, `show_sd_arrow`, `show_qq`.
#'
#' @return
#' - If `output = "draw"`: draws the plot and invisibly returns the ggplot object;
#' - If `output = "ggplot"`: returns a ggplot object (no drawing);
#' - If `output = "data"`: returns a list with
#'   \itemize{
#'     \item `t`: bootstrap draws;
#'     \item `t0`: point estimate;
#'     \item `sd`: sample standard deviation;
#'     \item `hist`: data frame for histogram layer;
#'     \item `dens`: data frame for density curve (`x`, `y`);
#'     \item `ci`: list with `lower`, `upper`, `y_lower`, `y_upper`.
#'   }
#'
#' @param object A `lavaan` object with stored bootstrap results, or an
#'   `sbt_std_boot` (from `standardizedSolution_boot()`) or `sbt_ustd_boot`
#'   (from `parameterEstimates_boot()`) object.
#' @param param Character. Name of the parameter to plot (as in `coef()` or `"lhs op rhs"`).
#' @param standardized Logical. Use standardized estimates. Ignored for `sbt_*` inputs.
#' @param bins Integer. Number of histogram bins when `show_hist = TRUE`. Default `35`.
#' @param show_ci Logical. Draw CI vertical dashed lines and labels. Default `TRUE`.
#' @param show_mean Logical. Draw **point estimate** vertical line and label. Default `TRUE`.
#' @param show_sd_arrow Logical. Draw ±SD double‐headed arrow and label. Default `TRUE`.
#' @param show_qq Logical. Add a side-by-side QQ plot (uses `patchwork` if available). Default `TRUE`.
#' @param show_hist Logical. Draw histogram (otherwise density only). Default `TRUE`.
#' @param text_size Numeric. Font size for labels. Default `3`.
#' @param ci_label_digits,mean_label_digits,sd_label_digits Integer digits for labels. Default `3`.
#' @param sd_clip_to_density Logical. Clip SD arrow within density curve span. Default `TRUE`.
#' @param sd_arrow_height Numeric in (0,1]. Relative height for the SD arrow. Default `0.58`.
#' @param trim_quantiles Length-2 numeric. Trim x-axis by quantiles (e.g., `c(.01,.99)`). `NULL` = no trimming.
#' @param ci_line_color Color for CI lines. Default `"#D62728"`.
#' @param dens_line_color Color for density curve. Default `"#8B0000CC"`.
#' @param bar_fill,bar_color Histogram fill/border colors. Defaults `"#5DADE233"` / `"#1B4F72"`.
#' @param enforce_serif Logical. Use serif base family in the theme. Default `TRUE`.
#' @param theme_override Optional `ggplot2::theme()` added on top of base theme.
#' @param dens_adjust Numeric ≥ 0. Bandwidth adjust for `stats::density()`. Default `1`.
#' @param dens_from,dens_to Optional numeric. Force `from`/`to` range for `density()`.
#' @param show_boot_mean Logical. Draw **bootstrap mean** vertical line. Default `TRUE`.
#' @param point_color Color for the **point estimate** line. Default `"#000000"`.      # <-- 新增，名称与代码一致
#' @param boot_mean_color Color for the **bootstrap mean** line. Default `"#AA3377"`.
#' @param output One of `"draw"`, `"ggplot"`, or `"data"`. What to return/draw. Default `"draw"`.
#' @param return Deprecated. Backward-compatible alias of `output`; if supplied, it overrides `output`.  # <-- 新增
#' @param ... Additional arguments passed to ggplot layers (e.g., `geom_histogram()`).
#'
#' @return If `output = "draw"`, draws the plot and (invisibly) returns the ggplot object;
#'   if `output = "ggplot"`, returns the ggplot (no drawing); if `output = "data"`,
#'   returns a list with `t`, `t0`, `sd`, `hist`, `dens`, and `ci` (class `"sbt_boot_plotdata"`).
#'
#' @references
#' Rousselet, G. A., Pernet, C. R., & Wilcox, R. R. (2021).
#' The percentile bootstrap: A primer with step-by-step instructions in R.
#' *Advances in Methods and Practices in Psychological Science*, 4(1), 1–10. \doi{10.1177/2515245920911881}
#'
#' @seealso [store_boot()], [standardizedSolution_boot()], [hist_qq_boot()]
#'
#' @importFrom stats density qqnorm quantile approx sd
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line labs theme_minimal
#' @importFrom ggplot2 theme element_text scale_x_continuous coord_cartesian
#' @importFrom ggplot2 annotate geom_segment geom_vline geom_point geom_smooth
#' @importFrom grid unit arrow
#' @importFrom rlang .data as_label
#' @importFrom patchwork wrap_plots plot_spacer

#' @export

gg_hist_qq_boot <- function(object,
                            param,
                            standardized = NULL,
                            bins = 35,
                            show_ci        = TRUE,
                            show_mean      = TRUE,   # point estimate line
                            show_sd_arrow  = TRUE,
                            show_qq        = TRUE,
                            show_hist      = TRUE,
                            text_size      = 3,
                            ci_label_digits   = 3,
                            mean_label_digits = 3,
                            sd_label_digits   = 3,
                            sd_clip_to_density = TRUE,
                            sd_arrow_height     = 0.58,
                            trim_quantiles = NULL,
                            ci_line_color   = "#D62728",
                            dens_line_color = "#8B0000CC",
                            bar_fill   = "#5DADE233",
                            bar_color  = "#1B4F72",
                            enforce_serif = TRUE,
                            theme_override = NULL,
                            dens_adjust = 1,
                            dens_from = NULL,
                            dens_to   = NULL,
                            ## --- keep only one copy of these:
                            show_boot_mean = TRUE,
                            point_color     = "#000000",   # point estimate line color
                            boot_mean_color = "#AA3377",
                            ## ---
                            output = c("draw", "ggplot", "data"),
                            return = NULL,
                            ...){
  if (is.null(return)) {
    output <- match.arg(output)
  } else {
    output <- match.arg(return, c("draw","ggplot","data"))
  }

  # --- Standardization inference
  if (is.null(standardized) &&
      !(inherits(object, "sbt_std_boot") || inherits(object, "sbt_ustd_boot"))) {
    stop("'standardized' must be TRUE or FALSE.")
  }
  if (inherits(object, "sbt_std_boot") && is.null(standardized))  standardized <- TRUE
  if (inherits(object, "sbt_ustd_boot") && is.null(standardized)) standardized <- FALSE

  # --- Bootstrap draws & point estimate
  bo <- param_find_boot(object = object, param = param, standardized = standardized)
  if (any(is.na(bo))) stop("Bootstrap estimates not found or not stored.")
  t  <- bo$t; t0 <- bo$t0
  if (!is.numeric(t) || length(t) < 2L || diff(range(t)) == 0) {
    stop("Bootstrap draws are degenerate (all identical or non-numeric).")
  }

  # --- Bootstrap mean and bias
  boot_mean <- mean(t, na.rm = TRUE)
  bias_hat  <- boot_mean - t0

  # --- Density & Data frames
  d_obj <- stats::density(
    t,
    adjust = dens_adjust,
    from = if (is.null(dens_from)) min(t) else dens_from,
    to   = if (is.null(dens_to))   max(t) else dens_to
  )
  if (all(!is.finite(d_obj$y)) || all(d_obj$y == 0)) {
    stop("Density is degenerate for the given samples.")
  }
  df_hist <- data.frame(t = t)
  df_dens <- data.frame(x = d_obj$x, y = d_obj$y)
  sdv     <- stats::sd(t, na.rm = TRUE)

  # --- CI retrieval (fallback to percentile)
  ci_list <- .ci_from_object_or_quantile(object, t, param)
  ci_lower <- ci_list$lower
  ci_upper <- ci_list$upper

  # --- Project CI to density heights
  ymax    <- max(df_dens$y)
  y_lower <- if (is.finite(ci_lower)) stats::approx(df_dens$x, df_dens$y, xout = ci_lower)$y else NA_real_
  y_upper <- if (is.finite(ci_upper)) stats::approx(df_dens$x, df_dens$y, xout = ci_upper)$y else NA_real_

  # --- Return data only
  boot_data <- list(
    t = t, t0 = t0, sd = sdv, param = param,
    hist = df_hist, dens = df_dens,
    ci = list(lower = ci_lower, upper = ci_upper, y_lower = y_lower, y_upper = y_upper),
    boot_mean = boot_mean,
    bias_hat  = bias_hat
  )
  class(boot_data) <- c("sbt_boot_plotdata", class(boot_data))
  if (identical(output, "data")) return(boot_data)

  # --- Base theme
  base_theme <- ggplot2::theme_minimal(base_family = if (enforce_serif) "serif" else NULL)

  # --- Main figure
  p <- ggplot2::ggplot(df_hist, ggplot2::aes(.data$t))
  if (isTRUE(show_hist)) {
    p <- p + ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = bar_fill, color = bar_color, linewidth = 0.6, ...
    )
  }
  p <- p +
    ggplot2::geom_line(data = df_dens, ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = 1.2, color = dens_line_color, ...) +
    ggplot2::labs(title = paste0("Histogram of ", param), x = param, y = "Density") +
    base_theme +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.title = ggplot2::element_text(size = 11),
                   axis.text  = ggplot2::element_text(size = 9))
  if (!is.null(theme_override)) p <- p + theme_override

  # --- X axis range and trimming
  x_lim <- range(t, na.rm = TRUE)
  if (!is.null(trim_quantiles)) {
    stopifnot(length(trim_quantiles) == 2, is.numeric(trim_quantiles),
              trim_quantiles[1] >= 0, trim_quantiles[2] <= 1,
              trim_quantiles[1] < trim_quantiles[2])
    q <- stats::quantile(t, probs = trim_quantiles, na.rm = TRUE)
    x_lim <- c(q[1], q[2])
  }
  pad_l_r <- .x_pad(x_lim, left = 0.02, right = 0.06)
  xr_eff  <- diff(x_lim)

  # --- Reserve bottom corridor for CI labels
  y_pad  <- 0.10 * ymax
  y_lab  <- -0.40 * y_pad

  p <- p +
    ggplot2::scale_x_continuous(limits = pad_l_r, expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_y_continuous(limits = c(-y_pad, NA_real_), expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::coord_cartesian(clip = "off")

  # --- CI vertical dashed lines + labels
  if (isTRUE(show_ci) && is.finite(ci_lower)) {
    p <- p +
      ggplot2::geom_segment(x = ci_lower, xend = ci_lower, y = 0, yend = y_lower,
                            linetype = "dashed", color = ci_line_color, linewidth = 1.0) +
      ggplot2::annotate("text",
                        x = ci_lower, y = y_lab,
                        label = paste0("L = ", round(ci_lower, ci_label_digits)),
                        size = text_size, hjust = 0.5, vjust = 1, fontface = 2, color = "black")
  }
  if (isTRUE(show_ci) && is.finite(ci_upper)) {
    p <- p +
      ggplot2::geom_segment(x = ci_upper, xend = ci_upper, y = 0, yend = y_upper,
                            linetype = "dashed", color = ci_line_color, linewidth = 1.0) +
      ggplot2::annotate("text",
                        x = ci_upper, y = y_lab,
                        label = paste0("U = ", round(ci_upper, ci_label_digits)),
                        size = text_size, hjust = 0.5, vjust = 1, fontface = 2, color = "black")
  }

  # --- Optional vertical lines (point estimate & bootstrap mean) + legend (bottom)
  lab_point <- sprintf("Point estimate = %.*f", mean_label_digits, t0)
  lab_boot  <- sprintf("Bootstrap mean = %.*f", mean_label_digits, boot_mean)

  labels <- character(0); xs <- numeric(0); cols <- character(0); ltys <- character(0)
  if (isTRUE(show_mean)) {
    labels <- c(labels, lab_point); xs <- c(xs, t0); cols <- c(cols, point_color); ltys <- c(ltys, "dashed")
  }
  if (isTRUE(show_boot_mean)) {
    labels <- c(labels, lab_boot); xs <- c(xs, boot_mean); cols <- c(cols, boot_mean_color); ltys <- c(ltys, "dashed")
  }

  if (length(labels) > 0) {
    lines_df <- data.frame(
      x    = xs,
      type = factor(labels, levels = labels),
      stringsAsFactors = FALSE
    )
    p <- p +
      ggplot2::geom_vline(
        data = lines_df,
        ggplot2::aes(xintercept = .data$x, colour = .data$type, linetype = .data$type),
        linewidth = 0.9, show.legend = TRUE
      ) +
      ggplot2::scale_colour_manual(values = stats::setNames(cols, labels), breaks = labels, name = NULL) +
      ggplot2::scale_linetype_manual(values = stats::setNames(ltys, labels), breaks = labels, name = NULL) +
      # 只保留 colour 图例的 override，关闭 linetype 图例以避免 duplicated override.aes
      ggplot2::guides(
        colour   = ggplot2::guide_legend(order = 1, nrow = 1, byrow = TRUE,
                                         override.aes = list(linewidth = 1.3)),
        linetype = "none"
      ) +
      ggplot2::theme(
        legend.position   = "bottom",
        legend.box        = "horizontal",
        legend.margin = ggplot2::margin(t = -2, r = 2, b = 2, l = 2),
        legend.key.width  = grid::unit(10, "pt"),
        legend.key.height = grid::unit(6, "pt"),
        legend.spacing.x  = grid::unit(4, "pt")
      )
  }

  # --- SD double-headed arrow
  if (isTRUE(show_sd_arrow)) {
    arrow_y <- sd_arrow_height * ymax
    if (isTRUE(sd_clip_to_density)) {
      span <- .dens_span_at(arrow_y, df_dens$x, df_dens$y)
      if (anyNA(span)) {
        for (k in seq(sd_arrow_height, 0.20, by = -0.05)) {
          arrow_y <- k * ymax
          span <- .dens_span_at(arrow_y, df_dens$x, df_dens$y)
          if (!anyNA(span)) break
        }
      }
      x_left  <- max(t0 - sdv, span[1], na.rm = TRUE)
      x_right <- min(t0 + sdv, span[2], na.rm = TRUE)
    } else {
      x_left  <- t0 - sdv
      x_right <- t0 + sdv
    }
    if (is.finite(x_left) && is.finite(x_right) && (x_left < x_right)) {
      p <- p +
        ggplot2::annotate("segment",
                          x = x_left, xend = x_right, y = arrow_y, yend = arrow_y,
                          linetype = "dashed", color = "grey40", linewidth = 0.9,
                          arrow = grid::arrow(ends = "both", type = "closed", length = grid::unit(5, "pt"))) +
        ggplot2::annotate("text",
                          x = x_right + 0.03 * xr_eff, y = arrow_y + 0.02 * ymax,
                          label = paste0("SD = ", round(sdv, sd_label_digits)),
                          size = text_size, hjust = 0, fontface = 2, color = "black")
    }
  }

  # --- Output or combine with QQ
  if (identical(output, "ggplot") && !isTRUE(show_qq)) return(p)
  if (!isTRUE(show_qq)) { print(p); return(invisible(p)) }

  # QQ (side-by-side)
  q <- stats::qqnorm(t, plot.it = FALSE)
  dfqq <- data.frame(theoretical = q$x, sample = q$y)
  pqq <- ggplot2::ggplot(dfqq, ggplot2::aes(.data$theoretical, .data$sample)) +
    ggplot2::geom_point(shape = 21, size = 2, color = "#1B4F72", fill = bar_fill) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = dens_line_color, formula = y ~ x) +
    ggplot2::labs(title = paste0("Normal QQ-Plot of ", param),
                  x = "Quantiles of Standard Normal", y = param) +
    base_theme +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.title = ggplot2::element_text(size = 11),
                   axis.text  = ggplot2::element_text(size = 9))
  if (!is.null(theme_override)) pqq <- pqq + theme_override

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combo <- patchwork::wrap_plots(p, patchwork::plot_spacer(), pqq,
                                   nrow = 1, widths = c(1, 0.05, 1))
    if (identical(output, "ggplot")) return(combo)
    print(combo); return(invisible(combo))
  } else {
    print(p); print(pqq); return(invisible(p))
  }
}


#' Scatterplot Matrix of Bootstrap Estimates (ggplot2/GGally)
#'
#' @description
#' Generates a scatterplot matrix of the bootstrap estimates for
#' two or more model parameters using `ggplot2` and `GGally`,
#' with optional smoothing lines, confidence ellipses, and
#' customizable density or histogram plots on the diagonal.
#'
#' @details
#' The function [gg_scatter_boot()] is an enhanced ggplot2-based
#' version of [scatter_boot()], designed to produce publication-quality
#' scatterplot matrices of bootstrap estimates.
#' It can be applied to the output of [standardizedSolution_boot()]
#' or [parameterEstimates_boot()].
#'
#' @param object A bootstrap result object of class
#'   `"sbt_std_boot"` or `"sbt_ustd_boot"`.
#' @param params Character vector of parameter names to plot (≥ 2).
#' @param standardized Logical; whether to use standardized estimates.
#' @param title Plot title.
#' @param point_size,point_alpha,point_color Aesthetics for points in lower panels.
#' @param show_smooth,smooth_method,smooth_se Control the regression smoother in lower panels.
#' @param show_ellipse,ellipse_level,ellipse_color Confidence ellipse options in lower panels.
#' @param diag_type One of `"both"`, `"density"`, `"hist"`, `"blank"` for diagonal panels.
#' @param hist_border_size Histogram outline width in diagonal `"hist"`/`"both"`.
#' @param bins,dens_fill,dens_color,hist_fill,hist_color Diagonal panel styles.
#' @param show_corr Logical; if TRUE, show correlation coefficients in upper panels.
#' @param corr_text_size,corr_digits Text size and digits for upper-panel correlation.
#' @param panel_border_color,panel_border_size Panel border color and width for all panels.
#' @param show_mean_diag,mean_line_color,mean_linewidth,mean_linetype
#'   Draw a vertical mean line in diagonal panels.
#' @param output Either `"draw"` (print) or `"ggplot"` (return ggplot object).
#'
#' @return Invisibly returns the `ggplot` object (or prints it if `output="draw"`).
#'
#' @seealso [scatter_boot()], [hist_qq_boot()],
#'   [standardizedSolution_boot()], [parameterEstimates_boot()]
#'
#' @export

gg_scatter_boot <- function(object,
                            params,
                            standardized = NULL,
                            title = "Bootstrap Estimates",
                            # point
                            point_size = 1.8,
                            point_alpha = 0.35,
                            point_color = "#5DADE233",
                            # smoother & ellipse
                            show_smooth  = TRUE,
                            smooth_method = "lm",
                            smooth_se     = FALSE,
                            show_ellipse  = TRUE,
                            ellipse_level = 0.95,
                            ellipse_color = "#8B0000CC",
                            # diagonal
                            diag_type = c("both", "density", "hist", "blank"),
                            hist_border_size = 0.2,
                            bins = 30,
                            dens_fill = "#5DADE233",
                            dens_color = "#8B0000CC",
                            hist_fill = "#5DADE233",
                            hist_color = "#1B4F72",
                            show_corr = TRUE,
                            corr_text_size = 7,
                            panel_border_color = "grey40",
                            panel_border_size  = 2,
                            corr_digits = 2,
                            show_mean_diag   = TRUE,
                            mean_line_color  = "black",
                            mean_linewidth   = 0.7,
                            mean_linetype    = "dashed",
                            output = c("draw", "ggplot")) {

  output   <- match.arg(output)
  diag_type <- match.arg(diag_type)

  if (!requireNamespace("GGally", quietly = TRUE))
    stop("Please install.packages('GGally') to use gg_scatter_boot().")

  # --- Standardization inference
  if (is.null(standardized) &&
      !(inherits(object, "sbt_std_boot") || inherits(object, "sbt_ustd_boot"))) {
    stop("'standardized' must be TRUE or FALSE.")
  }
  if (inherits(object, "sbt_std_boot") && is.null(standardized))  standardized <- TRUE
  if (inherits(object, "sbt_ustd_boot") && is.null(standardized)) standardized <- FALSE
  if (length(params) < 2) stop("Need to select two or more parameters.")

  # --- Collect bootstrap vectors
  boot_out_list <- lapply(params, function(p) {
    out <- param_find_boot(object = object, param = p, standardized = standardized)
    if (any(is.na(out))) stop("Bootstrap estimates not found or not stored for parameter: ", p)
    out$t
  })
  names(boot_out_list) <- params
  df <- as.data.frame(boot_out_list, check.names = FALSE)

  # --- Panel functions (no library(), use ggplot2:: fully qualified)
  lower_fun <- function(data, mapping, ...) {
    p <- ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha, color = point_color, ...)
    if (isTRUE(show_smooth)) {
      # 显式公式并静默默认提示
      p <- suppressMessages(
        p + ggplot2::geom_smooth(method = smooth_method, se = smooth_se,
                                 formula = y ~ x, linewidth = 0.8, color = dens_color, ...)
      )
    }
    if (isTRUE(show_ellipse)) {
      p <- p + ggplot2::stat_ellipse(level = ellipse_level, linewidth = 0.8,
                                     color = ellipse_color, linetype = "dashed", ...)
    }
    p +
      ggplot2::theme_minimal(base_family = "serif") +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA,
                                             linewidth = panel_border_size)
      )
  }

  upper_fun <- if (isTRUE(show_corr)) {
    function(data, mapping, ...) {
      x <- ggplot2::ggplot_build(ggplot2::ggplot(data, mapping) + ggplot2::geom_point())
      df_num <- ggplot2::layer_data(x$plot, 1)
      r <- suppressWarnings(stats::cor(df_num$x, df_num$y, use = "pairwise.complete.obs"))
      lbl <- if (is.finite(r)) sprintf("r = %.2f", round(r, corr_digits)) else "r = NA"
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = lbl,
                          size = corr_text_size, fontface = 2) +
        ggplot2::theme_minimal(base_family = "serif") +
        ggplot2::theme(
          axis.title = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA,
                                               linewidth = panel_border_size)
        )
    }
  } else {
    lower_fun
  }

  diag_fun <- switch(
    diag_type,
    "density" = function(data, mapping, ...) {
      x_var <- rlang::as_label(mapping$x)
      mu <- mean(data[[x_var]], na.rm = TRUE)
      ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_density(fill = dens_fill, color = dens_color, linewidth = 1.0, ...) +
        { if (isTRUE(show_mean_diag))
          ggplot2::geom_vline(xintercept = mu, color = mean_line_color,
                              linewidth = mean_linewidth, linetype = mean_linetype) } +
        ggplot2::theme_minimal(base_family = "serif") +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA,
                                               linewidth = panel_border_size)
        )
    },
    "hist" = function(data, mapping, ...) {
      x_var <- rlang::as_label(mapping$x)
      mu <- mean(data[[x_var]], na.rm = TRUE)
      ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = bins, fill = hist_fill,
                                colour = hist_color, size = hist_border_size, ...) +
        { if (isTRUE(show_mean_diag))
          ggplot2::geom_vline(xintercept = mu, color = mean_line_color,
                              linewidth = mean_linewidth, linetype = mean_linetype) } +
        ggplot2::theme_minimal(base_family = "serif") +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA,
                                               linewidth = panel_border_size)
        )
    },
    "both" = function(data, mapping, ...) {
      x_var <- rlang::as_label(mapping$x)
      mu <- mean(data[[x_var]], na.rm = TRUE)
      ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                                bins = bins, fill = hist_fill,
                                colour = hist_color, size = hist_border_size, ...) +
        ggplot2::geom_density(color = dens_color, linewidth = 1.0, ...) +
        { if (isTRUE(show_mean_diag))
          ggplot2::geom_vline(xintercept = mu, color = mean_line_color,
                              linewidth = mean_linewidth, linetype = mean_linetype) } +
        ggplot2::theme_minimal(base_family = "serif") +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA,
                                               linewidth = panel_border_size)
        )
    },
    "blank" = function(data, mapping, ...) {
      ggplot2::ggplot() +
        ggplot2::theme_void(base_family = "serif") +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA,
                                               linewidth = panel_border_size)
        )
    }
  )

  gp <- GGally::ggpairs(
    df,
    lower = list(continuous = lower_fun),
    upper = list(continuous = upper_fun),
    diag  = list(continuous = diag_fun),
    columnLabels = params,
    progress = FALSE,
    labeller = "label_parsed"
  ) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = title)

  if (identical(output, "ggplot")) return(gp)
  print(gp)
  invisible(gp)
}

# ---------- Internal helpers (not exported) ----------

# Robust CI extraction from sbt_* objects; fallback to percentile of t
#' @noRd
.ci_from_object_or_quantile <- function(object, t, param) {
  lower <- upper <- NA_real_
  has_ci_cols <- (inherits(object, "sbt_std_boot") || inherits(object, "sbt_ustd_boot")) &&
    all(c("boot.ci.lower","boot.ci.upper") %in% names(object))
  if (isTRUE(has_ci_cols)) {
    key   <- gsub("\\s+", "", param)
    rowid <- gsub("\\s+", "", paste0(object$lhs, object$op, object$rhs))
    hits  <- which(rowid == key | object$label == param)
    if (length(hits) > 1 && "group" %in% names(object)) {
      stop("Multiple matches for '", param, "'. Please disambiguate by group.")
    }
    if (length(hits) > 1) {
      stop("Multiple matches for '", param, "'. Please use an unambiguous name/label.")
    }
    if (length(hits) == 1) {
      lower <- object$boot.ci.lower[hits]
      upper <- object$boot.ci.upper[hits]
    }
  }
  if (!is.finite(lower) || !is.finite(upper)) {
    level <- attr(object, "level"); if (is.null(level) || !is.numeric(level)) level <- 0.95
    alpha <- (1 - level) / 2
    qs <- stats::quantile(t, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
    lower <- qs[1]; upper <- qs[2]
  }
  list(lower = lower, upper = upper)
}

# Find the horizontal span [xL, xR] of density curve at a given height y
#' @noRd
.dens_span_at <- function(y, x, f) {
  idx <- which(f >= y)
  if (length(idx) == 0L) return(c(NA_real_, NA_real_))
  c(min(x[idx]), max(x[idx]))
}

# Build padded x-limits around a given range
#' @noRd
.x_pad <- function(x_lim, left = 0.02, right = 0.06) {
  dx <- diff(x_lim)
  c(x_lim[1] - left * dx, x_lim[2] + right * dx)
}


# -------- param_find_boot (kept as in your package, minor tidy) --------

# Find and return the stored bootstrap estimates, if available.
#' @noRd
param_find_boot <- function(object,
                            param,
                            standardized) {
  boot_t0 <- NA
  boot_t <- NA
  is_sbt_std_boot <- inherits(object, "sbt_std_boot")
  if (is_sbt_std_boot) standardized <- TRUE
  is_sbt_ustd_boot <- inherits(object, "sbt_ustd_boot")
  if (is_sbt_ustd_boot) standardized <- FALSE

  if (standardized) {
    if (is_sbt_std_boot) {
      coef_names <- lavaan::lav_partable_labels(object)
      boot_i <- attr(object, "boot_est_std")
      if (is.null(boot_i)) stop("Bootstrap estimates not found in the object.")
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
        boot_t0 <- stats::setNames(lavaan::standardizedSolution(object, se = FALSE)[i, "est.std"], param)
      }
      boot_t <- boot_i[, param, drop = TRUE]
      out <- list(t0 = boot_t0, t = boot_t); return(out)
    }
  } else {
    if (is_sbt_ustd_boot) {
      coef_names <- lavaan::lav_partable_labels(object)
      boot_i <- attr(object, "boot_est_ustd")
      if (is.null(boot_i)) stop("Bootstrap estimates not found in the object.")
    } else {
      boot_i0 <- try(lavaan::lavInspect(object, "boot"), silent = TRUE)
      boot_i1 <- object@external$sbt_boot_ustd
      if (inherits(boot_i0, "try-error") && is.null(boot_i1)) {
        stop("Bootstrapping estimates not found.")
      }
      boot_i <- if (!inherits(boot_i0, "try-error")) boot_i0 else boot_i1
    }
    error_idx <- attr(boot_i, "error.idx")
    if (length(error_idx) != 0) boot_i <- boot_i[-error_idx, , drop = FALSE]
    if (param %in% colnames(boot_i)) {
      if (is_sbt_ustd_boot) {
        i <- match(param, coef_names)
        boot_t0 <- object[i, "est"]
      } else {
        boot_t0 <- lavaan::coef(object)[param]
      }
      boot_t <- boot_i[, param, drop = TRUE]
      out <- list(t0 = boot_t0, t = boot_t); return(out)
    } else {
      boot_i <- get_boot_def(object)
      if (param %in% colnames(boot_i)) {
        if (is_sbt_ustd_boot) {
          i <- match(param, object[, "label"])
          boot_t0 <- object[i, "est"]
        } else {
          boot_t0 <- lavaan::coef(object, type = "user")[param]
        }
        boot_t <- boot_i[, param, drop = TRUE]
        out <- list(t0 = boot_t0, t = boot_t); return(out)
      }
    }
  }
  list(t0 = boot_t0, t = boot_t)
}
