#' @title Diagnostic Plots of Bootstrap Estimates (ggplot2, Modular)
#'
#' @description
#' Produce diagnostic plots of bootstrap estimates using **ggplot2**.
#' Compared to [hist_qq_boot()] (base R graphics), this function provides
#' a more *modular* and *modern* syntax:
#' - Optional layers (histogram, mean line, CI lines, SD arrows, QQ plot, etc.)
#' - Ability to return only the data or the ggplot object for further customization
#'
#' @details
#' - For free (unstandardized) parameters, bootstrap draws are extracted directly
#'   from the `lavaan` object;
#' - For user-defined parameters or standardized solutions, [store_boot()] or
#'   [standardizedSolution_boot()] must be called in advance;
#' - Supports stored bootstrap objects of class `sbt_std_boot` (standardized)
#'   and `sbt_ustd_boot` (unstandardized).
#'
#' **Modular return:**
#' - `return = "draw"`: draw the plot (default) and return the ggplot object invisibly;
#' - `return = "ggplot"`: return the ggplot object without drawing; users can
#'   add external themes via `+ theme(...)`;
#' - `return = "data"`: return a list containing the data frames for histogram,
#'   density curve, and CI positions.
#'
#' **Axis control:**
#' - By default, the x-axis has a small margin beyond the observed range;
#'   use `trim_quantiles = c(.01, .99)` to crop extreme tails by quantiles.
#'
#' **Optional layers:**
#' - `show_hist`: histogram; `show_mean`: mean line and label;
#' - `show_ci`: percentile bootstrap CI short dashed lines and labels ("L=", "U=");
#' - `show_sd_arrow`: ±SD double-headed dashed arrows and labels;
#' - `show_qq`: side-by-side normal QQ plot (requires `patchwork` for layout).
#'
#' @return
#' - If `return = "draw"`: draws the plot and invisibly returns the ggplot object;
#' - If `return = "ggplot"`: returns a ggplot object (no drawing);
#' - If `return = "data"`: returns a list with:
#'   \itemize{
#'     \item `t`: bootstrap draws;
#'     \item `t0`: point estimate;
#'     \item `sd`: sample standard deviation;
#'     \item `hist`: data frame for histogram layer;
#'     \item `dens`: data frame for density curve (`x`, `y`);
#'     \item `ci`: list with `lower`, `upper`, `y_lower`, `y_upper`.
#'   }
#'
#' @param object
#' A `lavaan` object with stored bootstrap results,
#' an object of class `sbt_std_boot` (from [standardizedSolution_boot()]),
#' or `sbt_ustd_boot`.
#' @param param Character. Name of the parameter to plot
#'   (as in `coef()` or `"lhs op rhs"` form).
#' @param standardized Logical. Whether to use standardized bootstrap estimates.
#' Ignored for `sbt_std_boot` / `sbt_ustd_boot`; required for `lavaan` objects.
#' @param bins Integer. Number of histogram bins (if `show_hist = TRUE`). Default = 35.
#' @param show_ci Logical. Draw CI short dashed lines and labels. Default = TRUE.
#' @param show_mean Logical. Draw mean vertical dashed line and label. Default = TRUE.
#' @param show_sd_arrow Logical. Draw ±SD arrows and labels. Default = TRUE.
#' @param show_qq Logical. Add side-by-side QQ plot. Default = TRUE.
#' @param show_hist Logical. Draw histogram (otherwise density only). Default = FALSE.
#' @param text_size Numeric. Font size for labels. Default = 3.
#' @param ci_label_digits Integer. Digits for CI labels. Default = 3.
#' @param mean_label_digits Integer. Digits for mean label. Default = 3.
#' @param sd_label_digits Integer. Digits for SD labels. Default = 3.
#' @param trim_quantiles Length-2 numeric vector. Quantile range to trim x-axis
#' (e.g., `c(.005, .995)`); NULL = no trimming. Default = NULL.
#' @param ci_line_color Color for CI short dashed lines. Default = `"black"`.
#' @param dens_line_color Color for density curve. Default = `"#8B0000CC"`.
#' @param bar_fill Fill color for histogram (if `show_hist = TRUE`). Default = `"#5DADE233"`.
#' @param bar_color Border color for histogram. Default = `"#1B4F72"`.
#' @param return Character: `"draw"`, `"ggplot"`, or `"data"`. See Value.
#'
#' @references
#' Rousselet, G. A., Pernet, C. R., & Wilcox, R. R. (2021).
#' The percentile bootstrap: A primer with step-by-step instructions in R.
#' *Advances in Methods and Practices in Psychological Science*, 4(1), 1–10.
#' \doi{10.1177/2515245920911881}
#'
#' @seealso
#' [store_boot()], [standardizedSolution_boot()], [hist_qq_boot()]
#'
#' @importFrom stats density qqnorm quantile approx sd
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line labs theme_minimal
#' @importFrom ggplot2 theme element_text scale_x_continuous coord_cartesian
#' @importFrom ggplot2 annotate geom_segment geom_vline geom_point geom_smooth
#' @export

gg_hist_qq_boot <- function(object,
                            param,
                            standardized = NULL,
                            bins = 35,
                            show_ci        = TRUE,
                            show_mean      = TRUE,
                            show_sd_arrow  = TRUE,
                            show_qq        = FALSE,
                            show_hist = FALSE,
                            text_size      = 3,
                            ci_label_digits   = 3,
                            mean_label_digits = 3,
                            sd_label_digits   = 3,
                            trim_quantiles = NULL,          # e.g., c(.005,.995) 可裁尾
                            ci_line_color  = "black",     # 红色 CI 线
                            dens_line_color = "#8B0000CC",
                            bar_fill   = "#5DADE233",
                            bar_color  = "#1B4F72",
                            return = c("draw", "ggplot", "data")) {
  return <- match.arg(return)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  library(ggplot2)

  # --- standardized 推断
  if (is.null(standardized) &&
      !(inherits(object, "sbt_std_boot") || inherits(object, "sbt_ustd_boot"))) {
    stop("'standardized' must be TRUE or FALSE.")
  }
  if (inherits(object, "sbt_std_boot") && is.null(standardized))  standardized <- TRUE
  if (inherits(object, "sbt_ustd_boot") && is.null(standardized)) standardized <- FALSE

  # --- 取 bootstrap 样本与点估计
  bo <- param_find_boot(object = object, param = param, standardized = standardized)
  if (any(is.na(bo))) stop("Bootstrap estimates not found or not stored.")
  t  <- bo$t; t0 <- bo$t0
  if (diff(range(t)) == 0) stop("Identical estimates in all bootstrap samples.")

  # --- 数据
  d_obj   <- stats::density(t)
  df_hist <- data.frame(t = t)
  df_dens <- data.frame(x = d_obj$x, y = d_obj$y)
  sdv     <- stats::sd(t, na.rm = TRUE)
  x_min   <- min(t); x_max <- max(t)

  ## ---- obtain Bootstrap CI (robust matching) ----
  ci_lower <- ci_upper <- NA_real_

  has_ci_cols <- (inherits(object, "sbt_std_boot") || inherits(object, "sbt_ustd_boot")) &&
    all(c("boot.ci.lower","boot.ci.upper") %in% names(object))

  if (has_ci_cols) {
    key   <- gsub("\\s+", "", param)  # 去掉所有空白
    rowid <- gsub("\\s+", "", paste0(object$lhs, object$op, object$rhs))
    hits  <- which(rowid == key | object$label == param)  # 允许 label 命中

    # 如有多组，可以在这里再加一行 group 过滤，例如：
    # if ("group" %in% names(object)) hits <- hits[object$group[hits] == 1]

    if (length(hits) == 1) {
      ci_lower <- object$boot.ci.lower[hits]
      ci_upper <- object$boot.ci.upper[hits]
    }
  }

  # t: 你的 bootstrap 向量；df_dens: density 数据框
  if (is.na(ci_lower) || is.na(ci_upper)) {
    level <- attr(object, "level"); if (is.null(level) || !is.numeric(level)) level <- 0.95
    alpha <- (1 - level) / 2
    qs <- stats::quantile(t, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
    ci_lower <- qs[1]; ci_upper <- qs[2]
  }

  # 把 CI 投到密度曲线高度，供 segments/text 使用
  y_lower <- stats::approx(df_dens$x, df_dens$y, xout = ci_lower)$y
  y_upper <- stats::approx(df_dens$x, df_dens$y, xout = ci_upper)$y


  # 密度高度 & 尺度
  ymax <- max(df_dens$y)
  y_lower <- if (!is.na(ci_lower)) stats::approx(df_dens$x, df_dens$y, xout = ci_lower)$y else NA_real_
  y_upper <- if (!is.na(ci_upper)) stats::approx(df_dens$x, df_dens$y, xout = ci_upper)$y else NA_real_

  # 可只返回数据
  boot_data <- list(
    t = t, t0 = t0, sd = sdv, param = param,
    hist = df_hist, dens = df_dens,
    ci = list(lower = ci_lower, upper = ci_upper, y_lower = y_lower, y_upper = y_upper)
  )
  if (identical(return, "data")) return(boot_data)

  # --- 主图
  p <- ggplot(df_hist, aes(t))
  if (show_hist) {
  p <- p + geom_histogram(aes(y = after_stat(density)),
                          bins = bins, fill = bar_fill,
                          color = bar_color, linewidth = 1.1)}


  p <- p +
  geom_line(data = df_dens, aes(x = x, y = y),
            linewidth = 1.2, color = dens_line_color) +
  labs(title = paste0("Histogram of ", param), x = param, y = "Density") +
  theme_minimal(base_family = "serif") +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(size = 11),
        axis.text  = element_text(size = 9))

  # --- 横轴范围：仅少量边距；可选按分位数裁切
  x_lim <- range(t, na.rm = TRUE)
  if (!is.null(trim_quantiles)) {
    stopifnot(length(trim_quantiles) == 2, trim_quantiles[1] < trim_quantiles[2])
    q <- stats::quantile(t, probs = trim_quantiles, na.rm = TRUE)
    x_lim <- c(q[1], q[2])
  }
  xr_eff <- diff(x_lim)
  p <- p +
    coord_cartesian(xlim = c(x_lim[1] - 0.02 * xr_eff, x_lim[2] + 0.06 * xr_eff)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)))

  # --- CI 短虚线（红）+ 文本（靠近曲线，不压柱）
  if (show_ci && !is.na(ci_lower)) {
    p <- p +
      geom_segment(x = ci_lower, xend = ci_lower, y = 0, yend = y_lower,
                   linetype = "dashed", color = ci_line_color, linewidth = 1.0) +
      annotate("text",
               x = ci_lower - 0.003 * xr_eff,           # 向左微移
               y = y_lower + 0.10 * ymax,               # 曲线下方一点
               label = paste0("L = ", round(ci_lower, ci_label_digits)),
               size = text_size, hjust = 1, fontface = 2, color = "black")
  }
  if (show_ci && !is.na(ci_upper)) {
    p <- p +
      geom_segment(x = ci_upper, xend = ci_upper, y = 0, yend = y_upper,
                   linetype = "dashed", color = ci_line_color, linewidth = 1.0) +
      annotate("text",
               x = ci_upper + 0.015 * xr_eff,           # 向右微移
               y = y_upper + 0.05 * ymax,               # 曲线上方一点
               label = paste0("U = ", round(ci_upper, ci_label_digits)),
               size = text_size, hjust = 0, fontface = 2, color = "black")
  }

  # --- 均值虚线 + 文本（贴线但不重叠）
  if (show_mean) {
    p <- p +
      geom_vline(xintercept = t0, linetype = "dashed", color = "black", linewidth = 0.9) +
      annotate("text",
               x = t0 + 0.01 * xr_eff,
               y = 1.05 * ymax,
               label = paste0("Sample mean = ", round(t0, mean_label_digits)),
               size = text_size, hjust = 0, fontface = 2, color = "black")
  }

  # --- SD 双向虚线箭头 + 文本（曲线中上部）
  if (show_sd_arrow) {
    arrow_y <- 0.58 * ymax
    p <- p +
      annotate("segment",
               x = t0 - sdv, xend = t0 + sdv, y = arrow_y, yend = arrow_y,
               linetype = "dashed", color = "grey40", linewidth = 0.9,
               arrow = arrow(ends = "both", type = "closed", length = unit(5, "pt"))) +
      annotate("text",
               x = t0 + sdv + 0.015 * xr_eff,
               y = arrow_y + 0.02 * ymax,
               label = paste0("SD = ", round(sdv, sd_label_digits)),
               size = text_size, hjust = 0, fontface = 2, color = "black")
  }

  if (identical(return, "ggplot")) return(p)

  # --- 可选：QQ 图并排
  if (isTRUE(show_qq)) {
    q <- stats::qqnorm(t, plot.it = FALSE)
    dfqq <- data.frame(theoretical = q$x, sample = q$y)
    pqq <- ggplot(dfqq, aes(theoretical, sample)) +
      geom_point(shape = 21, size = 2, color = "#1B4F72", fill = bar_fill) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = dens_line_color) +
      labs(title = paste0("Normal QQ-Plot of ", param),
           x = "Quantiles of Standard Normal", y = param) +
      theme_minimal(base_family = "serif") +
      theme(plot.title = element_text(face = "bold"),
            axis.title = element_text(size = 11),
            axis.text  = element_text(size = 9))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      print(p + patchwork::plot_spacer() + pqq + patchwork::plot_layout(widths = c(1, .05, 1)))
    } else { print(p); print(pqq) }
  } else {
    print(p)
  }
  invisible(p)
}


#' Scatterplot Matrix of Bootstrap Estimates (ggplot2 Version)
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
#' Compared to [scatter_boot()], this function allows for:
#' \itemize{
#'   \item Custom point size, transparency, and color.
#'   \item Optional regression smoothing lines.
#'   \item Optional confidence ellipses.
#'   \item Choice between density plots, histograms, or blank panels on the diagonal.
#'   \item Customizable panel borders and correlation text size.
#' }
#'
#' @param object A bootstrap result object of class
#'   `"sbt_std_boot"` or `"sbt_ustd_boot"`, typically returned by
#'   [standardizedSolution_boot()] or [parameterEstimates_boot()].
#' @param params Character vector of parameter names to plot.
#'   At least two parameters must be provided. Names must match those
#'   returned by `coef()` in the fitted model.
#' @param standardized Logical; whether to use standardized estimates.
#'   If `NULL`, the function will infer this from the class of `object`.
#' @param title Title of the scatterplot matrix.
#'   Default is `"Bootstrap Estimates"`.
#' @param point_size Numeric; size of scatterplot points.
#' @param point_alpha Numeric (0–1); transparency of scatterplot points.
#' @param point_color Color for scatterplot points.
#' @param show_smooth Logical; whether to add a regression smoothing line
#'   to lower panels.
#' @param smooth_method Smoothing method passed to [geom_smooth()].
#' @param smooth_se Logical; whether to show the standard error band
#'   for the smoothing line.
#' @param show_ellipse Logical; whether to add a confidence ellipse
#'   to lower panels.
#' @param ellipse_level Confidence level for the ellipse.
#' @param ellipse_color Color for the ellipse.
#' @param diag_type Character; type of plot to show on the diagonal.
#'   One of `"density"`, `"hist"`, or `"blank"`.
#' @param bins Number of bins for histograms (if `diag_type = "hist"`).
#' @param dens_fill Fill color for density plots.
#' @param dens_color Outline color for density plots.
#' @param hist_fill Fill color for histograms.
#' @param hist_color Outline color for histograms.
#' @param show_corr Logical; whether to display correlation coefficients
#'   in the upper panels.
#' @param corr_text_size Numeric; text size for correlation coefficients.
#' @param panel_border_color Color for the panel borders.
#' @param panel_border_size Numeric; line width for panel borders.
#' @param corr_digits Number of decimal places for correlation coefficients.
#' @param return Either `"draw"` (print the plot) or `"ggplot"` (return
#'   the `ggplot` object without printing).
#'
#' @return Invisibly returns the `ggplot` object (or prints it if
#'   `return = "draw"`).
#'
#' @examples
#' \dontrun{
#' # Fit a model and store bootstrap estimates
#' fit <- sem(model_med, data = mydata)
#' fit <- store_boot(fit, R = 1000)
#'
#' # Standardized solution
#' std_boot <- standardizedSolution_boot(fit)
#'
#' # ggplot2-based scatterplot matrix
#' gg_scatter_boot(std_boot, c("a", "b", "ab"), standardized = TRUE)
#' }
#'
#' @seealso [scatter_boot()], [hist_qq_boot()],
#'   [standardizedSolution_boot()], [parameterEstimates_boot()]
#'
#' @export

gg_scatter_boot <- function(object,
                            params,
                            standardized = NULL,
                            title = "Bootstrap Estimates",
                            # 点的样式
                            point_size = 1.8,
                            point_alpha = 0.35,
                            point_color = "#5DADE233",
                            # 平滑与置信椭圆
                            show_smooth  = TRUE,
                            smooth_method = "lm",
                            smooth_se     = FALSE,
                            show_ellipse  = TRUE,
                            ellipse_level = 0.95,
                            ellipse_color = "#8B0000CC",
                            # 对角线：density 或 hist
                            diag_type = c("density", "hist", "blank"),
                            bins = 30,
                            dens_fill = "#5DADE233",
                            dens_color = "#8B0000CC",
                            hist_fill = "#5DADE233",
                            hist_color = "#1B4F72",
                            show_corr = TRUE,
                            corr_text_size = 6,                      # 放大 r 标签
                            panel_border_color = "grey40",           # 面板边框颜色
                            panel_border_size  = 2,                # 面板边框粗细
                            corr_digits = 2,
                            # 返回
                            return = c("draw", "ggplot")) {
  return <- match.arg(return)
  diag_type <- match.arg(diag_type)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install.packages('ggplot2')")
  if (!requireNamespace("GGally", quietly = TRUE))
    stop("Please install.packages('GGally')")

  # --- standardized 推断（与原逻辑一致）
  if (is.null(standardized) &&
      !(inherits(object, "sbt_std_boot") || inherits(object, "sbt_ustd_boot"))) {
    stop("'standardized' must be TRUE or FALSE.")
  }
  if (inherits(object, "sbt_std_boot") && is.null(standardized))  standardized <- TRUE
  if (inherits(object, "sbt_ustd_boot") && is.null(standardized)) standardized <- FALSE
  if (length(params) < 2) stop("Need to select two or more parameters.")

  # --- 收集每个参数的 bootstrap 向量
  boot_out_list <- lapply(params, function(p) {
    out <- param_find_boot(object = object, param = p, standardized = standardized)
    if (any(is.na(out))) stop("Bootstrap estimates not found or not stored for parameter: ", p)
    out$t
  })
  names(boot_out_list) <- params

  # 组装数据框
  df <- as.data.frame(boot_out_list, check.names = FALSE)

  # --- 自定义面板函数（GGally）
  library(GGally)
  library(ggplot2)

  # 下三角：散点 + 可选平滑 + 可选椭圆 + 边框
  lower_fun <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(size = point_size, alpha = point_alpha, color = point_color, ...)
    if (isTRUE(show_smooth)) {
      p <- p + geom_smooth(method = smooth_method, se = smooth_se,
                           linewidth = 0.8, color = dens_color, ...)
    }
    if (isTRUE(show_ellipse)) {
      p <- p + stat_ellipse(level = ellipse_level, linewidth = 0.8,
                            color = ellipse_color, linetype = "dashed", ...)
    }
    p +
      theme_minimal(base_family = "serif") +
      theme(
        panel.border = element_rect(color = panel_border_color, fill = NA,
                                    linewidth = panel_border_size)
      )
  }

  # 上三角：相关系数文本（放大） + 边框
  upper_fun <- if (isTRUE(show_corr)) {
    function(data, mapping, ...) {
      x <- ggplot2::ggplot_build(ggplot(data, mapping) + geom_point())
      df_num <- ggplot2::layer_data(x$plot, 1)
      r <- suppressWarnings(cor(df_num$x, df_num$y, use = "pairwise.complete.obs"))
      lbl <- if (is.finite(r)) sprintf("r = %.2f", round(r, corr_digits)) else "r = NA"
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = lbl,
                 size = corr_text_size, fontface = 2) +
        theme_minimal(base_family = "serif") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = panel_border_color, fill = NA,
                                      linewidth = panel_border_size)
        )
    }
  } else {
    lower_fun
  }

  # 对角线：密度/直方图/空白 + 边框
  diag_fun <- switch(
    diag_type,
    "density" = function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping) +
        geom_density(fill = dens_fill, color = dens_color, linewidth = 1.0, ...) +
        theme_minimal(base_family = "serif") +
        theme(
          panel.border = element_rect(color = panel_border_color, fill = NA,
                                      linewidth = panel_border_size)
        )
    },
    "hist" = function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping) +
        geom_histogram(bins = bins, fill = hist_fill, color = hist_color, ...) +
        theme_minimal(base_family = "serif") +
        theme(
          panel.border = element_rect(color = panel_border_color, fill = NA,
                                      linewidth = panel_border_size)
        )
    },
    "blank" = function(data, mapping, ...) {
      ggplot() +
        theme_void(base_family = "serif") +
        theme(
          panel.border = element_rect(color = panel_border_color, fill = NA,
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
    theme_minimal(base_family = "serif") +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = title)

  if (identical(return, "ggplot")) return(gp)

  print(gp)
  invisible(gp)
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
