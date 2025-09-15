# tests/testthat/test-jab_after_boot.R

testthat::test_that("jab_after_boot() errors cleanly when boot.idx missing", {
  testthat::skip_if_not_installed("lavaan")

  set.seed(1)
  HS <- lavaan::HolzingerSwineford1939
  model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
  '
  fit0 <- lavaan::cfa(model, data = HS, std.lv = TRUE)

  # 确认还没有 sbt_boot_ustd
  testthat::expect_true(is.null(fit0@external$sbt_boot_ustd))

  # 由于没有 keep.idx 的 bootstrap 索引，应报错
  testthat::expect_error(
    jab_after_boot(fit0, param = "x1~~x1", standardized = FALSE, verbose = FALSE),
    "sbt_boot_ustd|keep.idx"
  )
})

testthat::test_that("jab_after_boot() works for unstandardized (:= and free/var params)", {
  testthat::skip_if_not_installed("lavaan")

  # 允许本地通过环境变量控制 Bootstrap 次数（默认 200，避免过慢）
  R_boot <- as.integer(Sys.getenv("SBT_R", "200"))
  if (isTRUE(as.logical(Sys.getenv("SBT_SKIP_SLOW", "FALSE")))) {
    testthat::skip("slow integration test skipped by SBT_SKIP_SLOW")
  }

  set.seed(1234)
  n <- 200
  x <- runif(n) - 0.5
  m <- 0.4 * x + rnorm(n)
  y <- 0.3 * m + rnorm(n)
  dat <- data.frame(x, m, y)

  model <- '
    m ~ a * x
    y ~ b * m + cp * x
    ab := a * b
  '

  fit0 <- lavaan::sem(model, data = dat, fixed.x = FALSE)

  # 关键：keep.idx = TRUE，供 JAB 取 LOO 子分布
  testthat::expect_true(exists("store_boot", mode = "function"))
  fit2 <- store_boot(
    fit0,
    R        = R_boot,
    iseed    = 2345,
    keep.idx = TRUE
  )

  # ---- 分支 A：:= 定义参数（ab） ----
  # 小样本下会有“p-values not computed”警告，这里显式断言该警告
  testthat::expect_warning(
    res_ab <- jab_after_boot(
      fit2,
      param = "ab",
      standardized = FALSE,
      top_k = 5L,
      plot = FALSE,
      verbose = FALSE
    ),
    regexp = "not computed",
    ignore.case = TRUE
  )
  testthat::expect_type(res_ab, "list")
  testthat::expect_s3_class(res_ab$full_summary, "data.frame")
  testthat::expect_s3_class(res_ab$cases_summary, "data.frame")
  testthat::expect_true(is.matrix(res_ab$F))
  testthat::expect_type(res_ab$tstar, "double")
  testthat::expect_lte(nrow(res_ab$cases_summary), 5L)

  # ---- 分支 B：自由参数/方差列（自动选择存在的列名） ----
  ub <- suppressWarnings(parameterEstimates_boot(fit2))
  cols <- colnames(attr(ub, "boot_est_ustd"))
  # 你的实现里自由路径用 label（a,b,cp），方差用 lhs~~lhs
  cand <- intersect(c("a","b","cp"), cols)
  if (length(cand) == 0L) cand <- grep("~~", cols, value = TRUE)
  testthat::expect_true(length(cand) >= 1L)
  p_unstd <- cand[[1L]]

  testthat::expect_warning(
    res_path <- jab_after_boot(
      fit2,
      param = p_unstd,
      standardized = FALSE,
      top_k = 3L,
      plot = FALSE,
      verbose = FALSE
    ),
    regexp = "not computed",
    ignore.case = TRUE
  )
  testthat::expect_type(res_path$tstar, "double")
  testthat::expect_lte(nrow(res_path$cases_summary), 3L)

  # base 绘图分支：不应报错
  # base 绘图分支：同样会出现小样本警告，这里也显式断言
  testthat::expect_warning(
    jab_after_boot(
      fit2, param = p_unstd, standardized = FALSE,
      plot = TRUE, plot_engine = "base", verbose = FALSE
    ),
    regexp = "not computed",
    ignore.case = TRUE
  )
})

testthat::test_that("jab_after_boot() standardized branch + ggplot2 plotting", {
  testthat::skip_if_not_installed("lavaan")
  testthat::skip_if_not_installed("ggplot2")

  R_boot <- as.integer(Sys.getenv("SBT_R", "200"))
  if (isTRUE(as.logical(Sys.getenv("SBT_SKIP_SLOW", "FALSE")))) {
    testthat::skip("slow integration test skipped by SBT_SKIP_SLOW")
  }

  set.seed(5678)
  n <- 220
  x <- runif(n) - 0.5
  m <- 0.45 * x + rnorm(n)
  y <- 0.35 * m + rnorm(n)
  dat <- data.frame(x, m, y)

  model <- '
    m ~ a * x
    y ~ b * m + cp * x
    ab := a * b
  '

  fit0 <- lavaan::sem(model, data = dat, fixed.x = FALSE)
  fit2 <- store_boot(fit0, R = R_boot, iseed = 9876, keep.idx = TRUE)

  # ---- 自动挑一个标准化列名（优先回归，其次载荷，再次方差） ----
  std <- suppressWarnings(standardizedSolution_boot(fit2))
  std_cols <- colnames(attr(std, "boot_est_std"))

  # 回归：包含 "~" 且不包含 "=~" 或 "~~"
  reg_cols <- std_cols[grepl("~", std_cols) & !grepl("=~|~~", std_cols)]
  # 载荷：包含 "=~"
  load_cols <- std_cols[grepl("=~", std_cols)]
  # 方差：包含 "~~"
  var_cols <- std_cols[grepl("~~", std_cols)]

  p_std <- NA_character_
  if (length(reg_cols) > 0L) {
    p_std <- reg_cols[[1L]]
  } else if (length(load_cols) > 0L) {
    p_std <- load_cols[[1L]]
  } else if (length(var_cols) > 0L) {
    p_std <- var_cols[[1L]]
  }
  testthat::expect_true(is.character(p_std) && nchar(p_std) > 0L)

  # 标准化分支：同样断言小样本警告
  testthat::expect_warning(
    res_std <- jab_after_boot(
      fit2,
      param = p_std,
      standardized = TRUE,
      top_k = 4L,
      plot = FALSE,
      verbose = FALSE
    ),
    regexp = "not computed",
    ignore.case = TRUE
  )
  testthat::expect_true(res_std$standardized)
  testthat::expect_s3_class(res_std$full_summary, "data.frame")
  testthat::expect_lte(nrow(res_std$cases_summary), 4L)

  # ggplot2 分支：若 LOO 有效，应返回 ggplot 对象；若全 NA，允许为 NULL
  testthat::expect_warning(
    res_plot <- jab_after_boot(
      fit2,
      param = p_std,
      standardized = TRUE,
      plot = TRUE,
      plot_engine = "ggplot2",
      verbose = FALSE
    ),
    regexp = "not computed",
    ignore.case = TRUE
  )
  if (!is.null(res_plot$plot_obj)) {
    testthat::expect_s3_class(res_plot$plot_obj, "ggplot")
  }
})

testthat::test_that("jab_after_boot() gives informative error for unknown parameter", {
  testthat::skip_if_not_installed("lavaan")

  set.seed(2468)
  n <- 180
  x <- runif(n) - 0.5
  m <- 0.4 * x + rnorm(n)
  y <- 0.3 * m + rnorm(n)
  dat <- data.frame(x, m, y)

  model <- '
    m ~ a * x
    y ~ b * m + cp * x
    ab := a * b
  '
  fit0 <- lavaan::sem(model, data = dat, fixed.x = FALSE)
  fit2 <- store_boot(fit0, R = 120, iseed = 1111, keep.idx = TRUE)

  testthat::expect_error(
    suppressWarnings(
      jab_after_boot(
        fit2,
        param = "NOT_A_PARAM",
        standardized = TRUE,
        verbose = FALSE
      )
    ),
    "not found|Cannot locate|Could not resolve",
    ignore.case = TRUE
  )
})
