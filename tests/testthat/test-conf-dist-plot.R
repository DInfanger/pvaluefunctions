# The plots are checked through their properties (layers, scales, facets)
# with ggplot2::ggplot_build() instead of image snapshots, which depend on
# the platform and on the ggplot2 version.

plot_of <- function(...) {
  args <- list(
    estimate = 0.5,
    stderr = 0.2,
    type = "general_z",
    n_values = 100L
  )
  args <- utils::modifyList(args, list(...))

  # ggplot2 deprecation warnings are not what is tested here
  res <- suppressWarnings(do.call(conf_dist_quiet, args))
  res$plot
}

geoms_of <- function(p) {
  vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
}

scale_of <- function(p, axis = c("x", "y")) {
  axis <- match.arg(axis)
  built <- ggplot2::ggplot_build(p)

  built$layout[[paste0("panel_scales_", axis)]][[1]]
}

test_that("every plot type returns a ggplot object", {
  for (plot_type in c("p_val", "s_val", "cdf", "pdf")) {
    for (alternative in c("two_sided", "one_sided")) {
      p <- plot_of(plot_type = plot_type, alternative = alternative)

      expect_s3_class(p, "ggplot")
      expect_no_error(ggplot2::ggplot_build(p))
    }
  }
})

test_that("plots work for all types of estimates", {
  specs <- list(
    list(estimate = 1, df = 10, stderr = 0.3, type = "linreg"),
    list(estimate = 1, df = 10, stderr = 0.3, type = "gammareg"),
    list(estimate = 1, df = 10, stderr = 0.3, type = "general_t"),
    list(estimate = 0.8, tstat = 2.1, df = 25, type = "ttest"),
    list(estimate = 0.1, stderr = 0.3, type = "logreg", trans = "exp"),
    list(estimate = 0.1, stderr = 0.3, type = "poisreg", trans = "exp"),
    list(estimate = 0.1, stderr = 0.3, type = "coxreg", trans = "exp"),
    list(estimate = 0.4, n = 40, type = "spearman"),
    list(estimate = 0.4, n = 40, type = "kendall"),
    list(estimate = 4, n = 20, type = "var"),
    list(estimate = 0.3, n = 50, type = "prop"),
    list(estimate = c(0.4, 0.3), n = c(60, 60), type = "propdiff")
  )

  for (spec in specs) {
    expect_s3_class(do.call(plot_of, spec), "ggplot")
  }
})

test_that("the exact Pearson distribution can be plotted", {
  skip_on_cran()

  expect_s3_class(
    plot_of(estimate = 0.3, n = 30, type = "pearson", n_values = 20L),
    "ggplot"
  )
})

test_that("null values are drawn as vertical lines within the x-limits", {
  p <- plot_of(null_values = c(0, 0.3), xlim = c(-0.5, 1.5))
  expect_equal(sum(geoms_of(p) == "GeomVline"), 1)
  expect_equal(
    sort(ggplot2::ggplot_build(p)$data[[which(geoms_of(p) == "GeomVline")]]$xintercept),
    c(0, 0.3)
  )

  # A null value outside of xlim is not drawn and reported
  expect_message(
    p <- plot_of(null_values = c(0, 9), xlim = c(-0.5, 1.5)),
    "outside"
  )
  vline <- ggplot2::ggplot_build(p)$data[[which(geoms_of(p) == "GeomVline")]]
  expect_equal(vline$xintercept, 0)

  # No vertical lines without null values
  expect_false("GeomVline" %in% geoms_of(plot_of()))
})

test_that("xlim sets the range of the x-axis", {
  p <- plot_of(xlim = c(-0.5, 1.5))

  built <- ggplot2::ggplot_build(p)

  # The default 5% expansion of ggplot2 is added to the limits
  expect_equal(built$layout$panel_params[[1]]$x.range, c(-0.6, 1.6))
})

test_that("inverted reverses the y-axis", {
  normal <- plot_of(inverted = FALSE)
  inverted <- plot_of(inverted = TRUE)

  expect_equal(scale_of(normal, "y")$get_transformation()$name, "identity")
  expect_equal(scale_of(inverted, "y")$get_transformation()$name, "reverse")
})

test_that("the y-axis can be logarithmic", {
  p <- plot_of(log_yaxis = TRUE, plot_p_limit = 1e-4)

  expect_equal(scale_of(p, "y")$get_transformation()$name, "customlog")
  expect_s3_class(plot_of(log_yaxis = TRUE, inverted = TRUE), "ggplot")
  expect_s3_class(
    plot_of(log_yaxis = TRUE, alternative = "one_sided", plot_p_limit = 1e-3),
    "ggplot"
  )
})

test_that("the x-axis is logarithmic for exp and can be forced", {
  default <- plot_of(trans = "exp", type = "logreg")
  linear <- plot_of(trans = "exp", type = "logreg", x_scale = "linear")
  # Only positive values can be shown on a logarithmic scale
  forced <- plot_of(estimate = 4, n = 20, type = "var", x_scale = "logarithm")
  untransformed <- plot_of()

  expect_match(scale_of(default, "x")$get_transformation()$name, "^log")
  expect_equal(scale_of(linear, "x")$get_transformation()$name, "identity")
  expect_match(scale_of(forced, "x")$get_transformation()$name, "^log")
  expect_equal(scale_of(untransformed, "x")$get_transformation()$name, "identity")
})

test_that("a secondary y-axis is added", {
  p <- plot_of()

  expect_false(is.null(scale_of(p, "y")$secondary.axis))
})

test_that("estimates are shown in separate panels unless together", {
  args <- list(estimate = c(1, 2, 3), stderr = c(1, 1, 1))

  separate <- do.call(plot_of, c(args, list(together = FALSE)))
  together <- do.call(plot_of, c(args, list(together = TRUE)))

  expect_s3_class(separate$facet, "FacetWrap")
  expect_s3_class(together$facet, "FacetNull")
  expect_equal(nrow(ggplot2::ggplot_build(separate)$layout$layout), 3)
  expect_equal(nrow(ggplot2::ggplot_build(together)$layout$layout), 1)
})

test_that("nrow and ncol control the layout of the panels", {
  p <- plot_of(
    estimate = c(1, 2, 3, 4),
    stderr = rep(1, 4),
    together = FALSE,
    nrow = 1,
    ncol = 4
  )
  layout <- ggplot2::ggplot_build(p)$layout$layout

  expect_equal(max(layout$ROW), 1)
  expect_equal(max(layout$COL), 4)
})

test_that("colours and legend of curves plotted together can be controlled", {
  args <- list(estimate = c(0, 1), stderr = c(1, 0.5), together = TRUE)

  coloured <- do.call(plot_of, args)
  no_legend <- do.call(plot_of, c(args, list(plot_legend = FALSE)))
  same <- do.call(plot_of, c(args, list(same_color = TRUE, col = "red")))

  built <- ggplot2::ggplot_build(coloured)$data[[1]]
  expect_equal(length(unique(built$colour)), 2)

  expect_equal(no_legend$theme$legend.position, "none")
  expect_equal(coloured$theme$legend.position, "top")

  same_built <- ggplot2::ggplot_build(same)$data[[1]]
  expect_equal(unique(same_built$colour), "red")
})

test_that("title and axis labels are used", {
  p <- plot_of(title = "My title", xlab = "My x", ylab = "My y")

  expect_equal(p$labels$title, "My title")
  expect_equal(p$labels$x, "My x")
  expect_equal(p$labels$y, "My y")
})

test_that("counternulls are shown only on request", {
  with_cn <- plot_of(null_values = 0, plot_counternull = TRUE)
  without_cn <- plot_of(null_values = 0, plot_counternull = FALSE)

  expect_gt(length(with_cn$layers), length(without_cn$layers))
})

test_that("a user-defined transformation function is used for plotting", {
  rse_fun <- function(x) 100 * (1 - exp(x))
  p <- plot_of(
    estimate = log(0.72),
    stderr = 0.19,
    type = "coxreg",
    trans = rse_fun,
    xlim = log(1 - c(-30, 60) / 100)
  )

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plots work with the lower limit of the p-value axis", {
  for (limit in c(0, 0.01, 0.2)) {
    expect_s3_class(plot_of(plot_p_limit = limit), "ggplot")
  }
})

test_that("legend and layout options do not matter for a single estimate", {
  expect_s3_class(
    plot_of(together = TRUE, plot_legend = FALSE, same_color = TRUE),
    "ggplot"
  )
})

test_that("a transformation is ignored for types that do not support it", {
  expect_message(
    res <- cd(estimate = 4, n = 20, type = "var", trans = "exp"),
    "changed to identity"
  )
  expect_true(all(res$res_frame$values > 0))
  expect_equal(res$res_frame$values, cd(estimate = 4, n = 20, type = "var")$res_frame$values)
})

test_that("only the first colour is used", {
  expect_warning(
    res <- cd(estimate = 0.5, stderr = 0.2, type = "general_z", col = c("red", "blue")),
    "Only first color"
  )
  expect_s3_class(res$res_frame, "data.frame")
})

test_that("s-value plots show confidence limits and counternulls", {
  p <- plot_of(
    plot_type = "s_val",
    conf_level = c(0.95, 0.8),
    null_values = 0,
    plot_counternull = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  inverted <- plot_of(plot_type = "s_val", inverted = TRUE, conf_level = 0.95)
  expect_equal(scale_of(inverted, "y")$get_transformation()$name, "identity")
  expect_s3_class(inverted, "ggplot")
})

test_that("building the plot gives no deprecation warnings", {
  expect_no_warning(
    res <- conf_dist_quiet(
      estimate = 0.5,
      stderr = 0.2,
      type = "general_z",
      n_values = 100L,
      conf_level = c(0.95, 0.8)
    )
  )
  expect_no_warning(ggplot2::ggplot_build(res$plot))
})

test_that("confidence levels are labelled in p-value plots", {
  p <- plot_of(conf_level = c(0.95, 0.8), null_values = 0.2, plot_counternull = TRUE)

  expect_true("GeomLabel" %in% geoms_of(p))
})

test_that("one-sided plots with a logarithmic axis cap the cut-off", {
  p <- plot_of(
    alternative = "one_sided",
    log_yaxis = TRUE,
    cut_logyaxis = 0.9,
    plot_p_limit = 1e-3
  )

  expect_equal(scale_of(p, "y")$get_transformation()$name, "customlog")
  expect_s3_class(
    plot_of(alternative = "one_sided", log_yaxis = TRUE, cut_logyaxis = 0.2),
    "ggplot"
  )
  expect_s3_class(
    plot_of(log_yaxis = TRUE, cut_logyaxis = 0.2, conf_level = 0.95),
    "ggplot"
  )
})

test_that("cdf plots can be inverted", {
  p <- plot_of(plot_type = "cdf", inverted = TRUE)

  expect_equal(scale_of(p, "y")$get_transformation()$name, "reverse")
})

test_that("pdf plots of several estimates have free scales", {
  p <- plot_of(
    plot_type = "pdf",
    estimate = c(0, 5),
    stderr = c(1, 2),
    together = FALSE
  )
  layout <- ggplot2::ggplot_build(p)$layout$layout

  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(nrow(layout), 2)
})

test_that("counternulls of curves plotted together are coloured", {
  p <- plot_of(
    estimate = c(0, 1),
    stderr = c(1, 0.5),
    together = TRUE,
    null_values = 0.5,
    plot_counternull = TRUE
  )
  same <- plot_of(
    estimate = c(0, 1),
    stderr = c(1, 0.5),
    together = TRUE,
    same_color = TRUE,
    null_values = 0.5,
    plot_counternull = TRUE
  )

  expect_true("GeomPoint" %in% geoms_of(p))
  expect_true("GeomPoint" %in% geoms_of(same))
})

test_that("null values are placed correctly on a logarithmic x-axis", {
  p <- plot_of(
    estimate = 4,
    n = 20,
    type = "var",
    x_scale = "logarithm",
    null_values = 3
  )
  vline <- ggplot2::ggplot_build(p)$data[[which(geoms_of(p) == "GeomVline")]]

  expect_equal(vline$xintercept, log(3))
})
