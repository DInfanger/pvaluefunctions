# Plot of the results of conf_dist(). `res` contains the transformed and
# labelled frames, `text_frame` the labels of the confidence levels (NULL if
# there are none) and `n_est` the number of estimates.

plot_conf_dist <- function(
  res,
  text_frame,
  n_est,
  plot_type,
  alternative,
  conf_level,
  null_values,
  trans,
  trans_fun,
  x_scale,
  xlim,
  p_cutoff,
  log_yaxis,
  cut_logyaxis,
  inverted,
  together,
  same_color,
  plot_legend,
  col,
  nrow,
  ncol,
  plot_counternull,
  title,
  xlab,
  ylab,
  ylab_sec
) {
  # Create custom y-axis scale (mixed linear and logarithmic)

  # Transform cutoff for log-y-axis if applicable

  if (alternative %in% "one_sided") {
    if (cut_logyaxis > 0.5) {
      cut_logyaxis <- 0.5
    }
    cut_logyaxis_one <- cut_logyaxis
    cut_logyaxis <- cut_logyaxis_one * 2
  } else {
    cut_logyaxis_one <- cut_logyaxis / 2
  }

  # Labeller functions for custom log-scale

  lab_onesided <- Vectorize(function(x) {
    if (!is.na(x) && (x < cut_logyaxis_one) && (round((x %% 1) * 10) == 0)) {
      sprintf("%.5g", x)
    } else {
      sprintf("%.2f", x)
    }
  })

  lab_twosided <- Vectorize(function(x) {
    if (!is.na(x) && (x <= cut_logyaxis) && (round((x %% 1) * 10) == 0)) {
      sprintf("%.5g", x)
    } else {
      sprintf("%.1f", x)
    }
  })

  # Start plotting

  # Set variable to be plotted on y-axis depending on the type of the plot

  y_var <- switch(
    plot_type,
    p_val = "p_two",
    cdf = "conf_dist",
    pdf = "conf_dens",
    s_val = "s_val"
  )

  # Set label of the x-axis if not provided

  if (is.null(xlab)) {
    xlab <- "Estimate"
  }

  # Set label of the primary and secondary y-axis depending on the type of the plot

  if (is.null(ylab)) {
    ylab <- switch(
      plot_type,
      p_val = expression(paste(
        italic("P"),
        "-value (two-sided) / Significance level" ~ alpha,
        sep = ""
      )),
      cdf = "Confidence distribution",
      pdf = "Confidence density",
      s_val = expression(paste(
        "Surprisal in bits (two-sided ",
        ~ italic("P"),
        "-value)",
        sep = ""
      ))
    )
  }

  if (is.null(ylab_sec)) {
    ylab_sec <- switch(
      plot_type,
      p_val = expression(paste(
        italic("P"),
        "-value (one-sided) / Significance level" ~ alpha,
        sep = ""
      )),
      s_val = expression(paste(
        "Surprisal in bits (one-sided ",
        ~ italic("P"),
        "-value)",
        sep = ""
      ))
    )
  }

  # Create a ggplot2-object with no geoms

  p <- ggplot(
    res$res_frame,
    aes(x = .data$values, y = .data[[y_var]], group = .data$variable)
  ) +
    theme_bw()

  # Multiple estimates that are plotted in the same graph
  multi_together <- isTRUE(together) && n_est >= 2

  # If 2 or more estimates are plotted together, differentiate them by color (if user did not specify "same_color = TRUE")

  if (multi_together && !same_color) {
    p <- p + aes(colour = .data$variable)
  }

  if (alternative == "one_sided" && !multi_together) {
    # For only one one-sided p-value curve, set the colors to black and blue
    # (or to the specified color if "same_color = TRUE")
    curve_colours <- if (same_color) c(col, col) else c("black", "#08A9CF")

    p <- p +
      geom_line(aes(colour = .data$hypothesis), linewidth = 1.5) +
      scale_colour_manual(values = curve_colours) +
      theme(
        legend.position = "none"
      )
  } else if (
    alternative == "two_sided" && (!multi_together || isTRUE(same_color))
  ) {
    # For only one two-sided p-value curve, set the color to black
    p <- p + geom_line(linewidth = 1.5, colour = col)
  } else {
    # For 2 or more estimates plotted together: set the colors according to "Set1" palette
    p <- p +
      geom_line(linewidth = 1.5) +
      scale_colour_brewer(palette = "Set1", name = "") +
      theme(
        legend.position = if (isTRUE(plot_legend)) "top" else "none",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
      )
  }

  # Add the labels for the axes

  p <- p + xlab(xlab) + ylab(ylab)

  #-----------------------------------------------------------------------------
  # y-axis
  #-----------------------------------------------------------------------------

  # For p-value curves: Set the left and right y-axes, possibly with a logarithmic part

  if (plot_type %in% "p_val") {
    if (isTRUE(log_yaxis) && (p_cutoff < cut_logyaxis)) {
      lower_ylim_two <- round(10^(ceiling(round(log10(p_cutoff), 5))), 10)
      lower_ylim_one <- round(10^(ceiling(round(log10(p_cutoff / 2), 5))), 10)

      # Both axes need at least one power of ten below their cutoff
      if (
        lower_ylim_two <= cut_logyaxis && lower_ylim_one <= cut_logyaxis_one
      ) {
        # Split the breaks into two parts: i) below the cutoff i.e. the logarithmic part and ii) the linear part above the cutoff

        breaks_two <- c(
          10^(seq(log10(lower_ylim_two), log10(cut_logyaxis), by = 1)),
          seq(ceiling(cut_logyaxis / 0.1) * 0.1, 1, by = 0.1)
        )
        breaks_one <- c(
          10^(seq(log10(lower_ylim_one), log10(cut_logyaxis_one), by = 1)),
          seq(ceiling(cut_logyaxis_one / 0.05) * 0.05, 0.5, 0.05)
        )
      } else {
        breaks_two <- seq(ceiling(cut_logyaxis / 0.1) * 0.1, 1, by = 0.1)
        breaks_one <- seq(ceiling(cut_logyaxis_one / 0.05) * 0.05, 0.5, 0.05)
      }

      # Remove possible duplicates

      breaks_two <- unique(breaks_two)
      breaks_one <- unique(breaks_one)

      # The inverted axis uses the reversed limits and the reversed transformation
      y_limits <- if (isTRUE(inverted)) c(1, p_cutoff) else c(p_cutoff, 1)
      y_transform <- magnify_trans(
        rev = isTRUE(inverted),
        interval_low = cut_logyaxis,
        interval_high = 1,
        reducer = cut_logyaxis,
        reducer2 = 8
      )

      p <- p +
        scale_y_continuous(
          limits = y_limits,
          breaks = breaks_two,
          labels = lab_twosided,
          transform = y_transform,
          sec.axis = sec_axis(
            transform = ~ . * (1 / 2),
            name = ylab_sec,
            breaks = breaks_one,
            labels = lab_onesided
          )
        ) +
        # I() gives panel coordinates, which bypass the x-axis transformation
        # (-Inf or 0 would be NaN or -Inf on a logarithmic x-axis)
        annotate(
          "rect",
          xmin = I(0),
          xmax = I(1),
          ymin = p_cutoff,
          ymax = cut_logyaxis,
          alpha = 0.1,
          colour = "#E6E6E6"
        )
    } else {
      y_scale <- if (isTRUE(inverted)) scale_y_reverse else scale_y_continuous

      p <- p +
        y_scale(
          breaks = seq(0, 1, 0.1),
          sec.axis = sec_axis(
            ~ . * (1 / 2),
            name = ylab_sec,
            breaks = seq(0, 1, 0.1) / 2
          )
        )
    }
  }

  # For s-value curves: inverted the y-axis and transform it according to log2

  if (plot_type %in% "s_val") {
    s_val_max <- max(res$res_frame$s_val, na.rm = TRUE)

    # The inverted axis is a regular axis with the limits in ascending order
    # while the default axis is a reversed axis
    if (isTRUE(inverted)) {
      y_scale <- scale_y_continuous
      y_limits <- c(0, s_val_max)
    } else {
      y_scale <- scale_y_reverse
      y_limits <- c(s_val_max, 0)
    }

    p <- p +
      y_scale(
        limits = y_limits,
        breaks = scales::pretty_breaks(n = 10)(c(0, s_val_max)),
        sec.axis = sec_axis(
          transform = ~ . + log2(2),
          name = ylab_sec,
          breaks = scales::pretty_breaks(n = 10)(c(
            1,
            -log2(min(res$res_frame$p_two, na.rm = TRUE) / 2)
          ))
        )
      )
  }

  # Y-axis formatting For confidence distributions and densities

  if (plot_type %in% c("cdf", "pdf")) {
    # For the confidence distributions (plot_type == "cdf"), inverted y-axis if specified

    if (plot_type %in% c("cdf") && isTRUE(inverted)) {
      p <- p +
        scale_y_continuous(
          transform = "reverse",
          breaks = scales::pretty_breaks(n = 10)
        )
    } else {
      p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    }
  }

  #-----------------------------------------------------------------------------
  # x-axis
  #-----------------------------------------------------------------------------

  if (x_scale %in% "log") {
    # Plot x-axis on a log-scale (can't set "xlim" here, because then the gray area cannot be added!)
    p <- p +
      scale_x_continuous(
        transform = "log",
        breaks = scales::pretty_breaks(n = 10)
      )
  } else {
    p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  }

  # Set x-limits now

  p <- p + coord_cartesian(xlim = xlim, expand = TRUE)

  #-----------------------------------------------------------------------------
  # Horizontal lines at the significance levels if specified
  #-----------------------------------------------------------------------------

  if (plot_type %in% c("p_val", "s_val", "cdf") && !is.null(conf_level)) {
    hlines_tmp <- text_frame$p_value

    if (plot_type %in% "s_val") {
      hlines_tmp <- -log2(hlines_tmp)
    }

    p <- p + geom_hline(yintercept = hlines_tmp, linetype = 2)
  }

  #-----------------------------------------------------------------------------
  # Vertical lines at the null values if specified
  #-----------------------------------------------------------------------------

  if (!is.null(null_values)) {
    x_range <- ggplot_build(p)$layout$panel_params[[1]]$x.range

    if (trans %in% "exp" && x_scale %in% "log") {
      # If the x-axis was log-transformed, we need to backtransform the plotting limits because they are given on the log-scale
      plot_limits <- trans_fun(x_range)
    } else if (!trans %in% "exp" && x_scale %in% "log") {
      plot_limits <- exp(x_range)
    } else {
      plot_limits <- x_range
    }

    # Which null_values are outside of the plotting limits

    null_outside_plot <- which(
      (res$counternull_frame$null_value <= plot_limits[1]) |
        (res$counternull_frame$null_value >= plot_limits[2])
    )

    # Only add lines for those null values that are inside the plotting limits

    null_lines <- res$counternull_frame

    if (length(null_outside_plot) > 0) {
      # Print a message that shows which null values were outside of the x-axis limits.

      nulls_outside <- unique(
        res$counternull_frame$null_value[null_outside_plot]
      )
      cli::cli_inform(
        "The following null values are outside of the specified x-axis range (xlim) and are not shown: {.val {nulls_outside}}"
      )

      null_lines <- null_lines[-null_outside_plot, ]
    }

    if (nrow(null_lines) > 0L) {
      p <- p +
        geom_vline(
          data = null_lines,
          aes(xintercept = .data$null_value),
          linetype = 1,
          linewidth = 0.5
        )
    }
  }

  #-----------------------------------------------------------------------------
  # Facets for multiple estimates not plotted together
  #-----------------------------------------------------------------------------

  if (n_est >= 2 && isFALSE(together)) {
    p <- p +
      facet_wrap(
        vars(.data$variable),
        nrow = nrow,
        ncol = ncol,
        scales = if (plot_type %in% "pdf") "free" else "free_x"
      )
  }

  #-----------------------------------------------------------------------------
  # Add text boxes at the specified significance levels, if any
  #-----------------------------------------------------------------------------

  if (!plot_type %in% "pdf" && !is.null(conf_level)) {
    if (plot_type %in% "s_val") {
      text_frame$p_value <- -log2(text_frame$p_value)
    }

    # Boxes without border; `label.size` was deprecated in ggplot2 4.0.0
    border_args <- if (utils::packageVersion("ggplot2") >= "4.0.0") {
      list(linewidth = 0)
    } else {
      list(label.size = NA)
    }

    p <- p +
      do.call(
        geom_label,
        c(
          list(
            data = text_frame,
            mapping = aes(
              x = .data$theor_values,
              y = .data$p_value,
              label = .data$label
            ),
            inherit.aes = FALSE,
            parse = TRUE,
            size = 5.5,
            hjust = "inward"
          ),
          border_args
        )
      )
  }

  #-----------------------------------------------------------------------------
  # Add points for the counternull if specified
  #-----------------------------------------------------------------------------

  if (
    !is.null(null_values) &&
      isTRUE(plot_counternull) &&
      (plot_type %in% c("p_val", "s_val")) &&
      !all(is.na(res$res_frame$counternull))
  ) {
    if (multi_together && !same_color) {
      p <- p +
        geom_point(
          aes(x = .data$values, y = .data$counternull, colour = .data$variable),
          size = 4,
          shape = 21,
          fill = "white",
          stroke = 1.7
        ) +
        guides(colour = guide_legend(override.aes = list(pch = NA)))
    } else {
      p <- p +
        geom_point(
          aes(x = .data$values, y = .data$counternull),
          colour = col,
          size = 4,
          shape = 21,
          fill = "white",
          stroke = 1.7
        )
    }
  }

  #-----------------------------------------------------------------------------
  # Make the plot prettier by increasing font size
  #-----------------------------------------------------------------------------

  p <- p +
    theme(
      axis.title.y.left = element_text(
        colour = "black",
        size = 17,
        hjust = 0.5,
        margin = margin(0, 10, 0, 0)
      ),
      axis.title.y.right = element_text(
        colour = "black",
        size = 17,
        hjust = 0.5,
        margin = margin(0, 0, 0, 10)
      ),
      axis.title.x = element_text(colour = "black", size = 17),
      axis.text.x = element_text(colour = "black", size = 15),
      axis.text.y = element_text(colour = "black", size = 15),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text.x = element_text(size = 15)
    )

  #-----------------------------------------------------------------------------
  # Add title and labels if specified
  #-----------------------------------------------------------------------------

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }

  p
}
