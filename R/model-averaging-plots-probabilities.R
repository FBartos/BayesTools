#' Plot posterior model probabilities
#'
#' Draw aligned strips whose widths are the posterior probabilities of the
#' models. Each row uses colour to show the supplied state of each model.
#'
#' @param model_probabilities Named numeric vector of posterior model
#'   probabilities. The names must be unique model IDs and the probabilities
#'   must sum to one up to numerical tolerance.
#' @param cells Data frame with character columns `model` (model ID), `row`
#'   (row ID), `label` (row label), and `state`. The optional character column
#'   `group` provides group labels, with `""` for ungrouped rows. Supply
#'   exactly one cell per model and row, with the same label and group for
#'   every cell of a row. Rows appear from top to bottom in order of first
#'   occurrence in `cells`.
#' @param plot_type Either `"base"` or `"ggplot"`.
#' @param ... Reserved for future use.
#' @param group_labels Whether to label consecutive groups of rows in the left
#'   margin. A group that reappears later has its label repeated.
#' @param labels Named character vector of optional label overrides, keyed
#'   by row ID.
#' @param decreasing Whether to sort models from highest to lowest probability.
#'   Set to `FALSE` to sort from lowest to highest. Ties retain the order in
#'   `model_probabilities`.
#' @param col Named character vector of colours keyed by state. Built-in states
#'   are `"positive"` (blue), `"negative"` (red), `"included"` (black), and
#'   `"excluded"` (white). Supply a colour for each other state used in
#'   `cells`.
#' @param legend Whether to show a legend for visible states.
#' @param legend_position Place the legend outside the panel at `"bottom"` or
#'   `"top"`.
#' @param state_labels Named character vector of optional legend labels keyed
#'   by state. Custom states use their keys as labels by default.
#'
#' @details
#' The caller supplies the model states. For example, directional states can
#' describe the sign of a within-model estimate, while inclusion states can
#' describe components without a signed estimate. Cells in a model share their
#' probability interval across all rows. Models with zero probability have
#' no visible cells. To select or reorder rows, pass the corresponding cells
#' in the desired row order.
#'
#' @return With `plot_type = "base"`, draws a plot and returns `NULL`
#'   invisibly. With `plot_type = "ggplot"`, returns a ggplot object with an
#'   integer `"n_rows"` attribute giving the number of displayed rows.
#' @seealso [plot_models()]
#' @examples
#' probabilities <- c(null = 0.2, effect = 0.5, full = 0.3)
#' cells <- data.frame(
#'   model = rep(names(probabilities), 2),
#'   row = rep(c("effect", "heterogeneity"), each = 3),
#'   label = rep(c("Effect", "Heterogeneity"), each = 3),
#'   state = c("excluded", "positive", "positive",
#'             "excluded", "excluded", "included")
#' )
#' plot_model_probabilities(probabilities, cells, legend = TRUE)
#' @md
#' @export
plot_model_probabilities <- function(model_probabilities, cells, plot_type = "base", ..., group_labels = FALSE,
                                     labels = NULL, decreasing = TRUE, col = NULL, legend = FALSE,
                                     legend_position = "bottom", state_labels = NULL){

  rlang::check_dots_empty()
  check_bool(group_labels)
  check_bool(decreasing)
  check_bool(legend)
  check_char(legend_position, allow_values = c("bottom", "top"), allow_NA = FALSE)
  check_char(plot_type, allow_values = c("base", "ggplot"), allow_NA = FALSE)

  data <- .model_probability_plot_data(
    model_probabilities, cells, group_labels, labels, decreasing, col, state_labels
  )
  legend_at <- if(legend) legend_position else "none"
  if(plot_type == "ggplot"){
    return(.model_probability_plot_ggplot(data, legend_at))
  }
  .model_probability_plot_base(data, legend_at)
  invisible(NULL)
}

.model_probability_check_names <- function(x, name){

  keys <- names(x)
  if(length(x) && (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys))){
    stop("The '", name, "' argument must have unique, nonempty names.", call. = FALSE)
  }
}

.model_probability_plot_data <- function(model_probabilities, cells, group_labels, labels,
                                        decreasing, col, state_labels){

  check_real(model_probabilities, lower = 0, upper = 1, check_length = 0, allow_NA = FALSE)
  .model_probability_check_names(model_probabilities, "model_probabilities")
  if(abs(sum(model_probabilities) - 1) > sqrt(.Machine$double.eps)){
    stop("The 'model_probabilities' argument must sum to one.", call. = FALSE)
  }

  columns <- c("model", "row", "label", "state")
  if(!is.data.frame(cells) || !nrow(cells) ||
     !all(columns %in% names(cells)) || anyDuplicated(names(cells))){
    stop("The 'cells' argument must be a nonempty data frame with 'model', 'row', 'label', and 'state' columns.",
         call. = FALSE)
  }
  cells <- as.data.frame(cells)
  if(!"group" %in% names(cells)){
    cells$group <- ""
  }
  cells <- cells[c(columns, "group")]
  for(column in names(cells)){
    check_char(cells[[column]], paste0("cells$", column), check_length = 0, allow_NA = FALSE)
  }
  if(any(!nzchar(cells$row)) || any(!nzchar(cells$label)) || any(!nzchar(cells$state))){
    stop("The 'row', 'label', and 'state' columns in 'cells' must contain nonempty strings.", call. = FALSE)
  }
  rows <- cells[!duplicated(cells$row), c("row", "label", "group")]
  if(any(!cells$model %in% names(model_probabilities)) ||
     nrow(cells) != length(model_probabilities) * nrow(rows) || anyDuplicated(cells[c("row", "model")])){
    stop("The 'cells' argument must contain exactly one cell for each model and row.", call. = FALSE)
  }
  row_index <- match(cells$row, rows$row)
  if(any(cells$label != rows$label[row_index]) || any(cells$group != rows$group[row_index])){
    stop("The 'label' and 'group' in 'cells' must be constant within each row.", call. = FALSE)
  }

  check_char(labels, check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  .model_probability_check_names(labels, "labels")
  label_rows <- match(names(labels), rows$row)
  if(anyNA(label_rows) || any(!nzchar(labels))){
    stop("The 'labels' argument must contain nonempty labels keyed by displayed row IDs.", call. = FALSE)
  }
  if(length(labels)){
    rows$label[label_rows] <- unname(labels)
  }

  palette <- c(positive = "#316DA4", negative = "#BD4E48", included = "#252525", excluded = "#FFFFFF")
  legend_labels <- c(positive = "Positive", negative = "Negative", included = "Included", excluded = "Excluded")
  check_char(col, check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  .model_probability_check_names(col, "col")
  palette[names(col)] <- col
  missing_colors <- setdiff(unique(cells$state), names(palette))
  if(length(missing_colors)){
    stop("Supply a colour in 'col' for each state without a default: '", paste(missing_colors, collapse = "', '"),
         "'.", call. = FALSE)
  }
  rgba <- tryCatch(grDevices::col2rgb(palette, alpha = TRUE), error = function(e){
    stop("The 'col' argument must contain valid R colours.", call. = FALSE)
  })
  check_char(state_labels, check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  .model_probability_check_names(state_labels, "state_labels")
  if(any(!names(state_labels) %in% names(palette)) || any(!nzchar(state_labels))){
    stop("The 'state_labels' argument must contain nonempty labels keyed by states with a colour.", call. = FALSE)
  }
  custom_states <- setdiff(names(palette), names(legend_labels))
  legend_labels[custom_states] <- custom_states
  legend_labels[names(state_labels)] <- state_labels

  n_rows <- nrow(rows)
  rows$y <- rev(seq_len(n_rows))
  groups <- rle(if(group_labels) rows$group else rep("", n_rows))
  rows$run <- factor(rep(seq_along(groups$lengths), groups$lengths))
  runs <- data.frame(
    run = levels(rows$run),
    label = groups$values,
    y = rows$y[cumsum(groups$lengths)] + (groups$lengths - 1) / 2
  )
  row_rules <- rows[seq_len(n_rows - 1L), c("run", "y")]
  row_rules$y <- row_rules$y - 0.5
  line_colour <- "#939DA5"
  row_rules$colour <- ifelse(diff(as.integer(rows$run)) != 0L, line_colour, "#A6AFB685")

  model_probabilities <- model_probabilities[order(model_probabilities, decreasing = decreasing, method = "radix")]
  ends <- c(0, cumsum(model_probabilities))
  model_index <- match(cells$model, names(model_probabilities))
  cells$xmin <- ends[model_index]
  cells$xmax <- ends[model_index + 1L]
  cells$ymin <- rows$y[row_index] - 0.5
  cells$ymax <- rows$y[row_index] + 0.5
  cells$run <- rows$run[row_index]
  cells <- cells[order(row_index, model_index), ]
  cells <- cells[cells$xmax > cells$xmin & cells$xmin < 1, ]

  # Use pale dots where both neighbouring fills are dark.
  brightness <- drop(c(0.21, 0.72, 0.07) %*% rgba[1:3, , drop = FALSE])
  dark <- stats::setNames(255 + (brightness - 255) * rgba[4L, ] / 255 < 145, names(palette))
  left <- seq_len(nrow(cells) - 1L)
  right <- left + 1L
  adjacent <- cells$row[left] == cells$row[right]
  boundaries <- cells[right[adjacent], c("xmin", "ymin", "ymax", "run")]
  pale <- dark[cells$state[left[adjacent]]] & dark[cells$state[right[adjacent]]]
  boundaries$colour <- ifelse(pale, "#FFFFFF80", "#68747E80")

  visible_states <- names(palette)[names(palette) %in% cells$state]
  list(
    cells = cells, rows = rows, runs = runs, row_rules = row_rules, boundaries = boundaries,
    palette = palette[visible_states], legend_labels = legend_labels[visible_states],
    ticks = c(0, 0.25, 0.50, 0.75, 1), tick_labels = c("0", "0.25", "0.50", "0.75", "1"),
    axis_title = "Posterior Probability", line_colour = line_colour
  )
}

.model_probability_plot_base <- function(data, legend_at){

  # Advance to the next layout panel before measuring its figure dimensions.
  old_par <- graphics::par(mai = rep(0, 4L))
  on.exit(graphics::par(old_par))
  graphics::plot.new()

  rows <- data$rows
  runs <- data$runs[nzchar(data$runs$label), ]
  cells <- data$cells
  boundaries <- data$boundaries
  row_rules <- data$row_rules
  line_inches <- graphics::par("csi") * graphics::par("mex")
  figure_inches <- graphics::par("fin")
  axis_cex <- 0.8 * graphics::par("cex.axis")
  group_cex <- axis_cex * 0.9
  left_margin <- max(graphics::strwidth(rows$label, units = "inches", cex = axis_cex)) + 0.2
  if(nrow(runs)){
    group_width <- graphics::strwidth(runs$label, units = "inches", cex = group_cex, font = 2)
    left_margin <- left_margin + max(group_width) + 0.2
  }

  legend_cex <- 0.85
  margins <- c(3 * line_inches, left_margin, 0.15, 0.2)
  if(legend_at != "none"){
    legend_width <- max(graphics::strwidth(data$legend_labels, units = "inches", cex = legend_cex)) + 0.55
    legend_cols <- max(1L, min(4L, length(data$legend_labels), floor(figure_inches[1L] / legend_width)))
    legend_rows <- ceiling(length(data$legend_labels) / legend_cols)
    legend_side <- if(legend_at == "top") 3L else 1L
    margins[legend_side] <- margins[legend_side] + (legend_rows + 0.5) * line_inches
  }
  graphics::par(mai = margins)
  graphics::plot.window(c(0, 1), c(0.5, nrow(rows) + 0.5), xaxs = "i", yaxs = "i")
  graphics::rect(cells$xmin, cells$ymin, cells$xmax, cells$ymax, col = data$palette[cells$state], border = NA)
  if(nrow(boundaries)){
    graphics::segments(boundaries$xmin, boundaries$ymin, boundaries$xmin, boundaries$ymax,
                       col = boundaries$colour, lty = 3, lwd = 0.7)
  }
  if(nrow(row_rules)){
    graphics::segments(0, row_rules$y, 1, row_rules$y, lwd = 0.7, col = row_rules$colour)
  }
  graphics::box(col = data$line_colour, lwd = 0.7)
  graphics::axis(1L, at = data$ticks, labels = data$tick_labels, lwd = 0, lwd.ticks = 0.7,
                 col.ticks = data$line_colour, cex.axis = axis_cex, tcl = -0.2, mgp = c(2, 0.55, 0))
  graphics::title(xlab = data$axis_title, line = 1.85, cex.lab = 0.95 * graphics::par("cex.lab"))

  label_x <- graphics::grconvertX((left_margin - 0.1) / figure_inches[1L], "nfc", "user")
  graphics::text(label_x, rows$y, rows$label, adj = c(1, 0.5), xpd = NA, cex = axis_cex)
  if(nrow(runs)){
    group_x <- graphics::grconvertX(0.1 / figure_inches[1L], "nfc", "user")
    graphics::text(group_x, runs$y, runs$label, adj = c(0, 0.5), xpd = NA, cex = group_cex, font = 2)
  }
  if(legend_at != "none"){
    legend_x <- graphics::grconvertX(0.5, "nfc", "user")
    inset <- 0.03 / figure_inches[2L]
    legend_y <- graphics::grconvertY(if(legend_at == "top") 1 - inset else inset, "nfc", "user")
    graphics::legend(
      legend_x, legend_y, legend = unname(data$legend_labels), ncol = legend_cols,
      fill = unname(data$palette), border = data$line_colour, bty = "n", xpd = NA,
      xjust = 0.5, yjust = if(legend_at == "top") 1 else 0, cex = legend_cex
    )
  }
}

.model_probability_plot_ggplot <- function(data, legend_at){

  # Budget display columns for the widest labels, swatches, and spacing.
  # The output device can be chosen after this plot is constructed.
  key_widths <- nchar(data$legend_labels, type = "width") + 4L
  legend_cols <- max(1L, min(4L, sum(cumsum(sort(key_widths, decreasing = TRUE)) <= 64L)))
  p <- ggplot2::ggplot(data$cells) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
                   ymax = .data$ymax, fill = .data$state)
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$xmin, xend = .data$xmin, y = .data$ymin, yend = .data$ymax,
                   colour = .data$colour),
      data = data$boundaries, linewidth = 0.25, linetype = "dotted", show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = 1, y = .data$y, yend = .data$y, colour = .data$colour),
      data = data$row_rules, linewidth = 0.25, show.legend = FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_manual(
      values = data$palette, breaks = names(data$palette), labels = unname(data$legend_labels),
      guide = ggplot2::guide_legend(
        ncol = legend_cols, order = 1L,
        override.aes = list(colour = data$line_colour, linewidth = 0.25)
      )
    ) +
    ggplot2::scale_x_continuous(breaks = data$ticks, labels = data$tick_labels) +
    ggplot2::scale_y_continuous(breaks = data$rows$y, labels = data$rows$label) +
    ggplot2::coord_cartesian(xlim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = data$axis_title, y = NULL, fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(hjust = 1, margin = ggplot2::margin(r = 6)),
      axis.ticks.x = ggplot2::element_line(colour = data$line_colour, linewidth = 0.25),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(2.5, "pt"),
      panel.border = ggplot2::element_rect(fill = NA, colour = data$line_colour, linewidth = 0.25),
      panel.grid = ggplot2::element_blank(),
      legend.position = legend_at,
      legend.key = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.9, "lines")
    )

  if(any(nzchar(data$runs$label))){
    group_names <- stats::setNames(data$runs$label, data$runs$run)
    p <- p + ggplot2::facet_grid(
      rows = ggplot2::vars(run = .data$run), scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::labeller(run = ggplot2::as_labeller(group_names))
    ) +
      ggplot2::theme(
        panel.spacing.y = grid::unit(0, "pt"),
        strip.background = ggplot2::element_blank(),
        strip.placement = "outside",
        strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", hjust = 0)
      )
  }
  attr(p, "n_rows") <- nrow(data$rows)
  p
}
