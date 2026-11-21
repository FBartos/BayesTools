skip_if_not_test_profile(c("unit", "visual"))

.model_probability_plot_fixture <- function(){

  probabilities <- c(model_3 = 0.25, model_1 = 0.5, model_2 = 0.25, model_zero = 0)
  cells <- data.frame(
    model = rep(rev(names(probabilities)), 4L),
    row = rep(c("tau", "mu", "loc", "bias"), each = 4L),
    label = rep(c("Heterogeneity", "Effect", "Dose", "Publication bias"), each = 4L),
    group = rep(c("Components", "Components", "Location", "Components"), each = 4L),
    state = c("excluded", "included", "excluded", "included",
              "excluded", "negative", "positive", "positive",
              "excluded", "positive", "excluded", "positive",
              "negative", "adjustment", "excluded", "adjustment")
  )
  list(probabilities = probabilities, cells = cells, col = c(adjustment = "#6C4590"))
}

test_that("model strips preserve posterior masses and declared tie and row order", {
  skip_if_not_test_profile("unit")
  fixture <- .model_probability_plot_fixture()
  probabilities <- fixture$probabilities
  probabilities["model_1"] <- 0.5 - 1e-10
  data <- .model_probability_plot_data(
    probabilities, fixture$cells, TRUE, c(loc = "Dose (mg)"), TRUE,
    fixture$col, c(adjustment = "Selection model")
  )
  first_row <- data$cells[data$cells$row == "tau", ]
  models <- c("model_1", "model_3", "model_2")
  ends <- c(0, cumsum(unname(probabilities[models])))

  expect_identical(first_row$model, models)
  expect_identical(first_row$xmin, ends[1:3])
  expect_identical(first_row$xmax, ends[2:4])
  expect_identical(data$boundaries$xmin, rep(ends[2:3], 4L))
  expect_identical(data$rows$row, c("tau", "mu", "loc", "bias"))
  expect_identical(data$rows$label, c("Heterogeneity", "Effect", "Dose (mg)", "Publication bias"))
  expect_identical(data$runs$label, c("Components", "Location", "Components"))
  expect_identical(data$runs$y, c(3.5, 2, 1))
  expect_identical(unname(data$palette["adjustment"]), "#6C4590")
  expect_identical(unname(data$legend_labels["adjustment"]), "Selection model")

  one_row <- fixture$cells[fixture$cells$row == "bias", ]
  ascending <- .model_probability_plot_data(
    fixture$probabilities, one_row, FALSE, NULL, FALSE, fixture$col, NULL
  )
  expect_identical(ascending$cells$model, c("model_3", "model_2", "model_1"))
  expect_identical(ascending$cells$xmax, c(0.25, 0.5, 1))
  expect_identical(names(ascending$palette), c("excluded", "adjustment"))
})

test_that("the public model probability plot dispatches and preserves row count", {
  skip_if_not_test_profile("unit")
  fixture <- .model_probability_plot_fixture()
  cells <- fixture$cells[fixture$cells$row != "mu", ]
  p <- plot_model_probabilities(
    fixture$probabilities, cells, group_labels = TRUE, col = fixture$col,
    legend = TRUE, legend_position = "top", plot_type = "ggplot"
  )
  expect_s3_class(p, "ggplot")
  expect_identical(attr(p, "n_rows", exact = TRUE), 3L)
  expect_identical(p$theme$legend.position, "top")
  built <- ggplot2::ggplot_build(p)
  expect_identical(nrow(built$layout$layout), 3L)
  expect_identical(
    built$layout$panel_params[[1L]]$x$get_labels(),
    c("0", "0.25", "0.50", "0.75", "1")
  )
  expect_no_warning(ggplot2::ggplotGrob(p))

  one_cell <- data.frame(model = "m", row = "tau", label = "Heterogeneity", group = "", state = "included")
  p <- plot_model_probabilities(c(m = 1), one_cell, legend = TRUE, plot_type = "ggplot")
  expect_identical(attr(p, "n_rows", exact = TRUE), 1L)
  expect_no_warning(ggplot2::ggplotGrob(p))

  output <- tempfile(fileext = ".pdf")
  withr::defer(unlink(output))
  withr::local_pdf(output)
  margin <- graphics::par("mai")
  base <- withVisible(plot_model_probabilities(
    fixture$probabilities, cells, group_labels = TRUE, col = fixture$col, legend = TRUE
  ))
  expect_identical(base, list(value = NULL, visible = FALSE))
  expect_identical(graphics::par("mai"), margin)
})

test_that("states wholly beyond the probability axis are absent from the legend", {
  skip_if_not_test_profile("unit")
  probabilities <- c(a = 0.6, b = 0.40000001, c = 1e-9)
  cells <- data.frame(
    model = names(probabilities), row = "effect", label = "Effect", group = "",
    state = c("positive", "negative", "included")
  )
  p <- plot_model_probabilities(probabilities, cells, legend = TRUE, plot_type = "ggplot")
  expect_identical(p$scales$get_scales("fill")$breaks, c("positive", "negative"))
  expect_identical(p$data$xmax, cumsum(unname(probabilities[1:2])))
})

test_that("ggplot legend wraps long model states within the panel width", {
  skip_if_not_test_profile("unit")
  probabilities <- c(a = 0.4, b = 0.3, c = 0.2, d = 0.1)
  cells <- data.frame(
    model = names(probabilities), row = "bias", label = "Publication bias", group = "",
    state = c("selection", "symmetric", "positive", "negative")
  )
  col <- c(selection = "cornflowerblue", symmetric = "darkgreen")
  state_labels <- c(
    selection = "Selection model (four bands)", symmetric = "Symmetric selection model",
    positive = "Positive small-study effect", negative = "Negative small-study effect"
  )
  withr::local_pdf(file = NULL, width = 7, height = 3)
  for(side in c("bottom", "top")){
    p <- plot_model_probabilities(
      probabilities, cells, col = col, state_labels = state_labels,
      legend = TRUE, legend_position = side, plot_type = "ggplot"
    )
    plot_grid <- ggplot2::ggplotGrob(p)
    guide <- match(paste0("guide-box-", side), plot_grid$layout$name)
    if(is.na(guide)){
      guide <- match("guide-box", plot_grid$layout$name)
    }
    guide_width <- grid::convertWidth(grid::grobWidth(plot_grid$grobs[[guide]]), "inches", valueOnly = TRUE)
    panel <- plot_grid$layout[plot_grid$layout$name == "panel", ]
    margins <- grid::convertWidth(sum(plot_grid$widths[-seq.int(panel$l, panel$r)]), "inches", valueOnly = TRUE)
    expect_lt(guide_width, grDevices::dev.size("in")[1L] - margins)
  }
})

test_that("the public model probability plot rejects invalid inputs", {
  skip_if_not_test_profile("unit")
  fixture <- .model_probability_plot_fixture()
  probabilities <- fixture$probabilities
  cells <- fixture$cells
  expect_error(plot_model_probabilities(probabilities * 0.8, cells),
               "The 'model_probabilities' argument must sum to one.", fixed = TRUE)
  expect_error(plot_model_probabilities(unname(probabilities), cells),
               "The 'model_probabilities' argument must have unique, nonempty names.", fixed = TRUE)
  expect_error(plot_model_probabilities(probabilities, cells[-1L, ]),
               "The 'cells' argument must contain exactly one cell for each model and row.", fixed = TRUE)
  expect_error(plot_model_probabilities(probabilities, cells[c(2L, 2:nrow(cells)), ]),
               "The 'cells' argument must contain exactly one cell for each model and row.", fixed = TRUE)
  cells$label[2L] <- "Inconsistent"
  expect_error(plot_model_probabilities(probabilities, cells),
               "The 'label' and 'group' in 'cells' must be constant within each row.", fixed = TRUE)
  expect_error(plot_model_probabilities(probabilities, fixture$cells),
               "Supply a colour in 'col' for each state without a default: 'adjustment'.", fixed = TRUE)
  expect_error(plot_model_probabilities(probabilities, fixture$cells, col = c(adjustment = "invalid")),
               "The 'col' argument must contain valid R colours.", fixed = TRUE)
  expect_error(plot_model_probabilities(probabilities, fixture$cells, labels = c(eta = "Wrong ID"),
                                        col = fixture$col),
               "The 'labels' argument must contain nonempty labels keyed by displayed row IDs.", fixed = TRUE)
})

test_that("the model probability matrix renders in both graphics backends", {
  skip_if_not_test_profile("visual")
  skip_if_not_installed("vdiffr")
  fixture <- .model_probability_plot_fixture()

  vdiffr::expect_doppelganger("model probabilities ggplot", plot_model_probabilities(
    fixture$probabilities, fixture$cells, labels = c(loc = "Dose (mg)"), group_labels = TRUE,
    col = fixture$col, state_labels = c(adjustment = "Selection model"),
    legend = TRUE, legend_position = "top", plot_type = "ggplot"
  ))
  vdiffr::expect_doppelganger("model probabilities base", function(){
    plot_model_probabilities(
      fixture$probabilities, fixture$cells, labels = c(loc = "Dose (mg)"), group_labels = TRUE,
      col = fixture$col, state_labels = c(adjustment = "Selection model"), legend = TRUE
    )
  })
})
