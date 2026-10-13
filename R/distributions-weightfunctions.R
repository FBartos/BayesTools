# Step-plot vertex indices of a weightfunction with 'n_cuts' cut points
# (including the 0 and 1 end points): x indexes the cut points and y the bins.
.weightfunction_step_indices <- function(n_cuts){

  if(!is.numeric(n_cuts) || length(n_cuts) != 1L || !is.finite(n_cuts) ||
     n_cuts < 2 || n_cuts != as.integer(n_cuts)){
    stop("'n_cuts' must be an integer of at least 2.", call. = FALSE)
  }
  n_cuts <- as.integer(n_cuts)
  n_bins <- n_cuts - 1L
  if(n_cuts == 2L){
    return(list(x = c(1L, 2L), y = c(1L, 1L)))
  }

  list(
    x = c(1L, sort(rep(2L:(n_cuts - 1L), 2L)), n_cuts),
    y = sort(rep(seq_len(n_bins), 2L))
  )
}
# Plotted mean of a weight that has none (a heavy right tail on the log scale
# makes E[omega] infinite): the mean line is omitted for that bin, so NA breaks
# the line, and the y range comes from the finite values.
.weightfunction_plot_mean <- function(x_mean){

  x_mean[!is.finite(x_mean)] <- NA_real_
  x_mean
}
.weightfunction_plot_y_range <- function(x_mean, x_lCI, x_uCI){

  c(0, max(1, .weightfunction_plot_mean(x_mean), x_lCI, x_uCI, na.rm = TRUE))
}
