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
