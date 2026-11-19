skip_if_not_test_profile("unit")

test_that("centered shrinkage is exchangeable in the declared coefficient basis", {
  dat <- data.frame(f = factor(rep(LETTERS[1:4], c(1, 2, 3, 4))))
  coefficient <- c(.7, -.3, 1.1)
  coefficient_names <- paste0("mu_f[", seq_along(coefficient), "]")
  permutation <- c(3L, 1L, 4L, 2L)
  normalization <- sqrt(1 - 1 / nlevels(dat$f))

  for(family in c("ridge", "lasso", "horseshoe", "regularized horseshoe")){
    bind <- function(contrast, scale){
      slab <- identical(family, "regularized horseshoe")
      group <- prior_regularized(if(slab) "horseshoe" else family, scale,
        terms = "f", contrast = contrast, slab_scale = if(slab) 3 * scale)
      JAGS_formula(~ f, "mu", dat, list(intercept = prior("point", list(0)), sh = group))$prior_list
    }
    for(contrast in c("meandif", "orthonormal")){
      priors <- bind(contrast, .6)
      spec <- priors$mu_f$regularization
      C <- attr(priors$mu_f, "factor_design")
      effect <- drop(C %*% coefficient)
      rotated <- drop(qr.solve(C, effect[permutation]))
      # A relabeling genuinely mixes the Helmert coordinates.
      expect_false(isTRUE(all.equal(sort(abs(rotated)), sort(abs(coefficient)))))
      before <- stats::setNames(coefficient, coefficient_names)
      if(length(spec$local_names)) before[spec$local_names] <- 1.3
      if(length(spec$variance_names)) before[spec$variance_names] <- 1.7
      if(length(spec$slab_name)) before[spec$slab_name] <- 2.3
      after <- before
      after[coefficient_names] <- rotated
      expect_equal(JAGS_marglik_priors(before, priors), JAGS_marglik_priors(after, priors),
                   tolerance = 1e-10)
    }
    draw_levels <- function(contrast, scale){
      p <- bind(contrast, scale)$mu_f
      set.seed(619)
      unname(rng(p, 100))
    }
    # Both contrast matrices operate on identical raw-coordinate laws.
    # Matching the level laws requires explicitly converting the supplied widths.
    meandif <- draw_levels("meandif", .6)
    expect_equal(draw_levels("orthonormal", .6), normalization * meandif, tolerance = 1e-10)
    expect_equal(draw_levels("orthonormal", .6 / normalization), meandif, tolerance = 1e-10)
    expect_lt(max(abs(rowSums(meandif))), 1e-10)
  }
})
