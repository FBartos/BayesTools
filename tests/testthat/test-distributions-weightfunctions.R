context("Distributions - Weight functions")

test_that("Density function works", {

  ### based on the fact that marginals of cumulative Dirichlet distribution (beta distributions)
  expect_equal(mdone.sided(0.5, alpha = c(1, 1)), matrix(c(1, 0), nrow = 1, ncol = 2))

  ### visual verification
  test_dist <- function(x, samples, densities,i){
    hist(samples[,i], freq = FALSE, xlim = c(0, 1))
    lines(x, densities[,i])
  }
  x_seq     <- seq(0, 1, 0.01)


  set.seed(1)
  samples   <- t(apply(extraDistr::rdirichlet(10000, c(1, 1)), 1, cumsum))
  densities <- mdone.sided(x_seq, alpha = c(1, 1))

  expect_doppelganger("mdone.sided-1", function()test_dist(x_seq, samples, densities, 1))
  expect_true(all(densities[x_seq != 1,2] == 0) & is.infinite(densities[x_seq == 1,2]))


  set.seed(2)
  samples   <- t(apply(extraDistr::rdirichlet(10000, c(.3, 1, 2)), 1, cumsum))
  densities <- mdone.sided(x_seq, alpha = c(.3, 1, 2))

  expect_doppelganger("mdone.sided_2-1", function()test_dist(x_seq, samples, densities, 1))
  expect_doppelganger("mdone.sided_2-2", function()test_dist(x_seq, samples, densities, 2))
  expect_true(all(densities[x_seq != 1,3] == 0) & is.infinite(densities[x_seq == 1,3]))


  set.seed(3)
  samples   <- t(apply(extraDistr::rdirichlet(10000, c(4, 3, 5, 2)), 1, cumsum))
  densities <- mdone.sided(x_seq, alpha = c(4, 3, 5, 2))

  expect_doppelganger("mdone.sided_3-1", function()test_dist(x_seq, samples, densities, 1))
  expect_doppelganger("mdone.sided_3-2", function()test_dist(x_seq, samples, densities, 2))
  expect_doppelganger("mdone.sided_3-3", function()test_dist(x_seq, samples, densities, 3))
  expect_true(all(densities[x_seq != 1,4] == 0) & is.infinite(densities[x_seq == 1,4]))


  ### one-sided and two-sided should be equal
  densities1 <- mdone.sided(x_seq, alpha = c(1, 1))
  densities2 <- mdtwo.sided(x_seq, alpha = c(1, 1))
  expect_true(all.equal(densities1, densities2))

  densities1 <- mdone.sided(x_seq, alpha = c(3, 5, .4))
  densities2 <- mdtwo.sided(x_seq, alpha = c(3, 5, .4))
  expect_true(all.equal(densities1, densities2))


  ### the logarithm works
  expect_equal(mdone.sided(c(0, 1), alpha = c(1, 1), log = TRUE), matrix(c(0, 0, 0, Inf), ncol = 2, nrow = 2))

  densities1 <- mdone.sided(x_seq, alpha = c(.1, 11, 5, .4), log = TRUE)
  densities2 <- mdtwo.sided(x_seq, alpha = c(.1, 11, 5, .4), log = TRUE)
  expect_true(all.equal(densities1, densities2))


  ### the non-monotonic one-sided is not implemented since it requires a pdf of a product of beta distributed variables


  ### the fixed weight functions
  expect_equal(mdone.sided_fixed(c(.3, 0),    omega = c(.3, 1)), matrix(c(Inf, 0, 0, 0),         ncol = 2, nrow = 2))
  expect_equal(mdone.sided_fixed(c(0, .5, 1), omega = c(.5, 1)), matrix(c(0, Inf, 0, 0, 0, Inf), ncol = 2, nrow = 3))
  expect_equal(mdone.sided_fixed(c(0, .5, 1), omega = c(.5, 1), log = TRUE), log(matrix(c(0, Inf, 0, 0, 0, Inf), ncol = 2, nrow = 3)))
  expect_true(all.equal(mdone.sided_fixed(c(0, .5, 1), omega = c(.5, 1)), mdtwo.sided_fixed(c(0, .5, 1), omega = c(.5, 1))))
})

test_that("Random generator works", {

  ### visual verification that the samples corresponds to the marginals of cumulative Dirichlet distribution (beta distributions)
  test_dist <- function(alpha, beta, samples, i){
    hist(samples[,i], freq = FALSE, xlim = c(0, 1))
    lines(seq(0, 1, 0.01), stats::dbeta(seq(0, 1, 0.01), alpha, beta))
    return(invisible())
  }


  set.seed(1)
  samples   <- rone.sided(10000, alpha = c(1, 1))

  expect_doppelganger("rone.sided-1", function()test_dist(1, 1, samples, 1))
  expect_equal(samples[,ncol(samples)], rep(1, nrow(samples)))


  set.seed(2)
  samples   <- rone.sided(10000, alpha = c(.3, 1, 2))

  expect_doppelganger("rone.sided-2-1", function()test_dist(0.3, 3, samples, 1))
  expect_doppelganger("rone.sided-2-2", function()test_dist(1.3, 2, samples, 2))
  expect_equal(samples[,ncol(samples)], rep(1, nrow(samples)))


  set.seed(3)
  samples   <- rone.sided(10000, alpha = c(4, 3, 5, 2))

  expect_doppelganger("rone.sided-3-1", function()test_dist( 4, 10, samples, 1))
  expect_doppelganger("rone.sided-3-2", function()test_dist( 7,  7, samples, 2))
  expect_doppelganger("rone.sided-3-3", function()test_dist(12,  2, samples, 3))
  expect_equal(samples[,ncol(samples)], rep(1, nrow(samples)))


  ### one-sided and two-sided should be equal
  set.seed(4)
  samples1 <- rone.sided(10, alpha = c(1, 1))
  set.seed(4)
  samples2 <- rtwo.sided(10, alpha = c(1, 1))
  expect_true(all.equal(samples1, samples2))


  set.seed(5)
  samples1 <- rone.sided(10, alpha = c(3, 5, .4))
  set.seed(5)
  samples2 <- rtwo.sided(10, alpha = c(3, 5, .4))
  expect_true(all.equal(samples1, samples2))


  ### at least a visual check for the non-monotonic one-sided
  set.seed(6)
  samples   <- rone.sided(10000, alpha1 = c(1, 1), alpha2 = c(1, 1))
  expect_doppelganger("rone.sided-4-1", function()hist(samples[,1], xlim = c(0, 1), main = "one.sided(alpha1 = c(1, 1), alpha2 = c(1, 1))"))
  expect_doppelganger("rone.sided-4-2", function()hist(samples[,2], xlim = c(0, 1), main = "one.sided(alpha1 = c(1, 1), alpha2 = c(1, 1))"))
  expect_equal(samples[,ncol(samples)], rep(1, nrow(samples)))

  set.seed(7)
  samples   <- rone.sided(10000, alpha1 = c(1, 1, 1), alpha2 = c(1, 1))
  expect_doppelganger("rone.sided-5-1", function()hist(samples[,1], xlim = c(0, 1), main = "one.sided(alpha1 = c(1, 1, 1), alpha2 = c(1, 1))"))
  expect_doppelganger("rone.sided-5-2", function()hist(samples[,2], xlim = c(0, 1), main = "one.sided(alpha1 = c(1, 1, 1), alpha2 = c(1, 1))"))
  expect_doppelganger("rone.sided-5-3", function()hist(samples[,3], xlim = c(0, 1), main = "one.sided(alpha1 = c(1, 1, 1), alpha2 = c(1, 1))"))
  expect_equal(samples[,ncol(samples)], rep(1, nrow(samples)))


  ### the fixed weight functions
  expect_equal(rone.sided_fixed(3, omega = c(.3, 1)), cbind(rep(.3, 3), rep(1, 3)))
  expect_equal(rone.sided_fixed(5, omega = c(.5, 1, .8)), cbind(rep(.5, 5), rep(1, 5), rep(.8, 5)))
  expect_true(all.equal(rone.sided_fixed(5, omega = c(.5, 1)), rtwo.sided_fixed(5, omega = c(.5, 1))))
})

test_that("Distribution function works", {

  ### based on the fact that marginals of cumulative Dirichlet distribution (beta distributions)
  expect_equal(mpone.sided(0.5, alpha = c(1, 1)), matrix(c(0.5, 0), nrow = 1, ncol = 2))

  q_seq <- seq(0, 1, 0.10)

  expect_equal(mpone.sided(q_seq, alpha = c(1, 1)),
               cbind(
                 stats::pbeta(q_seq, 1, 1),
                 c(rep(0, length(q_seq) - 1), 1)
               ))
  expect_equal(mpone.sided(q_seq, alpha = c(.3, 1, 2)),
               cbind(
                 stats::pbeta(q_seq, 0.3, 3),
                 stats::pbeta(q_seq, 1.3, 2),
                 c(rep(0, length(q_seq) - 1), 1)
               ))
  expect_equal(mpone.sided(q_seq, alpha = c(4, 3, 5, 2)),
               cbind(
                 stats::pbeta(q_seq,  4, 10),
                 stats::pbeta(q_seq,  7,  7),
                 stats::pbeta(q_seq, 12,  2),
                 c(rep(0, length(q_seq) - 1), 1)
               ))


  ### one-sided and two-sided should be equal
  expect_equal(mqone.sided(q_seq, alpha = c(1, 1)),     mqtwo.sided(q_seq, alpha = c(1, 1)))
  expect_equal(mqone.sided(q_seq, alpha = c(3, 5, .4)), mqtwo.sided(q_seq, alpha = c(3, 5, .4)))


  ### the logarithm works
  expect_equal(mpone.sided(q_seq, alpha = c(1, 1), log.p = TRUE),
               cbind(
                 stats::pbeta(q_seq, 1, 1, log.p = TRUE),
                 log(c(rep(0, length(q_seq) - 1), 1))
               ))


  ### the lower-tail works
  expect_equal(mpone.sided(q_seq, alpha = c(1, 1), lower.tail = FALSE),
               cbind(
                 stats::pbeta(q_seq, 1, 1, lower.tail = FALSE),
                 1 - c(rep(0, length(q_seq) - 1), 1)
               ))


  ### the non-monotonic one-sided is not implemented since it requires a pdf of a product of beta distributed variables


  ### the fixed weight functions
  expect_equal(mpone.sided_fixed(c(0, .5, 1), omega = c(.5, 1)), matrix(c(0, 1, 1, 0, 0, 1), ncol = 2, nrow = 3))
  expect_equal(mpone.sided_fixed(c(0, .5, 1), omega = c(.5, 1), log.p = TRUE), log(matrix(c(0, 1, 1, 0, 0, 1), ncol = 2, nrow = 3)))
  expect_equal(mpone.sided_fixed(c(0, .5, 1), omega = c(.5, 1), lower.tail = FALSE), 1- matrix(c(0, 1, 1, 0, 0, 1), ncol = 2, nrow = 3))
  expect_equal(mpone.sided_fixed(c(0, .5, 1), omega = c(.5, 1)), mptwo.sided_fixed(c(0, .5, 1), omega = c(.5, 1)))
})

test_that("Quantile function works", {

  ### based on the fact that marginals of cumulative Dirichlet distribution (beta distributions)
  expect_equal(mqone.sided(0.5, alpha = c(1, 1)), matrix(c(0.5, 1), nrow = 1, ncol = 2))

  p_seq <- seq(0, 1, 0.10)

  expect_equal(mqone.sided(p_seq, alpha = c(1, 1)),
               cbind(
                 stats::qbeta(p_seq, 1, 1),
                 c(0, rep(1, length(p_seq) - 1))
               ))
  expect_equal(mqone.sided(p_seq, alpha = c(.3, 1, 2)),
               cbind(
                 stats::qbeta(p_seq, 0.3, 3),
                 stats::qbeta(p_seq, 1.3, 2),
                 c(0, rep(1, length(p_seq) - 1))
               ))
  expect_equal(mqone.sided(p_seq, alpha = c(4, 3, 5, 2)),
               cbind(
                 stats::qbeta(p_seq,  4, 10),
                 stats::qbeta(p_seq,  7,  7),
                 stats::qbeta(p_seq, 12,  2),
                 c(0, rep(1, length(p_seq) - 1))
               ))


  ### one-sided and two-sided should be equal
  expect_equal(mqone.sided(p_seq, alpha = c(1, 1)),     mqtwo.sided(p_seq, alpha = c(1, 1)))
  expect_equal(mqone.sided(p_seq, alpha = c(3, 5, .4)), mqtwo.sided(p_seq, alpha = c(3, 5, .4)))


  ### the logarithm works
  expect_equal(mqone.sided(log(p_seq), alpha = c(1, 1), log.p = TRUE),
               cbind(
                 stats::qbeta(log(p_seq), 1, 1, log.p = TRUE),
                 c(0, rep(1, length(p_seq) - 1))
               ))


  ### the lower-tail works
  expect_equal(mqone.sided(p_seq, alpha = c(1, 1), lower.tail = FALSE),
               cbind(
                 stats::qbeta(p_seq, 1, 1, lower.tail = FALSE),
                 rep(1, length(p_seq))
               ))


  ### the non-monotonic one-sided is not implemented since it requires a pdf of a product of beta distributed variables


  ### the fixed weight functions
  expect_equal(mqone.sided_fixed(c(0, .5, 1), omega = c(.5, 1)), matrix(c(0, .5, .5, 0, 1, 1), ncol = 2, nrow = 3))
  expect_equal(mqone.sided_fixed(log(c(0, .5, 1)), omega = c(.5, 1), log.p = TRUE), matrix(c(0, .5, .5, 0, 1, 1), ncol = 2, nrow = 3))
  expect_equal(mqone.sided_fixed(c(0, .5, 1), omega = c(.5, 1), lower.tail = FALSE), matrix(c(1, .5, .5, 1, 1, 1), ncol = 2, nrow = 3))
  expect_equal(mqone.sided_fixed(c(0, .5, 1), omega = c(.5, 1)), mqtwo.sided_fixed(c(0, .5, 1), omega = c(.5, 1)))
})
