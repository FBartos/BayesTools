skip_if_not_test_profile("unit")

# The orthonormal and mean-difference contrasts are the normalised Helmert
# basis: column j has -1 / sqrt(j (j + 1)) in rows 1..j and j / sqrt(j (j + 1))
# in row j + 1. The expectations below spell the entries as the same IEEE
# expressions, so identical() states bit-for-bit platform independence.

test_that("orthonormal contrasts are the closed-form Helmert basis for 2, 3 and 4 levels", {

  expected_2 <- cbind(c(-1, 1) / sqrt(2))
  expected_3 <- cbind(c(-1, 1, 0) / sqrt(2), c(-1, -1, 2) / sqrt(6))
  expected_4 <- cbind(
    c(-1, 1, 0, 0) / sqrt(2),
    c(-1, -1, 2, 0) / sqrt(6),
    c(-1, -1, -1, 3) / sqrt(12)
  )

  expect_identical(contr.orthonormal(2), expected_2)
  expect_identical(contr.orthonormal(3), expected_3)
  expect_identical(contr.orthonormal(4), expected_4)

  # the entries are single correctly rounded operations, not decimal literals
  expect_identical(contr.orthonormal(3)[1, 1], -1 / sqrt(2))
  expect_identical(contr.orthonormal(3)[3, 2], 2 / sqrt(6))
  expect_identical(contr.orthonormal(4)[4, 3], 3 / sqrt(12))

  # the entries below the sub-diagonal of the Helmert layout are exact zeros
  below <- outer(1:7, 1:6, function(i, j) i > j + 1)
  expect_true(all(contr.orthonormal(7)[below] == 0))
})

test_that("orthonormal contrasts have orthonormal, centred columns and imply the centring matrix", {

  for(a in 2:10){
    basis <- contr.orthonormal(a)

    expect_identical(dim(basis), c(a, a - 1L))
    expect_null(dimnames(basis))
    # absolute errors: the closed form is orthonormal and centred to rounding
    expect_lt(max(abs(crossprod(basis) - diag(a - 1L))), 1e-14)
    expect_lt(max(abs(colSums(basis))), 1e-14)
    expect_lt(max(abs(basis %*% t(basis) - (diag(a) - 1 / a))), 1e-14)
  }
})

test_that("mean-difference contrasts rescale the orthonormal basis", {

  for(a in 2:10){
    expect_identical(contr.meandif(a), contr.orthonormal(a) / sqrt(1 - 1 / a))
    expect_identical(dim(contr.meandif(a)), c(a, a - 1L))
    expect_null(dimnames(contr.meandif(a)))

    # B B' = (I - J / a) / (1 - 1 / a): every level has unit prior variance
    # relative to the grand mean under a unit-scale spherical prior
    meandif <- contr.meandif(a)
    expect_lt(max(abs(meandif %*% t(meandif) - (diag(a) - 1 / a) / (1 - 1 / a))), 1e-14)
  }

  # the second level of a two-level factor has coefficient 1 (to rounding)
  expect_lt(abs(contr.meandif(2)[2, 1] - 1), 1e-14)
})

test_that("contrasts = FALSE returns the identity and level names give the level count", {

  for(a in 2:5){
    expect_identical(contr.orthonormal(a, contrasts = FALSE), diag(a))
    expect_identical(contr.meandif(a, contrasts = FALSE), diag(a))

    levels <- letters[seq_len(a)]
    expect_identical(contr.orthonormal(levels), contr.orthonormal(a))
    expect_identical(contr.meandif(levels), contr.meandif(a))
    expect_identical(contr.orthonormal(levels, contrasts = FALSE), diag(a))
    expect_identical(contr.meandif(levels, contrasts = FALSE), diag(a))
  }

  expect_identical(contr.orthonormal(c(1, 2, 3)), contr.orthonormal(3))
  expect_identical(contr.meandif(factor(c("a", "b", "c", "d"))), contr.meandif(4))
})

test_that("contrast functions reject an insufficient or invalid level count", {

  for(contrast in list(contr.orthonormal, contr.meandif)){
    expect_error(contrast(1), "Not enough degrees of freedom to define contrasts.", fixed = TRUE)
    expect_error(contrast(0), "Not enough degrees of freedom to define contrasts.", fixed = TRUE)
    expect_error(contrast(3, contrasts = NA))
  }
})
