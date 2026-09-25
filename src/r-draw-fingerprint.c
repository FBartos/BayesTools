#include <math.h>

#include <Rinternals.h>
#include <R_ext/Arith.h>

/*
 * Value fingerprint of posterior draws, recorded in their metadata container
 * and recomputed when the metadata are read (R/draws-metadata.R).
 *
 * One pass over the values, without allocation, returns
 *   [0] the number of values,
 *   [1] the number of missing values (NA or NaN),
 *   [2] the sum of the observed values,
 *   [3] the sum of the observed values times their (1-based) positions,
 *   [4] the sum of the absolute observed values, and
 *   [5] the sum of the absolute products of [3],
 * where [4] and [5] scale the rounding-error bound of [2] and [3]. The sums
 * are accumulated in double precision in a fixed order: value i goes to
 * partial sum i mod 4, in increasing i, and the partial sums are combined as
 * (p0 + p1) + (p2 + p3).
 */

typedef struct {
  double sum, weighted, abs_sum, abs_weighted;
} bt_fingerprint_lane;

/* Adds value 'value' at 1-based position 'position' to a lane; a missing
 * value (NaN, which includes NA_real_) is only counted. */
#define BT_FINGERPRINT_ADD(lane, value, position, missing)  \
  do {                                                       \
    double bt_value_ = (value);                              \
    if(ISNAN(bt_value_)){                                    \
      (missing) += 1.0;                                      \
    }else{                                                   \
      double bt_product_ = bt_value_ * (position);           \
      (lane).sum          += bt_value_;                      \
      (lane).weighted     += bt_product_;                    \
      (lane).abs_sum      += fabs(bt_value_);                \
      (lane).abs_weighted += fabs(bt_product_);              \
    }                                                        \
  } while(0)

static double bt_as_double(int value)
{
  return value == NA_INTEGER ? NA_REAL : (double) value;
}

SEXP BayesTools_draw_fingerprint(SEXP x)
{
  bt_fingerprint_lane lane0 = {0.0, 0.0, 0.0, 0.0};
  bt_fingerprint_lane lane1 = {0.0, 0.0, 0.0, 0.0};
  bt_fingerprint_lane lane2 = {0.0, 0.0, 0.0, 0.0};
  bt_fingerprint_lane lane3 = {0.0, 0.0, 0.0, 0.0};
  double missing = 0.0;
  R_xlen_t n = XLENGTH(x);
  R_xlen_t blocks = n - n % 4;

  switch(TYPEOF(x)){
  case REALSXP: {
    const double *values = REAL_RO(x);
    R_xlen_t i = 0;
    for(; i < blocks; i += 4){
      BT_FINGERPRINT_ADD(lane0, values[i],     (double) (i + 1), missing);
      BT_FINGERPRINT_ADD(lane1, values[i + 1], (double) (i + 2), missing);
      BT_FINGERPRINT_ADD(lane2, values[i + 2], (double) (i + 3), missing);
      BT_FINGERPRINT_ADD(lane3, values[i + 3], (double) (i + 4), missing);
    }
    if(i < n)     BT_FINGERPRINT_ADD(lane0, values[i],     (double) (i + 1), missing);
    if(i + 1 < n) BT_FINGERPRINT_ADD(lane1, values[i + 1], (double) (i + 2), missing);
    if(i + 2 < n) BT_FINGERPRINT_ADD(lane2, values[i + 2], (double) (i + 3), missing);
    break;
  }
  case INTSXP:
  case LGLSXP: {
    const int *values = TYPEOF(x) == INTSXP ? INTEGER_RO(x) : LOGICAL_RO(x);
    R_xlen_t i = 0;
    for(; i < blocks; i += 4){
      BT_FINGERPRINT_ADD(lane0, bt_as_double(values[i]),     (double) (i + 1), missing);
      BT_FINGERPRINT_ADD(lane1, bt_as_double(values[i + 1]), (double) (i + 2), missing);
      BT_FINGERPRINT_ADD(lane2, bt_as_double(values[i + 2]), (double) (i + 3), missing);
      BT_FINGERPRINT_ADD(lane3, bt_as_double(values[i + 3]), (double) (i + 4), missing);
    }
    if(i < n)     BT_FINGERPRINT_ADD(lane0, bt_as_double(values[i]),     (double) (i + 1), missing);
    if(i + 1 < n) BT_FINGERPRINT_ADD(lane1, bt_as_double(values[i + 1]), (double) (i + 2), missing);
    if(i + 2 < n) BT_FINGERPRINT_ADD(lane2, bt_as_double(values[i + 2]), (double) (i + 3), missing);
    break;
  }
  default:
    Rf_error("Draw fingerprints require numeric or logical values.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, 6));
  double *result = REAL(out);
  result[0] = (double) n;
  result[1] = missing;
  result[2] = (lane0.sum + lane1.sum) + (lane2.sum + lane3.sum);
  result[3] = (lane0.weighted + lane1.weighted) + (lane2.weighted + lane3.weighted);
  result[4] = (lane0.abs_sum + lane1.abs_sum) + (lane2.abs_sum + lane3.abs_sum);
  result[5] = (lane0.abs_weighted + lane1.abs_weighted) +
    (lane2.abs_weighted + lane3.abs_weighted);
  UNPROTECT(1);
  return out;
}
