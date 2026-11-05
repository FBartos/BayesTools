#include <stdlib.h>

#include <Rinternals.h>

/*
 * The doubles the C library's strtod() reads from decimal text: JAGS reads the
 * numeric constants of a model with it, and it rounds correctly on the
 * platforms R supports. R reads decimal text (parse(), as.numeric()) with its
 * own algorithm, which accumulates the digits and the power of ten in long
 * double; where long double is double (Apple silicon), a 16- or 17-digit
 * constant can read as a neighbouring double. Text that is not entirely one
 * number is NA.
 */
SEXP BayesTools_decimal_read(SEXP text)
{
  if(TYPEOF(text) != STRSXP){
    Rf_error("'text' must be a character vector.");
  }
  R_xlen_t n = XLENGTH(text);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *value = REAL(out);
  for(R_xlen_t i = 0; i < n; ++i){
    SEXP element = STRING_ELT(text, i);
    if(element == NA_STRING){
      value[i] = NA_REAL;
      continue;
    }
    const char *start = CHAR(element);
    char *end = NULL;
    double parsed = strtod(start, &end);
    value[i] = (end == start || *end != '\0') ? NA_REAL : parsed;
  }
  UNPROTECT(1);
  return out;
}
