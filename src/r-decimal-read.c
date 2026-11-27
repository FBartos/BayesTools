#include <locale.h>
#include <stdlib.h>
#include <string.h>

#include <Rinternals.h>

/*
 * The doubles the C library's strtod() reads from decimal text: JAGS reads the
 * numeric constants of a model with it, and it rounds correctly on the
 * platforms R supports. R reads decimal text (parse(), as.numeric()) with its
 * own algorithm, which accumulates the digits and the power of ten in long
 * double; where long double is double (Apple silicon), a 16- or 17-digit
 * constant can read as a neighbouring double. Text that is not entirely one
 * number is NA.
 *
 * The text always has a decimal point. strtod() expects the decimal mark of
 * the session's LC_NUMERIC, which R keeps at "C" but a user can change: the
 * point is then replaced by that mark before reading, and text that already
 * holds the mark is not a number. A mark of more than one byte cannot be
 * substituted in place; such text reads as NA and the R callers fall back.
 */
SEXP BayesTools_decimal_read(SEXP text)
{
  if(TYPEOF(text) != STRSXP){
    Rf_error("'text' must be a character vector.");
  }
  const struct lconv *convention = localeconv();
  const char *mark = convention != NULL ? convention->decimal_point : NULL;
  const int point_is_mark = mark == NULL || mark[0] == '\0' || (mark[0] == '.' && mark[1] == '\0');
  const int single_byte_mark = !point_is_mark && mark[1] == '\0';

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
    if(point_is_mark){
      double parsed = strtod(start, &end);
      value[i] = (end == start || *end != '\0') ? NA_REAL : parsed;
      continue;
    }
    if(!single_byte_mark || strchr(start, mark[0]) != NULL){
      value[i] = NA_REAL;
      continue;
    }
    size_t length = strlen(start);
    char *copy = R_alloc(length + 1, sizeof(char));
    memcpy(copy, start, length + 1);
    for(size_t k = 0; k < length; ++k){
      if(copy[k] == '.') copy[k] = mark[0];
    }
    double parsed = strtod(copy, &end);
    value[i] = (end == copy || *end != '\0') ? NA_REAL : parsed;
  }
  UNPROTECT(1);
  return out;
}
