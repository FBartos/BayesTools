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
 * strtod() expects the decimal mark of the session's LC_NUMERIC, which R keeps
 * at "C" but a session or an embedding application can change; under a
 * decimal comma it stops at the point. The text always has a decimal point
 * (fits are refused under another mark, .bt_check_numeric_locale(), but
 * printed numbers and stored keys are read outside fits too), so it is read
 * under the "C" numeric locale: LC_NUMERIC is set to "C" for the reading and
 * restored.
 */
SEXP BayesTools_decimal_read(SEXP text)
{
  if(TYPEOF(text) != STRSXP){
    Rf_error("'text' must be a character vector.");
  }
  R_xlen_t n = XLENGTH(text);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *value = REAL(out);

  /* the copy is allocated before the locale changes; nothing after the
     change can raise an R error, so the session's locale is always restored */
  const char *current = setlocale(LC_NUMERIC, NULL);
  char *session = NULL;
  if(current != NULL && strcmp(current, "C") != 0){
    session = R_alloc(strlen(current) + 1, sizeof(char));
    strcpy(session, current);
    setlocale(LC_NUMERIC, "C");
  }
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
  if(session != NULL){
    setlocale(LC_NUMERIC, session);
  }
  UNPROTECT(1);
  return out;
}
