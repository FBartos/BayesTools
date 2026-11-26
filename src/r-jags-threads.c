#include <Rinternals.h>

#ifdef _WIN32
#include <windows.h>
#endif

/*
 * The OpenMP thread count of the JAGS library loaded in this process.
 *
 * The Windows build of JAGS 5 (libjags-5.dll) links an OpenBLAS and an OpenMP
 * runtime in statically. With the default thread count, every larger BLAS or
 * LAPACK call of a multivariate node (dmnorm, dwish, ...) starts a team of
 * threads, which costs milliseconds per call on a machine with many cores and
 * makes such models more than a thousand times slower than on one thread; the
 * draws were identical in every comparison. Neither OPENBLAS_NUM_THREADS nor
 * openblas_set_num_threads() reaches that copy, but the omp_set_num_threads()
 * the library exports does. That runtime is JAGS's own: R and other packages
 * keep theirs. The count applies to the calling thread, on which JAGS runs
 * (BayesTools compiles its models with 'n.threads = 1').
 *
 * Windows only: elsewhere JAGS uses the system BLAS, whose OpenMP runtime, if
 * any, is shared by the whole process and is left alone.
 *
 * 'libraries' names the candidate JAGS libraries in order; the first one that
 * is loaded is used. 'threads' is one integer: with NA the count is only read,
 * otherwise it is set to 'threads'. Returns the count before the call, or NA
 * when no candidate is loaded, when the loaded one lacks either thread
 * control, and on platforms other than Windows.
 */
SEXP BayesTools_jags_threads(SEXP libraries, SEXP threads)
{
  if(TYPEOF(libraries) != STRSXP){
    Rf_error("'libraries' must be a character vector.");
  }
  if(TYPEOF(threads) != INTSXP || XLENGTH(threads) != 1){
    Rf_error("'threads' must be one integer.");
  }
  int value = INTEGER(threads)[0];
  if(value != NA_INTEGER && value < 1){
    Rf_error("'threads' must be NA or a positive integer.");
  }
#ifdef _WIN32
  typedef int (*get_threads_fn)(void);
  typedef void (*set_threads_fn)(int);
  for(R_xlen_t i = 0; i < XLENGTH(libraries); ++i){
    SEXP name = STRING_ELT(libraries, i);
    if(name == NA_STRING) continue;
    HMODULE library = GetModuleHandleA(CHAR(name));
    if(library == NULL) continue;
    /* the cast through void (*)(void) keeps GCC's -Wcast-function-type quiet */
    get_threads_fn get_threads = (get_threads_fn) (void (*)(void))
      GetProcAddress(library, "omp_get_max_threads");
    set_threads_fn set_threads = (set_threads_fn) (void (*)(void))
      GetProcAddress(library, "omp_set_num_threads");
    if(get_threads == NULL || set_threads == NULL) return Rf_ScalarInteger(NA_INTEGER);
    int previous = get_threads();
    if(value != NA_INTEGER) set_threads(value);
    return Rf_ScalarInteger(previous);
  }
#endif
  return Rf_ScalarInteger(NA_INTEGER);
}
