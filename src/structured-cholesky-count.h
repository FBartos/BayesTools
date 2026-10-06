#ifndef BAYESTOOLS_STRUCTURED_CHOLESKY_COUNT_H
#define BAYESTOOLS_STRUCTURED_CHOLESKY_COUNT_H

#include <cstddef>

// No multiplication is performed until its result is bounded by the limit.
inline bool bt_structured_cholesky_checked_count(std::size_t n_draws,
                                                std::size_t n_columns,
                                                std::size_t limit,
                                                std::size_t *result)
{
  if(n_draws == 0 || n_columns == 0 || result == NULL ||
     n_draws > limit / n_columns){
    return false;
  }
  const std::size_t partial = n_draws * n_columns;
  if(partial > limit / n_columns){
    return false;
  }
  *result = partial * n_columns;
  return true;
}

#endif
