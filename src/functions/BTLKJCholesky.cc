#include "BTLKJCholesky.h"
#include "../lkj/BTLKJCore.h"

#include <cmath>
#include <vector>

namespace jags {
  namespace BayesTools {

    namespace {

      bool valid_K(double K_value)
      {
        return std::isfinite(K_value) && K_value >= 2.0 && std::floor(K_value) == K_value;
      }

      unsigned int as_K(double K_value)
      {
        return static_cast<unsigned int>(K_value);
      }

      bool check_lkj_args(std::vector<double const *> const &args,
                          std::vector<unsigned int> const &lengths)
      {
        if(lengths[1] != 1 || !valid_K(*args[1])){
          return false;
        }

        unsigned int K = as_K(*args[1]);
        return lengths[0] == bayestools::lkj::n_cpc(K) &&
          bayestools::lkj::valid_u(args[0], lengths[0]);
      }
    }

    BTLKJCholesky::BTLKJCholesky() : VectorFunction("bt_lkj_cholesky", 2)
    {
    }

    void BTLKJCholesky::evaluate(double *value,
                                 std::vector<double const *> const &args,
                                 std::vector<unsigned int> const &lengths) const
    {
      bayestools::lkj::fill_cholesky_from_u(value, args[0], as_K(*args[1]));
    }

    unsigned int BTLKJCholesky::length(std::vector<unsigned int> const &lengths,
                                       std::vector<double const *> const &args) const
    {
      if(lengths[1] != 1 || !valid_K(*args[1])){
        return 0;
      }
      unsigned int K = as_K(*args[1]);
      return K * K;
    }

    bool BTLKJCholesky::checkParameterLength(std::vector<unsigned int> const &lengths) const
    {
      return lengths[0] > 0 && lengths[1] == 1;
    }

    bool BTLKJCholesky::checkParameterValue(std::vector<double const *> const &args,
                                            std::vector<unsigned int> const &lengths) const
    {
      return check_lkj_args(args, lengths);
    }

    bool BTLKJCholesky::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return mask[1];
    }

    bool BTLKJCholesky::checkParameterFixed(std::vector<bool> const &mask) const
    {
      return mask[1];
    }

    bool BTLKJCholesky::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return false;
    }

    BTLKJCorr::BTLKJCorr() : VectorFunction("bt_lkj_corr", 2)
    {
    }

    void BTLKJCorr::evaluate(double *value,
                             std::vector<double const *> const &args,
                             std::vector<unsigned int> const &lengths) const
    {
      bayestools::lkj::fill_corr_from_u(value, args[0], as_K(*args[1]));
    }

    unsigned int BTLKJCorr::length(std::vector<unsigned int> const &lengths,
                                   std::vector<double const *> const &args) const
    {
      if(lengths[1] != 1 || !valid_K(*args[1])){
        return 0;
      }
      unsigned int K = as_K(*args[1]);
      return K * K;
    }

    bool BTLKJCorr::checkParameterLength(std::vector<unsigned int> const &lengths) const
    {
      return lengths[0] > 0 && lengths[1] == 1;
    }

    bool BTLKJCorr::checkParameterValue(std::vector<double const *> const &args,
                                        std::vector<unsigned int> const &lengths) const
    {
      return check_lkj_args(args, lengths);
    }

    bool BTLKJCorr::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return mask[1];
    }

    bool BTLKJCorr::checkParameterFixed(std::vector<bool> const &mask) const
    {
      return mask[1];
    }

    bool BTLKJCorr::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return false;
    }
  }
}
