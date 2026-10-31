#include "BTLKJCholesky.h"
#include "../lkj/BTLKJCore.h"

#include <cmath>
#include <limits>
#include <vector>

namespace jags {
  namespace BayesTools {

    namespace {

      unsigned int as_K(double K_value)
      {
        return static_cast<unsigned int>(K_value);
      }

      bool valid_geometry(std::vector<double const *> const &args,
                          std::vector<unsigned int> const &lengths)
      {
        if(lengths[1] != 1){
          return false;
        }
        const double K_value = *args[1];
        if(!std::isfinite(K_value) || K_value < 2.0 || std::floor(K_value) != K_value ||
           K_value * K_value > static_cast<double>(std::numeric_limits<unsigned int>::max())){
          return false;
        }

        return lengths[0] == bayestools::lkj::n_cpc(as_K(K_value));
      }

      bool finite_transform_u(double const *u, unsigned int length)
      {
        // JAGS can evaluate descendants before rejecting an out-of-support
        // stochastic proposal. The transform therefore accepts every finite
        // proposal; fill_cholesky_from_u() supplies the finite extension.
        for(unsigned int i = 0; i < length; ++i){
          if(!std::isfinite(u[i])){
            return false;
          }
        }

        return true;
      }

      bool check_lkj_args(std::vector<double const *> const &args,
                          std::vector<unsigned int> const &lengths)
      {
        return valid_geometry(args, lengths) &&
          finite_transform_u(args[0], lengths[0]);
      }
    }

    BTLKJCholesky::BTLKJCholesky() : VectorFunction("bt_lkj_cholesky", 2)
    {
    }

    void BTLKJCholesky::evaluate(double *value,
                                 std::vector<double const *> const &args,
                                 std::vector<unsigned int> const &lengths) const
    {
      if(!valid_geometry(args, lengths)){
        return;
      }
      bayestools::lkj::fill_cholesky_from_u(value, args[0], as_K(*args[1]));
    }

    unsigned int BTLKJCholesky::length(std::vector<unsigned int> const &lengths,
                                       std::vector<double const *> const &args) const
    {
      if(!valid_geometry(args, lengths)){
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

    // bt_lkj_cholesky_ordered(u, K, position): the lower Cholesky factor in
    // design order of the vine-order correlation matrix of u, where design
    // column d is vine variable position[d] (a fixed permutation of 1, ..., K).
    BTLKJCholeskyOrdered::BTLKJCholeskyOrdered() : VectorFunction("bt_lkj_cholesky_ordered", 3)
    {
    }

    void BTLKJCholeskyOrdered::evaluate(double *value,
                                        std::vector<double const *> const &args,
                                        std::vector<unsigned int> const &lengths) const
    {
      if(!valid_geometry(args, lengths) || lengths[2] != as_K(*args[1]) ||
         !bayestools::lkj::valid_position(args[2], lengths[2])){
        return;
      }
      bayestools::lkj::fill_cholesky_from_u_ordered(value, args[0], as_K(*args[1]), args[2]);
    }

    unsigned int BTLKJCholeskyOrdered::length(std::vector<unsigned int> const &lengths,
                                              std::vector<double const *> const &args) const
    {
      if(!valid_geometry(args, lengths)){
        return 0;
      }
      unsigned int K = as_K(*args[1]);
      return K * K;
    }

    bool BTLKJCholeskyOrdered::checkParameterLength(std::vector<unsigned int> const &lengths) const
    {
      return lengths[0] > 0 && lengths[1] == 1 && lengths[2] > 0;
    }

    bool BTLKJCholeskyOrdered::checkParameterValue(std::vector<double const *> const &args,
                                                   std::vector<unsigned int> const &lengths) const
    {
      return check_lkj_args(args, lengths) && lengths[2] == as_K(*args[1]) &&
        bayestools::lkj::valid_position(args[2], lengths[2]);
    }

    bool BTLKJCholeskyOrdered::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return mask[1] && mask[2];
    }

    bool BTLKJCholeskyOrdered::checkParameterFixed(std::vector<bool> const &mask) const
    {
      return mask[1] && mask[2];
    }

    bool BTLKJCholeskyOrdered::isDiscreteValued(std::vector<bool> const &mask) const
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
      if(!valid_geometry(args, lengths)){
        return;
      }
      bayestools::lkj::fill_corr_from_u(value, args[0], as_K(*args[1]));
    }

    unsigned int BTLKJCorr::length(std::vector<unsigned int> const &lengths,
                                   std::vector<double const *> const &args) const
    {
      if(!valid_geometry(args, lengths)){
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
