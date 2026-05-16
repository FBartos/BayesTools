#ifndef BTLKJCHOLESKY_H_
#define BTLKJCHOLESKY_H_

#include <function/VectorFunction.h>

namespace jags {
  namespace BayesTools {

    class BTLKJCholesky : public VectorFunction
    {
    public:
      BTLKJCholesky();
      void evaluate(double *value,
                    std::vector<double const *> const &args,
                    std::vector<unsigned int> const &lengths) const;
      unsigned int length(std::vector<unsigned int> const &lengths,
                          std::vector<double const *> const &args) const;
      bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
      bool checkParameterValue(std::vector<double const *> const &args,
                               std::vector<unsigned int> const &lengths) const;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const;
      bool checkParameterFixed(std::vector<bool> const &mask) const;
      bool isDiscreteValued(std::vector<bool> const &mask) const;
    };

    class BTLKJCorr : public VectorFunction
    {
    public:
      BTLKJCorr();
      void evaluate(double *value,
                    std::vector<double const *> const &args,
                    std::vector<unsigned int> const &lengths) const;
      unsigned int length(std::vector<unsigned int> const &lengths,
                          std::vector<double const *> const &args) const;
      bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
      bool checkParameterValue(std::vector<double const *> const &args,
                               std::vector<unsigned int> const &lengths) const;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const;
      bool checkParameterFixed(std::vector<bool> const &mask) const;
      bool isDiscreteValued(std::vector<bool> const &mask) const;
    };
  }
}

#endif
