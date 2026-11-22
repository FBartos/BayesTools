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
                    std::vector<unsigned long> const &lengths) const override;
      unsigned long length(std::vector<unsigned long> const &lengths,
                          std::vector<double const *> const &args) const override;
      bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
      bool checkParameterValue(std::vector<double const *> const &args,
                               std::vector<unsigned long> const &lengths) const override;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const override;
      bool checkParameterFixed(std::vector<bool> const &mask) const override;
      bool isDiscreteValued(std::vector<bool> const &mask) const override;
    };

    class BTLKJCholeskyOrdered : public VectorFunction
    {
    public:
      BTLKJCholeskyOrdered();
      void evaluate(double *value,
                    std::vector<double const *> const &args,
                    std::vector<unsigned long> const &lengths) const override;
      unsigned long length(std::vector<unsigned long> const &lengths,
                          std::vector<double const *> const &args) const override;
      bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
      bool checkParameterValue(std::vector<double const *> const &args,
                               std::vector<unsigned long> const &lengths) const override;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const override;
      bool checkParameterFixed(std::vector<bool> const &mask) const override;
      bool isDiscreteValued(std::vector<bool> const &mask) const override;
    };

    class BTLKJCorr : public VectorFunction
    {
    public:
      BTLKJCorr();
      void evaluate(double *value,
                    std::vector<double const *> const &args,
                    std::vector<unsigned long> const &lengths) const override;
      unsigned long length(std::vector<unsigned long> const &lengths,
                          std::vector<double const *> const &args) const override;
      bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
      bool checkParameterValue(std::vector<double const *> const &args,
                               std::vector<unsigned long> const &lengths) const override;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const override;
      bool checkParameterFixed(std::vector<bool> const &mask) const override;
      bool isDiscreteValued(std::vector<bool> const &mask) const override;
    };
  }
}

#endif
