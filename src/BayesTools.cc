#include <module/Module.h>
#include <version.h>

static_assert(JAGS_MAJOR >= 5, "BayesTools requires JAGS 5.0.0 or later");

#include "distributions/DBTLKJCPC.h"
#include "distributions/DBTInvGamma.h"
#include "distributions/DBTMoment.h"
#include "distributions/DBTInvMoment.h"
#include "functions/BTLKJCholesky.h"
#include "distributions/DBTSamplerMarker.h"
#include "samplers/BTGatedSlice.h"
#include "samplers/BTFactorSlice.h"
#include "samplers/BTFiniteReadout.h"

#include <cstdio>
#include <exception>

#include <Rinternals.h>

namespace jags {
  namespace BayesTools {

    class BayesToolsModule : public Module
    {
      FiniteReadout _finite_readout;
    public:
      BayesToolsModule();
      ~BayesToolsModule();
      void prepare_finite_readout() { _finite_readout.prepare(); }
    };

    BayesToolsModule::BayesToolsModule() : Module("BayesTools"), _finite_readout(*this)
    {
      insert(new DBTLKJCPC);
      insert(new DBTInvGamma);
      insert(new DBTMoment);
      insert(new DBTInvMoment);
      insert(new BTLKJCholesky);
      insert(new BTLKJCholeskyOrdered);
      insert(new BTLKJCorr);
      insert(new DBTSamplerMarker(gate_marker_name));
      insert(new DBTSamplerMarker(block_marker_name));
      // Process-lifetime sampler factories that this destructor never deletes.
      // JAGS tries a factory inserted later first, so the block factory is tried
      // before the gated one; a gated node stays out of blocks because the block
      // factory skips every node with a gate marker.
      insert(gated_slice_factory());
      insert(factor_slice_factory());
    }

    BayesToolsModule::~BayesToolsModule()
    {
      // Normal process shutdown does not invoke the R namespace unload hook.
      unload();

      std::vector<Function*> const &fvec = functions();
      for(unsigned int i = 0; i < fvec.size(); ++i){
        delete fvec[i];
      }
      std::vector<Distribution*> const &dvec = distributions();
      for(unsigned int i = 0; i < dvec.size(); ++i){
        delete dvec[i];
      }
    }
  }
}

jags::BayesTools::BayesToolsModule _BayesTools_module;

extern "C" SEXP BayesTools_prepare_finite_readout()
{
  // Rf_error() jumps out of the function, so it is called only after the
  // handler has finished and the exception is destroyed.
  bool failed = false;
  char message[512] = "";
  try {
    _BayesTools_module.prepare_finite_readout();
  } catch (std::exception const &error) {
    failed = true;
    std::snprintf(message, sizeof message, "%s", error.what());
  } catch (...) {
    failed = true;
    std::snprintf(message, sizeof message, "%s", "unknown exception");
  }
  if (failed) {
    Rf_error("Could not prepare the BayesTools finite readout: %s", message);
  }
  return R_NilValue;
}
