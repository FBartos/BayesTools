#include <module/Module.h>

#include "distributions/DBTLKJCPC.h"
#include "distributions/DBTInvGamma.h"
#include "distributions/DBTMoment.h"
#include "distributions/DBTInvMoment.h"
#include "functions/BTLKJCholesky.h"

namespace jags {
  namespace BayesTools {

    class BayesToolsModule : public Module
    {
    public:
      BayesToolsModule();
      ~BayesToolsModule();
    };

    BayesToolsModule::BayesToolsModule() : Module("BayesTools")
    {
      insert(new DBTLKJCPC);
      insert(new DBTInvGamma);
      insert(new DBTMoment);
      insert(new DBTInvMoment);
      insert(new BTLKJCholesky);
      insert(new BTLKJCorr);
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
