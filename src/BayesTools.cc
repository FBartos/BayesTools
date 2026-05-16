#include <module/Module.h>

#include "distributions/DBTLKJCPC.h"
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
      insert(new BTLKJCholesky);
      insert(new BTLKJCorr);
    }

    BayesToolsModule::~BayesToolsModule()
    {
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
