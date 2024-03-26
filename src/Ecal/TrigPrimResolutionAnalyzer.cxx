#include "Framework/Analyzer.h"

#include "Ecal/Event/EcalHit.h"
#include "Recon/Event/HgcrocTrigDigi.h"
#include "Recon/Event/HgcrocDigiCollection.h"

namespace ldmx::ecal {
class TrigPrimResolutionAnalyzer : public framework::Analyzer {
  // private member variables (if need be)
 public:
  TrigPrimResolutionAnalyzer(const std::string& name, framework::Process& process)
    : framework::Analyzer(name, process) {}
  virtual ~TrigPrimResolutionAnalyzer() = default;
  void configure(framework::config::Parameters& parameters) final;
  void onProcessStart() final;
  void analyze(const framework::Event& event) final;
};

void TrigPrimResolutionAnalyzer::configure(framework::config::Parameters& parameters) {
  // configure yourself by getting parameters from python (if need b)
}

void TrigPrimResolutionAnalyzer::onProcessStart() final {
  // initialize processing by making histograms and such
}

void TrigPrimResolutionAnalyzer::analyze(const framework::Event& event) final {
  // called once on each event, get objects and fill histograms
}

}

DECLARE_ANALYZER_NS(ldmx::ecal, TrigPrimResolutionAnalyzer);
