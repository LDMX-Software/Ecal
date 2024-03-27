#include "Framework/EventProcessor.h"

#include "Ecal/Event/EcalHit.h"
#include "Recon/Event/HgcrocTrigDigi.h"
#include "Recon/Event/HgcrocDigiCollection.h"
#include "DetDescr/EcalID.h"
#include "DetDescr/EcalTriggerID.h"

namespace ldmx::ecal {
class TrigPrimResolutionAnalyzer : public framework::Analyzer {
  std::string digi_collection_name_ = "EcalDigis";
  std::string digi_pass_name_ = "";
  std::string trig_collection_name_ = "ecalTrigDigis";
  std::string trig_pass_name_ = "";
  std::string hit_collection_name_ = "EcalRecHits";
  std::string hit_pass_name_ = "";
 public:
  TrigPrimResolutionAnalyzer(const std::string& name, framework::Process& process)
    : framework::Analyzer(name, process) {}
  virtual ~TrigPrimResolutionAnalyzer() = default;
  void configure(framework::config::Parameters& parameters) final;
  void onProcessStart() final;
  void analyze(const framework::Event& event) final;
};

void TrigPrimResolutionAnalyzer::configure(framework::config::Parameters& parameters) {
  digi_collection_name_ = parameters.getParameter("digi_collection_name", digi_collection_name_);
  digi_pass_name_ = parameters.getParameter("digi_pass_name", digi_pass_name_);
  trig_collection_name_ = parameters.getParameter("trig_collection_name", trig_collection_name_);
  trig_pass_name_ = parameters.getParameter("trig_pass_name", trig_pass_name_);
  hit_collection_name_ = parameters.getParameter("hit_collection_name", hit_collection_name_);
  hit_pass_name_ = parameters.getParameter("hit_pass_name", hit_pass_name_);
}

void TrigPrimResolutionAnalyzer::onProcessStart() {
  // initialize processing by making histograms and such
  // first, we get the directory for this processor in the histogram file
  getHistoDirectory();
  // then we can create histograms within it
  /*
  histograms_.create(
      "name",
      "xlabel", nbins, xstart, xend,
      "ylabel", nbins, ystart, yend
  );
  */
}

void TrigPrimResolutionAnalyzer::analyze(const framework::Event& event) {
  // called once on each event, get objects and fill histograms
  const auto& trigs = event.getCollection<ldmx::HgcrocTrigDigi>(trig_collection_name_, trig_pass_name_);
  // trigs are a std::vector<ldmx::HgcrocTrigDigi>
  const auto& digis = event.getObject<ldmx::HgcrocDigiCollection>(digi_collection_name_, digi_pass_name_);
  // digis are a ldmx::HgcrocDigiCollection
  const auto& hits  = event.getCollection<ldmx::EcalHit>(hit_collection_name_, hit_pass_name_);
  // hits are a std::vector<ldmx::EcalHit>
  
  /*
  histograms_.fill(xvalue, yvalue);
  */
}

}

DECLARE_ANALYZER_NS(ldmx::ecal, TrigPrimResolutionAnalyzer);
