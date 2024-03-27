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
  /*
   * Since I defined these member variables to have sensible values in the declaration above,
   * we DO NOT get their parameter values from the python config. Below, I have written out
   * the lines that we could use to get their values from the python config if changing them
   * is desirable (e.g. if you want to run a few different passes of digi and recon with different
   * parameters and then want to analyze the different passes).
   */
  /*
  digi_collection_name_ = parameters.getParameter("digi_collection_name");
  digi_pass_name_ = parameters.getParameter("digi_pass_name");
  trig_collection_name_ = parameters.getParameter("trig_collection_name");
  trig_pass_name_ = parameters.getParameter("trig_pass_name");
  hit_collection_name_ = parameters.getParameter("hit_collection_name");
  hit_pass_name_ = parameters.getParameter("hit_pass_name");
  */
}

void TrigPrimResolutionAnalyzer::onProcessStart() {
  // initialize processing by making histograms and such
  // first, we get the directory for this processor in the histogram file
  getHistoDirectory();
  // then we can create histograms within it
  histograms_.create(
      "total_trig_energy" /* name - as written in output ROOT file */,
      "Total of all Trig Digis" /* xlabel - axis label of histogram */,
      1000 /* number of bins */, 0 /* minimum value */, 10000 /* maximum value */
  ); // WARNING: I don't think this binning is good!! I just picked a random number!!
  /*
   * 2D Histograms are also possible
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

  int total{0};
  for (const auto& trig : trigs) {
    total += trig.linearPrimitive();
  }

  /*
   * after picking some dummy bins, I printout the values I'm going to fill
   * so I can see what an actual binning should be. For 10 events, the output was
   * 978
   * 949
   * 1070
   * 858
   * 846
   * 774
   * 853
   * 859
   * 887
   * 917
   * so I changed the binning to 1k bins ranging from 0 to 10k
  std::cout << total << std::endl;
   */
  histograms_.fill("total_trig_energy", total);
  
  /*
   * 2D fill example
  histograms_.fill("name", xvalue, yvalue);
  */
}

}

DECLARE_ANALYZER_NS(ldmx::ecal, TrigPrimResolutionAnalyzer);
