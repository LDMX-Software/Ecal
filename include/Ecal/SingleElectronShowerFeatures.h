#ifndef ECAL_SINGLEELECTRONSHOWERFEATURES_H
#define ECAL_SINGLEELECTRONSHOWERFEATURES_H

// LDMX
#include "DetDescr/EcalHexReadout.h"
#include "DetDescr/EcalID.h"
#include "Ecal/Event/EcalHit.h"
#include "Ecal/Event/EcalVetoResult.h"
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h"

// C++
#include <map>
#include <memory>

namespace ecal {

/**
 * @class SingleElectronShowerFeatures
 * @brief Calculates a variety of shower features assuming
 * a single primary electron in the event.
 */
class SingleElectronShowerFeatures : public framework::Producer {
 public:
  /// pair a cell's ID with the energy deposited in it
  typedef std::pair<ldmx::EcalID, float> CellEnergyPair;

  /// x,y coordinate pair
  typedef std::pair<float, float> XYCoords;

  SingleElectronShowerFeatures(const std::string& name, framework::Process& process)
      : Producer(name, process) {}

  /**
   * Configure the processor using the given user specified parameters.
   *
   * @param parameters Set of parameters used to configure this processor.
   */
  void configure(framework::config::Parameters& parameters) final override;

  /**
   * Calculate the single electron shower features.
   *
   * The features that are calculated here are
   *  - 
   *
   */
  void produce(framework::Event& event);

 private:

  std::vector<float> ecalLayerEdepRaw_;
  std::vector<float> ecalLayerEdepReadout_;
  std::vector<float> ecalLayerTime_;

  int nEcalLayers_{0};
  int backEcalStartingLayer_{0};
  int nReadoutHits_{0};
  int deepestLayerHit_{0};
  int doBdt_{0};

  double summedDet_{0};
  double summedTightIso_{0};
  double maxCellDep_{0};
  double showerRMS_{0};
  double xStd_{0};
  double yStd_{0};
  double avgLayerHit_{0};
  double stdLayerHit_{0};
  double ecalBackEnergy_{0};

  bool verbose_{false};

  /// name of rec pass to use
  std::string rec_pass_name_;
  /// name of rec collection to use
  std::string rec_coll_name_;
};

}  // namespace ecal

#endif
