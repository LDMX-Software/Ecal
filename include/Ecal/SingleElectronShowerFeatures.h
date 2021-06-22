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
  /// x,y coordinate pair
  typedef std::pair<float, float> XYCoords;

  /**
   * Normal blank constructor to register as a producer
   */
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
   *  - Total energy in ECal (summedDet)
   *  - Total isolated energy (summedTightIso)
   *  - Energy after layer 20 (ecalBackEnergy)
   *  - Maximum energy deposited in a single cell (maxCellDep)
   *  - Transverse RMS of all hits in shower (showerRMS)
   *  - STD of x-position of all hits (xStd)
   *  - STD of y-position of all hits (yStd)
   *  - Average layer hit (avgLayerHit)
   *  - STD of layers hit (stdLayerHit)
   *  - Num hits readout (nReadoutHits)
   *  - Deepest layer hit (deepestLayerHit)
   *  - Energy within electron containment region for 5 radii
   *  - Energy within photon containment region for 5 radii
   *  - N Hits, X and Y mean and STD for hits outside 
   *    either electron or photon containment region for 5 radii
   */
  void produce(framework::Event& event);

 private:

  bool verbose_{false};

  /// number of layers in the ECal
  int nEcalLayers_{0};
  /// name of rec pass to use
  std::string rec_pass_name_;
  /// name of rec collection to use
  std::string rec_coll_name_;
};

}  // namespace ecal

#endif
