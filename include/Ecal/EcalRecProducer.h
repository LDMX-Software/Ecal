/**
 * @file EcalRecProducer.h
 * @brief Class that performs basic ECal digitization
 * @author Owen Colegrove, UCSB
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef ECAL_ECALRECPRODUCER_H_
#define ECAL_ECALRECPRODUCER_H_

//----------------//
//   C++ StdLib   //
//----------------//
#include <memory>  //for smart pointers

//----------//
//   LDMX   //
//----------//
#include "DetDescr/DetectorID.h"
#include "DetDescr/EcalID.h"
#include "Framework/EventDef.h"
#include "Framework/EventProcessor.h"

namespace ecal {

/**
 * @class EcalRecProducer
 * @brief Performs basic ECal reconstruction
 *
 * Reconstruction is done from the EcalDigi samples.
 * Some hard-coded parameters are used for position and energy calculation.
 */
class EcalRecProducer : public framework::Producer {
 public:
  /**
   * Constructor
   */
  EcalRecProducer(const std::string& name, framework::Process& process);

  /**
   * Destructor
   */
  virtual ~EcalRecProducer();

  /**
   * Grabs configure parameters from the python config file.
   */
  virtual void configure(framework::config::Parameters&);

  /**
   * Produce EcalHits and put them into the event bus using the
   * EcalDigis as input.
   *
   * This function unfolds the digi samples taken by the HGC ROC
   * and reconstructs their energy using knowledge of how
   * the chip operates and the position using EcalGeometry.
   */
  virtual void produce(framework::Event& event);

 private:
  /** Digi Collection Name to use as input */
  std::string digiCollName_;

  /** Digi Pass Name to use as input */
  std::string digiPassName_;

  /// simhit collection name
  std::string simHitCollName_;

  /// simhit pass name
  std::string simHitPassName_;

  /// output hit collection name
  std::string recHitCollName_;

  /// store intermediate reconstruction values
  bool store_intermediate_values_;

  /// Energy [MeV] deposited by a MIP in Si 0.5mm thick
  double mip_si_energy_;

  /// Length of clock cycle [ns]
  double clock_cycle_;

  /// Number of electrons generated by average MIP in Si 0.5mm thick
  double charge_per_mip_;

  /**
   * Layer Weights to use for this reconstruction
   *
   * Layer weights account for the energy lost in the absorber directly
   * in front of the Silicon layer where the measured energy was deposited.
   * These are determined by calculating the average amount of energy lost
   * by a MIP passing through the extra material between sensitive layers.
   */
  std::vector<double> layerWeights_;

  /**
   * Second Order Energy Correction to use for this reconstruction
   *
   * This is a shift applied to all of the energies in order to have the
   * mean of the total energy deposited in the ECal be accurate.
   * This is less physically motivated than the layer weights and is more
   * of a calibration number.
   */
  double secondOrderEnergyCorrection_;
};
}  // namespace ecal

#endif
