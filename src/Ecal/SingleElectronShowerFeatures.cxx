
// LDMX
#include "DetDescr/SimSpecialID.h"
#include "DetDescr/EcalHexReadout.h"
#include "DetDescr/EcalID.h"
#include "Ecal/Event/EcalHit.h"
#include "Ecal/Event/EcalVetoResult.h"
#include "Framework/Configure/Parameters.h"
#include "Framework/EventProcessor.h"
#include "Recon/Event/EventConstants.h"
#include "SimCore/Event/SimParticle.h"
#include "SimCore/Event/SimTrackerHit.h"

/*~~~~~~~~~~~*/
/*   Tools   */
/*~~~~~~~~~~~*/
#include "Tools/AnalysisUtils.h"

// C++
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <fstream>
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


void SingleElectronShowerFeatures::configure(framework::config::Parameters &parameters) {
  nEcalLayers_ = parameters.getParameter<int>("num_ecal_layers");

  // Set the collection name as defined in the configuration
  rec_pass_name_ = parameters.getParameter<std::string>("rec_pass_name");
  rec_coll_name_ = parameters.getParameter<std::string>("rec_coll_name");
}

void SingleElectronShowerFeatures::produce(framework::Event &event) {
  // Get the Ecal Geometry
  const ldmx::EcalHexReadout &hexReadout = getCondition<ldmx::EcalHexReadout>(
      ldmx::EcalHexReadout::CONDITIONS_OBJECT_NAME);

  const auto& ele_trajectory{event.getCollection<std::pair<float,float>>(proj_ele_traj_, pass_)};
  const auto& photon_trajectory{event.getCollection<std::pair<float,float>>(proj_photon_traj_, pass_)};
  const auto& ele_radii{event.getCollection<double>(ele_roc_, pass_)};
  const auto& photon_radii{event.getCollection<double>(photon_roc_, pass_)};

  // Get the collection of reconstructed Ecal hits from the event.
  const auto& ecalRecHits = event.getCollection<ldmx::EcalHit>(
      rec_coll_name_, rec_pass_name_);

  /****************************************************************************
   * Calculate Energy-Weighted Global Centroid of Shower ~~ O(n)
   ***************************************************************************/

  std::map<ldmx::EcalID, float> cell_to_edep;
  float sumEdep = 0;
  float showerRMS = 0;
  auto wgtCentroidCoords = std::make_pair<float, float>(0., 0.);
  for (const ldmx::EcalHit &hit : ecalRecHits) {
    ldmx::EcalID id(hit.getID());
    auto cell_energy_pair = std::make_pair(id, hit.getEnergy());
    XYCoords centroidCoords = hexReadout.getCellCenterAbsolute(id);
    wgtCentroidCoords.first = wgtCentroidCoords.first +
                              centroidCoords.first * cell_energy_pair.second;
    wgtCentroidCoords.second = wgtCentroidCoords.second +
                               centroidCoords.second * cell_energy_pair.second;
    sumEdep += cell_energy_pair.second;
    cell_to_edep.emplace(id, hit.getEnergy());
  }
  wgtCentroidCoords.first = (sumEdep > 1E-6) ? wgtCentroidCoords.first / sumEdep
                                             : wgtCentroidCoords.first;
  wgtCentroidCoords.second = (sumEdep > 1E-6)
                                 ? wgtCentroidCoords.second / sumEdep
                                 : wgtCentroidCoords.second;

  // Find Nearest Cell to Centroid
  ldmx::EcalID globalCentroid;
  float maxDist = 1e6;
  for (const ldmx::EcalHit &hit : ecalRecHits) {
    XYCoords centroidCoords = hexReadout.getCellCenterAbsolute(ldmx::EcalID(hit.getID()));

    float deltaR =
        pow(pow((centroidCoords.first - wgtCentroidCoords.first), 2) +
                pow((centroidCoords.second - wgtCentroidCoords.second), 2),
            .5);
    showerRMS += deltaR * hit.getEnergy();
    if (deltaR < maxDist) {
      maxDist = deltaR;
      globalCentroid = ldmx::EcalID(hit.getID());
    }
  }
  if (sumEdep > 0) showerRMS = showerRMS / sumEdep;

  // flatten global centroid to zero layer
  globalCentroid = ldmx::EcalID(0, globalCentroid.module(), globalCentroid.cell());

  /****************************************************************************
   * Determine Map of Isolated Hits ~~ O(n)
   *  - needs global centroid and full map of cell IDs to energy deposited
   ***************************************************************************/

  bool doTight = true;

  std::map<ldmx::EcalID, float> isocell_to_edep;
  for (const ldmx::EcalHit &hit : ecalRecHits) {
    auto isolatedHit = std::make_pair(true, ldmx::EcalID());
    ldmx::EcalID id(hit.getID());
    ldmx::EcalID flatid(0, id.module(), id.cell());
    if (doTight) {
      // Disregard hits that are on the centroid.
      if (flatid == globalCentroid) 
        continue;

      // Skip hits that are on centroid inner ring
      if (hexReadout.isNN(globalCentroid, flatid))
        continue;
    }

    // Skip hits that have a readout neighbor
    // Get neighboring cell id's and try to look them up in the full cell map
    // (constant speed algo.)
    //  these ideas are only cell/module (must ignore layer)
    const std::vector<ldmx::EcalID>& cellNbrIds = hexReadout.getNN(id);

    for (const ldmx::EcalID& flat_neighbor : cellNbrIds) {
      // update neighbor ID to the current layer
      ldmx::EcalID neighbor(id.layer(), flat_neighbor.module(), 
          flat_neighbor.cell());
      // look in cell hit map to see if it is there
      if (cell_to_edep.find(neighbor) != cell_to_edep.end()) {
        isolatedHit = std::make_pair(false, neighbor);
        break;
      }
    }

    // skip if not isolated
    if (!isolatedHit.first)
      continue;

    // Insert isolated hit
    isocell_to_edep.emplace(id, hit.getEnergy());
  }

  /****************************************************************************
   * Containment Variable Calculation ~~ O(n)
   ***************************************************************************/

  float wavgLayerHit = 0;
  float xMean = 0;
  float yMean = 0;

  // Containment variables
  unsigned int nregions = 5;
  std::vector<int> outsideContainmentNHits(nregions, 0);
  std::vector<float> electronContainmentEnergy(nregions, 0.0),
      photonContainmentEnergy(nregions, 0.0),
      outsideContainmentEnergy(nregions, 0.0),
      outsideContainmentXmean(nregions, 0.0),
      outsideContainmentYmean(nregions, 0.0),
      outsideContainmentXstd(nregions, 0.0),
      outsideContainmentYstd(nregions, 0.0);
  std::vector<float> ecalLayerEdepRaw(nEcalLayers_, 0),
      ecalLayerEdepReadout(nEcalLayers_, 0), ecalLayerTime(nEcalLayers_, 0);
  int nReadoutHits{0}, deepestLayerHit{0};
  double ecalBackEnergy{0}, maxCellDep{0}, summedDet{0}, summedTightIso{0},
      xStd{0}, yStd{0}, avgLayerHit{0}, stdLayerHit{0};

  for (const ldmx::EcalHit &hit : ecalRecHits) {
    // Layer-wise quantities
    ldmx::EcalID id(hit.getID());
    ecalLayerEdepRaw[id.layer()] += hit.getEnergy();
    if (id.layer() >= 20) ecalBackEnergy += hit.getEnergy();
    if (maxCellDep < hit.getEnergy()) maxCellDep = hit.getEnergy();
    if (hit.getEnergy() > 0) {
      nReadoutHits++;
      ecalLayerEdepReadout[id.layer()] += hit.getEnergy();
      ecalLayerTime[id.layer()] += (hit.getEnergy()) * hit.getTime();
      XYCoords xy_pair = hexReadout.getCellCenterAbsolute(id);
      xMean += xy_pair.first * hit.getEnergy();
      yMean += xy_pair.second * hit.getEnergy();
      avgLayerHit += id.layer();
      wavgLayerHit += id.layer() * hit.getEnergy();
      if (deepestLayerHit < id.layer()) {
        deepestLayerHit = id.layer();
      }
      float distance_ele_trajectory =
          ele_trajectory.size()
              ? sqrt(
                    pow((xy_pair.first - ele_trajectory[id.layer()].first), 2) +
                    pow((xy_pair.second - ele_trajectory[id.layer()].second),
                        2))
              : -1.0;
      float distance_photon_trajectory =
          photon_trajectory.size()
              ? sqrt(
                    pow((xy_pair.first - photon_trajectory[id.layer()].first),
                        2) +
                    pow((xy_pair.second - photon_trajectory[id.layer()].second),
                        2))
              : -1.0;
      // Decide which region a hit goes into and add to sums
      for (unsigned int ireg = 0; ireg < nregions; ireg++) {
        if (distance_ele_trajectory >= ireg * ele_radii[id.layer()] &&
            distance_ele_trajectory < (ireg + 1) * ele_radii[id.layer()])
          electronContainmentEnergy[ireg] += hit.getEnergy();
        if (distance_photon_trajectory >= ireg * photon_radii[id.layer()] &&
            distance_photon_trajectory < (ireg + 1) * photon_radii[id.layer()])
          photonContainmentEnergy[ireg] += hit.getEnergy();
        if (distance_ele_trajectory > (ireg + 1) * ele_radii[id.layer()] &&
            distance_photon_trajectory >
                (ireg + 1) * photon_radii[id.layer()]) {
          outsideContainmentEnergy[ireg] += hit.getEnergy();
          outsideContainmentNHits[ireg] += 1;
          outsideContainmentXmean[ireg] += xy_pair.first * hit.getEnergy();
          outsideContainmentYmean[ireg] += xy_pair.second * hit.getEnergy();
        }
      }
    }
  }

  for (const auto &[id, energy] : isocell_to_edep) {
    if (energy > 0) summedTightIso += energy;
  }

  for (int iLayer = 0; iLayer < ecalLayerEdepReadout.size(); iLayer++) {
    ecalLayerTime[iLayer] /= ecalLayerEdepReadout[iLayer];
    summedDet += ecalLayerEdepReadout[iLayer];
  }

  if (nReadoutHits > 0) {
    avgLayerHit /= nReadoutHits;
    wavgLayerHit /= summedDet;
    xMean /= summedDet;
    yMean /= summedDet;
  } else {
    wavgLayerHit = 0;
    avgLayerHit = 0;
    xMean = 0;
    yMean = 0;
  }

  for (unsigned int ireg = 0; ireg < nregions; ireg++) {
    if (outsideContainmentEnergy[ireg] > 0) {
      outsideContainmentXmean[ireg] /= outsideContainmentEnergy[ireg];
      outsideContainmentYmean[ireg] /= outsideContainmentEnergy[ireg];
    }
  }

  /****************************************************************************
   * Second Loop to Calculate Deviations ~~ O(n)
   ***************************************************************************/

  for (const ldmx::EcalHit &hit : ecalRecHits) {
    ldmx::EcalID id(hit.getID());
    XYCoords xy_pair = hexReadout.getCellCenterAbsolute(id);
    if (hit.getEnergy() > 0) {
      xStd +=
          pow((xy_pair.first - xMean), 2) * hit.getEnergy();
      yStd +=
          pow((xy_pair.second - yMean), 2) * hit.getEnergy();
      stdLayerHit += pow((id.layer() - wavgLayerHit), 2) * hit.getEnergy();
    }
    float distance_ele_trajectory =
        ele_trajectory.size()
            ? sqrt(pow((xy_pair.first - ele_trajectory[id.layer()].first), 2) +
                   pow((xy_pair.second - ele_trajectory[id.layer()].second), 2))
            : -1.0;
    float distance_photon_trajectory =
        photon_trajectory.size()
            ? sqrt(pow((xy_pair.first - photon_trajectory[id.layer()].first),
                       2) +
                   pow((xy_pair.second - photon_trajectory[id.layer()].second),
                       2))
            : -1.0;
    for (unsigned int ireg = 0; ireg < nregions; ireg++) {
      if (distance_ele_trajectory > (ireg + 1) * ele_radii[id.layer()] &&
          distance_photon_trajectory > (ireg + 1) * photon_radii[id.layer()]) {
        outsideContainmentXstd[ireg] +=
            pow((xy_pair.first - outsideContainmentXmean[ireg]), 2) *
            hit.getEnergy();
        outsideContainmentYstd[ireg] +=
            pow((xy_pair.second - outsideContainmentYmean[ireg]), 2) *
            hit.getEnergy();
      }
    }
  }

  if (nReadoutHits > 0) {
    xStd = sqrt(xStd / summedDet);
    yStd = sqrt(yStd / summedDet);
    stdLayerHit = sqrt(stdLayerHit / summedDet);
  } else {
    xStd = 0;
    yStd = 0;
    stdLayerHit = 0;
  }

  for (unsigned int ireg = 0; ireg < nregions; ireg++) {
    if (outsideContainmentEnergy[ireg] > 0) {
      outsideContainmentXstd[ireg] =
          sqrt(outsideContainmentXstd[ireg] / outsideContainmentEnergy[ireg]);
      outsideContainmentYstd[ireg] =
          sqrt(outsideContainmentYstd[ireg] / outsideContainmentEnergy[ireg]);
    }
  }

  ldmx::EcalVetoResult result;
  result.setVariables(
      nReadoutHits, deepestLayerHit, summedDet, summedTightIso, maxCellDep,
      showerRMS, xStd, yStd, avgLayerHit, stdLayerHit, ecalBackEnergy,
      0/*nStraightTracks*/, 0/*nLinregTracks*/, 0/*firstNearPhLayer*/, 0/*epAng*/, 0/*epSep*/, 
      electronContainmentEnergy, photonContainmentEnergy,
      outsideContainmentEnergy, outsideContainmentNHits, outsideContainmentXstd,
      outsideContainmentYstd, ecalLayerEdepReadout, recoilP, recoilPos);

  event.add(getName(), result);
}

}  // namespace ecal

DECLARE_PRODUCER_NS(ecal, SingleElectronShowerFeatures);
