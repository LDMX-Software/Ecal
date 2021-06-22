#include "Ecal/SingleElectronShowerFeatures.h"

// LDMX
#include "DetDescr/SimSpecialID.h"
#include "Ecal/Event/EcalHit.h"
#include "Recon/Event/EventConstants.h"

/*~~~~~~~~~~~*/
/*   Tools   */
/*~~~~~~~~~~~*/
#include "Tools/AnalysisUtils.h"

// C++
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <fstream>

namespace ecal {

// arrays holding 68% containment radius per layer for different bins in
// momentum/angle
static const std::vector<double> radius68_thetalt10_plt500 = {
    4.045666158618167,  4.086393662224346,  4.359141107602775,
    4.666549994726691,  5.8569181911416015, 6.559716356124256,
    8.686967529043072,  10.063482736354674, 13.053528344041274,
    14.883496407943747, 18.246694748611368, 19.939799900443724,
    22.984795944506224, 25.14745829663406,  28.329169392203216,
    29.468032123356345, 34.03271241527079,  35.03747443690781,
    38.50748727211848,  39.41576583301171,  42.63622296033334,
    45.41123601592071,  48.618139095742876, 48.11801717451056,
    53.220539860213655, 58.87753380915155,  66.31550881539764,
    72.94685877928593,  85.95506228335348,  89.20607201266672,
    93.34370253818409,  96.59471226749734,  100.7323427930147,
    103.98335252232795};

static const std::vector<double> radius68_thetalt10_pgt500 = {
    4.081926458777424,  4.099431732299409,  4.262428482867968,
    4.362017581473145,  4.831341579961153,  4.998346041276382,
    6.2633736512415705, 6.588371889265881,  8.359969947444522,
    9.015085558044309,  11.262722588206483, 12.250305471269183,
    15.00547660437276,  16.187264014640103, 19.573764900578503,
    20.68072032434797,  24.13797140783321,  25.62942209291236,
    29.027596514735617, 30.215039667389316, 33.929540248019585,
    36.12911729771914,  39.184563500620946, 42.02062468386282,
    46.972125628650204, 47.78214816041894,  55.88428562462974,
    59.15520134927332,  63.31816666637158,  66.58908239101515,
    70.75204770811342,  74.022963432757,    78.18592874985525,
    81.45684447449884};

static const std::vector<double> radius68_theta10to20 = {
    4.0251896715647115, 4.071661598616328,  4.357690094817289,
    4.760224640141712,  6.002480766325418,  6.667318981016246,
    8.652513285172342,  9.72379373302137,   12.479492693251478,
    14.058548828317289, 17.544872909347912, 19.43616066939176,
    23.594162859513734, 25.197329065282954, 29.55995803074302,
    31.768946746958296, 35.79247330197688,  37.27810357669942,
    41.657281051476545, 42.628141392692626, 47.94208483539388,
    49.9289473559796,   54.604030254423975, 53.958762417361655,
    53.03339560920388,  57.026277390001425, 62.10810455035879,
    66.10098633115634,  71.1828134915137,   75.17569527231124,
    80.25752243266861,  84.25040421346615,  89.33223137382352,
    93.32511315462106};

static const std::vector<double> radius68_thetagt20 = {
    4.0754238481177705, 4.193693485630508,  5.14209420056253,
    6.114996249971468,  7.7376807326481645, 8.551663213602291,
    11.129110612057813, 13.106293737495639, 17.186617323282082,
    19.970887612094604, 25.04088272634407,  28.853696411302344,
    34.72538105333071,  40.21218694947545,  46.07344239520299,
    50.074953583805346, 62.944045771758645, 61.145621459396814,
    69.86940198299047,  74.82378572939959,  89.4528387422834,
    93.18228303096758,  92.51751129204555,  98.80228884380018,
    111.17537347472128, 120.89712563907408, 133.27021026999518,
    142.99196243434795, 155.36504706526904, 165.08679922962185,
    177.45988386054293, 187.18163602489574, 199.55472065581682,
    209.2764728201696};

void SingleElectronShowerFeatures::configure(framework::config::Parameters &parameters) {

  nEcalLayers_ = parameters.getParameter<int>("num_ecal_layers");

  // Set the collection name as defined in the configuration
  rec_pass_name_ = parameters.getParameter<std::string>("rec_pass_name");
  rec_coll_name_ = parameters.getParameter<std::string>("rec_coll_name");
}

void SingleElectronShowerFeatures::produce(framework::Event &event) {

  /****************************************************************************
   * Get Electron+Photon Tracks Upstream of ECAL
   *  - ignored if collection doesn't exist
   ***************************************************************************/

  std::vector<double> recoilP;
  std::vector<float> recoilPos;
  std::vector<double> recoilPAtTarget;
  std::vector<float> recoilPosAtTarget;

  if (event.exists("EcalScoringPlaneHits")) {
    //
    // Loop through all of the sim particles and find the recoil electron.
    //

    // Get the collection of simulated particles from the event
    auto particleMap{event.getMap<int, ldmx::SimParticle>("SimParticles")};

    // Search for the recoil electron
    auto [recoilTrackID, recoilElectron] = Analysis::getRecoil(particleMap);

    // Find ECAL SP hit for recoil electron
    auto ecalSpHits{
        event.getCollection<ldmx::SimTrackerHit>("EcalScoringPlaneHits")};
    float pmax = 0;
    for (ldmx::SimTrackerHit &spHit : ecalSpHits) {
      
      ldmx::SimSpecialID hit_id(spHit.getID());
      if (hit_id.plane() != 31 || spHit.getMomentum()[2] <= 0) continue;

      if (spHit.getTrackID() == recoilTrackID) {
        if (sqrt(pow(spHit.getMomentum()[0], 2) +
                 pow(spHit.getMomentum()[1], 2) +
                 pow(spHit.getMomentum()[2], 2)) > pmax) {
          recoilP = spHit.getMomentum();
          recoilPos = spHit.getPosition();
          pmax = sqrt(pow(recoilP[0], 2) + pow(recoilP[1], 2) +
                      pow(recoilP[2], 2));
        }
      }
    }

    // Find target SP hit for recoil electron
    if (event.exists("TargetScoringPlaneHits")) {
      std::vector<ldmx::SimTrackerHit> targetSpHits =
          event.getCollection<ldmx::SimTrackerHit>("TargetScoringPlaneHits");
      pmax = 0;
      for (ldmx::SimTrackerHit &spHit : targetSpHits) {
        ldmx::SimSpecialID hit_id(spHit.getID());
        if (hit_id.plane() != 1 || spHit.getMomentum()[2] <= 0) continue;

        if (spHit.getTrackID() == recoilTrackID) {
          if (sqrt(pow(spHit.getMomentum()[0], 2) +
                   pow(spHit.getMomentum()[1], 2) +
                   pow(spHit.getMomentum()[2], 2)) > pmax) {
            recoilPAtTarget = spHit.getMomentum();
            recoilPosAtTarget = spHit.getPosition();
            pmax =
                sqrt(pow(recoilPAtTarget[0], 2) + pow(recoilPAtTarget[1], 2) +
                     pow(recoilPAtTarget[2], 2));
          }
        }
      }
    }
  }

  // Get the Ecal Geometry
  const ldmx::EcalHexReadout &hexReadout = getCondition<ldmx::EcalHexReadout>(
      ldmx::EcalHexReadout::CONDITIONS_OBJECT_NAME);

  /****************************************************************************
   * Project Electron+Photon Trajectories into ECAL
   ***************************************************************************/

  std::vector<XYCoords> ele_trajectory, photon_trajectory;
  if (recoilP.size() > 0) {
    /* Calculate where trajectory intersects ECAL layers using position and momentum
     * at scoring plane 
     */
    auto getTrajectory[&hexReadout](std::vector<double> momentum,
                                    std::vector<float> position)
        ->std::vector<XYCoords> {
      std::vector<XYCoords> positions;
      for (int iLayer = 0; iLayer < nEcalLayers_; iLayer++) {
        float posX =
            position[0] + (momentum[0] / momentum[2]) *
                              (hexReadout.getZPosition(iLayer) - position[2]);
        float posY =
            position[1] + (momentum[1] / momentum[2]) *
                              (hexReadout.getZPosition(iLayer) - position[2]);
        positions.push_back(std::make_pair(posX, posY));
      }
      return positions;
    }

    ele_trajectory = getTrajectory(recoilP, recoilPos);
    std::vector<double> pvec = recoilPAtTarget.size()
                                   ? recoilPAtTarget
                                   : std::vector<double>{0.0, 0.0, 0.0};
    std::vector<float> posvec = recoilPosAtTarget.size()
                                    ? recoilPosAtTarget
                                    : std::vector<float>{0.0, 0.0, 0.0};
    photon_trajectory =
        getTrajectory({-pvec[0], -pvec[1], 4000.0 - pvec[2]}, posvec);
  }

  float recoilPMag =
      recoilP.size()
          ? sqrt(pow(recoilP[0], 2) + pow(recoilP[1], 2) + pow(recoilP[2], 2))
          : -1.0;
  float recoilTheta = recoilPMag > 0 ? recoilP[2] / recoilPMag : -1.0;

  std::vector<double> ele_radii = radius68_thetalt10_plt500;
  if (recoilTheta < 10 && recoilPMag >= 500)
    ele_radii = radius68_thetalt10_pgt500;
  else if (recoilTheta >= 10 && recoilTheta < 20)
    ele_radii = radius68_theta10to20;
  else if (recoilTheta >= 20)
    ele_radii = radius68_thetagt20;
  // Use default binning
  std::vector<double> photon_radii = radius68_thetalt10_plt500;

  // Get the collection of digitized Ecal hits from the event.
  const auto& ecalRecHits = event.getCollection<ldmx::EcalHit>(
      rec_coll_name_, rec_pass_name_);

  /****************************************************************************
   * Calculate Energy-Weighted Global Centroid of Shower ~~ O(n)
   ***************************************************************************/

  std::map<ldmx::EcalID, float> cell_to_edep;
  float sumEdep = 0;
  auto wgtCentroidCoords = std::make_pair<float, float>(0., 0.);
  for (const ldmx::EcalHit &hit : ecalRecHits) {
    ldmx::EcalID id(hit.getID());
    auto cell_energy_pair = std::make_pair(id, hit.getEnergy());
    XYCoords centroidCoords = hexReadout.getCellCentroidXYPair(id);
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
    XYCoords centroidCoords = hexReadout.getCellCentroidXYPair(hitID(hit));

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
  ldmx::EcalID globalCentroid = ldmx::EcalID(0, globalCentroid.module(), globalCentroid.cell());

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
      showerRMS{0}, xStd{0}, yStd{0}, avgLayerHit{0}, stdLayerHit{0};

  for (const ldmx::EcalHit &hit : ecalRecHits) {
    // Layer-wise quantities
    ldmx::EcalID id = hitID(hit);
    ecalLayerEdepRaw[id.layer()] += hit.getEnergy();
    if (id.layer() >= 20) ecalBackEnergy += hit.getEnergy();
    if (maxCellDep < hit.getEnergy()) maxCellDep = hit.getEnergy();
    if (hit.getEnergy() > 0) {
      nReadoutHits++;
      ecalLayerEdepReadout[id.layer()] += hit.getEnergy();
      ecalLayerTime[id.layer()] += (hit.getEnergy()) * hit.getTime();
      XYCoords xy_pair = hexReadout.getCellCentroidXYPair(id);
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
    if (energy > 0) summedTightIso_ += energy;
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
    avgLayerHit_ = 0;
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
    ldmx::EcalID id = hitID(hit);
    XYCoords xy_pair = hexReadout.getCellCentroidXYPair(id);
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

  if (nReadoutHits_ > 0) {
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

  result.setVariables(
      nReadoutHits, deepestLayerHit, summedDet, summedTightIso, maxCellDep,
      showerRMS, xStd, yStd, avgLayerHit, stdLayerHit, ecalBackEnergy,
      nStraightTracks, nLinregTracks, firstNearPhLayer, epAng, epSep, 
      electronContainmentEnergy, photonContainmentEnergy,
      outsideContainmentEnergy, outsideContainmentNHits, outsideContainmentXstd,
      outsideContainmentYstd, ecalLayerEdepReadout, recoilP, recoilPos);

  event.add(getName(), result);
}

}  // namespace ecal

DECLARE_PRODUCER_NS(ecal, SingleElectronShowerFeatures);
