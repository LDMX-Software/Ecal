#include "Framework/EventProcessor.h"

#include "SimCore/Event/SimTrackerHit.h"
#include "SimCore/Event/SimParticle.h"

namespace ecal {

class SimRecoilProjector : public framework::Producer {
 public:
  SimRecoilProjector(const std::string& name, framework::Process& p)
    : framework::Producer(name,p) {}
  virtual void configure(const framework::config::Parameters& ps) final override;
  virtual void produce(framework::Event& event) final override;
};

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

void SimRecoilProjector::configure(const framework::config::Parameters& ps) {
}

void SimRecoilProjector::produce(framework::Event& event) {
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
    auto getTrajectory = [&hexReadout,this](std::vector<double> momentum,
                                       std::vector<float> position) {
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
    };

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

  event.add(
}

}

DECLARE_PRODUCER_NS(ecal,SimRecoilProjector);
