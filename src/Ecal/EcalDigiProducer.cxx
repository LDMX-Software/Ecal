/**
 * @file EcalDigiProducer.cxx
 * @brief Class that performs basic ECal digitization
 * @author Cameron Bravo, SLAC National Accelerator Laboratory
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Ecal/EcalDigiProducer.h"
#include "Framework/RandomNumberSeedService.h"

namespace ecal {

EcalDigiProducer::EcalDigiProducer(const std::string& name,
                                   framework::Process& process)
    : Producer(name, process) {
  // noise generator by default uses a Gausian model for noise
  //  i.e. It assumes the noise is distributed around a mean (setPedestal)
  //  with a certain RMS (setNoise) and then calculates
  //  how many hits should be generated for a given number of empty
  //  channels and a minimum readout value (setNoiseThreshold)
  noiseGenerator_ = std::make_unique<ldmx::NoiseGenerator>();
}

EcalDigiProducer::~EcalDigiProducer() {}

void EcalDigiProducer::configure(framework::config::Parameters& ps) {
  // settings of readout chip
  //  used  in actual digitization
  auto hgcrocParams = ps.getParameter<framework::config::Parameters>("hgcroc");
  hgcroc_ = std::make_unique<ldmx::HgcrocEmulator>(hgcrocParams);
  clockCycle_ = hgcrocParams.getParameter<double>("clockCycle");
  nADCs_ = hgcrocParams.getParameter<int>("nADCs");
  iSOI_ = hgcrocParams.getParameter<int>("iSOI");
  noise_ = hgcrocParams.getParameter<bool>("noise");

  // collection names
  inputCollName_ = ps.getParameter<std::string>("inputCollName");
  inputPassName_ = ps.getParameter<std::string>("inputPassName");
  digiCollName_ = ps.getParameter<std::string>("digiCollName");

  // physical constants
  //  used to calculate unit conversions
  MeV_ = ps.getParameter<double>("MeV");

  // Time -> clock counts conversion
  //  time [ns] * ( 2^10 / max time in ns ) = clock counts
  ns_ = 1024. / clockCycle_;

  // Configure generator that will produce noise hits in empty channels
  double readoutThreshold = ps.getParameter<double>("avgReadoutThreshold");
  double gain = ps.getParameter<double>("avgGain");
  double pedestal = ps.getParameter<double>("avgPedestal");
  // rms noise in mV
  noiseGenerator_->setNoise(
      hgcrocParams.getParameter<double>("noiseRMS"));
  // mean noise amplitude (if using Gaussian Model for the noise) in mV
  noiseGenerator_->setPedestal(gain*pedestal);  
  // threshold for readout in mV
  noiseGenerator_->setNoiseThreshold(gain*readoutThreshold);
}

void EcalDigiProducer::produce(framework::Event& event) {
  // Need to handle seeding on the first event
  if (!noiseGenerator_->hasSeed()) {
    const auto& rseed = getCondition<framework::RandomNumberSeedService>(
        framework::RandomNumberSeedService::CONDITIONS_OBJECT_NAME);
    noiseGenerator_->seedGenerator(
        rseed.getSeed("EcalDigiProducer::NoiseGenerator"));
  }
  if (noiseInjector_.get() == nullptr) {
    const auto& rseed = getCondition<framework::RandomNumberSeedService>(
        framework::RandomNumberSeedService::CONDITIONS_OBJECT_NAME);
    noiseInjector_ = std::make_unique<TRandom3>(
        rseed.getSeed("EcalDigiProducer::NoiseInjector"));
  }
  if (!hgcroc_->hasSeed()) {
    const auto& rseed = getCondition<framework::RandomNumberSeedService>(
        framework::RandomNumberSeedService::CONDITIONS_OBJECT_NAME);
    hgcroc_->seedGenerator(rseed.getSeed("EcalDigiProducer::HgcrocEmulator"));
  }

  hgcroc_->condition(getCondition<conditions::DoubleTableCondition>("EcalHgcrocConditions"));

  // Empty collection to be filled
  ldmx::HgcrocDigiCollection ecalDigis;
  ecalDigis.setNumSamplesPerDigi(nADCs_);
  ecalDigis.setSampleOfInterestIndex(iSOI_);

  std::set<unsigned int>
      filledDetIDs;  // detector IDs that already have a hit in them

  /******************************************************************************************
   * HGCROC Emulation on Simulated Hits
   *****************************************************************************************/
  // std::cout << "Sim Hits" << std::endl;
  // get simulated ecal hits from Geant4
  //  the class EcalHitIO in the SimApplication module handles the translation
  //  from G4CalorimeterHits to SimCalorimeterHits this class ensures that only
  //  one SimCalorimeterHit is generated per cell, but multiple "contributions"
  //  are still handled within SimCalorimeterHit
  auto ecalSimHits{event.getCollection<ldmx::SimCalorimeterHit>(
      inputCollName_, inputPassName_)};

  /* debug printout
  std::cout << "Energy to Voltage Conversion: " << MeV_ << " mV/MeV" <<
  std::endl;
   */

  std::vector<std::pair<double,double>> pulses_at_chip;
  for (auto const& simHit : ecalSimHits) {
    std::vector<double> voltages, times;
    for (int iContrib = 0; iContrib < simHit.getNumberOfContribs();
         iContrib++) {
      /* debug printout
      std::cout << simHit.getContrib(iContrib).edep << " MeV" << std::endl;
       */
      /**
       * HACK ALERT
       * The shifting of the time should _not_ be done this sloppily.
       * In reality, each chip has a set time phase that it samples at (relative
       * to target), so the time shifting should be at the emulator level.
       */
      pulses_at_chip.emplace_back(
        simHit.getContrib(iContrib).edep * MeV_,
        simHit.getContrib(iContrib).time  // global time (t=0ns at target)
          - simHit.getPosition().at(2) /
                299.702547  // shift light-speed particle traveling along z
      );
    }

    unsigned int hitID = simHit.getID();
    filledDetIDs.insert(hitID);

    /* debug printout
    std::cout << hitID << " "
        << simHit.getEdep() << std::endl;
     */
    // container emulator uses to write out samples and
    // transfer samples into the digi collection
    std::vector<ldmx::HgcrocDigiCollection::Sample> digiToAdd;
    if (hgcroc_->digitize(hitID, pulses_at_chip, digiToAdd)) {
      ecalDigis.addDigi(hitID, digiToAdd);
    }
  }

  /******************************************************************************************
   * Noise Simulation on Empty Channels
   *****************************************************************************************/
  if (noise_) {
    // std::cout << "Noise Hits" << std::endl;
    // put noise into some empty channels

    // geometry constants
    //  These are used in the noise generation so that we can randomly
    //  distribute the noise uniformly throughout the ECal channels.
    const auto& hexGeom = getCondition<ldmx::EcalHexReadout>(
        ldmx::EcalHexReadout::CONDITIONS_OBJECT_NAME);
    int nEcalLayers = hexGeom.getNumLayers();
    int nModulesPerLayer = hexGeom.getNumModulesPerLayer();
    int nCellsPerModule = hexGeom.getNumCellsPerModule();
    int numEmptyChannels = nEcalLayers * nModulesPerLayer * nCellsPerModule -
                           ecalDigis.getNumDigis();
    // noise generator gives us a list of noise amplitudes [mV] that randomly
    // populate the empty channels and are above the readout threshold
    auto noiseHitAmplitudes{
        noiseGenerator_->generateNoiseHits(numEmptyChannels)};
    std::vector<std::pair<double,double>> fake_pulse(1,{0.,0.});
    for (double noiseHit : noiseHitAmplitudes) {
      // generate detector ID for noise hit
      // making sure that it is in an empty channel
      unsigned int noiseID;
      do {
        int layerID = noiseInjector_->Integer(nEcalLayers);
        int moduleID = noiseInjector_->Integer(nModulesPerLayer);
        int cellID = noiseInjector_->Integer(nCellsPerModule);
        auto detID = ldmx::EcalID(layerID, moduleID, cellID);
        noiseID = detID.raw();
      } while (filledDetIDs.find(noiseID) != filledDetIDs.end());
      filledDetIDs.insert(noiseID);

      // get a time for this noise hit
      fake_pulse[0].second = noiseInjector_->Uniform(clockCycle_);

      // noise generator gives the amplitude above the readout threshold
      //  we need to convert it to the amplitude above the pedestal
      double gain = hgcroc_->gain(noiseID);
      fake_pulse[0].first = noiseHit +
                            gain * hgcroc_->readoutThreshold(noiseID) -
                            gain * hgcroc_->pedestal(noiseID);

      std::vector<ldmx::HgcrocDigiCollection::Sample> digiToAdd;
      if (hgcroc_->digitize(noiseID, fake_pulse, digiToAdd)) {
        ecalDigis.addDigi(noiseID, digiToAdd);
      }
    }  // loop over noise amplitudes
  }    // if we should do the noise

  event.add(digiCollName_, ecalDigis);

  return;
}  // produce

}  // namespace ecal

DECLARE_PRODUCER_NS(ecal, EcalDigiProducer);
