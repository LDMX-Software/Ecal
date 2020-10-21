/**
 * @file EcalRecProducer.cxx
 * @brief Class that performs basic ECal reconstruction
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Ecal/EcalRecProducer.h"
#include "DetDescr/EcalHexReadout.h"

namespace ldmx {

    EcalRecProducer::EcalRecProducer(const std::string& name, Process& process) :
        Producer(name, process) {
    }

    EcalRecProducer::~EcalRecProducer() { }

    void EcalRecProducer::configure(Parameters& ps) {

        //collection names
        digiCollName_ = ps.getParameter<std::string>( "digiCollName" );
        digiPassName_ = ps.getParameter<std::string>( "digiPassName" );
        simHitCollName_  = ps.getParameter<std::string>("simHitCollName");
        simHitPassName_  = ps.getParameter<std::string>("simHitPassName");
        recHitCollName_ = ps.getParameter<std::string>("recHitCollName");

        layerWeights_ = ps.getParameter<std::vector<double>>( "layerWeights" );
        secondOrderEnergyCorrection_ = ps.getParameter<double>( "secondOrderEnergyCorrection" );

        mipSiEnergy_ = ps.getParameter<double>( "mipSiEnergy" );
        mV_          = ps.getParameter<double>( "mV" );

        auto hgcrocParams{ps.getParameter<Parameters>("hgcroc")};
        pedestal_    = hgcrocParams.getParameter<double>( "pedestal" );
        gain_        = hgcrocParams.getParameter<double>( "gain" );
        clockCycle_  = hgcrocParams.getParameter<double>( "clockCycle" );
        drainRate_   = hgcrocParams.getParameter<double>( "drainRate" );
        totMax_      = hgcrocParams.getParameter<double>( "totMax" );

    }

    void EcalRecProducer::produce(Event& event) {
        // Get the Ecal Geometry
        const EcalHexReadout& hexReadout = getCondition<EcalHexReadout>(EcalHexReadout::CONDITIONS_OBJECT_NAME);

        std::vector<EcalHit> ecalRecHits;
        auto ecalDigis = event.getObject<HgcrocDigiCollection>( digiCollName_ , digiPassName_ );
        int numDigiHits = ecalDigis.getNumDigis();
        //loop through digis
        for ( unsigned int iDigi = 0; iDigi < numDigiHits; iDigi++ ) {
            
            auto digi = ecalDigis.getDigi( iDigi );

            //ID from first digi sample
            //  assuming rest of samples have same ID
            EcalID id(digi.id());
            
            //ID to real space position
            double x,y,z;
            hexReadout.getCellAbsolutePosition( id , x , y , z );

            //TOA is the time of arrival with respect to the 25ns clock window
            //  TODO what to do if hit NOT in first clock cycle?
            double timeRelClock25 = digi.begin()->toa()*(clockCycle_/1024); //ns
            double hitTime = timeRelClock25;

            //get energy estimate from all digi samples
            double siEnergy(0.);

            /* debug printout
            std::cout << "Recon { "
                << "ID: " << id << ", "
                << "TOA: " << hitTime << "ns } ";
                */
            if ( digi.isTOT() ) {
                //TOT - number of clock ticks that pulse was over threshold
                //  this is related to the amplitude of the pulse approximately through a linear drain rate
                //  the amplitude of the pulse is related to the energy deposited

                int tdc = digi.tot();
                double tot = ( double(tdc)/4096. ) * totMax_;

                //convert the time over threshold into a total energy deposited in the silicon
                //  (time over threshold [ns]) * (rate of drain [mV/ns]) * (convert to energy [MeV/mV])
                siEnergy = tot * drainRate_ * mV_;

                /* debug printout
                std::cout << "TOT Mode -> "
                          << tdc << " TDC Counts -> " << tot << " ns -> "
                          << siEnergy << " MeV" << std::endl;
                 */
            } else {
                //ADC mode of readout
                //ADC - voltage measurement at a specific time of the pulse
                // Pulse Shape:
                //  p[0]/(1.0+exp(p[1](t-p[2]+p[3]-p[4])))/(1.0+exp(p[5]*(t-p[6]+p[3]-p[4])))
                //  p[0] = amplitude to be fit (TBD)
                //  p[1] = -0.345 shape parameter - rate of up slope
                //  p[2] = 70.6547 shape parameter - time of up slope relative to shape fit
                //  p[3] = 77.732 shape parameter - time of peak relative to shape fit
                //  p[4] = peak time to be fit (TBD)
                //  p[5] = 0.140068 shape parameter - rate of down slope
                //  p[6] = 87.7649 shape paramter - time of down slope relative to shape fit
                //These measurements can be used to fit the pulse shape if TOT is not available
                
                TH1F voltageMeasurements( "voltageMeasurements" , "voltageMeasurements" ,
                        10.*clockCycle_ , 0. , 10.*clockCycle_ );

                double maxMeas{0.};
                int numWholeClocks{0};
                for ( auto it = digi.begin(); it < digi.end(); it++) {
                    double voltage = (it->adc_t() - pedestal_)*gain_; //mV
                    if ( voltage > maxMeas ) maxMeas = voltage;
                    double time    = numWholeClocks*clockCycle_; //+ offestWithinClock; //ns
                    voltageMeasurements.Fill( time , voltage );
                }

                if ( false ) {
                    //fit the voltage measurements with the pulse function
                    //  would need to access the pulse function in HGCROC somehow
                    //voltageMeasurements.Fit( &pulseFunc_ , "QW" );
                    //get the silicon energy from the fitted voltage amplitude in mV
                    //siEnergy = (pulseFunc_.GetParameter( 0 ))*mV_;
                } else {
                    //just use the maximum measured voltage
                    siEnergy = (maxMeas)*mV_;
                }

                /* debug printout
                std::cout << "ADC Mode -> "
                          << siEnergy << " MeV" << std::endl;
                 */
            }
            
            //incorporate layer weights
            int layer = id.layer();
            double recHitEnergy = 
                ( (siEnergy / mipSiEnergy_ )*layerWeights_.at(layer)+siEnergy )*secondOrderEnergyCorrection_;

            //copy over information to rec hit structure in new collection
            EcalHit recHit;
            recHit.setID( id.raw() );
            recHit.setXPos( x );
            recHit.setYPos( y );
            recHit.setZPos( z );
            recHit.setAmplitude( siEnergy );
            recHit.setEnergy( recHitEnergy );
            recHit.setTime( hitTime );

            ecalRecHits.push_back( recHit );
        }

        if (event.exists( simHitCollName_, simHitPassName_ )) {
            //ecal sim hits exist ==> label which hits are real and which are pure noise
            auto ecalSimHits{event.getCollection<SimCalorimeterHit>( simHitCollName_, simHitPassName_ )};
            std::set<int> real_hits;
            for ( auto const& sim_hit : ecalSimHits ) real_hits.insert( sim_hit.getID() );
            for ( auto& hit : ecalRecHits ) hit.setNoise( real_hits.find(hit.getID()) == real_hits.end() );
        }

        //add collection to event bus
        event.add( recHitCollName_, ecalRecHits );
    }

}

DECLARE_PRODUCER_NS(ldmx, EcalRecProducer);