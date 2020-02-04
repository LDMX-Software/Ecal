/**
 * @file EcalRecProducer.cxx
 * @brief Class that performs basic ECal reconstruction
 * @author Tom Eichlersmith, University of Minnesota
 */

#include "Ecal/EcalRecProducer.h"

namespace ldmx {

    const double EcalRecProducer::MIP_SI_RESPONSE = 0.130; // MeV

    const double EcalRecProducer::DEFAULT_SECOND_ORDER_CORRECTION = 0.948;

    const std::vector<double> EcalRecProducer::DEFAULT_LAYER_WEIGHTS 
        = {1.641, 3.526, 5.184, 6.841,
        8.222, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775,
        8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 8.775, 12.642, 16.51,
        16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 16.51, 8.45}; 

    EcalRecProducer::EcalRecProducer(const std::string& name, Process& process) :
        Producer(name, process) {
    }

    EcalRecProducer::~EcalRecProducer() { }

    void EcalRecProducer::configure(const ParameterSet& ps) {

        digiCollName_ = ps.getString( "digiCollName" , "EcalDigis" );
        digiPassName_ = ps.getString( "digiPassName" , "" );

        layerWeights_ = ps.getVDouble( "layerWeights" , DEFAULT_LAYER_WEIGHTS );
        secondOrderEnergyCorrection_ = ps.getDouble( "secondOrderEnergyCorrection" , DEFAULT_SECOND_ORDER_CORRECTION );

        //if changing Ecal geometry, input geometry parameters here
        //  module radius - center-to-flat radius of a module in mm
        //  num cells wide - number of cells in center horizontal row
        //  gap between modules - air gap between adjacent flat edges of module
        //  z-position of front face of ECal [mm]
        //  list of positions of sensitive layers relative to front of ECal [mm]
        ecalHexReadout_ = std::make_unique<EcalHexReadout>();
    }

    void EcalRecProducer::produce(Event& event) {

        std::vector<EcalHit> ecalRecHits;
        EcalDigiCollection ecalDigis = event.getObject<EcalDigiCollection>( digiCollName_ , digiPassName_ );
        int numDigiHits = ecalDigis.getNumDigis();
        //loop through digis
        for ( unsigned int iDigi = 0; iDigi < numDigiHits; iDigi++ ) {
            
            //Right now, hard-coded to only use one sample in EcalDigiProducer
            //TODO: expand to multiple samples per digi
            EcalDigiSample sample = (ecalDigis.getDigi( iDigi )).at(0);
            int rawID = sample.rawID_;
            
            //ID to real space position
            double x,y,z;
            ecalHexReadout_->getCellAbsolutePosition( rawID , x , y , z );
            
            //get energy and time estimate from digi information
            double siEnergy, hitTime;
            
            //TODO: Energy estimate from N samples can (and should be) refined
            //TOA is the time of arrival with respect to the 25ns clock window
            double timeRelClock25 = sample.toa_*(25./1024); //ns
            hitTime = timeRelClock25;

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

            //TOT - number of clock ticks that pulse was over threshold
            //  this is related to the amplitude of the pulse through some convoluted relation using the pulse shape
            //  the amplitude of the pulse is related to the energy deposited
            //  TODO actually have mutliple samples instead of having adc_t_ count number of samples over threshold
            siEnergy = convertTOT( sample.adc_t_*1024 + sample.tot_ );
            
            //incorporate layer weights
            detID_.setRawValue( rawID );
            detID_.unpack();
            int layer = detID_.getLayerID();
            double recHitEnergy = 
                ( (siEnergy / MIP_SI_RESPONSE)*layerWeights_.at(layer)+siEnergy )*secondOrderEnergyCorrection_;

            //copy over information to rec hit structure in new collection
            EcalHit recHit;
            recHit.setID( rawID );
            recHit.setXPos( x );
            recHit.setYPos( y );
            recHit.setZPos( z );
            recHit.setAmplitude( siEnergy );
            recHit.setEnergy( recHitEnergy );
            recHit.setTime( hitTime );

            ecalRecHits.push_back( recHit );
        }

        //add collection to event bus
        event.add( "EcalRecHits", ecalRecHits );
    }

    double EcalRecProducer::convertTOT(const int tot) const {

        /**
         * Fit retrieved from a SimEDep vs TOT plot.
         * NOT physically motivated, will need to investigate further.
         *
         * Fit (for TOT > 3000):
         *    Function: expo ==> exp([0]+[1]*x)
         *    EDM=1.27359e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
         *    EXT PARAMETER                                   STEP         FIRST   
         *    NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
         *     1  Constant    -1.36729e+01   1.91947e-03   3.15252e-05  -6.27332e-01
         *     2  Slope        2.41246e-03   3.53378e-07   5.80385e-09  -1.99549e+03
         * Roughly flattens out (within uncertainty) to ~5e-2MeV when TOT < 3000
         */

        if ( tot > 3000 ) return exp( -1.36729e1 + 2.41246e-3*tot );
        
        return 5e-2;
    }
}

DECLARE_PRODUCER_NS(ldmx, EcalRecProducer);
