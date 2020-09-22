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
#include <memory> //for smart pointers

//----------//
//   LDMX   //
//----------//
#include "Framework/EventDef.h" 
#include "DetDescr/DetectorID.h"
#include "DetDescr/EcalID.h"
#include "DetDescr/EcalHexReadout.h"
#include "Framework/EventProcessor.h"

namespace ldmx {

    /**
     * @class EcalRecProducer
     * @brief Performs basic ECal reconstruction
     *
     * Reconstruction is done from the EcalDigi samples.
     * Some hard-coded parameters are used for position and energy calculation.
     */
    class EcalRecProducer : public Producer {

        public:

            /**
             * Constructor
             */
            EcalRecProducer(const std::string& name, Process& process);

            /**
             * Destructor
             */
            virtual ~EcalRecProducer();

            /**
             * Grabs configure parameters from the python config file.
             */
            virtual void configure(Parameters&);

            /**
             * Produce EcalHits and put them into the event bus using the
             * EcalDigis as input.
             *
             * This function unfolds the digi samples taken by the HGC ROC
             * and reconstructs their energy using knowledge of how
             * the chip operates and the position using EcalHexReadout.
             */
            virtual void produce(Event& event);

        private:

            /** Digi Collection Name to use as input */
            std::string digiCollName_;

            /** Digi Pass Name to use as input */
            std::string digiPassName_;

            /// Energy [MeV] deposited by a MIP in Si 0.5mm thick
            double mipSiEnergy_;

            /// Pedestal [ADC Counts] on below signal
            double pedestal_;

            /// Conversion from ADC counts to voltage [mV]
            double gain_;

            /// Length of clock cycle [ns]
            double clockCycle_;

            /// Rate that voltage drains off chip after saturation [mV/ns]
            double drainRate_;

            /// Maximum TOT measured by the chip [ns]
            double totMax_;

            /// Conversion from voltage [mV] to energy [MeV]
            double mV_;

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

            /**
             * Helper Instance of EcalHexReadout:
             *
             * performs real space postion <-> ID translation
             */
            std::unique_ptr<EcalHexReadout> ecalHexReadout_;

    };
}

#endif
