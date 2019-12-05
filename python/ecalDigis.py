#!/usr/bin/python

from LDMX.Framework import ldmxcfg

ecalDigis = ldmxcfg.Producer("ecalDigis","ldmx::EcalDigiProducer")

# Set the noise (in electrons) when the capacitance is 0.
ecalDigis.parameters["noiseIntercept"] = 900.

# Set the capacitative noise slope (electrons/pF)
ecalDigis.parameters["noiseSlope"] = 22.

# Set the capacitance per cell pad (pF)
ecalDigis.parameters["padCapacitance"] = 27.56

# set the readout threshold in multiples of RMS noise
ecalDigis.parameters["readoutThreshold"] = 4.
