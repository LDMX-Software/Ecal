# This config is only meant to be run inside of the ecal digi pipeline test
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process( 'test_ecal_digis' )
# Set the maximum number of events
p.maxEvents = 2000
# Import the Ecal conditions and geometry
from LDMX.Ecal import ecal_hardcoded_conditions
from LDMX.Ecal import digi, ecal_trig_digi
p.outputFiles = ['ecal_digi_pipeline_test.root']
p.histogramFile = 'ecal_digi_pipeline_test_histo.root'
# Geometry provider
from LDMX.Ecal import EcalGeometry
ecalDigis = digi.EcalDigiProducer()
# Turn of noise hits
ecalDigis.hgcroc.noise = False
p.sequence = [
    ldmxcfg.Producer('fakeSimHits','ecal::test::EcalFakeSimHits','Ecal'),
    ecalDigis,
    ecal_trig_digi.EcalTrigPrimDigiProducer(),
    digi.EcalRecProducer(),
    ldmxcfg.Analyzer('checkEcalHits','ecal::test::EcalCheckEnergyReconstruction','Ecal'),
]

