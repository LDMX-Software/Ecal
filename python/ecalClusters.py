"""Configuration for the EcalClusterProducer

Examples
--------
    from LDMX.Ecal.ecalClusters import ecalClusters
    p.sequence.append( ecalClusters )
"""

from LDMX.Framework import ldmxcfg

class EcalClusterProducer(ldmxcfg.Producer) :
    """Configure the clustering"""

    def __init__(self,name) :
        super().__init__(name,"ldmx::EcalClusterProducer","Ecal")

        self.cutoff = 10.
        self.seedThreshold = 100.0 #MeV

        # Pass name for ecal digis
        self.digisPassName = "recon"

        # Minimum number of hits in a cluster to save it to the collection
        self.nHitsMin = 3

        # Name of the algo to save to the root file 
        self.algoName = "MyClusterAlgo"

        # Name of the cluster collection to make
        self.clusterCollName = "ecalClusters"

        # Name of the cluster algo collection to make
        self.algoCollName = "ClusterAlgoResult"

        from LDMX.DetDescr import EcalHexReadout
        self.hexReadout = EcalHexReadout.EcalHexReadout()
        
ecalClusters = EcalClusterProducer("ecalClusters")