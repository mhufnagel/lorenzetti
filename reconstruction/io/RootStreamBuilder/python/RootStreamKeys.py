

__all__ = ['recordable']



def recordable( key ):

  keys = [
            # CaloCellMaker
            "Collection_EM1",
            "Collection_EM2",
            "Collection_EM3",
            "Collection_HAD1",
            "Collection_HAD2",
            "Collection_HAD3",
            "Hits",
            # CaloCellMerge
            "Cells",
            "TruthCells",
            # CrossTalk
            "XTCells",
            "XTClusters",
            "XTRings",
            # CaloClusterMaker
            "EventInfo",
            "Clusters",
            "TruthClusters",
            "Particles",
            "TruthParticles",
            "Seeds",
            # CaloRingsBuilder
            "Rings",
            "TruthRings",
            ]

  return key
