__all__ = ["RootStreamHITMaker"]

from GaugiKernel import Logger
from GaugiKernel.macros import *
from G4Kernel import treatPropertyValue


class RootStreamHITMaker( Logger ):

  __allow_keys = [
                  "InputEventKey",
                  "InputTruthKey",
                  "InputSeedsKey",
                  "InputHitsKey",
                  "OutputEventKey",
                  "OutputTruthKey",
                  "OutputSeedsKey",
                  "OutputHitsKey",
                  "OutputLevel", 
                  "NtupleName",
                  "OnlyRoI",
                  "EtaWindow",
                  "PhiWindow",
                  ]


  def __init__( self, name, **kw ): 
    
    Logger.__init__(self)
    import ROOT
    ROOT.gSystem.Load('liblorenzetti')
    from ROOT import RootStreamHITMaker
    self.__core = RootStreamHITMaker(name)

    for key, value in kw.items():
      self.setProperty( key,value )


  def core(self):
    return self.__core


  def setProperty( self, key, value ):
    if key in self.__allow_keys:
      setattr( self, '__' + key , value )
      self.core().setProperty( key, treatPropertyValue(value) )
    else:
      MSG_FATAL( self, "Property with name %s is not allow for %s object", key, self.__class__.__name__)

 
  def getProperty( self, key ):
    if key in self.__allow_keys:
      return getattr( self, '__' + key )
    else:
      MSG_FATAL( self, "Property with name %s is not allow for %s object", key, self.__class__.__name__)



