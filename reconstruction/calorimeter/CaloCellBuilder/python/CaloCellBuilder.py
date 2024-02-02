
__all__ = ["CaloCellBuilder"]

from GaugiKernel            import Logger
from GaugiKernel.constants  import *
from GaugiKernel.macros     import *
from GaugiKernel            import GeV
from G4Kernel               import treatPropertyValue
from RootStreamBuilder      import recordable
import os


#
# Calo cell builder
#
class CaloCellBuilder( Logger ):


  def __init__( self, name, detector,
                      HistogramPath        = "Expert", 
                      HitsKey              = "Hits",
                      DoCrosstalk          = False,
                      OutputLevel          = 2,
                      XTAmpCapacitive      = 4.2,
                      XTAmpInductive       = 2.3,
                      XTAmpResistive       = 1.0,
                      ):

    Logger.__init__(self)
    self.__recoAlgs = []
    self.HistogramPath = HistogramPath
    self.OutputLevel = OutputLevel
    self.HitsKey = HitsKey
    self.__detector = detector
    self.__doCrosstalk = DoCrosstalk
    self.__XTAmpCapacitive  = XTAmpCapacitive
    self.__XTAmpInductive   = XTAmpInductive
    self.__XTAmpResistive   = XTAmpResistive

    


  #
  # Configure 
  #
  def configure(self):

    MSG_INFO(self, "Configure CaloCellBuilder.")
    from CaloCellBuilder import CaloCellMaker, CaloCellMerge, PulseGenerator, OptimalFilter, CalibrationTool, CrossTalk

    collectionKeys = []
    
    for samp in self.__detector.samplings:

      MSG_INFO(self, "Create new CaloCellMaker and dump all cells into %s collection", samp.CollectionKey)
      
      calib = CalibrationTool( "CalibrationTool", 
                                OutputLevel  = self.OutputLevel)
      
      pulse = PulseGenerator( "PulseGenerator", 
                              NSamples        = samp.Samples, 
                              ShaperFile      = samp.Shaper,
                              OutputLevel     = self.OutputLevel,
                              SamplingRate    = 25.0,
                              Pedestal        = 0.0,
                              DeformationMean = 0.0, 
                              DeformationStd  = 0.0,
                              NoiseMean       = 0.0,
                              NoiseStd        = 0.0,#samp.Noise,
                              StartSamplingBC = samp.StartSamplingBC )

      of= OptimalFilter("OptimalFilter",
                        OutputLevel=self.OutputLevel)
      

      alg = CaloCellMaker("CaloCellMaker_" + samp.CollectionKey, 
                            # input key
                            EventKey                = recordable( "EventInfo" ), 
                            HitsKey                 = recordable( self.HitsKey ),
                            # output key
                            CollectionKey           = samp.CollectionKey, 
                            # Hits grid configuration
                            EtaBins                 = samp.sensitive().EtaBins,
                            PhiBins                 = samp.sensitive().PhiBins,
                            ZMin                    = samp.volume().ZMin,
                            ZMax                    = samp.volume().ZMax,
                            RMin                    = samp.volume().RMin,
                            RMax                    = samp.volume().RMax,
                            Sampling                = samp.Sampling,
                            Noise                   = samp.Noise,
                            OFCa                    = samp.OFWeightsEnergy,
                            OFCb                    = samp.OFWeightsTime,
                            Segment                 = samp.sensitive().Segment,
                            Detector                = samp.Detector,
                            # Bunch crossing configuration
                            BunchIdStart            = samp.BunchIdStart,
                            BunchIdEnd              = samp.BunchIdEnd,
                            BunchDuration           = 25, #ns
                            # monitoring configuration
                            HistogramPath           = self.HistogramPath + '/' + samp.name(),
                            OutputLevel             = self.OutputLevel,
                            DetailedHistograms      = False, # Use True when debug with only one thread
                            )
      
      alg.Calibration = calib
      alg.PulseGenerator = pulse # for all cell
      alg.Tools = [of] # for each cel
      self.__recoAlgs.append( alg )
      collectionKeys.append( samp.CollectionKey )



    MSG_INFO(self, "Create CaloCellMerge and dump all cell collections into %s container", "Cells")
    # Merge all collection into a container and split between truth and reco
    mergeAlg = CaloCellMerge( "CaloCellMerge" , 
                              CollectionKeys  = collectionKeys,
                              CellsKey        = recordable("Cells"),
                              TruthCellsKey   = recordable("TruthCells"),
                              OutputLevel     = self.OutputLevel )

    self.__recoAlgs.append( mergeAlg )

    if self.__doCrosstalk: 
      MSG_INFO(self, "Create CrossTalk effect into Cell Collection")
      
      ofxt = OptimalFilter("OptimalFilter_XT",
                            OutputLevel    = self.OutputLevel)
      
      xt = CrossTalk("CrossTalk",
                            SigmaNoiseCut   = 2,
                            CollectionKeys  = collectionKeys,
                            XTCellsKey      = recordable("XTCells"),
                            CellsKey        = recordable("Cells"),
                            HistogramPath   = "",
                            OutputLevel     = self.OutputLevel,
                            XTAmpCapacitive = self.__XTAmpCapacitive,
                            XTAmpInductive  = self.__XTAmpInductive,
                            XTAmpResistive  = self.__XTAmpResistive,
                            )
      xt.Tools = [ ofxt]
      self.__recoAlgs.append( xt )


  def merge( self, acc ):
    # configure
    self.configure()
    for reco in self.__recoAlgs:
      acc+=reco 















