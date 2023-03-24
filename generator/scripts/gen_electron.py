#!/usr/bin/env python3
from GaugiKernel import LoggingLevel, Logger
from GaugiKernel import GeV
import argparse
import sys,os


mainLogger = Logger.getModuleLogger("pythia")
parser = argparse.ArgumentParser(description = '', add_help = False)
parser = argparse.ArgumentParser()

#
# Mandatory arguments
#

parser.add_argument('--nov','--numberOfEvents', action='store', dest='numberOfEvents', 
                    required = False, type=int, default=1,
                    help = "The number of events to be generated.")

parser.add_argument('--eventNumber', action='store', dest='eventNumber', 
                    required = False, default=0, type=int,
                    help = "The list of numbers per event.")

parser.add_argument('--runNumber', action='store', dest='runNumber', 
                    required = False, type=int, default = 0,
                    help = "The run number.")

                    
parser.add_argument('-o','--outputFile', action='store', dest='outputFile', required = True,
                    help = "The event file generated by pythia.")


parser.add_argument('--energy_min', action='store', dest='energy_min', required = False, type=float, default=-1,
                    help = "Energy min in GeV.")

parser.add_argument('--energy_max', action='store', dest='energy_max', required = False, type=float, default=-1,
                    help = "Energy max in GeV.")

parser.add_argument('-e', '--energy', action='store', dest='energy', required = False, type=float, default=-1,
                    help = "Energy in GeV.")


parser.add_argument('--eta', action='store', dest='eta', required = False, type=float, default=0.00,
                    help = "Eta position.")

parser.add_argument('--phi', action='store', dest='phi', required = False, type=float, default=1.52,
                    help = "Phi position.")

parser.add_argument('--doEtaRanged', action='store', dest='doEtaRanged', required = False, type=bool, default=False,
                    help = "Enable eta range.")

parser.add_argument('--eta_min', action='store', dest='eta_min', required = False, type=float, default=-2.5,
                    help = "Minimum Eta.")

parser.add_argument('--eta_max', action='store', dest='eta_max', required = False, type=float, default=2.5,
                    help = "Maximum Eta.")

#
# Pileup simulation arguments
#

parser.add_argument('--pileupAvg', action='store', dest='pileupAvg', required = False, type=int, default=0,
                    help = "The pileup average (default is zero).")

parser.add_argument('--bc_id_start', action='store', dest='bc_id_start', required = False, type=int, default=-21,
                    help = "The bunch crossing id start.")

parser.add_argument('--bc_id_end', action='store', dest='bc_id_end', required = False, type=int, default=4,
                    help = "The bunch crossing id end.")

parser.add_argument('--bc_duration', action='store', dest='bc_duration', required = False, type=int, default=25,
                    help = "The bunch crossing duration (in nanoseconds).")


#
# Extra parameters
#

parser.add_argument('--outputLevel', action='store', dest='outputLevel', required = False, type=int, default=0,
                    help = "The output level messenger.")

parser.add_argument('-s','--seed', action='store', dest='seed', required = False, type=int, default=0,
                    help = "The pythia seed (zero is the clock system)")




if len(sys.argv)==1:
  parser.print_help()
  sys.exit(1)

args = parser.parse_args()


try:

  from evtgen import Pythia8 
  from filters import SingleParticle, Particle
  from GenKernel import EventTape


  tape = EventTape( "EventTape", OutputFile = args.outputFile, RunNumber=args.runNumber )
  
  
  # Create the seed
  electron = SingleParticle( "Electron",
                             Pythia8("Generator", 
                                   Seed=args.seed, 
                                   EventNumber = args.eventNumber),
                             Eta          = args.eta,
                             Phi          = args.phi,
                             EnergyMin    = args.energy_min*GeV,
                             EnergyMax    = args.energy_max*GeV,
                             Energy       = args.energy*GeV,
                             Particle     = Particle.Electron, 
                             DoRangedEta  = args.doEtaRanged,
                             EtaMin       = args.eta_min,
                             EtaMax       = args.eta_max )
  
  tape+=electron

  if args.pileupAvg > 0:

    mb_file   = os.environ['LZT_PATH']+'/generator/guns/data/minbias_config.cmnd'

    pileup = Pileup("Pileup",
                   Pythia8("Generator", File=mb_file, Seed=args.seed),
                   EtaMax         = 3.2,
                   Select         = 2,
                   PileupAvg      = args.pileupAvg,
                   BunchIdStart   = args.bc_id_start,
                   BunchIdEnd     = args.bc_id_end,
                   OutputLevel    = args.outputLevel,
                   DeltaEta       = 0.22,
                   DeltaPhi       = 0.22,
                  )

    tape+=pileup
  
 
  # Run!
  tape.run(args.numberOfEvents)

  sys.exit(0)
except  Exception as e:
  print(e)
  sys.exit(1)
