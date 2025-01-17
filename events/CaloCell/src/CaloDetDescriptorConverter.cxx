#include "CaloCell/CaloDetDescriptorConverter.h"
//#include "G4Kernel/macros.h"


using namespace xAOD;



bool CaloDetDescriptorConverter::convert( const CaloDetDescriptor *descriptor, CaloDetDescriptor_t &descriptor_t, int link )
{

  if(descriptor){
    descriptor_t.sampling    = (int)descriptor->sampling();
    descriptor_t.detector    = (int)descriptor->detector();
    descriptor_t.eta         = descriptor->eta();
    descriptor_t.phi         = descriptor->phi();
    descriptor_t.deta        = descriptor->deltaEta(); 
    descriptor_t.dphi        = descriptor->deltaPhi();
    descriptor_t.e           = descriptor->e();
    descriptor_t.tau         = descriptor->tau(); // estimated time of flight  (with respect to distance and collision clock)
    descriptor_t.bcid_start  = descriptor->bcid_start();
    descriptor_t.bcid_end    = descriptor->bcid_end();
    descriptor_t.bc_duration = descriptor->bc_duration();
    descriptor_t.hash        = descriptor->hash();
    descriptor_t.pulse       = descriptor->pulse();
    descriptor_t.cell_link   = link; // cross link to cell

    for ( int bcid = descriptor->bcid_start();  bcid <= descriptor->bcid_end(); ++bcid )
    {
      descriptor_t.edep_per_bunch.push_back( descriptor->edep(bcid) );
      descriptor_t.tof_per_bunch.push_back( descriptor->tof(bcid)   );
    }
    descriptor_t.edep   = descriptor->edep(); // truth energy into the bunch crossing zero
    descriptor_t.tof    = descriptor->tof(); // 'truth' time of flight in the cell at BC=0

    return true;
  }

  return false;
}



bool CaloDetDescriptorConverter::convert( const CaloDetDescriptor_t &descriptor_t, CaloDetDescriptor *&descriptor )
{

  descriptor = new xAOD::CaloDetDescriptor( descriptor_t.eta, 
                                            descriptor_t.phi,
                                            descriptor_t.deta,
                                            descriptor_t.dphi,
                                            descriptor_t.hash,
                                            (CaloSampling)descriptor_t.sampling, 
                                            (Detector)descriptor_t.detector,
                                            descriptor_t.bc_duration,
                                            descriptor_t.bcid_start,
                                            descriptor_t.bcid_end );

  descriptor->setE(descriptor_t.e); // estimated energy from OF
  descriptor->setTau(descriptor_t.tau); // 
  descriptor->setPulse( descriptor_t.pulse); // pulse from generator
  
  int pos=0;
  for ( int bcid = descriptor->bcid_start();  bcid <= descriptor->bcid_end(); ++bcid )
  {
    descriptor->edep( bcid, descriptor_t.edep_per_bunch.at(pos) ); // truth energy for each bunch crossing
    descriptor->tof ( bcid, descriptor_t.tof_per_bunch.at(pos)  ); //
    pos++;
  }
  return true;
}