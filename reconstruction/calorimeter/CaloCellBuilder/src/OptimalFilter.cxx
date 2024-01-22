#include "OptimalFilter.h"
#include "CaloCell/CaloDetDescriptor.h"

using namespace Gaugi;


OptimalFilter::OptimalFilter( std::string name ) : 
  IMsgService(name),
  AlgTool()
{
  // declareProperty( "WeightsEnergy"    , m_ofweightsEnergy={}  );
  // declareProperty( "WeightsTime"      , m_ofweightsTime={}    );
  // declareProperty( "NoiseStd"         , m_noiseStd=0.0        );
  declareProperty( "OutputLevel"      , m_outputLevel=1       );
}

//!=====================================================================

OptimalFilter::~OptimalFilter()
{}

//!=====================================================================

StatusCode OptimalFilter::initialize()
{
  setMsgLevel(m_outputLevel);
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode OptimalFilter::finalize()
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode OptimalFilter::execute( SG::EventContext &/*ctx*/, Gaugi::EDM *edm ) const
{

  auto *cell = static_cast<xAOD::CaloDetDescriptor*>(edm);
  auto pulse = cell->pulse();
  float energy=0.0;
  float tau=0.0, EneTau=0.0;
  auto ofweightsEnergy=cell->OFCa();
  auto ofweightsTime=cell->OFCb();
  float noise=cell->noise();

  // Energy estimation
  if( ofweightsEnergy.size() != pulse.size() ){
    MSG_ERROR( "The ofweightsEnergy size its different than the pulse size." );
    return StatusCode::FAILURE;
  }else{
    for( unsigned sample=0; sample < pulse.size(); ++sample) 
      energy += pulse[sample]*ofweightsEnergy[sample];
  }

  // Time estimation
  if( ofweightsTime.size() != pulse.size() ){
    MSG_ERROR( "The ofweightsTime size its different than the pulse size." );
    return StatusCode::FAILURE;
  }
  else{
    for( unsigned sample=0; sample < pulse.size(); ++sample){
      EneTau += pulse[sample]*ofweightsTime[sample];
    }
    if (energy > 2*noise){
      tau = EneTau/energy;
    }
    else{
      tau = 0.0;
    }
    MSG_DEBUG(" Cell hash: "<< cell->hash()<<", sampling(noise): "<< cell->sampling()<<"("<< noise<<"), Energy(truth): "<< cell->edep() <<", Energy(OF2): " << energy <<", time(Truth): "<< cell->tof() << ", EneTau: "<< EneTau <<"time(OF): "<< tau );
  }
  
  cell->setE(energy);
  cell->setTau(tau);
  return StatusCode::SUCCESS;
}


