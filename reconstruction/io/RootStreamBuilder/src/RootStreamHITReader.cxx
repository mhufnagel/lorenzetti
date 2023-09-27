

#include "CaloHit/CaloHitContainer.h"
#include "EventInfo/EventInfoContainer.h"
#include "TruthParticle/TruthParticleContainer.h"
#include "TruthParticle/ParticleSeedContainer.h"

#include "CaloHit/CaloHitConverter.h"
#include "EventInfo/EventInfoConverter.h"
#include "TruthParticle/TruthParticleConverter.h"
#include "TruthParticle/ParticleSeedConverter.h"

#include "RootStreamHITReader.h"
#include "GaugiKernel/EDM.h"


using namespace SG;
using namespace Gaugi;



RootStreamHITReader::RootStreamHITReader( std::string name ) : 
  IMsgService(name),
  Algorithm()
{
  declareProperty( "InputFile"          , m_inputFile=""                    );
  declareProperty( "EventKey"           , m_eventKey="EventInfo"            );
  declareProperty( "TruthKey"           , m_truthKey="Particles"            );
  declareProperty( "SeedsKey"           , m_seedsKey="Seeds"                );
  declareProperty( "HitsKey"            , m_hitsKey="Hits"                  );
  declareProperty( "OutputLevel"        , m_outputLevel=1                   );
  declareProperty( "NtupleName"         , m_ntupleName="CollectionTree"     );
}

//!=====================================================================

RootStreamHITReader::~RootStreamHITReader()
{}

//!=====================================================================

StatusCode RootStreamHITReader::initialize()
{
  CHECK_INIT();
  setMsgLevel(m_outputLevel);
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode RootStreamHITReader::finalize()
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode RootStreamHITReader::bookHistograms( EventContext &ctx ) const
{
  MSG_INFO("Reading file " << m_inputFile);
  auto store = ctx.getStoreGateSvc();
  TFile *file = new TFile(m_inputFile.c_str(), "read");
  store->decorate( "events", file );
  return StatusCode::SUCCESS; 
}

//!=====================================================================

StatusCode RootStreamHITReader::pre_execute( EventContext &/*ctx*/ ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode RootStreamHITReader::execute( EventContext &/*ctx*/, const G4Step * /*step*/ ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode RootStreamHITReader::execute( EventContext &ctx, int evt ) const
{
  return deserialize( evt, ctx );
}

//!=====================================================================

StatusCode RootStreamHITReader::post_execute( EventContext &/*ctx*/ ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode RootStreamHITReader::fillHistograms( EventContext &ctx ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================


StatusCode RootStreamHITReader::deserialize( int evt, EventContext &ctx ) const
{
  std::vector<xAOD::CaloHit_t           > *collection_hits       = nullptr;
  std::vector<xAOD::EventInfo_t         > *collection_event      = nullptr;
  std::vector<xAOD::TruthParticle_t     > *collection_truth      = nullptr;
  std::vector<xAOD::ParticleSeed_t      > *collection_seeds      = nullptr;

  MSG_DEBUG( "Link all branches..." );
  
  auto store = ctx.getStoreGateSvc();
  TFile *file = (TFile*)store->decorator("events");
  TTree *tree = (TTree*)file->Get(m_ntupleName.c_str());

  InitBranch( tree, ("EventInfoContainer_"+m_eventKey).c_str()     , &collection_event     );
  InitBranch( tree, ("TruthParticleContainer_"+m_truthKey).c_str() , &collection_truth     );
  InitBranch( tree, ("ParticleSeedContainer_"+m_seedsKey).c_str()  , &collection_seeds     );
  InitBranch( tree, ("CaloHitContainer_"+m_hitsKey).c_str()        , &collection_hits      );

  tree->GetEntry( evt );


  { // deserialize TruthParticle
    SG::WriteHandle<xAOD::TruthParticleContainer> container(m_truthKey, ctx);
    container.record( std::unique_ptr<xAOD::TruthParticleContainer>(new xAOD::TruthParticleContainer()));

    xAOD::TruthParticleConverter cnv;
    for( auto& par_t : *collection_truth)
    {
      xAOD::TruthParticle  *par=nullptr;
      cnv.convert(par_t, par);
      MSG_INFO( "Particle in eta = " << par->eta() << ", phi = " << par->phi() << ", pdgID = "<< par->pdgid());
      container->push_back(par);
    }
  }

  { // deserialize ParticleSeed
    SG::WriteHandle<xAOD::ParticleSeedContainer> container(m_seedsKey, ctx);
    container.record( std::unique_ptr<xAOD::ParticleSeedContainer>(new xAOD::ParticleSeedContainer()));

    xAOD::ParticleSeedConverter cnv;
    for( auto& par_t : *collection_seeds)
    {
      xAOD::ParticleSeed  *par=nullptr;
      cnv.convert(par_t, par);
      MSG_INFO( "Particle seeded in eta = " << par->eta() << ", phi = " << par->phi() << ", et_tot = "<< par->ettot());
      container->push_back(par);
    }
  }



  { // deserialize EventInfo

    SG::WriteHandle<xAOD::EventInfoContainer> container(m_eventKey, ctx);
    container.record( std::unique_ptr<xAOD::EventInfoContainer>(new xAOD::EventInfoContainer()));
    xAOD::EventInfo  *event=nullptr;
    xAOD::EventInfoConverter cnv;
    cnv.convert(  collection_event->at(0), event);
    MSG_INFO( "EventNumber = " << event->eventNumber() << ", Avgmu = " << event->avgmu());
    container->push_back(event);
  }
  

  {
    SG::WriteHandle<xAOD::CaloHitContainer> container(m_hitsKey, ctx);
    container.record( std::unique_ptr<xAOD::CaloHitContainer>(new xAOD::CaloHitContainer()));
    float etot=0;
    for( auto &hit_t : *collection_hits )
    {
      xAOD::CaloHit *hit = nullptr;
      xAOD::CaloHitConverter cnv;
      cnv.convert(hit_t, hit); // alloc memory
      container->push_back(hit);
      etot+=hit->edep();
    }

    MSG_DEBUG("Container hit size is " << container->size() << " and total energy " << etot << " MeV " );
  }



  delete collection_hits      ;
  delete collection_event     ;
  delete collection_truth     ;
  delete collection_seeds     ;
  return StatusCode::SUCCESS;
 
}

//!=====================================================================

template <class T>
void RootStreamHITReader::InitBranch(TTree* fChain, std::string branch_name, T* param) const
{
  std::string bname = branch_name;
  if (fChain->GetAlias(bname.c_str()))
     bname = std::string(fChain->GetAlias(bname.c_str()));

  if (!fChain->FindBranch(bname.c_str()) ) {
    MSG_WARNING( "unknown branch " << bname );
    return;
  }
  fChain->SetBranchStatus(bname.c_str(), 1.);
  fChain->SetBranchAddress(bname.c_str(), param);
}

//!=====================================================================

