#include "CrossTalk.h"

#include "CaloCell/CaloDetDescriptor.h"
#include "CaloCell/CaloCellContainer.h"
#include "CaloCell/CaloDetDescriptorCollection.h"
#include "CaloCell/CaloCellConverter.h"
#include "CaloCell/CaloDetDescriptorConverter.h"

#include "G4Kernel/CaloPhiRange.h"
#include "G4Kernel/constants.h"

#include "G4SystemOfUnits.hh"

#include <map>


using namespace Gaugi;


CrossTalk::CrossTalk( std::string name ) : 
  IMsgService(name),
  Algorithm()
  // AlgTool()
{
  declareProperty( "SigmaNoiseCut"    , m_sigmaNoiseCut=0                              );
  declareProperty( "CollectionKeys"   , m_collectionKeys={}                            ); // input
  declareProperty( "XTCellsKey"       , m_xtcellsKey="XTCells"                         ); // input
  declareProperty( "CellsKey"         , m_cellsKey="Cells"                             ); // input
  declareProperty( "HistogramPath"    , m_histPath="/CrossTalkSimulator"               );
  declareProperty( "OutputLevel"      , m_outputLevel=1                                );
  declareProperty( "XTAmpCapacitive"  , m_AmpXt_C=4.2                                  ); // in %
  declareProperty( "XTAmpInductive"   , m_AmpXt_L=2.3                                  ); // in %
  declareProperty( "XTAmpResistive"   , m_AmpXt_R=1.0                                  ); // in %
}

//!=====================================================================

// CrossTalk::~CrossTalk()
// {;}

//!=====================================================================

StatusCode CrossTalk::initialize()
{
  CHECK_INIT();
  setMsgLevel(m_outputLevel);

  // initialize tools
  for ( auto tool : m_toolHandles )
  {
    if (tool->initialize().isFailure() )
    {
      MSG_FATAL( "It's not possible to iniatialize " << tool->name() << " tool." );
    }
  }
  return StatusCode::SUCCESS;
}
//!=====================================================================

void CrossTalk::push_back( Gaugi::AlgTool* tool )
{
  m_toolHandles.push_back(tool);
}

//!=====================================================================

StatusCode CrossTalk::pre_execute( SG::EventContext &/*ctx*/ ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CrossTalk::execute( SG::EventContext &/*ctx*/ , const G4Step * /*step*/ ) const
{
  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CrossTalk::execute( SG::EventContext &ctx , int /*evt*/ ) const
{
  MSG_INFO("Executing CrossTalk module...");

  std::vector < float > samples_xtalk_ind    ;
  std::vector < float > samples_xtalk_cap    ;
  std::vector < float > samples_signal       ;
  std::vector < float > samples_signal_xtalk ;

  SG::ReadHandle<xAOD::CaloCellContainer> container   (m_cellsKey, ctx);

  MSG_INFO( "Creating reco XT cells containers with key " << m_xtcellsKey);
  SG::WriteHandle<xAOD::CaloCellContainer> xtContainer( m_xtcellsKey , ctx );
  xtContainer.record( std::unique_ptr<xAOD::CaloCellContainer>(new xAOD::CaloCellContainer()) );

  // create xt energy excess container for further corrections
  SG::WriteHandle<xAOD::CaloDetDescriptorCollection> xtEneExcess( "XTDescriptorEneExcess" , ctx );
  xtEneExcess.record( std::unique_ptr<xAOD::CaloDetDescriptorCollection>(new xAOD::CaloDetDescriptorCollection()) );
  
  MSG_DEBUG("Before execution: collection.size: "<< container.ptr()->size() << ", xtCollection.size(): "<< xtContainer->size());

  // loop over ordinary cell container
  for (const auto cell : **container.ptr() ){
    
    auto xtdescriptor = new xAOD::CaloDetDescriptor( *(cell->descriptor()) ); // get content to a new object pointer, to be "crosstalked"
    
    auto pulseBefore  = cell->descriptor()->pulse();
    // auto energyBefore = cell->descriptor()->e();
    // auto timeBefore   = cell->descriptor()->tau();

    // Step 1: check if we need to apply cx method for current cell. Only for cells higher than
    // min energy. Here, lets use the truth energy from the main bunch crossing.
    bool bCrossTalkConditions = ( !(xtdescriptor->edep() < m_sigmaNoiseCut*xtdescriptor->noise()) && !(xtdescriptor->pulse().size() == 0) && !((xtdescriptor->sampling() != 3) && (xtdescriptor->sampling()  != 12)) );

    // ------------------------------------------------------------------------------
    // If there IS xtalk conditions, apply XT model to cell 1st neighbors,
    // change the current cell pulse, then add cell to new XT cell container.
    // -------------------------------------------------------------------------------
    if (bCrossTalkConditions){

      MSG_DEBUG("bCrossTalkConditions=true: Sampling/Detector "<< xtdescriptor->sampling() <<"/"<< xtdescriptor->detector() <<", hash "<< xtdescriptor->hash() <<", nsamples " << xtdescriptor->pulse().size() << ", truthEne "<< xtdescriptor->edep() << ", ene " << xtdescriptor->e() << ", eta/phi "<< xtdescriptor->eta() << "/"<< xtdescriptor->phi() );

      // Step 2: build a 3x3 window around the central cell.
      //    Since this is a cell candidate, lets take all cells around this cells using a 3x3 window.
      //    First lets retrieve the full container in memory (not const objects inside of the collection)

      std::vector<const xAOD::CaloCell*> cells_around;

      // loop over ordinary cell container
      for (auto neighborCell : **container.ptr() ){
        const xAOD::CaloDetDescriptor *neighborDescriptor = neighborCell->descriptor();

        if ( neighborDescriptor->pulse().size() == 0) continue; // protection: if there is no pulseShape, skip that cell. 
        if ( xtdescriptor->sampling() != neighborDescriptor->sampling() ) continue;  // cells_around must belong to the same sampling of central_cell
        if ( xtdescriptor->hash() == neighborDescriptor->hash()) continue; // central_cell must not belong to cells_around
        
        // build a 3x3 window around the central cell
        float diffEta = std::abs( cell->eta() - neighborCell->eta() );
        float diffPhi = std::abs( CaloPhiRange::fix( cell->phi() - neighborCell->phi() ) );

        if( diffEta <= 3*cell->deltaEta()/2 && diffPhi <= 3*cell->deltaPhi()/2 ){
          cells_around.push_back(neighborCell); // cells that will have their charges 'leaked' to center cell, on XT computation.
        }
      }

      // Step 3: Loop over cells_around to extract xtalk effect from the central_cell surroundings.
      std::vector<float> final_xt_pulse(5);

      for (auto cellXT : cells_around){
    
        auto pulseCellXT = cellXT->descriptor()->pulse();
        // auto *excessEneDescriptor = new xAOD::CaloDetDescriptor( *(cellXT->descriptor()) ); // for cluster charge conservation correction

        std::vector<float> neighbor_xt_pulse;

        float distorted_sample_ind=0, distorted_sample_cap=0;
        
        // case 1: diagonal from central cell
        if (cellXT->eta() != cell->eta() && cellXT->phi() != cell->phi()){
          
          for (unsigned samp_index=0; samp_index<5; ++samp_index){
            distorted_sample_ind = XTalkTF( pulseCellXT[samp_index],samp_index,true, true);
            distorted_sample_cap = 0; // there is no capacitive cross-talk effect in the cell diagonal

            samples_xtalk_ind.push_back(distorted_sample_ind); // histogram
            samples_xtalk_cap.push_back(distorted_sample_cap);  // histogram

            neighbor_xt_pulse.push_back(distorted_sample_ind + distorted_sample_cap);
          }
        }
        else {
          // case 2: is inside central cross position        
          for (int samp_index=0; samp_index<5; samp_index++){
            distorted_sample_ind = XTalkTF( pulseCellXT[samp_index],samp_index,true, true);
            distorted_sample_cap = XTalkTF( pulseCellXT[samp_index],samp_index,true, false);

            samples_xtalk_ind.push_back(distorted_sample_ind); // histogram
            samples_xtalk_cap.push_back(distorted_sample_cap); // histogram

            neighbor_xt_pulse.push_back(distorted_sample_ind + distorted_sample_cap);
          }
        }
        // sum all xtalk effects around center cell
        for (int i=0; i<5; i++){
          final_xt_pulse[i] += neighbor_xt_pulse[i];
        }

        // *** Energy conservation correction: ***
        // Trye get the current neighbor cell descriptor from 'excess energy container.
        // Here, if that descriptor hash already exists in the container, then integrate the new computed energy to its 'Pulse'.
        xAOD::CaloDetDescriptor *excessEneDescriptor=nullptr;
        if ( xtEneExcess->retrieve( cellXT->descriptor()->hash(), excessEneDescriptor ) ){
          std::vector<float> integratedEne(5);
          std::vector<float> oldPulse = excessEneDescriptor->pulse();
          for (int i=0; i<5; i++){
            integratedEne[i] = oldPulse[i] + neighbor_xt_pulse[i];
          }

          MSG_DEBUG("(XTChargeConservation) Hash "<<cellXT->descriptor()->hash() << " exists on container, integrating its energy from "<< excessEneDescriptor->pulse() <<", to " << integratedEne);
          excessEneDescriptor->setPulse(integratedEne);
        }
        else{ // if hash do not exist in container, add new descriptor to it.
          auto *newExcessEneDescriptor = new xAOD::CaloDetDescriptor( *(cellXT->descriptor()) );
          newExcessEneDescriptor->setPulse( neighbor_xt_pulse );

          if (!xtEneExcess->insert(newExcessEneDescriptor->hash() , newExcessEneDescriptor)){
            MSG_FATAL("(XTChargeConservation) Descriptor is unique and it's not possible to insert new descriptor into collection!");
            return StatusCode::FAILURE;
          }
          MSG_DEBUG("(XTChargeConservation) New cell added with hash "<< newExcessEneDescriptor->hash() << " and signal " << newExcessEneDescriptor->pulse());
        }
        // ***************************************
      } // end-for in cells_around


      // Step 4: add total pulse distortion from neighbor cells into the central cell of the 3x3 window.
      auto centralCellPulse = xtdescriptor->pulse(); 

      for (int i=0; i<5; i++){
        samples_signal.push_back(centralCellPulse[i]); // add to fillHistograms
        centralCellPulse[i] = centralCellPulse[i] + final_xt_pulse[i]; 
        samples_signal_xtalk.push_back(centralCellPulse[i]); // add to fillHistograms
      }

      // Step 5: change pulse value of central cell of the 3x3 window with adjacent xtalk effects.
      xtdescriptor->setPulse(centralCellPulse);
      
      // // Step xxx: Correct for charge excess
      // xAOD::CaloDetDescriptor *excessEneDescriptor=nullptr;
      // if ( !xtEneExcess->retrieve( xtdescriptor->hash() , excessEneDescriptor) ){
      //   MSG_FATAL("Cannot find hash "<< xtdescriptor->hash()  << " on EnergyExcessContainer!");
      // }

      // Step 6: Call for Estimation Methods tool (or any other tool applied into cells, AFTER pulse generation.)
      // for ( auto tool : m_toolHandles )
      // {
      //   // digitalization
      //   if( tool->execute( ctx, xtdescriptor ).isFailure() ){
      //     MSG_ERROR( "It's not possible to execute the tool with name " << tool->name() );
      //     return StatusCode::FAILURE;
      //   }
      // }

      // auto pulseAfter   = xtdescriptor->pulse();//centralCellPulse;
      // auto energyAfter  = xtdescriptor->e();
      // auto timeAfter    = xtdescriptor->tau();

      // MSG_DEBUG(" e: "<< energyBefore <<" t: "<< timeBefore << "  Pulse before: " << pulseBefore[0] << "   "<< pulseBefore[1] << "   "<< pulseBefore[2] << "   "<< pulseBefore[3] << "   "<< pulseBefore[4]);
      // MSG_DEBUG(" e_xt: "<< energyAfter << " t_xt: "<< timeAfter <<"  Pulse after XT: " << pulseAfter[0] << "   "<< pulseAfter[1] << "   "<< pulseAfter[2] << "   "<< pulseAfter[3] << "   "<< pulseAfter[4]);
      MSG_DEBUG("Cell "<< cell->descriptor()->hash() <<", sampling "<< cell->descriptor()->sampling() <<", pulse() = "<< cell->descriptor()->pulse() << ", edep/tof= "<< cell->descriptor()->edep() <<"/"<<cell->descriptor()->tof()  <<", e/tau=" << cell->descriptor()->tau() << "/"<<cell->descriptor()->e() );
      // MSG_DEBUG("XTCell "<< xtdescriptor->hash() <<", sampling "<< xtdescriptor->sampling() <<", pulse() = "<< xtdescriptor->pulse() << ", tof= "<<xtdescriptor->tof() <<", tau=" << xtdescriptor->tau() );

      samples_xtalk_ind.clear();
      samples_xtalk_cap.clear();
      samples_signal.clear();
      samples_signal_xtalk.clear(); 
    }

    // ------------------------------------------------------------------------------
    // If there is NO xtalk conditions, add the cell normally into new XT Container.
    //  Look, here, the current descriptor hasn't been changed.
    // -------------------------------------------------------------------------------
    // Setup the caloCell to XT cell container (Copy from original CellContainer)
    auto xtcell = new xAOD::CaloCell();

    xtcell->setEta( xtdescriptor->eta() );
    xtcell->setPhi( xtdescriptor->phi() );
    xtcell->setDeltaEta( xtdescriptor->deltaEta() );
    xtcell->setDeltaPhi( xtdescriptor->deltaPhi() );
    xtcell->setE( xtdescriptor->e() ); // Estimated energy from OF/COF
    xtcell->setTau( xtdescriptor->tau() ); // Estimated time from OF/COF
    xtcell->setEt( xtdescriptor->e() / std::cosh( xtdescriptor->eta() ) );
    xtcell->setDescriptor( xtdescriptor );
    xtContainer->push_back( xtcell ); //add CaloCell to XTcontainer  

  }

  MSG_DEBUG("collection.size: "<< container.ptr()->size() << ", xtCollection.size(): "<< xtContainer->size());

  // Extra step i: Correct for charge excess
  // for (auto excessEneDescriptor : *xtEneExcess){
  for (auto xtcell : **xtContainer ){

    // xAOD::CaloCell *xtcell=nullptr;
    xAOD::CaloDetDescriptor *descriptor = const_cast<xAOD::CaloDetDescriptor*>(xtcell->descriptor());
    xAOD::CaloDetDescriptor *excessDescriptor = nullptr;

    if ( !xtEneExcess->retrieve( descriptor->hash() , excessDescriptor) ){
      // MSG_ERROR("Cannot find hash "<< descriptor->hash()  << " on EnergyExcess Collection!");
      continue;
    }

    std::vector<float> correctedPulse(5);
    std::vector<float> xtCellPulse = descriptor->pulse();
    std::vector<float> excessPulse = excessDescriptor->pulse();

    for (int i=0; i<5; i++){
      correctedPulse[i] = xtCellPulse[i] - excessPulse[i];
    }
    descriptor->setPulse(correctedPulse);

    // Extra Step ii: Call for Estimation Methods tool
    for ( auto tool : m_toolHandles )
    {
      // digitalization
      if( tool->execute( ctx, descriptor ).isFailure() ){
        MSG_ERROR( "It's not possible to execute the tool with name " << tool->name() );
        return StatusCode::FAILURE;
      }
    }

    const_cast<xAOD::CaloCell*>(xtcell)->setE(   descriptor->e() ); // Estimated energy from OF/COF
    const_cast<xAOD::CaloCell*>(xtcell)->setTau( descriptor->tau() ); // Estimated time from OF/COF
    const_cast<xAOD::CaloCell*>(xtcell)->setEt(  descriptor->e() / std::cosh( descriptor->eta() ) );

    MSG_DEBUG("(ChargeCorrection) Hash "<< descriptor->hash() <<" corrected pulse from "<< xtCellPulse << " to "<< descriptor->pulse()<<", e/tau = " << descriptor->e()<< "/"<<descriptor->tau());
    MSG_DEBUG("XTCell "<< descriptor->hash() <<", sampling "<< descriptor->sampling() <<", pulse() = "<< descriptor->pulse() << ", edep/tof= "<< descriptor->edep() <<"/"<< descriptor->tof() );
  }

  return StatusCode::SUCCESS;
  
}
//!=====================================================================

StatusCode CrossTalk::finalize()
{
  for ( auto tool : m_toolHandles )
  {
    if (tool->finalize().isFailure() )
    {
      MSG_ERROR( "It's not possible to finalize " << tool->name() << " tool." );
    }
  }

  return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CrossTalk::post_execute( SG::EventContext &/*ctx*/ ) const
{
 
  return StatusCode::SUCCESS;
}

float CrossTalk::XTalkTF(float sample, int samp_index, bool diagonal, bool inductive) const
{

  float BaseAmpXTc = m_AmpXt_C/100*sample ;
  float BaseAmpXTl = m_AmpXt_L/100*sample ;
  // float BaseAmpXTr = m_AmpXt_R*sample ;
  float XTcSamples = BaseAmpXTc * XTalk       (25*(samp_index+1) , false ); //+ delayPerCell[cell] + m_tau_0, false ) ) ;
  float XTlSamples = BaseAmpXTl * XTalk       (25*(samp_index+1) , false ); //+ delayPerCell[cell] + m_tau_0, false ) ) ;
  // float XTrSamples = BaseAmpXTr * CellFunction(25*(samp_index+1) , false ); //+ delayPerCell[cell] + m_tau_0, false ) ) ;

  if (diagonal && inductive){
    // ind_part = XTlSamples;
    return XTlSamples;
  }
  else{
    // return XTcSamples + XTlSamples;
    if (inductive){
      return XTlSamples;
      }
    else{
      return XTcSamples;
    }
  }
  
  // SampClusNoise.push_back( noise->Gaus(0, 2) ) ;

}

double CrossTalk::XTalk(double x, bool type) const
{
  TF1* XT_cellTF = new TF1("XT_cellTF","((exp(-x/[0])*x*x)/(2 *[0]*[0]*([0] - [1])) - (exp(-(x/[0]))*x*[1])/([0]*pow([0] - [1],2)) + exp(-(x/[0]))*[1]*[1]/pow([0] - [1],3) + (exp(-(x/[1]))*[1]*[1])/pow(-[0] + [1],3) + (1/(2*[2]*[0] *pow(([0] - [1]),3)))*exp(-x* (1/[0] + 1/[1]))* (-2 *exp(x *(1/[0] + 1/[1]))*[0] *pow(([0] - [1]),3) - 2 *exp(x/[0])*[0]*pow([1],3) + exp(x/[1]) *(x*x *pow(([0] - [1]),2) + 2*x*[0]*([0]*[0] - 3*[0]*[1] + 2*[1]*[1]) + 2*[0]*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))) + ((1 - (exp((-x + [2])/[0])*(x - [2])*([0] - 2*[1]))/pow(([0] - [1]),2) - (exp((-x + [2])/[0])*(x - [2])*(x- [2]))/(2*[0]*([0] - [1])) + (exp((-x + [2])/[1])*[1]*[1]*[1])/pow(([0] - [1]),3) - (exp((-x + [2])/[0])*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))/pow(([0] - [1]),3))* 0.5*( 1+sign(1, x -[2]) ) )/[2])*[3]*[4]*[3]*[0]*[0]",0., m_tmax2);
  XT_cellTF->SetParameter(0, m_taud);
  XT_cellTF->SetParameter(1, m_taupa);
  XT_cellTF->SetParameter(2, m_td);
  XT_cellTF->SetParameter(3, m_Rf);
  XT_cellTF->SetParameter(4, m_C1);

  double xt_cell = 0 ;

  if (type){
    xt_cell = XT_cellTF->Derivative(x) ;
  }
  else {
    xt_cell = XT_cellTF->Eval(x) ;
  }

  delete XT_cellTF ;

  return xt_cell ;
}// end of function


double CrossTalk::CellFunction(double x, bool type) const
{
  TF1* CellM = new TF1("CellM","[5]*((exp(-x/[0])*x*x)/(2 *[0]*[0]*([0] - [1])) - (exp(-(x/[0]))*x*[1])/([0]*pow([0] - [1],2)) + exp(-(x/[0]))*[1]*[1]/pow([0] - [1],3) + (exp(-(x/[1]))*[1]*[1])/pow(-[0] + [1],3) + (1/(2*[2]*[0] *pow(([0] - [1]),3)))*exp(-x* (1/[0] + 1/[1]))* (-2 *exp(x *(1/[0] + 1/[1]))*[0] *pow(([0] - [1]),3) - 2 *exp(x/[0])*[0]*pow([1],3) + exp(x/[1]) *(x*x *pow(([0] - [1]),2) + 2*x*[0]*([0]*[0] - 3*[0]*[1] + 2*[1]*[1]) + 2*[0]*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))) + ((1 - (exp((-x + [2])/[0])*(x - [2])*([0] - 2*[1]))/pow(([0] - [1]),2) - (exp((-x + [2])/[0])*(x - [2])*(x- [2]))/(2*[0]*([0] - [1])) + (exp((-x + [2])/[1])*[1]*[1]*[1])/pow(([0] - [1]),3) - (exp((-x + [2])/[0])*[0]*([0]*[0] - 3*[0]*[1] + 3*[1]*[1]))/pow(([0] - [1]),3))* 0.5*( 1+sign(1,x -[2]) ) )/[2])*[3]*[4]*[3]*[0]*[0]",0., m_tmax2);

  CellM->SetParameter(0, m_taud);
  CellM->SetParameter(1, m_taupa);
  CellM->SetParameter(2, m_td);
  CellM->SetParameter(3, m_Rf);
  CellM->SetParameter(4, m_C1);

  double cell = 0 ;

  if (type){
      cell = CellM->Derivative(x) ;        
  }
  else {
      cell = CellM->Eval(x) ;
  }

  delete CellM ;

  return cell ;
}// end of function


StatusCode CrossTalk::fillHistograms( SG::EventContext &/*ctx*/ ) const
{
  return StatusCode::SUCCESS;
}


StatusCode CrossTalk::bookHistograms( SG::EventContext &/*ctx*/ ) const
{
  return StatusCode::SUCCESS;
}

