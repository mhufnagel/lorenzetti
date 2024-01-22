#include "CalibrationTool.h"

using namespace Gaugi;


CalibrationTool::CalibrationTool( std::string name ) :
    IMsgService(name),
    AlgTool()
    {
        declareProperty( "OutputLevel"      , m_outputLevel=1   );
    }

//!=====================================================================

CalibrationTool::~CalibrationTool()
{;}

//!=====================================================================

StatusCode CalibrationTool::initialize(){
    setMsgLevel( (MSG::Level)m_outputLevel );
    MSG_DEBUG( " Initializing Calibration tool.");
    return StatusCode::SUCCESS;
}

//!=====================================================================

StatusCode CalibrationTool::finalize(){
    MSG_DEBUG( " Finalizing Calibration tool.");
    return StatusCode::SUCCESS;
}
//!=====================================================================

StatusCode CalibrationTool::execute(SG::EventContext &/*ctx*/, Gaugi::EDM *edm ) const{

    auto *descriptor    = static_cast<xAOD::CaloDetDescriptor*>(edm);
    auto *calibHelper   = descriptor->calibHelper();

    // Time Of Flight Calibration
    if (descriptor->tof() > 0.0){
        float oldTof = descriptor->tof();

        descriptor->setTofDelay( calibHelper->tofIP() ); //calibrate the TOF for impact point distance

        if (m_outputLevel == 2) MSG_INFO("Cell "<< descriptor->hash() <<", sampling "<< descriptor->sampling() << ": tof set from "<< oldTof << " to "<< descriptor->tof()<< ". RMin=" << calibHelper->RMin() << ", RMax=" << calibHelper->RMax() << ", yn="<< calibHelper->yn()<< ", zn="<< calibHelper->zn()<< ", rn="<< calibHelper->rn()<<", eta/theta="<< calibHelper->eta()<<"/"<<calibHelper->thetaDeg()<<", tof_ip="<<calibHelper->tofIP());
    }


    return StatusCode::SUCCESS;
}