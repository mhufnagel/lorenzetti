#ifndef CalibrationTool_h
#define CalibrationTool_h

#include "GaugiKernel/StatusCode.h"
#include "GaugiKernel/AlgTool.h"
#include "GaugiKernel/EDM.h"
#include "CaloCell/CaloDetDescriptor.h"
#include "CaloCell/CaloCalibHelper.h"


class CalibrationTool : public Gaugi::AlgTool
{
    public:
        /* Constructor */    
        CalibrationTool(std::string name);
        virtual ~CalibrationTool();

        virtual StatusCode initialize() override;
        virtual StatusCode finalize() override;
        virtual StatusCode execute( SG::EventContext &ctx, Gaugi::EDM * ) const override;


    private:
        
};

#endif