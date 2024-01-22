
#include "CaloCell/CaloCalibHelper.h"
// #include "CaloCell/CaloDetDescriptor.h"
// #include "CaloCell/CaloCell.h"
// #include "CaloCell/enumeration.h"
// #include "G4Kernel/constants.h"
#include <math.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


using namespace xAOD;

CaloCalibHelper::CaloCalibHelper(
                                    float Rmin, 
                                    float Rmax,
                                    float eta): 
    EDM(),
    m_Rmin(Rmin),
    m_Rmax(Rmax),
    m_eta(eta)
    // m_phi(phi),
    // m_deta(deta),
    // m_dphi(dphi)
{
    // Initialize geometric variables
    // Theta
    float arg       = -exp(abs(eta));
    float theta_rad = 2*atan(arg) +  M_PI;
    m_theta_deg     = theta_rad * 180/M_PI;

    // Projections
    m_yn = (m_Rmin*mm + m_Rmax*mm)/2;  
    m_rn = m_yn/sin(theta_rad);
    m_zn = m_rn*cos(theta_rad);

    // TOF Delay from IP
    m_tof_ip = m_rn/c_light;
}