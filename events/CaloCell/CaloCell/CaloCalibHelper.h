#ifndef CaloCalibHelper_h
#define CaloCalibHelper_h

/** simulator includes **/
#include "CaloCell/enumeration.h"
// #include "CaloCell/CaloDetDescriptor.h"
// #include "CaloCell/CaloCell.h"
#include "GaugiKernel/EDM.h"
#include "GaugiKernel/macros.h"


namespace xAOD{

    class CaloCalibHelper: public Gaugi::EDM
    {
        public:

            /** Constructor **/
            CaloCalibHelper()=default;

            /** Constructor **/
            CaloCalibHelper(float Rmin, float Rmax, float eta);

            /** Destructor **/
            ~CaloCalibHelper()=default;

            /* setters and getters */
            PRIMITIVE_SETTER_AND_GETTER(float, m_Rmin, setRMin,  RMin);
            PRIMITIVE_SETTER_AND_GETTER(float, m_Rmax, setRMax,  RMax);
            PRIMITIVE_SETTER_AND_GETTER(float, m_eta,  setEta,   eta);
            // PRIMITIVE_SETTER_AND_GETTER(float, m_phi,  setPhi,   phi);
            // PRIMITIVE_SETTER_AND_GETTER(float, m_deta, setDeta,  deta);
            // PRIMITIVE_SETTER_AND_GETTER(float, m_dphi, setDphi,  dphi);

            PRIMITIVE_SETTER_AND_GETTER( float, m_yn,        setyn,         yn );
            PRIMITIVE_SETTER_AND_GETTER( float, m_zn,        setzn,         zn );
            PRIMITIVE_SETTER_AND_GETTER( float, m_rn,        setrn,         rn );
            PRIMITIVE_SETTER_AND_GETTER( float, m_tof_ip,    setTofIP,      tofIP );
            PRIMITIVE_SETTER_AND_GETTER( float, m_theta_deg, setThetaDeg,   thetaDeg );
            // PRIMITIVE_SETTER_AND_GETTER( float, m_theta_rad, setThetaRad,   thetaRad );

        private:

            // Geometry from detector

            /* ! sampling 'y' coordinate start. Perpendicular to beam-axis (theta=90°), in mm */
            float m_Rmin;
            /* ! sampling 'y' coordinate end. Perpendicular to beam-axis (theta=90°) , in mm*/
            float m_Rmax;
            /*! cell eta center */
            float m_eta;
            // /*! cell phi center */
            // float m_phi;
            // /*! cell delta eta */
            // float m_deta;
            // /*! cell delta phi */
            // float m_dphi;


            // New Geometry variables calculated from detector

            /* ! sampling 'y' coordinate centroid. Projection of theta angle. */
            float m_yn;
            /* ! cell 'z' coordinate centroid. Projection of m_rn to z-axis. */
            float m_zn;
            /* ! cell centroid distance to impact point (0,0,0) in mm*/
            float m_rn;
            /* ! time of flight from impact point (0,0,0) to cell centroid at speed of light */
            float m_tof_ip;
            /* ! converted eta into theta angle in degrees */
            float m_theta_deg;
            /* ! converted eta into theta angle in rad */
            // float m_theta_rad;

    };
}
#endif
