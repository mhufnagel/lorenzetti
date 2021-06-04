#ifndef CaloHit_h
#define CaloHit_h

#include "CaloCell/enumeration.h"
#include "GaugiKernel/EDM.h"
#include "GaugiKernel/macros.h"
#include "G4Step.hh"
#include "globals.hh"



namespace xAOD{

  class CaloHit: public Gaugi::EDM
  {  
    public:

      CaloHit():EDM(){};

      /** Contructor **/
      CaloHit( 
               // eta/phi position at the detector, deta/dphi size into the detector, radius
               float eta, 
               float phi, 
               float deta, 
               float dphi,
               // Hash
               unsigned int hash,
               // cell identification
               CaloSampling sampling, 
               Detector detector,
               // bunch crossing information
               float bc_duration, 
               int bcid_start, 
               int bcid_end 
              );

      /** Destructor **/
      ~CaloHit()=default;
      
      /*! Fill the deposit energy into the cell */
      void Fill( const G4Step * );
      /** Zeroize the pulse/sample vectors **/
      void clear();


      /*
       * Cell identification
       */

      /*! Cell eta center */
      PRIMITIVE_SETTER_AND_GETTER( float, m_eta, setEta, eta );
      /*! Cell phi center */
      PRIMITIVE_SETTER_AND_GETTER( float, m_phi, setPhi, phi );
      /*! Cell delta eta */
      PRIMITIVE_SETTER_AND_GETTER( float, m_deta, setDeltaEta , deltaEta);
      /*! Cell delta phi */
      PRIMITIVE_SETTER_AND_GETTER( float, m_dphi, setDeltaPhi, deltaPhi );
      /*! Cell hash */
      PRIMITIVE_SETTER_AND_GETTER( unsigned int, m_hash, setHash, hash );
      /*! Cell sampling id */
      PRIMITIVE_SETTER_AND_GETTER( CaloSampling  , m_sampling , setSampling   , sampling  );
      /*! Cell layer id */
      PRIMITIVE_SETTER_AND_GETTER( Detector  , m_detector  , setDetector    , detector   );
     



      /*! Energy deposity from simulated hits **/
      float edep( int bc_id=0 ) const{
        if (m_edep.count(bc_id)){
          return m_edep.at(bc_id);
        }else{// return zero case bc not exist
          return 0;
        }
      }

      void edep( int bc_id, float e ){
        m_edep[bc_id] += e;
      }


      /*
       * Bunch crossing information
       */

      /*! Bunch crossing id start */
      PRIMITIVE_SETTER_AND_GETTER( int, m_bcid_start  , set_bcid_start  , bcid_start    );
      /*! Bunch crossing id end */
      PRIMITIVE_SETTER_AND_GETTER( int, m_bcid_end    , set_bcid_end    , bcid_end      );
      /* Time space (in ns) between two bunch crossings */
      PRIMITIVE_SETTER_AND_GETTER( float, m_bc_duration , set_bc_duration , bc_duration );
      /*! Time (in ns) for each bunch crossing */
      PRIMITIVE_SETTER_AND_GETTER( std::vector<float> , m_time , setTime , time   );



    private:
 
      int findIndex( float value) const ;
       

      /*! id sample */
      CaloSampling m_sampling;
      /*! id layer */
      Detector m_detector;

      
      /*! eta center */
      float m_eta;
      /*! phi center */
      float m_phi;
      /*! delta eta */
      float m_deta;
      /*! delta phi */
      float m_dphi;


      /*! bunch crossing start id */
      int m_bcid_start;
      /*! bunch crossing end id */
      int m_bcid_end;
      /*! bunch crossing space in ns between two bunchs */
      float m_bc_duration;


      /*! time (in ns) for each bunch between bcid_start and bcid_end */
      std::vector<float> m_time;
      /*! energy deposit between bcid_start and bcid_end */
      std::map< int, float> m_edep;
     
      /*! Access information unique ID number */
      unsigned int m_hash;

  };

}
#endif