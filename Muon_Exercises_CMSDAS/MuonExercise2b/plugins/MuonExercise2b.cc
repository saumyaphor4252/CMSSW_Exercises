// -*- C++ -*-
//
// Package:    MuonAnalysis/MuonExercise2
// Class:      MuonExercise2
// 
/**\class MuonExercise2 MuonExercise2.cc MuonAnalysis/MuonExercise2/plugins/MuonExercise2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Thu, 11 Jun 2020 14:59:13 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TProfile.h>
//#include "TRandom3.h"
//#include "RoccoR.cc"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuonExercise2b : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonExercise2b(const edm::ParameterSet&);
      ~MuonExercise2b();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

 //     RoccoR rc("//eos/uscms/store/user/cmsdas/2019/short_exercises/Muons/Exercises/MuonExercise2/plugin/rcdata.2016.v3"); 
 //     TRandom3 rnd(1234);
      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonsToken;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genMuonsToken;

      size_t nPatTotal,nPatPtCut,nPatEtaCut,nPatVectorSizeCut, nPatOppositeChargeCut,nPatMassCut,
             nGenTotal,nGenPtCut,nGenEtaCut,nGenVectorSizeCut,nGenOppositeChargeCut,nGenMassCut;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonExercise2b::MuonExercise2b(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   edm::InputTag muonLabel("slimmedMuons");
   edm::InputTag genMuonLabel ("packedGenParticles");

   muonsToken=consumes<pat::MuonCollection>(muonLabel); 
   genMuonsToken=consumes<pat::PackedGenParticleCollection>(genMuonLabel);

   nPatTotal=0; nPatPtCut=0; nPatEtaCut=0; nPatVectorSizeCut=0; nPatOppositeChargeCut=0; nPatMassCut=0; 
   nGenTotal=0; nGenPtCut=0; nGenEtaCut=0; nGenVectorSizeCut=0; nGenOppositeChargeCut=0; nGenMassCut=0;
}
MuonExercise2b::~MuonExercise2b()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise2b::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonsToken,muons);

   edm::Handle<pat::PackedGenParticleCollection> genMuons;
   iEvent.getByToken(genMuonsToken,genMuons);
   
   //Loop over PAT Muons
   std::vector<int> vpat,vgen;
   for(auto mu=muons->cbegin();mu!=muons->cend();++mu)
   {  
      nPatTotal++;
      if(! (mu->pt()>20) ) continue;
      nPatPtCut++;
      if(fabs(mu->eta())>2.4) continue;
      nPatEtaCut++;
      vpat.push_back( std::distance(muons->cbegin(),mu) );
   }
   if(vpat.size()==2) 
   {
      nPatVectorSizeCut++;
      if( muons->at(vpat[0]).charge()*muons->at(vpat[1]).charge()==-1 ) 
      {
         nPatOppositeChargeCut++;
         int idx_mum=-1,idx_mup=-1;
         if( muons->at(vpat[0]).charge()>0 ) { idx_mup=vpat[0]; idx_mum=vpat[1]; }
         else if( muons->at(vpat[0]).charge()<0 ) { idx_mup=vpat[1]; idx_mum=vpat[0]; }
         double MRec=( muons->at(idx_mup).p4()+ muons->at(idx_mum).p4() ).M();
         if( MRec>70 && MRec<110 ) nPatMassCut++;
       }
   }

   //Loop over Gen Muons
   for(auto mu_G=genMuons->cbegin();mu_G!=genMuons->cend();++mu_G)
   {
       const pat::PackedGenParticle& mcMuon=(*mu_G);
       if(! (fabs(mcMuon.pdgId())==13)) continue;
       nGenTotal++;
       if(! (mcMuon.pt()>20) ) continue;
       nGenPtCut++;
       if(fabs(mcMuon.eta())>2.4) continue;
       nGenEtaCut++;
       vgen.push_back( std::distance(genMuons->cbegin(),mu_G) );
   }
   if(vgen.size()==2) 
   {      
      nGenVectorSizeCut++;
      if( genMuons->at(vgen[0]).charge()*genMuons->at(vgen[1]).charge()==-1 ) 
      {
         nGenOppositeChargeCut++;
         int idx_mum_g=-1,idx_mup_g=-1;
         if( genMuons->at(vgen[0]).charge()>0 ) { idx_mup_g=vgen[0]; idx_mum_g=vgen[1]; }
         else if( genMuons->at(vgen[0]).charge()<0 ) { idx_mup_g=vgen[1]; idx_mum_g=vgen[0]; }
         double MGen=( genMuons->at(idx_mup_g).p4()+ genMuons->at(idx_mum_g).p4() ).M();
         if( MGen>70 && MGen<110 ) nGenMassCut++;
      }
   }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonExercise2b::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise2b::endJob() 
{  
   std::cout<<"Total PAT Muons: "<<nPatTotal<<std::endl;
   std::cout<<"Pat Muons After pt Cut: "<<nPatPtCut<<std::endl;
   std::cout<<"Pat Muons After eta Cut: "<<nPatEtaCut<<std::endl;
   std::cout<<"Pat Muons After vector size Cut: "<<nPatVectorSizeCut<<std::endl;
   std::cout<<"Pat Muons After opposite charge Cut: "<<nPatOppositeChargeCut<<std::endl;
   std::cout<<"Pat Muons After Mass range Cut: "<<nPatMassCut<<std::endl;

   std::cout<<"Total Gen Muons: "<<nGenTotal<<std::endl;
   std::cout<<"Gen Muons After pt Cut: "<<nGenPtCut<<std::endl;
   std::cout<<"Gen Muons After eta Cut: "<<nGenEtaCut<<std::endl;
   std::cout<<"Gen Muons After vector size Cut: "<<nGenVectorSizeCut<<std::endl;
   std::cout<<"Gen Muons After opposite charge Cut: "<<nGenOppositeChargeCut<<std::endl;
   std::cout<<"Gem Muons After Mass range Cut: "<<nGenMassCut<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonExercise2b::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise2b);
