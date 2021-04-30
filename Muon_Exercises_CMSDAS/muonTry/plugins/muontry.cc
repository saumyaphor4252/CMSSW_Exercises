// -*- C++ -*-
//
// Package:    MuonAnalysis/muontry
// Class:      muontry
// 
/**\class muontry muontry.cc MuonAnalysis/muontry/plugins/muontry.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Mon, 06 Jul 2020 08:21:35 GMT
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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class muontry : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit muontry(const edm::ParameterSet&);
      ~muontry();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonsToken;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genMuonsToken;
  
      size_t nPatTotal, nPatGlobal, nPatPt, nPatEta, nPatChargedIso, nPatPlusMuons, nPatMinusMuons,
             nGenTotal, nGenEta, nGenPt, nGenPlusMuons, nGenMinusMuons; 




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
muontry::muontry(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::InputTag muonLabel("slimmedMuons");
   edm::InputTag genMuonLabel ("packedGenParticles");

   muonsToken=consumes<pat::MuonCollection>(muonLabel);
   genMuonsToken=consumes<pat::PackedGenParticleCollection>(genMuonLabel);
   edm::Service<TFileService>fs;
 
   nPatTotal=0; nPatGlobal=0; nPatPt=0; nPatEta=0; nPatChargedIso=0; nPatPlusMuons=0; nPatMinusMuons=0;
   nGenTotal=0; nGenEta=0; nGenPt=0; nGenPlusMuons=0; nGenMinusMuons=0; 

}


muontry::~muontry()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
muontry::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonsToken,muons);

   edm::Handle<pat::PackedGenParticleCollection> genMuons;
   iEvent.getByToken(genMuonsToken,genMuons);

   for(auto mu=muons->cbegin();mu!=muons->cend();++mu)
   {
      nPatTotal++;
//      if(! mu->isGlobalMuon()) continue;
//      nPatGlobal++;
      if(! (mu->pt()>20) ) continue;
      nPatPt++;
      if(fabs(mu->eta())>2.4) continue;
      nPatEta++;
//      if(! (mu->chargedHadronIso()<0.15) ) continue;
//      nPatChargedIso++;
      if( mu->charge()>0 ) 
      nPatPlusMuons++;
      if( mu->charge()<0 ) 
      nPatMinusMuons++;
	  
   }
   
   for(auto mu_G=genMuons->cbegin();mu_G!=genMuons->cend();++mu_G)
   {
      const pat::PackedGenParticle& mcMuon=(*mu_G);
      if(! (fabs(mcMuon.pdgId())==13)) continue;
      nGenTotal++;
      //if(not mcMuon.isGlobalMuon()) continue;
      //nGenGlobal++;
      if(! (mcMuon.pt()>20) ) continue;
      nGenPt++;
      if(fabs(mcMuon.eta())>2.4) continue;
      nGenEta++;
      if( mcMuon.charge()>0 ) 
      nGenPlusMuons++;
      if( mcMuon.charge()<0 ) 
      nGenMinusMuons++;   
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
muontry::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
muontry::endJob() 
{
   std::cout<<"nPatTotal: "<<nPatTotal<<std::endl;
   std::cout<<"nPatGlobal: "<<nPatGlobal<<std::endl;
   std::cout<<"nPatPt: "<<nPatPt<<std::endl;
   std::cout<<"nPatEta: "<<nPatEta<<std::endl;
   std::cout<<"nPatChargedIso: "<<nPatChargedIso<<std::endl;
   std::cout<<"nPatPlusMuons: "<<nPatPlusMuons<<std::endl;
   std::cout<<"nPatMinusMuons: "<<nPatMinusMuons<<std::endl;
   std::cout<<"nGenTotal: "<<nGenTotal<<std::endl;
 //  std::cout<<"nGenGlobal: "<<nGenGlobal<<std::endl;
   std::cout<<"nGenPt: "<<nGenPt<<std::endl;
   std::cout<<"nGenEta: "<<nGenEta<<std::endl;
   std::cout<<"nGenPlusMuons: "<<nGenPlusMuons<<std::endl;
   std::cout<<"nGenMinusMuons: "<<nGenMinusMuons<<std::endl; 


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
muontry::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(muontry);
