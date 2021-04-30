// -*- C++ -*-
//
// Package:    TrackAnalysis/patTrackAnalyzer
// Class:      patTrackAnalyzer
// 
/**\class patTrackAnalyzer patTrackAnalyzer.cc TrackAnalysis/patTrackAnalyzer/plugins/patTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Mon, 25 May 2020 12:00:25 GMT
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


#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class patTrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit patTrackAnalyzer(const edm::ParameterSet&);
      ~patTrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfsToken_;
      TTree* treeEvent;
      
      std::vector<double> pt;
      std::vector<double> eta;
      std::vector<double> phi;
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
patTrackAnalyzer::patTrackAnalyzer(const edm::ParameterSet& iConfig)
:pfsToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfs")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");
}


patTrackAnalyzer::~patTrackAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
patTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfsToken_,pfs);

   std::cout<<"number of tracks: "<<pfs->size();

   for(pat::PackedCandidateCollection::const_iterator itTrack=pfs->begin(); itTrack!=pfs->end(); ++itTrack)
   {
     pt.push_back(itTrack->pt());
     eta.push_back(itTrack->eta());
     phi.push_back(itTrack->phi());
     
   }

   treeEvent->Fill();


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
patTrackAnalyzer::beginJob()
{
   treeEvent->Branch("pt", &pt);
   treeEvent->Branch("eta", &eta);
   treeEvent->Branch("phi", &phi);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
patTrackAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
patTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(patTrackAnalyzer);
