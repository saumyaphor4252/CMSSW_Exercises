// -*- C++ -*-
//
// Package:    TrackAnalysisAOD/TrackAnalyzer
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc TrackAnalysisAOD/TrackAnalyzer/plugins/TrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Mon, 25 May 2020 15:43:42 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
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

class TrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrackAnalyzer(const edm::ParameterSet&);
      ~TrackAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_; 
      TTree* treeEvent;

      int indexEvent_;
      int numTotal;
      int numLoose;
      int numTight; 
      int numHighPurity;

      std::vector<double> pt;
      std::vector<double> eta;
      std::vector<double> phi;
      std::vector<double> Loose_nvh;
      std::vector<double> Loose_numPixelHits;
      std::vector<double> Loose_numValidHits;
      std::vector<double> Loose_numTkLayers;
      std::vector<double> Tight_nvh;
      std::vector<double> Tight_numPixelHits;
      std::vector<double> Tight_numValidHits;
      std::vector<double> Tight_numTkLayers;
      std::vector<double> HighPurity_nvh;
      std::vector<double> HighPurity_numPixelHits;
      std::vector<double> HighPurity_numValidHits;
      std::vector<double> HighPurity_numTkLayers;
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
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)
: tracksToken_(consumes<edm::View<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("tracks", edm::InputTag("generalTracks")) ))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");
   indexEvent_=0;
   numTotal = 0;
   numLoose = 0;
   numTight = 0;
   numHighPurity = 0;

}


TrackAnalyzer::~TrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   pt.clear();
   eta.clear();
   phi.clear();
   Loose_nvh.clear();
   Loose_numPixelHits.clear(); 
   Loose_numValidHits.clear();
   Loose_numTkLayers.clear();
   Tight_nvh.clear();
   Tight_numPixelHits.clear();
   Tight_numValidHits.clear();
   Tight_numTkLayers.clear();
   HighPurity_nvh.clear();
   HighPurity_numPixelHits.clear();
   HighPurity_numValidHits.clear();
   HighPurity_numTkLayers.clear();
 
   edm::Handle<edm::View<reco::Track> > trackHandle;
   iEvent.getByToken(tracksToken_, trackHandle);
   if ( !trackHandle.isValid() ) return;

   const edm::View<reco::Track>& tracks = *trackHandle;

   size_t iTrack = 0;

   for ( auto track : tracks )
   {
     ++numTotal;
     pt.push_back(track.pt());
     eta.push_back(track.eta());
     phi.push_back(track.phi());

     if (track.quality(track.qualityByName("loose"))     )
         {   ++numLoose;
             Loose_nvh.push_back(track.numberOfValidHits());
             Loose_numPixelHits.push_back(track.hitPattern().numberOfValidPixelHits());
             Loose_numValidHits.push_back(track.hitPattern().numberOfValidHits());
             Loose_numTkLayers.push_back(track.hitPattern().trackerLayersWithMeasurement());
         }
     if (track.quality(track.qualityByName("tight"))     )
         {   ++numTight;
             Tight_nvh.push_back(track.numberOfValidHits());
             Tight_numPixelHits.push_back(track.hitPattern().numberOfValidPixelHits());
             Tight_numValidHits.push_back(track.hitPattern().numberOfValidHits());
             Tight_numTkLayers.push_back(track.hitPattern().trackerLayersWithMeasurement());
         }
     if (track.quality(track.qualityByName("highPurity")))
         {   ++numHighPurity;
             HighPurity_nvh.push_back(track.numberOfValidHits()); 
             HighPurity_numPixelHits.push_back(track.hitPattern().numberOfValidPixelHits());
             HighPurity_numValidHits.push_back(track.hitPattern().numberOfValidHits());
             HighPurity_numTkLayers.push_back(track.hitPattern().trackerLayersWithMeasurement());
         }
     iTrack++;
   }

   ++indexEvent_;
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
TrackAnalyzer::beginJob()
{
   treeEvent->Branch("pt", &pt);
   treeEvent->Branch("eta", &eta);
   treeEvent->Branch("phi", &phi);
   treeEvent->Branch("Loose_nvh",&Loose_nvh);
   treeEvent->Branch("Loose_numPixelHits",&Loose_numPixelHits);
   treeEvent->Branch("Loose_numValidHits",&Loose_numValidHits);
   treeEvent->Branch("Loose_numTkLayers", &Loose_numTkLayers);
   treeEvent->Branch("Tight_nvh", &Tight_nvh);
   treeEvent->Branch("Tight_numPixelHits", &Tight_numPixelHits);
   treeEvent->Branch("Tight_numValidHits", &Tight_numValidHits);
   treeEvent->Branch("Tight_numTkLayers", &Tight_numTkLayers);
   treeEvent->Branch("HighPurity_nvh", &HighPurity_nvh);
   treeEvent->Branch("HighPurity_numPixelHits", &HighPurity_numPixelHits);
   treeEvent->Branch("HighPurity_numValidHits", &HighPurity_numValidHits);
   treeEvent->Branch("HighPurity_numTkLayers", &HighPurity_numTkLayers);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalyzer::endJob() 
{
   std::cout << "Events " << indexEvent_<< std::endl
             << " numTotal: " << numTotal << std::endl
             << " numLoose: " << numLoose << std::endl
             << " numTight: " << numTight  << std::endl
             << " numHighPurity: " << numHighPurity << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
