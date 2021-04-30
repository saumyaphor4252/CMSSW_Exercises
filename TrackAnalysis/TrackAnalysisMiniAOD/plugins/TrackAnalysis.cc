// -*- C++ -*-
//
// Package:    TrackAnalysisMiniAOD/TrackAnalysis
// Class:      TrackAnalysis
// 
/**\class TrackAnalysis TrackAnalysis.cc TrackAnalysisMiniAOD/TrackAnalysis/plugins/TrackAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Tue, 26 May 2020 13:08:38 GMT
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

class TrackAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrackAnalysis(const edm::ParameterSet&);
      ~TrackAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > tracksToken_;
  //    edm::EDGetTokenT<edm::View<pat::PackedCandidate> > losttracksToken_;
 
      int indexEvent_;
      edm::Service<TFileService> fs;
      TTree* treeEvent;
      std::vector<double> track_pt;
      std::vector<double> track_eta;
      std::vector<double> track_phi;
      std::vector<double> track_normchi2;
      std::vector<double> track_nPixelHits;
      std::vector<double> track_nValidHits;
      std::vector<double> track_nTkLayers;
      std::vector<double> HighPurityTrack_pt;
      std::vector<double> HighPurityTrack_eta;
      std::vector<double> HighPurityTrack_phi;
      std::vector<double> HighPurityTrack_normchi2;
      std::vector<double> HighPurityTrack_nPixelHits;
      std::vector<double> HighPurityTrack_nValidHits;
      std::vector<double> HighPurityTrack_nTkLayers;
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
TrackAnalysis::TrackAnalysis(const edm::ParameterSet& iConfig)
: tracksToken_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("tracks")) )
//, losttracksToken_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("losttracks")) )
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");
   indexEvent_=0;
}


TrackAnalysis::~TrackAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   track_pt.clear();
   track_eta.clear();
   track_phi.clear();
   track_normchi2.clear();
   track_nPixelHits.clear();
   track_nValidHits.clear();
   track_nTkLayers.clear();
   HighPurityTrack_pt.clear();
   HighPurityTrack_eta.clear();
   HighPurityTrack_phi.clear();
   HighPurityTrack_normchi2.clear();
   HighPurityTrack_nPixelHits.clear();
   HighPurityTrack_nValidHits.clear();
   HighPurityTrack_nTkLayers.clear();

   edm::Handle<edm::View<pat::PackedCandidate> > trackHandle;
   iEvent.getByToken(tracksToken_, trackHandle);
   if ( !trackHandle.isValid() ) return;

   const edm::View<pat::PackedCandidate>& tracks = *trackHandle;

//   edm::Handle<edm::View<pat::PackedCandidate> > losttrackHandle;
//   iEvent.getByToken(losttracksToken_, losttrackHandle);
//   if ( !losttrackHandle.isValid() ) return;

   size_t iTrack = 0,iHPTrack = 0;
 
//   const edm::View<pat::PackedCandidate>& losttracks = *losttrackHandle;

   for ( auto track : tracks )
   {
        if (!track.hasTrackDetails() || track.charge() == 0 )  
            continue;

        track_pt.push_back(track.pt());
        track_eta.push_back(track.eta());
        track_phi.push_back(track.phi());
       
        track_normchi2.push_back(track.pseudoTrack().normalizedChi2());
        track_nPixelHits.push_back(track.numberOfPixelHits());
        track_nValidHits.push_back(track.pseudoTrack().hitPattern().numberOfValidHits());
        track_nTkLayers.push_back(track.pseudoTrack().hitPattern().trackerLayersWithMeasurement());

        if (track.trackHighPurity())
        {
           HighPurityTrack_pt.push_back(track.pt());
           HighPurityTrack_eta.push_back(track.eta());
           HighPurityTrack_phi.push_back(track.phi());

           HighPurityTrack_normchi2.push_back(track.pseudoTrack().normalizedChi2());
           HighPurityTrack_nPixelHits.push_back(track.numberOfPixelHits());
           HighPurityTrack_nValidHits.push_back(track.pseudoTrack().hitPattern().numberOfValidHits());
           HighPurityTrack_nTkLayers.push_back(track.pseudoTrack().hitPattern().trackerLayersWithMeasurement());

           iHPTrack++;
        }

        iTrack++;
   }

   ++indexEvent_;
   
   std::cout << "Event " << indexEvent_<< std::endl;
   if(indexEvent_ < 5){std::cout<<"Ratio: "<<iHPTrack/iTrack<<std::endl;}
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
TrackAnalysis::beginJob()
{
   treeEvent->Branch("Track_pt", &track_pt);
   treeEvent->Branch("Track_eta", &track_eta);
   treeEvent->Branch("Track_phi", &track_phi);
   treeEvent->Branch("Track_NormalizedChi2", &track_normchi2);
   treeEvent->Branch("Track_nPixelHits", &track_nPixelHits);
   treeEvent->Branch("Track_nValidHits", &track_nValidHits);
   treeEvent->Branch("Track_nTkLayers", &track_nTkLayers);
 
   treeEvent->Branch("HighPurityTrack_pt", &HighPurityTrack_pt);
   treeEvent->Branch("HighPurityTrack_eta", &HighPurityTrack_eta);
   treeEvent->Branch("HighPurityTrack_phi", &HighPurityTrack_phi);
   treeEvent->Branch("HighPurityTrack_NormalizedChi2", &HighPurityTrack_normchi2);
   treeEvent->Branch("HighPurityTrack_nPixelHits", &HighPurityTrack_nPixelHits);
   treeEvent->Branch("HighPurityTrack_nValidHits", &HighPurityTrack_nValidHits);
   treeEvent->Branch("HighPurityTrack_nTkLayers", &HighPurityTrack_nTkLayers);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalysis::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalysis);
