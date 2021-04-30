// -*- C++ -*-
//
// Package:    VertexAnalysis/VertexAnalyzerAOD
// Class:      VertexAnalyzerAOD
// 
/**\class VertexAnalyzerAOD VertexAnalyzerAOD.cc VertexAnalysis/VertexAnalyzerAOD/plugins/VertexAnalyzerAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Fri, 29 May 2020 12:00:52 GMT
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <TTree.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class VertexAnalyzerAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit VertexAnalyzerAOD(const edm::ParameterSet&);
      ~VertexAnalyzerAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;

      TTree* treeEvent;
    
      int number_of_vertices;
      int number_of_events;
 //     std::vector<int> v_nVertices;

      std::vector<double> x_coordinate; 
      std::vector<double> y_coordinate;
      std::vector<double> z_coordinate;
      std::vector<double> tracks_from_each_vertex;
      std::vector<int> vertices_per_event;
      std::vector<int> tracks_per_event;
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
VertexAnalyzerAOD::VertexAnalyzerAOD(const edm::ParameterSet& iConfig)
: verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");

   number_of_vertices=0;
   number_of_events=0;
   vertices_per_event.clear();
   tracks_per_event.clear();
}


VertexAnalyzerAOD::~VertexAnalyzerAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VertexAnalyzerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   x_coordinate.clear();
   y_coordinate.clear();
   z_coordinate.clear();
   tracks_from_each_vertex.clear();   

   //    Defining tokens
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(verticesToken_,vertices);
  // if ( !vertices.isValid() ) return;
   vertices_per_event.push_back(vertices->size());

   
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByToken(tracksToken_,tracks);
   tracks_per_event.push_back(tracks->size());

   //     LOOP over Vertices 
   //v_nVertices.push_back(vertices->size());
   for(reco::VertexCollection::const_iterator itVertex=vertices->begin(); itVertex!=vertices->end(); ++itVertex)
   {
       x_coordinate.push_back(itVertex->x());
       y_coordinate.push_back(itVertex->y());
       z_coordinate.push_back(itVertex->z());
       if(number_of_vertices<2){std::cout<<itVertex->x()<<std::endl;}     
       tracks_from_each_vertex.push_back(itVertex->tracksSize());
       number_of_vertices++;
   }
   
   treeEvent->Fill();
   number_of_events++;  
   
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
VertexAnalyzerAOD::beginJob()
{
   treeEvent->Branch("X_Coordinate",&x_coordinate);
   treeEvent->Branch("Y_Cooridnate", &y_coordinate);
   treeEvent->Branch("Z_Cooridnate", &z_coordinate);
   treeEvent->Branch("Tracks_per_vertex",&tracks_from_each_vertex);
   treeEvent->Branch("Tracks",&tracks_per_event);
   treeEvent->Branch("Vertices",&vertices_per_event);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexAnalyzerAOD::endJob() 
{ 
   std::cout<<"Total events: "<<number_of_events<<std::endl;
   std::cout<<"Tota; vertices in all events: "<<number_of_vertices<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexAnalyzerAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexAnalyzerAOD);
