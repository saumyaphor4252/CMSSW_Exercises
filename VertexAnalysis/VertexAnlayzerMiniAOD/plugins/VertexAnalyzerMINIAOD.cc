// -*- C++ -*-
//
// Package:    VertexAnalysis/VertexAnalyzerMINIAOD
// Class:      VertexAnalyzerMINIAOD
// 
/**\class VertexAnalyzerMINIAOD VertexAnalyzerMINIAOD.cc VertexAnalysis/VertexAnalyzerMINIAOD/plugins/VertexAnalyzerMINIAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Sat, 30 May 2020 10:41:42 GMT
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
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
 #include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include <TTree.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class VertexAnalyzerMINIAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit VertexAnalyzerMINIAOD(const edm::ParameterSet&);
      ~VertexAnalyzerMINIAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
 /*     edm::EDGetTokenT<reco::VertexCompositePtrCandidate> KshortsToken_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> LambdasToken_;
   */       edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate> > KshortsToken_;
  
    TTree* treeEvent;

      int Kshort_number_of_vertices;
      int number_of_events;
      int Lambda_number_of_vertices;
      std::vector<double> lambda_mass;
      std::vector<double> kshort_mass;

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
VertexAnalyzerMINIAOD::VertexAnalyzerMINIAOD(const edm::ParameterSet& iConfig)
: KshortsToken_(consumes<edm::View<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("KshortVertices")) )
//: KshortsToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortVertices"))),
// LambdasToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaVertices")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");

 //  Lambda_number_of_vertices=0;
   Kshort_number_of_vertices=0;
   number_of_events=0;

}


VertexAnalyzerMINIAOD::~VertexAnalyzerMINIAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VertexAnalyzerMINIAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   kshort_mass.clear();
 //  lambda_mass.clear();

 //  edm::Handle<reco::VertexCompositePtrCandidateCollection> kshort_vertices;
    edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > kshort_vertices;
   iEvent.getByToken(KshortsToken_,kshort_vertices);
   if ( !kshort_vertices.isValid() ) return;
  
  // edm::Handle<reco::VertexCompositePtrCandidateCollection> lambda_vertices;
 //  iEvent.getByToken(LambdasToken_,lambda_vertices);
//   if ( !lambda_vertices.isValid() ) return;
 
 // for(reco::VertexCompositePtrCandidateCollection::const_iterator itKshort=kshort_vertices->begin(); itKshort!=kshort_vertices->end(); ++itKshort)
      const edm::View<reco::VertexCompositePtrCandidate>& KshortVertices = *kshort_vertices;
   size_t iTrack = 0;
   for ( auto Kshort : KshortVertices )
    {
      kshort_mass.push_back(Kshort.mass());
      Kshort_number_of_vertices++;
   }
/*
  for(reco::VertexCompositePtrCandidateCollection::const_iterator itLambda=lambda_vertices->begin(); itLambda!=lambda_vertices->end(); ++itLambda)
   {
       lambda_mass.push_back(itLambda->mass());
       Lambda_number_of_vertices++;
   }
*/
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
VertexAnalyzerMINIAOD::beginJob()
{
   treeEvent->Branch("Kshort",&kshort_mass);
 //  treeEvent->Branch("Lambda", &lambda_mass);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexAnalyzerMINIAOD::endJob() 
{
   std::cout<<"Total events: "<<number_of_events<<std::endl;
   std::cout<<"Total Kshort vertices in all events: "<<Kshort_number_of_vertices<<std::endl;
 //  std::cout<<"Total Lambda vertices in all events: "<<Lambda_number_of_vertices<<std::endl;
}

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexAnalyzerMINIAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexAnalyzerMINIAOD);
