// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Sun, 10 May 2020 09:00:49 GMT
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

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/ValueMap.h"
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

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Electron> > electronsToken_;
      edm::Service<TFileService> fs;
      TTree* treeEvent;

      // all variables for the output tree
      std::vector<double> charge;
      std::vector<double> energy;
      std::vector<double> et;
      std::vector<double> et2;
      std::vector<double> eta;
      std::vector<double> mt;
      std::vector<double> mtSqr;
      std::vector<double> theta;
      std::vector<double> phi;
      std::vector<double> p;
      std::vector<double> pt;
      std::vector<double> px;
      std::vector<double> py;
      std::vector<double> pz;
      
      // VID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapLooseToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMediumToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
     
      // All Electron filters and variables
      std::vector<bool> passEleIdLoose_;
      std::vector<bool> passEleIdMedium_;
      std::vector<bool> passEleIdTight_;
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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
:electronsToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
 eleIdMapLooseToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapLoose"))),
  eleIdMapMediumToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapMedium"))),
  eleIdMapTightToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapTight")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   charge.clear();
   energy.clear();
   et.clear();
   et2.clear();
   eta.clear();
   mt.clear();
   mtSqr.clear();
   theta.clear();
   phi.clear();
   p.clear();
   pt.clear();
   px.clear();
   py.clear();
   pz.clear();
   passEleIdLoose_.clear();
   passEleIdMedium_.clear();
   passEleIdTight_.clear();
  
    using namespace edm;

   //Electrons Token
   edm::Handle<edm::View<pat::Electron> > electronHandle;
   iEvent.getByToken(electronsToken_,electronHandle);

   // Get the electron ID data from the event stream.
   // Note: this implies that the VID ID modules have been run upstream.
   edm::Handle<edm::ValueMap<bool> > loose_ele_id_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_ele_id_decisions;   
   edm::Handle<edm::ValueMap<bool> > tight_ele_id_decisions; 
   iEvent.getByToken(eleIdMapLooseToken_ ,loose_ele_id_decisions);
   iEvent.getByToken(eleIdMapMediumToken_ ,medium_ele_id_decisions);
   iEvent.getByToken(eleIdMapTightToken_ ,tight_ele_id_decisions);

   //const edm::View<pat::Electron>& electrons = *electronHandle;
   //size_t iEle = 0;
   //for ( auto ele : electrons ) 
   for (size_t i = 0; i < electronHandle->size(); ++i)
   {   
        const auto el = electronHandle->ptrAt(i);

        charge.push_back(el->charge());
        energy.push_back(el->energy());
        et.push_back(el->et());
        et2.push_back(el->et2());
        eta.push_back(el->eta());
        mt.push_back(el->mt());
        mtSqr.push_back(el->mtSqr());
        theta.push_back(el->theta());
        phi.push_back(el->phi());
        p.push_back(el->p());
        pt.push_back(el->pt());
        px.push_back(el->px());
        py.push_back(el->py());
        pz.push_back(el->pz());
        //
        // Look up and save the ID decisions
        // 
        bool isPassEleIdLoose  = (*loose_ele_id_decisions)[el];
        bool isPassEleIdMedium  = (*medium_ele_id_decisions)[el];
        bool isPassEleIdTight  = (*tight_ele_id_decisions)[el];
        passEleIdLoose_.push_back  ( (int)isPassEleIdLoose  );
        passEleIdMedium_.push_back  ( (int)isPassEleIdMedium  );
        passEleIdTight_.push_back  ( (int)isPassEleIdTight  ); 
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
DemoAnalyzer::beginJob()
{
   treeEvent->Branch("charge", &charge);
   treeEvent->Branch("energy", &energy);
   treeEvent->Branch("et", &et);
   treeEvent->Branch("et2", &et2);
   treeEvent->Branch("eta", &eta);
   treeEvent->Branch("mt", &mt);
   treeEvent->Branch("mtSqr", &mtSqr);
   treeEvent->Branch("theta", &theta);
   treeEvent->Branch("phi", &phi);
   treeEvent->Branch("p", &p);
   treeEvent->Branch("pt", &pt);
   treeEvent->Branch("px", &px);
   treeEvent->Branch("py", &py);
   treeEvent->Branch("pz", &pz);
   treeEvent->Branch("passEleIdLoose"  ,  &passEleIdLoose_ );
   treeEvent->Branch("passEleIdMedium"  ,  &passEleIdMedium_ );
   treeEvent->Branch("passEleIdTight"  ,  &passEleIdTight_ );
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
