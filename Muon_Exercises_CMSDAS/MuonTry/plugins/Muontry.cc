// -*- C++ -*-
//
// Package:    MuonAnalysis/Muontry
// Class:      Muontry
// 
/**\class Muontry Muontry.cc MuonAnalysis/Muontry/plugins/Muontry.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Tue, 14 Jul 2020 05:15:22 GMT
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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Muontry : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
      explicit Muontry(const edm::ParameterSet&);
      ~Muontry();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;

  TH1F* h_pt_loose;
  TH1F* h_pt_medium;
  TH1F* h_pt_tight;
  TH1F* h_pt_soft;
  TH1F* h_pt_highpt;

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
Muontry::Muontry(const edm::ParameterSet& iConfig):
  muonsToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonTag"))),
  verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexTag")))
{
   //now do what ever initialization is needed
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    h_pt_loose = fs->make<TH1F>("pt_loose",";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_medium = fs->make<TH1F>("pt_medium", ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_tight = fs->make<TH1F>("pt_tight", ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_soft  = fs->make<TH1F>("pt_soft" , ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_highpt = fs->make<TH1F>("pt_highpt", ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
}


Muontry::~Muontry()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Muontry::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  //pat::Muons 
   edm::Handle<edm::View<pat::Muon> > muons;
   iEvent.getByToken(muonsToken_, muons);
   if(!muons.isValid()) 
   {
      throw cms::Exception("Muon collection not valid!");
   }

   //vertices 
   edm::Handle<std::vector<reco::Vertex>> vertices;
   iEvent.getByToken(verticesToken_, vertices);
   if(!vertices.isValid())   
   {
      throw cms::Exception("Vertex collection not valid!");
   }

  // Let's check that we have at least one good vertex! 
  std::vector<reco::Vertex>::const_iterator firstGoodVertex = vertices->end();
  for (std::vector<reco::Vertex>::const_iterator it=vertices->begin(); it!=firstGoodVertex; ++it) 
  {
     if (!it->isFake() && it->ndof()>4 && it->position().Rho()<2. && std::abs(it->position().Z())<24.) 
     {
        if(firstGoodVertex == vertices->end()) firstGoodVertex = it;
        break;
     }
  }
             	                                           
  // Require a good vertex
  if(firstGoodVertex == vertices->end()) return;
  edm::View<pat::Muon>::const_iterator muend = muons->end();
  for (edm::View<pat::Muon>::const_iterator it=muons->begin(); it!=muend; ++it) 
  {
    // Require muon to have a silicon track (i.e. global- or tracker-muon) 
    if (!it->innerTrack().isNonnull()) continue;
    // Require that muon be within the tracker volume and have pt > 20 GeV
    if (std::abs(it->eta())>2.5) continue;
    if (it->pt()<20) continue;
   
    if(it->isLooseMuon()) 
    {
      h_pt_loose->Fill(std::min(it->pt(),99.9));
    }
    if(it->isMediumMuon()) 
    {
      h_pt_medium->Fill(std::min(it->pt(),99.9));
    }
    if(it->isTightMuon(*firstGoodVertex)) 
    {
      h_pt_tight->Fill(std::min(it->pt(),99.9));
    }
    if(it->isSoftMuon(*firstGoodVertex)) 
    {
      h_pt_soft->Fill(std::min(it->pt(),99.9));
    }
    if(it->isHighPtMuon(*firstGoodVertex)) 
    {
      h_pt_highpt->Fill(std::min(it->pt(),99.9));
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
Muontry::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Muontry::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Muontry::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Muontry);
