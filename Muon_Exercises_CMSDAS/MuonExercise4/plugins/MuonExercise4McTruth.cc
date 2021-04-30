// -*- C++ -*-
//
// Package:    CMSDASExercises/MuonExercise4McTruth
// Class:      MuonExercise4McTruth
// 
/**\class MuonExercise4McTruth MuonExercise4McTruth.cc CMSDASExercises/MuonExercise4McTruth/plugins/MuonExercise4McTruth.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Daniele Trocino
//         Created:  Tue, 29 Aug 2017 13:13:55 GMT
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" 

#include "TH1.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuonExercise4McTruth : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MuonExercise4McTruth(const edm::ParameterSet&);
  ~MuonExercise4McTruth();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  // Token declarations
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genCollToken;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonCollToken;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertexCollToken;

  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;

  // TFileService
  edm::Service<TFileService> fs;

  std::vector<std::string> triggerList;
  std::map<TString, size_t> triggerIdxList; 

  // Trigger matching 
  bool matchTriggerObject(const pat::Muon& mu,
			  const std::map<TString, size_t> & trigList,
			  const edm::TriggerNames& trigNames, 
			  edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjs,
			  edm::Handle<edm::TriggerResults> trigBits);

  // Histos
  TH1D *h_pt; 
  TH1D *h_eta;
  TH1D *h_nvtx;

  TH1D *h_rec_pt; 
  TH1D *h_rec_eta; 
  TH1D *h_rec_nvtx; 

  TH1D *h_rec_loose_pt; 
  TH1D *h_rec_loose_eta; 
  TH1D *h_rec_loose_nvtx; 

  TH1D *h_rec_medium_pt; 
  TH1D *h_rec_medium_eta; 
  TH1D *h_rec_medium_nvtx; 

  TH1D *h_rec_tight_pt; 
  TH1D *h_rec_tight_eta; 
  TH1D *h_rec_tight_nvtx; 

  TH1D *h_rec_tight_iso_pt; 
  TH1D *h_rec_tight_iso_eta; 
  TH1D *h_rec_tight_iso_nvtx; 

  TH1D *h_rec_tight_iso_hlt_pt; 
  TH1D *h_rec_tight_iso_hlt_eta; 
  TH1D *h_rec_tight_iso_hlt_nvtx; 

  // Efficiencies
  TGraphAsymmErrors *gae_rec_pt; 
  TGraphAsymmErrors *gae_rec_eta; 
  TGraphAsymmErrors *gae_rec_nvtx; 

  TGraphAsymmErrors *gae_rec_loose_pt; 
  TGraphAsymmErrors *gae_rec_loose_eta; 
  TGraphAsymmErrors *gae_rec_loose_nvtx; 

  TGraphAsymmErrors *gae_rec_medium_pt; 
  TGraphAsymmErrors *gae_rec_medium_eta; 
  TGraphAsymmErrors *gae_rec_medium_nvtx; 

  TGraphAsymmErrors *gae_rec_tight_pt; 
  TGraphAsymmErrors *gae_rec_tight_eta; 
  TGraphAsymmErrors *gae_rec_tight_nvtx; 

  TGraphAsymmErrors *gae_rec_tight_iso_pt; 
  TGraphAsymmErrors *gae_rec_tight_iso_eta; 
  TGraphAsymmErrors *gae_rec_tight_iso_nvtx; 

  TGraphAsymmErrors *gae_rec_tight_iso_hlt_pt; 
  TGraphAsymmErrors *gae_rec_tight_iso_hlt_eta; 
  TGraphAsymmErrors *gae_rec_tight_iso_hlt_nvtx; 
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
MuonExercise4McTruth::MuonExercise4McTruth(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  usesResource("TFileService");

  edm::InputTag genPartTag("prunedGenParticles");
  edm::InputTag muonTag("slimmedMuons");
  edm::InputTag vertexTag("offlineSlimmedPrimaryVertices");
  edm::InputTag triggerTag("TriggerResults", "", "HLT");
  edm::InputTag trigObjTag("selectedPatTrigger");

  genCollToken = consumes<reco::GenParticleCollection>(genPartTag); 
  muonCollToken = consumes<pat::MuonCollection>(muonTag);
  vertexCollToken = consumes<std::vector<reco::Vertex> >(vertexTag);
  trigResultsToken = consumes<edm::TriggerResults>(triggerTag);
  trigObjCollToken = consumes<pat::TriggerObjectStandAloneCollection>(trigObjTag);

  triggerList = iConfig.getParameter<std::vector<std::string> >("TriggerList");

  // Histograms 
  h_pt   = fs->make<TH1D>("pt"  , ";GEN muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_eta  = fs->make<TH1D>("eta" , ";GEN muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_nvtx = fs->make<TH1D>("nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 

  h_rec_pt   = fs->make<TH1D>("rec_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_rec_eta  = fs->make<TH1D>("rec_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_rec_nvtx = fs->make<TH1D>("rec_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 
        		           
  h_rec_loose_pt   = fs->make<TH1D>("rec_loose_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_rec_loose_eta  = fs->make<TH1D>("rec_loose_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_rec_loose_nvtx = fs->make<TH1D>("rec_loose_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 
              		           	       
  h_rec_medium_pt   = fs->make<TH1D>("rec_medium_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_rec_medium_eta  = fs->make<TH1D>("rec_medium_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_rec_medium_nvtx = fs->make<TH1D>("rec_medium_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 
              		           	       
  h_rec_tight_pt   = fs->make<TH1D>("rec_tight_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_rec_tight_eta  = fs->make<TH1D>("rec_tight_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_rec_tight_nvtx = fs->make<TH1D>("rec_tight_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 
              		           	       
  h_rec_tight_iso_pt   = fs->make<TH1D>("rec_tight_iso_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_rec_tight_iso_eta  = fs->make<TH1D>("rec_tight_iso_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_rec_tight_iso_nvtx = fs->make<TH1D>("rec_tight_iso_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 
              		           	       
  h_rec_tight_iso_hlt_pt   = fs->make<TH1D>("rec_tight_iso_hlt_pt"  , ";PAT muon p_{T} [GeV];Events",  20,  20.0, 120.0); 
  h_rec_tight_iso_hlt_eta  = fs->make<TH1D>("rec_tight_iso_hlt_eta" , ";PAT muon #eta;Events"       ,  24,  -2.4,   2.4); 
  h_rec_tight_iso_hlt_nvtx = fs->make<TH1D>("rec_tight_iso_hlt_nvtx", ";N. vertices;Events"         ,  60,   0.5,  60.5); 

  // Efficiencies
  unsigned int nbinpt   = h_rec_pt->GetNbinsX(); 
  unsigned int nbineta  = h_rec_eta->GetNbinsX(); 
  unsigned int nbinnvtx = h_rec_nvtx->GetNbinsX(); 
  
  gae_rec_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_pt->SetName("eff_rec_pt");
  gae_rec_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_eta->SetName("eff_rec_eta");
  gae_rec_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_nvtx->SetName("eff_rec_nvtx");

  gae_rec_loose_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_loose_pt->SetName("eff_rec_loose_pt");
  gae_rec_loose_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_loose_eta->SetName("eff_rec_loose_eta");
  gae_rec_loose_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_loose_nvtx->SetName("eff_rec_loose_nvtx");

  gae_rec_medium_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_medium_pt->SetName("eff_rec_medium_pt");  
  gae_rec_medium_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_medium_eta->SetName("eff_rec_medium_eta");
  gae_rec_medium_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_medium_nvtx->SetName("eff_rec_medium_nvtx");

  gae_rec_tight_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_tight_pt->SetName("eff_rec_tight_pt");
  gae_rec_tight_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_tight_eta->SetName("eff_rec_tight_eta");
  gae_rec_tight_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_tight_nvtx->SetName("eff_rec_tight_nvtx");

  gae_rec_tight_iso_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_tight_iso_pt->SetName("eff_rec_tight_iso_pt");
  gae_rec_tight_iso_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_tight_iso_eta->SetName("eff_rec_tight_iso_eta");
  gae_rec_tight_iso_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_tight_iso_nvtx->SetName("eff_rec_tight_iso_nvtx");

  gae_rec_tight_iso_hlt_pt   = fs->make<TGraphAsymmErrors>(nbinpt  ); gae_rec_tight_iso_hlt_pt->SetName("eff_rec_tight_iso_hlt_pt");
  gae_rec_tight_iso_hlt_eta  = fs->make<TGraphAsymmErrors>(nbineta ); gae_rec_tight_iso_hlt_eta->SetName("eff_rec_tight_iso_hlt_eta");
  gae_rec_tight_iso_hlt_nvtx = fs->make<TGraphAsymmErrors>(nbinnvtx); gae_rec_tight_iso_hlt_nvtx->SetName("eff_rec_tight_iso_hlt_nvtx");

}


MuonExercise4McTruth::~MuonExercise4McTruth()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise4McTruth::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(vertexCollToken, vertices);
  if(!vertices.isValid())
  {
    throw cms::Exception("Vertex collection not valid!"); 
  } 

  int nGoodVtx = 0;
  const reco::Vertex *goodVtx = nullptr; 
  for(std::vector<reco::Vertex>::const_iterator it=vertices->begin(), endVtx = vertices->end(); it!=endVtx; ++it)
  {
    if(!it->isFake() && it->ndof()>4 && it->position().Rho()<2. && std::abs(it->position().Z())<24.) 
    {
      nGoodVtx++;
      if(goodVtx==nullptr) goodVtx = &(*it); 
    }
  }
  // Require a good vertex
  if(nGoodVtx==0) return; 

  // Retrieve the GenParticle collection and loop over it 
  edm::Handle<reco::GenParticleCollection> genColl; 
  iEvent.getByToken(genCollToken, genColl); 
  if(!genColl.isValid())
  {
    throw cms::Exception("GenParticle collection not valid!"); 
  }

  std::vector<const reco::GenParticle*> genmus; 
  for(auto&& genPart : *(genColl.product()))
  { 
     // Check if it's a muon from Drell-Yan process
     if(genPart.isPromptFinalState() && std::abs(genPart.pdgId()) == 13) 
     {  
        // Only muons within acceptance and pt>20 
        if(genPart.pt()>20. && std::abs(genPart.eta())<2.4)	
        {
           // Fill the collection of selected gen muons 
           genmus.push_back(&genPart); 

	   // Fill the gen-muon plots 
           h_pt->Fill(genPart.pt());
           h_eta->Fill(genPart.eta()); 
           h_nvtx->Fill(nGoodVtx);
        }
     }
  }

  if(genmus.size()==0) return;
  std::vector<const reco::GenParticle*>::const_iterator gmbeg = genmus.begin(), gmend = genmus.end();

  // --- Now prepare all the trigger objects --- 
  // Retrieve TriggerResults 
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(trigResultsToken, triggerBits);
  if(!triggerBits.isValid())
  {
    throw cms::Exception("TriggerResults collection not valid!"); 
  } 
  
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerBits);

  // Find trigger indexes (only once)
  if(triggerIdxList.size()==0)
  { 
     for(auto&& t : triggerList)
     { 
         for(size_t i=0, n=triggerBits->size(); i<n; ++i)
         {
       	     if(TString(triggerNames.triggerName(i)).Contains((t+"_v").c_str()))
             { 
                 triggerIdxList[(t+"_v").c_str()] = i;
	         break; 
             } 
         }
     }
  }

  // Retrieve TriggerObjects 
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigObjCollToken, triggerObjects);
  if(!triggerObjects.isValid())
  {
    throw cms::Exception("TriggerObjectStandAloneCollection collection not valid!"); 
  } 

  // Retrieve the pat::Muon collection and loop over it 
  edm::Handle<pat::MuonCollection> muonColl; 
  iEvent.getByToken(muonCollToken, muonColl); 
  if(!muonColl.isValid())
  {
    throw cms::Exception("Muon collection not valid!"); 
  }

  for(auto&& muon : *(muonColl.product())) 
  {
    // Let's skip muons that are standalone only
    if( muon.isStandAloneMuon() && !(muon.isTrackerMuon() || muon.isGlobalMuon()) ) continue;

    // Check if it is matched to a GenParticle
    const reco::GenParticle *gmu = muon.genParticle(); 
    if( gmu==0 ) continue;

    // Check if the GenParticle is among the ones we selected 
    if(std::find(gmbeg, gmend, gmu)==gmend) continue;  

    // Fill the plots here 
    // (Let's not cut on pat::Muon pt and |eta| (already cut at GEN level))
    h_rec_pt->Fill(muon.pt());
    h_rec_eta->Fill(muon.eta());
    h_rec_nvtx->Fill(nGoodVtx);
 
    // -- ID --
    //  Fill ID plots 
    //  - Loose 
    if(muon.isLooseMuon())
    {
       h_rec_loose_pt->Fill(muon.pt());
       h_rec_loose_eta->Fill(muon.eta());
       h_rec_loose_nvtx->Fill(nGoodVtx);
    }

    //  - Medium 
    if(muon.isMediumMuon())
    {
       h_rec_medium_pt->Fill(muon.pt());
       h_rec_medium_eta->Fill(muon.eta());
       h_rec_medium_nvtx->Fill(nGoodVtx);
    }

    //  - Tight
    if( muon.isTightMuon(*goodVtx)==false ) continue; // if it's not tight, nothing else to do 
    h_rec_tight_pt->Fill(muon.pt());
    h_rec_tight_eta->Fill(muon.eta());
    h_rec_tight_nvtx->Fill(nGoodVtx);

    // Now isolation, only for tight muons
    if(muon.isIsolationValid()==false) continue; 
    const reco::MuonPFIsolation &pfR04 = muon.pfIsolationR04();

    // Calculate PF combined relative isolation with Delta-beta correction 
    double corriso=pfR04.sumChargedHadronPt+std::max(0.0,double(pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt)) ;

    if(corriso/muon.pt()>0.12) continue; // not isolated, nothing else to do
    h_rec_tight_iso_pt->Fill(muon.pt());
    h_rec_tight_iso_eta->Fill(muon.eta());
    h_rec_tight_iso_nvtx->Fill(nGoodVtx);

    // Finally, let's see if the isolated tight muon fired a trigger
    bool passTrigger = matchTriggerObject(muon, triggerIdxList, triggerNames, triggerObjects, triggerBits);
    if(passTrigger==false) continue; 
    h_rec_tight_iso_hlt_pt->Fill(muon.pt());
    h_rec_tight_iso_hlt_eta->Fill(muon.eta());
    h_rec_tight_iso_hlt_nvtx->Fill(nGoodVtx);

  } // end for(auto&& muon : *(muonColl.product()))
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonExercise4McTruth::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise4McTruth::endJob() 
{
  // Compute efficiencies
  try{gae_rec_pt  ->Divide(h_rec_pt  , h_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_eta ->Divide(h_rec_eta , h_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_nvtx->Divide(h_rec_nvtx, h_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  

  try{gae_rec_loose_pt  ->Divide(h_rec_loose_pt  , h_rec_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_loose_eta ->Divide(h_rec_loose_eta , h_rec_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_loose_nvtx->Divide(h_rec_loose_nvtx, h_rec_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  

  try{gae_rec_medium_pt  ->Divide(h_rec_medium_pt  , h_rec_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_medium_eta ->Divide(h_rec_medium_eta , h_rec_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_medium_nvtx->Divide(h_rec_medium_nvtx, h_rec_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  

  try{gae_rec_tight_pt  ->Divide(h_rec_tight_pt  , h_rec_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_tight_eta ->Divide(h_rec_tight_eta , h_rec_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_tight_nvtx->Divide(h_rec_tight_nvtx, h_rec_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  

  try{gae_rec_tight_iso_pt  ->Divide(h_rec_tight_iso_pt  , h_rec_tight_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_tight_iso_eta ->Divide(h_rec_tight_iso_eta , h_rec_tight_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_tight_iso_nvtx->Divide(h_rec_tight_iso_nvtx, h_rec_tight_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  

  try{gae_rec_tight_iso_hlt_pt  ->Divide(h_rec_tight_iso_hlt_pt  , h_rec_tight_iso_pt  , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_tight_iso_hlt_eta ->Divide(h_rec_tight_iso_hlt_eta , h_rec_tight_iso_eta , "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {}  
  try{gae_rec_tight_iso_hlt_nvtx->Divide(h_rec_tight_iso_hlt_nvtx, h_rec_tight_iso_nvtx, "cl=0.683 b(1,1) mode");} catch(cms::Exception& ex) {} 

  // Labels 
  gae_rec_pt  ->GetXaxis()->SetTitle(h_rec_pt  ->GetXaxis()->GetTitle()); 
  gae_rec_eta ->GetXaxis()->SetTitle(h_rec_eta ->GetXaxis()->GetTitle()); 
  gae_rec_nvtx->GetXaxis()->SetTitle(h_rec_nvtx->GetXaxis()->GetTitle()); 

  gae_rec_loose_pt  ->GetXaxis()->SetTitle(h_rec_loose_pt  ->GetXaxis()->GetTitle()); 
  gae_rec_loose_eta ->GetXaxis()->SetTitle(h_rec_loose_eta ->GetXaxis()->GetTitle()); 
  gae_rec_loose_nvtx->GetXaxis()->SetTitle(h_rec_loose_nvtx->GetXaxis()->GetTitle()); 

  gae_rec_medium_pt  ->GetXaxis()->SetTitle(h_rec_medium_pt  ->GetXaxis()->GetTitle()); 
  gae_rec_medium_eta ->GetXaxis()->SetTitle(h_rec_medium_eta ->GetXaxis()->GetTitle()); 
  gae_rec_medium_nvtx->GetXaxis()->SetTitle(h_rec_medium_nvtx->GetXaxis()->GetTitle()); 

  gae_rec_tight_pt  ->GetXaxis()->SetTitle(h_rec_tight_pt  ->GetXaxis()->GetTitle()); 
  gae_rec_tight_eta ->GetXaxis()->SetTitle(h_rec_tight_eta ->GetXaxis()->GetTitle()); 
  gae_rec_tight_nvtx->GetXaxis()->SetTitle(h_rec_tight_nvtx->GetXaxis()->GetTitle()); 

  gae_rec_tight_iso_pt  ->GetXaxis()->SetTitle(h_rec_tight_iso_pt  ->GetXaxis()->GetTitle()); 
  gae_rec_tight_iso_eta ->GetXaxis()->SetTitle(h_rec_tight_iso_eta ->GetXaxis()->GetTitle()); 
  gae_rec_tight_iso_nvtx->GetXaxis()->SetTitle(h_rec_tight_iso_nvtx->GetXaxis()->GetTitle()); 

  gae_rec_tight_iso_hlt_pt  ->GetXaxis()->SetTitle(h_rec_tight_iso_hlt_pt  ->GetXaxis()->GetTitle()); 
  gae_rec_tight_iso_hlt_eta ->GetXaxis()->SetTitle(h_rec_tight_iso_hlt_eta ->GetXaxis()->GetTitle()); 
  gae_rec_tight_iso_hlt_nvtx->GetXaxis()->SetTitle(h_rec_tight_iso_hlt_nvtx->GetXaxis()->GetTitle()); 

  gae_rec_pt  ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_eta ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_nvtx->GetYaxis()->SetTitle("Efficiency"); 

  gae_rec_loose_pt  ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_loose_eta ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_loose_nvtx->GetYaxis()->SetTitle("Efficiency"); 

  gae_rec_medium_pt  ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_medium_eta ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_medium_nvtx->GetYaxis()->SetTitle("Efficiency"); 

  gae_rec_tight_pt  ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_tight_eta ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_tight_nvtx->GetYaxis()->SetTitle("Efficiency"); 

  gae_rec_tight_iso_pt  ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_tight_iso_eta ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_tight_iso_nvtx->GetYaxis()->SetTitle("Efficiency"); 

  gae_rec_tight_iso_hlt_pt  ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_tight_iso_hlt_eta ->GetYaxis()->SetTitle("Efficiency"); 
  gae_rec_tight_iso_hlt_nvtx->GetYaxis()->SetTitle("Efficiency"); 

}

bool MuonExercise4McTruth::matchTriggerObject(const pat::Muon& mu,
					      const std::map<TString, size_t>& trigList,
					      const edm::TriggerNames& trigNames, 
					      edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjs,
					      edm::Handle<edm::TriggerResults> trigBits)
{
   double dR = 0.1;

   for(auto&& t : trigList)
   {
      if(trigBits->accept(t.second)==false) continue; 
      for(pat::TriggerObjectStandAlone obj : *trigObjs)
      {
         obj.unpackPathNames(trigNames);
         if(obj.hasPathName((t.first+"*").Data(), true, true))
         { 
            if(reco::deltaR(obj, *(mu.innerTrack())) < dR) 
            {
	       return true; 
            }
         }
      } // end for(pat::TriggerObjectStandAlone obj : *trigObjs)
   } // end for(auto&& t : trigList)

  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise4McTruth);
