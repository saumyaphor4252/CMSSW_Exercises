// -*- C++ -*-
//
// Package:    MuonAnalysis/MuonExercise1
// Class:      MuonExercise1
// 
/**\class MuonExercise1 MuonExercise1.cc MuonAnalysis/MuonExercise1/plugins/MuonExercise1.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Fri, 05 Jun 2020 06:01:23 GMT
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

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuonExercise1 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonExercise1(const edm::ParameterSet&);
      ~MuonExercise1();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonsToken_; 
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genmuonsToken_;
      edm::EDGetTokenT<l1extra::L1MuonParticleCollection> l1muonsToken_;
      //edm::EDGetTokenT<BXVector<l1t::Muon>> l1muonsToken_;
      edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjsToken_;

      TH1D *h_pt;
      TH1D *h_eta;
      TH1D *h_phi;
      TH1D *h_nmuonsPerEvent;
      TH1D *h_leading_pt;
      TH1D *h_leading_eta;
      TH1D *h_subleading_pt;
      TH1D *h_subleading_eta;
      TH1D *h_genpt;
      TH1D *h_geneta;
      TH1D *h_genphi;
      TH1D *h_nGenmuonsPerEvent;
      TH1D *h_l1pt;
      TH1D *h_l1eta;
      TH1D *h_l1phi;
      TH1D *h_hltpt;
      TH1D *h_hlteta;
      TH1D *h_hltphi;
      TH1D* h_GenMatch_pt;
      TH1D* h_GenMatch_eta;
      TH1D* h_GenMatch_phi;
      TH1D* h_GenUnmatch_pt;
      TH1D* h_GenUnmatch_eta;
      TH1D* h_GenUnmatch_phi;
      TH1D* h_L1Match_pt;
      TH1D* h_L1Match_eta;                
      TH1D* h_L1Match_phi;            
      TH1D* h_L1Unmatch_pt;
      TH1D* h_L1Unmatch_eta;
      TH1D* h_L1Unmatch_phi;
      TH1D* h_HLTMatch_pt;
      TH1D* h_HLTMatch_eta;
      TH1D* h_HLTMatch_phi;
      TH1D* h_HLTUnmatch_pt;
      TH1D* h_HLTUnmatch_eta;
      TH1D* h_HLTUnmatch_phi;

      size_t nglobalmuons=0,nrecomuons=0,nleadingmuons=0,nsubleadingmuons=0,ngenmuons=0,nl1muons=0,nhltmuons=0,nhlt=0;
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
MuonExercise1::MuonExercise1(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::InputTag muonLabel("slimmedMuons");
   edm::InputTag GenmuonLabel("packedGenParticles");
   edm::InputTag L1muonLabel("l1extraParticles");
   //edm::InputTag L1muonLabel("gmtStage2Digis:Muon");
   edm::InputTag TriggerLabel("TriggerResults","","HLT");
   edm::InputTag TrigObjLabel("selectedPatTrigger");
   //edm::InputTag TrigObjLabel("slimmedPatTrigger");

   muonsToken_ = consumes<pat::MuonCollection>(muonLabel);
   genmuonsToken_ = consumes<pat::PackedGenParticleCollection>(GenmuonLabel);
   l1muonsToken_ = consumes<l1extra::L1MuonParticleCollection>(L1muonLabel);
   //l1muonsToken_ = consumes<BXVector<l1t::Muon>>(L1muonLabel);
   trigResultsToken_ = consumes<edm::TriggerResults>(TriggerLabel);
   trigObjsToken_ = consumes<pat::TriggerObjectStandAloneCollection>(TrigObjLabel);

   edm::Service<TFileService> fs;
   h_pt = fs->make<TH1D>("pt","PAT Pt",100,0.0,100.0);
   h_eta = fs->make<TH1D>("eta","PAT Eta",100,-3.0,3.0);
   h_phi = fs->make<TH1D>("phi","PAT Phi",100,-4.0,4.0);
   h_nmuonsPerEvent=fs->make<TH1D>("nmuonsPerEvent","Muon Multiplicity",10,0,10.0);
   h_leading_pt=fs->make<TH1D>("leading_pt","Leading pT",100,0,100);
   h_leading_eta=fs->make<TH1D>("leading_eta","Leading Eta",100,-3.0,3.0);
   h_subleading_pt=fs->make<TH1D>("subleading_pt","Subleading pT",100,0,100);
   h_subleading_eta=fs->make<TH1D>("subleading_eta","Subleading Eta",100,-3.0,3.0);
   h_genpt = fs->make<TH1D>("genpt","GEN Pt",100,0.0,100.0);
   h_geneta = fs->make<TH1D>("geneta","GEN Eta",100,-3.0,3.0);
   h_genphi = fs->make<TH1D>("genphi","GEN Phi",100,-4,4);
   h_nGenmuonsPerEvent=fs->make<TH1D>("nGenmuonsPerEvent","GenMuon Multiplicity",10,0,10.0);
   h_l1pt = fs->make<TH1D>("l1pt","L1 Pt",100,0.0,200.0);
   h_l1eta = fs->make<TH1D>("l1eta","L1 Eta",100,-3.0,3.0);
   h_l1phi = fs->make<TH1D>("l1phi","L1 Phi",100,-4.0,4.0);
   h_hltpt = fs->make<TH1D>("hltpt","HLT Pt",100,0.0,200.0);
   h_hlteta = fs->make<TH1D>("hlteta","HLT Eta",100,-3.0,3.0);
   h_hltphi = fs->make<TH1D>("hltphi","HLT Phi",100,-4.0,4.0);
   h_GenMatch_pt=fs->make<TH1D>("h_GenMatch_pt","h_GenMatch_pt",100,0,100);
   h_GenMatch_eta=fs->make<TH1D>("h_GenMatch_eta","h_GenMatch_eta",100,-3.0,3.0);
   h_GenMatch_phi=fs->make<TH1D>("h_GenMatch_phi","h_GenMatch_phi",100,-4.0,4.0);
   h_GenUnmatch_pt=fs->make<TH1D>("h_GenUnmatch_pt","h_GenUnmatch_pt",100,0,100);
   h_GenUnmatch_eta=fs->make<TH1D>("h_GenUnmatch_eta","h_GenUnmatch_eta",100,-3.0,3.0);
   h_GenUnmatch_phi=fs->make<TH1D>("h_GenUnmatch_phi","h_GenUnmatch_phi",100,-4.0,4.0);
   h_L1Match_pt=fs->make<TH1D>("h_L1Match_pt","h_L1Match_pt",100,0,100);
   h_L1Match_eta=fs->make<TH1D>("h_L1Match_eta","h_L1Match_eta",100,-3.0,3.0);                  
   h_L1Match_phi=fs->make<TH1D>("h_L1Match_phi","h_L1Match_phi",100,-4.0,4.0);              
   h_L1Unmatch_pt=fs->make<TH1D>("h_L1Unmatch_pt","h_L1Unmatch_pt",100,0,100);
   h_L1Unmatch_eta=fs->make<TH1D>("h_L1Unmatch_eta","h_L1Unmatch_eta",100,-3.0,3.0);
   h_L1Unmatch_phi=fs->make<TH1D>("h_L1Unmatch_phi","h_L1Unmatch_phi",100,-3.0,3.0);
   h_HLTMatch_pt=fs->make<TH1D>("h_HLTMatch_pt","h_HLTMatch_pt",100,0,100);
   h_HLTMatch_eta=fs->make<TH1D>("h_HLTMatch_eta","h_HLTMatch_eta",100,-3.0,3.0);
   h_HLTMatch_phi=fs->make<TH1D>("h_HLTMatch_phi","h_HLTMatch_phi",100,-4.0,4.0);
   h_HLTUnmatch_pt=fs->make<TH1D>("h_HLTUnmatch_pt","h_HLTUnmatch_pt",100,0,100);
   h_HLTUnmatch_eta=fs->make<TH1D>("h_HLTUnmatch_eta","h_HLTUnmatch_eta",100,-3.0,3.0);
   h_HLTUnmatch_phi=fs->make<TH1D>("h_HLTUnmatch_phi","h_HLTUnmatch_phi",100,-4.0,4.0);

}


MuonExercise1::~MuonExercise1()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   //
   // PAT Muons
   //
   edm::Handle<vector<pat::Muon> > muons;
   iEvent.getByToken(muonsToken_,muons);
 
   h_nmuonsPerEvent->Fill(muons->size());
 
   size_t leading_muon=0;   
   for(auto it=muons->begin();it!=muons->end();++it)
   {
     nrecomuons++;leading_muon++;
     if(!it->isGlobalMuon()) continue;
     nglobalmuons++;
     reco::TrackRef mu_g(it->globalTrack());
     reco::TrackRef mu_i(it->innerTrack());
     reco::TrackRef mu_o(it->outerTrack());
     h_pt->Fill(it->pt());
     h_eta->Fill(it->eta());
     h_phi->Fill(it->phi()); 
    
     // Leading Muons 
     if(leading_muon==1)
     {  
        nleadingmuons++;
        h_leading_pt->Fill(it->pt());
        h_leading_eta->Fill(it->eta());
     }
     // Sub Leading Muons
     if(leading_muon==2)
     {
        nsubleadingmuons++;  
        h_subleading_pt->Fill(it->pt());
        h_subleading_eta->Fill(it->eta());
     }
    
   }

   //
   // GEN Muons
   //
   //if(!iEvent.isRealData())
   edm::Handle<pat::PackedGenParticleCollection> genmuons;
   iEvent.getByToken(genmuonsToken_,genmuons);
 
   int genmuonsPerEvent=0;
   for(auto it=genmuons->begin();it!=genmuons->end();++it)
   {
     const pat::PackedGenParticle& mcMuon=(*it);
     if(abs(mcMuon.pdgId())!=13) continue;
     if(fabs(mcMuon.eta())>2.4 || mcMuon.pt()<1.5) continue;
     ngenmuons++;genmuonsPerEvent++;
     h_genpt->Fill(mcMuon.pt());
     h_geneta->Fill(mcMuon.eta());
     h_genphi->Fill(mcMuon.phi());
   }
   h_nGenmuonsPerEvent->Fill(genmuonsPerEvent);

   //
   // L1 Trigger Muons
   //
   edm::Handle<l1extra::L1MuonParticleCollection> l1muons;
   //edm::Handle<BXVector<l1t::Muon>> l1muons;
   iEvent.getByToken(l1muonsToken_,l1muons);
   
   for(auto it=l1muons->begin();it!=l1muons->end();++it)
   {
      const L1MuGMTExtendedCand l1muon = (*it).gmtMuonCand();
      //if(l1muons.empty()) continue;
      nl1muons++;
      float pt=(*it).pt();
      float eta=(*it).eta();
      float phi=(*it).phi();
      h_l1pt->Fill(pt);
      h_l1eta->Fill(eta);
      h_l1phi->Fill(phi);
   }
  
   //
   // HLT Muons
   //
   edm::Handle<edm::TriggerResults> triggerBits;
   iEvent.getByToken(trigResultsToken_,triggerBits);
   const edm::TriggerNames& names=iEvent.triggerNames(*triggerBits);

   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(trigObjsToken_,triggerObjects);
   
   for(auto it=triggerObjects->cbegin();it!=triggerObjects->cend();++it)
   { if((*it).hasTriggerObjectType(trigger::TriggerMuon))
     nhlt++;
   }

   for(pat::TriggerObjectStandAlone hltMuon:*triggerObjects)
   {
      hltMuon.unpackPathNames(names);
      if(!hltMuon.hasTriggerObjectType(trigger::TriggerMuon)) continue;
      nhltmuons++;
      h_hltpt->Fill(hltMuon.pt());
      h_hlteta->Fill(hltMuon.eta());
      h_hltphi->Fill(hltMuon.phi());
   }

   //
   // Matching
   //
   for(auto it1=muons->cbegin();it1!=muons->cend();++it1)
   {
       if(! (it1->isGlobalMuon()) ) continue;
       int idx_G = -1;
       int idx_L1 = -1;
       int idx_HLT = -1;

       // match gen
       for(auto it2=genmuons->cbegin();it2!=genmuons->cend(); ++it2)
       {   
          const pat::PackedGenParticle& genParticle= (*it2);
          if(abs(genParticle.pdgId())!=13) continue;
          if(fabs(genParticle.eta())>2.4) continue;
          double dR=deltaR(*it2,*(it1->innerTrack()));
          if(dR<0.1)
          {
           idx_G=std::distance(genmuons->cbegin(), it2);
          }
       }
       // match trigger
       int n=0;
       for(pat::TriggerObjectStandAlone obj: *triggerObjects)
       {
          obj.unpackPathNames(names);
          if(not (obj.hasTriggerObjectType(trigger::TriggerMuon) || obj.hasTriggerObjectType(trigger::TriggerL1Mu) )) continue;
          double dR = deltaR(obj, *(it1->innerTrack()));
   
          // L1
          if(obj.hasTriggerObjectType(trigger::TriggerL1Mu) && dR<0.25)
          {
           idx_L1 = n;
          }
          // HLT
          if(obj.hasTriggerObjectType(trigger::TriggerMuon) && dR<0.1)
          {
           idx_HLT = n;
          }
          n++;
       }

       if(idx_G>-1) 
       {   h_GenMatch_pt->Fill(it1->pt());
           h_GenMatch_eta->Fill(it1->eta());
           h_GenMatch_phi->Fill(it1->phi());
       }
       else if(idx_G==-1)
       {   h_GenUnmatch_pt->Fill(it1->pt());
           h_GenUnmatch_eta->Fill(it1->eta());
           h_GenUnmatch_phi->Fill(it1->phi());
       }

       if(idx_L1>-1)
       {   h_L1Match_pt->Fill(it1->pt()); 
           h_L1Match_eta->Fill(it1->eta()); 
           h_L1Match_phi->Fill(it1->phi());                            
       }
       else if(idx_L1==-1)
       {   h_L1Unmatch_pt->Fill(it1->pt());                      
           h_L1Unmatch_eta->Fill(it1->eta());                      
           h_L1Unmatch_phi->Fill(it1->phi());                      
       }
  
       if(idx_HLT>-1)
       {   h_HLTMatch_pt->Fill(it1->pt()); 
           h_HLTMatch_eta->Fill(it1->eta()); 
           h_HLTMatch_phi->Fill(it1->phi()); 
       }
       else if(idx_HLT==-1)
       {   h_HLTUnmatch_pt->Fill(it1->pt());
           h_HLTUnmatch_eta->Fill(it1->eta());
           h_HLTUnmatch_phi->Fill(it1->phi());
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
MuonExercise1::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise1::endJob() 
{
   std::cout<<" Total Reco Muons: "<<nrecomuons<<std::endl;
   std::cout<<" Total Global Muons: "<<nglobalmuons<<std::endl;    
   std::cout<<" Total Leading Global Muons: "<<nleadingmuons<<std::endl;
   std::cout<<" Total Sunleading Global Muons: "<<nsubleadingmuons<<std::endl;
   std::cout<<" Total Gen Muons: "<<ngenmuons<<std::endl;
   std::cout<<" Total L1 Trigger Muons: "<<nl1muons<<std::endl;
   std::cout<<" Total HLT Muons: "<<nhltmuons<<std::endl;
   std::cout<<" Total HLT Muons new: "<<nhlt<<std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonExercise1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise1);
