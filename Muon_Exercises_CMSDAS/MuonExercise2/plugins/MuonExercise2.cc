// -*- C++ -*-
//
// Package:    MuonAnalysis/MuonExercise2
// Class:      MuonExercise2
// 
/**\class MuonExercise2 MuonExercise2.cc MuonAnalysis/MuonExercise2/plugins/MuonExercise2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Thu, 11 Jun 2020 14:59:13 GMT
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

class MuonExercise2 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonExercise2(const edm::ParameterSet&);
      ~MuonExercise2();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

 //     RoccoR rc("//eos/uscms/store/user/cmsdas/2019/short_exercises/Muons/Exercises/MuonExercise2/plugin/rcdata.2016.v3"); 
 //     TRandom3 rnd(1234);
      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::MuonCollection> muonsToken;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genMuonsToken;

        // Histograms
        TH1F* h_RecDiMuonM;
        TH1F* h_GenDiMuonM;
        TH1F* h_MassResolution;
        TH1F* h_MupTResolution;

        // Profile Histograms Declaration 
        //---------For Scale Estimation----------
        TProfile* prof_DiMuonMVsMuPlusPhi;
        TProfile* prof_DiMuonMVsMuMinusPhi;
        TProfile* prof_DiMuonMVsMuonEta;

        //---------For REsolution Study----------
        TProfile* prof_pTResolutionVsEta;
        TProfile* prof_pTResolutionVsMupT;
      
        size_t nPatTotal, nPatGlobalCut, nPatPtCut, nPatEtaCut, nPatHadronChargedIsoCut, nVectorSizeCut, nOppositeChargeCut, nGenTotal, nGenEtaCut,
               nGenPtCut, nGenMatchedMuonPairs, nRecMassCut, nRecProfileEntries, nResolutionHistogramEntries, nResolutionTProfileEntries;


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
MuonExercise2::MuonExercise2(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::InputTag muonLabel("slimmedMuons");
   edm::InputTag genMuonLabel ("packedGenParticles");

   muonsToken=consumes<pat::MuonCollection>(muonLabel); 
   genMuonsToken=consumes<pat::PackedGenParticleCollection>(genMuonLabel);

   edm::Service<TFileService>fs;
   
   h_RecDiMuonM=fs->make<TH1F>("h_RecDiMuonM",";m_{#mu^{+}#mu^{-}};",80,70,110);
   h_GenDiMuonM=fs->make<TH1F>("h_GemDiMuonM",";m_{#mu^{+}#mu^{-}};",80,70,110);
   h_MassResolution=fs->make<TH1F>("h_MassResolution","Mass Resolution",80,-0.15,0.15);
   h_MupTResolution=fs->make<TH1F>("h_MupTResolution","Muon p_{T} Resolution;#Delta p_{T}/p_{T};",80,-0.2,0.2);

   prof_DiMuonMVsMuPlusPhi=fs->make<TProfile>("prof_DiMuonMVsMuPlusPhi","m_{#mu^{+}#mu^{-}} Vs #mu^{+}#phi; RecoMuon{+} #phi[rad];Z Peak position[GeV/c^2]",16,-3.14,3.14,88,93);
   prof_DiMuonMVsMuMinusPhi=fs->make<TProfile>("prof_DiMuonMVsMuMinusPhi","m_{#mu^{+}#mu^{-}} Vs #mu^{-}#phi; RecoMuon{-} #phi[rad];Z Peak position[GeV/c^2]",16,-3.14,3.14,88,93);
   prof_DiMuonMVsMuonEta=fs->make<TProfile>("prof_DiMuonMVsMuEta","m_{#mu^{+}#mu^{-}} Vs Muon #eta; RecoMuon #eta;Z Peak position[GeV/c^2]",50,-2.4,2.4,88,93);
   
   prof_pTResolutionVsEta=fs->make<TProfile>("prof_pTResolutionVsEta","p_{T} Resolution;Gen Muon #eta;#Delta p_{T}/p_{T}",50,-2.4,2.4,0,1);
   prof_pTResolutionVsMupT=fs->make<TProfile>("prof_pTResolutionVsMupT","p_{T} Resolution;Gen Muon p_{T} [GEV];#Delta p_{T}/p_{T}",25,20,100,0,1);

   nPatTotal=0; nPatGlobalCut=0; nPatPtCut=0; nPatEtaCut=0; nPatHadronChargedIsoCut=0; nVectorSizeCut=0; nOppositeChargeCut=0; nGenTotal=0; nGenEtaCut=0;
   nGenPtCut=0; nGenMatchedMuonPairs=0; nRecMassCut=0; nRecProfileEntries=0; nResolutionHistogramEntries=0; nResolutionTProfileEntries=0;
}
MuonExercise2::~MuonExercise2()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonsToken,muons);

   edm::Handle<pat::PackedGenParticleCollection> genMuons;
   iEvent.getByToken(genMuonsToken,genMuons);
   
   //Loop over muons store the indices of muons passing selection criteria in a vector v with charge of first being + and second -
   std::vector<int> v;
   for(auto mu=muons->cbegin();mu!=muons->cend();++mu)
   {  
      nPatTotal++; 
      if(! mu->isGlobalMuon()) continue;
      nPatGlobalCut++;
      if(! (mu->pt()>20) ) continue;
      nPatPtCut++;
      if(fabs(mu->eta())>2.4) continue;
      nPatEtaCut++;
      if(! (mu->chargedHadronIso()<0.15) ) continue;
      nPatHadronChargedIsoCut++;      
      
      v.push_back( std::distance(muons->cbegin(),mu) );
   }

   if(v.size()==2) 
   {
       nVectorSizeCut++;
       if( muons->at(v[0]).charge()*muons->at(v[1]).charge()==-1 ) 
       { 
           nOppositeChargeCut++;
           int idx_mum=-1,idx_mup=-1;
           if( muons->at(v[0]).charge()>0 ) { idx_mup=v[0]; idx_mum=v[1]; }
           else if( muons->at(v[0]).charge()<0 ) { idx_mup=v[1]; idx_mum=v[0]; }  
    
           int idxmup_G=-1,idxmum_G=-1;
           for(auto mu_G=genMuons->cbegin();mu_G!=genMuons->cend();++mu_G)
           {   
               const pat::PackedGenParticle& mcMuon=(*mu_G);
               if(! (fabs(mcMuon.pdgId())==13)) continue;
               nGenTotal++;
               if(fabs(mcMuon.eta())>2.4) continue;
               nGenEtaCut++;
               if(! (mcMuon.pt()>1.5) ) continue;
               nGenPtCut++;
               if( deltaR(mcMuon,*(muons->at(idx_mup).innerTrack()))<0.1 && mcMuon.charge()>0) idxmup_G=std::distance(genMuons->cbegin(),mu_G);
               if( deltaR(mcMuon,*(muons->at(idx_mum).innerTrack()))<0.1 && mcMuon.charge()<0) idxmum_G=std::distance(genMuons->cbegin(),mu_G);
           }
    
           //Get Gen Dimuon Invariant Mass
           if(idxmup_G>-1 && idxmum_G>-1)
           {
               nGenMatchedMuonPairs++;

               // Get Reconstructed Dimuon Invariant Mass
               double MRec=( muons->at(idx_mup).p4()+ muons->at(idx_mum).p4() ).M();
               if( MRec>70 && MRec<110 ) 
               {
               nRecMassCut++;
               h_RecDiMuonM->Fill(MRec);
          
               //Fill TProfile Histograms
               prof_DiMuonMVsMuPlusPhi->Fill(muons->at(idx_mup).phi(),MRec,1);
               prof_DiMuonMVsMuMinusPhi->Fill(muons->at(idx_mum).phi(),MRec,1);
               prof_DiMuonMVsMuonEta->Fill(muons->at(idx_mup).eta(),MRec,1);
               nRecProfileEntries++;
   
               //Fill Resolution Histograms
               double MGen=(genMuons->at(idxmup_G).p4()+genMuons->at(idxmum_G).p4()).M();
               h_GenDiMuonM->Fill(MGen);
               h_MassResolution->Fill((MRec-MGen)/MGen);
               h_MupTResolution->Fill((1/muons->at(idx_mup).pt()-1/genMuons->at(idxmup_G).pt())/(1/genMuons->at(idxmup_G).pt()));
               h_MupTResolution->Fill((1/muons->at(idx_mum).pt()-1/genMuons->at(idxmum_G).pt())/(1/genMuons->at(idxmum_G).pt()));
               nResolutionHistogramEntries++;

               // Fill Resoltuion TProfiles
               prof_pTResolutionVsEta->Fill(genMuons->at(idxmup_G).eta(),(1/muons->at(idx_mup).pt()-1/genMuons->at(idxmup_G).pt())/(1/genMuons->at(idxmup_G).pt())); 
               prof_pTResolutionVsEta->Fill(genMuons->at(idxmum_G).eta(),(1/muons->at(idx_mum).pt()-1/genMuons->at(idxmum_G).pt())/(1/genMuons->at(idxmum_G).pt()));
               prof_pTResolutionVsMupT->Fill(genMuons->at(idxmup_G).pt(),(1/muons->at(idx_mup).pt()-1/genMuons->at(idxmup_G).pt())/(1/genMuons->at(idxmup_G).pt()));
               prof_pTResolutionVsMupT->Fill(genMuons->at(idxmum_G).pt(),(1/muons->at(idx_mum).pt()-1/genMuons->at(idxmum_G).pt())/(1/genMuons->at(idxmum_G).pt()));
               nResolutionTProfileEntries++;
               }
           }
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
MuonExercise2::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise2::endJob() 
{
  std::cout<<"nPatTotal: "<<nPatTotal<<std::endl;
  std::cout<<"nPatGlobalCut: "<<nPatGlobalCut<<std::endl;
  std::cout<<"nPatPtCut: "<<nPatPtCut<<std::endl;
  std::cout<<"nPatEtaCut: "<<nPatEtaCut<<std::endl;
  std::cout<<"nPatHadronChargedIsoCut:"<<nPatHadronChargedIsoCut<<std::endl; 
  std::cout<<"nVectorSizeCut: "<<nVectorSizeCut<<std::endl;
  std::cout<<"nOppositeChargeCut: "<<nOppositeChargeCut<<std::endl;
  std::cout<<"nGenTotal: "<<nGenTotal<<std::endl;
  std::cout<<"nGenEtaCut: "<<nGenEtaCut<<std::endl;
  std::cout<<"nGenPtCUt: "<<nGenPtCut<<std::endl;
  std::cout<<"nGenMatchedMuonPairs: "<<nGenMatchedMuonPairs<<std::endl;
  std::cout<<"nRecMassCut: "<<nRecMassCut<<std::endl;
  std::cout<<"nRecProfileEntries: " <<nRecProfileEntries<<std::endl;
  std::cout<<"nResolutionHistogramEntries:" <<nResolutionHistogramEntries<<std::endl;
  std::cout<<"nResolutionTProfileEntries: " <<nResolutionTProfileEntries<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonExercise2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise2);
