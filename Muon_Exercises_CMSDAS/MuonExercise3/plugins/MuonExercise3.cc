// -*- C++ -*-
//
// Package:    MuonAnalysis/MuonExercise3
// Class:      MuonExercise3
// 
/**\class MuonExercise3 MuonExercise3.cc MuonAnalysis/MuonExercise3/plugins/MuonExercise3.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Sat, 13 Jun 2020 12:49:02 GMT
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/VectorUtil.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "MSETools.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuonExercise3 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonExercise3(const edm::ParameterSet&);
      ~MuonExercise3();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genMuonsToken_;
      edm::EDGetTokenT<reco::VertexCollection> verticesToken_;

      TH1F* h_pt[4];
      TH1F* h_eta[4];
      TH1F* h_nchi2[4];
      TH1F* h_nhits[4];
      TH1F* h_nstations[4];
      TH1F* h_npixels[4];
      TH1F* h_nlayers[4];
      TH1F* h_d0[4];
      TH1F* h_chiso[4];
      TH1F* h_nhiso[4];
      TH1F* h_emiso[4];
      TH1F* h_puiso[4];
      TH1F* h_coriso[4];

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
MuonExercise3::MuonExercise3(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::InputTag muonlabel("slimmedMuons");
   edm::InputTag genMuonlabel("packedGenParticles");
   edm::InputTag vertexlabel("offlineSlimmedPrimaryVertices");

   
   muonsToken_ = consumes<edm::View<pat::Muon>>(muonlabel);
   genMuonsToken_=consumes<pat::PackedGenParticleCollection>(genMuonlabel);
   verticesToken_=consumes<reco::VertexCollection>(vertexlabel);
 
   edm::Service<TFileService>fs;
   for(unsigned int idx=0;idx<mse::MuonParentage::NMU_PAR_TYPES;idx++)
   {
      h_pt[idx] =fs->make<TH1F>(Form("h_pt_%s",mse::enumNames[idx]),Form("h_pt_%s",mse::enumNames[idx]),20,0,100);
      h_eta[idx] =fs->make<TH1F>(Form("h_eta_%s",mse::enumNames[idx]),Form("h_eta_%s",mse::enumNames[idx]),25,-2.5,2.5);
      h_nchi2[idx] =fs->make<TH1F>(Form("h_nchi2_%s",mse::enumNames[idx]),Form("h_nchi2_%s",mse::enumNames[idx]),55,0,11);
      h_nhits[idx] =fs->make<TH1F>(Form("h_nhits_%s",mse::enumNames[idx]),Form("h_nhits_%s",mse::enumNames[idx]),40,-0.5,39.5);
      h_nstations[idx] =fs->make<TH1F>(Form("h_nstations_%s",mse::enumNames[idx]),Form("h_nstations_%s",mse::enumNames[idx]),6,-0.5,5.5);
      h_npixels[idx] =fs->make<TH1F>(Form("h_npixels_%s",mse::enumNames[idx]),Form("h_npixels_%s",mse::enumNames[idx]),6,-0.5,5.5);
      h_nlayers[idx] =fs->make<TH1F>(Form("h_nlayers_%s",mse::enumNames[idx]),Form("h_nlayers_%s",mse::enumNames[idx]),25,-0.5,5.5);
      h_d0[idx] =fs->make<TH1F>(Form("h_d0_%s",mse::enumNames[idx]),Form("h_d0_%s",mse::enumNames[idx]),40,-0.1,0.11);
      h_chiso[idx] =fs->make<TH1F>(Form("h_chiso_%s",mse::enumNames[idx]),Form("h_chiso_%s",mse::enumNames[idx]),20,0,10);
      h_nhiso[idx] =fs->make<TH1F>(Form("h_nhiso_%s",mse::enumNames[idx]),Form("h_nhiso_%s",mse::enumNames[idx]),20,0,10);
      h_emiso[idx] =fs->make<TH1F>(Form("h_emiso_%s",mse::enumNames[idx]),Form("h_emiso_%s",mse::enumNames[idx]),20,0,10);
      h_puiso[idx] =fs->make<TH1F>(Form("h_puiso_%s",mse::enumNames[idx]),Form("h_puiso_%s",mse::enumNames[idx]),20,0,10);
      h_coriso[idx] =fs->make<TH1F>(Form("h_coriso_%s",mse::enumNames[idx]),Form("h_coriso_%s",mse::enumNames[idx]),30,0,0.3);
   }   
}


MuonExercise3::~MuonExercise3()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise3::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<edm::View<pat::Muon>> muons;
   iEvent.getByToken(muonsToken_,muons);

   edm::Handle<pat::PackedGenParticleCollection> genMuons;
   iEvent.getByToken(genMuonsToken_,genMuons);
   std::vector<pat::PackedGenParticle> gpcol=(*genMuons);

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(verticesToken_,vertices);
   reco::VertexCollection::const_iterator firstGoodVertex=vertices->end();
   for(reco::VertexCollection::const_iterator it=vertices->begin();it!=firstGoodVertex;it++)
   {
      mse::isGoodVertex(*it);
      firstGoodVertex=it;
      break;
   }

   //require a good vertex
   if(firstGoodVertex==vertices->end()) return;
   
   //
   // Loop over Muons
   //
   if( (!(muons.isValid())) || (!(genMuons.isValid())) ) return;
   edm::View<pat::Muon>::const_iterator muon_end=muons->end();
   for(edm::View<pat::Muon>::const_iterator it=muons->begin();it!=muon_end;it++)
   {
      //globalmuon
      if(!(it->globalTrack().isNonnull())) continue;
      
      //muon should have a silicon trackand be matched to the PV(dz<0.2)
      if(!(it->innerTrack().isNonnull())) continue;
      if(fabs(it->innerTrack()->dz(firstGoodVertex->position()))>0.2) continue;     
 
      //muon should be within tracker volume with pt>20GeV
      if(fabs(it->eta())>2.5) continue;
      if(it->pt()<20) continue;
      
      mse::MuonParentage parentage=mse::MuonParentage::NOT_A_MUON;
      int momid=0;
      const pat::PackedGenParticle matchedPackedGenParticle=mse::getMatchedGenParticle(*it,gpcol);
      if(abs(matchedPackedGenParticle.pdgId())==13)
      {
        const reco::GenParticle momgp=mse::getMotherPacked(matchedPackedGenParticle);
        if(momgp.pdgId()!=0)
        {
          momid=momgp.pdgId();
          parentage=mse::getParentType(momgp);
        }
      }

      h_pt[parentage]->Fill(std::min(it->pt(),99.9));
      h_eta[parentage]->Fill(std::max(std::min(it->eta(),2.49),-2.49));
      h_nchi2[parentage]->Fill(std::min(it->globalTrack()->normalizedChi2(),10.99));
      h_nhits[parentage]->Fill(std::min(it->globalTrack()->hitPattern().numberOfValidMuonHits(),39));
      h_nstations[parentage]->Fill(std::min(it->numberOfMatchedStations(),5));
      h_npixels[parentage]->Fill(std::min(it->innerTrack()->hitPattern().numberOfValidPixelHits(),5));
      h_nlayers[parentage]->Fill(std::min(it->innerTrack()->hitPattern().trackerLayersWithMeasurement(),14));
      h_d0[parentage]->Fill(std::min(std::max(it->muonBestTrack()->dxy(firstGoodVertex->position()),-0.299),-0.299));

      if(it->isIsolationValid())
      {  reco::MuonPFIsolation pfR04=it->pfIsolationR04();
         h_chiso[parentage]->Fill(std::min(pfR04.sumChargedHadronPt,(float)9.9));
         h_nhiso[parentage]->Fill(std::min(pfR04.sumNeutralHadronEt,(float)9.9));
         h_emiso[parentage]->Fill(std::min(pfR04.sumPhotonEt,(float)9.9));
         h_puiso[parentage]->Fill(std::min(pfR04.sumPUPt,(float)9.9));
         double coriso = pfR04.sumChargedHadronPt+std::max(0.0,pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt);
         h_coriso[parentage]->Fill(std::min(coriso/it->pt(),0.299));
      }
    if(parentage<mse::MuonParentage::LF)
    printf("pt,eta,pdgId,momid,parentage: %4.2f,%4.2f, %d,%d,%d \n",it->pt(),it->eta(),it->pdgId(),momid,parentage);
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
MuonExercise3::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise3::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonExercise3::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise3);
