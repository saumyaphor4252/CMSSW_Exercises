// -*- C++ -*-
//
// Package:    cmsdas/ZprimeAnalyser
// Class:      ZprimeAnalyser
// 
/**\class ZprimeAnalyser ZprimeAnalyser.cc cmsdas/ZprimeAnalyser/plugins/ZprimeAnalyser.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Jason Lee
//         Created:  Thu, 13 Aug 2015 16:40:39 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std;

class ZprimeAnalyser : public edm::EDAnalyzer {
public:
  explicit ZprimeAnalyser(const edm::ParameterSet&);
  ~ZprimeAnalyser();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  edm::EDGetTokenT<edm::View<pat::Muon> >     muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > elecToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<pat::MET> >      metToken_;
  edm::EDGetTokenT<reco::VertexCollection >   vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> mcLabel_;

  TTree * ttree_;
  int b_nmuon, b_njet, b_step, b_channel;
  float b_MET;
  float b_lep1_pt, b_lep1_eta, b_lep1_phi;
  float b_lep2_pt, b_lep2_eta, b_lep2_phi;
  float b_reco_diMu_pt, b_reco_diMu_eta, b_reco_diMu_phi, b_reco_diMu_m;
  float b_gen_diMu_pt, b_gen_diMu_eta, b_gen_diMu_phi, b_gen_diMu_m;

  vector<int> *b_nSiTrackerLayers, *b_nPixelHits, *b_nMuonHits, *b_nMatches;
  vector<float> *b_relPtError, *b_dxy, *b_relTrkIso, *b_ptBestTrack;
  bool b_isMedium, b_isTight;
  
  bool runOnMC_;
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
ZprimeAnalyser::ZprimeAnalyser(const edm::ParameterSet& iConfig)
{
  muonToken_ = consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"));
  elecToken_ = consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("elecs"));
  jetToken_  = consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"));
  metToken_  = consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("mets"));     
  vtxToken_  = consumes<reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("vertices"));
  mcLabel_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("mcLabel"));
  
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("tree", "tree");
  ttree_->Branch("nSiTrackerLayers", &b_nSiTrackerLayers);
  ttree_->Branch("nPixelHits", &b_nPixelHits);
  ttree_->Branch("nMuonHits", &b_nMuonHits);
  ttree_->Branch("nMatches", &b_nMatches);
  ttree_->Branch("relPtError", &b_relPtError);
  ttree_->Branch("dxy", &b_dxy);
  ttree_->Branch("relTrkIso", &b_relTrkIso);
  ttree_->Branch("ptBestTrack", &b_ptBestTrack);

  ttree_->Branch("nmuon", &b_nmuon, "nmuon/I");
  ttree_->Branch("njet", &b_njet, "njet/I");

  ttree_->Branch("reco_diMu_pt", &b_reco_diMu_pt, "reco_diMu_pt/F");
  ttree_->Branch("reco_diMu_eta", &b_reco_diMu_eta, "reco_diMu_eta/F");
  ttree_->Branch("reco_diMu_phi", &b_reco_diMu_phi, "reco_diMu_phi/F");
  ttree_->Branch("reco_diMu_m", &b_reco_diMu_m, "reco_diMu_m/F");

  ttree_->Branch("gen_diMu_pt", &b_gen_diMu_pt, "gen_diMu_pt/F");
  ttree_->Branch("gen_diMu_eta", &b_gen_diMu_eta, "gen_diMu_eta/F");
  ttree_->Branch("gen_diMu_phi", &b_gen_diMu_phi, "gen_diMu_phi/F");
  ttree_->Branch("gen_diMu_m", &b_gen_diMu_m, "gen_diMu_m/F");
  
  ttree_->Branch("MET", &b_MET, "MET/F");
  ttree_->Branch("channel", &b_channel, "channel/I");
  ttree_->Branch("step", &b_step, "step/I");

  ttree_->Branch("isMedium", &b_isMedium, "isMedium/B");
  ttree_->Branch("isTight", &b_isTight, "isTight/B");

  ttree_->Branch("mu1_pt", &b_lep1_pt, "mu1_pt/F");
  ttree_->Branch("mu1_eta", &b_lep1_eta, "mu1_eta/F");
  ttree_->Branch("mu1_phi", &b_lep1_phi, "mu1_phi/F");

  ttree_->Branch("mu2_pt", &b_lep2_pt, "mu2_pt/F");
  ttree_->Branch("mu2_eta", &b_lep2_eta, "mu2_eta/F");
  ttree_->Branch("mu2_phi", &b_lep2_phi, "mu2_phi/F");


}
ZprimeAnalyser::~ZprimeAnalyser(){}

// ------------ method called for each event  ------------
void
ZprimeAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_nmuon = -1; b_njet = -1; b_step = 0; b_channel = -1;
  b_MET = -1;
  b_lep1_pt = -9; b_lep1_eta = -9; b_lep1_phi = -9;
  b_lep2_pt = -9; b_lep2_eta = -9; b_lep2_phi = -9;
  b_reco_diMu_pt = -9; b_reco_diMu_eta = -9; b_reco_diMu_phi = -9; b_reco_diMu_m = -9;
  b_gen_diMu_pt = -9; b_gen_diMu_eta = -9; b_gen_diMu_phi = -9; b_gen_diMu_m = -9;
  b_isMedium = 0; b_isTight = 0;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &pv = vertices->front();

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<edm::View<pat::Electron> > elecs;
  iEvent.getByToken(elecToken_, elecs);

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByToken(metToken_, mets);

  if (!iEvent.isRealData()){
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(mcLabel_,genParticles);
    for (const reco::GenParticle & g : *genParticles){
      if (fabs(g.pdgId()) == 32 && fabs(g.daughter(0)->pdgId()) == 13 ){
	b_gen_diMu_pt  = g.pt();
	b_gen_diMu_m   = g.mass();
	b_gen_diMu_eta = g.eta();
	b_gen_diMu_phi = g.phi();
      }
    }
  }

  b_nSiTrackerLayers = new vector<int>();
  b_nPixelHits = new vector<int>();
  b_nMuonHits = new vector<int>();
  b_nMatches = new vector<int>();
  b_relPtError = new vector<float>();
  b_dxy = new vector<float>();
  b_relTrkIso = new vector<float>();
  b_ptBestTrack = new vector<float>();

  b_nmuon = muons->size();
  
  vector<pat::Muon> selectedmuons; selectedmuons.clear();  
  for (auto mu : *muons) {    
    if (not mu.isGlobalMuon())
      continue;
    if (not mu.isTrackerMuon())
      continue;
    b_nSiTrackerLayers->push_back(mu.globalTrack()->hitPattern().trackerLayersWithMeasurement());
    b_nPixelHits->push_back(mu.globalTrack()->hitPattern().numberOfValidPixelHits());
    b_nMuonHits->push_back(mu.globalTrack()->hitPattern().numberOfMuonHits());
    b_nMatches->push_back(mu.numberOfChambers());

    float dxy = abs(mu.muonBestTrack()->dxy(pv.position()));
    if (dxy > 0.3) dxy = 0.3;
    b_dxy->push_back(dxy);
    
    float relTrkIso = mu.isolationR03().sumPt / mu.innerTrack()->pt();
    if (relTrkIso > 0.3) relTrkIso = 0.3;
    b_relTrkIso->push_back(relTrkIso);

    if (mu.pt() > 100.){
      reco::TrackRef cktTrack = (muon::tevOptimized(mu)).first;
      math::XYZTLorentzVectorD momentum = mu.p4() ;
      momentum *= cktTrack->pt()/momentum.pt();
      mu.setP4(momentum);// adding best TeV momentum
      float relPtError = cktTrack->ptError()/cktTrack->pt();
      if (relPtError > 0.5) relPtError = 0.5;
      b_relPtError->push_back(relPtError);
      if (cktTrack->ptError()/cktTrack->pt() > 0.3)
    	continue;
    }
    const reco::TrackRef& tunePTrack = mu.muonBestTrack();
    b_ptBestTrack->push_back(tunePTrack->pt());
    
    if (tunePTrack->pt() <= 45.)
      continue;
    if (abs(mu.muonBestTrack()->dxy(pv.position())) > 0.2)
      continue;
    if (abs(mu.dB()) > 0.2)
      continue;
    if (mu.isolationR03().sumPt / mu.innerTrack()->pt() > 0.10 )
      continue;
    if (mu.globalTrack()->hitPattern().trackerLayersWithMeasurement() <= 5)
      continue;
    if (mu.globalTrack()->hitPattern().numberOfValidPixelHits() == 0)
      continue;
    if (mu.globalTrack()->hitPattern().numberOfValidMuonHits() == 0)
      continue;
    if (mu.numberOfMatchedStations() < 2)
      continue;
    selectedmuons.push_back(mu);
  }
  vector<pat::Electron> selectedelecs; selectedelecs.clear();  
  for (auto el : *elecs) {
    //cout << "el.electronID() "<< el.electronID("egmGsfElectronIDs:heepElectronID-HEEPV60") << endl;
    if (!el.electronID("heepElectronID-HEEPV50-CSA14-25ns") ) continue;
    selectedelecs.push_back(el);
  }

  if (selectedmuons.size() == 2){
    if (selectedmuons[0].charge() * selectedmuons[1].charge() < 0){
      TLorentzVector mu1, mu2, reco_diMu;
      mu1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.105658);
      mu2.SetPtEtaPhiM(selectedmuons[1].pt(), selectedmuons[1].eta(), selectedmuons[1].phi(), 0.105658);
      reco_diMu = mu1 + mu2;
      b_reco_diMu_pt  = reco_diMu.Pt();
      b_reco_diMu_m   = reco_diMu.M();
      b_reco_diMu_eta = reco_diMu.Eta();
      b_reco_diMu_phi = reco_diMu.Phi();
    }
  }

  ttree_->Fill();
}
void ZprimeAnalyser::beginJob(){}
void ZprimeAnalyser::endJob(){}
void ZprimeAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZprimeAnalyser);
