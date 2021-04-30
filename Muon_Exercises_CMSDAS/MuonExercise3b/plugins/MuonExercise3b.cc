// -*- C++ -*-
//
// Package:    MuonAnalysis/MuonExercise3b
// Class:      MuonExercise3b
// 
/**\class MuonExercise3b MuonExercise3b.cc MuonAnalysis/MuonExercise3b/plugins/MuonExercise3b.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saumya Saumya
//         Created:  Sun, 14 Jun 2020 11:16:01 GMT
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

class MuonExercise3b : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonExercise3b(const edm::ParameterSet&);
      ~MuonExercise3b();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      enum MuonParentage {PROMPT,HF,LF,OTHER,N_MUPAR_TYPES};
      const char* enumNames[MuonParentage::N_MUPAR_TYPES]={"prompt","hf","lf","other"};


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      MuonParentage getParentType(const reco::GenParticle& prt);

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;

      bool isQCD;

      // Instantiate plots for each of the variables you want to study
      // Since we want to study these variables for all the 4 categories above
      // (prompt, heavy flavor, light flavor, other)
      // a good idea is to instantiate ARRAYS of histograms
      // Each array contains a histogram per category 
      TH1F* h_pt       [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_eta      [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_isglobal [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_istracker[MuonParentage::N_MUPAR_TYPES];
      TH1F* h_ispf     [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_nchi2    [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_nhits    [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_nstations[MuonParentage::N_MUPAR_TYPES];
      TH1F* h_npixels  [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_nlayers  [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_d0       [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_chiso    [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_nhiso    [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_emiso    [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_puiso    [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_combiso  [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_corriso  [MuonParentage::N_MUPAR_TYPES];

 
      TH1F* h_pt_loose [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_pt_medium[MuonParentage::N_MUPAR_TYPES];
      TH1F* h_pt_tight [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_pt_soft  [MuonParentage::N_MUPAR_TYPES];
      TH1F* h_pt_highpt[MuonParentage::N_MUPAR_TYPES];

      TH1F* h_nvtx               [MuonParentage::N_MUPAR_TYPES]; 
      TH1F* h_nvtx_tight         [MuonParentage::N_MUPAR_TYPES]; 
      TH1F* h_nvtx_tight_iso     [MuonParentage::N_MUPAR_TYPES]; 
      TH1F* h_nvtx_tight_iso_corr[MuonParentage::N_MUPAR_TYPES]; 

      TH1F* h_nvtx_tight_iso_eff     [MuonParentage::N_MUPAR_TYPES]; 
      TH1F* h_nvtx_tight_iso_corr_eff[MuonParentage::N_MUPAR_TYPES]; 
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
MuonExercise3b::MuonExercise3b(const edm::ParameterSet& iConfig):
 muonsToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonTag"))),
 genParticlesToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("GenPartTag"))),
 verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VertexTag")))
{
   //now do what ever initialization is needed

   // Flag to identify the QCD sample (the reason will be clear later)
   // In your configuration file, set it to 'True' only when running on the QCD samples 
   isQCD = iConfig.getUntrackedParameter<bool>("IsQCD");

   //TFileService to save histograms to an output file
   usesResource("TFileService");
   edm::Service<TFileService> fs;

  for (size_t idx=0; idx<MuonParentage::N_MUPAR_TYPES; idx++) 
  {
    h_pt       [idx] = fs->make<TH1F>(Form("pt_%s"       , enumNames[idx]), ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_eta      [idx] = fs->make<TH1F>(Form("eta_%s"      , enumNames[idx]), ";Muon #eta;Events", 25, -2.5 ,   2.5 );
    h_isglobal [idx] = fs->make<TH1F>(Form("isglobal_%s" , enumNames[idx]), ";isGlobalMuon();Events",  2, -0.5 ,   1.5 );
    h_istracker[idx] = fs->make<TH1F>(Form("istracker_%s", enumNames[idx]), ";isTrackerMuon();Events",  2, -0.5 ,   1.5 );
    h_nchi2    [idx] = fs->make<TH1F>(Form("nchi2_%s"    , enumNames[idx]), ";Glb-track #chi^{2};Events", 55,  0.  ,  11.  );
    h_nhits    [idx] = fs->make<TH1F>(Form("nhits_%s"    , enumNames[idx]), ";Muon hits;Events", 40, -0.5 ,  39.5 );
    h_nstations[idx] = fs->make<TH1F>(Form("nstations_%s", enumNames[idx]), ";Muon stations;Events",  6, -0.5 ,   5.5 );
    h_npixels  [idx] = fs->make<TH1F>(Form("npixels_%s"  , enumNames[idx]), ";Pixel hits;Events",  6, -0.5 ,   5.5 );
    h_nlayers  [idx] = fs->make<TH1F>(Form("nlayers_%s"  , enumNames[idx]), ";Tracker layers;Events", 25, -0.5 ,  24.5 );
    h_d0       [idx] = fs->make<TH1F>(Form("d0_%s"       , enumNames[idx]), ";Track d_{0} [cm];Events", 40, -0.10,   0.11);
    h_chiso    [idx] = fs->make<TH1F>(Form("chiso_%s"    , enumNames[idx]), ";Charged hadron isolation;Events", 20,  0.  ,  10.  );
    h_nhiso    [idx] = fs->make<TH1F>(Form("nhiso_%s"    , enumNames[idx]), ";Neutral hadron isolation;Events", 20,  0.  ,  10.  );
    h_emiso    [idx] = fs->make<TH1F>(Form("emiso_%s"    , enumNames[idx]), ";Photon isolation ;Events", 20,  0.  ,  10.  );
    h_puiso    [idx] = fs->make<TH1F>(Form("puiso_%s"    , enumNames[idx]), ";PU charged hadron isolation;Events", 20,  0.  ,  10.  );
    h_combiso  [idx] = fs->make<TH1F>(Form("combiso_%s"  , enumNames[idx]), ";Combined relative isolation;Events", 30,  0.  ,   0.3 );
    h_corriso  [idx] = fs->make<TH1F>(Form("corriso_%s"  , enumNames[idx]), ";#Delta#beta combined relative isolation;Events", 30,  0.  ,   0.3 );

    h_pt_loose [idx] = fs->make<TH1F>(Form("pt_loose_%s" , enumNames[idx]), ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_medium[idx] = fs->make<TH1F>(Form("pt_medium_%s", enumNames[idx]), ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_tight [idx] = fs->make<TH1F>(Form("pt_tight_%s" , enumNames[idx]), ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_soft  [idx] = fs->make<TH1F>(Form("pt_soft_%s"  , enumNames[idx]), ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
    h_pt_highpt[idx] = fs->make<TH1F>(Form("pt_highpt_%s", enumNames[idx]), ";Muon p_{T} [GeV];Events", 20,  0.  , 100.  );
	
    h_nvtx               [idx] = fs->make<TH1F>(Form("nvtx_%s"               , enumNames[idx]), ";Reconstructed vertices;Events", 60,  0.5, 60.5);
    h_nvtx_tight         [idx] = fs->make<TH1F>(Form("nvtx_tight_%s"         , enumNames[idx]), ";Reconstructed vertices;Events", 60,  0.5, 60.5);
    h_nvtx_tight_iso     [idx] = fs->make<TH1F>(Form("nvtx_tight_iso_%s"     , enumNames[idx]), ";Reconstructed vertices;Events", 60,  0.5, 60.5);
    h_nvtx_tight_iso_corr[idx] = fs->make<TH1F>(Form("nvtx_tight_iso_corr_%s", enumNames[idx]), ";Reconstructed vertices;Events", 60,  0.5, 60.5);

    h_nvtx_tight_iso_eff     [idx] = fs->make<TH1F>(Form("nvtx_tight_iso_eff_%s"     , enumNames[idx]), ";Reconstructed vertices;Efficiency", 60,  0.5, 60.5);
    h_nvtx_tight_iso_corr_eff[idx] = fs->make<TH1F>(Form("nvtx_tight_iso_corr_eff_%s", enumNames[idx]), ";Reconstructed vertices;Efficiency", 60,  0.5, 60.5);
  }
}


MuonExercise3b::~MuonExercise3b()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonExercise3b::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  // Retrieve the GenParticle collection from the event 
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);
  std::vector<reco::GenParticle> genColl = (*genParticles);
  // Check that the collection is valid 
  if(!genParticles.isValid())
  {
    throw cms::Exception("GenParticle collection not valid!");
  }

  // Retrieve the pat::Muon Collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonsToken_, muons);
  if(!muons.isValid()) 
  {
    throw cms::Exception("Muon collection not valid!");
  }

  // Retrieve the vertex Collection 
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(verticesToken_, vertices);
  if(!vertices.isValid()) 
  {
    throw cms::Exception("Vertex collection not valid!");
  }

  // Check that we have at least one good vertex
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
  int nvtx=vertices->size();

  edm::View<pat::Muon>::const_iterator muend = muons->end();
  for (edm::View<pat::Muon>::const_iterator it=muons->begin(); it!=muend; ++it) 
  {
    // Muon should have a silicon track (i.e. global- or tracker-muon) 
    if (!it->globalTrack().isNonnull()) continue;
    if (!it->innerTrack().isNonnull()) continue;
    if (fabs(it->innerTrack()->dz(firstGoodVertex->position()))>0.2) continue;   
 
    // Muon sholud be within the tracker volume and have pt > 20 GeV
    if (std::abs(it->eta())>2.5) continue;
    if (it->pt()<20) continue;

    
    //Check the origin of the muon: prompt, HF, LF, other?
    MuonParentage parentage = MuonParentage::OTHER;

    // For this, we need to find the gen particle associated with our reconstructed muon
    // Since we are using pat::Muons, the gen-matching is already done for us!
    // No need to loop over the GenParticle collection and perform a geometrical matching 
    const reco::GenParticle* gp = it->genParticle();
   
    // Check if the pat::Muon has a GenParticle associated
    // In what cases there is no gen-matching? 
    if(gp!=0)
    {
        // Function determines the muon origin for you
        parentage = getParentType(*gp);
    }
    else
    {
        if(!isQCD) continue;
        // If there is no genParticle() and it's a Drell-Yan or top-top sample, stop here!
        // Proceed with the classification only when running on the QCD sample! 
        // In all the other samples, muons from light flavor decays are NOT saved
        // in the GenParticle collection. Therefore, light-flavor decays would be
        // classified as "other", which is not correct.
        // QCD samples, on the other hand, have gen-particles also for light-flavor decays
    }

    // Fill plots
    h_pt[parentage]->Fill(std::min(it->pt(),99.9));
    h_eta[parentage]->Fill(std::max(std::min(it->eta(), 2.49), -2.49));
    h_isglobal[parentage]->Fill(it->isGlobalMuon());
    h_istracker[parentage]->Fill(it->isTrackerMuon());
    h_nchi2[parentage]->Fill(std::min(it->globalTrack()->normalizedChi2(), 10.99));
    h_nhits[parentage]->Fill(std::min(it->globalTrack()->hitPattern().numberOfValidMuonHits(), 39) );
    h_nstations[parentage]->Fill(std::min(it->numberOfMatchedStations(), 5));
    h_npixels[parentage]->Fill(std::min(it->innerTrack()->hitPattern().numberOfValidPixelHits(), 5));
    h_nlayers[parentage]->Fill(std::min(it->innerTrack()->hitPattern().trackerLayersWithMeasurement(), 14));
    h_d0[parentage]->Fill(std::min(std::max(it->muonBestTrack()->dxy(firstGoodVertex->position()), -0.299), 0.299));
   
    // Wanring: before using the isolation information stored in the Muon object, check if it's valid 
    if (it->isIsolationValid()) 
    {
        // Plot all the isolation variables declared above, then build the "PF combined relative isolation", with and without Delta-beta correction 
        const reco::MuonPFIsolation &pfR04 = it->pfIsolationR04();
        h_chiso[parentage]->Fill(std::min(pfR04.sumChargedHadronPt, (float)9.9));
        h_nhiso[parentage]->Fill(std::min(pfR04.sumNeutralHadronEt, (float)9.9));
        h_emiso[parentage]->Fill(std::min(pfR04.sumPhotonEt, (float)9.9));
        h_puiso[parentage]->Fill(std::min(pfR04.sumPUPt, (float)9.9));
        double combiso = pfR04.sumChargedHadronPt + std::max(0., double(pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt));
        double corriso = pfR04.sumChargedHadronPt + std::max(0., double(pfR04.sumNeutralHadronEt+pfR04.sumPhotonEt-0.5*pfR04.sumPUPt));
        h_combiso[parentage]->Fill(std::min(combiso/it->pt(), 0.299));
        h_corriso[parentage]->Fill(std::min(corriso/it->pt(), 0.299));


        // Pileup dependence: plot the number of reconstructed vertices
        // 1) before applying an isolation cut,
        // 2) after  applying an isolation cut without Delta-beta correction, and 
        // 3) after  applying an isolation cut with Delta-beta correction 
        h_nvtx[parentage]->Fill(std::min(nvtx,60));
        if(it->isTightMuon(*firstGoodVertex)) 
        {
            h_nvtx_tight[parentage]->Fill(std::min(nvtx,60));
            if(combiso/it->pt()<0.15) 
	    {
              h_nvtx_tight_iso[parentage]->Fill(std::min(nvtx,60));
            }
            if(corriso/it->pt()<0.15) 
            {
              h_nvtx_tight_iso_corr[parentage]->Fill(std::min(nvtx,60));
            }
        }
    }

    // Fill some histogram after requiring that the muon passes
    // each of the standard Muon POG ID's (Loose, Medium, Tight, Soft, HighPt)
    if(it->isLooseMuon()) 
    {
      h_pt_loose[parentage]->Fill(std::min(it->pt(),99.9));
    }
    if(it->isMediumMuon()) 
    {
      h_pt_medium[parentage]->Fill(std::min(it->pt(),99.9));
    }
    if(it->isTightMuon(*firstGoodVertex)) 
    {
      h_pt_tight[parentage]->Fill(std::min(it->pt(),99.9));
    }
    if(it->isSoftMuon(*firstGoodVertex)) 
    {
      h_pt_soft[parentage]->Fill(std::min(it->pt(),99.9));
    }
    if(it->isHighPtMuon(*firstGoodVertex)) 
    {
      h_pt_highpt[parentage]->Fill(std::min(it->pt(),99.9));
    }
 } // end for (edm::View<pat::Muon>::const_iterator it = muons->begin(); it != muon_end; it++)

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
MuonExercise3b::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonExercise3b::endJob() 
{
    // Normalize efficiency plots:
    // divide the number of reconstructed vertices after the isolation and before the isolation cut!
    // Then check the efficiency dependence as a function of the number of vertices 
    for (size_t idx=0; idx<MuonParentage::N_MUPAR_TYPES; idx++) 
    {  
        //  - iso/tight    
        h_nvtx_tight_iso_eff[idx]->Divide(h_nvtx_tight_iso[idx], h_nvtx_tight[idx], 1, 1, "B");

        //  - corrected iso/tight 
        h_nvtx_tight_iso_corr_eff[idx]->Divide(h_nvtx_tight_iso_corr[idx], h_nvtx_tight[idx], 1, 1, "B");
    }
}

// ------------ method to determine whether a nonprompt muon is from HF  ------------
MuonExercise3b::MuonParentage MuonExercise3b::getParentType (const reco::GenParticle& prt)
{
  // Prompt: mother ID = 13, 15 (with prompt tau), 23, 24, 25 
  if(prt.isPromptFinalState() || prt.isDirectPromptTauDecayProductFinalState()) 
  {
    return MuonParentage::PROMPT;
  }

  auto mom = &prt;
  bool sameid = true;
  while(mom->numberOfMothers()>0 && sameid) 
  {
    for(size_t im=0; im<mom->numberOfMothers(); ++im) 
    {
       mom = dynamic_cast<const reco::GenParticle*>(mom->mother(im));
       if(mom->pdgId()!=prt.pdgId() && std::abs(mom->pdgId())!=15) 
       {  // if it's a tau, it's not a prompt tau -> check tau's mother
          sameid = false;
          break;
       }
    }
  }

  // Classification
  unsigned int pdgId = abs(mom->pdgId());

  // - Prompt -- send warning: why didn't it pass isPromptFinalState() or isDirectPromptTauDecayProductFinalState()??
  if (pdgId == 13 || pdgId == 15 || pdgId == 23 || pdgId == 24 || pdgId == 25) 
  {
    std::cout << "[MuonExercise3::getParentType] WARNING: mother's PDG ID is "
              << pdgId << ", yet daughter did not pass isPromptFinalState() nor "
              << "isDirectPromptTauDecayProductFinalState()" << std::endl;

   return MuonParentage::PROMPT;
  }

  // - From heavy-flavor hadron decay
  bool ishf = false;
  int thirdDigit = (pdgId/100) % 10;
  if ( thirdDigit == 5 || thirdDigit == 4 ) ishf = true ; // should catch all B and C mesons
  int fourthDigit = (pdgId/1000) % 10;
  if ( fourthDigit == 5 || fourthDigit == 4 ) ishf = true ; // should catch all B and C baryons
  if (ishf) return MuonParentage::HF;

  // - From light-flavor hadron decay
  return MuonParentage::LF;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonExercise3b::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonExercise3b);
