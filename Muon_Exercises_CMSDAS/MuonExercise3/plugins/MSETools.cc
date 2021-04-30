#include "MSETools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/VectorUtil.h"

const char* mse::enumNames[4]={"prompt","hf","lf","other"};

bool mse::isBHadron(int pdgId)
{
    bool retval=false;
    int thirdDigit=(abs(pdgId)/100)%10;
    if(thirdDigit==5) retval=true; // should catch all B mesons
    int fourthDigit=(abs(pdgId)/1000)%10;
    if(fourthDigit==5) retval=true; //should catch all B mesons
    return retval;
}

bool mse::decaysToB(const reco::GenParticle& gp)
{
    for(size_t di=0;di<gp.numberOfDaughters();di++)
    {
        const reco::Candidate* dau=gp.daughter(di);
        if(isBHadron(dau->pdgId())) return true; 
    }
    return false;
}

bool mse::isCHadron(int pdgId)
{
    bool retval=false;
    int thirdDigit=(abs(pdgId)/100)%10;
    if(thirdDigit==4) retval=true; // should catch all C mesons
    int fourthDigit=(abs(pdgId)/1000)%10;
    if(fourthDigit==4) retval=true; //should catch all C mesons
    if(mse::isBHadron(pdgId)) retval=false;
    return retval;
}

bool mse::decaysToC(const reco::GenParticle& gp)    
{
    for(size_t di=0;di<gp.numberOfDaughters();di++)
    {
        const reco::Candidate* dau=gp.daughter(di);
        if(isCHadron(dau->pdgId())) return true; 
    }
    return false;
}

const pat::PackedGenParticle mse::getMatchedGenParticle(const pat::Muon& muon,const std::vector<pat::PackedGenParticle>& gpcol)
{
    LorentzVector mup4(muon.p4());
    double minDR = 999;
    std::vector<pat::PackedGenParticle>::const_iterator gen_end=gpcol.end();
    pat::PackedGenParticle matched_gen;
    for(std::vector<pat::PackedGenParticle>::const_iterator it=gpcol.begin();it!=gen_end;it++)
    {
        LorentzVector gp4(it->p4());
        double dr=ROOT::Math::VectorUtil::DeltaR(mup4,gp4);
        if(dr<0.1 && abs(it->pdgId())==13) return (*it);
        if(dr<minDR)
        {
         minDR=dr;
         matched_gen=*it;
        }
    }
    if(minDR<0.1) return matched_gen;
    return pat::PackedGenParticle();
}

const reco::GenParticle mse::getMotherPacked(const pat::PackedGenParticle& pgp)
{
    if(pgp.numberOfMothers()>0)
    {  const reco::GenParticle* firstMother=(const reco::GenParticle*)pgp.mother(0);
       if(firstMother!=0)
       { if(firstMother->pdgId()!=pgp.pdgId()) return (*firstMother);
         const reco::GenParticle newMother= mse::getMother(*firstMother);
         return newMother;
       }
       else return reco::GenParticle();
    }
    return reco::GenParticle();
}

const reco::GenParticle mse::getMother(const reco::GenParticle& gp)
{
    const reco::GenParticle* mom = &gp;
    while(mom->numberOfMothers()>0)
    {
        for(unsigned int idx=0;idx<mom->numberOfMothers();idx++)
        {
            mom=dynamic_cast<const reco::GenParticle*>(mom->mother(idx));
            if(mom->pdgId()!=gp.pdgId())
            return (*mom);
        }
    }
    return (*mom);
}  

mse::MuonParentage mse::getParentType(const reco::GenParticle& gp)
{
    unsigned int pdgId=abs(gp.pdgId());
    if(pdgId==15 || pdgId==23 || pdgId==24 || pdgId==25) return mse::MuonParentage::PROMPT;
    if(mse::isBHadron(pdgId) || mse::isCHadron(pdgId)) return mse::MuonParentage::HF;
    return mse::MuonParentage::LF;
}


bool mse::isGoodVertex(const reco::Vertex& vtx)
{
    if(vtx.isFake()) return false;
    if(vtx.ndof()<4) return false;
    if(vtx.position().Rho()>2.0) return false;
    if(fabs(vtx.position().Z())>24.0) return false;
    return true;
} 
