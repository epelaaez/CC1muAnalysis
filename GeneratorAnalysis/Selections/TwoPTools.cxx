
#ifndef TwoPTools_cxx
#define TwoPTools_cxx

#include <TMath.h>
#include <TVector3.h>

#include "TwoPTools.h"

using namespace std;

TwoPTools::TwoPTools(TVector3 MuonVector, TVector3 LeadingProtonVector, TVector3 RecoilProtonVector) {
    // Swap lead and recoil protons to make sure lead has higher magnitude
    if (RecoilProtonVector.Mag() > LeadingProtonVector.Mag()) { 
        std::swap(LeadingProtonVector, RecoilProtonVector); 
    }

    // Momenta and theta variables
    fMuonMomentum = MuonVector.Mag();
    fMuonCosTheta = MuonVector.CosTheta();
    fLeadingProtonMomentum = LeadingProtonVector.Mag();
    fLeadingProtonCosTheta = LeadingProtonVector.CosTheta();
    fRecoilProtonMomentum = RecoilProtonVector.Mag();
    fRecoilProtonCosTheta = RecoilProtonVector.CosTheta();

    // Opening angle variables
    TVector3 TotalProtonVector = LeadingProtonVector + RecoilProtonVector;
    fCosOpeningAngleProtons = std::cos(LeadingProtonVector.Angle(RecoilProtonVector));
    fCosOpeningAngleMuonTotalProton = std::cos(MuonVector.Angle(TotalProtonVector));

    // Transverse vectors
    TVector3 MuonVectorTrans(MuonVector.X(), MuonVector.Y(), 0);
    TVector3 TotalProtonVectorTrans(TotalProtonVector.X(), TotalProtonVector.Y(), 0);
    TVector3 DeltaPVectorTrans = MuonVectorTrans + TotalProtonVectorTrans;

    // Transverse momentum variable
    fTransverseMomentum = DeltaPVectorTrans.Mag();

    // Transverse momentum angular orientation with respect to transverse muon
    fDeltaAlphaT = TMath::ACos(
        (-MuonVectorTrans).Dot(DeltaPVectorTrans) / 
        (MuonVectorTrans.Mag() * DeltaPVectorTrans.Mag())
    ) * (180. / TMath::Pi()) ;
    if (fDeltaAlphaT > 180.) { fDeltaAlphaT -= 180.; }
	if (fDeltaAlphaT < 0.) { fDeltaAlphaT += 180.; }
}

double TwoPTools::ReturnMuonMomentum() {
    return fMuonMomentum;
}

double TwoPTools::ReturnMuonCosTheta() {
    return fMuonCosTheta;
}

double TwoPTools::ReturnLeadingProtonMomentum() {
    return fLeadingProtonMomentum;
}

double TwoPTools::ReturnLeadingProtonCosTheta() {
    return fLeadingProtonCosTheta;
}

double TwoPTools::ReturnRecoilProtonMomentum() {
    return fRecoilProtonMomentum;
}

double TwoPTools::ReturnRecoilProtonCosTheta() {
    return fRecoilProtonCosTheta;
}

double TwoPTools::ReturnCosOpeningAngleProtons() {
    return fCosOpeningAngleProtons;
}

double TwoPTools::ReturnCosOpeningAngleMuonTotalProton() {
    return fCosOpeningAngleMuonTotalProton;
}

double TwoPTools::ReturnTransverseMomentum() {
    return fTransverseMomentum;
}

double TwoPTools::ReturnDeltaAlphaT() {
    return fDeltaAlphaT;
}

#endif
