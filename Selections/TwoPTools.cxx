
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
    fLeadingProtonMomentum = LeadingProtonVector.Mag();
    fLeadingProtonCosTheta = LeadingProtonVector.CosTheta();
    fRecoilProtonMomentum = RecoilProtonVector.Mag();
    fRecoilProtonCosTheta = RecoilProtonVector.CosTheta();

    // Opening angle variables
    // fCosOpeningAngleProtons = (LeadingProtonVector.Dot(RecoilProtonVector)) / (fLeadingProtonMomentum * fRecoilProtonMomentum);
    TVector3 TotalProtonVector = LeadingProtonVector + RecoilProtonVector;
    double fTotalProtonMomentum = TotalProtonVector.Mag();
    // fCosOpeningAngleMuonTotalProton = (TotalProtonVector.Dot(MuonVector)) / (fMuonMomentum * fTotalProtonMomentum);

    fCosOpeningAngleProtons = std::cos(LeadingProtonVector.Angle(RecoilProtonVector));
    fCosOpeningAngleMuonTotalProton = std::cos(MuonVector.Angle(TotalProtonVector));

    // Transverse momentum variable
    fTransverseMomentum = (TotalProtonVector + MuonVector).Perp();
}

double TwoPTools::ReturnMuonMomentum() {
    return fMuonMomentum;
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

#endif
