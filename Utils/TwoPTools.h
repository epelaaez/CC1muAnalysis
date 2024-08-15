#ifndef TwoPTools_h
#define TwoPTools_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class TwoPTools {
    private: 
        double fMuonMomentum;
        double fMuonCosTheta;
        double fLeadingProtonMomentum;
        double fLeadingProtonCosTheta;
        double fRecoilProtonMomentum;
        double fRecoilProtonCosTheta;
        double fCosOpeningAngleProtons;
        double fCosOpeningAngleMuonTotalProton;
        double fTransverseMomentum;
        double fDeltaAlphaT;
        
 	//Rick's new added variables
 	double fInvariantMass;
	double fDoubleTransverseMissingMomentum;   //the double-transverse missing momentum/momentum imbalance
	double fCosOpeningAngleLProtonMuon; //LProton = Leading Proton
	double fCosOpeningAngleRProtonMuon; //RProton = Recoil Proton

    public:
        // Default constructor
        TwoPTools(TVector3 MuonVector, TVector3 LeadingProtonVector, TVector3 RecoilProtonVector);

        // Default destructor
        ~TwoPTools() = default;

        // Getter functions
        double ReturnMuonMomentum();
        double ReturnMuonCosTheta();
        double ReturnLeadingProtonMomentum();
        double ReturnLeadingProtonCosTheta();
        double ReturnRecoilProtonMomentum();
        double ReturnRecoilProtonCosTheta();
        double ReturnCosOpeningAngleProtons();
        double ReturnCosOpeningAngleMuonTotalProton();
        double ReturnTransverseMomentum();
        double ReturnDeltaAlphaT();
	
	//Additional Getter Functions
	double ReturnInvariantMass();
	double ReturnDoubleTransverseMissingMomentum();
	double ReturnCosOpeningAngleLProtonMuon();
	double ReturnCosOpeningAngleRProtonMuon();

};

#endif
