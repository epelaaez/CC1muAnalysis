#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"
#include "./Selections/TwoPTools.h"
#include "./Utils/Tools.h"
#include "Constants.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

using namespace std;
using namespace Constants;

// Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h);

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

    //----------------------------------------//	

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
    double A = 40.; // so that we can have xsecs per nucleus

    int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
    std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

    //----------------------------------------//	

    Tools tools;

    //----------------------------------------//	

    // Output file
    TString Directory = "/pnfs/sbnd/persistent/users/epelaez/HighSamples/";
    TString FileNameAndPath = Directory+"FlatTree/FlatTreeAnalyzerOutput_"+fOutputFile+".root";
    TFile* file = new TFile(FileNameAndPath,"recreate");

    std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
    std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
    
    //----------------------------------------//

    // Plot declaration

    // Single differential
    TH1D* TrueMuonCosThetaPlot[NInte];
    TH1D* TrueLeadingProtonCosThetaPlot[NInte];
    TH1D* TrueRecoilProtonCosThetaPlot[NInte];
    TH1D* TrueLeadingProtonMomentumPlot[NInte];
    TH1D* TrueRecoilProtonMomentumPlot[NInte];
    TH1D* TrueMuonMomentumPlot[NInte];
    TH1D* TrueCosOpeningAngleProtonsPlot[NInte];
    TH1D* TrueCosOpeningAngleMuonTotalProtonPlot[NInte];
    TH1D* TrueTransverseMomentumPlot[NInte];
    TH1D* TrueDeltaAlphaTPlot[NInte];

    TH1D* TrueNoFSILeadingProtonCosThetaPlot[NInte];
    TH1D* TrueNoFSIRecoilProtonCosThetaPlot[NInte];
    TH1D* TrueNoFSILeadingProtonMomentumPlot[NInte];
    TH1D* TrueNoFSIRecoilProtonMomentumPlot[NInte];
    TH1D* TrueNoFSIMuonMomentumPlot[NInte];
    TH1D* TrueNoFSICosOpeningAngleProtonsPlot[NInte];
    TH1D* TrueNoFSICosOpeningAngleMuonTotalProtonPlot[NInte];
    TH1D* TrueNoFSITransverseMomentumPlot[NInte];
    TH1D* TrueNoFSIDeltaAlphaTPlot[NInte];

    // Double differential
    TH1D* TrueSerialTransverseMomentum_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialDeltaAlphaT_InMuonCosThetaPlot[NInte];

    // Loop over the interaction processes

    for (int inte = 0; inte < NInte; inte++) {

    //--------------------------------------------------//

    // Final state
    TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#mu})",10,-1.,1.);
    TrueLeadingProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueLeadingProtonCosThetaPlot",";cos(#theta_{#vec{p}_{L}})",10,-1.,1.);
    TrueRecoilProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueRecoilProtonCosThetaPlot",";cos(#theta_{#vec{p}_{R}})",10,-1.,1.);
    TrueLeadingProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueLeadingProtonMomentumPlot",";|#vec{p}_{L}|",10,0.3,1);
    TrueRecoilProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueRecoilProtonMomentumPlot",";|#vec{p}_{R}|",10,0.3,1);
    TrueMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonMomentumPlot",";|#vec{p}_{#mu}|",10,0.1,1.2);
    TrueCosOpeningAngleProtonsPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleProtonsPlot",";cos(#theta_{#vec{p}_{L},#vec{p}_{R}})",10,-1.,1.);
    TrueCosOpeningAngleMuonTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleMuonTotalProtonPlot",";cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})",10,-1.,1.);
    TrueTransverseMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueTransverseMomentumPlot",";#delta P_{T}",10,0.,1.);
    TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",";#delta #alpha_{T}",10,0.,180.);

    // Before final state interactions
    TrueNoFSILeadingProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSILeadingProtonCosThetaPlot",";cos(#theta_{#vec{p}_{L}})",10,-1.,1.);
    TrueNoFSIRecoilProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIRecoilProtonCosThetaPlot",";cos(#theta_{#vec{p}_{R}})",10,-1.,1.);
    TrueNoFSILeadingProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSILeadingProtonMomentumPlot",";|#vec{p}_{L}|",10,0.3,1);
    TrueNoFSIRecoilProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIRecoilProtonMomentumPlot",";|#vec{p}_{R}|",10,0.3,1);
    TrueNoFSIMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIMuonMomentumPlot",";|#vec{p}_{#mu}|",10,0.1,1.2);
    TrueNoFSICosOpeningAngleProtonsPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSICosOpeningAngleProtonsPlot",";cos(#theta_{#vec{p}_{L},#vec{p}_{R}})",10,-1.,1.);
    TrueNoFSICosOpeningAngleMuonTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSICosOpeningAngleMuonTotalProtonPlot",";cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})",10,-1.,1.);
    TrueNoFSITransverseMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSITransverseMomentumPlot",";#delta P_{T}",10,0.,1.);
    TrueNoFSIDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIDeltaAlphaTPlot",";#delta #alpha_{T}",10,0.,180.);

    // Double differential final state
    TrueSerialTransverseMomentum_InMuonCosThetaPlot[inte] = new TH1D(
        InteractionLabels[inte]+"TrueSerialTransverseMomentum_InMuonCosThetaPlot",
        LabelXAxisTwoDTransverseMomentumInMuonCosTheta,
        tools.Return2DNBins(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices),
        &tools.Return2DBinIndices(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices)[0]
    );
    TrueSerialDeltaAlphaT_InMuonCosThetaPlot[inte] = new TH1D(
        InteractionLabels[inte]+"TrueSerialDeltaAlphaT_InMuonCosThetaPlot",
        LabelXAxisTwoDDeltaAlphaTInMuonCosTheta,
        tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),
        &tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]
    );

    //--------------------------------------------------//

    } // End of the loop over the interaction processes							

    //----------------------------------------//

    // Counters

    int CounterEventsPassedSelection = 0;
    int CounterQEEventsPassedSelection = 0;
    int CounterMECEventsPassedSelection = 0;
    int CounterRESEventsPassedSelection = 0;
    int CounterDISEventsPassedSelection = 0;
    int CounterCOHEventsPassedSelection = 0;

    // Counters for pre-FSI RES events
    int RESMode[24];
    for (int i = 0; i < 24; i++) {
        RESMode[i] = 0;
    }
    int NoFSICounterRESEventsPassedSelection = 0;

    //----------------------------------------//
    
    // Loop over the events

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        //----------------------------------------//	
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
        if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

        //----------------------------------------//	
            
        double weight = fScaleFactor*Units*A*Weight;
        if (fOutputFile.Contains("GiBUU")) { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 105 samples

        //----------------------------------------//	

        // Signal definition
        if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
        if (cc != 1) { continue; } // make sure that we have only CC interactions		

        // CC2p0pi event selection

        // Loop over final state particles
        int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0;
        int MuonTagging = 0, ElectronTagging = 0, PhotonTagging = 0;
        vector <int> ProtonID; ProtonID.clear();
        vector <int> MuonID; MuonID.clear();

        for (int i = 0; i < nfsp; i++) {
            double pf = TMath::Sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            if (pdg[i] == 13 && (pf > 0.1 && pf < 1.2)) {
                MuonTagging++;
                MuonID.push_back(i);
            }
            if (pdg[i] == 2212 && (pf > 0.3 && pf < 1.)) {
                ProtonTagging++;
                ProtonID.push_back(i);
            }
            if (fabs(pdg[i]) == 211 && pf > 0.07)  {
                ChargedPionTagging++;
            }
            if (pdg[i] == 111)  {
                NeutralPionTagging++;
            }
            if (fabs(pdg[i]) == 11)  {
                ElectronTagging++;
            }
            if (fabs(pdg[i]) == 22)  {
                PhotonTagging++;
            }
        } // End of the loop over the final state particles

        // Loop over pre-FSI particles
        int NoFSIProtonTagging = 0, NoFSIChargedPionTagging = 0, NoFSINeutralPionTagging = 0;
        int NoFSIMuonTagging = 0, NoFSIElectronTagging = 0, NoFSIPhotonTagging = 0;
        vector <int> NoFSIProtonID; NoFSIProtonID.clear();
        vector <int> NoFSIMuonID; NoFSIMuonID.clear();

        for (int i = 0; i < nvertp; i++) {
            double NoFSIpf = TMath::Sqrt(px_vert[i]*px_vert[i] + py_vert[i]*py_vert[i] + pz_vert[i]*pz_vert[i]);
            if (pdg_vert[i] == 13 && (NoFSIpf > 0.1 && NoFSIpf < 1.2)) {
                NoFSIMuonTagging++;
                NoFSIMuonID.push_back(i);
            }
            if (pdg_vert[i] == 2212 && (NoFSIpf > 0.3 && NoFSIpf < 1.)) {
                NoFSIProtonTagging++;
                NoFSIProtonID.push_back(i);
            }
            if (fabs(pdg_vert[i]) == 211 && NoFSIpf > 0.07)  {
                NoFSIChargedPionTagging++;
            }
            if (pdg_vert[i] == 111)  {
                NoFSINeutralPionTagging++;
            }
            if (fabs(pdg_vert[i]) == 11)  {
                NoFSIElectronTagging++;
            }
            if (fabs(pdg_vert[i]) == 22)  {
                NoFSIPhotonTagging++;
            }
        } // End of loop over pre-FSI particles


        // Check if signal definition for final state particles is satisfied
        if (
            ProtonTagging == 2 &&
            ChargedPionTagging == 0 &&
            NeutralPionTagging == 0 && 
            MuonTagging == 1
        ) {
            //----------------------------------------//	

            // https://arxiv.org/pdf/2106.15809.pdf

            CounterEventsPassedSelection++;
            
            // Classify the events based on the interaction type

            int genie_mode = -1.;
            if (TMath::Abs(Mode) == 1) { CounterQEEventsPassedSelection++; genie_mode = 1; } // QE
            else if (TMath::Abs(Mode) == 2) { CounterMECEventsPassedSelection++; genie_mode = 2; } // MEC
            else if (
                TMath::Abs(Mode) == 10 ||
                TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
                TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
                ) { CounterRESEventsPassedSelection++; genie_mode = 3; } // RES
            else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { CounterDISEventsPassedSelection++; genie_mode = 4; } // DIS
            else if (TMath::Abs(Mode) == 16) { CounterCOHEventsPassedSelection++; genie_mode = 5;} // COH
            else { continue; } 

            // Feb 8 2022: Only case that is not covered is 15 = diffractive

            //----------------------------------------//

            // Create momentum vectors and helper
            TVector3 Muon(px[MuonID[0]], py[MuonID[0]], pz[MuonID[0]]);
            TVector3 LeadingProton(px[ProtonID[0]], py[ProtonID[0]], pz[ProtonID[0]]);
            TVector3 RecoilProton(px[ProtonID[1]], py[ProtonID[1]], pz[ProtonID[1]]);
            TwoPTools Helper(Muon, LeadingProton, RecoilProton);

            // Get variables of interest
            double LeadingProtonCosTheta = Helper.ReturnLeadingProtonCosTheta();
            double RecoilProtonCosTheta = Helper.ReturnRecoilProtonCosTheta();
            double LeadingProtonMomentum = Helper.ReturnLeadingProtonMomentum();
            double RecoilProtonMomentum = Helper.ReturnRecoilProtonMomentum();
            double MuonMomentum = Helper.ReturnMuonMomentum();
            double CosOpeningAngleProtons = Helper.ReturnCosOpeningAngleProtons();
            double CosOpeningAngleMuonTotalProton = Helper.ReturnCosOpeningAngleMuonTotalProton();
            double TransverseMomentum = Helper.ReturnTransverseMomentum();
            double DeltaAlphaT = Helper.ReturnDeltaAlphaT();

            //----------------------------------------//

            // Filling in the histo regardless of interaction mode

            TrueMuonCosThetaPlot[0]->Fill(CosLep,weight);
            TrueLeadingProtonCosThetaPlot[0]->Fill(LeadingProtonCosTheta,weight);
            TrueRecoilProtonCosThetaPlot[0]->Fill(RecoilProtonCosTheta,weight);
            TrueLeadingProtonMomentumPlot[0]->Fill(LeadingProtonMomentum,weight);
            TrueRecoilProtonMomentumPlot[0]->Fill(RecoilProtonMomentum,weight);
            TrueMuonMomentumPlot[0]->Fill(MuonMomentum,weight);
            TrueCosOpeningAngleProtonsPlot[0]->Fill(CosOpeningAngleProtons,weight);
            TrueCosOpeningAngleMuonTotalProtonPlot[0]->Fill(CosOpeningAngleMuonTotalProton,weight);
            TrueTransverseMomentumPlot[0]->Fill(TransverseMomentum,weight);
            TrueDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);

            //----------------------------------------//

            // Filling in the histo based on the interaction mode

            TrueMuonCosThetaPlot[genie_mode]->Fill(CosLep,weight);
            TrueLeadingProtonCosThetaPlot[genie_mode]->Fill(LeadingProtonCosTheta,weight);
            TrueRecoilProtonCosThetaPlot[genie_mode]->Fill(RecoilProtonCosTheta,weight);
            TrueLeadingProtonMomentumPlot[genie_mode]->Fill(LeadingProtonMomentum,weight);
            TrueRecoilProtonMomentumPlot[genie_mode]->Fill(RecoilProtonMomentum,weight);
            TrueMuonMomentumPlot[genie_mode]->Fill(MuonMomentum,weight);
            TrueCosOpeningAngleProtonsPlot[genie_mode]->Fill(CosOpeningAngleProtons,weight);
            TrueCosOpeningAngleMuonTotalProtonPlot[genie_mode]->Fill(CosOpeningAngleMuonTotalProton,weight);
            TrueTransverseMomentumPlot[genie_mode]->Fill(TransverseMomentum,weight);
            TrueDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);

            //----------------------------------------//

            // Check for underflow/overflow for double differential plots
            if (CosLep < TwoDArrayNBinsMuonCosTheta[0]) { CosLep = (TwoDArrayNBinsMuonCosTheta[0] + TwoDArrayNBinsMuonCosTheta[1]) / 2.; }
            if (CosLep > TwoDArrayNBinsMuonCosTheta[TwoDNBinsMuonCosTheta]) { CosLep = (TwoDArrayNBinsMuonCosTheta[TwoDNBinsMuonCosTheta] + TwoDArrayNBinsMuonCosTheta[TwoDNBinsMuonCosTheta - 1]) / 2.; }

            if (TransverseMomentum < TwoDArrayTransverseMomentum[0]) { TransverseMomentum = (TwoDArrayTransverseMomentum[0] + TwoDArrayTransverseMomentum[1]) / 2.; }
            if (TransverseMomentum > TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum]) { TransverseMomentum = (TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum] + TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum - 1]) / 2.; }

            if (DeltaAlphaT < TwoDArrayDeltaAlphaT[0]) { DeltaAlphaT = (TwoDArrayDeltaAlphaT[0] + TwoDArrayDeltaAlphaT[1]) / 2.; }
            if (DeltaAlphaT > TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT]) { DeltaAlphaT = (TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT] + TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT - 1]) / 2.; }

            // Double differential indices
            int MuonCosThetaTwoDIndex = tools.ReturnIndex(CosLep, TwoDArrayNBinsMuonCosTheta);
            int SerialTransverseMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                TransverseMomentum
            );
            int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                DeltaAlphaT
            );

            // Fill in the histograms for double differential plots
            TrueSerialTransverseMomentum_InMuonCosThetaPlot[0]->Fill(SerialTransverseMomentumInMuonCosThetaIndex,weight);
            TrueSerialDeltaAlphaT_InMuonCosThetaPlot[0]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);

            TrueSerialTransverseMomentum_InMuonCosThetaPlot[genie_mode]->Fill(SerialTransverseMomentumInMuonCosThetaIndex,weight);
            TrueSerialDeltaAlphaT_InMuonCosThetaPlot[genie_mode]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);

            //----------------------------------------//
        } // End final state particles check

        // Check if signal definition for pre-FSI particles is satisfied
        else if (
            NoFSIProtonTagging == 2 &&
            NoFSIChargedPionTagging == 0 &&
            NoFSINeutralPionTagging == 0 && 
            NoFSIMuonTagging == 1
        ) {
            // Classify events based on interaction type
            int NoFSIgenie_mode = -1.;
            if (TMath::Abs(Mode) == 1) { NoFSIgenie_mode = 1; } // QE
            else if (TMath::Abs(Mode) == 2) { NoFSIgenie_mode = 2; } // MEC
            else if (
                TMath::Abs(Mode) == 10 ||
                TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
                TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
                ) { 
                    NoFSIgenie_mode = 3; 
                    if (fOutputFile.Contains("GiBUU")) { RESMode[Mode]++; NoFSICounterRESEventsPassedSelection++; }
                } // RES
            else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { NoFSIgenie_mode = 4; } // DIS
            else if (TMath::Abs(Mode) == 16) { NoFSIgenie_mode = 5;} // COH
            else { continue; } 

            // Create momentum vectors and helper
            TVector3 NoFSIMuon(px_vert[NoFSIMuonID[0]], py_vert[NoFSIMuonID[0]], pz_vert[NoFSIMuonID[0]]);
            TVector3 NoFSILeadingProton(px_vert[NoFSIProtonID[0]], py_vert[NoFSIProtonID[0]], pz_vert[NoFSIProtonID[0]]);
            TVector3 NoFSIRecoilProton(px_vert[NoFSIProtonID[1]], py_vert[NoFSIProtonID[1]], pz_vert[NoFSIProtonID[1]]);
            TwoPTools NoFSIHelper(NoFSIMuon, NoFSILeadingProton, NoFSIRecoilProton);

            // Filling in the histo regardless of interaction mode
            TrueNoFSILeadingProtonCosThetaPlot[0]->Fill(NoFSIHelper.ReturnLeadingProtonCosTheta(),weight);
            TrueNoFSIRecoilProtonCosThetaPlot[0]->Fill(NoFSIHelper.ReturnRecoilProtonCosTheta(),weight);
            TrueNoFSILeadingProtonMomentumPlot[0]->Fill(NoFSIHelper.ReturnLeadingProtonMomentum(),weight);
            TrueNoFSIRecoilProtonMomentumPlot[0]->Fill(NoFSIHelper.ReturnRecoilProtonMomentum(),weight);
            TrueNoFSIMuonMomentumPlot[0]->Fill(NoFSIHelper.ReturnMuonMomentum(),weight);
            TrueNoFSICosOpeningAngleProtonsPlot[0]->Fill(NoFSIHelper.ReturnCosOpeningAngleProtons(),weight);
            TrueNoFSICosOpeningAngleMuonTotalProtonPlot[0]->Fill(NoFSIHelper.ReturnCosOpeningAngleMuonTotalProton(),weight);
            TrueNoFSITransverseMomentumPlot[0]->Fill(NoFSIHelper.ReturnTransverseMomentum(),weight);
            TrueNoFSIDeltaAlphaTPlot[0]->Fill(NoFSIHelper.ReturnDeltaAlphaT(),weight);

            // Filling in the histo based on the interaction mode
            TrueNoFSILeadingProtonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnLeadingProtonCosTheta(),weight);
            TrueNoFSIRecoilProtonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnRecoilProtonCosTheta(),weight);
            TrueNoFSILeadingProtonMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnLeadingProtonMomentum(),weight);
            TrueNoFSIRecoilProtonMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnRecoilProtonMomentum(),weight);
            TrueNoFSIMuonMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnMuonMomentum(),weight);
            TrueNoFSICosOpeningAngleProtonsPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnCosOpeningAngleProtons(),weight);
            TrueNoFSICosOpeningAngleMuonTotalProtonPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnCosOpeningAngleMuonTotalProton(),weight);
            TrueNoFSITransverseMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnTransverseMomentum(),weight);
            TrueNoFSIDeltaAlphaTPlot[NoFSIgenie_mode]->Fill(NoFSIHelper.ReturnDeltaAlphaT(),weight);
        } // End pre-FSI particles check
        else { continue; }
    } // End of the loop over the events

    //----------------------------------------//	

    std::cout << "Percentage of events passing the selection cuts = " << 
    double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting QE events = " << 
    double(CounterQEEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting MEC events = " << 
    double(CounterMECEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting RES events = " << 
    double(CounterRESEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting DIS events = " << 
    double(CounterDISEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting COH events = " << 
    double(CounterCOHEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	

    if (fOutputFile.Contains("GiBUU")) {
        std::cout << "Mode percentage for RES events" << std::endl;
        std::cout << "10: " << double(RESMode[10]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "11: " << double(RESMode[11]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "12: " << double(RESMode[12]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "13: " << double(RESMode[13]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "17: " << double(RESMode[17]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "22: " << double(RESMode[22]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "23: " << double(RESMode[23]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
    }

    //----------------------------------------//	
    //----------------------------------------//	

    // Division by bin width to get the cross sections	
    // Loop over the interaction processes

    for (int inte = 0; inte < NInte; inte++) {

        //----------------------------------------//
    
        Reweight(TrueMuonCosThetaPlot[inte]);
        Reweight(TrueLeadingProtonCosThetaPlot[inte]);
        Reweight(TrueRecoilProtonCosThetaPlot[inte]);
        Reweight(TrueLeadingProtonMomentumPlot[inte]);
        Reweight(TrueRecoilProtonMomentumPlot[inte]);
        Reweight(TrueMuonMomentumPlot[inte]);
        Reweight(TrueCosOpeningAngleProtonsPlot[inte]);
        Reweight(TrueCosOpeningAngleMuonTotalProtonPlot[inte]);
        Reweight(TrueTransverseMomentumPlot[inte]);
        Reweight(TrueDeltaAlphaTPlot[inte]);

        Reweight(TrueNoFSILeadingProtonCosThetaPlot[inte]);
        Reweight(TrueNoFSIRecoilProtonCosThetaPlot[inte]);
        Reweight(TrueNoFSILeadingProtonMomentumPlot[inte]);
        Reweight(TrueNoFSIRecoilProtonMomentumPlot[inte]);
        Reweight(TrueNoFSIMuonMomentumPlot[inte]);
        Reweight(TrueNoFSICosOpeningAngleProtonsPlot[inte]);
        Reweight(TrueNoFSICosOpeningAngleMuonTotalProtonPlot[inte]);
        Reweight(TrueNoFSITransverseMomentumPlot[inte]);
        Reweight(TrueNoFSIDeltaAlphaTPlot[inte]);

        //----------------------------------------//

    } // End of the loop over the interaction processes		

    //----------------------------------------//		
        
    file->cd();
    file->Write();
    fFile->Close();

    std::cout << std::endl;
    std::cout << "File " << FileNameAndPath +" has been created " << std::endl; 
    std::cout << std::endl;

    std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

    //----------------------------------------//		

} // End of the program

//----------------------------------------//		

void Reweight(TH1D* h) {
    int NBins = h->GetXaxis()->GetNbins();

    for (int i = 0; i < NBins; i++) {

        double CurrentEntry = h->GetBinContent(i+1);
        double NewEntry = CurrentEntry / h->GetBinWidth(i+1);

        double CurrentError = h->GetBinError(i+1);
        double NewError = CurrentError / h->GetBinWidth(i+1);

        h->SetBinContent(i+1,NewEntry); 
        h->SetBinError(i+1,NewError); 
        //h->SetBinError(i+1,0.000001); 

    }
}
//----------------------------------------//		
