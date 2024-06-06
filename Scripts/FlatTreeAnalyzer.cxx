#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"
#include "TwoPTools.h"

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

//Function to divide by the bin width and to get xsecs
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

	// Output file
	TString Directory = "/pnfs/sbnd/persistent/users/epelaez/HighSamples/";
	TString FileNameAndPath = Directory+"FlatTree/FlatTreeAnalyzerOutput_"+fOutputFile+".root";

	// TString FileNameAndPath = "OutputFilesHighStats/FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	// TString FileNameAndPath = "OutputFiles/FlatTreeAnalyzerOutput_"+fOutputFile+".root";
	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//

	// Plot declaration

	TH1D* TrueMuonCosThetaPlot[NInte];
	TH1D* TrueLeadingProtonCosThetaPlot[NInte];
	TH1D* TrueRecoilProtonCosThetaPlot[NInte];
	TH1D* TrueLeadingProtonMomentumPlot[NInte];
	TH1D* TrueRecoilProtonMomentumPlot[NInte];
	TH1D* TrueMuonMomentumPlot[NInte];
	TH1D* TrueCosOpeningAngleProtonsPlot[NInte];
	TH1D* TrueCosOpeningAngleMuonTotalProtonPlot[NInte];
	TH1D* TrueTransverseMomentumPlot[NInte];

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#mu})",10,-1.,1.);
		TrueLeadingProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueLeadingProtonCosThetaPlot",";cos(#theta_{#vec{p}_{L}})",10,-1.,1.);
		TrueRecoilProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueRecoilProtonCosThetaPlot",";cos(#theta_{#vec{p}_{R}})",10,-1.,1.);
		TrueLeadingProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueLeadingProtonMomentumPlot",";|#vec{p}_{L}|",10,0.3,1);
		TrueRecoilProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueRecoilProtonMomentumPlot",";|#vec{p}_{R}|",10,0.3,1);
		TrueMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonMomentumPlot",";|#vec{p}_{#mu}|",10,0.1,1.2);
		TrueCosOpeningAngleProtonsPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleProtonsPlot",";cos(#theta_{#vec{p}_{L},#vec{p}_{R}})",10,-1.,1.);
		TrueCosOpeningAngleMuonTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleMuonTotalProtonPlot",";cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})",10,-1.,1.);
		TrueTransverseMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueTransverseMomentumPlot",";#delta P_{T}",10,0.,1.);

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

	//----------------------------------------//
	
	// Loop over the events

	for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  //----------------------------------------//	
	
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
	  if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

	  //----------------------------------------//	
		
	  double weight = fScaleFactor*Units*A*Weight;
		if (fOutputFile == "GiBUU") { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 105 samples

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  if (cc != 1) { continue; } // make sure that we have only CC interactions		

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0;
		int ElectronTagging = 0, PhotonTagging = 0;
	  vector <int> ProtonID; ProtonID.clear();
	  vector <int> MuonID; MuonID.clear();

	  // Example selection with CC1p0pi (units in GeV/c)
	  // Loop over final state particles

		// Modified to CC2p0pi
	  for (int i = 0; i < nfsp; i++) {
		
	    double pf = TMath::Sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);

	    if (pdg[i] == 13 && (pf > 0.1 && pf < 1.2)) {

	      MuonTagging ++;
	      MuonID.push_back(i);

	    }

	    if (pdg[i] == 2212 && (pf > 0.3 && pf < 1.) ) {

	      ProtonTagging ++;
	      ProtonID.push_back(i);

	    }

	    if (fabs(pdg[i]) == 211 && pf > 0.07)  {

	      ChargedPionTagging ++;

	    }

	    if (pdg[i] == 111)  {

	      NeutralPionTagging ++;

	    }

	    if (fabs(pdg[i]) == 11)  {

	      ElectronTagging ++;

	    }

	    if (fabs(pdg[i]) == 22)  {

	      PhotonTagging ++;

	    }


	  } // End of the loop over the final state particles

	  // If the signal definition is not satisfied, continue

	  if (
 		ProtonTagging != 2 || 
		ChargedPionTagging != 0 || 
		NeutralPionTagging != 0 || 
		MuonTagging != 1
	) { continue; }

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

	  //----------------------------------------//

	  // filling in the histo regardless of interaction mode

		TrueMuonCosThetaPlot[0]->Fill(CosLep,weight);
		TrueLeadingProtonCosThetaPlot[0]->Fill(Helper.ReturnLeadingProtonCosTheta(),weight);
		TrueRecoilProtonCosThetaPlot[0]->Fill(Helper.ReturnRecoilProtonCosTheta(),weight);
		TrueLeadingProtonMomentumPlot[0]->Fill(Helper.ReturnLeadingProtonMomentum(),weight);
		TrueRecoilProtonMomentumPlot[0]->Fill(Helper.ReturnRecoilProtonMomentum(),weight);
		TrueMuonMomentumPlot[0]->Fill(Helper.ReturnMuonMomentum(),weight);
		TrueCosOpeningAngleProtonsPlot[0]->Fill(Helper.ReturnCosOpeningAngleProtons(),weight);
		TrueCosOpeningAngleMuonTotalProtonPlot[0]->Fill(Helper.ReturnCosOpeningAngleMuonTotalProton(),weight);
		TrueTransverseMomentumPlot[0]->Fill(Helper.ReturnTransverseMomentum(),weight);
 
	  //----------------------------------------//

	  // filling in the histo based on the interaction mode

	  TrueMuonCosThetaPlot[genie_mode]->Fill(CosLep,weight);
		TrueLeadingProtonCosThetaPlot[genie_mode]->Fill(Helper.ReturnLeadingProtonCosTheta(),weight);
		TrueRecoilProtonCosThetaPlot[genie_mode]->Fill(Helper.ReturnRecoilProtonCosTheta(),weight);
		TrueLeadingProtonMomentumPlot[genie_mode]->Fill(Helper.ReturnLeadingProtonMomentum(),weight);
		TrueRecoilProtonMomentumPlot[genie_mode]->Fill(Helper.ReturnRecoilProtonMomentum(),weight);
		TrueMuonMomentumPlot[genie_mode]->Fill(Helper.ReturnMuonMomentum(),weight);
		TrueCosOpeningAngleProtonsPlot[genie_mode]->Fill(Helper.ReturnCosOpeningAngleProtons(),weight);
		TrueCosOpeningAngleMuonTotalProtonPlot[genie_mode]->Fill(Helper.ReturnCosOpeningAngleMuonTotalProton(),weight);
		TrueTransverseMomentumPlot[genie_mode]->Fill(Helper.ReturnTransverseMomentum(),weight);

	  //----------------------------------------//
	
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
