// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/ISyst.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionSystematics(int SystIndex) {
    std::string SystName = std::get<0>(SystsVector.at(SystIndex));
    int SystNUniv = std::get<1>(SystsVector.at(SystIndex));

    std::cout << "==============================" << std::endl;
    std::cout << "Systematic name: " << SystName << ", number of universes: " << SystNUniv <<  std::endl;
    std::cout << "==============================" << std::endl;

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // Get integrated flux
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root");
    TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
    double IntegratedFlux = (HistoFlux->Integral() * TargetPOT / POTPerSpill / Nominal_UB_XY_Surface);    

    // Some useful variables for later.
    // const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const std::string TargetFile = "/pnfs/sbnd/persistent/users/apapadop/CAF_Files/*.flat.caf.root";

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // We now create overlaid plots for several reconstructed variables and three lines:
    //     1. all selected reconstructed events
    //     2. reco signal events
    //     3. reco background events

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematics.root";
    TFile* SaveFile = new TFile(RootFilePath, "RECREATE");

    // Vectors to fill with variables and variable information to plot
    std::vector<Var> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;
    
    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Muon angle
    Vars.push_back(kMuonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("MuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

    // Leading proton angle
    Vars.push_back(kLeadingProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("LeadingProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L}})");

    // Recoil proton angle
    Vars.push_back(kRecoilProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("RecoilProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{R}})");

    // Opening angle between protons
    Vars.push_back(kCosOpeningAngleProtons); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleProtons"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");

    // Opening angle between muon and total proton
    Vars.push_back(kCosOpeningAngleMuonTotalProton); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleMuonTotalProton"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");

    // Delta alpha transverse
    Vars.push_back(kDeltaAlphaT); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("DeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Transverse momentum
    Vars.push_back(kTransverseMomentum); VarBins.push_back(bTransverseMomentumBins);
    PlotNames.push_back("TransverseMomentum"); VarLabels.push_back("#delta P_{T}");

    // Muon momentum 
    Vars.push_back(kMuonMomentum); VarBins.push_back(bMuonMomentumBins);
    PlotNames.push_back("MuonMomentum"); VarLabels.push_back("|#vec{p}_{#mu}|");

    // Leading proton momentum 
    Vars.push_back(kLeadingProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("LeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back(kRecoilProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("RecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    ////////////////////////////////
    // Double differential variables
    ////////////////////////////////

    // Serial transverse momentum in muon cos theta
    Vars.push_back(kTransverseMomentumInMuonCosTheta); VarBins.push_back(bTransverseMomentumInMuonCosTheta);
    PlotNames.push_back("SerialTransverseMomentum_InMuonCosTheta"); VarLabels.push_back("#delta P_{T} (bin #)");

    // Delta alpha transverse in muon cos theta
    Vars.push_back(kDeltaAlphaTInMuonCosTheta); VarBins.push_back(bDeltaAlphaTInMuonCosTheta);
    PlotNames.push_back("SerialDeltaAlphaT_InMuonCosTheta"); VarLabels.push_back("#delta #alpha_{T} (bin #)");

    // Opening angle between protons in muon cos theta
    Vars.push_back(kCosOpeningAngleProtonsInMuonCosTheta); VarBins.push_back(bCosOpeningAngleProtonsInMuonCosTheta);
    PlotNames.push_back("SerialCosOpeningAngleProtons_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)");
    
    // Opening angle between muon and protons in muon cos theta
    Vars.push_back(kCosOpeningAngleMuonTotalProtonInMuonCosTheta); VarBins.push_back(bCosOpeningAngleMuonTotalProtonInMuonCosTheta);
    PlotNames.push_back("SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)");

    // Construct all spectra
    std::vector<std::tuple<
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>
    >> Spectra;
    for (std::size_t iVar = 0; iVar < Vars.size(); iVar++) {
        // Create shift depending on number of universes, 6/10 for multisigma
	// and everything else is treated as nuniv
     	ISyst* syst = new SBNWeightSyst(SystName);
	std::vector<SystShifts> Shifts;
	SystShifts SigP1Shift(syst, +1);

	if (SystNUniv == 6 || SystNUniv == 10) {
    	    // Add +1 sigma shift
	    Shifts.push_back(SigP1Shift);
	} else {
	    // Add random Gaussian shifts
    	    for (int i = 0; i < SystNUniv; i++) {
		SystShifts RandomShift(syst, gRandom->Gaus(0,1));
    		Shifts.push_back(RandomShift);
    	    }
	}

	// Create reco spectrum with shift
    	auto RecoSpectra = std::make_unique<EnsembleSpectrum>(
	    NuLoader,
	    HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)),
	    kNoSpillCut,
	    kRecoIsSignal,
 	    Shifts
	);

	// Create reco true signal spectrum with shift
    	auto RecoTrueSpectra = std::make_unique<EnsembleSpectrum>(
	    NuLoader,
     	    HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)),
            kNoSpillCut,
            kRecoIsTrueReco,
            Shifts
        );

        // Create reco background spectrum with shift
        auto RecoBkgSpectra = std::make_unique<EnsembleSpectrum>(
            NuLoader,
    	    HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)),
	    kNoSpillCut,
	    kRecoIsBackground,
	    Shifts
	);
        
	// Add everything to main vector
        Spectra.push_back({std::move(RecoSpectra), std::move(RecoTrueSpectra), std::move(RecoBkgSpectra)});
    }

    // Load spectra
    NuLoader.Go();

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto& [RecoSpectra, RecoTrueSpectra, RecoBkgSpectra] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1D* RecoHisto = RecoSpectra->Nominal().ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSpectra->Nominal().ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = RecoBkgSpectra->Nominal().ToTH1(TargetPOT);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoHisto,"Reconstructed","l");
        RecoHisto->SetLineColor(kBlue+2);
        RecoHisto->SetLineWidth(4);

        // Style histograms
        RecoHisto->GetXaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetXaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetXaxis()->SetNdivisions(8);
        RecoHisto->GetXaxis()->SetLabelSize(TextSize);
        RecoHisto->GetXaxis()->SetTitleSize(TextSize);
        RecoHisto->GetXaxis()->SetTitleOffset(1.1);
        RecoHisto->GetXaxis()->CenterTitle();
        RecoHisto->GetXaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());

        RecoHisto->GetYaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetYaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetYaxis()->SetNdivisions(6);
        RecoHisto->GetYaxis()->SetLabelSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleOffset(1.3);
        RecoHisto->GetYaxis()->SetTickSize(0);
        RecoHisto->GetYaxis()->CenterTitle();

        double imax = RecoHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"True","l");
        RecoTrueHisto->SetLineColor(kRed+1);
        RecoTrueHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"Background","l");
        RecoBkgHisto->SetLineColor(kOrange+7);
        RecoBkgHisto->SetLineWidth(4);

	// Create covariance matrix
	std::string CovName = "Cov" + SystName;
	TH2* CovMatrix = new TH2D(
	    CovName.c_str(),
	    CovName.c_str(),
	    VarBins.at(i).NBins(),
	    VarBins.at(i).Min(),
	    VarBins.at(i).Max(),
	    VarBins.at(i).NBins(),
	    VarBins.at(i).Min(),
	    VarBins.at(i).Max()
	);

	// Loop over all universes
	int NUniv = (SystNUniv == 6 || SystNUniv == 10) ? 1 : SystNUniv;
    	for (int iUniv = 0; iUniv < NUniv; iUniv++) {
	    TH1* UnivRecoSpectrum = RecoSpectra->Universe(iUniv).ToTH1(TargetPOT);
    	    TH1* UnivRecoTrueSpectrum = RecoTrueSpectra->Universe(iUniv).ToTH1(TargetPOT);
            TH1* UnivRecoBkgSpectrum = RecoBkgSpectra->Universe(iUniv).ToTH1(TargetPOT);

    	    for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
		double XEventRateCV = (RecoHisto->GetBinContent(x) / (IntegratedFlux * NTargets)) * Units;
    		double XEventRateVar = (UnivRecoSpectrum->GetBinContent(x) / (IntegratedFlux * NTargets)) * Units;
		for (int y = 1; y <= x; y++) {
		    double YEventRateCV = (RecoHisto->GetBinContent(y) / (IntegratedFlux * NTargets)) * Units;
	    	    double YEventRateVar = (UnivRecoSpectrum->GetBinContent(y) / (IntegratedFlux * NTargets)) * Units; 

	    	    CovMatrix->Fill(
			RecoHisto->GetXaxis()->GetBinCenter(x),
	    		RecoHisto->GetXaxis()->GetBinCenter(y),
			((XEventRateVar - XEventRateCV) * (YEventRateVar - YEventRateCV)) / NUniv
		    );

	    	    CovMatrix->Fill(
			RecoHisto->GetXaxis()->GetBinCenter(y),
	    		RecoHisto->GetXaxis()->GetBinCenter(x),
			((XEventRateVar - XEventRateCV) * (YEventRateVar - YEventRateCV)) / NUniv
	 	    );
		}
	    }
	    // Save syst univ spectrum
    	    TString UnivString = TString(std::to_string(iUniv));
	    SaveFile->WriteObject(UnivRecoSpectrum, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco");
            SaveFile->WriteObject(UnivRecoTrueSpectrum, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_true");
            SaveFile->WriteObject(UnivRecoBkgSpectrum, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_bkg");
	}
	// Create directory for this sytematic if it does not exist yet
	std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/Uncertainties/"+SystName);
	    
	// Plot and save cov matrix
     	CovMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
 	CovMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());

 	PlotCanvas->cd();
	CovMatrix->Draw("colz text");
	PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/Cov"+PlotNames[i]+".png");
            
	// Get all error bands
	TGraphAsymmErrors* RecoErrorBand = RecoSpectra->ErrorBand(TargetPOT);
	TGraphAsymmErrors* RecoTrueErrorBand = RecoTrueSpectra->ErrorBand(TargetPOT);
	TGraphAsymmErrors* RecoBkgErrorBand = RecoBkgSpectra->ErrorBand(TargetPOT);

	PlotCanvas->cd();
	RecoHisto->Draw("hist");
	ana::DrawErrorBand(RecoHisto, RecoErrorBand);
	RecoTrueHisto->Draw("hist same");
	ana::DrawErrorBand(RecoTrueHisto, RecoTrueErrorBand);
	RecoBkgHisto->Draw("hist same");
	ana::DrawErrorBand(RecoBkgHisto, RecoBkgErrorBand);
	leg->Draw();

	// Save as png
	PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/"+PlotNames[i]+".png");
	SaveFile->WriteObject(CovMatrix, (TString)SystName+PlotNames[i]+"_cov");
	delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}
