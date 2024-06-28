// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1D.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionEfficiency() {
    // Some useful variables for later.
    const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const double TargetPOT(6.6e20);

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // Create the binning schemes for the Vars we wish to plot.
    const Binning bPrimaryEnergy = Binning::Simple(1, 0, 3.0); // one bin
    const Binning bAngleBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bDeltaAlphaBins = Binning::Simple(20, 0.0, 180.0);
    const Binning bTransverseMomentumBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bMuonMomentumBins = Binning::Simple(20, 0.1, 1.2);
    const Binning bProtonMomentumBins = Binning::Simple(20, 0.3, 1.0);

    // Double differential bins
    Tools tools; // tools for double differential bins

    const Binning bTransverseMomentumInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices)
    );
    const Binning bDeltaAlphaTInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)
    );
    const Binning bCosOpeningAngleProtonsInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices)
    );
    const Binning bCosOpeningAngleMuonTotalProtonInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices)
    );

    // We will create efficiency plots using true variables and definining signal efficiency
    // as the number of reconstructed the events that pass our signal definition and are true
    // signal events over the total true signal events; these two histograms are plotted

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Root file to store objects in
    TString RootFilePath = "/pnfs/sbnd/persistent/users/epelaez/CAFAnaOutput/SelectionEfficiency.root";
    TFile* SaveFile = new TFile(RootFilePath, "RECREATE");

    // Vectors to fill with variables and variable information to plot
    std::vector<TruthVar> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Muon angle
    Vars.push_back(kTruthMuonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

    // Leading proton angle
    Vars.push_back(kTruthLeadingProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthLeadingProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L}})");

    // Recoil proton angle
    Vars.push_back(kTruthRecoilProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthRecoilProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{R}})");

    // Opening angle between protons
    Vars.push_back(kTruthCosOpeningAngleProtons); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthCosOpeningAngleProtons"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");

    // Opening angle between muon and total proton
    Vars.push_back(kTruthCosOpeningAngleMuonTotalProton); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthCosOpeningAngleMuonTotalProton"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");

    // Delta alpha transverse
    Vars.push_back(kTruthDeltaAlphaT); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("TruthDeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Transverse momentum
    Vars.push_back(kTruthTransverseMomentum); VarBins.push_back(bTransverseMomentumBins);
    PlotNames.push_back("TruthTransverseMomentum"); VarLabels.push_back("#delta P_{T}");

    // Muon momentum 
    Vars.push_back(kTruthMuonMomentum); VarBins.push_back(bMuonMomentumBins);
    PlotNames.push_back("TruthMuonMomentum"); VarLabels.push_back("|#vec{p}_{#mu}|");

    // Leading proton momentum 
    Vars.push_back(kTruthLeadingProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("TruthLeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back(kTruthRecoilProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("TruthRecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    // Construct all spectra
    std::vector<std::tuple<std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto TrueSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kTruthIsSignal, kNoSpillCut);
        auto RecoTrueSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kTruthIsSignal, kNoSpillCut, kRecoIsSignal);
        Spectra.push_back({std::move(TrueSignals), std::move(RecoTrueSignals)});
    }

    NuLoader.Go();

    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto& [TrueSignals, RecoTrueSignals] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1* TrueHisto = TrueSignals->ToTH1(TargetPOT);
        TH1* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);

        TrueHisto->GetYaxis()->SetTitle("Signal efficiency");
        RecoTrueHisto->GetYaxis()->SetTitle("Signal efficiency");

        TEfficiency* Eff = new TEfficiency(*RecoTrueHisto, *TrueHisto);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        Eff->SetMarkerStyle(21);
        Eff->SetMarkerColor(kBlack);
        Eff->SetTitle("All events");
        Eff->Draw("AP");
        gPad->Update();
        Eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(VarBins.at(i).Min(),VarBins.at(i).Max());

        // Save as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Efficiency/"+PlotNames[i]+".png");

        // Save to root file
        SaveFile->WriteObject(Eff, PlotNames[i]+"_eff");

        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}