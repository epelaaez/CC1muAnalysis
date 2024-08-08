// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionResolution() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    // Load root file with histograms
    TString HistoFile = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> File(TFile::Open(HistoFile));

    ////////////
    // Variables
    ////////////

    // Vectors to fill with variable pairs and information to plot
    std::vector<std::tuple<Var, Var, TruthVar>> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

    // Dummy variable
    Vars.push_back({kEventCount, kEventCount, kTrueEventCount}); VarBins.push_back(bEventCount); 
    PlotNames.push_back("EventCount"); VarLabels.push_back("single bin");

    // Muon angle
    Vars.push_back({kMuonCosTheta, kRecoTruthMuonCosTheta, kTruthMuonCosTheta}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("MuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

    // Leading proton angle
    Vars.push_back({kLeadingProtonCosTheta, kRecoTruthLeadingProtonCosTheta, kTruthLeadingProtonCosTheta}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("LeadingProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L}})");

    // Recoil proton angle
    Vars.push_back({kRecoilProtonCosTheta, kRecoTruthRecoilProtonCosTheta, kTruthRecoilProtonCosTheta}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("RecoilProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{R}})");

    // Opening angle between protons
    Vars.push_back({kCosOpeningAngleProtons, kRecoTruthCosOpeningAngleProtons, kTruthCosOpeningAngleProtons}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleProtons"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");

    // Opening angle between muon and total proton
    Vars.push_back({kCosOpeningAngleMuonTotalProton, kRecoTruthCosOpeningAngleMuonTotalProton, kTruthCosOpeningAngleMuonTotalProton}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleMuonTotalProton"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");

    // Delta alpha transverse
    Vars.push_back({kDeltaAlphaT, kRecoTruthDeltaAlphaT, kTruthDeltaAlphaT}); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("DeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Transverse momentum
    Vars.push_back({kTransverseMomentum, kRecoTruthTransverseMomentum, kTruthTransverseMomentum}); VarBins.push_back(bTransverseMomentumBins);
    PlotNames.push_back("TransverseMomentum"); VarLabels.push_back("#delta P_{T}");

    // Muon momentum 
    Vars.push_back({kMuonMomentum, kRecoTruthMuonMomentum, kTruthMuonMomentum}); VarBins.push_back(bMuonMomentumBins);
    PlotNames.push_back("MuonMomentum"); VarLabels.push_back("|#vec{p}_{#mu}|");

    // Leading proton momentum 
    Vars.push_back({kLeadingProtonMomentum, kRecoTruthLeadingProtonMomentum, kTruthLeadingProtonMomentum}); VarBins.push_back(bLeadingProtonMomentumBins);
    PlotNames.push_back("LeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back({kRecoilProtonMomentum, kRecoTruthRecoilProtonMomentum, kTruthRecoilProtonMomentum}); VarBins.push_back(bRecoilProtonMomentumBins);
    PlotNames.push_back("RecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    ////////////////////////////////
    // Double differential variables
    ////////////////////////////////

    // // Serial transverse momentum in muon cos theta
    // Vars.push_back({kTransverseMomentumInMuonCosTheta, kRecoTruthTransverseMomentumInMuonCosTheta, kTruthTransverseMomentumInMuonCosTheta}); VarBins.push_back(bTransverseMomentumInMuonCosTheta);
    // PlotNames.push_back("SerialTransverseMomentum_InMuonCosTheta"); VarLabels.push_back("#delta P_{T} (bin #)");

    // // Delta alpha transverse in muon cos theta
    // Vars.push_back({kDeltaAlphaTInMuonCosTheta, kRecoTruthDeltaAlphaTInMuonCosTheta, kTruthDeltaAlphaTInMuonCosTheta}); VarBins.push_back(bDeltaAlphaTInMuonCosTheta);
    // PlotNames.push_back("SerialDeltaAlphaT_InMuonCosTheta"); VarLabels.push_back("#delta #alpha_{T} (bin #)");

    // // Opening angle between protons in muon cos theta
    // Vars.push_back({kCosOpeningAngleProtonsInMuonCosTheta, kRecoTruthCosOpeningAngleProtonsInMuonCosTheta, kTruthCosOpeningAngleProtonsInMuonCosTheta}); VarBins.push_back(bCosOpeningAngleProtonsInMuonCosTheta);
    // PlotNames.push_back("SerialCosOpeningAngleProtons_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)");
    
    // // Opening angle between muon and protons in muon cos theta
    // Vars.push_back({kCosOpeningAngleMuonTotalProtonInMuonCosTheta, kRecoTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta, kTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta}); VarBins.push_back(bCosOpeningAngleMuonTotalProtonInMuonCosTheta);
    // PlotNames.push_back("SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)");


    const int NVars = PlotNames.size();

    // Construct spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;
    for (int iVar = 0; iVar < NVars; ++iVar) {
        Var kRecoVar = std::get<0>(Vars.at(iVar)); Var kRecoTruthVar = std::get<1>(Vars.at(iVar));

        double MaxDiff = TMath::Abs(VarBins.at(iVar).Max() - VarBins.at(iVar).Min());
        const Binning bResolutionBin = Binning::Simple(41, -MaxDiff, MaxDiff);

        auto ResolutionValues = std::make_unique<Spectrum>(VarLabels.at(iVar), bResolutionBin, NuLoader, kRecoVar - kRecoTruthVar, kNoSpillCut, kRecoIsTrueReco);
        auto ActualResolutionValues = std::make_unique<Spectrum>(VarLabels.at(iVar), bActualResolution, NuLoader, (kRecoVar - kRecoTruthVar) / kRecoTruthVar, kNoSpillCut, kRecoIsTrueReco);
        Spectra.push_back({std::move(ResolutionValues), std::move(ActualResolutionValues)});
    }
    NuLoader.Go();

    // Loop over variables
    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto& [ResolutionValues, ActualResolutionValues] = Spectra.at(iVar);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1* ResolutionHist = ResolutionValues->ToTH1(TargetPOT);
        TH1* ActualResolutionHist = ActualResolutionValues->ToTH1(TargetPOT);

        // Manage under/overflow bins
        ResolutionHist->SetBinContent(ResolutionHist->GetNbinsX(), ResolutionHist->GetBinContent(ResolutionHist->GetNbinsX()) + ResolutionHist->GetBinContent(ResolutionHist->GetNbinsX() + 1));
        // ActualResolutionHist->SetBinContent(ActualResolutionHist->GetNbinsX(), ActualResolutionHist->GetBinContent(ActualResolutionHist->GetNbinsX()) + ActualResolutionHist->GetBinContent(ActualResolutionHist->GetNbinsX() + 1));

        ResolutionHist->SetBinContent(1, ResolutionHist->GetBinContent(0) + ResolutionHist->GetBinContent(1));
        // ActualResolutionHist->SetBinContent(1, ActualResolutionHist->GetBinContent(0) + ActualResolutionHist->GetBinContent(1));

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        ResolutionHist->GetXaxis()->SetTitleFont(FontStyle);
        ResolutionHist->GetXaxis()->SetLabelFont(FontStyle);
        ResolutionHist->GetXaxis()->SetNdivisions(8);
        ResolutionHist->GetXaxis()->SetLabelSize(TextSize);
        ResolutionHist->GetXaxis()->SetTitleSize(TextSize);
        ResolutionHist->GetXaxis()->SetTitleOffset(1.1);
        ResolutionHist->GetXaxis()->CenterTitle();
        ResolutionHist->GetXaxis()->SetTitle((VarLabels.at(iVar) + " difference").c_str());

        ResolutionHist->GetYaxis()->SetTitleFont(FontStyle);
        ResolutionHist->GetYaxis()->SetLabelFont(FontStyle);
        ResolutionHist->GetYaxis()->SetNdivisions(6);
        ResolutionHist->GetYaxis()->SetLabelSize(TextSize);
        ResolutionHist->GetYaxis()->SetTitleSize(TextSize);
        ResolutionHist->GetYaxis()->SetTitleOffset(1.3);
        ResolutionHist->GetYaxis()->SetTickSize(0);
        ResolutionHist->GetYaxis()->CenterTitle();
        ResolutionHist->GetYaxis()->SetTitle("# events");

        ActualResolutionHist->GetXaxis()->SetTitleFont(FontStyle);
        ActualResolutionHist->GetXaxis()->SetLabelFont(FontStyle);
        ActualResolutionHist->GetXaxis()->SetNdivisions(8);
        ActualResolutionHist->GetXaxis()->SetLabelSize(TextSize);
        ActualResolutionHist->GetXaxis()->SetTitleSize(TextSize);
        ActualResolutionHist->GetXaxis()->SetTitleOffset(1.1);
        ActualResolutionHist->GetXaxis()->CenterTitle();
        ActualResolutionHist->GetXaxis()->SetTitle((VarLabels.at(iVar) + " range resolution").c_str());

        ActualResolutionHist->GetYaxis()->SetTitleFont(FontStyle);
        ActualResolutionHist->GetYaxis()->SetLabelFont(FontStyle);
        ActualResolutionHist->GetYaxis()->SetNdivisions(6);
        ActualResolutionHist->GetYaxis()->SetLabelSize(TextSize);
        ActualResolutionHist->GetYaxis()->SetTitleSize(TextSize);
        ActualResolutionHist->GetYaxis()->SetTitleOffset(1.3);
        ActualResolutionHist->GetYaxis()->SetTickSize(0);
        ActualResolutionHist->GetYaxis()->CenterTitle();
        ActualResolutionHist->GetYaxis()->SetTitle("# events");

        PlotCanvas->cd();
        ResolutionHist->Draw("hist e1");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Resolution/Diff"+PlotNames[iVar]+".png");

        PlotCanvas->cd();
        ActualResolutionHist->Draw("hist e1");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Resolution/Res"+PlotNames[iVar]+".png");

        delete PlotCanvas;
    }
}