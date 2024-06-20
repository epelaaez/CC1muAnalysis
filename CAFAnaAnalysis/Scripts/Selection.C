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
#include "TFile.h"
#include "TH1D.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

using namespace std;
using namespace ana;

void Selection()
{
    // Some useful variables for later.
    const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const double TargetPOT(6.6e20);

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // Create the binning schemes for the Vars we wish to plot.
    const Binning bPrimaryEnergy = Binning::Simple(20, 0, 3.0);
    const Binning bAngleBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bDeltaAlphaBins = Binning::Simple(20, 0.0, 180.0);

    // We now create overlaid plots for several reconstructed variables and three lines:
    //     1. all selected reconstructed events
    //     2. reco signal events
    //     3. reco background events

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Variables to plot
    std::vector<Var> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

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

    // Spectrum with all true signal events
    Spectrum TrueSignals("TrueSignals", bPrimaryEnergy, NuLoader, kTrueEnergy, kTruthIsSignal, kNoSpillCut, kNoCut);

    // Construct all spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kNoSpillCut, kRecoIsSignal); 
        auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kNoSpillCut, kRecoIsTrueReco); 
        auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kNoSpillCut, kRecoIsBackground); 
        Spectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(RecoBkgSignals)});
    }
    NuLoader.Go();

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto& [RecoSignals, RecoTrueSignals, RecoBkgSignals] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1D* RecoHisto = RecoSignals->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = RecoBkgSignals->ToTH1(TargetPOT);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TLegend* leg = new TLegend(0.13,0.88,0.50,0.99);
        leg->SetBorderSize(0);
        leg->SetNColumns(2);
        leg->SetMargin(0.2);
        leg->SetFillColor(0);

        TLegendEntry* legReco = leg->AddEntry(RecoHisto,"all reco","l");
        legReco->SetTextColor(602); // blue
        RecoHisto->SetLineColor(602);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"signal","l");
        legRecoTrue->SetTextColor(797); // orange
        RecoTrueHisto->SetLineColor(797);  

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"bkg","l");
        legRecoBkg->SetTextColor(417); // green
        RecoBkgHisto->SetLineColor(417);

        PlotCanvas->cd();
        RecoHisto->Draw("hist same");
        RecoTrueHisto->Draw("hist same");
        RecoBkgHisto->Draw("hist same");
        leg->Draw();

        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/"+PlotNames[i]+".png");
        delete PlotCanvas;

        if (i == (Vars.size() - 1)) {
            // We can get efficienty and purity data with last histograms
            TH1D* TrueHisto = TrueSignals.ToTH1(TargetPOT);
            double TrueEvents = TrueHisto->Integral();
            double RecoEvents = RecoHisto->Integral();
            double BkgEvents  = RecoBkgHisto->Integral();

            std::cout << std::endl;
            std::cout << "True events: " << TrueEvents << std::endl;
            std::cout << "Reconstructed events: " << RecoEvents << std::endl;
            std::cout << "Background events: " << BkgEvents << std::endl;
            std::cout << "Efficiency (reco / total): " << RecoEvents / TrueEvents << std::endl;
            std::cout << "Purity (1 - bkg / reco): " << 1 - (BkgEvents / RecoEvents) << std::endl;
        }
    }
}