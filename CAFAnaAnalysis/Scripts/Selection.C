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
using namespace Constants;

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

    // Serial transverse momentum in muon cos theta
    Vars.push_back(kTransverseMomentumInMuonCosTheta); VarBins.push_back(bTransverseMomentumInMuonCosTheta);
    PlotNames.push_back("TransverseMomentumInMuonCosTheta"); VarLabels.push_back("#delta P_{T} (bin #)");

    // Delta alpha transverse in muon cos theta
    Vars.push_back(kDeltaAlphaTInMuonCosTheta); VarBins.push_back(bDeltaAlphaTInMuonCosTheta);
    PlotNames.push_back("DeltaAlphaTInMuonCosTheta"); VarLabels.push_back("#delta #alpha_{T} (bin #)");

    // Opening angle between protons in muon cos theta
    Vars.push_back(kCosOpeningAngleProtonsInMuonCosTheta); VarBins.push_back(bCosOpeningAngleProtonsInMuonCosTheta);
    PlotNames.push_back("CosOpeningAngleProtonsInMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)");
    
    // Opening angle between muon and protons in muon cos theta
    Vars.push_back(kCosOpeningAngleMuonTotalProtonInMuonCosTheta); VarBins.push_back(bCosOpeningAngleMuonTotalProtonInMuonCosTheta);
    PlotNames.push_back("CosOpeningAngleMuonTotalProtonInMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)");

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

    // We now create spectra that will help us get the efficiency and purity data for each of the cuts
    // These spectra are going to use the primary energy as a variable

    // Spectrum with all events
    Spectrum sAllEvents("AllEvents", bPrimaryEnergy, NuLoader, kTrueEnergy, kValidEnergyTruthCut, kNoSpillCut);
    // Spectrum with all reco events
    Spectrum sAllRecoEvents("AllRecoEvents", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kValidEnergyCut);
    // Spectrum with all true signal events
    Spectrum sAllTrueEvents("AllTrueEvents", bPrimaryEnergy, NuLoader, kTrueEnergy, kTruthIsSignalAndEnergy, kNoSpillCut);
    // Spectrum with all true signal events that were reconstructed
    Spectrum sAllTrueRecoEvents("AllTrueRecoEvents", bPrimaryEnergy, NuLoader, kTrueEnergy, kTruthIsSignalAndEnergy, kNoSpillCut, kNoCut);
    // Spectrum with first cut (cosmic)
    Spectrum sFirstCut("FirstCut", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kFirstCut);
    Spectrum sFirstCutTrue("FirstCutTrue", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kFirstCutTrue);
    // Spectrum with second cut (cosmic and vertex FV)
    Spectrum sSecondCut("SecondCut", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kSecondCut);
    Spectrum sSecondCutTrue("SecondCutTrue", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kSecondCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, and one muon)
    Spectrum sThirdCut("ThirdCut", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kThirdCut);
    Spectrum sThirdCutTrue("ThirdCutTrue", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kThirdCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, one muon, and two protons)
    Spectrum sFourthCut("FourthCut", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kFourthCut);
    Spectrum sFourthCutTrue("FourthCutTrue", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kFourthCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, and no charged pions)
    Spectrum sFifthCut("FifthCut", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kFifthCut);
    Spectrum sFifthCutTrue("FifthCutTrue", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kFifthCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, no charged pions, and no neutral pions)
    Spectrum sSixthCut("SixthCut", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kSixthCut);
    Spectrum sSixthCutTrue("SixthCutTrue", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kSixthCutTrue);
    // Spectrum with overall signal definition to sanity check it matches
    Spectrum sRecoSignal("RecoSignal", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kRecoIsSignal); 
    Spectrum sRecoTrueSignal("RecoTrueSignal", bPrimaryEnergy, NuLoader, kPrimaryEnergy, kNoSpillCut, kRecoIsTrueReco); 

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
    }

    // Get histograms for all cuts
    TH1D* AllEventsHisto = sAllEvents.ToTH1(TargetPOT);
    TH1D* AllRecoEventsHisto = sAllRecoEvents.ToTH1(TargetPOT);
    TH1D* AllTrueEventsHisto = sAllTrueEvents.ToTH1(TargetPOT);
    TH1D* AllTrueRecoEventsHisto = sAllTrueRecoEvents.ToTH1(TargetPOT);

    TH1D* FirstCutHisto = sFirstCut.ToTH1(TargetPOT);
    TH1D* FirstCutTrueHisto = sFirstCutTrue.ToTH1(TargetPOT);

    TH1D* SecondCutHisto = sSecondCut.ToTH1(TargetPOT);
    TH1D* SecondCutTrueHisto = sSecondCutTrue.ToTH1(TargetPOT);

    TH1D* ThirdCutHisto = sThirdCut.ToTH1(TargetPOT);
    TH1D* ThirdCutTrueHisto = sThirdCutTrue.ToTH1(TargetPOT);

    TH1D* FourthCutHisto = sFourthCut.ToTH1(TargetPOT);
    TH1D* FourthCutTrueHisto = sFourthCutTrue.ToTH1(TargetPOT);

    TH1D* FifthCutHisto = sFifthCut.ToTH1(TargetPOT);
    TH1D* FifthCutTrueHisto = sFifthCutTrue.ToTH1(TargetPOT);

    TH1D* SixthCutHisto = sSixthCut.ToTH1(TargetPOT);
    TH1D* SixthCutTrueHisto = sSixthCutTrue.ToTH1(TargetPOT);

    TH1D* RecoSignalHisto = sRecoSignal.ToTH1(TargetPOT);
    TH1D* RecoTrueSignalHisto = sRecoTrueSignal.ToTH1(TargetPOT);

    // Get integrals for all cuts
    double AllEventsInt = AllEventsHisto->Integral("width");
    double AllRecoEventsInt = AllRecoEventsHisto->Integral("width");
    double AllTrueEventsInt = AllTrueEventsHisto->Integral("width");
    double AllTrueRecoEventsInt = AllTrueRecoEventsHisto->Integral("width");

    double FirstCutInt = FirstCutHisto->Integral("width");
    double FirstCutTrueInt = FirstCutTrueHisto->Integral("width");

    double SecondCutInt = SecondCutHisto->Integral("width");
    double SecondCutTrueInt = SecondCutTrueHisto->Integral("width");

    double ThirdCutInt = ThirdCutHisto->Integral("width");
    double ThirdCutTrueInt = ThirdCutTrueHisto->Integral("width");

    double FourthCutInt = FourthCutHisto->Integral("width");
    double FourthCutTrueInt = FourthCutTrueHisto->Integral("width");

    double FifthCutInt = FifthCutHisto->Integral("width");
    double FifthCutTrueInt = FifthCutTrueHisto->Integral("width");

    double SixthCutInt = SixthCutHisto->Integral("width");
    double SixthCutTrueInt = SixthCutTrueHisto->Integral("width");

    double RecoSignalInt = RecoSignalHisto->Integral("width");
    double RecoTrueSignalInt = RecoTrueSignalHisto->Integral("width");

    // Print results
    std::cout << "================================" << std::endl;
    std::cout << "All events: " << AllEventsInt << std::endl;
    std::cout << "Reconstructed events: " << AllRecoEventsInt << std::endl;
    std::cout << "True signal events: " << AllTrueEventsInt << std::endl;
    std::cout << "True signal events that were reconstructed: " << AllTrueRecoEventsInt << std::endl;
    std::cout << "Cuts: " << std::endl;
    std::cout << "    Cosmic cut: " << FirstCutInt << ". G.E: " <<  (FirstCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (FirstCutTrueInt / AllTrueRecoEventsInt) * 100. << ". Purity: " << (FirstCutTrueInt / FirstCutInt) * 100. << std::endl;
    std::cout << "    Vertex in FV cut: " << SecondCutInt << ". G.E: " <<  (SecondCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (SecondCutTrueInt / AllTrueRecoEventsInt) * 100. << ". Purity: " << (SecondCutTrueInt / SecondCutInt) * 100. << std::endl;
    std::cout << "    One muon cut: " << ThirdCutInt << ". G.E: " <<  (ThirdCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (ThirdCutTrueInt / AllTrueRecoEventsInt) * 100. << ". Purity: " << (ThirdCutTrueInt / ThirdCutInt) * 100. << std::endl;
    std::cout << "    Two protons cut: " << FourthCutInt << ". G.E: " <<  (FourthCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (FourthCutTrueInt / AllTrueRecoEventsInt) * 100. << ". Purity: " << (FourthCutTrueInt / FourthCutInt) * 100. << std::endl;
    std::cout << "    No charged pions cut: " << FifthCutInt << ". G.E: " <<  (FifthCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (FifthCutTrueInt / AllTrueRecoEventsInt) * 100. << ". Purity: " << (FifthCutTrueInt / FifthCutInt) * 100. << std::endl;
    std::cout << "    No neutral pions cut: " << SixthCutInt << ". G.E: " <<  (SixthCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (SixthCutTrueInt / AllTrueRecoEventsInt) * 100. << ". Purity: " << (SixthCutTrueInt / SixthCutInt) * 100. << std::endl;
    std::cout << "Reconstructed events satisfying signal definition: " << RecoSignalInt << ". Purity: " << (RecoTrueSignalInt / RecoSignalInt) * 100. << std::endl;
}