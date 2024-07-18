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

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionMigrationMatrix() {
    // Some useful variables for later.
    // const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const std::string TargetFile = "/pnfs/sbnd/persistent/users/apapadop/CAF_Files/*.flat.caf.root";

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // We now create migration matrices for several variables, where these represent the ratio
    // between the events with a reconstructed value in bin i over the total true signal events
    // with true value in bin i

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/Matrix.root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    // Vectors to fill with variable pairs and information to plot
    std::vector<std::tuple<Var, Var, TruthVar>> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

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
    Vars.push_back({kLeadingProtonMomentum, kRecoTruthLeadingProtonMomentum, kTruthLeadingProtonMomentum}); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("LeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back({kRecoilProtonMomentum, kRecoTruthRecoilProtonMomentum, kTruthRecoilProtonMomentum}); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("RecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    ////////////////////////////////
    // Double differential variables
    ////////////////////////////////

    // Serial transverse momentum in muon cos theta
    Vars.push_back({kTransverseMomentumInMuonCosTheta, kRecoTruthTransverseMomentumInMuonCosTheta, kTruthTransverseMomentumInMuonCosTheta}); VarBins.push_back(bTransverseMomentumInMuonCosTheta);
    PlotNames.push_back("SerialTransverseMomentum_InMuonCosTheta"); VarLabels.push_back("#delta P_{T} (bin #)");

    // Delta alpha transverse in muon cos theta
    Vars.push_back({kDeltaAlphaTInMuonCosTheta, kRecoTruthDeltaAlphaTInMuonCosTheta, kTruthDeltaAlphaTInMuonCosTheta}); VarBins.push_back(bDeltaAlphaTInMuonCosTheta);
    PlotNames.push_back("SerialDeltaAlphaT_InMuonCosTheta"); VarLabels.push_back("#delta #alpha_{T} (bin #)");

    // Opening angle between protons in muon cos theta
    Vars.push_back({kCosOpeningAngleProtonsInMuonCosTheta, kRecoTruthCosOpeningAngleProtonsInMuonCosTheta, kTruthCosOpeningAngleProtonsInMuonCosTheta}); VarBins.push_back(bCosOpeningAngleProtonsInMuonCosTheta);
    PlotNames.push_back("SerialCosOpeningAngleProtons_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)");
    
    // Opening angle between muon and protons in muon cos theta
    Vars.push_back({kCosOpeningAngleMuonTotalProtonInMuonCosTheta, kRecoTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta, kTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta}); VarBins.push_back(bCosOpeningAngleMuonTotalProtonInMuonCosTheta);
    PlotNames.push_back("SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)");

    // Construct spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>, 
        std::unique_ptr<Spectrum>,
        std::vector<std::unique_ptr<Spectrum>>
    >> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        const std::vector<double>& BinEdges = VarBins.at(i).Edges();
        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;
        Var kCurrentVar = std::get<1>(Vars.at(i));

        for (int j = 0; j < VarBins.at(i).NBins(); j++) {
            double BinMin = BinEdges.at(j);
            double BinMax = (j == VarBins.at(i).NBins() - 1) ?  VarBins.at(i).Max() : BinEdges.at(j + 1);

            const Cut TempCut([=](const caf::SRSliceProxy* slc) {
                return (
                    kRecoIsTrueReco(slc) && 
                    kCurrentVar(slc) >= BinMin &&
                    kCurrentVar(slc) < BinMax
                );
            });
            auto RecoBinValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, TempCut);
            InnerSpectra.push_back(std::move(RecoBinValues));
        }
        auto RecoTruthValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<1>(Vars.at(i)), kNoSpillCut, kRecoIsTrueReco);
        auto TruthValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<2>(Vars.at(i)), kTruthIsSignal, kNoSpillCut);
        Spectra.push_back({std::move(TruthValues), std::move(RecoTruthValues), std::move(InnerSpectra)});
    }

    NuLoader.Go();

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        double RecoTotalEvents = 0.;
        auto& [TruthValues, RecoTruthValues, InnerSpectra] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1* TruthValuesHist = TruthValues->ToTH1(TargetPOT);
        TH1* RecoTruthValuesHist = RecoTruthValues->ToTH1(TargetPOT);

        // Manage under/overflow bins
        TruthValuesHist->SetBinContent(TruthValuesHist->GetNbinsX(), TruthValuesHist->GetBinContent(TruthValuesHist->GetNbinsX()) + TruthValuesHist->GetBinContent(TruthValuesHist->GetNbinsX() + 1));
        RecoTruthValuesHist->SetBinContent(RecoTruthValuesHist->GetNbinsX(), RecoTruthValuesHist->GetBinContent(RecoTruthValuesHist->GetNbinsX()) + RecoTruthValuesHist->GetBinContent(RecoTruthValuesHist->GetNbinsX() + 1));

        TruthValuesHist->SetBinContent(1, TruthValuesHist->GetBinContent(0) + TruthValuesHist->GetBinContent(1));
        RecoTruthValuesHist->SetBinContent(1, RecoTruthValuesHist->GetBinContent(0) + RecoTruthValuesHist->GetBinContent(1));

        TH2* MigrationMatrix = new TH2D(
            "Migration",
            "Migration",
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max()
        );
        TH2* ResponseMatrix = new TH2D(
            "Response",
            "Response",
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max()
        );

        std::vector<TH1*> RecoValuesHistos;
        for (int j = 0; j < VarBins.at(i).NBins(); j++) {
            RecoValuesHistos.push_back(InnerSpectra.at(j)->ToTH1(TargetPOT));
        }

        // Debugging
        std::cout << PlotNames.at(i) << std::endl;

        for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
            double RecoTruthCounts = RecoTruthValuesHist->GetBinContent(x);
            double TruthCounts = TruthValuesHist->GetBinContent(x);
            TH1* RecoValuesHist = RecoValuesHistos.at(x - 1); // -1 because ROOT labels bins starting from 1

            // Manage under/overflow bins
            RecoValuesHist->SetBinContent(RecoValuesHist->GetNbinsX(), RecoValuesHist->GetBinContent(RecoValuesHist->GetNbinsX()) + RecoValuesHist->GetBinContent(RecoValuesHist->GetNbinsX() + 1));
            RecoValuesHist->SetBinContent(1, RecoValuesHist->GetBinContent(0) + RecoValuesHist->GetBinContent(1));

            for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                double RecoCounts = RecoValuesHist->GetBinContent(y);

                // Migration matrix
                double MigrationRatio = RecoCounts / RecoTruthCounts;
                if (RecoTruthCounts == 0.) MigrationRatio = 0.0;

                // Response matrix
                double ResponseRatio = RecoCounts / TruthCounts;
                if (TruthCounts == 0.) ResponseRatio = 0.0;

                // Debugging
                std::cout << "Reco low bin: " << y << ". Counts: " << RecoCounts << std::endl;
                std::cout << "True low bin: " << x << ". Counts: " << TruthCounts << std::endl;
                std::cout << "Reco truth counts: " << RecoTruthCounts << std::endl;
                std::cout << "Migration ratio: " << MigrationRatio << std::endl;
                std::cout << "Response ratio: " << ResponseRatio << std::endl;
                std::cout << std::endl;

                MigrationMatrix->Fill(
                    TruthValuesHist->GetXaxis()->GetBinCenter(x),
                    RecoValuesHist->GetXaxis()->GetBinCenter(y),
                    MigrationRatio
                );
                ResponseMatrix->Fill(
                    TruthValuesHist->GetXaxis()->GetBinCenter(x),
                    RecoValuesHist->GetXaxis()->GetBinCenter(y),
                    ResponseRatio
                );
            }
            RecoTotalEvents += RecoValuesHist->Integral();
        }
        // Sanity check, these should be the same
        std::cout << "Total true var events: " << RecoTruthValuesHist->Integral() << std::endl;
        std::cout << "Total reco var events: " << RecoTotalEvents << std::endl;

        MigrationMatrix->GetXaxis()->SetTitle(("True " + VarLabels.at(i)).c_str());
        MigrationMatrix->GetYaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());
        double MigrationMin = MigrationMatrix->GetMinimum();
        double MigrationMax = MigrationMatrix->GetMaximum();
        MigrationMatrix->GetZaxis()->SetRangeUser(MigrationMin,MigrationMax);
        MigrationMatrix->GetZaxis()->CenterTitle();
        MigrationMatrix->GetZaxis()->SetTitleFont(FontStyle);
        MigrationMatrix->GetZaxis()->SetTitleSize(TextSize);
        MigrationMatrix->GetZaxis()->SetLabelFont(FontStyle);
        MigrationMatrix->GetZaxis()->SetLabelSize(TextSize);
        MigrationMatrix->GetZaxis()->SetNdivisions(5);

        ResponseMatrix->GetXaxis()->SetTitle(("True " + VarLabels.at(i)).c_str());
        ResponseMatrix->GetYaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());
        double ResponseMin = ResponseMatrix->GetMinimum();
        double ResponseMax = ResponseMatrix->GetMaximum();
        ResponseMatrix->GetZaxis()->SetRangeUser(ResponseMin,ResponseMax);
        ResponseMatrix->GetZaxis()->CenterTitle();
        ResponseMatrix->GetZaxis()->SetTitleFont(FontStyle);
        ResponseMatrix->GetZaxis()->SetTitleSize(TextSize);
        ResponseMatrix->GetZaxis()->SetLabelFont(FontStyle);
        ResponseMatrix->GetZaxis()->SetLabelSize(TextSize);
        ResponseMatrix->GetZaxis()->SetNdivisions(5);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        MigrationMatrix->Draw("colz");

        // Save migration matrix as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Matrices/Migration"+PlotNames[i]+".png");

        PlotCanvas->cd();
        ResponseMatrix->Draw("colz");

        // Save response matrix as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Matrices/Response"+PlotNames[i]+".png");
        
        // Save both to root file
        SaveFile->WriteObject(MigrationMatrix, PlotNames[i]+"_migration");
        SaveFile->WriteObject(ResponseMatrix, PlotNames[i]+"_response");
        SaveFile->WriteObject(TruthValuesHist, PlotNames[i]+"_true");

        delete PlotCanvas;
        delete MigrationMatrix;
        delete ResponseMatrix;
    }
    // Close file
    SaveFile->Close();
}
