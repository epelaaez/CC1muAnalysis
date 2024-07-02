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

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionMigrationMatrix() {
    // Some useful variables for later.
    const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const double TargetPOT(6.6e20);

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
    TString RootFilePath = "/pnfs/sbnd/persistent/users/epelaez/CAFAnaOutput/SelectionMigrationMatrix.root";
    TFile* SaveFile = new TFile(RootFilePath, "RECREATE");

    // Vectors to fill with variable pairs and information to plot
    std::vector<std::pair<Var, Var>> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

    // Muon angle
    Vars.push_back({kMuonCosTheta, kRecoTruthMuonCosTheta}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("MuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

    // Leading proton angle
    Vars.push_back({kLeadingProtonCosTheta, kRecoTruthLeadingProtonCosTheta}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("LeadingProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L}})");

    // Recoil proton angle
    Vars.push_back({kRecoilProtonCosTheta, kRecoTruthRecoilProtonCosTheta}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("RecoilProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{R}})");

    // Opening angle between protons
    Vars.push_back({kCosOpeningAngleProtons, kRecoTruthCosOpeningAngleProtons}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleProtons"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");

    // Opening angle between muon and total proton
    Vars.push_back({kCosOpeningAngleMuonTotalProton, kRecoTruthCosOpeningAngleMuonTotalProton}); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleMuonTotalProton"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");

    // Delta alpha transverse
    Vars.push_back({kDeltaAlphaT, kRecoTruthDeltaAlphaT}); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("DeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Transverse momentum
    Vars.push_back({kTransverseMomentum, kRecoTruthTransverseMomentum}); VarBins.push_back(bTransverseMomentumBins);
    PlotNames.push_back("TransverseMomentum"); VarLabels.push_back("#delta P_{T}");

    // Muon momentum 
    Vars.push_back({kMuonMomentum, kRecoTruthMuonMomentum}); VarBins.push_back(bMuonMomentumBins);
    PlotNames.push_back("MuonMomentum"); VarLabels.push_back("|#vec{p}_{#mu}|");

    // Leading proton momentum 
    Vars.push_back({kLeadingProtonMomentum, kRecoTruthLeadingProtonMomentum}); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("LeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back({kRecoilProtonMomentum, kRecoTruthRecoilProtonMomentum}); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("RecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    // Construct spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>, 
        std::vector<std::unique_ptr<Spectrum>>
    >> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        const std::vector<double>& BinEdges = VarBins.at(i).Edges();
        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;
        Var kCurrentVar = Vars.at(i).second;

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
            auto RecoBinValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i).first, kNoSpillCut, TempCut);
            InnerSpectra.push_back(std::move(RecoBinValues));
        }
        auto TruthValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i).second, kNoSpillCut, kRecoIsTrueReco);
        Spectra.push_back({std::move(TruthValues), std::move(InnerSpectra)});
    }

    NuLoader.Go();

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        double RecoTotalEvents = 0.;
        auto& [TruthValues, InnerSpectra] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1* TruthValuesHist = TruthValues->ToTH1(TargetPOT);
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

        std::vector<TH1*> RecoValuesHistos;
        for (int j = 0; j < VarBins.at(i).NBins(); j++) {
            RecoValuesHistos.push_back(InnerSpectra.at(j)->ToTH1(TargetPOT));
        }

        // Debugging
        std::cout << PlotNames.at(i) << std::endl;

        for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
            double TruthCounts = TruthValuesHist->GetBinContent(x);
            TH1* RecoValuesHist = RecoValuesHistos.at(x - 1); // -1 because ROOT lables bins starting from 1
            for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                double RecoCounts = RecoValuesHist->GetBinContent(y);
                double Ratio = RecoCounts / TruthCounts;
                if (TruthCounts == 0.) Ratio = 0.0;

                // Debugging
                std::cout << "Reco low bin: " << y << ". Counts: " << RecoCounts << std::endl;
                std::cout << "True low bin: " << x << ". Counts: " << TruthCounts << std::endl;
                std::cout << "Ratio: " << Ratio << std::endl;
                std::cout << std::endl;

                MigrationMatrix->Fill(
                    TruthValuesHist->GetXaxis()->GetBinCenter(x),
                    RecoValuesHist->GetXaxis()->GetBinCenter(y),
                    Ratio
                );
            }
            RecoTotalEvents += RecoValuesHist->Integral();
        }

        MigrationMatrix->GetXaxis()->SetTitle(("True " + VarLabels.at(i)).c_str());
        MigrationMatrix->GetYaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        MigrationMatrix->Draw("colz text");

        // Sanity check, these should be the same
        std::cout << "Total true var events: " << TruthValuesHist->Integral() << std::endl;
        std::cout << "Total reco var events: " << RecoTotalEvents << std::endl;

        // Save as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Migration/Migration"+PlotNames[i]+".png");
        
        // Save to root file
        SaveFile->WriteObject(MigrationMatrix, PlotNames[i]+"_migration");

        delete PlotCanvas;
        delete MigrationMatrix;
    }
    // Close file
    SaveFile->Close();
}