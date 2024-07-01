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

    // Create the binning schemes for the Vars we wish to plot.
    const Binning bPrimaryEnergy = Binning::Simple(1, 0, 3.0); // one bin
    const Binning bAngleBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bDeltaAlphaBins = Binning::Simple(20, 0.0, 180.0);
    const Binning bTransverseMomentumBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bMuonMomentumBins = Binning::Simple(20, 0.1, 1.2);
    const Binning bProtonMomentumBins = Binning::Simple(20, 0.3, 1.0);

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

    // Delta alpha transverse
    Vars.push_back({kDeltaAlphaT, kRecoTruthDeltaAlphaT}); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("DeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Construct spectra
    std::vector<std::tuple<std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto RecoValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i).first, kNoSpillCut, kRecoIsTrueReco);
        auto TruthValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i).second, kNoSpillCut, kRecoIsTrueReco);
        Spectra.push_back({std::move(RecoValues), std::move(TruthValues)});
    }

    NuLoader.Go();

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto& [RecoValues, TruthValues] = Spectra.at(i);


        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1* RecoValuesHist = RecoValues->ToTH1(TargetPOT);
        TH1* TruthValuesHist = TruthValues->ToTH1(TargetPOT);
        TH2* MigrationMatrix = new TH2F(
            "Migration",
            "Migration",
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max()
        );

        for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
            double TruthCounts = TruthValuesHist->GetBinContent(x);
            for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                double RecoCounts = RecoValuesHist->GetBinContent(y);
                double Ratio = RecoCounts / TruthCounts;

                std::cout << RecoValuesHist->GetXaxis()->GetBinLowEdge(y) << "  " << RecoCounts << std::endl;
                std::cout << TruthValuesHist->GetXaxis()->GetBinLowEdge(x) << "  " << TruthCounts << std::endl;
                std::cout << std::endl;

                if (TruthCounts == 0.) Ratio = 0.;
                MigrationMatrix->Fill(
                    TruthValuesHist->GetXaxis()->GetBinLowEdge(x),
                    RecoValuesHist->GetXaxis()->GetBinLowEdge(y),
                    Ratio
                );
            }
        }

        MigrationMatrix->GetXaxis()->SetTitle(("True " + VarLabels.at(i)).c_str());
        MigrationMatrix->GetYaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        MigrationMatrix->Draw("colz text");

        // Save as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Migration/Migration"+PlotNames[i]+".png");

        PlotCanvas->cd();

        TLegend* leg = new TLegend(0.2,0.73,0.45,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoValuesHist,"Reco","l");
        RecoValuesHist->SetLineColor(kBlue+2);
        RecoValuesHist->SetLineWidth(4);

        TLegendEntry* legTrue = leg->AddEntry(TruthValuesHist,"True","l");
        TruthValuesHist->SetLineColor(kRed+1);
        TruthValuesHist->SetLineWidth(4);

        double TrueMax = TruthValuesHist->GetMaximum();
        double RecoMax = RecoValuesHist->GetMaximum();
        double YAxisRange = 1.15 * ((TrueMax < RecoMax) ? RecoMax:  TrueMax);

        TruthValuesHist->GetYaxis()->SetRangeUser(0., YAxisRange);
        RecoValuesHist->GetYaxis()->SetRangeUser(0., YAxisRange);

        TruthValuesHist->Draw("hist");
        RecoValuesHist->Draw("hist same");
        leg->Draw();

        // Save as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Migration/TrueAndReco"+PlotNames[i]+".png");

        // Save to root file
        SaveFile->WriteObject(MigrationMatrix, PlotNames[i]+"_migration");

        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}