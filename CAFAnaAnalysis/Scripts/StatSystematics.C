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

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace Constants;

void StatSystematics() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;

    // Get integrated flux
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root");
    TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
    double IntegratedFlux = (HistoFlux->Integral() * TargetPOT / POTPerSpill / Nominal_UB_XY_Surface);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/Uncertainties/Statistical");

    // Load root file with histograms
    TString HistoFile = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> File(TFile::Open(HistoFile));

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsStats.root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    // Plot names
    std::vector<TString> PlotNames;
    PlotNames.push_back("MuonCosTheta");
    PlotNames.push_back("LeadingProtonCosTheta");
    PlotNames.push_back("RecoilProtonCosTheta");
    PlotNames.push_back("CosOpeningAngleProtons");
    PlotNames.push_back("CosOpeningAngleMuonTotalProton");
    PlotNames.push_back("DeltaAlphaT");
    PlotNames.push_back("TransverseMomentum");
    PlotNames.push_back("MuonMomentum");
    PlotNames.push_back("LeadingProtonMomentum");
    PlotNames.push_back("RecoilProtonMomentum");
    PlotNames.push_back("SerialTransverseMomentum_InMuonCosTheta");
    PlotNames.push_back("SerialDeltaAlphaT_InMuonCosTheta");
    PlotNames.push_back("SerialCosOpeningAngleProtons_InMuonCosTheta");
    PlotNames.push_back("SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta");

    const int NVars = PlotNames.size();

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
    
    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.15);
    PlotCanvas->SetRightMargin(0.15);
    PlotCanvas->SetBottomMargin(0.16);

    for (int iVar = 0; iVar < NVars; iVar++) {
        // Load true plots
        TH1D* RecoHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_reco"));
        TH1D* BkgHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_bkg"));

        int n = RecoHist->GetXaxis()->GetNbins();
        int max = RecoHist->GetXaxis()->GetXmax();
        int min = RecoHist->GetXaxis()->GetXmin();

        // Create covariance matrices
        std::string RecoCovName = "CovStatReco";
        TH2* RecoCovMatrix = new TH2D(
            (RecoCovName + (std::string)PlotNames[iVar]).c_str(),
            RecoCovName.c_str(),
            n, min, max,
            n, min, max
        );

        std::string BkgCovName = "CovStatBkg";
        TH2* BkgCovMatrix = new TH2D(
            (BkgCovName + (std::string)PlotNames[iVar]).c_str(),
            BkgCovName.c_str(),
            n, min, max,
            n, min, max
        );

        for (int x = 1; x < n + 1; x++) {
            for (int y = 1; y <= x; y++) {
                if (x == y) {
                    double RecoBinContent = RecoHist->GetBinContent(x);
                    double RecoValue = TMath::Sqrt(RecoBinContent) * (Units / (IntegratedFlux * NTargets));
                    RecoCovMatrix->SetBinContent(x, y, RecoValue);
                    RecoCovMatrix->SetBinContent(y, x, RecoValue);

                    double BkgBinContent = BkgHist->GetBinContent(x);
                    double BkgValue = TMath::Sqrt(BkgBinContent) * (Units / (IntegratedFlux * NTargets));
                    BkgCovMatrix->SetBinContent(x, y, BkgValue);
                    BkgCovMatrix->SetBinContent(y, x, BkgValue);
                } else {
                    RecoCovMatrix->SetBinContent(x, y, 0.);
                    RecoCovMatrix->SetBinContent(y, x, 0.);
                    BkgCovMatrix->SetBinContent(x, y, 0.);
                    BkgCovMatrix->SetBinContent(y, x, 0.);
                }
            }
        }

        // Plot reco cov matrix
        double RecoCovMin = RecoCovMatrix->GetMinimum();
        double RecoCovMax = RecoCovMatrix->GetMaximum();
        RecoCovMatrix->GetZaxis()->SetRangeUser(RecoCovMin,RecoCovMax);
        RecoCovMatrix->GetZaxis()->CenterTitle();
        RecoCovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        RecoCovMatrix->GetZaxis()->SetTitleSize(TextSize);
        RecoCovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        RecoCovMatrix->GetZaxis()->SetLabelSize(TextSize);
        RecoCovMatrix->GetZaxis()->SetNdivisions(5);

        std::string RecoAxisTitle = (std::string)RecoHist->GetXaxis()->GetTitle();
        std::string AxisTitle = "bin i " + RecoAxisTitle.substr(5, RecoAxisTitle.size() - 1);
        RecoCovMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        RecoCovMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        RecoCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/CovReco"+PlotNames[iVar]+".png");

        // Plot bkg cov matrix
        double BkgCovMin = BkgCovMatrix->GetMinimum();
        double BkgCovMax = BkgCovMatrix->GetMaximum();
        BkgCovMatrix->GetZaxis()->SetRangeUser(BkgCovMin,BkgCovMax);
        BkgCovMatrix->GetZaxis()->CenterTitle();
        BkgCovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        BkgCovMatrix->GetZaxis()->SetTitleSize(TextSize);
        BkgCovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        BkgCovMatrix->GetZaxis()->SetLabelSize(TextSize);
        BkgCovMatrix->GetZaxis()->SetNdivisions(5);

        BkgCovMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        BkgCovMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        BkgCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/CovBkg"+PlotNames[iVar]+".png");

        SaveFile->WriteObject(RecoCovMatrix, PlotNames[iVar]+"reco_cov");
        SaveFile->WriteObject(BkgCovMatrix, PlotNames[iVar]+"bkg_cov");
    }
}