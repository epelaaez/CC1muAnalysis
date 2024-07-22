// Root includes.
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>

// std includes.
#include <vector>

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

// Generator analysis includes.
#include "../../Utils/Util.C"

using namespace std;
using namespace Constants;

void TotalCovMatrices() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    gStyle->SetOptStat(0);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

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

    // Get cross-section and flux systematics
    std::vector<std::tuple<std::string, int>> SystsVector(XSecSystsVector);
    // SystsVector.insert(SystsVector.end(), FluxSystsVector.begin(), FluxSystsVector.end());
    // Flux systematics not working yet

    // Vector with all cross section systematic files
    std::vector<std::unique_ptr<TFile>> CovFiles;
    for (int iSyst = 0; iSyst < (int) SystsVector.size(); iSyst++) {
        std::string SystName = std::get<0>(SystsVector.at(iSyst));
        TString FilePath =  "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
        std::unique_ptr<TFile> File(TFile::Open(FilePath));
        CovFiles.push_back(std::move(File));
    }

    // File with stat systematics
    std::unique_ptr<TFile> StatsFile(TFile::Open("/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsStats.root"));
    CovFiles.push_back(std::move(StatsFile));

    // File to store total cov matrices
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/TotalCovMatrices.root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    TCanvas* PlotCanvas = new TCanvas("Cov","Cov",205,34,1124,768);
    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.15);
    PlotCanvas->SetRightMargin(0.15);
    PlotCanvas->SetBottomMargin(0.16);

    // Loop over all vars
    for (int iVar = 0; iVar < (int) PlotNames.size(); iVar++) {
        // Get first matrix to get correct dimensions
        TH2D* FirstCovHist = (TH2D*)(CovFiles[0]->Get<TH2D>(PlotNames[iVar]+"_cov"));
        int n = FirstCovHist->GetXaxis()->GetNbins();
        int max = FirstCovHist->GetXaxis()->GetXmax();
        int min = FirstCovHist->GetXaxis()->GetXmin();
        TMatrixD TotalCovMatrix(n, n); H2M(FirstCovHist, TotalCovMatrix, kTRUE);

        // Add cross section uncertainties
        for (int iCov = 1; iCov < (int) CovFiles.size(); iCov++) {
            TH2D* CovHist = (TH2D*)(CovFiles[iCov]->Get<TH2D>(PlotNames[iVar]+"_cov"));
            TMatrixD CovMatrix(n, n); H2M(CovHist, CovMatrix, kTRUE);
            TotalCovMatrix += CovMatrix;
        }
        
        TH2D* TotalCovHist = new TH2D(
            "TotalCov" + PlotNames[iVar],
            "TotalCov" + PlotNames[iVar],
            n, min, max, 
            n, min, max
        );
        // std::cout << TotalCovMatrix.Determinant() << std::endl;
        M2H(TotalCovMatrix, TotalCovHist);

        double CovMin = TotalCovHist->GetMinimum();
        double CovMax = TotalCovHist->GetMaximum();
        TotalCovHist->GetZaxis()->SetRangeUser(CovMin,CovMax);
        TotalCovHist->GetZaxis()->CenterTitle();
        TotalCovHist->GetZaxis()->SetTitleFont(FontStyle);
        TotalCovHist->GetZaxis()->SetTitleSize(TextSize);
        TotalCovHist->GetZaxis()->SetLabelFont(FontStyle);
        TotalCovHist->GetZaxis()->SetLabelSize(TextSize);
        TotalCovHist->GetZaxis()->SetNdivisions(5);
        TotalCovHist->GetXaxis()->SetTitle(FirstCovHist->GetXaxis()->GetTitle());
        TotalCovHist->GetYaxis()->SetTitle(FirstCovHist->GetYaxis()->GetTitle());

        PlotCanvas->cd();
        TotalCovHist->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/TotalCov"+PlotNames[iVar]+".png");

        SaveFile->WriteObject(TotalCovHist, PlotNames[iVar]);
    }
}