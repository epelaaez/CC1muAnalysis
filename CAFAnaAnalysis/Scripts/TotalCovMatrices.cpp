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
    gStyle->SetPaintTextFormat("4.1f");

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Plot names
    std::vector<TString> PlotNames;
    PlotNames.push_back("EventCount");
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

    // Vector with all cross section systematic files
    std::vector<std::unique_ptr<TFile>> CovFiles;

    // Add xsec systematics
    for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
        std::string SystName = std::get<0>(XSecSystsVector.at(iSyst));
        TString FilePath =  "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
        std::unique_ptr<TFile> File(TFile::Open(FilePath));
        CovFiles.push_back(std::move(File));
    }

    // Add flux systematics
    for (int iSyst = 0; iSyst < (int) FluxSystsVector.size(); iSyst++) {
        std::string SystName = std::get<0>(FluxSystsVector.at(iSyst));
        TString FilePath =  "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
        std::unique_ptr<TFile> File(TFile::Open(FilePath));
        CovFiles.push_back(std::move(File));
    }

    // Add stat systematics
    std::unique_ptr<TFile> StatsFile(TFile::Open("/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsStats.root"));
    CovFiles.push_back(std::move(StatsFile));

    // Add POT systematics
    std::unique_ptr<TFile> POTFile(TFile::Open("/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsPOT.root"));
    CovFiles.push_back(std::move(POTFile));

    // Add NTargets systematics
    std::unique_ptr<TFile> NTargetsFile(TFile::Open("/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsNTargets.root"));
    CovFiles.push_back(std::move(NTargetsFile));

    // Add Detector systematics
    std::unique_ptr<TFile> DetectorFile(TFile::Open("/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsDetector.root"));
    CovFiles.push_back(std::move(DetectorFile));

    // Add Reinteraction systematics
    std::unique_ptr<TFile> ReinteractionFile(TFile::Open("/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematicsReinteraction.root"));
    CovFiles.push_back(std::move(ReinteractionFile));

    // File with reco histograms
    TString HistoFile = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> File(TFile::Open(HistoFile));

    // File to store total cov matrices
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/TotalCovMatrices.root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    // Get integrated flux
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root");
    TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
    double IntegratedFlux = (HistoFlux->Integral() * TargetPOT / POTPerSpill / Nominal_UB_XY_Surface);

    TCanvas* PlotCanvas = new TCanvas("Cov","Cov",205,34,1124,768);

    // Loop over all vars
    for (int iVar = 0; iVar < (int) PlotNames.size(); iVar++) {
        // Margins for matrices
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        // Get first matrix to get correct dimensions
        TH2D* FirstCovHist = (TH2D*)(CovFiles[0]->Get<TH2D>(PlotNames[iVar]+"_cov"));
        int n = FirstCovHist->GetXaxis()->GetNbins();
        double max = FirstCovHist->GetXaxis()->GetXmax();
        double min = FirstCovHist->GetXaxis()->GetXmin();
        TMatrixD TotalCovMatrix(n, n); H2M(FirstCovHist, TotalCovMatrix, kTRUE);

        // Add all uncertainties
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

        // Get bin by bin uncertainties
        if (PlotNames[iVar] == "EventCount") {
            double Total = 0.;

            // Margins for single bin plot
            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(3);
            leg->SetTextSize(TextSize*0.8);
            leg->SetTextFont(FontStyle);

            // Add xsec uncertainties
            TMatrixD XSecFracCov(1, 1);
            for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
                TH2D* FracCovHist = (TH2D*)(CovFiles[iSyst]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
                TMatrixD FracCovMatrix(1, 1); H2M(FracCovHist, FracCovMatrix, kTRUE);
                XSecFracCov += FracCovMatrix;
            }
            TH1D* XSecHisto = new TH1D("XSec", ";;Uncertainty [%]", 1, 0, 1);
            XSecHisto->SetBinContent(1, TMath::Sqrt(XSecFracCov(0, 0)) * 100);
            Total += TMath::Power(TMath::Sqrt(XSecFracCov(0, 0)) * 100, 2);
            TLegendEntry* legXSec = leg->AddEntry(XSecHisto,"XSec","l");
            XSecHisto->SetLineColor(kBlue+2);
            XSecHisto->SetLineWidth(4);
            XSecHisto->SetMarkerSize(1.5);
            XSecHisto->SetMarkerColor(kBlue+2);

            // Style 
            XSecHisto->GetXaxis()->SetNdivisions(5);
            XSecHisto->GetXaxis()->SetLabelSize(0);
            XSecHisto->GetXaxis()->SetTitleSize(0);

            XSecHisto->GetYaxis()->SetTitleFont(FontStyle);
            XSecHisto->GetYaxis()->SetLabelFont(FontStyle);
            XSecHisto->GetYaxis()->SetLabelSize(TextSize);
            XSecHisto->GetYaxis()->SetTitleSize(TextSize);
            XSecHisto->GetYaxis()->SetNdivisions(6);
            XSecHisto->GetYaxis()->SetTitleOffset(1.);
            XSecHisto->GetYaxis()->SetTickSize(0);
            XSecHisto->GetYaxis()->CenterTitle();
            XSecHisto->GetYaxis()->SetRangeUser(0.,35);
            XSecHisto->Draw("hist text0");

            // Add flux uncertainties
            TMatrixD FluxFracCov(1, 1);
            for (int iSyst = 0; iSyst < (int) FluxSystsVector.size(); iSyst++) {
                TH2D* FracCovHist = (TH2D*)(CovFiles[iSyst + (int)XSecSystsVector.size()]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
                TMatrixD FracCovMatrix(1, 1); H2M(FracCovHist, FracCovMatrix, kTRUE);
                FluxFracCov += FracCovMatrix;
            }
            TH1D* FluxHisto = new TH1D("Flux", ";;Uncertainty [%]", 1, 0, 1);
            FluxHisto->SetBinContent(1, TMath::Sqrt(FluxFracCov(0, 0)) * 100);
            Total +=  TMath::Power(TMath::Sqrt(FluxFracCov(0, 0)) * 100, 2);
            TLegendEntry* legFlux = leg->AddEntry(FluxHisto,"Flux","l");
            FluxHisto->SetLineColor(kRed+1);
            FluxHisto->SetLineWidth(4);
            FluxHisto->SetMarkerSize(1.5);
            FluxHisto->SetMarkerColor(kRed+1);
            FluxHisto->Draw("hist text0 same");

            int offset = XSecSystsVector.size() + FluxSystsVector.size();

            // Add Stat systematics
            TMatrixD StatFracCov(1, 1);
            TH2D* StatFracCovHist = (TH2D*)(CovFiles[offset]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            H2M(StatFracCovHist, StatFracCov, kTRUE);
            TH1D* StatsHisto = new TH1D("Stats", ";;Uncertainty [%]", 1, 0, 1);
            StatsHisto->SetBinContent(1, TMath::Sqrt(StatFracCov(0, 0)) * 100);
            Total += TMath::Power(TMath::Sqrt(StatFracCov(0, 0)) * 100, 2);
            TLegendEntry* legStats= leg->AddEntry(StatsHisto,"Stat","l");
            StatsHisto->SetLineColor(kMagenta+1);
            StatsHisto->SetLineWidth(4);
            StatsHisto->SetMarkerSize(1.5);
            StatsHisto->SetMarkerColor(kMagenta+1);
            StatsHisto->Draw("hist text0 same");

            // Add POT systematics
            TMatrixD POTFracCov(1, 1);
            TH2D* POTFracCovHist = (TH2D*)(CovFiles[offset + 1]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            H2M(POTFracCovHist, POTFracCov, kTRUE);
            TH1D* POTHisto = new TH1D("POT", ";;Uncertainty [%]", 1, 0, 1);
            POTHisto->SetBinContent(1, TMath::Sqrt(POTFracCov(0, 0)) * 100);
            Total += TMath::Power(TMath::Sqrt(POTFracCov(0, 0)) * 100, 2);
            TLegendEntry* legPOT= leg->AddEntry(POTHisto,"POT","l");
            POTHisto->SetLineColor(kGreen);
            POTHisto->SetLineWidth(4);
            POTHisto->SetMarkerSize(1.5);
            POTHisto->SetMarkerColor(kGreen);
            POTHisto->Draw("hist text0 same");

            // Add NTargets systematics
            TMatrixD NTargetsFracCov(1, 1);
            TH2D* NTargetsFracCovHist = (TH2D*)(CovFiles[offset + 2]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            H2M(NTargetsFracCovHist, NTargetsFracCov, kTRUE);
            TH1D* NTargetsHisto = new TH1D("NTargets", ";;Uncertainty [%]", 1, 0, 1);
            NTargetsHisto->SetBinContent(1, TMath::Sqrt(NTargetsFracCov(0, 0)) * 100);
            Total += TMath::Power(TMath::Sqrt(NTargetsFracCov(0, 0)) * 100, 2);
            TLegendEntry* legNTargets= leg->AddEntry(NTargetsHisto,"NTargets","l");
            NTargetsHisto->SetLineColor(kYellow+1);
            NTargetsHisto->SetLineWidth(4);
            NTargetsHisto->SetMarkerSize(1.5);
            NTargetsHisto->SetMarkerColor(kYellow+1);
            NTargetsHisto->Draw("hist text0 same");

            // Add Detector systematics
            TMatrixD DetectorFracCov(1, 1);
            TH2D* DetectorFracCovHist = (TH2D*)(CovFiles[offset + 3]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            H2M(DetectorFracCovHist, DetectorFracCov, kTRUE);
            TH1D* DetectorHisto = new TH1D("Detector", ";;Uncertainty [%]", 1, 0, 1);
            DetectorHisto->SetBinContent(1, TMath::Sqrt(DetectorFracCov(0, 0)) * 100);
            Total += TMath::Power(TMath::Sqrt(DetectorFracCov(0, 0)) * 100, 2);
            TLegendEntry* legDetector= leg->AddEntry(DetectorHisto,"Detector","l");
            DetectorHisto->SetLineColor(kCyan-3);
            DetectorHisto->SetLineWidth(4);
            DetectorHisto->SetMarkerSize(1.5);
            DetectorHisto->SetMarkerColor(kCyan-3);
            DetectorHisto->Draw("hist text0 same");

            // Add Reinteraction systematics
            TMatrixD ReinteractionFracCov(1, 1);
            TH2D* ReinteractionFracCovHist = (TH2D*)(CovFiles[offset + 4]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            H2M(ReinteractionFracCovHist, ReinteractionFracCov, kTRUE);
            TH1D* ReinteractionHisto = new TH1D("Reinteraction", ";;Uncertainty [%]", 1, 0, 1);
            ReinteractionHisto->SetBinContent(1, TMath::Sqrt(ReinteractionFracCov(0, 0)) * 100);
            Total += TMath::Power(TMath::Sqrt(ReinteractionFracCov(0, 0)) * 100, 2);
            TLegendEntry* legReinteraction= leg->AddEntry(ReinteractionHisto,"Reinteraction","l");
            ReinteractionHisto->SetLineColor(kTeal + 3);
            ReinteractionHisto->SetLineWidth(4);
            ReinteractionHisto->SetMarkerSize(1.5);
            ReinteractionHisto->SetMarkerColor(kTeal+3);
            ReinteractionHisto->Draw("hist text0 same");

            // Total 
            TH1D* RecoHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_reco"));
            double BinContentScaled = RecoHist->GetBinContent(1) * (Units / (IntegratedFlux * NTargets));
            double TotalFromCovMatrix = TMath::Sqrt(TotalCovMatrix(0, 0) / (BinContentScaled * BinContentScaled));
            TH1D* TotalHisto = new TH1D("Total", ";;Uncertainty [%]", 1, 0, 1);
            TotalHisto->SetBinContent(1, TotalFromCovMatrix * 100);
            TLegendEntry* legTotal= leg->AddEntry(TotalHisto,"Total","l");
            TotalHisto->SetLineColor(kBlack);
            TotalHisto->SetLineWidth(4);
            TotalHisto->SetMarkerSize(1.5);
            TotalHisto->SetMarkerColor(kBlack);
            TotalHisto->Draw("hist text0 same");

            // These should be the same
            std::cout << TotalFromCovMatrix * 100 << std::endl;
            std::cout << TMath::Sqrt(Total) << std::endl;

            // Save plot
            leg->Draw();
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/SingleBinUnc.png");
        }
    }
}