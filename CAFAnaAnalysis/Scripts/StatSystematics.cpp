// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"
#include "Helpers.cpp"

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

    const int NVars = PlotNames.size();

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

    for (int iVar = 0; iVar < NVars; iVar++) {
        // Load true plots
        TH1D* RecoHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_reco"));
        TH1D* RecoTrueHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_reco_true"));
        TH1D* BkgHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_bkg"));

        // Legend for plot with error bands
        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoHist,"Reconstructed","l");
        RecoHist->SetLineColor(kBlue+2);
        RecoHist->SetLineWidth(4);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHist,"True","l");
        RecoTrueHist->SetLineColor(kRed+1);
        RecoTrueHist->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(BkgHist,"Background","l");
        BkgHist->SetLineColor(kOrange+7);
        BkgHist->SetLineWidth(4);

        int n = RecoHist->GetXaxis()->GetNbins();
        double max = RecoHist->GetXaxis()->GetXmax();
        double min = RecoHist->GetXaxis()->GetXmin();

        // Get graph with error bands
        TGraphAsymmErrors* RecoErrorBand = new TGraphAsymmErrors;
        TGraphAsymmErrors* RecoTrueErrorBand = new TGraphAsymmErrors;
        TGraphAsymmErrors* BkgErrorBand = new TGraphAsymmErrors;
        for (int binIdx = 0; binIdx < n + 2; ++binIdx) {
            const double recoxnom = RecoHist->GetXaxis()->GetBinCenter(binIdx);
            const double recoynom = RecoHist->GetBinContent(binIdx);
            RecoErrorBand->SetPoint(binIdx, recoxnom, recoynom);

            const double bkgxnom = BkgHist->GetXaxis()->GetBinCenter(binIdx);
            const double bkgynom = BkgHist->GetBinContent(binIdx);
            BkgErrorBand->SetPoint(binIdx, bkgxnom, bkgynom);

            const double recotruexnom = RecoTrueHist->GetXaxis()->GetBinCenter(binIdx);
            const double recotrueynom = RecoTrueHist->GetBinContent(binIdx);
            RecoTrueErrorBand->SetPoint(binIdx, recotruexnom, recotrueynom);

            const double dx = RecoHist->GetXaxis()->GetBinWidth(binIdx);

            std::vector<double> recoys;
            recoys.push_back(RecoHist->GetBinContent(binIdx)+TMath::Sqrt(RecoHist->GetBinContent(binIdx)));
            const double recoy0 = SelectionHelpers::FindQuantile(.5-0.6827/2, recoys);
            const double recoy1 = SelectionHelpers::FindQuantile(.5+0.6827/2, recoys);

            std::vector<double> bkgys;
            bkgys.push_back(BkgHist->GetBinContent(binIdx)+TMath::Sqrt(BkgHist->GetBinContent(binIdx)));
            const double bkgy0 = SelectionHelpers::FindQuantile(.5-0.6827/2, bkgys);
            const double bkgy1 = SelectionHelpers::FindQuantile(.5+0.6827/2, bkgys);

            std::vector<double> recotrueys;
            recotrueys.push_back(RecoTrueHist->GetBinContent(binIdx)+TMath::Sqrt(RecoTrueHist->GetBinContent(binIdx)));
            const double recotruey0 = SelectionHelpers::FindQuantile(.5-0.6827/2, recotrueys);
            const double recotruey1 = SelectionHelpers::FindQuantile(.5+0.6827/2, recotrueys);

            RecoErrorBand->SetPointError(binIdx, dx/2, dx/2, std::max(recoxnom-recoy0, 0.), std::max(recoy1-recoynom, 0.));
            BkgErrorBand->SetPointError(binIdx, dx/2, dx/2, std::max(bkgxnom-bkgy0, 0.), std::max(bkgy1-bkgynom, 0.));
            RecoTrueErrorBand->SetPointError(binIdx, dx/2, dx/2, std::max(recotruexnom-recotruey0, 0.), std::max(recotruey1-recotrueynom, 0.));
        }
        PlotCanvas->cd();
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        double imax = RecoHist->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoHist->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoTrueHist->GetYaxis()->SetRangeUser(0.,YAxisRange);
        BkgHist->GetYaxis()->SetRangeUser(0.,YAxisRange);

        RecoHist->Draw("hist");
        SelectionHelpers::DrawErrorBand(RecoHist, RecoErrorBand, -1, 1);
        RecoTrueHist->Draw("hist same");
        SelectionHelpers::DrawErrorBand(RecoTrueHist, RecoTrueErrorBand, -1, 1);
        BkgHist->Draw("hist same");
        SelectionHelpers::DrawErrorBand(BkgHist, BkgErrorBand, -1, 1);
        leg->Draw();
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/"+PlotNames[iVar]+".png");

        // Set margins for matrices
        PlotCanvas->cd();
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

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

        // Create fractional covariance matrices
        std::string RecoFracCovName = "FracCovStatReco";
        TH2* RecoFracCovMatrix = new TH2D(
            (RecoFracCovName + (std::string)PlotNames[iVar]).c_str(),
            RecoFracCovName.c_str(),
            n, min, max,
            n, min, max
        );

        std::string BkgFracCovName = "FracCovStatBkg";
        TH2* BkgFracCovMatrix = new TH2D(
            (BkgFracCovName + (std::string)PlotNames[iVar]).c_str(),
            BkgFracCovName.c_str(),
            n, min, max,
            n, min, max
        );

        // Create correlation matrices
        std::string RecoCorrName = "CorrStatReco";
        TH2* RecoCorrMatrix = new TH2D(
            (RecoCorrName + (std::string)PlotNames[iVar]).c_str(),
            RecoCorrName.c_str(),
            n, min, max,
            n, min, max
        );

        std::string BkgCorrName = "CorrStatBkg";
        TH2* BkgCorrMatrix = new TH2D(
            (BkgCorrName + (std::string)PlotNames[iVar]).c_str(),
            BkgCorrName.c_str(),
            n, min, max,
            n, min, max
        );

        // Create matrices 
        for (int x = 1; x < n + 1; x++) {
            for (int y = 1; y <= x; y++) {
                if (x == y) {
                    double RecoBinContent = RecoHist->GetBinContent(x);
                    double RecoValue = TMath::Sqrt(RecoBinContent) * (Units / (IntegratedFlux * NTargets));
                    RecoCovMatrix->SetBinContent(x, y, TMath::Max(RecoValue, 1e-8));

                    double RecoFracValue = (RecoBinContent == 0.) ? 0. : RecoValue / (RecoBinContent * RecoBinContent);
                    RecoFracCovMatrix->SetBinContent(x, y, TMath::Max(RecoFracValue, 1e-8));
                    RecoCorrMatrix->SetBinContent(x, y, 1.);

                    double BkgBinContent = BkgHist->GetBinContent(x);
                    double BkgValue = TMath::Sqrt(BkgBinContent) * (Units / (IntegratedFlux * NTargets));
                    BkgCovMatrix->SetBinContent(x, y, TMath::Max(BkgValue, 1e-8));

                    double BkgFracValue = (BkgBinContent == 0.) ? 0. : BkgValue / (BkgBinContent * BkgBinContent);
                    BkgFracCovMatrix->SetBinContent(x, y, TMath::Max(BkgFracValue, 1e-8));
                    BkgCorrMatrix->SetBinContent(x, y, 1.);
                } else {
                    RecoCovMatrix->SetBinContent(x, y, 0.);
                    RecoCovMatrix->SetBinContent(y, x, 0.);
                    BkgCovMatrix->SetBinContent(x, y, 0.);
                    BkgCovMatrix->SetBinContent(y, x, 0.);

                    RecoFracCovMatrix->SetBinContent(x, y, 0.);
                    RecoFracCovMatrix->SetBinContent(y, x, 0.);
                    BkgFracCovMatrix->SetBinContent(x, y, 0.);
                    BkgFracCovMatrix->SetBinContent(y, x, 0.);

                    RecoCorrMatrix->SetBinContent(x, y, 0.);
                    RecoCorrMatrix->SetBinContent(y, x, 0.);
                    BkgCorrMatrix->SetBinContent(x, y, 0.);
                    BkgCorrMatrix->SetBinContent(y, x, 0.);
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

        // Plot reco frac cov matrix
        double RecoFracCovMin = RecoFracCovMatrix->GetMinimum();
        double RecoFracCovMax = RecoFracCovMatrix->GetMaximum();
        RecoFracCovMatrix->GetZaxis()->SetRangeUser(RecoFracCovMin,RecoFracCovMax);
        RecoFracCovMatrix->GetZaxis()->CenterTitle();
        RecoFracCovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        RecoFracCovMatrix->GetZaxis()->SetTitleSize(TextSize);
        RecoFracCovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        RecoFracCovMatrix->GetZaxis()->SetLabelSize(TextSize);
        RecoFracCovMatrix->GetZaxis()->SetNdivisions(5);

        RecoFracCovMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        RecoFracCovMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        RecoFracCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/FracCovReco"+PlotNames[iVar]+".png");

        // Plot bkg frac cov matrix
        double BkgFracCovMin = BkgFracCovMatrix->GetMinimum();
        double BkgFracCovMax = BkgFracCovMatrix->GetMaximum();
        BkgFracCovMatrix->GetZaxis()->SetRangeUser(BkgFracCovMin,BkgFracCovMax);
        BkgFracCovMatrix->GetZaxis()->CenterTitle();
        BkgFracCovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        BkgFracCovMatrix->GetZaxis()->SetTitleSize(TextSize);
        BkgFracCovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        BkgFracCovMatrix->GetZaxis()->SetLabelSize(TextSize);
        BkgFracCovMatrix->GetZaxis()->SetNdivisions(5);

        BkgFracCovMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        BkgFracCovMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        BkgFracCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/FracCovBkg"+PlotNames[iVar]+".png");

        // Plot reco corr matrix
        RecoCorrMatrix->GetZaxis()->SetRangeUser(-1,1);
        RecoCorrMatrix->GetZaxis()->CenterTitle();
        RecoCorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        RecoCorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        RecoCorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        RecoCorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        RecoCorrMatrix->GetZaxis()->SetNdivisions(5);

        RecoCorrMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        RecoCorrMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        RecoCorrMatrix->Draw("colz text");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/CorrReco"+PlotNames[iVar]+".png");

        // Plot bkg corr matrix
        BkgCorrMatrix->GetZaxis()->SetRangeUser(-1,1);
        BkgCorrMatrix->GetZaxis()->CenterTitle();
        BkgCorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        BkgCorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        BkgCorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        BkgCorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        BkgCorrMatrix->GetZaxis()->SetNdivisions(5);

        BkgCorrMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        BkgCorrMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        BkgCorrMatrix->Draw("colz text");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/CorrBkg"+PlotNames[iVar]+".png");

        SaveFile->WriteObject(RecoCovMatrix, PlotNames[iVar]+"_cov");
        SaveFile->WriteObject(BkgCovMatrix, PlotNames[iVar]+"bkg_cov");
    }
}