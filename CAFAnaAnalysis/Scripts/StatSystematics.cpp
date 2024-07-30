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

using namespace std;
using namespace Constants;

// https://github.com/SBNSoftware/sbnana/blob/245aecf3f422e54820ffcdc53db326e5cd8e8fe0/sbnana/CAFAna/Core/EnsembleSpectrum.cxx#L196C3-L224C4
void DrawErrorBand(TH1* nom, TGraphAsymmErrors* band, int bandCol, double alpha) {
    if(bandCol == -1) bandCol = nom->GetLineColor()-10; // hopefully a lighter version

    // Check if this pad has already been drawn in
    const bool same = gPad && !gPad->GetListOfPrimitives()->IsEmpty();

    nom->Draw(same ? "hist same" : "hist");

    band->SetFillColorAlpha(bandCol, alpha);
    band->Draw("e2 same");

    nom->Draw("hist same");

    // If we are the first spectrum to draw, scale the y-axis appropriately to
    // fit the error band as well as the nominal
    if(!same){
      double maxY = 0;
      // Don't consider underflow or overflow bins when determining maximum
      for(int i = 1; i < band->GetN()-1; ++i){
        maxY = std::max(maxY, band->GetY()[i] + band->GetErrorYhigh(i));
      }

      // Use non-zero lower value so that SetLogy() still works
      nom->GetYaxis()->SetRangeUser(1e-10, 1.1 * maxY);
    }
    gPad->RedrawAxis();
}

// https://github.com/SBNSoftware/sbnana/blob/245aecf3f422e54820ffcdc53db326e5cd8e8fe0/sbnana/CAFAna/Core/Utilities.cxx#L974-L1000
double FindQuantile(double frac, std::vector<double>& xs) {
    // This turns out to be a much more fraught issue than you would naively
    // expect. This algorithm is equivalent to R-6 here:
    // https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample

    // In principle we could use std::nth_element(). Probably doesn't matter
    // much in practice since this is only for plotting.
    std::sort(xs.begin(), xs.end());

    const int N = xs.size();
    // The index we would ideally be sampling at
    const double h = frac*(N+1);
    // The indices on either side where we have to actually evaluate
    const unsigned int h0 = std::floor(h);
    const unsigned int h1 = std::ceil(h);
    if(h0 == 0) return xs[0]; // Don't underflow indexing
    if(h1 > xs.size()) return xs.back(); // Don't overflow indexing
    // The values at those indices
    const double x0 = xs[h0-1]; // wikipedia is using 1-based indexing
    const double x1 = xs[h1-1];

    if(h0 == h1) return x0;

    // Linear interpolation
    return (h1-h)*x0 + (h-h0)*x1;
}

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
        int max = RecoHist->GetXaxis()->GetXmax();
        int min = RecoHist->GetXaxis()->GetXmin();

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
            const double recoy0 = FindQuantile(.5-0.6827/2, recoys);
            const double recoy1 = FindQuantile(.5+0.6827/2, recoys);

            std::vector<double> bkgys;
            bkgys.push_back(BkgHist->GetBinContent(binIdx)+TMath::Sqrt(BkgHist->GetBinContent(binIdx)));
            const double bkgy0 = FindQuantile(.5-0.6827/2, bkgys);
            const double bkgy1 = FindQuantile(.5+0.6827/2, bkgys);

            std::vector<double> recotrueys;
            recotrueys.push_back(RecoTrueHist->GetBinContent(binIdx)+TMath::Sqrt(RecoTrueHist->GetBinContent(binIdx)));
            const double recotruey0 = FindQuantile(.5-0.6827/2, recotrueys);
            const double recotruey1 = FindQuantile(.5+0.6827/2, recotrueys);

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
        DrawErrorBand(RecoHist, RecoErrorBand, -1, 1);
        RecoTrueHist->Draw("hist same");
        DrawErrorBand(RecoTrueHist, RecoTrueErrorBand, -1, 1);
        BkgHist->Draw("hist same");
        DrawErrorBand(BkgHist, BkgErrorBand, -1, 1);
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
        double RecoCorrMin = RecoCorrMatrix->GetMinimum();
        double RecoCorrMax = RecoCorrMatrix->GetMaximum();
        RecoCorrMatrix->GetZaxis()->SetRangeUser(RecoCorrMin,RecoCorrMax);
        RecoCorrMatrix->GetZaxis()->CenterTitle();
        RecoCorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        RecoCorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        RecoCorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        RecoCorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        RecoCorrMatrix->GetZaxis()->SetNdivisions(5);

        RecoCorrMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        RecoCorrMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        RecoCorrMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/CorrReco"+PlotNames[iVar]+".png");

        // Plot bkg corr matrix
        double BkgCorrMin = BkgCorrMatrix->GetMinimum();
        double BkgCorrMax = BkgCorrMatrix->GetMaximum();
        BkgCorrMatrix->GetZaxis()->SetRangeUser(BkgCorrMin,BkgCorrMax);
        BkgCorrMatrix->GetZaxis()->CenterTitle();
        BkgCorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        BkgCorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        BkgCorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        BkgCorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        BkgCorrMatrix->GetZaxis()->SetNdivisions(5);

        BkgCorrMatrix->GetXaxis()->SetTitle(AxisTitle.c_str());
        BkgCorrMatrix->GetYaxis()->SetTitle(AxisTitle.c_str());

        PlotCanvas->cd();
        BkgCorrMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/Statistical/CorrBkg"+PlotNames[iVar]+".png");

        SaveFile->WriteObject(RecoCovMatrix, PlotNames[iVar]+"_cov");
        SaveFile->WriteObject(BkgCovMatrix, PlotNames[iVar]+"bkg_cov");
    }
}