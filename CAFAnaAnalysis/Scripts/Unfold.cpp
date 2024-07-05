// ROOT includes.
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

// Utils includes.
#include "../../Utils/Tools.cxx"
#include "../../Utils/Util.C"
#include "../../Utils/WienerSVD.C"

using namespace std;
using namespace Constants;

void ReweightXSec(TH1D* h, double SF = 1.) {
	int NBins = h->GetXaxis()->GetNbins();
	for (int i = 0; i < NBins; i++) {
		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
		h->SetBinError(i+1,NewError); 
	}
}

void Unfold() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;	

    // Load root file(s) with histograms
    TString SelectionRootFilePath = "/pnfs/sbnd/persistent/users/epelaez/CAFAnaOutput/Selection.root";
    TString MatrixRootFilePath = "/pnfs/sbnd/persistent/users/epelaez/CAFAnaOutput/Matrix.root";
    std::unique_ptr<TFile> SelectionFile(TFile::Open(SelectionRootFilePath));
    std::unique_ptr<TFile> MatrixFile(TFile::Open(MatrixRootFilePath));

    // Flux file
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); // make sure file is in path
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));

    // Integrated flux
    double IntegratedFlux = (HistoFlux->Integral() * TargetPOT / POTPerSpill / Nominal_UB_XY_Surface);

    // Plots to unfold
    std::vector<TString> PlotNames; std::vector<TString> XLabels; std::vector<TString> YLabels;

    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Muon angle
    PlotNames.push_back("MuonCosTheta"); 
    XLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");
    YLabels.push_back("#frac{dcos(#theta_{#vec{p}_{#mu}})}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TransverseMomentum"); 
    XLabels.push_back("#delta P_{T}");
    YLabels.push_back("#frac{d#delta P_{T}}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
    
    const int NPlots = PlotNames.size();

    for (int iPlot = 0; iPlot < NPlots; iPlot++) {
        // Load necessary plots
        TH2D* ResponseHist = (TH2D*)(MatrixFile->Get<TH2D>(PlotNames[iPlot]+"_response")); // response matrix
        TH1D* TruePlot = (TH1D*)(MatrixFile->Get<TH1D>(PlotNames[iPlot]+"_true")); // all true generated events
        TH1D* RecoPlot = (TH1D*)(SelectionFile->Get<TH1D>(PlotNames[iPlot]+"_reco")); // reco events
        TH1D* BkgPlot = (TH1D*)(SelectionFile->Get<TH1D>(PlotNames[iPlot]+"_bkg")); // bkg events
        RecoPlot->Add(BkgPlot, -1); // subtract background from reco events

        int n = TruePlot->GetNbinsX();
        int m = RecoPlot->GetNbinsX();

        double Nuedges[n+1];
        for (int i = 0; i < n+1; i++) { Nuedges[i] = TruePlot->GetBinLowEdge(i+1); }

        // Create objects to store matrices/vectors from Wiener SVD
        TMatrixD AddSmear(n,n);
        TVectorD WF(n);
        TMatrixD UnfoldCov(n,n);
        TMatrixD CovRotation(n,n);

        // Convert histograms to matrices/vectors
        TVectorD SignalVector(n); H2V(TruePlot, SignalVector);
        TVectorD MeasureVector(m); H2V(RecoPlot, MeasureVector);
        TMatrixD ResponseMatrix(m, n); H2M(ResponseHist, ResponseMatrix, kFALSE);
        TMatrixD CovarianceMatrix(m, m); CovarianceMatrix.UnitMatrix();

        TVectorD unfold = WienerSVD(
            ResponseMatrix,
            SignalVector,
            MeasureVector,
            CovarianceMatrix,
            2,
            0.5,
            AddSmear,
            WF,
            UnfoldCov,
            CovRotation
        );

        // Add smear to signal
        TH1D* UnfoldedSpectrum = new TH1D("Unfolded"+PlotNames[iPlot],";"+XLabels[iPlot]+";"+YLabels[iPlot],n,Nuedges);
        V2H(unfold, UnfoldedSpectrum); 
        
        ReweightXSec(UnfoldedSpectrum);
        UnfoldedSpectrum->Scale(Units / (IntegratedFlux * NTargets));

        ReweightXSec(TruePlot);
        TruePlot->Scale(Units / (IntegratedFlux * NTargets));

        // Declare canvas and legend
        TCanvas* PlotCanvas = new TCanvas(PlotNames[iPlot],PlotNames[iPlot],205,34,1124,768);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TLegend* leg = new TLegend(0.2,0.73,0.55,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legRecoTrue = leg->AddEntry(TruePlot,"True","l");
        TruePlot->SetLineColor(kRed+1);
        TruePlot->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(UnfoldedSpectrum,"Unfolded","l");
        UnfoldedSpectrum->SetLineColor(kOrange+7);
        UnfoldedSpectrum->SetLineWidth(4);

        double imax = TMath::Max(UnfoldedSpectrum->GetMaximum(),TruePlot->GetMaximum());
        double YAxisRange = 1.15*imax;
        UnfoldedSpectrum->GetYaxis()->SetRangeUser(0.,YAxisRange);
        TruePlot->GetYaxis()->SetRangeUser(0.,YAxisRange);	

        PlotCanvas->cd();
        UnfoldedSpectrum->Draw("hist same");
        TruePlot->Draw("hist same");
        leg->Draw();

        // Save histogram
        TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Unfolded/"+PlotNames[iPlot]+".png");
        delete PlotCanvas;
    }
}