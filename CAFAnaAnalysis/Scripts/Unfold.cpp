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

    // Plots to unfold
    std::vector<TString> PlotNames; std::vector<TString> XLabels;

    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Muon angle
    PlotNames.push_back("MuonCosTheta"); XLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

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
        TMatrixD ResponseMatrix(m, n); H2M(ResponseHist, ResponseMatrix, kTRUE); // X axis: True, Y axis: Reco
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
        TH1D* TrueUnf = new TH1D("TrueUnf_"+PlotNames[iPlot],";"+XLabels[iPlot]+";Events",n,Nuedges);
        TVectorD AcTrueUnfold = AddSmear * SignalVector;
        V2H(AcTrueUnfold, TrueUnf); ReweightXSec(TrueUnf);

        // Declare canvas and legend
        TCanvas* PlotCanvas = new TCanvas(PlotNames[iPlot],PlotNames[iPlot],205,34,1124,768);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        TrueUnf->Draw("hist");

        // Save histogram
        TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Unfolded/"+PlotNames[iPlot]+".png");
        delete PlotCanvas;
    }
}