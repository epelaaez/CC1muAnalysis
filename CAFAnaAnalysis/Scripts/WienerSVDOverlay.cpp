// Root includes.
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TGraphAsymmErrors.h>

// std includes.
#include <vector>

// Helpers includes.
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Util.C"
#include "../../Utils/Constants.h"
#include "../../Utils/Tools.cxx"

using namespace std;
using namespace Constants;

void WienerSVDOverlay() {
    int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;

	// Set defaults
	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(6);

	// Load tools
    Tools tools;

	// Dir to save plots
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/BuildEventGenerators/CC1muAnalysis";

	// File with unfolded spectrum
	TString UnfoldedFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Unfolded.root";
    std::unique_ptr<TFile> UnfoldedFile(TFile::Open(UnfoldedFilePath));

	// File with true signal CAF histogram
	TString MatrixRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Matrix.root";
    std::unique_ptr<TFile> MatrixFile(TFile::Open(MatrixRootFilePath));

	// Directory with alternative generators
	TString AltGeneratorsPath = "/pnfs/sbnd/persistent/users/" + (TString)UserName + "/HighSamples/FlatTree/";
	std::vector<TString> AltGenNames; std::vector<TString> AltGenLabels; std::vector<int> Colors;
    
	AltGenNames.push_back(AltGeneratorsPath+"FlatTreeAnalyzerOutput_GENIE_AR23.root"); AltGenLabels.push_back("GENIE AR23"); Colors.push_back(kBlue+8);
	AltGenNames.push_back(AltGeneratorsPath+"FlatTreeAnalyzerOutput_NuWro.root"); AltGenLabels.push_back("NuWro"); Colors.push_back(kRed+1);
	AltGenNames.push_back(AltGeneratorsPath+"FlatTreeAnalyzerOutput_NEUT.root"); AltGenLabels.push_back("NEUT"); Colors.push_back(kOrange+7);
	AltGenNames.push_back(AltGeneratorsPath+"FlatTreeAnalyzerOutput_GiBUU.root"); AltGenLabels.push_back("GiBUU"); Colors.push_back(kGreen+1);

	int NAltGen = AltGenNames.size();
	std::vector<TFile*> AltGenFiles; AltGenFiles.resize(NAltGen);
	for (int iAltGen = 0; iAltGen < NAltGen; ++iAltGen) {
        AltGenFiles[iAltGen] = new TFile(AltGenNames[iAltGen],"readonly");
    }

	TCanvas* PlotCanvas = new TCanvas("WienerSVDOverlay","WienerSVDOverlay",205,34,1124,768);
	PlotCanvas->SetTopMargin(0.13);
	PlotCanvas->SetLeftMargin(0.17);
	PlotCanvas->SetRightMargin(0.05);
	PlotCanvas->SetBottomMargin(0.16);

	// Loop over variables, start at 1 to skip event count
	for (int iVar = 1; iVar < PlotNames.size(); ++iVar) {
		/////////////////////////
		// Get unfolded spectrum
		/////////////////////////

		TH1D* UnfSpectrum = (TH1D*) UnfoldedFile->Get(PlotNames[iVar]+"_unf_spectrum");
        int n = UnfSpectrum->GetNbinsX();
		double edges[n+1]; for (int i = 0; i < n+1; i++) { edges[i] = UnfSpectrum->GetBinLowEdge(i+1); }

		//////////////////////
		// Get smearing matrix
		//////////////////////

		TH2D* SmearingMatrixHisto = (TH2D*) UnfoldedFile->Get(PlotNames[iVar]+"_smear");
		TMatrixD SmearingMatrix(n, n); H2M(SmearingMatrixHisto, SmearingMatrix, kTRUE);

		////////////////////
		// Get true from CAF
		////////////////////

		TH1D* TrueHisto = (TH1D*) MatrixFile->Get(PlotNames[iVar]+"_true");
		TrueHisto->Scale(Units / (IntegratedFlux * NTargets));
		TVectorD TrueVector(n); H2V(TrueHisto, TrueVector);
		TVectorD SmearedTrueVector = SmearingMatrix * TrueVector;

		TH1D* SmearedTrueHisto = new TH1D("SmearedTrue"+PlotNames[iVar],";"+(TString)VarLabels[iVar]+";"+(TString)YLabels[iVar], n, edges);
		V2H(SmearedTrueVector, SmearedTrueHisto); tools.Reweight(SmearedTrueHisto);

		/////////////////////////////
		// Get alternative generators
		/////////////////////////////

		TString GenPlotName = "True" + PlotNames[iVar] + "Plot";
		std::vector<TH1D*> AltGenHistos; AltGenHistos.resize(NAltGen);

		for (int iAltGen = 0; iAltGen < NAltGen; ++iAltGen) {
			TH1D* GenHisto = (TH1D*)(AltGenFiles[iAltGen]->Get(GenPlotName)); tools.Unweight(GenHisto);

        	TVectorD GenVector(n); H2V(GenHisto, GenVector);
			TVectorD SmearedGenVector = SmearingMatrix * GenVector;

			TH1D* SmearedGenHisto = new TH1D("Smeared"+AltGenNames[iAltGen]+PlotNames[iVar],";"+(TString)VarLabels[iVar]+";"+(TString)YLabels[iVar], n, edges);
			V2H(SmearedGenVector, SmearedGenHisto); tools.Reweight(SmearedGenHisto);

			AltGenHistos[iAltGen] = SmearedGenHisto;
		}

		if (PlotNames[iVar].Contains("Serial")) {
			auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[iVar]+"Plot"];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            int StartIndex = 0;

			// Loop over slices
			for (int iSlice = 0; iSlice < NSlices; ++iSlice) {
				TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
				leg->SetBorderSize(0);
				leg->SetNColumns(3);
				leg->SetTextSize(TextSize*0.8);
				leg->SetTextFont(FontStyle);

				////////////////////
				// Slice information
				////////////////////

				TString SlicePlotName = PlotNames[iVar] + "_" + TString(std::to_string(iSlice));
				double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 
				int SliceNBins = SerialVectorBins.at(iSlice); std::vector<double> SerialSliceBinning;
				for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
					double value = SerialVectorRanges.at(StartIndex + iBin);
					SerialSliceBinning.push_back(value);
				}

				//////////////////////////
				// Get error band for data
				//////////////////////////

				TGraphAsymmErrors* ErrorBand = (TGraphAsymmErrors*)UnfoldedFile->Get(SlicePlotName+"_band");

				////////////////
				// Slice histos
				////////////////

				TH1D* SlicedSmearedTrueHisto = tools.GetHistoBins(
					SmearedTrueHisto,
					SerialVectorLowBin.at(iSlice),
					SerialVectorHighBin.at(iSlice),
					SliceWidth,
					SerialSliceBinning,
					"SmearedSignal"
				);
				
				TH1D* SlicedUnfSpectrum = tools.GetHistoBins(
					UnfSpectrum,
					SerialVectorLowBin.at(iSlice),
					SerialVectorHighBin.at(iSlice),
					SliceWidth,
					SerialSliceBinning,
					"UnfoldedSpectrum"
				);

				std::vector<TH1D*> SlicedAltGenHistos; SlicedAltGenHistos.resize(NAltGen);
				for (int iAltGen = 0; iAltGen < NAltGen; ++iAltGen) {
					SlicedAltGenHistos[iAltGen] = tools.GetHistoBins(
						AltGenHistos[iAltGen],
						SerialVectorLowBin.at(iSlice),
						SerialVectorHighBin.at(iSlice),
						SliceWidth,
						SerialSliceBinning,
						"Unfolded" + AltGenNames[iAltGen]
					);
				}
				
				//////////////////
				// Plot everything
				//////////////////

				SlicedSmearedTrueHisto->GetXaxis()->SetTitleFont(FontStyle);
				SlicedSmearedTrueHisto->GetXaxis()->SetLabelFont(FontStyle);
				SlicedSmearedTrueHisto->GetXaxis()->SetLabelSize(TextSize);
				SlicedSmearedTrueHisto->GetXaxis()->SetTitleSize(TextSize);
				SlicedSmearedTrueHisto->GetXaxis()->SetTitleOffset(1.1);
				SlicedSmearedTrueHisto->GetXaxis()->CenterTitle();

				SlicedSmearedTrueHisto->GetYaxis()->SetTitleFont(FontStyle);
				SlicedSmearedTrueHisto->GetYaxis()->SetLabelFont(FontStyle);
				SlicedSmearedTrueHisto->GetYaxis()->SetLabelSize(TextSize*0.8);
				SlicedSmearedTrueHisto->GetYaxis()->SetTitleSize(TextSize*0.9);
				SlicedSmearedTrueHisto->GetYaxis()->SetNdivisions(6);
				SlicedSmearedTrueHisto->GetYaxis()->SetTitleOffset(1.3);
				SlicedSmearedTrueHisto->GetYaxis()->SetTickSize(0);
				SlicedSmearedTrueHisto->GetYaxis()->CenterTitle();

				double Max = 0;
				for (int i = 1; i < ErrorBand->GetN() - 1; ++i){
					Max = std::max(Max, ErrorBand->GetY()[i] + ErrorBand->GetErrorYhigh(i));
				}
				for (int i = 0; i < NAltGen; ++i) {
					Max = std::max(Max, SlicedAltGenHistos[i]->GetMaximum());
				}
				SlicedSmearedTrueHisto->GetYaxis()->SetRangeUser(0., 1.3 * Max);

				TLegendEntry* legSmearedTrue = leg->AddEntry(SlicedSmearedTrueHisto,"True","l");
				SlicedSmearedTrueHisto->SetLineColor(kBlue+2);
				SlicedSmearedTrueHisto->SetLineWidth(4);
				SlicedSmearedTrueHisto->Draw("hist");

				for (int iAltGen = 0; iAltGen < NAltGen; ++iAltGen) {
					TLegendEntry* legRecoTrue = leg->AddEntry(SlicedAltGenHistos[iAltGen],AltGenLabels[iAltGen],"l");
					SlicedAltGenHistos[iAltGen]->SetLineColor(Colors[iAltGen]);
					SlicedAltGenHistos[iAltGen]->SetLineWidth(4);
					SlicedAltGenHistos[iAltGen]->Draw("hist same");
				}

				TLegendEntry* legUnfSpectrum = leg->AddEntry(SlicedUnfSpectrum,"Unfolded data","ep");
				SlicedUnfSpectrum->SetLineColor(kBlack);
				SlicedUnfSpectrum->SetMarkerColor(kBlack);
				SlicedUnfSpectrum->SetMarkerStyle(20);
				SlicedUnfSpectrum->SetMarkerSize(1.);
				SlicedUnfSpectrum->Draw("e1x0 same");
				ErrorBand->Draw("e1 same");
				
				// Slice label
                TLatex *textSlice = new TLatex();
                TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iVar]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

				leg->Draw();
				PlotCanvas->SaveAs(dir+"/Figs/CAFAna/WienerSVDOverlay/"+SlicePlotName+".png");

				delete leg;

            	StartIndex += (SliceNBins + 1);
			}
		} else {
			TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
			leg->SetBorderSize(0);
			leg->SetNColumns(3);
			leg->SetTextSize(TextSize*0.8);
			leg->SetTextFont(FontStyle);

			//////////////////////////
			// Get error band for data
			//////////////////////////

			TGraphAsymmErrors* ErrorBand = (TGraphAsymmErrors*)UnfoldedFile->Get(PlotNames[iVar]+"_band");

			//////////////////
			// Plot everything
			//////////////////

			SmearedTrueHisto->GetXaxis()->SetTitleFont(FontStyle);
			SmearedTrueHisto->GetXaxis()->SetLabelFont(FontStyle);
			SmearedTrueHisto->GetXaxis()->SetLabelSize(TextSize);
			SmearedTrueHisto->GetXaxis()->SetTitleSize(TextSize);
			SmearedTrueHisto->GetXaxis()->SetTitleOffset(1.1);
			SmearedTrueHisto->GetXaxis()->CenterTitle();

			SmearedTrueHisto->GetYaxis()->SetTitleFont(FontStyle);
			SmearedTrueHisto->GetYaxis()->SetLabelFont(FontStyle);
			SmearedTrueHisto->GetYaxis()->SetLabelSize(TextSize);
			SmearedTrueHisto->GetYaxis()->SetTitleSize(TextSize);
			SmearedTrueHisto->GetYaxis()->SetNdivisions(6);
			SmearedTrueHisto->GetYaxis()->SetTitleOffset(1.2);
			SmearedTrueHisto->GetYaxis()->SetTickSize(0);
			SmearedTrueHisto->GetYaxis()->CenterTitle();

			double Max = 0;
			for (int i = 1; i < ErrorBand->GetN() - 1; ++i){
				Max = std::max(Max, ErrorBand->GetY()[i] + ErrorBand->GetErrorYhigh(i));
			}
			for (int i = 0; i < NAltGen; ++i) {
				Max = std::max(Max, AltGenHistos[i]->GetMaximum());
			}
			SmearedTrueHisto->GetYaxis()->SetRangeUser(0., 1.3 * Max);

			TLegendEntry* legSmearedTrue = leg->AddEntry(SmearedTrueHisto,"True","l");
			SmearedTrueHisto->SetLineColor(kBlue+2);
			SmearedTrueHisto->SetLineWidth(4);
			SmearedTrueHisto->Draw("hist");

			for (int iAltGen = 0; iAltGen < NAltGen; ++iAltGen) {
				TLegendEntry* legRecoTrue = leg->AddEntry(AltGenHistos[iAltGen],AltGenLabels[iAltGen],"l");
				AltGenHistos[iAltGen]->SetLineColor(Colors[iAltGen]);
				AltGenHistos[iAltGen]->SetLineWidth(4);
				AltGenHistos[iAltGen]->Draw("hist same");
			}

			TLegendEntry* legUnfSpectrum = leg->AddEntry(UnfSpectrum,"Unfolded data","ep");
			UnfSpectrum->SetLineColor(kBlack);
			UnfSpectrum->SetMarkerColor(kBlack);
			UnfSpectrum->SetMarkerStyle(20);
			UnfSpectrum->SetMarkerSize(1.);
			UnfSpectrum->Draw("e1x0  same");
			ErrorBand->Draw("e1 same");

			leg->Draw();
			PlotCanvas->SaveAs(dir+"/Figs/CAFAna/WienerSVDOverlay/"+PlotNames[iVar]+".png");

			delete leg;
		}
	}
}
