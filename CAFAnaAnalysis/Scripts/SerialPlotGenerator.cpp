// ROOT includes.
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>

// std includes.
#include <vector>

// Generator analysis includes.
#include "../../GeneratorAnalysis/Utils/Tools.cxx"
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace Constants;

void SerialPlotGenerator() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;	

    Tools tools;

    // Load root file with histograms
    TString RootFilePath = "/pnfs/sbnd/persistent/users/epelaez/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));

    // Double differential plots to deserialize
    std::vector<TString> PlotNames; std::vector<TString> XAxisLabel; std::vector<TString> YAxisLabel;

    PlotNames.push_back("SerialTransverseMomentum_InMuonCosTheta");
    XAxisLabel.push_back("#delta P_{T}");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("SerialDeltaAlphaT_InMuonCosTheta");
    XAxisLabel.push_back("#delta #alpha_{T}");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("SerialCosOpeningAngleProtons_InMuonCosTheta");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");


    const int NPlots = PlotNames.size();
    for (int iPlot = 0; iPlot < NPlots; iPlot++) {
        // Load histogram with serial plot
        TH1D* TrueHisto(File->Get<TH1D>(PlotNames.at(iPlot)));
        
        // Flatten out double differential plots
        auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[iPlot]+"Plot"];
        auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
        int StartIndex = 0;

        // Create vector to store deserialize plots
        std::vector<TH1D*> Histos;
        Histos.resize(NSlices);

        // Loop over slices
        for (int iSlice = 0; iSlice < NSlices; iSlice++) {
            TString SlicePlotName = PlotNames[iPlot] + "_" + TString(std::to_string(iSlice));
            double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 

            // Get number of bins
            int SliceNBins = SerialVectorBins.at(iSlice);
            std::vector<double> SerialSliceBinning;

            for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                double value = SerialVectorRanges.at(StartIndex + iBin);
                SerialSliceBinning.push_back(value);
            } // End of the number of bins and the bin ranges declaration

            // Declare canvas and legend
            TString CanvasName = "Canvas_" + SlicePlotName;
            TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1124,768);

            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            TLegend* leg = new TLegend(0.2,0.7,0.75,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(2);
            leg->SetTextSize(TextSize*0.8);
            leg->SetTextFont(FontStyle);
            
            Histos[iSlice]= tools.GetHistoBins(
                    TrueHisto,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "SBND"
                );
            Histos[iSlice]->SetLineWidth(4);
            Histos[iSlice]->SetLineColor(602); // blue

            Histos[iSlice]->GetXaxis()->SetTitleFont(FontStyle);
            Histos[iSlice]->GetXaxis()->SetLabelFont(FontStyle);
            Histos[iSlice]->GetXaxis()->SetNdivisions(8);
            Histos[iSlice]->GetXaxis()->SetLabelSize(TextSize);
            Histos[iSlice]->GetXaxis()->SetTitle(XAxisLabel.at(iPlot));
            Histos[iSlice]->GetXaxis()->SetTitleSize(TextSize);
            Histos[iSlice]->GetXaxis()->SetTitleOffset(1.1);
            Histos[iSlice]->GetXaxis()->CenterTitle();

            Histos[iSlice]->GetYaxis()->SetTitleFont(FontStyle);
            Histos[iSlice]->GetYaxis()->SetLabelFont(FontStyle);
            Histos[iSlice]->GetYaxis()->SetNdivisions(6);
            Histos[iSlice]->GetYaxis()->SetLabelSize(TextSize);
            Histos[iSlice]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
            Histos[iSlice]->GetYaxis()->SetTitleSize(TextSize);
            Histos[iSlice]->GetYaxis()->SetTitleOffset(1.3);
            Histos[iSlice]->GetYaxis()->SetTickSize(0);
            Histos[iSlice]->GetYaxis()->CenterTitle();

            PlotCanvas->cd();
            Histos[iSlice]->Draw("hist same");
            leg->AddEntry(Histos[iSlice],"SBND","l");

            PlotCanvas->cd();
            leg->Draw();

            TLatex *textSlice = new TLatex();
            textSlice->SetTextFont(FontStyle);
            textSlice->SetTextSize(TextSize);
            TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iPlot]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
            textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

            TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Serial/"+SlicePlotName+".png");
            delete PlotCanvas;

            // Update starting index to move to next slice
            StartIndex += (SliceNBins + 1);
        } // End loop over slices
    } // End loop over plots
}