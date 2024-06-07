#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "Constants.h"

using namespace std;
using namespace Constants;

TString to_string_with_precision(double a_value, const int n);

void GeneratorInteBreakDown() {

    //------------------------------//

    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;			

    TString OutFilePath = "/pnfs/sbnd/persistent/users/epelaez/HighSamples/FlatTree/";

    //------------------------------//

    // Event generators

    std::vector<TString> Names; std::vector<TString> Labels; 
    std::vector<TString> Process{"","QE","MEC","RES","DIS"};
    std::vector<int> Colors{kBlack,kBlue,kRed+1,kOrange+7,kGreen+1,kMagenta+1};

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23.root"); 
    Labels.push_back("GENIE_G18");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_G18.root"); 
    Labels.push_back("GENIE_AR23");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro.root"); 
    Labels.push_back("NuWro");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT.root"); 
    Labels.push_back("NEUT");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU.root"); 
    Labels.push_back("GiBUU");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_NoFSI.root"); 
    Labels.push_back("GiBUU_NoFSI");

    const int NSamples = Names.size();
    const int NColors = Colors.size();
    const int NProcesses = Process.size();

    // Sanity check
    if (NColors < NProcesses) { cout << "Give me some more colors for all the processes!" << endl; return; }

    std::vector<TFile*> Files; Files.resize(NSamples);

    //------------------------------//

    // Plots to overlay

    std::vector<TString> PlotNames;
    std::vector<TString> YAxisLabel;

    //------------------------------//

    // Pre FSI
    PlotNames.push_back("TrueNoFSILeadingProtonCosThetaPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIRecoilProtonCosThetaPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSILeadingProtonMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{L}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIRecoilProtonMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{R}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIMuonMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{#mu}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSICosOpeningAngleProtonsPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSICosOpeningAngleMuonTotalProtonPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSITransverseMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // Post FSI
    PlotNames.push_back("TrueMuonCosThetaPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos#theta_{#mu}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueLeadingProtonCosThetaPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueRecoilProtonCosThetaPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueLeadingProtonMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{L}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueRecoilProtonMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{R}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueMuonMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{#mu}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueCosOpeningAngleProtonsPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueCosOpeningAngleMuonTotalProtonPlot");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueTransverseMomentumPlot");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    //------------------------------//

    const int NPlots = PlotNames.size();
    const int NLabels = YAxisLabel.size();

    // sanity check
    if ( NPlots != NLabels) { cout << "Inconsistent number of plots and labels! Aborting !" << endl; return; }

    //------------------------------//	

    // Loop over the samples to open the files and the TTree

    for (int iSample = 0; iSample < NSamples; iSample++) {

        Files[iSample] = new TFile(Names[iSample],"readonly");

    } // End of the loop over the samples

    //------------------------------//

    // Loop over the plots to be compared

    for (int iPlot = 0; iPlot < NPlots; iPlot++) {

        for (int iSample = 0; iSample < NSamples; iSample++) {	

            TString LabelCopy = Labels[iSample];
            TString CanvasName = "ThreeDKI_"+LabelCopy.ReplaceAll(" ","_")+"_InteBreakDown_" + PlotNames[iPlot];
            TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
            PlotCanvas->cd();
            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);		
            PlotCanvas->Draw();	

            TLegend* leg = new TLegend(0.13,0.88,0.95,0.99);
            leg->SetBorderSize(0);
            leg->SetNColumns(3);
            leg->SetTextSize(TextSize);	
            leg->SetTextFont(FontStyle);						
            leg->SetMargin(0.2);				
            leg->SetFillColor(0);				

            // Loop over the interaction processes

            std::vector<TH1D*> Histos; Histos.resize(NProcesses);

            for (int iprocess = 0; iprocess < NProcesses; iprocess++) {	

                Histos[iprocess] = (TH1D*)(Files[iSample]->Get(Process[iprocess]+PlotNames[iPlot]));

                Histos[iprocess]->SetLineWidth(4);
                Histos[iprocess]->SetLineColor( Colors.at(iprocess) );	

                Histos[iprocess]->GetXaxis()->SetTitleFont(FontStyle);
                Histos[iprocess]->GetXaxis()->SetLabelFont(FontStyle);
                Histos[iprocess]->GetXaxis()->SetNdivisions(8);
                Histos[iprocess]->GetXaxis()->SetLabelSize(TextSize);
                Histos[iprocess]->GetXaxis()->SetTitleSize(TextSize);	
                Histos[iprocess]->GetXaxis()->SetTitleOffset(1.1);					
                Histos[iprocess]->GetXaxis()->CenterTitle();						

                Histos[iprocess]->GetYaxis()->SetTitleFont(FontStyle);
                Histos[iprocess]->GetYaxis()->SetLabelFont(FontStyle);
                Histos[iprocess]->GetYaxis()->SetNdivisions(6);
                Histos[iprocess]->GetYaxis()->SetLabelSize(TextSize);
                // Histos[iprocess]->GetYaxis()->SetTitle("Cross Section [10^{-38} cm^{2}/Ar]");
                Histos[iprocess]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
                Histos[iprocess]->GetYaxis()->SetTitleSize(TextSize);
                Histos[iprocess]->GetYaxis()->SetTitleOffset(1.25);
                Histos[iprocess]->GetYaxis()->SetTickSize(0);
                Histos[iprocess]->GetYaxis()->CenterTitle();	
                Histos[iprocess]->GetYaxis()->SetRangeUser(0.,1.15*Histos[0]->GetMaximum());

                Histos[iprocess]->Draw("hist same");
                Histos[0]->Draw("hist same");	

                double frac = Histos[iprocess]->Integral("width")/Histos[0]->Integral("width") * 100.;
                TString LegLabel = Process[iprocess] + " (" + to_string_with_precision(frac,1) + "%)";
                if (iprocess == 0) { LegLabel = "Total (" + to_string_with_precision(frac,1) + "%)"; }
                TLegendEntry* legColor = leg->AddEntry(Histos[iprocess],LegLabel,"l");
                legColor->SetTextColor( Colors.at(iprocess) ); 

            } // End of the loop over the processes

            PlotCanvas->cd();
            leg->Draw();

            TLatex *textSlice = new TLatex();
            textSlice->SetTextFont(FontStyle);
            textSlice->SetTextSize(TextSize);
            TString PlotNameDuplicate = PlotNames[iPlot];
            TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("True","") ;
            textSlice->DrawLatexNDC(0.2, 0.81, Labels[iSample] + "      " + LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

            gPad->RedrawAxis();
            TString SaveDirectory = (PlotNames[iPlot].Contains("NoFSI")) ? "/PreFSI" : "/PostFSI";
            PlotCanvas->SaveAs("./Figs/InteBreakDown/"+SaveDirectory+"/InteBreakDown_"+Labels[iSample]+"_"+PlotNames[iPlot]+".png");
            delete PlotCanvas;
        } // End of the loop over the samples grabing the plots	
    } // End of the loop over the plots
} // End of the program

TString to_string_with_precision(double a_value, const int n = 3) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());
}
