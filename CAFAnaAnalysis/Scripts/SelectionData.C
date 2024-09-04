// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionData() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader("/exp/sbnd/data/users/munjung/SBND/data/run_14500_14503.flat.caf.root");

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Data.root";
    TFile* SaveFile = new TFile(RootFilePath, "RECREATE");

    std::vector<std::unique_ptr<Spectrum>> Spectra;

    // Variables after cutting 
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto RecoSignal = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal);
        Spectra.push_back(std::move(RecoSignal));
    }

    // Variables without any cut
    std::vector<Var> NoCutVars; std::vector<Binning> NoCutBins; std::vector<std::string> NoCutLabels; std::vector<TString> NoCutNames;
    NoCutVars.push_back(kVertexX); NoCutBins.push_back(bVertexX); NoCutLabels.push_back("#vec{v}_{x}"); NoCutNames.push_back("VertexX");
    NoCutVars.push_back(kVertexY); NoCutBins.push_back(bVertexY); NoCutLabels.push_back("#vec{v}_{y}"); NoCutNames.push_back("VertexY");
    NoCutVars.push_back(kVertexZ); NoCutBins.push_back(bVertexZ); NoCutLabels.push_back("#vec{v}_{z}"); NoCutNames.push_back("VertexZ");
    for (std::size_t i = 0; i < NoCutVars.size(); i++) {
        auto RecoSignal = std::make_unique<Spectrum>(NoCutLabels.at(i), NoCutBins.at(i), NuLoader, NoCutVars.at(i), kNoSpillCut, kNoInvalidVariables);
        Spectra.push_back(std::move(RecoSignal));
    }
    
    NuLoader.Go();

    // Variables after cutting
    for (std::size_t iVar = 0; iVar < Vars.size(); iVar++) {
        auto& RecoSignal = Spectra[iVar];

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TH1D* RecoHist = RecoSignal->ToTH1(TargetPOT);

        RecoHist->SetBinContent(RecoHist->GetNbinsX(), RecoHist->GetBinContent(RecoHist->GetNbinsX()) + RecoHist->GetBinContent(RecoHist->GetNbinsX() + 1));
        RecoHist->SetBinContent(1, RecoHist->GetBinContent(0) + RecoHist->GetBinContent(1));

        double Area = RecoHist->Integral();
        RecoHist->Scale(1 / Area);
        RecoHist->SetLineColor(kBlue+2);
        RecoHist->SetLineWidth(4);

        RecoHist->GetXaxis()->SetTitleFont(FontStyle);
        RecoHist->GetXaxis()->SetLabelFont(FontStyle);
        RecoHist->GetXaxis()->SetNdivisions(8);
        RecoHist->GetXaxis()->SetLabelSize(TextSize);
        RecoHist->GetXaxis()->SetTitleSize(TextSize);
        RecoHist->GetXaxis()->SetTitleOffset(1.1);
        RecoHist->GetXaxis()->CenterTitle();
        RecoHist->GetXaxis()->SetTitle(("Reco " + VarLabels.at(iVar)).c_str());

        RecoHist->GetYaxis()->SetTitleFont(FontStyle);
        RecoHist->GetYaxis()->SetLabelFont(FontStyle);
        RecoHist->GetYaxis()->SetNdivisions(6);
        RecoHist->GetYaxis()->SetLabelSize(TextSize);
        RecoHist->GetYaxis()->SetTitleSize(TextSize);
        RecoHist->GetYaxis()->SetTitleOffset(1.3);
        RecoHist->GetYaxis()->SetTickSize(0);
        RecoHist->GetYaxis()->CenterTitle();

        RecoHist->GetYaxis()->SetRangeUser(0., RecoHist->GetMaximum() * 1.3);
        PlotCanvas->cd();
        RecoHist->Draw("hist");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/DataCounts/"+PlotNames[iVar]+".png");

        delete PlotCanvas;
    }

    // Variables without any cutting
    for (std::size_t iVar = Vars.size(); iVar < Vars.size() + NoCutVars.size(); iVar++) {
        auto& RecoSignal = Spectra[iVar];

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TH1D* RecoHist = RecoSignal->ToTH1(TargetPOT);

        RecoHist->SetBinContent(RecoHist->GetNbinsX(), RecoHist->GetBinContent(RecoHist->GetNbinsX()) + RecoHist->GetBinContent(RecoHist->GetNbinsX() + 1));
        RecoHist->SetBinContent(1, RecoHist->GetBinContent(0) + RecoHist->GetBinContent(1));

        double Area = RecoHist->Integral();
        RecoHist->Scale(1 / Area);
        RecoHist->SetLineColor(kBlue+2);
        RecoHist->SetLineWidth(4);

        RecoHist->GetXaxis()->SetTitleFont(FontStyle);
        RecoHist->GetXaxis()->SetLabelFont(FontStyle);
        RecoHist->GetXaxis()->SetNdivisions(8);
        RecoHist->GetXaxis()->SetLabelSize(TextSize);
        RecoHist->GetXaxis()->SetTitleSize(TextSize);
        RecoHist->GetXaxis()->SetTitleOffset(1.1);
        RecoHist->GetXaxis()->CenterTitle();
        RecoHist->GetXaxis()->SetTitle(("Reco " + NoCutLabels.at(iVar - Vars.size())).c_str());

        RecoHist->GetYaxis()->SetTitleFont(FontStyle);
        RecoHist->GetYaxis()->SetLabelFont(FontStyle);
        RecoHist->GetYaxis()->SetNdivisions(6);
        RecoHist->GetYaxis()->SetLabelSize(TextSize);
        RecoHist->GetYaxis()->SetTitleSize(TextSize);
        RecoHist->GetYaxis()->SetTitleOffset(1.3);
        RecoHist->GetYaxis()->SetTickSize(0);
        RecoHist->GetYaxis()->CenterTitle();

        RecoHist->GetYaxis()->SetRangeUser(0., RecoHist->GetMaximum() * 1.3);
        PlotCanvas->cd();
        RecoHist->Draw("hist");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/DataCounts/NoCut/"+NoCutNames[iVar - Vars.size()]+".png");

        delete PlotCanvas;
    }

    SaveFile->Close();
}