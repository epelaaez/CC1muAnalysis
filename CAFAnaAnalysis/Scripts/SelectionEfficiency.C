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
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionEfficiency() {
    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // We will create efficiency plots using true variables and definining signal efficiency
    // as the number of reconstructed the events that pass our signal definition and are true
    // signal events over the total true signal events; these two histograms are plotted

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionEfficiency.root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    // Vectors to fill with variables and variable information to plot
    std::vector<TruthVar> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Muon angle
    Vars.push_back(kTruthMuonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

    // Leading proton angle
    Vars.push_back(kTruthLeadingProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthLeadingProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L}})");

    // Recoil proton angle
    Vars.push_back(kTruthRecoilProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthRecoilProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{R}})");

    // Opening angle between protons
    Vars.push_back(kTruthCosOpeningAngleProtons); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthCosOpeningAngleProtons"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");

    // Opening angle between muon and total proton
    Vars.push_back(kTruthCosOpeningAngleMuonTotalProton); VarBins.push_back(bAngleBins);
    PlotNames.push_back("TruthCosOpeningAngleMuonTotalProton"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");

    // Delta alpha transverse
    Vars.push_back(kTruthDeltaAlphaT); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("TruthDeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Transverse momentum
    Vars.push_back(kTruthTransverseMomentum); VarBins.push_back(bTransverseMomentumBins);
    PlotNames.push_back("TruthTransverseMomentum"); VarLabels.push_back("#delta P_{T}");

    // Muon momentum 
    Vars.push_back(kTruthMuonMomentum); VarBins.push_back(bMuonMomentumBins);
    PlotNames.push_back("TruthMuonMomentum"); VarLabels.push_back("|#vec{p}_{#mu}|");

    // Leading proton momentum 
    Vars.push_back(kTruthLeadingProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("TruthLeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back(kTruthRecoilProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("TruthRecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    ////////////////////////////////
    // Double differential variables
    ////////////////////////////////

    // Serial transverse momentum in muon cos theta
    Vars.push_back(kTruthTransverseMomentumInMuonCosTheta); VarBins.push_back(bTransverseMomentumInMuonCosTheta);
    PlotNames.push_back("TrueSerialTransverseMomentum_InMuonCosTheta"); VarLabels.push_back("#delta P_{T} (bin #)");

    // Delta alpha transverse in muon cos theta
    Vars.push_back(kTruthDeltaAlphaTInMuonCosTheta); VarBins.push_back(bDeltaAlphaTInMuonCosTheta);
    PlotNames.push_back("TrueSerialDeltaAlphaT_InMuonCosTheta"); VarLabels.push_back("#delta #alpha_{T} (bin #)");

    // Opening angle between protons in muon cos theta
    Vars.push_back(kTruthCosOpeningAngleProtonsInMuonCosTheta); VarBins.push_back(bCosOpeningAngleProtonsInMuonCosTheta);
    PlotNames.push_back("TrueSerialCosOpeningAngleProtons_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)");
    
    // Opening angle between muon and protons in muon cos theta
    Vars.push_back(kTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta); VarBins.push_back(bCosOpeningAngleMuonTotalProtonInMuonCosTheta);
    PlotNames.push_back("TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)");

    // Construct all spectra
    std::vector<std::tuple<std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto TrueSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kTruthIsSignal, kNoSpillCut);
        auto RecoTrueSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, Vars.at(i), kTruthIsSignal, kNoSpillCut, kRecoIsSignal);
        Spectra.push_back({std::move(TrueSignals), std::move(RecoTrueSignals)});
    }

    NuLoader.Go();

    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto& [TrueSignals, RecoTrueSignals] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1D* TrueHisto = TrueSignals->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);

        // Manage under/overflow bins
        TrueHisto->SetBinContent(TrueHisto->GetNbinsX(), TrueHisto->GetBinContent(TrueHisto->GetNbinsX()) + TrueHisto->GetBinContent(TrueHisto->GetNbinsX() + 1));
        RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));

        TrueHisto->SetBinContent(1, TrueHisto->GetBinContent(0) + TrueHisto->GetBinContent(1));
        RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));

        // Change y axis title so efficiency plot inherits it
        TrueHisto->GetYaxis()->SetTitle("Signal efficiency");
        RecoTrueHisto->GetYaxis()->SetTitle("Signal efficiency");

        // Canvas margins
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        if (PlotNames[i].Contains("Serial")) {
            // Flatten out double differential plots
            auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator[PlotNames[i]+"Plot"];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            int StartIndex = 0;

            // Loop over slices
            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                // Slice name
                TString SlicePlotName = PlotNames[i] + "_" + TString(std::to_string(iSlice));

                // Get slice width
                double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 

                // Get number of bins
                int SliceNBins = SerialVectorBins.at(iSlice);
                std::vector<double> SerialSliceBinning;

                for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                    double value = SerialVectorRanges.at(StartIndex + iBin);
                    SerialSliceBinning.push_back(value);
                } // End of the number of bins and the bin ranges declaration

                // Slice true and reco true histos
                TH1D* SlicedTrueHisto = tools.GetHistoBins(
                    TrueHisto,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "True"
                );
                TH1D* SlicedRecoTrueHisto = tools.GetHistoBins(
                    RecoTrueHisto,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "RecoTrue"
                );

                // Create efficiency plot
                TEfficiency* Eff = new TEfficiency(*SlicedRecoTrueHisto, *SlicedTrueHisto);
                Eff->SetTitle((";True " + VarLabels.at(i) + ";").c_str());

                PlotCanvas->cd();
                Eff->SetMarkerStyle(21);
                Eff->SetMarkerColor(kBlack);
                Eff->Draw("AP");
                gPad->Update();
                Eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(SerialSliceBinning.at(0),SerialSliceBinning.at(SerialSliceBinning.size() - 1));

                // Slice label
                TLatex *textSlice = new TLatex();
                TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel[PlotNames[i]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

                // Save as png
                PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Efficiency/"+SlicePlotName+".png");

                // Save to root file
                SaveFile->WriteObject(Eff, SlicePlotName+"_eff");
            }
        } else {
            // Create efficiency plot
            TEfficiency* Eff = new TEfficiency(*RecoTrueHisto, *TrueHisto);
            Eff->SetTitle((";True " + VarLabels.at(i) + ";").c_str());

            PlotCanvas->cd();
            Eff->SetMarkerStyle(21);
            Eff->SetMarkerColor(kBlack);
            Eff->SetTitle("All events");
            Eff->Draw("AP");
            gPad->Update();
            Eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(VarBins.at(i).Min(),VarBins.at(i).Max());

            // Save as png
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Efficiency/"+PlotNames[i]+".png");

            // Save to root file
            SaveFile->WriteObject(Eff, PlotNames[i]+"_eff");
        }
        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}
