
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

void SelectionStubs() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    std::vector<Var> StubVars = {kPStubMult, kPStubProtonMult, kPStubMuonMult};
    std::vector<TString> StubVarNames = {"NStubs", "ProtonMult", "MuonMult"};
    std::vector<TString> StubVarLabel = {"# stubs", "# proton stubs", "# muon stubs"};
    std::vector<Cut> Cuts = {kRecoIsSignal, kRecoIsTrueReco, kRecoIsBackground};
    std::vector<TString> CutNames = {"Reco", "RecoTrue", "RecoBkg"};

    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;

    for (std::size_t i = 0; i < Cuts.size(); ++i) {
        auto NStubs = std::make_unique<Spectrum>("NStubs" + (std::string)CutNames[i], bStubMult, NuLoader, kPStubMult, kNoSpillCut, Cuts[i]); 
        auto ProtonStubs = std::make_unique<Spectrum>("ProtonStubs" + (std::string)CutNames[i], bStubMult, NuLoader, kPStubProtonMult, kNoSpillCut, Cuts[i]); 
        auto MuonStubs = std::make_unique<Spectrum>("MuonStubs" + (std::string)CutNames[i], bStubMult, NuLoader, kPStubMuonMult, kNoSpillCut, Cuts[i]); 
        Spectra.push_back({std::move(NStubs), std::move(ProtonStubs), std::move(MuonStubs)});
    }
    for (std::size_t i = 0; i < StubVars.size(); ++i) {
        auto RecoStubs = std::make_unique<Spectrum>("Reco" + (std::string)StubVarNames[i], bStubMult, NuLoader, StubVars[i], kNoSpillCut, kRecoIsSignal);
        auto RecoSignalStubs = std::make_unique<Spectrum>("RecoSignal" + (std::string)StubVarNames[i], bStubMult, NuLoader, StubVars[i], kNoSpillCut, kRecoIsTrueReco);
        auto RecoBkgStubs = std::make_unique<Spectrum>("RecoBkg" + (std::string)StubVarNames[i], bStubMult, NuLoader, StubVars[i], kNoSpillCut, kRecoIsBackground);
        Spectra.push_back({std::move(RecoStubs), std::move(RecoSignalStubs), std::move(RecoBkgStubs)});
    }
    NuLoader.Go();

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);  
    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.17);
    PlotCanvas->SetRightMargin(0.05);
    PlotCanvas->SetBottomMargin(0.16);

    for (std::size_t i = 0; i < Cuts.size(); ++i) {
        auto& [NStubs, ProtonStubs, MuonStubs] = Spectra[i];

        TH1D* NStubsHisto = NStubs->ToTH1(TargetPOT);
        TH1D* ProtonStubsHisto = ProtonStubs->ToTH1(TargetPOT);
        TH1D* MuonStubsHisto = MuonStubs->ToTH1(TargetPOT);  

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(NStubsHisto,"Total stubs","l");
        NStubsHisto->SetLineColor(kBlue+2);
        NStubsHisto->SetLineWidth(4);

        // Style histograms
        NStubsHisto->GetXaxis()->SetTitleFont(FontStyle);
        NStubsHisto->GetXaxis()->SetLabelFont(FontStyle);
        NStubsHisto->GetXaxis()->SetNdivisions(8);
        NStubsHisto->GetXaxis()->SetLabelSize(TextSize);
        NStubsHisto->GetXaxis()->SetTitleSize(TextSize);
        NStubsHisto->GetXaxis()->SetTitleOffset(1.1);
        NStubsHisto->GetXaxis()->CenterTitle();
        NStubsHisto->GetXaxis()->SetTitle("# stubs");

        NStubsHisto->GetYaxis()->SetTitleFont(FontStyle);
        NStubsHisto->GetYaxis()->SetLabelFont(FontStyle);
        NStubsHisto->GetYaxis()->SetNdivisions(6);
        NStubsHisto->GetYaxis()->SetLabelSize(TextSize);
        NStubsHisto->GetYaxis()->SetTitleSize(TextSize);
        NStubsHisto->GetYaxis()->SetTitleOffset(1.3);
        NStubsHisto->GetYaxis()->SetTickSize(0);
        NStubsHisto->GetYaxis()->CenterTitle();
        NStubsHisto->GetYaxis()->SetTitle("# events");

        double imax = MuonStubsHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        NStubsHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

        TLegendEntry* legRecoTrue = leg->AddEntry(ProtonStubsHisto,"Proton stubs","l");
        ProtonStubsHisto->SetLineColor(kRed+1);
        ProtonStubsHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(MuonStubsHisto,"Muon stubs","l");
        MuonStubsHisto->SetLineColor(kOrange+7);
        MuonStubsHisto->SetLineWidth(4);

        PlotCanvas->cd();
        NStubsHisto->Draw("hist");
        ProtonStubsHisto->Draw("hist same");
        MuonStubsHisto->Draw("hist same");
        leg->Draw();

        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Stubs/Multiplicity" + CutNames[i] + ".png");
    }

    for (std::size_t i = 0; i < StubVars.size(); ++i) {
        auto& [RecoStubs, RecoSignalStubs, RecoBkgStubs] = Spectra[i + Cuts.size()];

        TH1D* RecoStubsHisto = RecoStubs->ToTH1(TargetPOT);
        TH1D* RecoSignalStubsHisto = RecoSignalStubs->ToTH1(TargetPOT);
        TH1D* RecoBkgStubsHisto = RecoBkgStubs->ToTH1(TargetPOT);  

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoStubsHisto,"Reconstructed","l");
        RecoStubsHisto->SetLineColor(kBlue+2);
        RecoStubsHisto->SetLineWidth(4);

        // Style histograms
        RecoStubsHisto->GetXaxis()->SetTitleFont(FontStyle);
        RecoStubsHisto->GetXaxis()->SetLabelFont(FontStyle);
        RecoStubsHisto->GetXaxis()->SetNdivisions(8);
        RecoStubsHisto->GetXaxis()->SetLabelSize(TextSize);
        RecoStubsHisto->GetXaxis()->SetTitleSize(TextSize);
        RecoStubsHisto->GetXaxis()->SetTitleOffset(1.1);
        RecoStubsHisto->GetXaxis()->CenterTitle();
        RecoStubsHisto->GetXaxis()->SetTitle(StubVarLabel[i]);

        RecoStubsHisto->GetYaxis()->SetTitleFont(FontStyle);
        RecoStubsHisto->GetYaxis()->SetLabelFont(FontStyle);
        RecoStubsHisto->GetYaxis()->SetNdivisions(6);
        RecoStubsHisto->GetYaxis()->SetLabelSize(TextSize);
        RecoStubsHisto->GetYaxis()->SetTitleSize(TextSize);
        RecoStubsHisto->GetYaxis()->SetTitleOffset(1.3);
        RecoStubsHisto->GetYaxis()->SetTickSize(0);
        RecoStubsHisto->GetYaxis()->CenterTitle();
        RecoStubsHisto->GetYaxis()->SetTitle("# events");

        double imax = RecoStubsHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoStubsHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoSignalStubsHisto,"Signal","l");
        RecoSignalStubsHisto->SetLineColor(kRed+1);
        RecoSignalStubsHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgStubsHisto,"Background","l");
        RecoBkgStubsHisto->SetLineColor(kOrange+7);
        RecoBkgStubsHisto->SetLineWidth(4);

        PlotCanvas->cd();
        RecoStubsHisto->Draw("hist");
        RecoSignalStubsHisto->Draw("hist same");
        RecoBkgStubsHisto->Draw("hist same");
        leg->Draw();

        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Stubs/Reco" + StubVarNames[i] + ".png");
    }
}