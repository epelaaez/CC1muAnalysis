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
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

// std includes.
#include <filesystem>
#include <vector>
#include <string>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Tools.cxx"
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionDetectorSystematics() {
    Tools tools;

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // Target files for each sample
    const std::string CVPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04_3/cv";
    const std::vector<std::string> CVInputFiles = tools.GetInputFiles(CVPath);
    SpectrumLoader CVLoader(CVInputFiles);

    const std::string NoDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04_2/no_diffusion";
    const std::vector<std::string> NoDifInputFiles = tools.GetInputFiles(NoDifPath);
    SpectrumLoader NoDifLoader(NoDifInputFiles);

    const std::string NoLonDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04_2/no_longitudinal_diffusion";
    const std::vector<std::string> NoLonDifInputFiles = tools.GetInputFiles(NoLonDifPath);
    SpectrumLoader NoLonDifLoader(NoLonDifInputFiles);

    const std::string NoTraDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04_2/no_transverse_diffusion";
    const std::vector<std::string> NoTraDifInputFiles = tools.GetInputFiles(NoTraDifPath);
    SpectrumLoader NoTraDifLoader(NoTraDifInputFiles);

    const std::string SCENoDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04_2/with_sce_no_diffusion";
    const std::vector<std::string> SCENoDifInputFiles = tools.GetInputFiles(SCENoDifPath);
    SpectrumLoader SCENoDifLoader(SCENoDifInputFiles);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/Uncertainties/Detector");

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsDetector.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    const int NVars = Vars.size();

    // Construct all spectra
    auto CVSignal = std::make_unique<Spectrum>("event count", bEventCount, CVLoader, kEventCount, kSpillPrintFile, kNoCut); 
    auto NoDifSignals = std::make_unique<Spectrum>("event count", bEventCount, NoDifLoader, kEventCount, kNoSpillCut, kNoCut); 
    auto NoLonDifSignals = std::make_unique<Spectrum>("event count", bEventCount, NoLonDifLoader, kEventCount, kNoSpillCut, kNoCut); 
    auto NoTraDifSignals = std::make_unique<Spectrum>("event count", bEventCount, NoTraDifLoader, kEventCount, kNoSpillCut, kNoCut); 
    auto SCENoDifSignals = std::make_unique<Spectrum>("event count", bEventCount, SCENoDifLoader, kEventCount, kNoSpillCut, kNoCut); 

    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;
    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto CVRecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), CVLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        auto NoDifRecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NoDifLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        auto NoLonDifRecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NoLonDifLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        auto NoTraDifRecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NoTraDifLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        auto SCENoDifRecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), SCENoDifLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        Spectra.push_back({std::move(CVRecoSignals), std::move(NoDifRecoSignals), std::move(NoLonDifRecoSignals), std::move(NoTraDifRecoSignals), std::move(SCENoDifRecoSignals)});
    }
    CVLoader.Go(); NoDifLoader.Go(); NoLonDifLoader.Go(); NoTraDifLoader.Go(); SCENoDifLoader.Go();

    TH1D* CVCountHisto = CVSignal->ToTH1(1.23677e17);
    TH1D* NoDifCountHisto = NoDifSignals->ToTH1(8.78885e16);
    TH1D* NoLonDifCountHisto = NoLonDifSignals->ToTH1(8.50138e16);
    TH1D* NoTraDifCountHisto = NoTraDifSignals->ToTH1(8.65587e16);
    TH1D* SCENoDifCountHisto = SCENoDifSignals->ToTH1(8.50441e16);

    std::cout << CVCountHisto->Integral() << std::endl;
    std::cout << NoDifCountHisto->Integral() << std::endl;
    std::cout << NoLonDifCountHisto->Integral() << std::endl;
    std::cout << NoTraDifCountHisto->Integral() << std::endl;
    std::cout << SCENoDifCountHisto->Integral() << std::endl;

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.17);
    PlotCanvas->SetRightMargin(0.05);
    PlotCanvas->SetBottomMargin(0.16);

    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto& [CVSpectra, NoDifSpectra, NoLonDifSpectra, NoTraDifSpectra, SCENoDifSpectra] = Spectra.at(iVar);

        // Get plots
        TH1D* CVHisto = CVSpectra->ToTH1(1.23677e17);
        TH1D* NoDifHisto = NoDifSpectra->ToTH1(8.78885e16);
        TH1D* NoLonDifHisto = NoLonDifSpectra->ToTH1(8.50138e16);
        TH1D* NoTraDifHisto = NoTraDifSpectra->ToTH1(8.65587e16);
        TH1D* SCENoDifHisto = SCENoDifSpectra->ToTH1(8.50441e16);

        std::cout << CVHisto->Integral() << std::endl;
        std::cout << NoDifHisto->Integral() << std::endl;
        std::cout << NoLonDifHisto->Integral() << std::endl;
        std::cout << NoTraDifHisto->Integral() << std::endl;
        std::cout << SCENoDifHisto->Integral() << std::endl;

        TLegend* leg = new TLegend(0.2,0.73,0.8,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legCV = leg->AddEntry(CVHisto,"CV","l");
        CVHisto->SetLineColor(kBlue+2);
        CVHisto->SetLineWidth(4);

        CVHisto->GetXaxis()->SetTitleFont(FontStyle);
        CVHisto->GetXaxis()->SetLabelFont(FontStyle);
        CVHisto->GetXaxis()->SetNdivisions(8);
        CVHisto->GetXaxis()->SetLabelSize(TextSize);
        CVHisto->GetXaxis()->SetTitleSize(TextSize);
        CVHisto->GetXaxis()->SetTitleOffset(1.1);
        CVHisto->GetXaxis()->CenterTitle();
        CVHisto->GetXaxis()->SetTitle(("Reco " + VarLabels.at(iVar)).c_str());

        CVHisto->GetYaxis()->SetTitleFont(FontStyle);
        CVHisto->GetYaxis()->SetLabelFont(FontStyle);
        CVHisto->GetYaxis()->SetNdivisions(6);
        CVHisto->GetYaxis()->SetLabelSize(TextSize);
        CVHisto->GetYaxis()->SetTitleSize(TextSize);
        CVHisto->GetYaxis()->SetTitleOffset(1.3);
        CVHisto->GetYaxis()->SetTickSize(0);
        CVHisto->GetYaxis()->CenterTitle();

        TLegendEntry* legNoDif = leg->AddEntry(NoDifHisto,"No diff.","l");
        NoDifHisto->SetLineColor(kRed+1);
        NoDifHisto->SetLineWidth(4);

        TLegendEntry* legNoLonDif = leg->AddEntry(NoLonDifHisto,"No lon. diff.","l");
        NoLonDifHisto->SetLineColor(kOrange+7);
        NoLonDifHisto->SetLineWidth(4);

        TLegendEntry* legNoTraDif = leg->AddEntry(NoTraDifHisto,"No tran. diff.","l");
        NoTraDifHisto->SetLineColor(kMagenta);
        NoTraDifHisto->SetLineWidth(4);

        TLegendEntry* legSCENoDif = leg->AddEntry(SCENoDifHisto,"SCE no diff.","l");
        SCENoDifHisto->SetLineColor(kGreen+1);
        SCENoDifHisto->SetLineWidth(4);

        double Max = CVHisto->GetMaximum();
        Max = std::max(Max, NoDifHisto->GetMaximum());
        Max = std::max(Max, NoLonDifHisto->GetMaximum());
        Max = std::max(Max, NoTraDifHisto->GetMaximum());
        Max = std::max(Max, SCENoDifHisto->GetMaximum());
        CVHisto->GetYaxis()->SetRangeUser(0., 1.3*Max);

        PlotCanvas->cd();
        CVHisto->Draw("hist");
        NoDifHisto->Draw("hist same");
        NoLonDifHisto->Draw("hist same");
        NoTraDifHisto->Draw("hist same");
        SCENoDifHisto->Draw("hist same");
        leg->Draw();

        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Detector/"+PlotNames[iVar]+".png");

        // // Scale 15%
        // ModifiedCVHisto->Scale(1.15);
        // ModifiedRecoTrueHisto->Scale(1.15);
        // ModifiedRecoBkgHisto->Scale(1.15);

        // SelectionHelpers::DrawHistosWithErrorBands(
        //     CVHisto,
        //     RecoTrueHisto,
        //     RecoBkgHisto,
        //     {ModifiedCVHisto},
        //     {ModifiedRecoTrueHisto},
        //     {ModifiedRecoBkgHisto},
        //     dir,
        //     "Detector",
        //     PlotNames[iVar]
        // );

        // // Scale histograms for cov matrices
        // CVHisto->Scale(Units / (IntegratedFlux * NTargets));
        // RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
        // RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

        // // Scale modified histograms for cov matrices
        // ModifiedCVHisto->Scale(Units / (IntegratedFlux * NTargets));
        // ModifiedRecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
        // ModifiedRecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

        // // Plot cov, frac cov, and corr matrices
        // SelectionHelpers::DrawMatrices(
        //     CVHisto,
        //     RecoTrueHisto,
        //     RecoBkgHisto,
        //     {ModifiedCVHisto},
        //     {ModifiedRecoTrueHisto},
        //     {ModifiedRecoBkgHisto},
        //     dir,
        //     "Detector",
        //     VarLabels[iVar],
        //     PlotNames[iVar],
        //     SaveFile
        // );
    }
    // Close file
    SaveFile->Close();
}
