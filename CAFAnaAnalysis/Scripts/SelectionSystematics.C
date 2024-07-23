// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/ISyst.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionSystematics(std::string SystName, int SystNUniv) {
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Systematic with name " << SystName << ", and number of universes " << SystNUniv <<  std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // Get integrated flux
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root");
    TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
    double IntegratedFlux = (HistoFlux->Integral() * TargetPOT / POTPerSpill / Nominal_UB_XY_Surface);    

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // We now create overlaid plots for several reconstructed variables and three lines:
    //     1. all selected reconstructed events
    //     2. reco signal events
    //     3. reco background events

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/Uncertainties/"+SystName);

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/epelaez/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    // Vectors to fill with variables and variable information to plot
    std::vector<Var> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;
    
    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Muon angle
    Vars.push_back(kMuonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("MuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu}})");

    // Leading proton angle
    Vars.push_back(kLeadingProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("LeadingProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L}})");

    // Recoil proton angle
    Vars.push_back(kRecoilProtonCosTheta); VarBins.push_back(bAngleBins);
    PlotNames.push_back("RecoilProtonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{R}})");

    // Opening angle between protons
    Vars.push_back(kCosOpeningAngleProtons); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleProtons"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");

    // Opening angle between muon and total proton
    Vars.push_back(kCosOpeningAngleMuonTotalProton); VarBins.push_back(bAngleBins);
    PlotNames.push_back("CosOpeningAngleMuonTotalProton"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");

    // Delta alpha transverse
    Vars.push_back(kDeltaAlphaT); VarBins.push_back(bDeltaAlphaBins);
    PlotNames.push_back("DeltaAlphaT"); VarLabels.push_back("#delta #alpha_{T}");

    // Transverse momentum
    Vars.push_back(kTransverseMomentum); VarBins.push_back(bTransverseMomentumBins);
    PlotNames.push_back("TransverseMomentum"); VarLabels.push_back("#delta P_{T}");

    // Muon momentum 
    Vars.push_back(kMuonMomentum); VarBins.push_back(bMuonMomentumBins);
    PlotNames.push_back("MuonMomentum"); VarLabels.push_back("|#vec{p}_{#mu}|");

    // Leading proton momentum 
    Vars.push_back(kLeadingProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("LeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back(kRecoilProtonMomentum); VarBins.push_back(bProtonMomentumBins);
    PlotNames.push_back("RecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    ////////////////////////////////
    // Double differential variables
    ////////////////////////////////

    // Serial transverse momentum in muon cos theta
    Vars.push_back(kTransverseMomentumInMuonCosTheta); VarBins.push_back(bTransverseMomentumInMuonCosTheta);
    PlotNames.push_back("SerialTransverseMomentum_InMuonCosTheta"); VarLabels.push_back("#delta P_{T} (bin #)");

    // Delta alpha transverse in muon cos theta
    Vars.push_back(kDeltaAlphaTInMuonCosTheta); VarBins.push_back(bDeltaAlphaTInMuonCosTheta);
    PlotNames.push_back("SerialDeltaAlphaT_InMuonCosTheta"); VarLabels.push_back("#delta #alpha_{T} (bin #)");

    // Opening angle between protons in muon cos theta
    Vars.push_back(kCosOpeningAngleProtonsInMuonCosTheta); VarBins.push_back(bCosOpeningAngleProtonsInMuonCosTheta);
    PlotNames.push_back("SerialCosOpeningAngleProtons_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)");
    
    // Opening angle between muon and protons in muon cos theta
    Vars.push_back(kCosOpeningAngleMuonTotalProtonInMuonCosTheta); VarBins.push_back(bCosOpeningAngleMuonTotalProtonInMuonCosTheta);
    PlotNames.push_back("SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta"); VarLabels.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)");
    
    // Create shift depending on number of universes
    ISyst* syst = new SBNWeightSyst(SystName);
    std::vector<SystShifts> Shifts;

    if (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) {
    	// Add +1 sigma shift
        SystShifts SigP1Shift(syst, +1);
	    Shifts.push_back(SigP1Shift);
    } else {
	    // Add random Gaussian shifts
        for (int i = 0; i < SystNUniv; i++) {
            SystShifts RandomShift(syst, gRandom->Gaus(0,1));
            Shifts.push_back(RandomShift);
        }
    }

    // We now have the option to either load all the spectra from a previous run or 
    // run the spectra in this run
    const bool ConstructSpectra = true;

    // Where we store spectra if we are going to construct them    
    std::vector<std::tuple<
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>
    >> Spectra;

    // Where we load histograms if we do not construct spectra
    std::vector<std::vector<std::tuple<TH1*, TH1*, TH1*>>> LoadedHistos;

    if (ConstructSpectra) {
        // Construct all spectra
        for (std::size_t iVar = 0; iVar < Vars.size(); iVar++) {
            // Create reco spectrum with shift
            auto RecoSpectra = std::make_unique<EnsembleSpectrum>(
                NuLoader,
                HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)),
                kNoSpillCut,
                kRecoIsSignal,
                Shifts
            );
            // Create reco true signal spectrum with shift
            auto RecoTrueSpectra = std::make_unique<EnsembleSpectrum>(
                NuLoader,
                HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)),
                kNoSpillCut,
                kRecoIsTrueReco,
                Shifts
            );
            // Create reco background spectrum with shift
            auto RecoBkgSpectra = std::make_unique<EnsembleSpectrum>(
                NuLoader,
                HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)),
                kNoSpillCut,
                kRecoIsBackground,
                Shifts
            );
            // Add everything to main vector
            Spectra.push_back({std::move(RecoSpectra), std::move(RecoTrueSpectra), std::move(RecoBkgSpectra)});
        }
        // Load spectra
        NuLoader.Go();
    } else {
        // Load previously constructed histograms from file
        for (std::size_t i = 0; i < Vars.size(); i++) {
            std::vector<std::tuple<TH1*, TH1*, TH1*>> VarHistos;

            // Nominal plots
            TH1* RecoHisto = (TH1*)(SaveFile->Get<TH1>(PlotNames[i]+"_reco"));
            TH1* RecoTrueHisto = (TH1*)(SaveFile->Get<TH1>(PlotNames[i]+"_reco_true"));
            TH1* RecoBkgHisto = (TH1*)(SaveFile->Get<TH1>(PlotNames[i]+"_reco_bkg"));
            VarHistos.push_back({std::move(RecoHisto), std::move(RecoTrueHisto), std::move(RecoBkgHisto)});

            // Var univ plots
            int NUniv = (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) ? 1 : SystNUniv;
            for (int iUniv = 0; iUniv < NUniv; iUniv++) {
                TString UnivString = TString(std::to_string(iUniv));
                TH1* UnivRecoHisto = (TH1D*)(SaveFile->Get<TH1>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco"));
                TH1* UnivRecoTrueHisto = (TH1D*)(SaveFile->Get<TH1>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_true"));
                TH1* UnivRecoBkgHisto = (TH1D*)(SaveFile->Get<TH1>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_bkg"));
                VarHistos.push_back({std::move(UnivRecoHisto), std::move(UnivRecoTrueHisto), std::move(UnivRecoBkgHisto)});
            }
            LoadedHistos.push_back(std::move(VarHistos));

            // Resize Spectra so everything compiles correctly
            Spectra.resize(Vars.size());
        }
    }

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        // Get histograms
        TH1* RecoHisto; TH1* RecoTrueHisto; TH1* RecoBkgHisto;
        auto& [RecoSpectra, RecoTrueSpectra, RecoBkgSpectra] = Spectra.at(i);
        if (ConstructSpectra) {
            RecoHisto = RecoSpectra->Nominal().ToTH1(TargetPOT);
            RecoTrueHisto = RecoTrueSpectra->Nominal().ToTH1(TargetPOT);
            RecoBkgHisto = RecoBkgSpectra->Nominal().ToTH1(TargetPOT);

            // Manage under/overflow bins
            RecoHisto->SetBinContent(RecoHisto->GetNbinsX(), RecoHisto->GetBinContent(RecoHisto->GetNbinsX()) + RecoHisto->GetBinContent(RecoHisto->GetNbinsX() + 1));
            RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));
            RecoBkgHisto->SetBinContent(RecoBkgHisto->GetNbinsX(), RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX()) + RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX() + 1));

            RecoHisto->SetBinContent(1, RecoHisto->GetBinContent(0) + RecoHisto->GetBinContent(1));
            RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));
            RecoBkgHisto->SetBinContent(1, RecoBkgHisto->GetBinContent(0) + RecoBkgHisto->GetBinContent(1));
        } else {
            RecoHisto = std::get<0>(LoadedHistos.at(i)[0]);
            RecoTrueHisto = std::get<1>(LoadedHistos.at(i)[0]);
            RecoBkgHisto = std::get<2>(LoadedHistos.at(i)[0]);
        }

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoHisto,"Reconstructed","l");
        RecoHisto->SetLineColor(kBlue+2);
        RecoHisto->SetLineWidth(4);

        // Style histograms
        RecoHisto->GetXaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetXaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetXaxis()->SetNdivisions(8);
        RecoHisto->GetXaxis()->SetLabelSize(TextSize);
        RecoHisto->GetXaxis()->SetTitleSize(TextSize);
        RecoHisto->GetXaxis()->SetTitleOffset(1.1);
        RecoHisto->GetXaxis()->CenterTitle();
        RecoHisto->GetXaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());

        RecoHisto->GetYaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetYaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetYaxis()->SetNdivisions(6);
        RecoHisto->GetYaxis()->SetLabelSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleOffset(1.3);
        RecoHisto->GetYaxis()->SetTickSize(0);
        RecoHisto->GetYaxis()->CenterTitle();

        double imax = RecoHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"True","l");
        RecoTrueHisto->SetLineColor(kRed+1);
        RecoTrueHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"Background","l");
        RecoBkgHisto->SetLineColor(kOrange+7);
        RecoBkgHisto->SetLineWidth(4);

        // Create covariance matrix
        std::string CovName = "Cov" + SystName;
        TH2* CovMatrix = new TH2D(
            (CovName + (std::string)PlotNames[i]).c_str(),
            CovName.c_str(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max()
        );

        // Create fractional covariance matrix
        std::string FracCovName = "FracCov" + SystName;
        TH2* FracCovMatrix = new TH2D(
            (FracCovName + (std::string)PlotNames[i]).c_str(),
            FracCovName.c_str(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max()
        );

        // Create correlation matrix
        std::string CorrName = "Corr" + SystName;
        TH2* CorrMatrix = new TH2D(
            (CorrName + (std::string)PlotNames[i]).c_str(),
            CorrName.c_str(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max(),
            VarBins.at(i).NBins(),
            VarBins.at(i).Min(),
            VarBins.at(i).Max()
        );

        // Loop over all universes
        int NUniv = (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) ? 1 : SystNUniv;
        for (int iUniv = 0; iUniv < NUniv; iUniv++) {
            // Get histograms
            TH1* UnivRecoHisto; TH1* UnivRecoTrueHisto; TH1* UnivRecoBkgHisto;
            if (ConstructSpectra) {
                UnivRecoHisto = RecoSpectra->Universe(iUniv).ToTH1(TargetPOT);
                UnivRecoTrueHisto = RecoTrueSpectra->Universe(iUniv).ToTH1(TargetPOT);
                UnivRecoBkgHisto = RecoBkgSpectra->Universe(iUniv).ToTH1(TargetPOT);

                // Manage under/overflow bins
                UnivRecoHisto->SetBinContent(UnivRecoHisto->GetNbinsX(), UnivRecoHisto->GetBinContent(UnivRecoHisto->GetNbinsX()) + UnivRecoHisto->GetBinContent(UnivRecoHisto->GetNbinsX() + 1));
                UnivRecoTrueHisto->SetBinContent(UnivRecoTrueHisto->GetNbinsX(), UnivRecoTrueHisto->GetBinContent(UnivRecoTrueHisto->GetNbinsX()) + UnivRecoTrueHisto->GetBinContent(UnivRecoTrueHisto->GetNbinsX() + 1));
                UnivRecoBkgHisto->SetBinContent(UnivRecoBkgHisto->GetNbinsX(), UnivRecoBkgHisto->GetBinContent(UnivRecoBkgHisto->GetNbinsX()) + UnivRecoBkgHisto->GetBinContent(UnivRecoBkgHisto->GetNbinsX() + 1));

                UnivRecoHisto->SetBinContent(1, UnivRecoHisto->GetBinContent(0) + UnivRecoHisto->GetBinContent(1));
                UnivRecoTrueHisto->SetBinContent(1, UnivRecoTrueHisto->GetBinContent(0) + UnivRecoTrueHisto->GetBinContent(1));
                UnivRecoBkgHisto->SetBinContent(1, UnivRecoBkgHisto->GetBinContent(0) + UnivRecoBkgHisto->GetBinContent(1));
            } else {
                UnivRecoHisto = std::get<0>(LoadedHistos.at(i)[iUniv + 1]);
                UnivRecoTrueHisto = std::get<1>(LoadedHistos.at(i)[iUniv + 1]);
                UnivRecoBkgHisto = std::get<2>(LoadedHistos.at(i)[iUniv + 1]);
            }

            for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
                double XEventRateCV = (RecoHisto->GetBinContent(x) / (IntegratedFlux * NTargets)) * Units;
                double XEventRateVar = (UnivRecoHisto->GetBinContent(x) / (IntegratedFlux * NTargets)) * Units;
                for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                    double YEventRateCV = (RecoHisto->GetBinContent(y) / (IntegratedFlux * NTargets)) * Units;
                    double YEventRateVar = (UnivRecoHisto->GetBinContent(y) / (IntegratedFlux * NTargets)) * Units; 
                    double Value = ((XEventRateVar - XEventRateCV) * (YEventRateVar - YEventRateCV)) / NUniv;

                    // Fill covariance matrix
                    CovMatrix->Fill(
                        RecoHisto->GetXaxis()->GetBinCenter(x),
                        RecoHisto->GetXaxis()->GetBinCenter(y),
                        TMath::Max(Value, 1e-8)
                    );
                }
	        }

            // Save syst univ spectrum
            TString UnivString = TString(std::to_string(iUniv));
            if (ConstructSpectra) {
                SaveFile->WriteObject(UnivRecoHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco");
                SaveFile->WriteObject(UnivRecoTrueHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_true");
                SaveFile->WriteObject(UnivRecoBkgHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_bkg");
            }
        }

        // Create fractional covariance and correlation matrices
        for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
            double XEventRateCV = (RecoHisto->GetBinContent(x) / (IntegratedFlux * NTargets)) * Units;
            for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                double YEventRateCV = (RecoHisto->GetBinContent(y) / (IntegratedFlux * NTargets)) * Units;
                double CovBinValue = CovMatrix->GetBinContent(x,y);
                double XBinValue = CovMatrix->GetBinContent(x,x);
                double YBinValue = CovMatrix->GetBinContent(y,y);

                // Fill frac cov matrix
                double FracValue = (XBinValue == 0. || YEventRateCV == 0.) ? 0. : CovBinValue / (XEventRateCV * YEventRateCV);
                FracCovMatrix->SetBinContent(x, y, TMath::Max(FracValue, 1e-8));

                // Fill corr matrix
                double CorrValue = (XBinValue == 0. || YBinValue == 0.) ? 0. : CovBinValue / (TMath::Sqrt(XBinValue) * TMath::Sqrt(YBinValue));
                CorrMatrix->SetBinContent(x, y, TMath::Max(CorrValue, 1e-8));
            }
        }
            
        // Plot cov matrix
        double CovMin = CovMatrix->GetMinimum();
        double CovMax = CovMatrix->GetMaximum();
        CovMatrix->GetZaxis()->SetRangeUser(CovMin,CovMax); // set the ranges accordingly, for frac cov should be [0,100], for corr matrices [-1,1]
        CovMatrix->GetZaxis()->CenterTitle();
        CovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        CovMatrix->GetZaxis()->SetTitleSize(TextSize);
        CovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        CovMatrix->GetZaxis()->SetLabelSize(TextSize);
        CovMatrix->GetZaxis()->SetNdivisions(5);

        CovMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
        CovMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        CovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/Cov"+PlotNames[i]+".png");

        // Plot frac cov matrix
        double FracCovMin = FracCovMatrix->GetMinimum();
        double FracCovMax = FracCovMatrix->GetMaximum();
        FracCovMatrix->GetZaxis()->SetRangeUser(FracCovMin,FracCovMax); // set the ranges accordingly, for frac cov should be [0,100], for corr matrices [-1,1]
        FracCovMatrix->GetZaxis()->CenterTitle();
        FracCovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        FracCovMatrix->GetZaxis()->SetTitleSize(TextSize);
        FracCovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        FracCovMatrix->GetZaxis()->SetLabelSize(TextSize);
        FracCovMatrix->GetZaxis()->SetNdivisions(5);

        FracCovMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
        FracCovMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        FracCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/FracCov"+PlotNames[i]+".png");

        // Plot correlation matrix
        CorrMatrix->GetZaxis()->SetRangeUser(-1,1);
        CorrMatrix->GetZaxis()->CenterTitle();
        CorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        CorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        CorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        CorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        CorrMatrix->GetZaxis()->SetNdivisions(5);

        CorrMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
        CorrMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        CorrMatrix->Draw("colz text");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/Corr"+PlotNames[i]+".png");

        // Save objects
        SaveFile->WriteObject(CovMatrix, PlotNames[i]+"_cov");
        SaveFile->WriteObject(FracCovMatrix, PlotNames[i]+"_fraccov");
        SaveFile->WriteObject(CorrMatrix, PlotNames[i]+"_corr");

        if (ConstructSpectra) {
            // Plot histograms with error bands, only when constructing spectra
            TGraphAsymmErrors* RecoErrorBand = RecoSpectra->ErrorBand(TargetPOT);
            TGraphAsymmErrors* RecoTrueErrorBand = RecoTrueSpectra->ErrorBand(TargetPOT);
            TGraphAsymmErrors* RecoBkgErrorBand = RecoBkgSpectra->ErrorBand(TargetPOT);

            PlotCanvas->cd();

            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            RecoHisto->Draw("hist");
            ana::DrawErrorBand(RecoHisto, RecoErrorBand);
            RecoTrueHisto->Draw("hist same");
            ana::DrawErrorBand(RecoTrueHisto, RecoTrueErrorBand);
            RecoBkgHisto->Draw("hist same");
            ana::DrawErrorBand(RecoBkgHisto, RecoBkgErrorBand);
            leg->Draw();
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/"+PlotNames[i]+".png");

            // Save objects
            SaveFile->WriteObject(RecoHisto, PlotNames[i]+"_reco");
            SaveFile->WriteObject(RecoTrueHisto, PlotNames[i]+"_reco_true");
            SaveFile->WriteObject(RecoBkgHisto, PlotNames[i]+"_reco_bkg");
        }
        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}
