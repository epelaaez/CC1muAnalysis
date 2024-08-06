// SBNAna includes.
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"

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
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Constants.h"
#include "../../Utils/Tools.cxx"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionMCStatSystematics() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // Initialize tools
    Tools tools;

    // Get integrated flux
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root");
    TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
    double IntegratedFlux = (HistoFlux->Integral() * TargetPOT / POTPerSpill / Nominal_UB_XY_Surface);

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/Uncertainties/MCStat");

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsMCStat.root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");

    // Vectors to fill with variables and variable information to plot
    std::vector<Var> Vars; std::vector<Binning> VarBins;
    std::vector<TString> PlotNames; std::vector<std::string> VarLabels;

    ////////////////////////////////
    // Single differential variables
    ////////////////////////////////

    // Event count 
    Vars.push_back(kEventCount); VarBins.push_back(bEventCount);
    PlotNames.push_back("EventCount"); VarLabels.push_back("single bin");

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
    Vars.push_back(kLeadingProtonMomentum); VarBins.push_back(bLeadingProtonMomentumBins);
    PlotNames.push_back("LeadingProtonMomentum"); VarLabels.push_back("|#vec{p}_{L}|");

    // Recoil proton momentum 
    Vars.push_back(kRecoilProtonMomentum); VarBins.push_back(bRecoilProtonMomentumBins);
    PlotNames.push_back("RecoilProtonMomentum"); VarLabels.push_back("|#vec{p}_{R}|");

    //////////////////////////////
    // Double differential variables
    //////////////////////////////

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

    const int NVars = Vars.size();

    // Create ISyst that computes weight for each event
    class MCStatSyst : public ISyst {
        public:
            MCStatSyst(int UnivIndex) : ISyst("MCStat", "MCStat"), UnivIndex(UnivIndex) {}

            // Use time to create random seed
            void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const {
                int seed = (int) slc->truth.time * 1000;
                weight *= tools.PoissonRandomNumber(seed);
            }

            void Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const {
                int seed = (int) nu->time * 1000;
                weight *= tools.PoissonRandomNumber(seed);
            }

            const Tools tools;
            int UnivIndex;
    };

    std::vector<SystShifts> Shifts; 
    for (int i = 0; i < 100; ++i) {
        ISyst* Syst = new MCStatSyst(i);
        SystShifts Shift(Syst, +1);
        Shifts.push_back(Shift);
    }

    // Construct all spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>
    >> Spectra;
    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, Vars.at(iVar), kNoSpillCut, kRecoIsSignal); 
        auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(iVar), VarBins.at(iVar), NuLoader, Vars.at(iVar), kNoSpillCut, kRecoIsTrueReco); 
        auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, Vars.at(iVar), kNoSpillCut, kRecoIsBackground); 

        auto UnivRecoSignals = std::make_unique<EnsembleSpectrum>(
            NuLoader,
            HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)), 
            kNoSpillCut,
            kRecoIsSignal, 
            Shifts
        );
        auto UnivRecoTrueSignals = std::make_unique<EnsembleSpectrum>(
            NuLoader,
            HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)), 
            kNoSpillCut,
            kRecoIsTrueReco, 
            Shifts
        ); 
        auto UnivRecoBkgSignals = std::make_unique<EnsembleSpectrum>(
            NuLoader,
            HistAxis(VarLabels.at(iVar), VarBins.at(iVar), Vars.at(iVar)), 
            kNoSpillCut,
            kRecoIsBackground, 
            Shifts
        ); 
        Spectra.push_back({
            std::move(RecoSignals), 
            std::move(RecoTrueSignals), 
            std::move(RecoBkgSignals),
            std::move(UnivRecoSignals), 
            std::move(UnivRecoTrueSignals), 
            std::move(UnivRecoBkgSignals)
        });
    }
    NuLoader.Go();

    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto& [
            RecoSpectra, 
            RecoTrueSpectra, 
            RecoBkgSpectra,
            UnivRecoSpectra, 
            UnivRecoTrueSpectra, 
            UnivRecoBkgSpectra
        ] = Spectra.at(iVar);

        // Get nominal plots
        TH1D* RecoHisto = RecoSpectra->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSpectra->ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = RecoBkgSpectra->ToTH1(TargetPOT);

        // Get universe plots
        std::vector<TH1D*> UnivRecoHistos;
        std::vector<TH1D*> UnivRecoTrueHistos;
        std::vector<TH1D*> UnivRecoBkgHistos;
        for (int i = 0; i < 100; ++i) {
            TH1D* ModifiedRecoHisto = UnivRecoSpectra->Universe(i).ToTH1(TargetPOT);
            UnivRecoHistos.push_back(std::move(ModifiedRecoHisto));
            TH1D* ModifiedRecoTrueHisto = UnivRecoTrueSpectra->Universe(i).ToTH1(TargetPOT);
            UnivRecoTrueHistos.push_back(std::move(ModifiedRecoTrueHisto));
            TH1D* ModifiedRecoBkgHisto = UnivRecoBkgSpectra->Universe(i).ToTH1(TargetPOT);
            UnivRecoBkgHistos.push_back(std::move(ModifiedRecoBkgHisto));
        }

        // Plot histograms with error band
        SelectionHelpers::DrawHistosWithErrorBands(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            UnivRecoHistos,
            UnivRecoTrueHistos,
            UnivRecoBkgHistos,
            dir, 
            "MCStat",
            PlotNames[iVar]
        );

        // Scale histograms for cov matrices
        RecoHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

        // Scale modified histograms for cov matrices
        for (int i = 0; i < 100; ++i) {
            UnivRecoHistos[i]->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoTrueHistos[i]->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoBkgHistos[i]->Scale(Units / (IntegratedFlux * NTargets));
        }

        // Plot cov, frac cov, and corr matrices
        SelectionHelpers::DrawMatrices(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            UnivRecoHistos,
            UnivRecoTrueHistos,
            UnivRecoBkgHistos,
            dir,
            "MCStat",
            VarLabels[iVar],
            PlotNames[iVar],
            SaveFile
        );
    }
    // Close file
    SaveFile->Close();
}