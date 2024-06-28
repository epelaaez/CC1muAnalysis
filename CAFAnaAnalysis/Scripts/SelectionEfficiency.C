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

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionEfficiency() {
    // Some useful variables for later.
    const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const double TargetPOT(6.6e20);

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // Create the binning schemes for the Vars we wish to plot.
    const Binning bPrimaryEnergy = Binning::Simple(1, 0, 3.0); // one bin
    const Binning bAngleBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bDeltaAlphaBins = Binning::Simple(20, 0.0, 180.0);
    const Binning bTransverseMomentumBins = Binning::Simple(20, 0.0, 1.0);
    const Binning bMuonMomentumBins = Binning::Simple(20, 0.1, 1.2);
    const Binning bProtonMomentumBins = Binning::Simple(20, 0.3, 1.0);

    // Double differential bins
    Tools tools; // tools for double differential bins

    const Binning bTransverseMomentumInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices)
    );
    const Binning bDeltaAlphaTInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)
    );
    const Binning bCosOpeningAngleProtonsInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices)
    );
    const Binning bCosOpeningAngleMuonTotalProtonInMuonCosTheta = Binning::Custom(
        tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices)
    );
}