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
#include "TROOT.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>
#include <string>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Generator analysis includes.
#include "../../GeneratorAnalysis/Scripts/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void RunAllSystematics() {
    gROOT->ProcessLine(".L ./Scripts/SelectionSystematics.C");

    int nSysts = SystsVector.size();
    for (int i = 93; i < nSysts; i++) {
        gROOT->ProcessLine(("SelectionSystematics(" + std::to_string(i) + ")").c_str()); 
    }
}
