// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"

// std includes.
#include <vector>

// Definitions for Vars and Cuts.
#include "Definitions.h"

using namespace std;
using namespace ana;

void Selection()
{
    // Some useful variables for later.
    const std::string TargetFile = "/exp/sbnd/data/users/munjung/SBND/2023B/cnnid/cnnid.flat.caf.root";
    const double TargetPOT(6.6e20);

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(TargetFile);

    // Side note: I find it to be good practice to add a letter to each object to denote the type.
    // It helps me keep track of the many different objects and makes the code more readable (opinion).

    // Create the binning schemes for the Vars we wish to plot.
    const Binning bPrimaryEnergy = Binning::Simple(20, 0, 3.0);
    const Binning bTrackLength = Binning::Simple(40, 0, 200);

    // Create a simple Spectrum showing the primary neutrino energy of only slices with a nu mu CC interaction.
    Spectrum sNuEnergy( "Neutrino Energy [GeV]",     // A label for the Spectrum.
                bPrimaryEnergy,              // Use 20 bins from 0.0 to 3.0 GeV
                NuLoader,                    // Associate this Spectrum with the NuLoader object (and its target CAF)
                kPrimaryEnergy,              // The Var to plot.
                kNoSpillCut,                 // The SpillCut to use (none in this case).
                kIsNuMuCC );                 // The Cut to use (only true nu mu CC slices).

    // Let's do the same as above, but add a SpillCut to see how many nu mu CC slices survive a simple CRT hit veto.
    Spectrum sNuEnergyNoCRT( "Neutrino Energy [GeV]",     // A label for the Spectrum.
                bPrimaryEnergy,              // Use 20 bins from 0.0 to 3.0 GeV
                NuLoader,                    // Associate this Spectrum with the NuLoader object (and its target CAF)
                kPrimaryEnergy,              // The Var to plot.
                kCRTHitVeto,                 // Use a simple CRT Hit SpillCut.
                kIsNuMuCC );                 // The Cut to use (only true nu mu CC slices).

    // Create a simple Spectrum showing the length of all tracks in all slices with a nu mu CC interaction.
    Spectrum sTrackLen( "Track Length [cm]",         // A label for the Spectrum.
                bTrackLength,                // Use 20 bins from 0.0 to 3.0 GeV
                NuLoader,                    // Associate this Spectrum with the NuLoader object (and its target CAF)
                kAllTrkLen,                  // The Var to plot.
                kNoSpillCut,                 // The SpillCut to use (none in this case).
                kIsNuMuCC );                 // The Cut to use (only true nu mu CC slices).

    // Create a simple Spectrum showing the primary neutrino energy of MC events that match signal definition.
    Spectrum sNuEnergySignal( "True Neutrino Energy [GeV]",
                bPrimaryEnergy, // Use 20 bins from 0.0 to 3.0 GeV
                NuLoader,       // Associate this Spectrum with the NuLoader object (and its target CAF)
                kTrueEnergy,    // The TruthVar to plot
                kIsSignal,      // The TruthCut to use (signal definition)
                kNoSpillCut     // The SpillCut to use (none)
            );  

    // Create a simple Spectrum showing the primary neutrino energy of slices truth matched to signal definition.
    Spectrum sNuEnergySignalTruthMatched( "True Neutrino Energy [GeV]",
                bPrimaryEnergy, // Use 20 bins from 0.0 to 3.0 GeV
                NuLoader,       // Associate this Spectrum with the NuLoader object (and its target CAF)
                kTrueEnergy,    // The TruthVar to plot
                kIsSignal,      // The TruthCut to use (signal definition)
                kNoSpillCut,    // The SpillCut to use (none)
                kIsNuMuCC       // The Cut to use (true nu mu CC slices)
            );  

    // Now that each Spectrum is defined, use the Go() method to populate the Spectrum objects.
    NuLoader.Go();

    // Write out the Spectrum objects to a TCanvas
    std::vector<Spectrum> Spectra{
        sNuEnergy,
        sNuEnergyNoCRT,
        sTrackLen,
    };
    std::vector<TString> PlotNames{
        "sNuEnergy",
        "sNuEnergyNoCRT",
        "sTrackLen",
    };
    const int nSpectra = Spectra.size();

    TString dir = "/exp/sbnd/app/users/epelaez/CC1muAnalysis";
    for (int i = 0; i < nSpectra; i++) {
        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1D* Histogram = Spectra[i].ToTH1(TargetPOT);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        Histogram->Draw("hist same");

        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/"+PlotNames[i]+".png");
        delete PlotCanvas;
    }

    // We want to plot some Spectrum objects overlaid with each other
    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
    TH1D* Histo1 = sNuEnergySignal.ToTH1(TargetPOT);
    TH1D* Histo2 = sNuEnergySignalTruthMatched.ToTH1(TargetPOT);

    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.17);
    PlotCanvas->SetRightMargin(0.05);
    PlotCanvas->SetBottomMargin(0.16);

    TLegend* leg = new TLegend(0.13,0.88,0.95,0.99);
    leg->SetBorderSize(0);
    leg->SetNColumns(3);
    leg->SetMargin(0.2);
    leg->SetFillColor(0);

    TLegendEntry* legColor1 = leg->AddEntry(Histo1,"all","l");
    legColor1->SetTextColor(602); // blue
    Histo1->SetLineColor(602);
    TLegendEntry* legColor2 = leg->AddEntry(Histo2,"truth matched","l");
    legColor2->SetTextColor(797); // orange
    Histo2->SetLineColor(797);  

    PlotCanvas->cd();
    Histo1->Draw("hist same");
    Histo2->Draw("hist same");
    leg->Draw();
    PlotCanvas->SaveAs(dir+"/Figs/CAFAna/EnergyTruthMatch.png");
    delete PlotCanvas;
}