#ifndef Constants_h
#define Constants_h

#include "TString.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <tuple>

using namespace std;

namespace Constants {
    const double Units = 1E38;

    const double TargetPOT(6.6e20);
    const double NTargets = 1.05E30; // Argon nuclei, not nucleons

    double Nominal_UB_XY_Surface = 175.*180.*2.*2.; // cm2
	double POTPerSpill = 4997.*5e8;

    // Variables for double differential analysis
    static const int TwoDNBinsMuonCosTheta = 2; 
    std::vector<double> TwoDArrayNBinsMuonCosTheta{-1.0,0.5,1.0};

    static const int TwoDNBinsTransverseMomentum = 11;
    std::vector<double> TwoDArrayTransverseMomentum{0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2};

    static const int TwoDNBinsDeltaAlphaT = 9;
    std::vector<double> TwoDArrayDeltaAlphaT{0.,20.,40.,60.,80.,100.,120.,140.,160.,180.};

    std::vector<std::vector<double>> TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices{
        {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2},
        {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2},
    };
    static const TString LabelXAxisTwoDTransverseMomentumInMuonCosTheta = ";#deltap_{T} [bin #]";

    std::vector<std::vector<double>> TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices{
        {0.,20.,40.,60.,80.,100.,120.,140.,160.,180.},
        {0.,20.,40.,60.,80.,100.,120.,140.,160.,180.},
    };
    static const TString LabelXAxisTwoDDeltaAlphaTInMuonCosTheta = ";#delta #alpha_{T} [bin #]";

    std::vector<std::vector<double>> TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices{
        {-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.},
        {-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.},
    };
    static const TString LabelXAxisTwoDCosOpeningAngleProtonsInMuonCosTheta = ";cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) [bin #]";

    std::vector<std::vector<double>> TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices{
        {-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.},
        {-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.},
    };
    static const TString LabelXAxisTwoDCosOpeningMuonTotalProtonInMuonCosTheta = ";cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) [bin #]";

    static std::map<TString,TString> LatexLabel = {
        { "MuonCosThetaPlot",  "All events" },
        { "LeadingProtonCosThetaPlot",  "All events" },	
        { "RecoilProtonCosThetaPlot",  "All events" },	
        { "LeadingProtonMomentumPlot",  "All events" },	
        { "RecoilProtonMomentumPlot",  "All events" },	
        { "MuonMomentumPlot",  "All events" },	
        { "CosOpeningAngleProtonsPlot",  "All events" },	
        { "CosOpeningAngleMuonTotalProtonPlot",  "All events" },	
        { "TransverseMomentumPlot",  "All events" },	
    };

    static std::map<TString, std::tuple<vector<double>, vector<vector<double>>>> PlotNameToDiscriminator = {
        {"TrueSerialTransverseMomentum_InMuonCosThetaPlot", {TwoDArrayNBinsMuonCosTheta, TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices}},
        {"TrueSerialDeltaAlphaT_InMuonCosThetaPlot", {TwoDArrayNBinsMuonCosTheta, TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices}},
        {"TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot", {TwoDArrayNBinsMuonCosTheta, TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices}},
        {"TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot", {TwoDArrayNBinsMuonCosTheta, TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices}}
    };

    static std::map<TString, TString> PlotNameToSliceLabel = {
        {"TrueSerialTransverseMomentum_InMuonCosThetaPlot", "cos(#theta_{#vec{p}_{#mu}})"},
        {"TrueSerialDeltaAlphaT_InMuonCosThetaPlot", "cos(#theta_{#vec{p}_{#mu}})"},
        {"TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot", "cos(#theta_{#vec{p}_{#mu}})"},
        {"TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot", "cos(#theta_{#vec{p}_{#mu}})"}
    };
}

#endif