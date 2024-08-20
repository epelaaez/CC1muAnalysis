#ifndef Constants_h
#define Constants_h

#include "TString.h"
#include "TMath.h"
#include "TH1D.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <tuple>

namespace Constants {
    // User to access
    const std::string UserName = "rickyoon";

    const double Units = 1E38;

    const double TargetPOT = 6.79e20;
    const double NTargets = 1.05E30; // Argon nuclei, not nucleons
    // const double NTargets = 4.6712e31; // Nucleons

    double Nominal_UB_XY_Surface = 175. * 180. * 2. * 2.; // cm2
	double POTPerSpill = 5e12;

    // Integrated flux
    TFile* FluxFile = TFile::Open("MCC9_FluxHist_volTPCActive.root"); // make sure file is in path
	TH1D* HistoFlux = (TH1D*)(FluxFile->Get("hEnumu_cv"));
    double IntegratedFlux = (HistoFlux->Integral() * (TargetPOT / POTPerSpill / Nominal_UB_XY_Surface));
    // double IntegratedFlux = 1.65974e13; // from Henry Lay

    // Binning for single differential analysis
    static const int NBinsEventCount = 1;
    static const std::vector<double> ArrayNBinsEventCount{0., 1.};

    static const int NBinsAngle = 10;
    static const std::vector<double> ArrayNBinsAngle{-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.};

    static const int NBinsDeltaAlphaT = 6;
    static const std::vector<double> ArrayNBinsDeltaAlphaT{0.,30.,60.,90.,120.,150.,180.};

    static const int NBinsTransverseMomentum = 7;
    static const std::vector<double> ArrayNBinsTransverseMomentum{0.,0.1,0.2,0.3,0.4,0.5,0.6,1.};

    //Additional Binning for single differential
    static const int NBinsInvariantMass = 20;
    static const std::vector<double> ArrayNBinsInvariantMass{1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.};

    static const int NBinsCosAngleLPMu = 20;
    static const std::vector<double> ArrayNBinsCosAngleLPMu = {-1.,-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

    static const int NBinsCosAngleRPMu = 20;
    static const std::vector<double> ArrayNBinsCosAngleRPMu = {-1.,-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

    // 0.183 spacing
    static const int NBinsMuonMomentum = 6;
    static const std::vector<double> ArrayNBinsMuonMomentum{0.1, 0.283, 0.466, 0.649, 0.832, 1.015, 1.2};

    // 0.116 spacing
    static const int NBinsLeadingProtonMomentum = 6;
    static const std::vector<double> ArrayNBinsLeadingProtonMomentum{0.3, 0.416, 0.532, 0.648, 0.764, 0.880, 1.};

    static const int NBinsRecoilProtonMomentum = 5;
    static const std::vector<double> ArrayNBinsRecoilProtonMomentum{0.3, 0.416, 0.532, 0.648, 0.764, 1.};

    // Variables for double differential analysis
    static const int TwoDNBinsMuonCosTheta = 2; 
    std::vector<double> TwoDArrayNBinsMuonCosTheta{-1.0,0.5,1.0};

    static const int TwoDNBinsTransverseMomentum = 5;
    std::vector<double> TwoDArrayTransverseMomentum{0.,0.2,0.4,0.6,0.8,1.};

    static const int TwoDNBinsDeltaAlphaT = 6;
    std::vector<double> TwoDArrayDeltaAlphaT{0.,30.,60.,90.,120.,150.,180.};

    std::vector<std::vector<double>> TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices{
        {0.,0.1,0.2,0.3,0.4,0.5,0.6,1.},
        {0.,0.1,0.2,0.3,0.4,0.5,0.6,1.},
    };
    static const TString LabelXAxisTwoDTransverseMomentumInMuonCosTheta = ";#deltap_{T} [bin #]";

    std::vector<std::vector<double>> TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices{
        {0.,30.,60.,90.,120.,150.,180.},
        {0.,30.,60.,90.,120.,150.,180.},
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
        { "InvariantMassPlot", "All events" },
        { "CosOpeningAngleLProtonMuonPlot", "All events" },
	{ "CosOpeningAngleRProtonMuonPlot", "All events" }	
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

    /////////////
    // Plot names
    /////////////

    static const std::vector<TString> PlotNames = {
        "EventCount",
        "MuonCosTheta",
        "LeadingProtonCosTheta",
        "RecoilProtonCosTheta",
        "CosOpeningAngleProtons",
        "CosOpeningAngleMuonTotalProton",
        "DeltaAlphaT",
        "TransverseMomentum",
        "MuonMomentum",
        "LeadingProtonMomentum",
        "RecoilProtonMomentum",
        "SerialTransverseMomentum_InMuonCosTheta",
        "SerialDeltaAlphaT_InMuonCosTheta",
        "SerialCosOpeningAngleProtons_InMuonCosTheta",
        "SerialCosOpeningAngleMuonTotalProton_InMuonCosTheta",
        "InvariantMass",
        "CosOpeningAngleLProtonMuon",
        "CosOpeningAngleRProtonMuon"
    };

    static const std::vector<std::string> VarLabels = {
        "single bin",
        "cos(#theta_{#vec{p}_{#mu}})",
        "cos(#theta_{#vec{p}_{L}})",
        "cos(#theta_{#vec{p}_{R}})",
        "cos(#theta_{#vec{p}_{L},#vec{p}_{R}})",
        "cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})",
        "#delta #alpha_{T}",
        "#delta P_{T}",
        "|#vec{p}_{#mu}|",
        "|#vec{p}_{L}|",
        "|#vec{p}_{R}|",
        "#delta P_{T} (bin #)",
        "#delta #alpha_{T} (bin #)",
        "cos(#theta_{#vec{p}_{L},#vec{p}_{R}}) (bin #)",
        "cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}}) (bin #)",
 	"W",
 	"cos(#theta_{p^{L}, p^{#mu}})",
	"cos(#theta_{p^{R}, p^{#mu}})"
    };

    static const std::vector<std::string> YLabels = {
        "#frac{d#sigma}{# events} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{dcos(#theta_{#vec{p}_{L}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{dcos(#theta_{#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{d|#vec{p}_{#mu}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{d|#vec{p}_{L}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d#sigma}{d|#vec{p}_{R}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d^{2}#sigma}{dcos(#theta_{#vec{p}_{#mu}}) d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d^{2}#sigma}{dcos(#theta_{#vec{p}_{#mu}}) d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d^{2}#sigma}{dcos(#theta_{#vec{p}_{#mu}}) dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
        "#frac{d^{2}#sigma}{dcos(#theta_{#vec{p}_{#mu}}) dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
	"#frac{d#sigma}{dW} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
	"#frac{d#sigma}{d#cos(#theta_{p^{L}, p^{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]",
	"#frac{d#sigma}{d#cos(#theta_{p^{R}, p^{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]"
    };

    ///////////////
    // Systematics
    ///////////////

    // If using most recent files
    // const std::vector<std::tuple<std::string, int>> XSecSystsVector = {
    //     {"GENIEReWeight_SBND_v1_multisigma_MaCCQE", 6}, // 0
    //     {"GENIEReWeight_SBND_v1_multisigma_MaNCEL", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_EtaNCEL", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_MaCCRES", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_MvCCRES", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_MaNCRES", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_MvNCRES", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpCC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpCC2pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpNC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpNC2pi", 6}, // 10
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnCC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnCC2pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnNC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnNC2pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpCC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpCC2pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpNC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpNC2pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnCC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnCC2pi", 6}, // 20
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnNC1pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnNC2pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_RDecBR1gamma", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_RDecBR1eta", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_Theta_Delta2Npi", 10},
    //     {"GENIEReWeight_SBND_v1_multisigma_AhtBY", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_BhtBY", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_CV1uBY", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_CV2uBY", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FormZone", 6}, // 30
    //     {"GENIEReWeight_SBND_v1_multisigma_MFP_pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrCEx_pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrInel_pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrAbs_pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrPiProd_pi", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_MFP_N", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrCEx_N", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrInel_N", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrAbs_N", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_FrPiProd_N", 6}, // 40
    //     {"GENIEReWeight_SBND_v1_multisigma_CCQEPauliSupViaKF", 6},
    //     {"GENIEReWeight_SBND_v1_multisigma_CCQEMomDistroFGtoSF", 10},
    //     {"GENIEReWeight_SBND_v1_multisim_MaCCQE", 100}, 
    //     {"GENIEReWeight_SBND_v1_multisim_MaNCEL", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_EtaNCEL", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_MaCCRES", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_MvCCRES", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_MaNCRES", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_MvNCRES", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpCC1pi", 100}, // 50
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpCC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpNC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpNC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnCC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnCC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnNC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnNC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpCC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpCC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpNC1pi", 100}, // 60
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpNC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnCC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnCC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnNC1pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnNC2pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_RDecBR1gamma", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_RDecBR1eta", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_AhtBY", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_BhtBY", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_CV1uBY", 100}, // 70
    //     {"GENIEReWeight_SBND_v1_multisim_CV2uBY", 100},
    //     // {"GENIEReWeight_SBND_v1_multisim_FormZone", 100}, 
    //     {"GENIEReWeight_SBND_v1_multisim_MFP_pi", 100}, 
    //     {"GENIEReWeight_SBND_v1_multisim_FrCEx_pi", 100}, 
    //     {"GENIEReWeight_SBND_v1_multisim_FrInel_pi", 100}, 
    //     {"GENIEReWeight_SBND_v1_multisim_FrAbs_pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_FrPiProd_pi", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_MFP_N", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_FrCEx_N", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_FrInel_N", 100}, // 80
    //     {"GENIEReWeight_SBND_v1_multisim_FrAbs_N", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_FrPiProd_N", 100},
    //     {"GENIEReWeight_SBND_v1_multisim_CCQEPauliSupViaKF", 100},
    //     {"MINERvAE2p2h_ICARUS_v1_E2p2h_A_nu", 6},
    //     {"MINERvAE2p2h_ICARUS_v1_E2p2h_B_nu", 6},
    //     {"MINERvAE2p2h_ICARUS_v1_E2p2h_A_nubar", 6},
    //     {"MINERvAE2p2h_ICARUS_v1_E2p2h_B_nubar", 6},
    //     // {"MINERvAq0q3Weighting_SBND_v1_Mnv2p2hGaussEnhancement", 4},
    //     {"MiscInteractionSysts_SBND_v1_C12ToAr40_2p2hScaling_nu", 6},
    //     {"MiscInteractionSysts_SBND_v1_C12ToAr40_2p2hScaling_nubar", 6}, // 90
    //     // {"MiscInteractionSysts_SBND_v1_nuenuebar_xsec_ratio", 2},
    //     // {"MiscInteractionSysts_SBND_v1_nuenumu_xsec_ratio", 2},
    //     {"MiscInteractionSysts_SBND_v1_SPPLowQ2Suppression", 10},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_CC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_CC_3Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_CC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_CC_3Pi", 6},
    //     // {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_np_CC_1Pi", 7},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_NC_1Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_NC_2Pi", 6}, // 100
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_NC_3Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_NC_1Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_NC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_NC_3Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_CC_1Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_CC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_CC_3Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_CC_1Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_CC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_CC_3Pi", 6}, // 110
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_NC_1Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_NC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_NC_3Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_NC_1Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_NC_2Pi", 6},
    //     {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_NC_3Pi", 6}	
    // };

    // If using not so recent files
    const std::vector<std::tuple<std::string, int>> XSecSystsVector = {
        {"AhtBY_multisigma_Genie", 6}, // 0
        {"BhtBY_multisigma_Genie", 6},
        {"CV1uBY_multisigma_Genie", 6},
        {"CV2uBY_multisigma_Genie", 6},
        {"EtaNCEL_multisigma_Genie", 6},
        {"FormZone_multisigma_Genie", 6},
        {"FrAbs_N_multisigma_Genie", 6},
        {"FrAbs_pi_multisigma_Genie", 6},
        {"FrCEx_N_multisigma_Genie", 6},
        {"FrCEx_pi_multisigma_Genie", 6},
        {"FrInel_N_multisigma_Genie", 6}, // 10
        {"FrInel_pi_multisigma_Genie", 6},
        {"FrPiProd_N_multisigma_Genie", 6},
        {"FrPiProd_pi_multisigma_Genie", 6},
        {"MFP_N_multisigma_Genie", 6},
        {"MFP_pi_multisigma_Genie", 6},
        {"MaCCQE_multisigma_Genie", 6},
        {"MaCCRES_multisigma_Genie", 6},
        {"MaNCEL_multisigma_Genie", 6},
        {"MaNCRES_multisigma_Genie", 6},
        {"MvCCRES_multisigma_Genie", 6}, // 20
        {"MvNCRES_multisigma_Genie", 6},
        {"NonRESBGvbarnCC1pi_multisigma_Genie", 6},
        {"NonRESBGvbarnCC2pi_multisigma_Genie", 6},
        {"NonRESBGvbarnNC1pi_multisigma_Genie", 6},
        {"NonRESBGvbarnNC2pi_multisigma_Genie", 6},
        {"NonRESBGvbarpCC1pi_multisigma_Genie", 6},
        {"NonRESBGvbarpCC2pi_multisigma_Genie", 6},
        {"NonRESBGvbarpNC1pi_multisigma_Genie", 6},
        {"NonRESBGvbarpNC2pi_multisigma_Genie", 6},
        {"NonRESBGvnCC1pi_multisigma_Genie", 6}, // 30
        {"NonRESBGvnCC2pi_multisigma_Genie", 6}, 
        {"NonRESBGvnNC1pi_multisigma_Genie", 6},
        {"NonRESBGvnNC2pi_multisigma_Genie", 6},
        {"NonRESBGvpCC1pi_multisigma_Genie", 6},
        {"NonRESBGvpCC2pi_multisigma_Genie", 6},
        {"NonRESBGvpNC1pi_multisigma_Genie", 6},
        {"NonRESBGvpNC2pi_multisigma_Genie", 6},
        // {"multisim_Genie", 1000}
    };

    const std::vector<std::tuple<std::string, int>> FluxSystsVector = {
        {"expskin_Flux", 1000},
        {"horncurrent_Flux", 1000},
        {"kminus_Flux", 1000},
        {"kplus_Flux", 1000},
        {"kzero_Flux", 1000},
        {"nucleoninexsec_Flux", 1000},
        {"nucleonqexsec_Flux", 1000},
        {"nucleontotxsec_Flux", 1000},
        {"piminus_Flux", 1000},
        {"pioninexsec_Flux", 1000},
        {"pionqexsec_Flux", 1000},
        {"piontotxsec_Flux", 1000},
        {"piplus_Flux", 1000}
    };

    ////////////
    // Fake data
    ////////////

    const std::vector<TString> FakeDataNames = {
        "TwiceMEC",
        "TwiceQE",
        "Combined"
    };
}

#endif
