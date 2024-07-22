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

    static const int TwoDNBinsTransverseMomentum = 10;
    std::vector<double> TwoDArrayTransverseMomentum{0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.};

    static const int TwoDNBinsDeltaAlphaT = 9;
    std::vector<double> TwoDArrayDeltaAlphaT{0.,20.,40.,60.,80.,100.,120.,140.,160.,180.};

    std::vector<std::vector<double>> TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices{
        {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.},
        {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.},
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

    ///////////////
    // Systematics
    ///////////////

    const std::vector<std::tuple<std::string, int>> XSecSystsVector = {
        {"GENIEReWeight_SBND_v1_multisigma_MaCCQE", 6}, // 0
        {"GENIEReWeight_SBND_v1_multisigma_MaNCEL", 6},
        {"GENIEReWeight_SBND_v1_multisigma_EtaNCEL", 6},
        {"GENIEReWeight_SBND_v1_multisigma_MaCCRES", 6},
        {"GENIEReWeight_SBND_v1_multisigma_MvCCRES", 6},
        {"GENIEReWeight_SBND_v1_multisigma_MaNCRES", 6},
        {"GENIEReWeight_SBND_v1_multisigma_MvNCRES", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpCC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpCC2pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpNC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvpNC2pi", 6}, // 10
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnCC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnCC2pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnNC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvnNC2pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpCC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpCC2pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpNC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarpNC2pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnCC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnCC2pi", 6}, // 20
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnNC1pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_NonRESBGvbarnNC2pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_RDecBR1gamma", 6},
        {"GENIEReWeight_SBND_v1_multisigma_RDecBR1eta", 6},
        {"GENIEReWeight_SBND_v1_multisigma_Theta_Delta2Npi", 10},
        {"GENIEReWeight_SBND_v1_multisigma_AhtBY", 6},
        {"GENIEReWeight_SBND_v1_multisigma_BhtBY", 6},
        {"GENIEReWeight_SBND_v1_multisigma_CV1uBY", 6},
        {"GENIEReWeight_SBND_v1_multisigma_CV2uBY", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FormZone", 6}, // 30
        {"GENIEReWeight_SBND_v1_multisigma_MFP_pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrCEx_pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrInel_pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrAbs_pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrPiProd_pi", 6},
        {"GENIEReWeight_SBND_v1_multisigma_MFP_N", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrCEx_N", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrInel_N", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrAbs_N", 6},
        {"GENIEReWeight_SBND_v1_multisigma_FrPiProd_N", 6}, // 40
        {"GENIEReWeight_SBND_v1_multisigma_CCQEPauliSupViaKF", 6},
        {"GENIEReWeight_SBND_v1_multisigma_CCQEMomDistroFGtoSF", 10},
        {"GENIEReWeight_SBND_v1_multisim_MaCCQE", 100}, 
        {"GENIEReWeight_SBND_v1_multisim_MaNCEL", 100},
        {"GENIEReWeight_SBND_v1_multisim_EtaNCEL", 100},
        {"GENIEReWeight_SBND_v1_multisim_MaCCRES", 100},
        {"GENIEReWeight_SBND_v1_multisim_MvCCRES", 100},
        {"GENIEReWeight_SBND_v1_multisim_MaNCRES", 100},
        {"GENIEReWeight_SBND_v1_multisim_MvNCRES", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpCC1pi", 100}, // 50
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpCC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpNC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvpNC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnCC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnCC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnNC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvnNC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpCC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpCC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpNC1pi", 100}, // 60
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarpNC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnCC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnCC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnNC1pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_NonRESBGvbarnNC2pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_RDecBR1gamma", 100},
        {"GENIEReWeight_SBND_v1_multisim_RDecBR1eta", 100},
        {"GENIEReWeight_SBND_v1_multisim_AhtBY", 100},
        {"GENIEReWeight_SBND_v1_multisim_BhtBY", 100},
        {"GENIEReWeight_SBND_v1_multisim_CV1uBY", 100}, // 70
        {"GENIEReWeight_SBND_v1_multisim_CV2uBY", 100},
        // {"GENIEReWeight_SBND_v1_multisim_FormZone", 100}, 
        {"GENIEReWeight_SBND_v1_multisim_MFP_pi", 100}, 
        {"GENIEReWeight_SBND_v1_multisim_FrCEx_pi", 100}, 
        {"GENIEReWeight_SBND_v1_multisim_FrInel_pi", 100}, 
        {"GENIEReWeight_SBND_v1_multisim_FrAbs_pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_FrPiProd_pi", 100},
        {"GENIEReWeight_SBND_v1_multisim_MFP_N", 100},
        {"GENIEReWeight_SBND_v1_multisim_FrCEx_N", 100},
        {"GENIEReWeight_SBND_v1_multisim_FrInel_N", 100}, // 80
        {"GENIEReWeight_SBND_v1_multisim_FrAbs_N", 100},
        {"GENIEReWeight_SBND_v1_multisim_FrPiProd_N", 100},
        {"GENIEReWeight_SBND_v1_multisim_CCQEPauliSupViaKF", 100},
        {"MINERvAE2p2h_ICARUS_v1_E2p2h_A_nu", 6},
        {"MINERvAE2p2h_ICARUS_v1_E2p2h_B_nu", 6},
        {"MINERvAE2p2h_ICARUS_v1_E2p2h_A_nubar", 6},
        {"MINERvAE2p2h_ICARUS_v1_E2p2h_B_nubar", 6},
        // {"MINERvAq0q3Weighting_SBND_v1_Mnv2p2hGaussEnhancement", 4},
        {"MiscInteractionSysts_SBND_v1_C12ToAr40_2p2hScaling_nu", 6},
        {"MiscInteractionSysts_SBND_v1_C12ToAr40_2p2hScaling_nubar", 6}, // 90
        // {"MiscInteractionSysts_SBND_v1_nuenuebar_xsec_ratio", 2},
        // {"MiscInteractionSysts_SBND_v1_nuenumu_xsec_ratio", 2},
        {"MiscInteractionSysts_SBND_v1_SPPLowQ2Suppression", 10},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_CC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_CC_3Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_CC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_CC_3Pi", 6},
        // {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_np_CC_1Pi", 7},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_NC_1Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_NC_2Pi", 6}, // 100
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_n_NC_3Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_NC_1Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_NC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nu_p_NC_3Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_CC_1Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_CC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_CC_3Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_CC_1Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_CC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_CC_3Pi", 6}, // 110
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_NC_1Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_NC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_n_NC_3Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_NC_1Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_NC_2Pi", 6},
        {"NOvAStyleNonResPionNorm_SBND_v1_NR_nubar_p_NC_3Pi", 6}	
    };

    const std::vector<std::tuple<std::string, int>> FluxSystsVector = {
        {"expskin_Flux", 100},
        {"horncurrent_Flux", 100},
        {"kminus_Flux", 100},
        {"kplus_Flux", 100},
        {"kzero_Flux", 100},
        {"nucleoninexsec_Flux", 100},
        {"nucleonqexsec_Flux", 100},
        {"nucleontotxsec_Flux", 100},
        {"piminus_Flux", 100},
        {"pioninexsec_Flux", 100},
        {"pionqexsec_Flux", 100},
        {"piontotxsec_Flux", 100},
        {"piplus_Flux", 100}
    };
}

#endif