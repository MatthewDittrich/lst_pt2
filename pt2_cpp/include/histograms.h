#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include <TH1D.h>
#include <vector>
#include <string>

class HistogramManager {
public:
    HistogramManager() = default;
    ~HistogramManager() = default;

    void init();
    void write();

    // --- The 11 Valid Connections ---
    std::vector<std::string> catNames = {
        "L1F_to_L2F", "L1F_to_L2T", "L1T_to_L2F", "L1T_to_L2T", "L1T_to_E1",
        "L2F_to_L3P", "L2T_to_L3P", "L2T_to_E1",  "L3P_to_E1",  "E1_to_E2", "E2_to_E3P"
    };
    std::vector<std::string> catTitles = {
        "Bar L1 Flat -> Bar L2 Flat", "Bar L1 Flat -> Bar L2 Tilt", 
        "Bar L1 Tilt -> Bar L2 Flat", "Bar L1 Tilt -> Bar L2 Tilt", "Bar L1 Tilt -> Endcap L1",
        "Bar L2 Flat -> Bar L3+",     "Bar L2 Tilt -> Bar L3+",     "Bar L2 Tilt -> Endcap L1", 
        "Bar L3+ -> Endcap L1",       "Endcap L1 -> Endcap L2",     "Endcap L2 -> Endcap L3+"
    };

    // =========================================================================
    // GENERAL VARIABLES (Change all of these to [11])
    // =========================================================================
    TH1D* real_pt2_deltaPT[11];
    TH1D* real_pt2_deltaETA[11];
    TH1D* real_pt2_deltaPHI[11];
    TH1D* real_pt2_deltaR[11];
    TH1D* real_pt2_deltaAngle[11];

    TH1D* fake_pt2_deltaPT[11];
    TH1D* fake_pt2_deltaETA[11];
    TH1D* fake_pt2_deltaPHI[11];
    TH1D* fake_pt2_deltaR[11];
    TH1D* fake_pt2_deltaAngle[11];

    TH1D* real_unused_pt2_deltaPT[11];
    TH1D* real_unused_pt2_deltaETA[11];
    TH1D* real_unused_pt2_deltaPHI[11];
    TH1D* real_unused_pt2_deltaR[11];
    TH1D* real_unused_pt2_deltaAngle[11];

    TH1D* fake_unused_pt2_deltaPT[11];
    TH1D* fake_unused_pt2_deltaETA[11];
    TH1D* fake_unused_pt2_deltaPHI[11];
    TH1D* fake_unused_pt2_deltaR[11];
    TH1D* fake_unused_pt2_deltaAngle[11];

    // LST specific variables
    TH1D* real_pt2_LSTdPhi[11];
    TH1D* fake_pt2_LSTdPhi[11];
    TH1D* real_unused_pt2_LSTdPhi[11];
    TH1D* fake_unused_pt2_LSTdPhi[11];

    TH1D* real_pt2_LSTdBeta[11];
    TH1D* fake_pt2_LSTdBeta[11];
    TH1D* real_unused_pt2_LSTdBeta[11];
    TH1D* fake_unused_pt2_LSTdBeta[11];

    TH1D* real_pt2_LSTbetaOut[11];
    TH1D* fake_pt2_LSTbetaOut[11];
    TH1D* real_unused_pt2_LSTbetaOut[11];
    TH1D* fake_unused_pt2_LSTbetaOut[11];

    TH1D* real_pt2_LSTKinZRes[11];
    TH1D* fake_pt2_LSTKinZRes[11];
    TH1D* real_unused_pt2_LSTKinZRes[11];
    TH1D* fake_unused_pt2_LSTKinZRes[11];

    TH1D* real_pt2_LSTOrgZRes[11];
    TH1D* fake_pt2_LSTOrgZRes[11];
    TH1D* real_unused_pt2_LSTOrgZRes[11];
    TH1D* fake_unused_pt2_LSTOrgZRes[11];

    // MD0 Components
    TH1D* real_pt2_MD0_dXY[11];
    TH1D* real_pt2_MD0_dZ[11];
    TH1D* fake_pt2_MD0_dXY[11];
    TH1D* fake_pt2_MD0_dZ[11];

    TH1D* real_unused_pt2_MD0_dXY[11];
    TH1D* real_unused_pt2_MD0_dZ[11];
    TH1D* fake_unused_pt2_MD0_dXY[11];
    TH1D* fake_unused_pt2_MD0_dZ[11];

    // MD1 Components
    TH1D* real_pt2_MD1_dXY[11];
    TH1D* real_pt2_MD1_dZ[11];
    TH1D* fake_pt2_MD1_dXY[11];
    TH1D* fake_pt2_MD1_dZ[11];

    TH1D* real_unused_pt2_MD1_dXY[11];
    TH1D* real_unused_pt2_MD1_dZ[11];
    TH1D* fake_unused_pt2_MD1_dXY[11];
    TH1D* fake_unused_pt2_MD1_dZ[11];

    // Separated R-Z Simple Pointing
    TH1D* real_pt2_MD0_rz_simple[11];
    TH1D* real_pt2_MD1_rz_simple[11];
    TH1D* fake_pt2_MD0_rz_simple[11];
    TH1D* fake_pt2_MD1_rz_simple[11];
    
    TH1D* real_unused_pt2_MD0_rz_simple[11];
    TH1D* real_unused_pt2_MD1_rz_simple[11];
    TH1D* fake_unused_pt2_MD0_rz_simple[11];
    TH1D* fake_unused_pt2_MD1_rz_simple[11];
};

#endif
