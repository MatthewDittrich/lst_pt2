#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include <TH1F.h>
#include <TH2D.h>

class HistogramManager {
public:
    HistogramManager() = default;
    ~HistogramManager() = default;

    void init();
    void write();

    std::vector<std::string> catNames = {
        "Barrel_L1_Flat", "Barrel_L1_Tilted",
        "Barrel_L2_Flat", "Barrel_L2_Tilted",
        "Barrel_L3plus",  "Endcap_L1",
        "Endcap_L2",      "Endcap_L3plus"
    };
    std::vector<std::string> catTitles = {
        "Barrel L1 Flat", "Barrel L1 Tilted",
        "Barrel L2 Flat", "Barrel L2 Tilted",
        "Barrel L3+",     "Endcap L1",
        "Endcap L2",      "Endcap L3+"
    };

    // Kinematics & Angles
    TH1D* real_pt2_deltaPT[8];
    TH1D* real_pt2_deltaETA[8];
    TH1D* real_pt2_deltaPHI[8];
    TH1D* real_pt2_deltaR[8];
    TH1D* real_pt2_deltaAngle[8];

    TH1D* fake_pt2_deltaPT[8];
    TH1D* fake_pt2_deltaETA[8];
    TH1D* fake_pt2_deltaPHI[8];
    TH1D* fake_pt2_deltaR[8];
    TH1D* fake_pt2_deltaAngle[8];

    TH1D* real_unused_pt2_deltaPT[8];
    TH1D* real_unused_pt2_deltaETA[8];
    TH1D* real_unused_pt2_deltaPHI[8];
    TH1D* real_unused_pt2_deltaR[8];
    TH1D* real_unused_pt2_deltaAngle[8];

    TH1D* fake_unused_pt2_deltaPT[8];
    TH1D* fake_unused_pt2_deltaETA[8];
    TH1D* fake_unused_pt2_deltaPHI[8];
    TH1D* fake_unused_pt2_deltaR[8];
    TH1D* fake_unused_pt2_deltaAngle[8];

    // LST specific variables
    TH1D* real_pt2_LSTdPhi[8];
    TH1D* fake_pt2_LSTdPhi[8];
    TH1D* real_unused_pt2_LSTdPhi[8];
    TH1D* fake_unused_pt2_LSTdPhi[8];

    TH1D* real_pt2_LSTdBeta[8];
    TH1D* fake_pt2_LSTdBeta[8];
    TH1D* real_unused_pt2_LSTdBeta[8];
    TH1D* fake_unused_pt2_LSTdBeta[8];

    TH1D* real_pt2_LSTbetaOut[8];
    TH1D* fake_pt2_LSTbetaOut[8];
    TH1D* real_unused_pt2_LSTbetaOut[8];
    TH1D* fake_unused_pt2_LSTbetaOut[8];

    TH1D* real_pt2_LSTKinZRes[8];
    TH1D* fake_pt2_LSTKinZRes[8];
    TH1D* real_unused_pt2_LSTKinZRes[8];
    TH1D* fake_unused_pt2_LSTKinZRes[8];

    TH1D* real_pt2_LSTOrgZRes[8];
    TH1D* fake_pt2_LSTOrgZRes[8];
    TH1D* real_unused_pt2_LSTOrgZRes[8];
    TH1D* fake_unused_pt2_LSTOrgZRes[8];

    // MD0 Components
    TH1D* real_pt2_MD0_dXY[8];
    TH1D* real_pt2_MD0_dZ[8];
    TH1D* fake_pt2_MD0_dXY[8];
    TH1D* fake_pt2_MD0_dZ[8];

    TH1D* real_unused_pt2_MD0_dXY[8];
    TH1D* real_unused_pt2_MD0_dZ[8];
    TH1D* fake_unused_pt2_MD0_dXY[8];
    TH1D* fake_unused_pt2_MD0_dZ[8];

    // MD1 Components
    TH1D* real_pt2_MD1_dXY[8];
    TH1D* real_pt2_MD1_dZ[8];
    TH1D* fake_pt2_MD1_dXY[8];
    TH1D* fake_pt2_MD1_dZ[8];

    TH1D* real_unused_pt2_MD1_dXY[8];
    TH1D* real_unused_pt2_MD1_dZ[8];
    TH1D* fake_unused_pt2_MD1_dXY[8];
    TH1D* fake_unused_pt2_MD1_dZ[8];

    // Separated R-Z Simple Pointing
    TH1D* real_pt2_MD0_rz_simple[8];
    TH1D* real_pt2_MD1_rz_simple[8];
    TH1D* fake_pt2_MD0_rz_simple[8];
    TH1D* fake_pt2_MD1_rz_simple[8];
    
    TH1D* real_unused_pt2_MD0_rz_simple[8];
    TH1D* real_unused_pt2_MD1_rz_simple[8];
    TH1D* fake_unused_pt2_MD0_rz_simple[8];
    TH1D* fake_unused_pt2_MD1_rz_simple[8];

    // --- 2D Cut Study Histograms ---
    TH2D* h2_real_MD0_dXY_vs_dZ;
    TH2D* h2_fake_MD0_dXY_vs_dZ;
    TH2D* h2_real_MD1_dXY_vs_dZ;
    TH2D* h2_fake_MD1_dXY_vs_dZ;

    TH2D* h2_real_unused_MD0_dXY_vs_dZ;
    TH2D* h2_fake_unused_MD0_dXY_vs_dZ;
    TH2D* h2_real_unused_MD1_dXY_vs_dZ;
    TH2D* h2_fake_unused_MD1_dXY_vs_dZ;
};

#endif
