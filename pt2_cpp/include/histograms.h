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

    // Delta pT2 Hists
    TH1D* real_pt2_deltaPT  = nullptr;
    TH1D* real_pt2_deltaETA = nullptr;
    TH1D* real_pt2_deltaPHI = nullptr;

    TH1D* real_unused_pt2_deltaPT  = nullptr;
    TH1D* real_unused_pt2_deltaETA = nullptr;
    TH1D* real_unused_pt2_deltaPHI = nullptr;

    TH1D* fake_pt2_deltaPT  = nullptr;
    TH1D* fake_pt2_deltaETA = nullptr;
    TH1D* fake_pt2_deltaPHI = nullptr;

    TH1D* fake_unused_pt2_deltaPT  = nullptr;
    TH1D* fake_unused_pt2_deltaETA = nullptr;
    TH1D* fake_unused_pt2_deltaPHI = nullptr;

    // --- 3D Helical Components (MD0) ---
    TH1D* real_pt2_MD0_dXY = nullptr;
    TH1D* real_pt2_MD0_dZ  = nullptr;
    TH1D* fake_pt2_MD0_dXY = nullptr;
    TH1D* fake_pt2_MD0_dZ  = nullptr;

    TH1D* real_unused_pt2_MD0_dXY = nullptr;
    TH1D* real_unused_pt2_MD0_dZ  = nullptr;
    TH1D* fake_unused_pt2_MD0_dXY = nullptr;
    TH1D* fake_unused_pt2_MD0_dZ  = nullptr;

    // --- 3D Helical Components (MD1) ---
    TH1D* real_pt2_MD1_dXY = nullptr;
    TH1D* real_pt2_MD1_dZ  = nullptr;
    TH1D* fake_pt2_MD1_dXY = nullptr;
    TH1D* fake_pt2_MD1_dZ  = nullptr;

    TH1D* real_unused_pt2_MD1_dXY = nullptr;
    TH1D* real_unused_pt2_MD1_dZ  = nullptr;
    TH1D* fake_unused_pt2_MD1_dXY = nullptr;
    TH1D* fake_unused_pt2_MD1_dZ  = nullptr;

    // --- Separated R-Z Simple Pointing ---
    TH1D* real_pt2_MD0_rz_simple = nullptr;
    TH1D* real_pt2_MD1_rz_simple = nullptr;
    TH1D* fake_pt2_MD0_rz_simple = nullptr;
    TH1D* fake_pt2_MD1_rz_simple = nullptr;
    TH1D* real_unused_pt2_MD0_rz_simple = nullptr;
    TH1D* real_unused_pt2_MD1_rz_simple = nullptr;
    TH1D* fake_unused_pt2_MD0_rz_simple = nullptr;
    TH1D* fake_unused_pt2_MD1_rz_simple = nullptr;

    // --- 2D Cut Study Histograms (X = dZ, Y = dXY) ---
    TH2D* h2_real_MD0_dXY_vs_dZ = nullptr;
    TH2D* h2_fake_MD0_dXY_vs_dZ = nullptr;
    TH2D* h2_real_MD1_dXY_vs_dZ = nullptr;
    TH2D* h2_fake_MD1_dXY_vs_dZ = nullptr;

    TH2D* h2_real_unused_MD0_dXY_vs_dZ = nullptr;
    TH2D* h2_fake_unused_MD0_dXY_vs_dZ = nullptr;
    TH2D* h2_real_unused_MD1_dXY_vs_dZ = nullptr;
    TH2D* h2_fake_unused_MD1_dXY_vs_dZ = nullptr;

    // Delta R Hists
    TH1D* real_pt2_deltaR = nullptr;
    TH1D* real_unused_pt2_deltaR = nullptr;
    TH1D* fake_pt2_deltaR = nullptr;
    TH1D* fake_unused_pt2_deltaR = nullptr;

    // Directional Angle Hists
    TH1D* real_pt2_deltaAngle = nullptr;
    TH1D* real_unused_pt2_deltaAngle = nullptr;
    TH1D* fake_pt2_deltaAngle = nullptr;
    TH1D* fake_unused_pt2_deltaAngle = nullptr;
};

#endif
