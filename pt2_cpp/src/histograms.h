#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include <TH1F.h>

class HistogramManager {
public:
    HistogramManager() = default;
    ~HistogramManager() = default;

    void init();
    void write();

    // Delta pT2 Hists
    TH1F* real_pt2_deltaPT  = nullptr;
    TH1F* real_pt2_deltaETA = nullptr;
    TH1F* real_pt2_deltaPHI = nullptr;

    TH1F* real_unused_pt2_deltaPT  = nullptr;
    TH1F* real_unused_pt2_deltaETA = nullptr;
    TH1F* real_unused_pt2_deltaPHI = nullptr;

    TH1F* fake_pt2_deltaPT  = nullptr;
    TH1F* fake_pt2_deltaETA = nullptr;
    TH1F* fake_pt2_deltaPHI = nullptr;

    TH1F* fake_unused_pt2_deltaPT  = nullptr;
    TH1F* fake_unused_pt2_deltaETA = nullptr;
    TH1F* fake_unused_pt2_deltaPHI = nullptr;
    //------------extrapolation--------------
    TH1F* real_pt2_dist3D;
    TH1F* fake_pt2_dist3D;
    TH1F* real_unused_pt2_dist3D;
    TH1F* fake_unused_pt2_dist3D;

<<<<<<< Updated upstream
=======
    TH1D* fake_unused_pt2_deltaPT  = nullptr;
    TH1D* fake_unused_pt2_deltaETA = nullptr;
    TH1D* fake_unused_pt2_deltaPHI = nullptr;

    //------------extrapolation--------------
    TH1D* real_pt2_dist3D = nullptr;
    TH1D* fake_pt2_dist3D = nullptr;
    TH1D* real_unused_pt2_dist3D = nullptr;
    TH1D* fake_unused_pt2_dist3D = nullptr;

    TH1D* real_pt2_rz_reg = nullptr;
    TH1D* fake_pt2_rz_reg = nullptr;
    TH1D* real_unused_pt2_rz_reg = nullptr;
    TH1D* fake_unused_pt2_rz_reg = nullptr;

    // 2. Simple Pointing Histograms
    TH1D* real_pt2_rz_simple = nullptr;
    TH1D* fake_pt2_rz_simple = nullptr;
    TH1D* real_unused_pt2_rz_simple = nullptr;
    TH1D* fake_unused_pt2_rz_simple = nullptr;

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
>>>>>>> Stashed changes
};

#endif

