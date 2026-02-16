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
};

#endif

