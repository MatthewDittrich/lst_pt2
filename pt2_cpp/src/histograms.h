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
};

#endif

