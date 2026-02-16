#include "histograms.h"
#include <cmath>

void HistogramManager::init() {

    real_pt2_deltaPT = new TH1D(
        "real_pt2_deltaPT",
        "Real pT2 #Delta pT;#Delta pT;Entries",
        180, -10.0, 10.0
    );

    real_pt2_deltaETA = new TH1D(
        "real_pt2_deltaETA",
        "Real pT2 #Delta #eta;#Delta #eta;Entries",
        180, -1.0, 1.0
    );

    real_pt2_deltaPHI = new TH1D(
        "real_pt2_deltaPHI",
        "Real pT2 #Delta #phi;#Delta #phi;Entries",
        180, -1.0, 1.0
    );

    real_unused_pt2_deltaPT = new TH1D(
        "real_unused_pt2_deltaPT",
        "Real Unused pT2 #Delta pT;#Delta pT;Entries",
        180, -10.0, 10.0
    );

    real_unused_pt2_deltaETA = new TH1D(
        "real_unused_pt2_deltaETA",
        "Real Unused pT2 #Delta #eta;#Delta #eta;Entries",
        180, -1.0, 1.0
    );

    real_unused_pt2_deltaPHI = new TH1D(
        "real_unused_pt2_deltaPHI",
        "Real Unused pT2 #Delta #phi;#Delta #phi;Entries",
        180, -1.0, 1.0
    );

    fake_pt2_deltaPT = new TH1D(
        "fake_pt2_deltaPT",
        "Fake pT2 #Delta pT;#Delta pT;Entries",
        180, -10.0, 10.0
    );

    fake_pt2_deltaETA = new TH1D(
        "fake_pt2_deltaETA",
        "Fake pT2 #Delta #eta;#Delta #eta;Entries",
        180, -1.0, 1.0
    );

    fake_pt2_deltaPHI = new TH1D(
        "fake_pt2_deltaPHI",
        "Fake pT2 #Delta #phi;#Delta #phi;Entries",
        180, -1.0, 1.0
    );

    fake_unused_pt2_deltaPT = new TH1D(
        "fake_unused_pt2_deltaPT",
        "Fake Unused pT2 #Delta pT;#Delta pT;Entries",
        180, -10.0, 10.0
    );

    fake_unused_pt2_deltaETA = new TH1D(
        "fake_unused_pt2_deltaETA",
        "Fake Unused pT2 #Delta #eta;#Delta #eta;Entries",
        180, -1.0, 1.0
    );

    fake_unused_pt2_deltaPHI = new TH1D(
        "fake_unused_pt2_deltaPHI",
        "Fake Unused pT2 #Delta #phi;#Delta #phi;Entries",
        180, -1.0, 1.0
    );

}
