#include "histograms.h"
#include <cmath>

void HistogramManager::init() {

    real_pt2_deltaPT = new TH1F(
        "real_pt2_deltaPT",
        "Real pT2 #Delta pT;#Delta pT;Entries",
        180, -0.6, 0.6
    );

    real_pt2_deltaETA = new TH1F(
        "real_pt2_deltaETA",
        "Real pT2 #Delta #eta;#Delta #eta;Entries",
        180, -0.1, 0.1
    );

    real_pt2_deltaPHI = new TH1F(
        "real_pt2_deltaPHI",
        "Real pT2 #Delta #phi;#Delta #phi;Entries",
        180, -0.25, 0.25
    );

    real_unused_pt2_deltaPT = new TH1F(
        "real_unused_pt2_deltaPT",
        "Real Unused pT2 #Delta pT;#Delta pT;Entries",
        180, -0.6, 0.6
    );

    real_unused_pt2_deltaETA = new TH1F(
        "real_unused_pt2_deltaETA",
        "Real Unused pT2 #Delta #eta;#Delta #eta;Entries",
        180, -0.1, 0.1
    );

    real_unused_pt2_deltaPHI = new TH1F(
        "real_unused_pt2_deltaPHI",
        "Real Unused pT2 #Delta #phi;#Delta #phi;Entries",
        180, -0.25, 0.25
    );

    fake_pt2_deltaPT = new TH1F(
        "fake_pt2_deltaPT",
        "Fake pT2 #Delta pT;#Delta pT;Entries",
        180, -0.6, 0.6
    );

    fake_pt2_deltaETA = new TH1F(
        "fake_pt2_deltaETA",
        "Fake pT2 #Delta #eta;#Delta #eta;Entries",
        180, -0.1, 0.1
    );

    fake_pt2_deltaPHI = new TH1F(
        "fake_pt2_deltaPHI",
        "Fake pT2 #Delta #phi;#Delta #phi;Entries",
        180, -0.25, 0.25
    );

    fake_unused_pt2_deltaPT = new TH1F(
        "fake_unused_pt2_deltaPT",
        "Fake Unused pT2 #Delta pT;#Delta pT;Entries",
        180, -0.6, 0.6
    );

    fake_unused_pt2_deltaETA = new TH1F(
        "fake_unused_pt2_deltaETA",
        "Fake Unused pT2 #Delta #eta;#Delta #eta;Entries",
        180, -0.1, 0.1
    );

    fake_unused_pt2_deltaPHI = new TH1F(
        "fake_unused_pt2_deltaPHI",
        "Fake Unused pT2 #Delta #phi;#Delta #phi;Entries",
        180, -0.25, 0.25
    );
    //------------Extrapolation--------------------
    real_pt2_dist3D = new TH1F(
        "real_pt2_dist3D", 
         "Real pT2 Distance;Dist [cm];Entries", 
         100, 0, 15.0
    );
    
    fake_pt2_dist3D = new TH1F(
        "fake_pt2_dist3D", 
        "Fake pT2 Distance;Dist [cm];Entries", 
        100, 0, 15.0
        );

    real_unused_pt2_dist3D = new TH1F(
        "real_unused_pt2_dist3D", 
        "Real Unused pT2 Distance;Dist [cm];Entries", 
        100, 0, 15.0
        );

    fake_unused_pt2_dist3D = new TH1F(
        "fake_unused_pt2_dist3D", 
        "Fake Unused pT2 Distance;Dist [cm];Entries", 
        100, 0, 15.0
        );

}
