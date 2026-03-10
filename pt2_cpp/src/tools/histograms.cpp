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

    //------------Extrapolation--------------------

    real_pt2_dist3D = new TH1D(
        "real_pt2_dist3D",
        "Real pT2 Distance;Dist [cm];Entries",
        100, 0, 15.0
    );

    fake_pt2_dist3D = new TH1D(
        "fake_pt2_dist3D",
        "Fake pT2 Distance;Dist [cm];Entries",
        100, 0, 15.0
    );

    real_unused_pt2_dist3D = new TH1D(
        "real_unused_pt2_dist3D",
        "Real Unused pT2 Distance;Dist [cm];Entries",
        100, 0, 15.0
    );

   fake_unused_pt2_dist3D = new TH1D(
       "fake_unused_pt2_dist3D",
       "Fake Unused pT2 Distance;Dist [cm];Entries",
       100, 0, 15.0
    );
   
    real_pt2_rz_simple = new TH1D(
        "real_pt2_rz_simple", 
        "Simple Pointing R-Z Residual;[cm];Entries", 
        200, -10, 10
    );
    
    fake_pt2_rz_simple = new TH1D(
        "fake_pt2_rz_simple", 
        "Simple Pointing R-Z Residual;[cm];Entries", 
        200, -10, 10
    );
     
    real_unused_pt2_rz_simple = new TH1D(
        "real_unused_pt2_rz_simple", 
        "Real Unused Simple R-Z;[cm];Entries", 
        200, -10, 10
    );
     
    fake_unused_pt2_rz_simple = new TH1D(
        "fake_unused_pt2_rz_simple", 
        "Fake Unused Simple R-Z;[cm];Entries",
        200, -10, 10
    );
   
    //Delta R

    real_pt2_deltaR = new TH1D(
        "real_pt2_deltaR", 
        "Real pT2 #Delta R; #Delta R; Entries", 
        100, 0, 1.0
    );
    
    fake_pt2_deltaR = new TH1D(
        "fake_pt2_deltaR", 
        "Fake pT2 #Delta R; #Delta R; Entries", 
        100, 0, 1.0
    );
    
    real_unused_pt2_deltaR = new TH1D(
        "real_unused_pt2_deltaR", 
        "Real Unused #Delta R; #Delta R; Entries", 
        100, 0, 1.0
    );
    
    fake_unused_pt2_deltaR = new TH1D(
        "fake_unused_pt2_deltaR", 
        "Fake Unused #Delta R; #Delta R; Entries", 
        100, 0, 1.0
    );
    
    //Delta angle

    real_pt2_deltaAngle = new TH1D(
        "real_pt2_deltaAngle", 
        "Real pT2 Directional #Delta#alpha; #Delta#alpha [rad]; Entries", 
        100, 0, 0.5
    );
    
    fake_pt2_deltaAngle = new TH1D(
        "fake_pt2_deltaAngle", 
        "Fake pT2 Directional #Delta#alpha; #Delta#alpha [rad]; Entries", 
        100, 0, 0.5
    );
    
    real_unused_pt2_deltaAngle = new TH1D(
        "real_unused_pt2_deltaAngle", 
        "Real Unused #Delta#alpha; [rad]; Entries", 
        100, 0, 0.5
    );
    
    fake_unused_pt2_deltaAngle = new TH1D(
        "fake_unused_pt2_deltaAngle", 
        "Fake Unused #Delta#alpha; [rad]; Entries", 
        100, 0, 0.5
    );
}
