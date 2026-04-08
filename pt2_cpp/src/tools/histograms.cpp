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

   /* real_pt2_dist3D = new TH1D(
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
   */
    // =========================================================================
    // MD0 Helical Components
    // =========================================================================

    real_pt2_MD0_dXY = new TH1D(
        "real_pt2_MD0_dXY",
        "Real MD0 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    real_pt2_MD0_dZ = new TH1D(
        "real_pt2_MD0_dZ",
        "Real MD0 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    fake_pt2_MD0_dXY = new TH1D(
        "fake_pt2_MD0_dXY",
        "Fake MD0 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    fake_pt2_MD0_dZ = new TH1D(
        "fake_pt2_MD0_dZ",
        "Fake MD0 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    real_unused_pt2_MD0_dXY = new TH1D(
        "real_unused_pt2_MD0_dXY",
        "Real Unused MD0 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    real_unused_pt2_MD0_dZ = new TH1D(
        "real_unused_pt2_MD0_dZ",
        "Real Unused MD0 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    fake_unused_pt2_MD0_dXY = new TH1D(
        "fake_unused_pt2_MD0_dXY",
        "Fake Unused MD0 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    fake_unused_pt2_MD0_dZ = new TH1D(
        "fake_unused_pt2_MD0_dZ",
        "Fake Unused MD0 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    // =========================================================================
    // MD1 Helical Components
    // =========================================================================

    real_pt2_MD1_dXY = new TH1D(
        "real_pt2_MD1_dXY",
        "Real MD1 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    real_pt2_MD1_dZ = new TH1D(
        "real_pt2_MD1_dZ",
        "Real MD1 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    fake_pt2_MD1_dXY = new TH1D(
        "fake_pt2_MD1_dXY",
        "Fake MD1 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    fake_pt2_MD1_dZ = new TH1D(
        "fake_pt2_MD1_dZ",
        "Fake MD1 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    real_unused_pt2_MD1_dXY = new TH1D(
        "real_unused_pt2_MD1_dXY",
        "Real Unused MD1 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    real_unused_pt2_MD1_dZ = new TH1D(
        "real_unused_pt2_MD1_dZ",
        "Real Unused MD1 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    fake_unused_pt2_MD1_dXY = new TH1D(
        "fake_unused_pt2_MD1_dXY",
        "Fake Unused MD1 #DeltaXY;#DeltaXY [cm];Entries",
        180, 0, 5.0
    );

    fake_unused_pt2_MD1_dZ = new TH1D(
        "fake_unused_pt2_MD1_dZ",
        "Fake Unused MD1 #DeltaZ;#DeltaZ [cm];Entries",
        180, 0, 10.0
    );

    real_pt2_MD0_rz_simple = new TH1D(
        "real_pt2_MD0_rz_simple", 
        "Real MD0 R-Z Residual;[cm];Entries", 
        180, -10, 10
    );
    
    real_pt2_MD1_rz_simple = new TH1D(
        "real_pt2_MD1_rz_simple", 
        "Real MD1 R-Z Residual;[cm];Entries", 
        180, -10, 10
    );
    
    fake_pt2_MD0_rz_simple = new TH1D(
        "fake_pt2_MD0_rz_simple", 
        "Fake MD0 R-Z Residual;[cm];Entries", 
        180, -10, 10
    );
    
    fake_pt2_MD1_rz_simple = new TH1D(
        "fake_pt2_MD1_rz_simple", 
        "Fake MD1 R-Z Residual;[cm];Entries", 
        180, -10, 10
    );

    real_unused_pt2_MD0_rz_simple = new TH1D(
        "real_unused_pt2_MD0_rz_simple",
        "Real Unused MD0 R-Z Residual;[cm];Entries",
        180, -10, 10
    );

    real_unused_pt2_MD1_rz_simple = new TH1D(
        "real_unused_pt2_MD1_rz_simple",
        "Real Unused MD1 R-Z Residual;[cm];Entries",
        180, -10, 10
    );

    fake_unused_pt2_MD0_rz_simple = new TH1D(
        "fake_unused_pt2_MD0_rz_simple",
        "Fake Unused MD0 R-Z Residual;[cm];Entries",
        180, -10, 10
    );

    fake_unused_pt2_MD1_rz_simple = new TH1D(
        "fake_unused_pt2_MD1_rz_simple",
        "Fake Unused MD1 R-Z Residual;[cm];Entries",
        180, -10, 10
    );

    // Cut Study: dXY distribution for events that pass |dZ| < 2.5 cm
    // Add this inside HistogramManager::init()
    
    h2_real_MD0_dXY_vs_dZ = new TH2D(
        "h2_real_MD0_dXY_vs_dZ", 
        "Real MD0; |#DeltaZ| [cm]; dXY [cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );
    
    h2_fake_MD0_dXY_vs_dZ = new TH2D(
        "h2_fake_MD0_dXY_vs_dZ", 
        "Fake MD0; |#DeltaZ| [cm]; dXY [cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );
    
    h2_real_MD1_dXY_vs_dZ = new TH2D(
        "h2_real_MD1_dXY_vs_dZ", 
        "Real MD1; |#DeltaZ| [cm]; dXY [cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );
    
    h2_fake_MD1_dXY_vs_dZ = new TH2D(
        "h2_fake_MD1_dXY_vs_dZ", 
        "Fake MD1; |#DeltaZ|[cm]; dXY [cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );

    h2_real_unused_MD0_dXY_vs_dZ = new TH2D(
        "h2_real_unused_MD0_dXY_vs_dZ", 
        "Real Unused MD0; |#DeltaZ| [cm]; dXY[cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );
    
    h2_fake_unused_MD0_dXY_vs_dZ = new TH2D(
        "h2_fake_unused_MD0_dXY_vs_dZ", 
        "Fake Unused MD0; |#DeltaZ| [cm]; dXY [cm]",
        180, 0, 10.0, 180, 0, 10.0
    );
    
    h2_real_unused_MD1_dXY_vs_dZ = new TH2D(
        "h2_real_unused_MD1_dXY_vs_dZ", 
        "Real Unused MD1; |#DeltaZ| [cm]; dXY [cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );
    
    h2_fake_unused_MD1_dXY_vs_dZ = new TH2D(
        "h2_fake_unused_MD1_dXY_vs_dZ", 
        "Fake Unused MD1; |#DeltaZ| [cm]; dXY [cm]", 
        180, 0, 10.0, 180, 0, 10.0
    );

    //Delta R

    real_pt2_deltaR = new TH1D(
        "real_pt2_deltaR", 
        "Real pT2 #Delta R; #Delta R; Entries", 
        180, 0, 1.0
    );
    
    fake_pt2_deltaR = new TH1D(
        "fake_pt2_deltaR", 
        "Fake pT2 #Delta R; #Delta R; Entries", 
        180, 0, 1.0
    );
    
    real_unused_pt2_deltaR = new TH1D(
        "real_unused_pt2_deltaR", 
        "Real Unused #Delta R; #Delta R; Entries", 
        180, 0, 1.0
    );
    
    fake_unused_pt2_deltaR = new TH1D(
        "fake_unused_pt2_deltaR", 
        "Fake Unused #Delta R; #Delta R; Entries", 
        180, 0, 1.0
    );
    
    //Delta angle

    real_pt2_deltaAngle = new TH1D(
        "real_pt2_deltaAngle", 
        "Real pT2 Directional #Delta#alpha; #Delta#alpha [rad]; Entries", 
        180, 0, 0.5
    );
    
    fake_pt2_deltaAngle = new TH1D(
        "fake_pt2_deltaAngle", 
        "Fake pT2 Directional #Delta#alpha; #Delta#alpha [rad]; Entries", 
        180, 0, 0.5
    );
    
    real_unused_pt2_deltaAngle = new TH1D(
        "real_unused_pt2_deltaAngle", 
        "Real Unused #Delta#alpha; [rad]; Entries", 
        180, 0, 0.5
    );
    
    fake_unused_pt2_deltaAngle = new TH1D(
        "fake_unused_pt2_deltaAngle", 
        "Fake Unused #Delta#alpha; [rad]; Entries", 
        180, 0, 0.5
    );
}
