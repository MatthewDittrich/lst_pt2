#include "plot_recipes.h"
#include "histograms.h"

#include <vector>

std::vector<PlotRecipe> getPt2Recipes(const HistogramManager& hists) {
    return {
        {
            .title = "All pT2: Delta p_{T}",
            .xAxis = "#Delta p_{T}",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPT",
            .hists = { hists.real_pt2_deltaPT, hists.fake_pt2_deltaPT },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "All pT2: Delta #eta",
            .xAxis = "#Delta #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaETA",
            .hists = { hists.real_pt2_deltaETA, hists.fake_pt2_deltaETA },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "All pT2: Delta #phi",
            .xAxis = "#Delta #phi",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPHI",
            .hists = { hists.real_pt2_deltaPHI, hists.fake_pt2_deltaPHI },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: Delta p_{T}",
            .xAxis = "#Delta p_{T}",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPT",
            .hists = { hists.real_unused_pt2_deltaPT, hists.fake_unused_pt2_deltaPT },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: Delta #eta",
            .xAxis = "#Delta #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaETA",
            .hists = { hists.real_unused_pt2_deltaETA, hists.fake_unused_pt2_deltaETA },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: Delta #phi",
            .xAxis = "#Delta #phi",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPHI",
            .hists = { hists.real_unused_pt2_deltaPHI, hists.fake_unused_pt2_deltaPHI },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // --- NEW: 3D Helical Distance ---

     /*   {
            .title = "All pT2: 3D Extrapolation Distance",
            .xAxis = "Distance [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_dist3D",
            .hists = { hists.real_pt2_dist3D, hists.fake_pt2_dist3D },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: 3D Extrapolation Distance",
            .xAxis = "Distance [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_dist3D",
            .hists = { hists.real_unused_pt2_dist3D, hists.fake_unused_pt2_dist3D },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
*/
        // =========================================================================
        // MD0 Helical Components (All)
        // =========================================================================
        {
            .title = "All pT2: MD0 Transverse Distance (dXY)",
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dxy",
            .hists = { hists.real_pt2_MD0_dXY, hists.fake_pt2_MD0_dXY },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "All pT2: MD0 Longitudinal Distance (dZ)",
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dz",
            .hists = { hists.real_pt2_MD0_dZ, hists.fake_pt2_MD0_dZ },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // =========================================================================
        // MD0 Helical Components (Unused)
        // =========================================================================
        {
            .title = "Unused pT2: MD0 Transverse Distance (dXY)",
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dxy",
            .hists = { hists.real_unused_pt2_MD0_dXY, hists.fake_unused_pt2_MD0_dXY },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: MD0 Longitudinal Distance (dZ)",
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dz",
            .hists = { hists.real_unused_pt2_MD0_dZ, hists.fake_unused_pt2_MD0_dZ },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // =========================================================================
        // MD1 Helical Components (All)
        // =========================================================================
        {
            .title = "All pT2: MD1 Transverse Distance (dXY)",
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dxy",
            .hists = { hists.real_pt2_MD1_dXY, hists.fake_pt2_MD1_dXY },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "All pT2: MD1 Longitudinal Distance (dZ)",
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dz",
            .hists = { hists.real_pt2_MD1_dZ, hists.fake_pt2_MD1_dZ },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // =========================================================================
        // MD1 Helical Components (Unused)
        // =========================================================================
        {
            .title = "Unused pT2: MD1 Transverse Distance (dXY)",
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dxy",
            .hists = { hists.real_unused_pt2_MD1_dXY, hists.fake_unused_pt2_MD1_dXY },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: MD1 Longitudinal Distance (dZ)",
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dz",
            .hists = { hists.real_unused_pt2_MD1_dZ, hists.fake_unused_pt2_MD1_dZ },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        // =========================================================================
        // R-Z Simple Pointing (Separated by MD0 and MD1)
        // =========================================================================

        // MD0 - All
        {
            .title = "All pT2: MD0 R-Z Simple Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_rz_simple",
            .hists = { hists.real_pt2_MD0_rz_simple, hists.fake_pt2_MD0_rz_simple },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // MD1 - All
        {
            .title = "All pT2: MD1 R-Z Simple Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_rz_simple",
            .hists = { hists.real_pt2_MD1_rz_simple, hists.fake_pt2_MD1_rz_simple },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // MD0 - Unused Only
        {
            .title = "Unused pT2: MD0 R-Z Simple Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_rz_simple",
            .hists = { hists.real_unused_pt2_MD0_rz_simple, hists.fake_unused_pt2_MD0_rz_simple },
            .legend = {"Real", "Fake"},
            .printYields = true
        },

        // MD1 - Unused Only
        {
            .title = "Unused pT2: MD1 R-Z Simple Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_rz_simple",
            .hists = { hists.real_unused_pt2_MD1_rz_simple, hists.fake_unused_pt2_MD1_rz_simple },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
              
        {          
            .title = "pT2 Delta R", 
            .xAxis = "#Delta R", 
            .yAxis = "Entries",
            .filename = "pt2_all_deltaR", 
            .hists = {hists.real_pt2_deltaR, hists.fake_pt2_deltaR},
            .legend = {"Real", "Fake"}, 
            .printYields = true
        },
        {
            .title = "Unused pT2: Delta R",
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaR",
            .hists = { hists.real_unused_pt2_deltaR, hists.fake_unused_pt2_deltaR },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "pT2 Directional Delta Angle", 
            .xAxis = "#Delta#alpha [rad]", 
            .yAxis = "Entries",
            .filename = "pt2_all_deltaAngle", 
            .hists = {hists.real_pt2_deltaAngle, hists.fake_pt2_deltaAngle},
            .legend = {"Real", "Fake"}, 
            .printYields = true
        },
        {
            .title = "Unused pT2: Directional Delta Angle",
            .xAxis = "#Delta#alpha [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaAngle",
            .hists = { hists.real_unused_pt2_deltaAngle, hists.fake_unused_pt2_deltaAngle },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        // --- NEW: LST dPhi ---
        {
            .title = "All pT2: LST #Delta#phi", 
            .xAxis = "LST #Delta#phi [rad]", 
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdPhi", 
            .hists = {hists.real_pt2_LSTdPhi, hists.fake_pt2_LSTdPhi},
            .legend = {"Real", "Fake"}, 
            .printYields = true
        },
        {
            .title = "Unused pT2: LST #Delta#phi",
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdPhi",
            .hists = { hists.real_unused_pt2_LSTdPhi, hists.fake_unused_pt2_LSTdPhi },
            .legend = {"Real", "Fake"},
            .printYields = true
        }, 
        // --- LST dBeta ---
        {
            .title = "All pT2: LST #Delta#beta",
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdBeta",
            .hists = {hists.real_pt2_LSTdBeta, hists.fake_pt2_LSTdBeta},
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: LST #Delta#beta",
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdBeta",
            .hists = { hists.real_unused_pt2_LSTdBeta, hists.fake_unused_pt2_LSTdBeta },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        // --- LST betaOut ---
        {
            .title = "All pT2: LST #beta_{Out}",
            .xAxis = "LST #beta_{Out} [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTbetaOut",
            .hists = {hists.real_pt2_LSTbetaOut, hists.fake_pt2_LSTbetaOut},
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: LST #beta_{Out}",
            .xAxis = "LST #beta_{Out} [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTbetaOut",
            .hists = { hists.real_unused_pt2_LSTbetaOut, hists.fake_unused_pt2_LSTbetaOut },
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        // --- LST Z-Residuals ---
        // --- NEW: LST Z-Residuals ---
        {
            .title = "All pT2: LST Geometric Z-Residual (Beamspot Cut)", 
            .xAxis = "Actual Z - Origin Predicted Z [cm]", 
            .yAxis = "Entries",
            .filename = "pt2_all_LSTOrgZRes", 
            .hists = {hists.real_pt2_LSTOrgZRes, hists.fake_pt2_LSTOrgZRes},
            .legend = {"Real", "Fake"}, 
            .printYields = true
        },
        {
            .title = "All pT2: LST Kinematic Z-Residual (Momentum Cut)", 
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]", 
            .yAxis = "Entries",
            .filename = "pt2_all_LSTKinZRes", 
            .hists = {hists.real_pt2_LSTKinZRes, hists.fake_pt2_LSTKinZRes},
            .legend = {"Real", "Fake"}, 
            .printYields = true
        },
        {
            .title = "Unused pT2: LST Geometric Z-Residual (Beamspot Cut)",
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTOrgZRes",
            .hists = {hists.real_unused_pt2_LSTOrgZRes, hists.fake_unused_pt2_LSTOrgZRes},
            .legend = {"Real", "Fake"},
            .printYields = true
        },
        {
            .title = "Unused pT2: LST Kinematic Z-Residual (Momentum Cut)",
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTKinZRes",
            .hists = {hists.real_unused_pt2_LSTKinZRes, hists.fake_unused_pt2_LSTKinZRes},
            .legend = {"Real", "Fake"},
            .printYields = true
        }

    };
}

