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
        {
            .title = "All pT2: 3D Extrapolation Distance",
            .xAxis = "Distance [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_dist3D",
            .hists = { hists.real_pt2_dist3D, hists.fake_pt2_dist3D },
            .legend = {"Real", "Fake"},
          //  .normalize = true,
            .printYields = true
        },
        {
            .title = "Unused pT2: 3D Extrapolation Distance",
            .xAxis = "Distance [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_dist3D",
            .hists = { hists.real_unused_pt2_dist3D, hists.fake_unused_pt2_dist3D },
            .legend = {"Real", "Fake"},
          //  .normalize = true,
            .printYields = true
        },

/*        // --- NEW: R-Z Regression Residuals ---
        {
            .title = "All pT2: R-Z Regression Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_rz_reg",
            .hists = { hists.real_pt2_rz_reg, hists.fake_pt2_rz_reg },
            .legend = {"Real", "Fake"},
          //  .normalize = true,
            .printYields = true
        },
        {
            .title = "Unused pT2: R-Z Regression Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_rz_reg",
            .hists = { hists.real_unused_pt2_rz_reg, hists.fake_unused_pt2_rz_reg },
            .legend = {"Real", "Fake"},
          //  .normalize = true,
            .printYields = true
        },
*/
        // --- NEW: R-Z Simple Pointing Residuals ---
        {
            .title = "All pT2: R-Z Simple Pointing Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_rz_simple",
            .hists = { hists.real_pt2_rz_simple, hists.fake_pt2_rz_simple },
            .legend = {"Real", "Fake"},
          //  .normalize = true,
            .printYields = true
        },
        {
            .title = "Unused pT2: R-Z Simple Pointing Residual",
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_rz_simple",
            .hists = { hists.real_unused_pt2_rz_simple, hists.fake_unused_pt2_rz_simple },
            .legend = {"Real", "Fake"},
           // .normalize = true,
            .printYields = true
        },
        {
            .title = "pT2 Delta R", 
            .xAxis = "#Delta R", 
            .yAxis = "Normalized",
            .filename = "pt2_all_deltaR", 
            .hists = {hists.real_pt2_deltaR, hists.fake_pt2_deltaR},
            .legend = {"Real", "Fake"}, 
            //.normalize = true, 
            .printYields = true
        },
        
        {
            .title = "Unused pT2: Delta R",
            .xAxis = "#Delta R",
            .yAxis = "Normalized",
            .filename = "pt2_unused_deltaR",
            .hists = { hists.real_unused_pt2_deltaR, hists.fake_unused_pt2_deltaR },
            .legend = {"Real", "Fake"},
           // .normalize = true,
            .printYields = true
        },
        
        {
            .title = "pT2 Directional Delta Angle", 
            .xAxis = "#Delta#alpha [rad]", 
            .yAxis = "Normalized",
            .filename = "pt2_all_deltaAngle", 
            .hists = {hists.real_pt2_deltaAngle, hists.fake_pt2_deltaAngle},
            .legend = {"Real", "Fake"}, 
            //.normalize = true, 
            .printYields = true
        },
        
        {
            .title = "Unused pT2: Directional Delta Angle",
            .xAxis = "#Delta#alpha [rad]",
            .yAxis = "Normalized",
            .filename = "pt2_unused_deltaAngle",
            .hists = { hists.real_unused_pt2_deltaAngle, hists.fake_unused_pt2_deltaAngle },
            .legend = {"Real", "Fake"},
            //.normalize = true,
            .printYields = true
        },


     };
}

