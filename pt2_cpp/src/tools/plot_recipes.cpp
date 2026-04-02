#include "plot_recipes.h"
#include "histograms.h"
#include <vector>
#include <string>

std::vector<PlotRecipe> getPt2Recipes(const HistogramManager& hists) {
    std::vector<PlotRecipe> recipes;

    for (int i = 0; i < 8; ++i) {
        std::string sfx = "_" + hists.catNames[i];
        std::string ttl = " (" + hists.catTitles[i] + ")";

        // =====================================================================
        // ALL pT2s - KINEMATICS & LST VARIABLES
        // =====================================================================

        recipes.push_back({
            .title = "pT2 Delta p_{T}" + ttl,
            .xAxis = "#Delta p_{T} [GeV]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPT" + sfx,
            .hists = {hists.real_pt2_deltaPT[i], hists.fake_pt2_deltaPT[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPHI" + sfx,
            .hists = {hists.real_pt2_deltaPHI[i], hists.fake_pt2_deltaPHI[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaR" + sfx,
            .hists = {hists.real_pt2_deltaR[i], hists.fake_pt2_deltaR[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdBeta" + sfx,
            .hists = {hists.real_pt2_LSTdBeta[i], hists.fake_pt2_LSTdBeta[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTKinZRes" + sfx,
            .hists = {hists.real_pt2_LSTKinZRes[i], hists.fake_pt2_LSTKinZRes[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTOrgZRes" + sfx,
            .hists = {hists.real_pt2_LSTOrgZRes[i], hists.fake_pt2_LSTOrgZRes[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // =====================================================================
        // ALL pT2s - MD COMPONENTS
        // =====================================================================

        recipes.push_back({
            .title = "MD0 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dxy" + sfx,
            .hists = {hists.real_pt2_MD0_dXY[i], hists.fake_pt2_MD0_dXY[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dz" + sfx,
            .hists = {hists.real_pt2_MD0_dZ[i], hists.fake_pt2_MD0_dZ[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dxy" + sfx,
            .hists = {hists.real_pt2_MD1_dXY[i], hists.fake_pt2_MD1_dXY[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dz" + sfx,
            .hists = {hists.real_pt2_MD1_dZ[i], hists.fake_pt2_MD1_dZ[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_rz_simple" + sfx,
            .hists = {hists.real_pt2_MD0_rz_simple[i], hists.fake_pt2_MD0_rz_simple[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_rz_simple" + sfx,
            .hists = {hists.real_pt2_MD1_rz_simple[i], hists.fake_pt2_MD1_rz_simple[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // =====================================================================
        // UNUSED pT2s - KINEMATICS & LST VARIABLES
        // =====================================================================

        recipes.push_back({
            .title = "Unused pT2 Delta p_{T}" + ttl,
            .xAxis = "#Delta p_{T} [GeV]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPT" + sfx,
            .hists = {hists.real_unused_pt2_deltaPT[i], hists.fake_unused_pt2_deltaPT[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPHI" + sfx,
            .hists = {hists.real_unused_pt2_deltaPHI[i], hists.fake_unused_pt2_deltaPHI[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdBeta" + sfx,
            .hists = {hists.real_unused_pt2_LSTdBeta[i], hists.fake_unused_pt2_LSTdBeta[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTKinZRes" + sfx,
            .hists = {hists.real_unused_pt2_LSTKinZRes[i], hists.fake_unused_pt2_LSTKinZRes[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTOrgZRes" + sfx,
            .hists = {hists.real_unused_pt2_LSTOrgZRes[i], hists.fake_unused_pt2_LSTOrgZRes[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // =====================================================================
        // UNUSED pT2s - MD COMPONENTS
        // =====================================================================

        recipes.push_back({
            .title = "Unused MD0 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dxy" + sfx,
            .hists = {hists.real_unused_pt2_MD0_dXY[i], hists.fake_unused_pt2_MD0_dXY[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dz" + sfx,
            .hists = {hists.real_unused_pt2_MD0_dZ[i], hists.fake_unused_pt2_MD0_dZ[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dxy" + sfx,
            .hists = {hists.real_unused_pt2_MD1_dXY[i], hists.fake_unused_pt2_MD1_dXY[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dz" + sfx,
            .hists = {hists.real_unused_pt2_MD1_dZ[i], hists.fake_unused_pt2_MD1_dZ[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_rz_simple" + sfx,
            .hists = {hists.real_unused_pt2_MD0_rz_simple[i], hists.fake_unused_pt2_MD0_rz_simple[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_rz_simple" + sfx,
            .hists = {hists.real_unused_pt2_MD1_rz_simple[i], hists.fake_unused_pt2_MD1_rz_simple[i]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });
    }

    return recipes;
}
