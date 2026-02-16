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
        }
    };
}
