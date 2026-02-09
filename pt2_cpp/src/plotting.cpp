#include "plotting.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>
#include <filesystem>
#include <algorithm>

namespace fs = std::filesystem;

void Plotting::plotPt2Histograms(const HistogramManager& hists,
                                 const std::string& outputDir)
{
    // Ensure the output directory exists
    if (!fs::exists(outputDir)) {
        fs::create_directories(outputDir);
    }

    struct PlotInfo {
        TH1F* hist1;
        TH1F* hist2;
        std::string title;
        std::string filename;
    };

    // Define plots
    PlotInfo plots[] = {
        {hists.real_pt2_deltaPT,        hists.fake_pt2_deltaPT,
         "Delta p_{T} (Real vs Fake)",   "deltaPT_case1"},

        {hists.real_pt2_deltaETA,       hists.fake_pt2_deltaETA,
         "Delta #eta (Real vs Fake)",   "deltaETA_case1"},

        {hists.real_pt2_deltaPHI,       hists.fake_pt2_deltaPHI,
         "Delta #phi (Real vs Fake)",   "deltaPHI_case1"},

        {hists.real_unused_pt2_deltaPT,  hists.fake_unused_pt2_deltaPT,
         "Delta p_{T} Unused (Real vs Fake)", "deltaPT_case2"},

        {hists.real_unused_pt2_deltaETA, hists.fake_unused_pt2_deltaETA,
         "Delta #eta Unused (Real vs Fake)",  "deltaETA_case2"},

        {hists.real_unused_pt2_deltaPHI, hists.fake_unused_pt2_deltaPHI,
         "Delta #phi Unused (Real vs Fake)",  "deltaPHI_case2"},
        //----------------Extrapolation----------------------

        {hists.real_pt2_dist3D, hists.fake_pt2_dist3D,
         "Helical Extrapolation Distance (Real vs Fake)", "extrapDist_case1"},

        {hists.real_unused_pt2_dist3D, hists.fake_unused_pt2_dist3D,
         "Helical Extrapolation Distance Unused (Real vs Fake)", "extrapDist_case2"},
        
        //----------------------------------------------------
    };

    // Loop over all plots
    for (const auto& plot : plots) {

        if (!plot.hist1 || !plot.hist2) {
            std::cerr << "Warning: Missing histogram for "
                      << plot.title << std::endl;
            continue;
        }

        // ---------------------------------------------------------
        // Compute shared y-axis max (Real vs Fake)
        // ---------------------------------------------------------
        double ymax = std::max(plot.hist1->GetMaximum(),
                               plot.hist2->GetMaximum());

        // Safety fallback
        if (ymax <= 0.0) ymax = 1.0;

        // ---------------------------------------------------------
        // Make TWO versions: linear-Y and log-Y
        // ---------------------------------------------------------
        for (bool useLogY : {false, true}) {

            TCanvas canvas("canvas", plot.title.c_str(), 800, 600);
            canvas.cd();
            canvas.SetLogy(useLogY);

            // Scale y-axis appropriately
            double ymaxScaled = useLogY ? ymax * 10.0 : ymax * 1.2;

            plot.hist1->SetMaximum(ymaxScaled);
            plot.hist2->SetMaximum(ymaxScaled);

            // Log plots must not start at zero
            if (useLogY) {
                plot.hist1->SetMinimum(0.1);
                plot.hist2->SetMinimum(0.1);
            } else {
                plot.hist1->SetMinimum(0.0);
                plot.hist2->SetMinimum(0.0);
            }

            // Style
            plot.hist1->SetLineColor(kBlue);
            plot.hist1->SetLineWidth(2);
            plot.hist1->SetStats(0);
            plot.hist1->SetTitle(plot.title.c_str());

            plot.hist2->SetLineColor(kRed);
            plot.hist2->SetLineWidth(2);
            plot.hist2->SetStats(0);

            // Draw
            plot.hist1->Draw("HIST");
            plot.hist2->Draw("HIST SAME");

            // Legend
            TLegend legend(0.7, 0.75, 0.9, 0.9);
            legend.AddEntry(plot.hist1, "Real", "l");
            legend.AddEntry(plot.hist2, "Fake", "l");
            legend.Draw();

            // Output names
            std::string suffix = useLogY ? "_logY" : "_linY";
            std::string pdfPath = outputDir + "/" + plot.filename + suffix + ".pdf";
            std::string pngPath = outputDir + "/" + plot.filename + suffix + ".png";

            canvas.SaveAs(pdfPath.c_str());
            canvas.SaveAs(pngPath.c_str());
        }
    }

    std::cout << "All plots saved to " << outputDir << std::endl;
}

