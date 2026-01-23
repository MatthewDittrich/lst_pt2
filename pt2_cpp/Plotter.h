#ifndef PLOTTER_H
#define PLOTTER_H

#include <string>
#include <vector>

// ROOT Includes
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <Rtypes.h> // For Color_t

namespace Plotter {

    /**
     * Creates and saves a stacked histogram comparing Real vs Fake objects.
     * Used for: Simple Used/Unused comparisons.
     */
    void createAndSaveSimpleStack(
        const std::string& canvasName, 
        const std::string& stackTitle, 
        const std::string& saveFileName, 
        const std::string& legendHeader, 
        TH1D* h_real, 
        TH1D* h_fake, 
        const std::string& real_label, 
        const std::string& fake_label, 
        Color_t real_color, 
        Color_t fake_color, 
        bool setLogY = false
    );

    /**
     * Creates and saves a detailed stacked histogram showing matching categories.
     * Used for: Detailed Used/Unused comparisons (Matched to Both, Unused Only, etc).
     */
    void createAndSaveDetailedStack(
        const std::string& canvasName, 
        const std::string& stackTitle, 
        const std::string& saveFileName, 
        const std::string& legendHeader, 
        bool setLogY, 
        TH1D* h_fake, 
        TH1D* h_real_unmatched, 
        TH1D* h_real_matchedToUnused, 
        TH1D* h_real_matchedToBoth, 
        TH1D* h_real_matchedToUsed
    );

    /**
     * Creates and saves a standard 1D histogram with a statistics box for entries.
     * Used for: Delta plots, Extrapolation distances, Layer plots.
     */
    void createAndSaveSimplePlot(
        const std::string& canvasName, 
        const std::string& saveFileName, 
        TH1D* histo, 
        const std::string& plotTitle, 
        const std::string& xLabel, 
        const std::string& yLabel
    );

}

#endif // PLOTTER_H
