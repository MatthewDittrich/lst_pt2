#ifndef PLOTTING_H
#define PLOTTING_H

#include <string>
#include "histograms.h"

class Plotting {
public:
    Plotting() = default;
    ~Plotting() = default;

    // Plot all pt2 histograms and save to outputDir
    static void plotPt2Histograms(const HistogramManager& hists, const std::string& outputDir);
};

#endif

