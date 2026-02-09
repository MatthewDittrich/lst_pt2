#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <getopt.h>

#include "gator.h"
#include "histograms.h"
#include "plotting.h"
#include "rootReader.h"
#include "tools.h"
#include "pt2.h"
#include "extrapolation.h" // Added header for extrapolation

int main(int argc, char** argv) {

    // 1. Default Arguments
    bool makePlots = false;
    bool writeRoot = false;
    bool lowPT = false;
    std::string inputFile = "/blue/avery/matthew.dittrich/Permanent/LSTNtuple.root";
    std::string outputDir = "output";

    // 2. Command Line Arguments
    int opt;
    while ((opt = getopt(argc, argv, "pi:o:r:k")) != -1) {
        switch (opt) {
            case 'p': makePlots = true; break;
            case 'r': writeRoot = true; break;
            case 'k': lowPT = true;     break;
            case 'i': inputFile = optarg; break;
            case 'o': outputDir = optarg; break;
            default:
                std::cerr << "Usage: " << argv[0] << "\n"
                    << "[-p] Make Plots\n"
                    << "[-k] Run Low pT\n"
                    << "[-r] Make Root File \n"
                    << "[-i] Input File Path \n"
                    << "[-o] Output Directory \n";
                return 1;
        }    
    }

    // 3. Setup and Configuration
    print_gator();
    std::cout << "\n=== Configuration ===\n";
    std::cout << "Input file:        " << inputFile << "\n";
    std::cout << "Output directory:  " << outputDir << "\n";
    std::cout << "Make plots:        " << (makePlots ? "yes" : "no") << "\n";
    std::cout << "=====================\n\n";

    std::filesystem::create_directories(outputDir);

    HistogramManager hists;
    hists.init();

    if (lowPT || writeRoot) {
        throw std::runtime_error("Error: options --lowPT and --writeRoot are not currently set up.");
    }

    // 4. Load Pixel Map
    SuperbinToDetIdMap superbinToDetIds;
    std::string pixelMapFile = "pLS_map_ElCheapo.txt";
    try {
        loadSuperbinDetIdMap(pixelMapFile, superbinToDetIds);
    } catch (const std::exception& e) {
        std::cerr << "Error loading superbin map: " << e.what() << std::endl;
        return 1;
    }

    // 5. Setup Root Reader
    rootReader reader;
    if (!reader.Init(inputFile, "tree")) {
        std::cerr << "Error: Could not find Root File or TTree! \n" << std::endl;
        return 1;
    } 

    DetIdToLSMap detidToLS;
    pT2Collection pt2s;

    // 6. Main Event Looper
    for (Long64_t ievt = 0; ievt < reader.GetEntries(); ++ievt) {
        std::cout << "Processing Event: " << ievt << " / " << reader.GetEntries() << std::endl;
        reader.GetEntry(ievt);

        // Reset per-event variables
        reader.pls_origin_z.clear();
        reader.pls_superbin.clear();
        detidToLS.clear();
        pt2s.clear();
        reader.pls_isUsed.clear();
        reader.ls_isUsed.clear();

        UsedMask usedMask = buildUsedMask(reader);
        reader.ls_isUsed  = std::move(usedMask.ls_isUsed);
        reader.pls_isUsed = std::move(usedMask.pls_isUsed);

        size_t nLS = reader.ls_pt->size();
        size_t nPLS = reader.pls_pt->size();

        reader.pls_origin_z.reserve(nPLS);
        reader.pls_superbin.reserve(nPLS);

        // Map LS to DetIds
        for (size_t k = 0; k < nLS; ++k) {
            std::vector<int> detIds = getDetIdsForLS(reader, k);
            for (int detId : detIds) {
                detidToLS[detId].push_back(k);
            }
        }

        // Build pT2 Candidates
        for (size_t j = 0; j < nPLS; ++j) {
            reader.pls_origin_z.push_back(CalculatePlsZ(reader, j));
            reader.pls_superbin.push_back(CalculateSuperbin(reader, j));
            buildPt2sForPLS(j, reader, superbinToDetIds, detidToLS, pt2s);
        }

        // --- pT2 Object Analysis Loop ---
        for (auto& pt2 : pt2s) {
            size_t plsIdx = pt2.pls_idx;
            size_t lsIdx = pt2.ls_idx;
            
            // Calculate Physics observables
            pt2.delta_pt  = deltaPt(reader.pls_pt->at(plsIdx), reader.ls_pt->at(lsIdx));
            pt2.delta_eta = deltaEta(reader.pls_eta->at(plsIdx), reader.ls_eta->at(lsIdx));
            pt2.delta_phi = deltaPhi(reader.pls_phi->at(plsIdx), reader.ls_phi->at(lsIdx));
            pt2.is_real   = pt2TruthFinder(reader, plsIdx, lsIdx);
            pt2.is_used   = pt2UsedCalculator(reader, plsIdx, lsIdx);

            // Calculate Helical Extrapolation Distances
            std::pair<double, double> dists = extrapolation::extrapolatePlsHelicallyAndGetDistance(plsIdx, lsIdx, reader);

            // Fill Histograms: Case 1 (All Real vs Fake)
            if (pt2.is_real) {
                hists.real_pt2_deltaPT->Fill(pt2.delta_pt);
                hists.real_pt2_deltaETA->Fill(pt2.delta_eta);
                hists.real_pt2_deltaPHI->Fill(pt2.delta_phi);
                if (dists.first >= 0) {
                    hists.real_pt2_dist3D->Fill(dists.first);
                    hists.real_pt2_dist3D->Fill(dists.second);
                }

                // Fill Histograms: Case 2 (Unused Only)
                if (!pt2.is_used) {
                    hists.real_unused_pt2_deltaPT->Fill(pt2.delta_pt);
                    hists.real_unused_pt2_deltaETA->Fill(pt2.delta_eta);
                    hists.real_unused_pt2_deltaPHI->Fill(pt2.delta_phi);
                    if (dists.first >= 0) {
                        hists.real_unused_pt2_dist3D->Fill(dists.first);
                        hists.real_unused_pt2_dist3D->Fill(dists.second);
                    }
                }
            } 
            else {
                // Fake - Case 1 (All)
                hists.fake_pt2_deltaPT->Fill(pt2.delta_pt);
                hists.fake_pt2_deltaETA->Fill(pt2.delta_eta);
                hists.fake_pt2_deltaPHI->Fill(pt2.delta_phi);
                if (dists.first >= 0) {
                    hists.fake_pt2_dist3D->Fill(dists.first);
                    hists.fake_pt2_dist3D->Fill(dists.second);
                }

                // Fake - Case 2 (Unused Only)
                if (!pt2.is_used) {
                    hists.fake_unused_pt2_deltaPT->Fill(pt2.delta_pt);
                    hists.fake_unused_pt2_deltaETA->Fill(pt2.delta_eta);
                    hists.fake_unused_pt2_deltaPHI->Fill(pt2.delta_phi);
                    if (dists.first >= 0) {
                        hists.fake_unused_pt2_dist3D->Fill(dists.first);
                        hists.fake_unused_pt2_dist3D->Fill(dists.second);
                    }
                }
            }
        } // End pT2 Loop
    } // End Event Looper

    // 7. Output and Plotting
    if (makePlots) {
        Plotting::plotPt2Histograms(hists, outputDir);
    }

    return 0;
}
