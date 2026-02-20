#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <getopt.h>
#include <cmath>

#include "gator.h"
#include "histograms.h"
#include "plotting.h"
#include "plot_recipes.h"
#include "rootReader.h"
#include "tools.h"
#include "pt2.h"
#include "extrapolation.h"

int main(int argc, char** argv) {

    // Default Arguments
    bool makePlots = false;
    bool writeRoot = false;
    bool lowPT = false;
    int nEvents = -1;
    std::string inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple.root";
    std::string outputDir = "output";

    // Command Line Arguments
    int opt;
    while ((opt = getopt(argc, argv, "pi:o:r:k:n:")) != -1) {
        switch (opt) {
            case 'p':
                makePlots = true;
                break;
            case 'r':
                writeRoot = true;
                break;
            case 'k':
                lowPT = true;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputDir = optarg;
                break;
            case 'n':
                nEvents = std::stoi(optarg);
                break;
            default:
                std::cerr << "Usage: " << argv[0] << "\n"
                    << "[-p] Make Plots\n"
                    << "[-k] Run Low pT\n"
                    << "[-r] Make Root File \n"
                    << "[-i] Input File Path \n"
                    << "[-o] Output Directory \n"
                    << "[-n] Number of Events \n";
                return 1;
        }    
    }

    // Print Gator
    print_gator();

    // Print Configuration
    std::cout << "\n=== Configuration ===\n";
    std::cout << "Input file:        " << inputFile << "\n";
    std::cout << "Output directory:  " << outputDir << "\n";
    std::cout << "Make plots:        " << (makePlots ? "yes" : "no") << "\n";
    std::cout << "Use Low pT:        " << (lowPT ? "yes" : "no") << "\n";
    std::cout << "Write ROOT file:   " << (writeRoot ? "yes" : "no") << "\n";
    if (nEvents > 0) std::cout << "Number of events:  " << nEvents << "\n";
    std::cout << "=====================\n\n";

    // Create the output directory
    std::filesystem::create_directories(outputDir);

    // Initialize Histograms
    HistogramManager hists;
    hists.init();

    // Options that are not ready
    if (lowPT || writeRoot) {
        throw std::runtime_error(
            "Error: options --lowPT and --writeRoot are not currently set up."
        );
    }

    // Get the pixelMap loaded
    SuperbinToDetIdMap superbinToDetIds;
    std::string pixelMapFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/PIXEL_MAPS/Standard_0p8GeV_Maps/pLS_map_ElCheapo.txt";
    try {
        loadSuperbinDetIdMap(pixelMapFile, superbinToDetIds);
    } catch (const std::exception& e) {
        std::cerr << "Error loading superbin map: " << e.what() << std::endl;
        return 1;
    }

    // Initialize per-event variables
    DetIdToLSMap detidToLS;
    pT2Collection pt2s;

    // Setup the Root Reader
    rootReader reader;
    if (!reader.Init(inputFile, "tree")) {
        std::cerr << "Error: Could not find Root File or TTree! \n" << std::endl;
        return 1;
    } 

    Long64_t totalEntries = reader.GetEntries();
    Long64_t entriesToProcess = (nEvents > 0 && nEvents < totalEntries) ? nEvents : totalEntries;

    print_creature();

    // Main Looper
    for (Long64_t ievt = 0; ievt < entriesToProcess; ++ievt) {

        reader.GetEntry(ievt);

        // Update progress bar every 5 events
        if (ievt % 2 == 0 || ievt == entriesToProcess)
            printProgressBar(ievt, entriesToProcess);

        // Clear Variables from the last event
        reader.pls_origin_z.clear();
        reader.pls_superbin.clear();
        detidToLS.clear();
        pt2s.clear();
        reader.pls_isUsed.clear();
        reader.ls_isUsed.clear();

        // Get used masks for the LS and pLS
        UsedMask usedMask = buildUsedMask(reader);
        reader.ls_isUsed  = std::move(usedMask.ls_isUsed);
        reader.pls_isUsed = std::move(usedMask.pls_isUsed);

        // Get Size of pLS and LS
        size_t nLS = reader.ls_pt->size();
        size_t nPLS = reader.pls_pt->size();

        // Reserve New Calculated Variables
        reader.pls_origin_z.reserve(nPLS);
        reader.pls_superbin.reserve(nPLS);

        // LS Loop
        for (size_t k = 0; k < nLS; ++k) {
            std::vector<int> detIds = getDetIdsForLS(reader, k);
                for (int detId : detIds) {
                    detidToLS[detId].push_back(k);
                }
        }

        // pLS Loop
        for (size_t j = 0; j < nPLS; ++j) {
            reader.pls_origin_z.push_back(CalculatePlsZ(reader, j));
            reader.pls_superbin.push_back(CalculateSuperbin(reader, j));
            buildPt2sForPLS(j, reader, superbinToDetIds, detidToLS, pt2s);
        }

        // pT2 Loop
        for (auto& pt2 : pt2s) {
            size_t plsIdx = pt2.pls_idx;
            size_t lsIdx = pt2.ls_idx;
            
            pt2.delta_pt  = deltaPt(reader.pls_pt->at(plsIdx), reader.ls_pt->at(lsIdx));
            pt2.delta_eta = deltaEta(reader.pls_eta->at(plsIdx), reader.ls_eta->at(lsIdx));
            pt2.delta_phi = deltaPhi(reader.pls_phi->at(plsIdx), reader.ls_phi->at(lsIdx));
            pt2.is_real = pt2TruthFinder(reader, plsIdx, lsIdx);
            pt2.is_used = pt2UsedCalculator(reader, plsIdx, lsIdx);

            // New Physics Calculations
            float dR = std::sqrt(pt2.delta_eta * pt2.delta_eta + pt2.delta_phi * pt2.delta_phi);
            std::pair<double, double> dists = extrapolation::extrapolatePlsHelicallyAndGetDistance(plsIdx, lsIdx, reader);
            std::pair<double, double> rz_simple = extrapolation::extrapolateSimplePointingInRZ(plsIdx, lsIdx, reader);
            double dAngle = extrapolation::calculateDeltaAngle(plsIdx, lsIdx, reader);

            // Fill some histos
            if (pt2.is_real) {
                hists.real_pt2_deltaPT->Fill(pt2.delta_pt);
                hists.real_pt2_deltaETA->Fill(pt2.delta_eta);
                hists.real_pt2_deltaPHI->Fill(pt2.delta_phi);
                hists.real_pt2_deltaR->Fill(dR);

                if (dAngle > -1.0) hists.real_pt2_deltaAngle->Fill(dAngle);

                if (rz_simple.first > -900) {
                    hists.real_pt2_rz_simple->Fill(rz_simple.first);
                    hists.real_pt2_rz_simple->Fill(rz_simple.second);
                }

                if (dists.first >= 0) {
                    hists.real_pt2_dist3D->Fill(dists.first);
                    hists.real_pt2_dist3D->Fill(dists.second);
                }

                if (!pt2.is_used) {
                    hists.real_unused_pt2_deltaPT->Fill(pt2.delta_pt);
                    hists.real_unused_pt2_deltaETA->Fill(pt2.delta_eta);
                    hists.real_unused_pt2_deltaPHI->Fill(pt2.delta_phi);
                    hists.real_unused_pt2_deltaR->Fill(dR);

                    if (dAngle > -1.0) hists.real_unused_pt2_deltaAngle->Fill(dAngle);

                    if (rz_simple.first > -900) {
                        hists.real_unused_pt2_rz_simple->Fill(rz_simple.first);
                        hists.real_unused_pt2_rz_simple->Fill(rz_simple.second);
                    }

                    if (dists.first >= 0) {
                        hists.real_unused_pt2_dist3D->Fill(dists.first);
                        hists.real_unused_pt2_dist3D->Fill(dists.second);
                    }
                }
            } 
            else {
                hists.fake_pt2_deltaPT->Fill(pt2.delta_pt);
                hists.fake_pt2_deltaETA->Fill(pt2.delta_eta);
                hists.fake_pt2_deltaPHI->Fill(pt2.delta_phi);
                hists.fake_pt2_deltaR->Fill(dR);

                if (dAngle > -1.0) hists.fake_pt2_deltaAngle->Fill(dAngle);

                if (rz_simple.first > -900) {
                    hists.fake_pt2_rz_simple->Fill(rz_simple.first);
                    hists.fake_pt2_rz_simple->Fill(rz_simple.second);
                }

                if (dists.first >= 0) {
                    hists.fake_pt2_dist3D->Fill(dists.first);
                    hists.fake_pt2_dist3D->Fill(dists.second);
                }

                if (!pt2.is_used) {
                    hists.fake_unused_pt2_deltaPT->Fill(pt2.delta_pt);
                    hists.fake_unused_pt2_deltaETA->Fill(pt2.delta_eta);
                    hists.fake_unused_pt2_deltaPHI->Fill(pt2.delta_phi);
                    hists.fake_unused_pt2_deltaR->Fill(dR);

                    if (dAngle > -1.0) hists.fake_unused_pt2_deltaAngle->Fill(dAngle);

                    if (rz_simple.first > -900) {
                        hists.fake_unused_pt2_rz_simple->Fill(rz_simple.first);
                        hists.fake_unused_pt2_rz_simple->Fill(rz_simple.second);
                    }

                    if (dists.first >= 0) {
                        hists.fake_unused_pt2_dist3D->Fill(dists.first);
                        hists.fake_unused_pt2_dist3D->Fill(dists.second);
                    }
                }
            }
        }
    }

    auto recipes = getPt2Recipes(hists);
    Plotting plotter; 
    plotter.plotRecipes(recipes, outputDir);

    return 0;
}
