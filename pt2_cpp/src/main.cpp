#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <getopt.h>
#include <cmath>
#include <iomanip>
#include <memory>

#include "gator.h"
#include "histograms.h"
#include "plotting.h"
#include "plot_recipes.h"
#include "rootReader.h"
#include "tools.h"
#include "pt2.h"
#include "extrapolation.h"
#include "root_writer.h"
#include "model_inference.h"

int main(int argc, char** argv) {

    // Default Arguments
    bool makePlots = false;
    bool writeRoot = false;
    bool lowPT = false;
    int nEvents = -1;
    std::string inputFile; 
    std::string outputDir;

    // Command Line Arguments
    int opt;
    while ((opt = getopt(argc, argv, "prki:o:n:")) != -1) {
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
                    << "[-r] Make Root File for Training\n"
                    << "[-i] Input File Path \n"
                    << "[-o] Output Directory \n"
                    << "[-n] Number of Events \n";
                return 1;
        }    
    }

    // Code should do something
    if ((!makePlots) && (!writeRoot)){
        std::cerr << "Error: the code should have some sort of output!" << std::endl;
        return 1;
    }
    // Set input and output file paths
    if (outputDir.empty()){
        outputDir = "output";
    }
    if (inputFile.empty()){ 
        if (lowPT){
            inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple_LowPT.root";
        }
        else {
            inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple.root";
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

    // ML MODEL
    const std::string MODEL_PATH = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/NN_MODEL/model.onnx";
    const std::string MEAN_PATH = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/NN_MODEL/mean.npy";
    const std::string STD_PATH = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/NN_MODEL/std.npy";
    // ----- CUT SCORES -----
    //const float CUT_SCORE = 0.9976133704185486f; // 90% Real Efficiency
    //const float CUT_SCORE = 0.9933087825775146f; // 95% Real Efficiency
    //const float CUT_SCORE = 0.9894251227378845f; // 96% Real Efficiency
    //const float CUT_SCORE = 0.9768215417861938f; // 97% Real Efficiency
    const float CUT_SCORE = 0.9252402186393738f; // 98% Real Efficiency
    //const float CUT_SCORE = 0.7104159593582153f; // 99% Real Efficiency
    //const float CUT_SCORE = 0.011754416860640049f; // 99.9% Real Efficiency
    // Load the Model
    Pt2Scorer scorer(MODEL_PATH, MEAN_PATH, STD_PATH);
    std::cout << "ML model loaded.  Active cut score = " << CUT_SCORE << "\n\n";

    // Initialize Histograms
    HistogramManager hists;
    hists.init();

    // Get the Correct Pixel Map Directory
    std::string pixelMapFileDir;
    if (lowPT){
        pixelMapFileDir = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/PIXEL_MAPS/Pixel_Maps_0p6GeV/";
    }
    else{
        pixelMapFileDir = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/PIXEL_MAPS/Pixel_Maps_0p8GeV/";
    }

    // Load superbin --> detID
    SuperbinToDetIdMap superbinToDetIds_POS;
    SuperbinToDetIdMap superbinToDetIds_NEG;
    SuperbinToDetIdMap superbinToDetIds_NON;
    try {
        loadSuperbinDetIdMap(pixelMapFileDir, 
                            superbinToDetIds_POS, 
                            superbinToDetIds_NEG,
                            superbinToDetIds_NON);
    }
    catch (const std::exception& e) {
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

        // Update progress bar every N events
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

        // TODO: Does this double count pT2 candidates?
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
            reader.pls_superbin.push_back(CalculateSuperbin(reader, j, lowPT));
            buildPt2sForPLS(j, reader, superbinToDetIds_POS, superbinToDetIds_NEG, superbinToDetIds_NON, detidToLS, pt2s);
        }

        // Create a new writer for this event
        std::ostringstream filename;
        filename << outputDir << "/event_" << ievt << ".root";
        std::unique_ptr<Pt2Writer> writer = std::make_unique<Pt2Writer>(outputDir);
        writer->createNewFile(filename.str());

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
            std::vector<double> heli = extrapolation::extrapolatePlsHelicallyAndGetDistance(plsIdx, lsIdx, reader);
            std::pair<double, double> rz_simple = extrapolation::extrapolateSimplePointingInRZ(plsIdx, lsIdx, reader);
            double dAngle = extrapolation::calculateDeltaAngle(plsIdx, lsIdx, reader);
            // ML Score
            std::array<float, 21> features = {{
                reader.ls_pt->at(lsIdx),
                reader.ls_eta->at(lsIdx),
                std::sin(reader.ls_phi->at(lsIdx)),
                std::cos(reader.ls_phi->at(lsIdx)),
                reader.pls_pt->at(plsIdx),
                reader.pls_eta->at(plsIdx),
                std::sin(reader.pls_phi->at(plsIdx)),
                std::cos(reader.pls_phi->at(plsIdx)),
                static_cast<float>(reader.pls_charge->at(plsIdx)),
                static_cast<float>(reader.pls_nhit->at(plsIdx)),
                pt2.delta_pt,
                pt2.delta_eta,
                pt2.delta_phi,
                dR,
                static_cast<float>(heli[0]),
                static_cast<float>(heli[1]),
                static_cast<float>(heli[2]),
                static_cast<float>(heli[3]),
                static_cast<float>(rz_simple.first),
                static_cast<float>(rz_simple.second),
                std::log(std::abs(static_cast<float>(heli[0])) + 1e-6f),
            }};
            float mlscore = scorer.score(features);

            if (
                mlscore >= CUT_SCORE &&
                true
            ){
                if (writeRoot) {
                    writer->fill(
                        // Truth
                        pt2.is_real,
                        // LS
                        reader.ls_pt->at(lsIdx),
                        reader.ls_eta->at(lsIdx),
                        reader.ls_phi->at(lsIdx),
                        // M_0
                        reader.md_pt->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_eta->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_phi->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_anchor_x->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_anchor_y->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_anchor_z->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_other_x->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_other_y->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_other_z->at(reader.ls_mdIdx0->at(lsIdx)),
                        reader.md_layer->at(reader.ls_mdIdx0->at(lsIdx)),
                        // MD_1
                        reader.md_pt->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_eta->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_phi->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_anchor_x->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_anchor_y->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_anchor_z->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_other_x->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_other_y->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_other_z->at(reader.ls_mdIdx1->at(lsIdx)),
                        reader.md_layer->at(reader.ls_mdIdx1->at(lsIdx)),
                        // PLS
                        reader.pls_pt->at(plsIdx),
                        reader.pls_phi->at(plsIdx),
                        reader.pls_eta->at(plsIdx),
                        reader.pls_hit0_x->at(plsIdx),
                        reader.pls_hit0_y->at(plsIdx),
                        reader.pls_hit0_z->at(plsIdx),
                        reader.pls_hit1_x->at(plsIdx),
                        reader.pls_hit1_y->at(plsIdx),
                        reader.pls_hit1_z->at(plsIdx),
                        reader.pls_hit2_x->at(plsIdx),
                        reader.pls_hit2_y->at(plsIdx),
                        reader.pls_hit2_z->at(plsIdx),
                        reader.pls_hit3_x->at(plsIdx),
                        reader.pls_hit3_y->at(plsIdx),
                        reader.pls_hit3_z->at(plsIdx),
                        reader.pls_charge->at(plsIdx),
                        reader.pls_nhit->at(plsIdx),
                        // Calculated
                        pt2.delta_pt,
                        pt2.delta_eta,
                        pt2.delta_phi,
                        dR,
                        heli[0],
                        heli[1],
                        heli[2],
                        heli[3],
                        rz_simple.first,
                        rz_simple.second
                    );
                }                
                if (pt2.is_real) {
                    hists.real_pt2_deltaPT->Fill(pt2.delta_pt);
                    hists.real_pt2_deltaETA->Fill(pt2.delta_eta);
                    hists.real_pt2_deltaPHI->Fill(pt2.delta_phi);
                    hists.real_pt2_deltaR->Fill(dR);

                    if (dAngle > -1.0) hists.real_pt2_deltaAngle->Fill(dAngle);

                    // Fill Separated 3D components for Real
                    if (heli[0] >= 0) {
                        hists.real_pt2_MD0_dXY->Fill(heli[0]); 
                        hists.real_pt2_MD0_dZ->Fill(heli[1]);
                        hists.real_pt2_MD1_dXY->Fill(heli[2]); 
                        hists.real_pt2_MD1_dZ->Fill(heli[3]);
                        //Fill 2d histograms
                        hists.h2_real_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                        hists.h2_real_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);
                    }

                    if (rz_simple.first > -900) {
                        hists.real_pt2_MD0_rz_simple->Fill(rz_simple.first);
                        hists.real_pt2_MD1_rz_simple->Fill(rz_simple.second);
                    }

                    if (!pt2.is_used) {
                        hists.real_unused_pt2_deltaPT->Fill(pt2.delta_pt);
                        hists.real_unused_pt2_deltaETA->Fill(pt2.delta_eta);
                        hists.real_unused_pt2_deltaPHI->Fill(pt2.delta_phi);
                        hists.real_unused_pt2_deltaR->Fill(dR);
                        if (dAngle > -1.0) hists.real_unused_pt2_deltaAngle->Fill(dAngle);

                        if (heli[0] >= 0) {
                            hists.real_unused_pt2_MD0_dXY->Fill(heli[0]); 
                            hists.real_unused_pt2_MD0_dZ->Fill(heli[1]);
                            hists.real_unused_pt2_MD1_dXY->Fill(heli[2]); 
                            hists.real_unused_pt2_MD1_dZ->Fill(heli[3]);

                            //Fill 2d histograms
                            hists.h2_real_unused_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                            hists.h2_real_unused_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);
                        }

                        if (rz_simple.first > -900) {
                            hists.real_unused_pt2_MD0_rz_simple->Fill(rz_simple.first);
                            hists.real_unused_pt2_MD1_rz_simple->Fill(rz_simple.second);
                        }
                    }
                } 
                else {
                    hists.fake_pt2_deltaPT->Fill(pt2.delta_pt);
                    hists.fake_pt2_deltaETA->Fill(pt2.delta_eta);
                    hists.fake_pt2_deltaPHI->Fill(pt2.delta_phi);
                    hists.fake_pt2_deltaR->Fill(dR);
                    if (dAngle > -1.0) hists.fake_pt2_deltaAngle->Fill(dAngle);

                    // Fill Separated 3D components for Fake
                    if (heli[0] >= 0) {
                        hists.fake_pt2_MD0_dXY->Fill(heli[0]); 
                        hists.fake_pt2_MD0_dZ->Fill(heli[1]);
                        hists.fake_pt2_MD1_dXY->Fill(heli[2]); 
                        hists.fake_pt2_MD1_dZ->Fill(heli[3]);

                        //Fill 2d histograms
                        hists.h2_fake_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                        hists.h2_fake_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);
                    }

                    if (rz_simple.first > -900) {
                        hists.fake_pt2_MD0_rz_simple->Fill(rz_simple.first);
                        hists.fake_pt2_MD1_rz_simple->Fill(rz_simple.second);
                    }

                    if (!pt2.is_used) {
                        hists.fake_unused_pt2_deltaPT->Fill(pt2.delta_pt);
                        hists.fake_unused_pt2_deltaETA->Fill(pt2.delta_eta);
                        hists.fake_unused_pt2_deltaPHI->Fill(pt2.delta_phi);
                        hists.fake_unused_pt2_deltaR->Fill(dR);
                        if (dAngle > -1.0) hists.fake_unused_pt2_deltaAngle->Fill(dAngle);

                        if (heli[0] >= 0) {
                            hists.fake_unused_pt2_MD0_dXY->Fill(heli[0]); 
                            hists.fake_unused_pt2_MD0_dZ->Fill(heli[1]);
                            hists.fake_unused_pt2_MD1_dXY->Fill(heli[2]); 
                            hists.fake_unused_pt2_MD1_dZ->Fill(heli[3]);
                        
                            //Fill 2d histograms
                            hists.h2_fake_unused_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                            hists.h2_fake_unused_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);
                        }

                        if (rz_simple.first > -900) {
                            hists.fake_unused_pt2_MD0_rz_simple->Fill(rz_simple.first);
                            hists.fake_unused_pt2_MD1_rz_simple->Fill(rz_simple.second);
                        }
                    }
                }      
            }
        } 
        writer->close();
    }

    if (makePlots) {
        auto recipes = getPt2Recipes(hists);
        Plotting plotter; 
        plotter.plotRecipes(recipes, outputDir);
    }

    return 0;
}
