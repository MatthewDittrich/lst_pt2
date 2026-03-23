#include <iostream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <getopt.h>
#include <cmath>
#include <iomanip>

#include "gator.h"
#include "histograms.h"
#include "plotting.h"
#include "plot_recipes.h"
#include "rootReader.h"
#include "tools.h"
#include "pt2.h"
#include "extrapolation.h"
#include "extra_cuts.h"

int main(int argc, char** argv) {

    // Default Arguments
    bool makePlots = false;
    bool writeRoot = false;
    bool lowPT = false;
    int nEvents = -1;
    double targetPercent = 90.0; 
    double myCutZ0 = 4.4169;  
    double myCutZ1 = 5.7503;
    std::string inputFile; 
    std::string outputDir;

    // Command Line Arguments
    int opt;
    while ((opt = getopt(argc, argv, "prki:o:n:e:")) != -1) {
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
            case 'e': 
                targetPercent = std::stod(optarg); 
                break; 
            default:
                std::cerr << "Usage: " << argv[0] << "\n"
                    << "[-p] Make Plots\n"
                    << "[-k] Run Low pT\n"
                    << "[-r] Make Root File \n"
                    << "[-i] Input File Path \n"
                    << "[-o] Output Directory \n"
                    << "[-n] Number of Events \n"
                    << "[-e] Percentage of real events \n";
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

    // Initialize Histograms
    HistogramManager hists;
    hists.init();

    // Options that are not ready
    if (writeRoot) {
        throw std::runtime_error(
            "Error: options --writeRoot is not currently set up."
        );
    }

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
            
            double lst_dPhi = extra_cuts::calculateLSTDPhi(plsIdx, lsIdx, reader);
            std::vector<double> betas = extra_cuts::calculateLSTdBeta(plsIdx, lsIdx, reader);
            double dBeta = betas[2];
            double betaOut = betas[1]; 
            double lst_zResGeo = extra_cuts::calculateLSTOriginZResidual(plsIdx, lsIdx, reader);
            double lst_zResKin = extra_cuts::calculateLSTKinematicZResidual(plsIdx, lsIdx, reader);
           

            if(heli[1] > myCutZ0 ||  heli[3] > myCutZ1 ){continue;}
           // if(heli[0] > 1.3920 ||  heli[2] > 2.2394 ){continue;}
            //if(pt2.delta_phi < -0.2991 || pt2.delta_phi > 0.3039){continue;}
           // if(rz_simple.first < -1.7022  || rz_simple.first > 0.5405 || rz_simple.second < -2.9496 || rz_simple.second > 0.8334){continue;}
           // if(pt2.delta_pt < -1.2848 || pt2.delta_pt > 0.9443){continue;}
            if (pt2.is_real) {
                hists.real_pt2_deltaPT->Fill(pt2.delta_pt);
                hists.real_pt2_deltaETA->Fill(pt2.delta_eta);
                hists.real_pt2_deltaPHI->Fill(pt2.delta_phi);
                hists.real_pt2_deltaR->Fill(dR);

                if (dAngle > -1.0) hists.real_pt2_deltaAngle->Fill(dAngle);
                if (lst_dPhi > -100.0) hists.real_pt2_LSTdPhi->Fill(lst_dPhi);
                if (dBeta > -100.0) hists.real_pt2_LSTdBeta->Fill(dBeta);
                if (betaOut > -100.0) hists.real_pt2_LSTbetaOut->Fill(betaOut);
                if (lst_zResGeo > -100.0) hists.real_pt2_LSTOrgZRes->Fill(lst_zResGeo);
                if (lst_zResKin > -100.0) hists.real_pt2_LSTKinZRes->Fill(lst_zResKin);

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
                    if (lst_dPhi > -100.0) hists.real_unused_pt2_LSTdPhi->Fill(lst_dPhi);
                    if (dBeta > -100.0) hists.real_unused_pt2_LSTdBeta->Fill(dBeta);
                    if (betaOut > -100.0) hists.real_unused_pt2_LSTbetaOut->Fill(betaOut);
                    if (lst_zResGeo > -100.0) hists.real_unused_pt2_LSTOrgZRes->Fill(lst_zResGeo);
                    if (lst_zResKin > -100.0) hists.real_unused_pt2_LSTKinZRes->Fill(lst_zResKin);

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
                if (lst_dPhi > -100.0) hists.fake_pt2_LSTdPhi->Fill(lst_dPhi);
                if (dBeta > -100.0) hists.fake_pt2_LSTdBeta->Fill(dBeta);
                if (betaOut > -100.0) hists.fake_pt2_LSTbetaOut->Fill(betaOut);
                if (lst_zResGeo > -100.0) hists.fake_pt2_LSTOrgZRes->Fill(lst_zResGeo);
                if (lst_zResKin > -100.0) hists.fake_pt2_LSTKinZRes->Fill(lst_zResKin);

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
                    if (lst_dPhi > -100.0) hists.fake_unused_pt2_LSTdPhi->Fill(lst_dPhi);
                    if (dBeta > -100.0) hists.fake_unused_pt2_LSTdBeta->Fill(dBeta);
                    if (betaOut > -100.0) hists.fake_unused_pt2_LSTbetaOut->Fill(betaOut);
                    if (lst_zResGeo > -100.0) hists.fake_unused_pt2_LSTOrgZRes->Fill(lst_zResGeo);
                    if (lst_zResKin > -100.0) hists.fake_unused_pt2_LSTKinZRes->Fill(lst_zResKin);

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

    // --- AUTOMATED CUT CALCULATION ---
    double q[1];
    double p[1] = { targetPercent / 100.0 }; // e.g., 0.90
    double idealCut0 = 0, idealCut1 = 0;
    
    // 2-Sided Bounds Variables
    double idealCutPhiMin = 0, idealCutPhiMax = 0;
    double idealCutRZ0Min = 0, idealCutRZ0Max = 0;
    double idealCutRZ1Min = 0, idealCutRZ1Max = 0;
    double idealCutPtMin = 0,  idealCutPtMax = 0; // <--- NEW: Delta pT Bounds

    // For 2-sided variables
    double tail = (1.0 - (targetPercent / 100.0)) / 2.0; 
    double p_sym[2] = { tail, 1.0 - tail }; 
    double q_sym[2]; 

    if (hists.real_pt2_MD0_dZ->GetEntries() > 0) {
        
        // --- 1-Sided Cuts (dZ) ---
        hists.real_pt2_MD0_dZ->GetQuantiles(1, q, p);
        idealCut0 = q[0];

        hists.real_pt2_MD1_dZ->GetQuantiles(1, q, p);
        idealCut1 = q[0];

        // --- 2-Sided Cuts (Delta Phi) ---
        if (hists.real_pt2_deltaPHI->GetEntries() > 0) {
            hists.real_pt2_deltaPHI->GetQuantiles(2, q_sym, p_sym);
            idealCutPhiMin = q_sym[0]; 
            idealCutPhiMax = q_sym[1]; 
        }
        
        // --- 2-Sided Cuts (MD0 R-Z Simple) ---
        if (hists.real_pt2_MD0_rz_simple->GetEntries() > 0) {
            hists.real_pt2_MD0_rz_simple->GetQuantiles(2, q_sym, p_sym);
            idealCutRZ0Min = q_sym[0]; 
            idealCutRZ0Max = q_sym[1]; 
        }

        // --- 2-Sided Cuts (MD1 R-Z Simple) ---
        if (hists.real_pt2_MD1_rz_simple->GetEntries() > 0) {
            hists.real_pt2_MD1_rz_simple->GetQuantiles(2, q_sym, p_sym);
            idealCutRZ1Min = q_sym[0]; 
            idealCutRZ1Max = q_sym[1]; 
        }

        // --- 2-Sided Cuts (Delta pT) ---
        if (hists.real_pt2_deltaPT->GetEntries() > 0) {
            hists.real_pt2_deltaPT->GetQuantiles(2, q_sym, p_sym);
            idealCutPtMin = q_sym[0]; 
            idealCutPtMax = q_sym[1]; 
        }

        std::cout << "\n=================================================" << std::endl;
        std::cout << "Target Efficiency set to: " << targetPercent << "%" << std::endl;
        std::cout << "IDEAL CUTS TO REACH THIS EFFICIENCY:" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "MD0 dZ Cut:    < " << idealCut0 << " cm" << std::endl;
        std::cout << "MD1 dZ Cut:    < " << idealCut1 << " cm" << std::endl;
        std::cout << "Delta Phi Cut: " << idealCutPhiMin << "  to  " << idealCutPhiMax << " rad" << std::endl;
        std::cout << "MD0 R-Z Cut:   " << idealCutRZ0Min << "  to  " << idealCutRZ0Max << " cm" << std::endl;
        std::cout << "MD1 R-Z Cut:   " << idealCutRZ1Min << "  to  " << idealCutRZ1Max << " cm" << std::endl;
        std::cout << "Delta pT Cut:  " << idealCutPtMin  << "  to  " << idealCutPtMax  << " GeV" << std::endl;
        std::cout << "=================================================\n" << std::endl;
    }

    auto recipes = getPt2Recipes(hists);
    Plotting plotter; 
    plotter.plotRecipes(recipes, outputDir);


    if (writeRoot) {
    std::string outFileName = outputDir + "/cut_study_hists.root";
    TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

    hists.h2_real_MD0_dXY_vs_dZ->Write();
    hists.h2_fake_MD0_dXY_vs_dZ->Write();
    hists.h2_real_MD1_dXY_vs_dZ->Write();
    hists.h2_fake_MD1_dXY_vs_dZ->Write();

    hists.h2_real_unused_MD0_dXY_vs_dZ->Write();
    hists.h2_fake_unused_MD0_dXY_vs_dZ->Write();
    hists.h2_real_unused_MD1_dXY_vs_dZ->Write();
    hists.h2_fake_unused_MD1_dXY_vs_dZ->Write();

    outFile->Close();
    std::cout << "Saved 2D cut study histograms to: " << outFileName << std::endl;
    }

    return 0;
}
