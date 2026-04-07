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
    double myCutZ0 = 0.4896;  
    double myCutZ1 = 0.9304;
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
       
    // =========================================================================
    //                        SETUP BDT OUTPUT TREE
    // =========================================================================
    std::string bdtFileName = outputDir + "/bdt_training_data.root";
    TFile* bdtFile = new TFile(bdtFileName.c_str(), "RECREATE");
    TTree* bdtTree = new TTree("pt2_features", "LST Variables for BDT Training");

    // Variables to hold the data for each row
    double b_lst_dPhi, b_betaIn, b_betaOut, b_dBeta;
    double b_zResGeo, b_zResKin;
    double b_delta_pt, b_delta_eta, b_delta_phi;
    int b_is_real;
    double b_heli_dXY0, b_heli_dZ0, b_heli_dXY1, b_heli_dZ1;
    double b_rz_simple0, b_rz_simple1;
    int b_event_id;
    int b_md0_cat, b_md1_cat;
    
    bdtTree->Branch("event_id", &b_event_id);
    bdtTree->Branch("md0_cat", &b_md0_cat); 
    bdtTree->Branch("md1_cat", &b_md1_cat); 
    // Create the branches (columns) in the TTree
    bdtTree->Branch("lst_dPhi", &b_lst_dPhi);
    bdtTree->Branch("betaIn", &b_betaIn);
    bdtTree->Branch("betaOut", &b_betaOut);
    bdtTree->Branch("dBeta", &b_dBeta);
    bdtTree->Branch("zResGeo", &b_zResGeo);
    bdtTree->Branch("zResKin", &b_zResKin);
    bdtTree->Branch("delta_pt", &b_delta_pt);
    bdtTree->Branch("delta_eta", &b_delta_eta);
    bdtTree->Branch("delta_phi", &b_delta_phi);
    bdtTree->Branch("heli_dXY0", &b_heli_dXY0);
    bdtTree->Branch("heli_dZ0", &b_heli_dZ0);
    bdtTree->Branch("heli_dXY1", &b_heli_dXY1);
    bdtTree->Branch("heli_dZ1", &b_heli_dZ1);
    bdtTree->Branch("rz_simple0", &b_rz_simple0);
    bdtTree->Branch("rz_simple1", &b_rz_simple1);

    
    // The "Target" variable: 1 for Real, 0 for Fake
    bdtTree->Branch("is_real", &b_is_real);
    
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
           
            // --- CATEGORIZE MD0 and MD1 ---
            int md0_idx = reader.ls_mdIdx0->at(lsIdx);
            int md1_idx = reader.ls_mdIdx1->at(lsIdx);
            uint32_t detId0 = reader.md_detId->at(md0_idx);
            uint32_t detId1 = reader.md_detId->at(md1_idx);
            
            int cat0 = extra_cuts::getCategoryFromDetId(detId0);
            int cat1 = extra_cuts::getCategoryFromDetId(detId1);
            int comboIdx = extra_cuts::getConnectionIndex(cat0, cat1); 

            pt2.delta_pt  = deltaPt(reader.pls_pt->at(plsIdx), reader.ls_pt->at(lsIdx));
            pt2.delta_eta = deltaEta(reader.pls_eta->at(plsIdx), reader.ls_eta->at(lsIdx));
            pt2.delta_phi = deltaPhi(reader.pls_phi->at(plsIdx), reader.ls_phi->at(lsIdx));
            pt2.is_real = pt2TruthFinder(reader, plsIdx, lsIdx);
            pt2.is_used = pt2UsedCalculator(reader, plsIdx, lsIdx);
            if((reader.pls_pt->at(plsIdx)) > 0.8){continue;}

            // New Physics Calculations
            float dR = std::sqrt(pt2.delta_eta * pt2.delta_eta + pt2.delta_phi * pt2.delta_phi);
            std::vector<double> heli = extrapolation::extrapolatePlsHelicallyAndGetDistance(plsIdx, lsIdx, reader);
            std::pair<double, double> rz_simple = extrapolation::extrapolateSimplePointingInRZ(plsIdx, lsIdx, reader);
            double dAngle = extrapolation::calculateDeltaAngle(plsIdx, lsIdx, reader);
            
            double lst_dPhi = extra_cuts::calculateLSTDPhi(plsIdx, lsIdx, reader);
            std::vector<double> betas = extra_cuts::calculateLSTdBeta(plsIdx, lsIdx, reader);
            double dBeta = betas[2];
            double betaOut = betas[1]; 
            double betaIn = betas[0];
            double lst_zResGeo = extra_cuts::calculateLSTOriginZResidual(plsIdx, lsIdx, reader);
            double lst_zResKin = extra_cuts::calculateLSTKinematicZResidual(plsIdx, lsIdx, reader);
            

            //---FILLING ROOT FILE FOR BDT -------
            if (lst_dPhi > -100.0 && dBeta > -100.0) {
                b_lst_dPhi  = lst_dPhi;
                b_betaIn    = betaIn;
                b_betaOut   = betaOut;
                b_dBeta     = dBeta;
                b_zResGeo   = lst_zResGeo;
                b_zResKin   = lst_zResKin;
                b_delta_pt  = pt2.delta_pt;
                b_delta_eta = pt2.delta_eta;
                b_delta_phi = pt2.delta_phi;
                b_heli_dXY0 = heli[0];
                b_heli_dZ0  = heli[1];
                b_heli_dXY1 = heli[2];
                b_heli_dZ1  = heli[3];
                b_rz_simple0 = rz_simple.first;
                b_rz_simple1 = rz_simple.second;
                b_event_id = ievt;
                // The Machine Learning label!
                b_is_real   = pt2.is_real ? 1 : 0;
                b_md0_cat   = cat0; // Save categories to Python!
                b_md1_cat   = cat1; 
                bdtTree->Fill(); 
            }
            //if(heli[1] > myCutZ0 ||  heli[3] > myCutZ1 ){continue;}
            //if(heli[0] > 2.3896 ||  heli[2] > 3.4234 ){continue;}
            //if(pt2.delta_phi < -0.3493 || pt2.delta_phi > 0.3457){continue;}
            //if(rz_simple.first < -2.8875  || rz_simple.first > 1.2586 || rz_simple.second < -4.7368 || rz_simple.second > 1.8688){continue;}
            //if(pt2.delta_pt < -0.6123 || pt2.delta_pt > 0.1846){continue;}
            //if(dBeta < -0.0445 || dBeta > 0.0393){continue;}
            //if(lst_zResKin < -3.5895 || lst_zResKin > 3.7145){continue;}
            if (comboIdx >= 0) {
            if (pt2.is_real) {
                hists.real_pt2_deltaPT[comboIdx]->Fill(pt2.delta_pt);
                hists.real_pt2_deltaETA[comboIdx]->Fill(pt2.delta_eta);
                hists.real_pt2_deltaPHI[comboIdx]->Fill(pt2.delta_phi);
                hists.real_pt2_deltaR[comboIdx]->Fill(dR);

                if (dAngle > -1.0) hists.real_pt2_deltaAngle[comboIdx]->Fill(dAngle);
                if (lst_dPhi > -100.0) hists.real_pt2_LSTdPhi[comboIdx]->Fill(lst_dPhi);
                if (dBeta > -100.0) hists.real_pt2_LSTdBeta[comboIdx]->Fill(dBeta);
                if (betaOut > -100.0) hists.real_pt2_LSTbetaOut[comboIdx]->Fill(betaOut);
                if (lst_zResGeo > -100.0) hists.real_pt2_LSTOrgZRes[comboIdx]->Fill(lst_zResGeo);
                if (lst_zResKin > -100.0) hists.real_pt2_LSTKinZRes[comboIdx]->Fill(lst_zResKin);

                // Fill Separated 3D components for Real
                if (heli[0] >= 0) {
                    hists.real_pt2_MD0_dXY[comboIdx]->Fill(heli[0]); 
                    hists.real_pt2_MD0_dZ[comboIdx]->Fill(heli[1]);
                    hists.real_pt2_MD1_dXY[comboIdx]->Fill(heli[2]); 
                    hists.real_pt2_MD1_dZ[comboIdx]->Fill(heli[3]);
                    //Fill 2d histograms
                    //hists.h2_real_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                    //hists.h2_real_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);

                }

                if (rz_simple.first > -900) {
                    hists.real_pt2_MD0_rz_simple[comboIdx]->Fill(rz_simple.first);
                    hists.real_pt2_MD1_rz_simple[comboIdx]->Fill(rz_simple.second);
                }

                if (!pt2.is_used) {
                    hists.real_unused_pt2_deltaPT[comboIdx]->Fill(pt2.delta_pt);
                    hists.real_unused_pt2_deltaETA[comboIdx]->Fill(pt2.delta_eta);
                    hists.real_unused_pt2_deltaPHI[comboIdx]->Fill(pt2.delta_phi);
                    hists.real_unused_pt2_deltaR[comboIdx]->Fill(dR);
                    if (dAngle > -1.0) hists.real_unused_pt2_deltaAngle[comboIdx]->Fill(dAngle);
                    if (lst_dPhi > -100.0) hists.real_unused_pt2_LSTdPhi[comboIdx]->Fill(lst_dPhi);
                    if (dBeta > -100.0) hists.real_unused_pt2_LSTdBeta[comboIdx]->Fill(dBeta);
                    if (betaOut > -100.0) hists.real_unused_pt2_LSTbetaOut[comboIdx]->Fill(betaOut);
                    if (lst_zResGeo > -100.0) hists.real_unused_pt2_LSTOrgZRes[comboIdx]->Fill(lst_zResGeo);
                    if (lst_zResKin > -100.0) hists.real_unused_pt2_LSTKinZRes[comboIdx]->Fill(lst_zResKin);

                    if (heli[0] >= 0) {
                        hists.real_unused_pt2_MD0_dXY[comboIdx]->Fill(heli[0]); 
                        hists.real_unused_pt2_MD0_dZ[comboIdx]->Fill(heli[1]);
                        hists.real_unused_pt2_MD1_dXY[comboIdx]->Fill(heli[2]); 
                        hists.real_unused_pt2_MD1_dZ[comboIdx]->Fill(heli[3]);

                        //Fill 2d histograms
                        //hists.h2_real_unused_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                        //hists.h2_real_unused_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);

                    }

                    if (rz_simple.first > -900) {
                        hists.real_unused_pt2_MD0_rz_simple[comboIdx]->Fill(rz_simple.first);
                        hists.real_unused_pt2_MD1_rz_simple[comboIdx]->Fill(rz_simple.second);
                    }
                }
            } 
            else {
                hists.fake_pt2_deltaPT[comboIdx]->Fill(pt2.delta_pt);
                hists.fake_pt2_deltaETA[comboIdx]->Fill(pt2.delta_eta);
                hists.fake_pt2_deltaPHI[comboIdx]->Fill(pt2.delta_phi);
                hists.fake_pt2_deltaR[comboIdx]->Fill(dR);

                if (dAngle > -1.0) hists.fake_pt2_deltaAngle[comboIdx]->Fill(dAngle);
                if (lst_dPhi > -100.0) hists.fake_pt2_LSTdPhi[comboIdx]->Fill(lst_dPhi);
                if (dBeta > -100.0) hists.fake_pt2_LSTdBeta[comboIdx]->Fill(dBeta);
                if (betaOut > -100.0) hists.fake_pt2_LSTbetaOut[comboIdx]->Fill(betaOut);
                if (lst_zResGeo > -100.0) hists.fake_pt2_LSTOrgZRes[comboIdx]->Fill(lst_zResGeo);
                if (lst_zResKin > -100.0) hists.fake_pt2_LSTKinZRes[comboIdx]->Fill(lst_zResKin);

                // Fill Separated 3D components for Fake
                if (heli[0] >= 0) {
                    hists.fake_pt2_MD0_dXY[comboIdx]->Fill(heli[0]); 
                    hists.fake_pt2_MD0_dZ[comboIdx]->Fill(heli[1]);
                    hists.fake_pt2_MD1_dXY[comboIdx]->Fill(heli[2]); 
                    hists.fake_pt2_MD1_dZ[comboIdx]->Fill(heli[3]);

                    //Fill 2d histograms
                    //hists.h2_fake_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                    //hists.h2_fake_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);

                }

                if (rz_simple.first > -900) {
                    hists.fake_pt2_MD0_rz_simple[comboIdx]->Fill(rz_simple.first);
                    hists.fake_pt2_MD1_rz_simple[comboIdx]->Fill(rz_simple.second);
                }

                if (!pt2.is_used) {
                    hists.fake_unused_pt2_deltaPT[comboIdx]->Fill(pt2.delta_pt);
                    hists.fake_unused_pt2_deltaETA[comboIdx]->Fill(pt2.delta_eta);
                    hists.fake_unused_pt2_deltaPHI[comboIdx]->Fill(pt2.delta_phi);
                    hists.fake_unused_pt2_deltaR[comboIdx]->Fill(dR);

                    if (dAngle > -1.0) hists.fake_unused_pt2_deltaAngle[comboIdx]->Fill(dAngle);
                    if (lst_dPhi > -100.0) hists.fake_unused_pt2_LSTdPhi[comboIdx]->Fill(lst_dPhi);
                    if (dBeta > -100.0) hists.fake_unused_pt2_LSTdBeta[comboIdx]->Fill(dBeta);
                    if (betaOut > -100.0) hists.fake_unused_pt2_LSTbetaOut[comboIdx]->Fill(betaOut);
                    if (lst_zResGeo > -100.0) hists.fake_unused_pt2_LSTOrgZRes[comboIdx]->Fill(lst_zResGeo);
                    if (lst_zResKin > -100.0) hists.fake_unused_pt2_LSTKinZRes[comboIdx]->Fill(lst_zResKin);

                    if (heli[0] >= 0) {
                        hists.fake_unused_pt2_MD0_dXY[comboIdx]->Fill(heli[0]); 
                        hists.fake_unused_pt2_MD0_dZ[comboIdx]->Fill(heli[1]);
                        hists.fake_unused_pt2_MD1_dXY[comboIdx]->Fill(heli[2]); 
                        hists.fake_unused_pt2_MD1_dZ[comboIdx]->Fill(heli[3]);
                        
                         //Fill 2d histograms
                       // hists.h2_fake_unused_MD0_dXY_vs_dZ->Fill(std::abs(heli[1]), heli[0]);
                        //hists.h2_fake_unused_MD1_dXY_vs_dZ->Fill(std::abs(heli[3]), heli[2]);
      
                    }

                    if (rz_simple.first > -900) {
                        hists.fake_unused_pt2_MD0_rz_simple[comboIdx]->Fill(rz_simple.first);
                        hists.fake_unused_pt2_MD1_rz_simple[comboIdx]->Fill(rz_simple.second);
                    }
                }
            }      
        } 
    }
    }
    // --- AUTOMATED CUT CALCULATION ---
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "IDEAL CUTS PER CATEGORY (Target Efficiency: " << targetPercent << "%)" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    for (int i = 0; i < 11; ++i) {
        // Skip category if it has no entries to avoid division by zero / ROOT errors
        if (hists.real_pt2_deltaPT[i]->GetEntries() == 0) continue;

        // 1-Sided configuration (upper bound only)
        double q_1[1];
        double p_1[1] = { targetPercent / 100.0 };

        // 2-Sided configuration (symmetric tails)
        double tail = (1.0 - (targetPercent / 100.0)) / 2.0;
        double q_2[2];
        double p_2[2] = { tail, 1.0 - tail };

        // --- Extrapolation Variables (1-Sided) ---
        double cut_dZ0 = 0, cut_dZ1 = 0;
        double cut_dXY0 = 0, cut_dXY1 = 0;

        // --- Simple R-Z & Kinematic Variables (2-Sided) ---
        double cut_RZ0Min = 0, cut_RZ0Max = 0;
        double cut_RZ1Min = 0, cut_RZ1Max = 0;
        double cut_dPhiMin = 0, cut_dPhiMax = 0;
        double cut_dPtMin = 0,  cut_dPtMax = 0;

        // --- LST Variables (2-Sided) ---
        double cut_LSTdPhiMin = 0,    cut_LSTdPhiMax = 0;
        double cut_LSTdBetaMin = 0,   cut_LSTdBetaMax = 0;
        double cut_LSTbetaOutMin = 0, cut_LSTbetaOutMax = 0;
        double cut_LSTOrgZMin = 0,    cut_LSTOrgZMax = 0;
        double cut_LSTKinZMin = 0,    cut_LSTKinZMax = 0;

        // --- Calculate 1-Sided Quantiles ---
        if (hists.real_pt2_MD0_dZ[i]->GetEntries() > 0) {
            hists.real_pt2_MD0_dZ[i]->GetQuantiles(1, q_1, p_1);
            cut_dZ0 = q_1[0];
        }
        if (hists.real_pt2_MD1_dZ[i]->GetEntries() > 0) {
            hists.real_pt2_MD1_dZ[i]->GetQuantiles(1, q_1, p_1);
            cut_dZ1 = q_1[0];
        }
        if (hists.real_pt2_MD0_dXY[i]->GetEntries() > 0) {
            hists.real_pt2_MD0_dXY[i]->GetQuantiles(1, q_1, p_1);
            cut_dXY0 = q_1[0];
        }
        if (hists.real_pt2_MD1_dXY[i]->GetEntries() > 0) {
            hists.real_pt2_MD1_dXY[i]->GetQuantiles(1, q_1, p_1);
            cut_dXY1 = q_1[0];
        }

        // --- Calculate 2-Sided Quantiles (Standard) ---
        if (hists.real_pt2_MD0_rz_simple[i]->GetEntries() > 0) {
            hists.real_pt2_MD0_rz_simple[i]->GetQuantiles(2, q_2, p_2);
            cut_RZ0Min = q_2[0]; cut_RZ0Max = q_2[1];
        }
        if (hists.real_pt2_MD1_rz_simple[i]->GetEntries() > 0) {
            hists.real_pt2_MD1_rz_simple[i]->GetQuantiles(2, q_2, p_2);
            cut_RZ1Min = q_2[0]; cut_RZ1Max = q_2[1];
        }
        if (hists.real_pt2_deltaPHI[i]->GetEntries() > 0) {
            hists.real_pt2_deltaPHI[i]->GetQuantiles(2, q_2, p_2);
            cut_dPhiMin = q_2[0]; cut_dPhiMax = q_2[1];
        }
        if (hists.real_pt2_deltaPT[i]->GetEntries() > 0) {
            hists.real_pt2_deltaPT[i]->GetQuantiles(2, q_2, p_2);
            cut_dPtMin = q_2[0]; cut_dPtMax = q_2[1];
        }

        // --- Calculate 2-Sided Quantiles (LST) ---
        if (hists.real_pt2_LSTdPhi[i]->GetEntries() > 0) {
            hists.real_pt2_LSTdPhi[i]->GetQuantiles(2, q_2, p_2);
            cut_LSTdPhiMin = q_2[0]; cut_LSTdPhiMax = q_2[1];
        }
        if (hists.real_pt2_LSTdBeta[i]->GetEntries() > 0) {
            hists.real_pt2_LSTdBeta[i]->GetQuantiles(2, q_2, p_2);
            cut_LSTdBetaMin = q_2[0]; cut_LSTdBetaMax = q_2[1];
        }
        if (hists.real_pt2_LSTbetaOut[i]->GetEntries() > 0) {
            hists.real_pt2_LSTbetaOut[i]->GetQuantiles(2, q_2, p_2);
            cut_LSTbetaOutMin = q_2[0]; cut_LSTbetaOutMax = q_2[1];
        }
        if (hists.real_pt2_LSTOrgZRes[i]->GetEntries() > 0) {
            hists.real_pt2_LSTOrgZRes[i]->GetQuantiles(2, q_2, p_2);
            cut_LSTOrgZMin = q_2[0]; cut_LSTOrgZMax = q_2[1];
        }
        if (hists.real_pt2_LSTKinZRes[i]->GetEntries() > 0) {
            hists.real_pt2_LSTKinZRes[i]->GetQuantiles(2, q_2, p_2);
            cut_LSTKinZMin = q_2[0]; cut_LSTKinZMax = q_2[1];
        }

        // --- Print Formatting ---
        std::cout << "\n>>> CATEGORY: " << hists.catTitles[i] << " <<<" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        
        std::cout << "  [Helical Extrapolation Cuts]" << std::endl;
        std::cout << "    MD0 dZ Cut:          < " << cut_dZ0 << " cm" << std::endl;
        std::cout << "    MD1 dZ Cut:          < " << cut_dZ1 << " cm" << std::endl;
        std::cout << "    MD0 dXY Cut:         < " << cut_dXY0 << " cm" << std::endl;
        std::cout << "    MD1 dXY Cut:         < " << cut_dXY1 << " cm" << std::endl;
        
        std::cout << "  [Simple Pointing Cuts]" << std::endl;
        std::cout << "    MD0 R-Z Res:           " << std::setw(8) << cut_RZ0Min << " to " << cut_RZ0Max << " cm" << std::endl;
        std::cout << "    MD1 R-Z Res:           " << std::setw(8) << cut_RZ1Min << " to " << cut_RZ1Max << " cm" << std::endl;
        
        std::cout << "  [Kinematic Cuts]" << std::endl;
        std::cout << "    Delta Phi:             " << std::setw(8) << cut_dPhiMin << " to " << cut_dPhiMax << " rad" << std::endl;
        std::cout << "    Delta pT:              " << std::setw(8) << cut_dPtMin << " to " << cut_dPtMax << " GeV" << std::endl;
        
        std::cout << "  [LST Component Cuts]" << std::endl;
        std::cout << "    LST Delta Phi:         " << std::setw(8) << cut_LSTdPhiMin << " to " << cut_LSTdPhiMax << " rad" << std::endl;
        std::cout << "    LST Delta Beta:        " << std::setw(8) << cut_LSTdBetaMin << " to " << cut_LSTdBetaMax << " rad" << std::endl;
        std::cout << "    LST Beta Out:          " << std::setw(8) << cut_LSTbetaOutMin << " to " << cut_LSTbetaOutMax << " rad" << std::endl;
        std::cout << "    LST Geometric Z-Res:   " << std::setw(8) << cut_LSTOrgZMin << " to " << cut_LSTOrgZMax << " cm" << std::endl;
        std::cout << "    LST Kinematic Z-Res:   " << std::setw(8) << cut_LSTKinZMin << " to " << cut_LSTKinZMax << " cm" << std::endl;
    }
    std::cout << "\n" << std::string(80, '=') << "\n" << std::endl;


    auto recipes = getPt2Recipes(hists);
    Plotting plotter; 
    plotter.plotRecipes(recipes, outputDir);

    bdtFile->cd();      
    bdtTree->Write();   
    bdtFile->Close();

    return 0;
}
