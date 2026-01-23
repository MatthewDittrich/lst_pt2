#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>

// ROOT Includes
#include <TStyle.h>
#include <TROOT.h>

// Project Includes
#include "LSTReader.h"
#include "PhysicsAlgo.h"
#include "HistManager.h"
#include "Plotter.h"
#include "Gator.h" // Ensure Gator.h and Gator.cpp exist as discussed

int main(int argc, char** argv) {
    // -------------------------------------------------------
    // 1. Configuration & Initialization
    // -------------------------------------------------------
    
    // Print the Gator!
    print_coding_gator();

    // Allow filename argument, default to LSTNtuple.root
    std::string input_filename = (argc > 1) ? argv[1] : "LSTNtuple.root";
    std::string tree_name = "tree";

    std::cout << "--- Starting LST Analysis ---" << std::endl;
    std::cout << "Opening file: " << input_filename << std::endl;

    // Initialize Reader
    LSTReader reader;
    if (!reader.Init(input_filename, tree_name)) {
        std::cerr << "Error: Could not initialize reader (File or Tree not found)." << std::endl;
        return 1;
    }

    // Initialize Histograms
    HistManager hists;
    hists.Init();

    // -------------------------------------------------------
    // 2. Event Loop
    // -------------------------------------------------------
    Long64_t nEntries = reader.GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;

    long long ideal_pt2_count = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        if (i % 1000 == 0) std::cout << "Event " << i << " / " << nEntries << "\r" << std::flush;
        
        // Load data for this event
        reader.GetEntry(i);

        // -------------------------------------------------------------
        // STEP 0: Index Finding (Required First)
        // -------------------------------------------------------------
        // We use vector<bool> for O(1) instant access instead of std::set
        std::vector<bool> ls_is_used(reader.ls_pt->size(), false);
        std::vector<bool> pls_is_used(reader.pls_pt->size(), false);
        
        PhysicsAlgo::findUsedIndices(ls_is_used, pls_is_used, reader);

        // Data structures for stats and physics mapping
        std::map<int, PlsMatchStatus> simIdx_to_pls_status;
        std::map<int, LsMatchStatus> simIdx_to_ls_status;
        std::map<int, std::vector<int>> ls_simIdx_map;

        // -------------------------------------------------------------
        // LOOP 1: PLS Pre-Calculation
        // (We MUST do this first to get stats for LS Detailed Plots)
        // This replaces 'getPlsMatchStatus'
        // -------------------------------------------------------------
        for (size_t j = 0; j < reader.pls_pt->size(); ++j) {
            if (reader.pls_isFake->at(j) == 0) {
                int sim = reader.pls_simIdx->at(j);
                if (sim != -999) {
                    if (pls_is_used[j]) simIdx_to_pls_status[sim].used_count++;
                    else simIdx_to_pls_status[sim].unused_count++;
                }
            }
        }

        // -------------------------------------------------------------
        // LOOP 2: LS "Mega-Loop"
        // (Calc Stats + Build Map + Fill LS Plots)
        // This replaces 'getLsMatchStatus' AND the map building loop
        // -------------------------------------------------------------
        for (size_t j = 0; j < reader.ls_pt->size(); ++j) {
            bool is_used = ls_is_used[j];
            bool is_fake = (reader.ls_isFake->at(j) != 0);
            int sim = reader.ls_simIdx->at(j);

            // A. Calc Stats & Build Map
            if (!is_fake && sim != -999) {
                // Update Status
                if (is_used) simIdx_to_ls_status[sim].used_count++;
                else simIdx_to_ls_status[sim].unused_count++;
                
                // Build Map for pT2 logic (Optimization!)
                ls_simIdx_map[sim].push_back(j);
            }

            // B. Fill LS Histograms (Uses PLS Stats from Loop 1)
            hists.FillSimpleHistograms(
                is_used, is_fake, 
                reader.ls_pt->at(j), reader.ls_eta->at(j), reader.ls_phi->at(j),
                hists.h_pt_real_ls_used,   hists.h_pt_fake_ls_used,
                hists.h_pt_real_ls_unused, hists.h_pt_fake_ls_unused,
                hists.h_eta_real_ls_used,  hists.h_eta_fake_ls_used,
                hists.h_eta_real_ls_unused,hists.h_eta_fake_ls_unused,
                hists.h_phi_real_ls_used,  hists.h_phi_fake_ls_used,
                hists.h_phi_real_ls_unused,hists.h_phi_fake_ls_unused
            );

            hists.FillDetailedLsHistograms(
                is_used, is_fake,
                reader.ls_pt->at(j), reader.ls_eta->at(j), reader.ls_phi->at(j),
                sim, simIdx_to_pls_status
            );
        }

        // -------------------------------------------------------------
        // LOOP 3: PLS "Physics Loop"
        // (Fill PLS Plots + Process pT2 Physics)
        // -------------------------------------------------------------
        for (size_t j = 0; j < reader.pls_pt->size(); ++j) {
            bool is_used = pls_is_used[j];
            bool is_fake = (reader.pls_isFake->at(j) != 0);

            // A. Fill PLS Histograms (Uses LS Stats from Loop 2)
            hists.FillSimpleHistograms(
                is_used, is_fake, 
                reader.pls_pt->at(j), reader.pls_eta->at(j), reader.pls_phi->at(j),
                hists.h_pt_real_pls_used,   hists.h_pt_fake_pls_used,
                hists.h_pt_real_pls_unused, hists.h_pt_fake_pls_unused,
                hists.h_eta_real_pls_used,  hists.h_eta_fake_pls_used,
                hists.h_eta_real_pls_unused,hists.h_eta_fake_pls_unused,
                hists.h_phi_real_pls_used,  hists.h_phi_fake_pls_used,
                hists.h_phi_real_pls_unused,hists.h_phi_fake_pls_unused
            );

            hists.FillDetailedPlsHistograms(
                is_used, is_fake,
                reader.pls_pt->at(j), reader.pls_eta->at(j), reader.pls_phi->at(j),
                reader.pls_simIdx->at(j), simIdx_to_ls_status
            );

            // B. Process pT2 Physics
            // We only need to check physics if it's NOT used in a track
            if (!is_used) {
                // We call the Single-Index version of runPT2Physics here
                PhysicsAlgo::runPT2Physics(j, reader, ls_is_used, ls_simIdx_map, hists, ideal_pt2_count);
            }
        }
    }
    
    std::cout << "\nFinished processing entries." << std::endl;
    std::cout << "Total ideal pt2 objects found: " << ideal_pt2_count << std::endl;

    // -------------------------------------------------------
    // 3. Plotting & Saving
    // -------------------------------------------------------
    gStyle->SetOptStat(0);
    std::cout << "Generating plots..." << std::endl;

    // --- Simple LS Plots ---
    Plotter::createAndSaveSimpleStack("c_pt_ls_used", "pT of Used LS;p_{T} [GeV];LS Objects", "compare_pt_ls_used.png", "Used LS", hists.h_pt_real_ls_used, hists.h_pt_fake_ls_used, "Real, Used", "Fake, Used", kGreen+2, kOrange-3, true);
    Plotter::createAndSaveSimpleStack("c_eta_ls_used", "Eta of Used LS;#eta;LS Objects", "compare_eta_ls_used.png", "Used LS", hists.h_eta_real_ls_used, hists.h_eta_fake_ls_used, "Real, Used", "Fake, Used", kGreen+2, kOrange-3);
    Plotter::createAndSaveSimpleStack("c_phi_ls_used", "Phi of Used LS;#phi;LS Objects", "compare_phi_ls_used.png", "Used LS", hists.h_phi_real_ls_used, hists.h_phi_fake_ls_used, "Real, Used", "Fake, Used", kGreen+2, kOrange-3);

    Plotter::createAndSaveSimpleStack("c_pt_ls_unused", "pT of Unused LS;p_{T} [GeV];LS Objects", "compare_pt_ls_unused.png", "Unused LS", hists.h_pt_real_ls_unused, hists.h_pt_fake_ls_unused, "Real, Unused", "Fake, Unused", kGreen-7, kOrange-9, true);
    Plotter::createAndSaveSimpleStack("c_eta_ls_unused", "Eta of Unused LS;#eta;LS Objects", "compare_eta_ls_unused.png", "Unused LS", hists.h_eta_real_ls_unused, hists.h_eta_fake_ls_unused, "Real, Unused", "Fake, Unused", kGreen-7, kOrange-9);
    Plotter::createAndSaveSimpleStack("c_phi_ls_unused", "Phi of Unused LS;#phi;LS Objects", "compare_phi_ls_unused.png", "Unused LS", hists.h_phi_real_ls_unused, hists.h_phi_fake_ls_unused, "Real, Unused", "Fake, Unused", kGreen-7, kOrange-9);

    // --- Simple PLS Plots ---
    Plotter::createAndSaveSimpleStack("c_pt_pls_used", "pT of Used PLS;p_{T} [GeV];PLS Objects", "compare_pt_pls_used.png", "Used PLS", hists.h_pt_real_pls_used, hists.h_pt_fake_pls_used, "Real, Used", "Fake, Used", kAzure+1, kRed-7, true);
    Plotter::createAndSaveSimpleStack("c_eta_pls_used", "Eta of Used PLS;#eta;PLS Objects", "compare_eta_pls_used.png", "Used PLS", hists.h_eta_real_pls_used, hists.h_eta_fake_pls_used, "Real, Used", "Fake, Used", kAzure+1, kRed-7);
    Plotter::createAndSaveSimpleStack("c_phi_pls_used", "Phi of Used PLS;#phi;PLS Objects", "compare_phi_pls_used.png", "Used PLS", hists.h_phi_real_pls_used, hists.h_phi_fake_pls_used, "Real, Used", "Fake, Used", kAzure+1, kRed-7);

    Plotter::createAndSaveSimpleStack("c_pt_pls_unused", "pT of Unused PLS;p_{T} [GeV];PLS Objects", "compare_pt_pls_unused.png", "Unused PLS", hists.h_pt_real_pls_unused, hists.h_pt_fake_pls_unused, "Real, Unused", "Fake, Unused", kAzure-9, kRed-9, true);
    Plotter::createAndSaveSimpleStack("c_eta_pls_unused", "Eta of Unused PLS;#eta;PLS Objects", "compare_eta_pls_unused.png", "Unused PLS", hists.h_eta_real_pls_unused, hists.h_eta_fake_pls_unused, "Real, Unused", "Fake, Unused", kAzure-9, kRed-9);
    Plotter::createAndSaveSimpleStack("c_phi_pls_unused", "Phi of Unused PLS;#phi;PLS Objects", "compare_phi_pls_unused.png", "Unused PLS", hists.h_phi_real_pls_unused, hists.h_phi_fake_pls_unused, "Real, Unused", "Fake, Unused", kAzure-9, kRed-9);

    // --- Detailed LS Plots ---
    Plotter::createAndSaveDetailedStack("c_pt_ls_used_detailed", "pT of Used LS (Detailed);p_{T} [GeV];LS Objects", "compare_pt_ls_used_detailed.png", "Used Real LS Matched To PLS:", true, hists.h_pt_fake_ls_used_detailed, hists.h_pt_real_ls_used_unmatched, hists.h_pt_real_ls_used_matchedToUnusedOnly, hists.h_pt_real_ls_used_matchedToBoth, hists.h_pt_real_ls_used_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_eta_ls_used_detailed", "Eta of Used LS (Detailed);#eta;LS Objects", "compare_eta_ls_used_detailed.png", "Used Real LS Matched To PLS:", false, hists.h_eta_fake_ls_used_detailed, hists.h_eta_real_ls_used_unmatched, hists.h_eta_real_ls_used_matchedToUnusedOnly, hists.h_eta_real_ls_used_matchedToBoth, hists.h_eta_real_ls_used_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_phi_ls_used_detailed", "Phi of Used LS (Detailed);#phi;LS Objects", "compare_phi_ls_used_detailed.png", "Used Real LS Matched To PLS:", false, hists.h_phi_fake_ls_used_detailed, hists.h_phi_real_ls_used_unmatched, hists.h_phi_real_ls_used_matchedToUnusedOnly, hists.h_phi_real_ls_used_matchedToBoth, hists.h_phi_real_ls_used_matchedToUsedOnly);
    
    Plotter::createAndSaveDetailedStack("c_pt_ls_unused_detailed", "pT of Unused LS (Detailed);p_{T} [GeV];LS Objects", "compare_pt_ls_unused_detailed.png", "Unused Real LS Matched To PLS:", true, hists.h_pt_fake_ls_unused_detailed, hists.h_pt_real_ls_unused_unmatched, hists.h_pt_real_ls_unused_matchedToUnusedOnly, hists.h_pt_real_ls_unused_matchedToBoth, hists.h_pt_real_ls_unused_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_eta_ls_unused_detailed", "Eta of Unused LS (Detailed);#eta;LS Objects", "compare_eta_ls_unused_detailed.png", "Unused Real LS Matched To PLS:", false, hists.h_eta_fake_ls_unused_detailed, hists.h_eta_real_ls_unused_unmatched, hists.h_eta_real_ls_unused_matchedToUnusedOnly, hists.h_eta_real_ls_unused_matchedToBoth, hists.h_eta_real_ls_unused_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_phi_ls_unused_detailed", "Phi of Unused LS (Detailed);#phi;LS Objects", "compare_phi_ls_unused_detailed.png", "Unused Real LS Matched To PLS:", false, hists.h_phi_fake_ls_unused_detailed, hists.h_phi_real_ls_unused_unmatched, hists.h_phi_real_ls_unused_matchedToUnusedOnly, hists.h_phi_real_ls_unused_matchedToBoth, hists.h_phi_real_ls_unused_matchedToUsedOnly);

    // --- Detailed PLS Plots ---
    Plotter::createAndSaveDetailedStack("c_pt_pls_used_detailed", "pT of Used PLS (Detailed);p_{T} [GeV];PLS Objects", "compare_pt_pls_used_detailed.png", "Used Real PLS Matched To LS:", true, hists.h_pt_fake_pls_used_detailed, hists.h_pt_real_pls_used_unmatched, hists.h_pt_real_pls_used_matchedToUnusedOnly, hists.h_pt_real_pls_used_matchedToBoth, hists.h_pt_real_pls_used_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_eta_pls_used_detailed", "Eta of Used PLS (Detailed);#eta;PLS Objects", "compare_eta_pls_used_detailed.png", "Used Real PLS Matched To LS:", false, hists.h_eta_fake_pls_used_detailed, hists.h_eta_real_pls_used_unmatched, hists.h_eta_real_pls_used_matchedToUnusedOnly, hists.h_eta_real_pls_used_matchedToBoth, hists.h_eta_real_pls_used_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_phi_pls_used_detailed", "Phi of Used PLS (Detailed);#phi;PLS Objects", "compare_phi_pls_used_detailed.png", "Used Real PLS Matched To LS:", false, hists.h_phi_fake_pls_used_detailed, hists.h_phi_real_pls_used_unmatched, hists.h_phi_real_pls_used_matchedToUnusedOnly, hists.h_phi_real_pls_used_matchedToBoth, hists.h_phi_real_pls_used_matchedToUsedOnly);
    
    Plotter::createAndSaveDetailedStack("c_pt_pls_unused_detailed", "pT of Unused PLS (Detailed);p_{T} [GeV];PLS Objects", "compare_pt_pls_unused_detailed.png", "Unused Real PLS Matched To LS:", true, hists.h_pt_fake_pls_unused_detailed, hists.h_pt_real_pls_unused_unmatched, hists.h_pt_real_pls_unused_matchedToUnusedOnly, hists.h_pt_real_pls_unused_matchedToBoth, hists.h_pt_real_pls_unused_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_eta_pls_unused_detailed", "Eta of Unused PLS (Detailed);#eta;PLS Objects", "compare_eta_pls_unused_detailed.png", "Unused Real PLS Matched To LS:", false, hists.h_eta_fake_pls_unused_detailed, hists.h_eta_real_pls_unused_unmatched, hists.h_eta_real_pls_unused_matchedToUnusedOnly, hists.h_eta_real_pls_unused_matchedToBoth, hists.h_eta_real_pls_unused_matchedToUsedOnly);
    Plotter::createAndSaveDetailedStack("c_phi_pls_unused_detailed", "Phi of Unused PLS (Detailed);#phi;PLS Objects", "compare_phi_pls_unused_detailed.png", "Unused Real PLS Matched To LS:", false, hists.h_phi_fake_pls_unused_detailed, hists.h_phi_real_pls_unused_unmatched, hists.h_phi_real_pls_unused_matchedToUnusedOnly, hists.h_phi_real_pls_unused_matchedToBoth, hists.h_phi_real_pls_unused_matchedToUsedOnly);

    // --- Delta Plots ---
    Plotter::createAndSaveSimplePlot("c_delta_pt", "delta_pt.png", hists.h_delta_pt, "#Delta p_{T} (LS - PLS)", "#Delta p_{T} [GeV]", "pT2 Objects");
    Plotter::createAndSaveSimplePlot("c_delta_eta", "delta_eta.png", hists.h_delta_eta, "#Delta #eta (LS - PLS)", "#Delta #eta", "pT2 Objects");
    Plotter::createAndSaveSimplePlot("c_delta_phi", "delta_phi.png", hists.h_delta_phi, "#Delta #phi (LS - PLS)", "#Delta #phi [rad]", "pT2 Objects");
    Plotter::createAndSaveSimplePlot("c_delta_R", "delta_R.png", hists.h_delta_R, "#Delta R (LS - PLS)", "#Delta R", "pT2 Objects");

    Plotter::createAndSaveSimplePlot("c_delta_pt_cand", "delta_pt_cand.png", hists.h_delta_pt_cand, "#Delta p_{T} (LS - PLS)", "#Delta p_{T} [GeV]", "pT2 Objects");
    Plotter::createAndSaveSimplePlot("c_delta_eta_cand", "delta_eta_cand.png", hists.h_delta_eta_cand, "#Delta #eta (LS - PLS)", "#Delta #eta", "pT2 Objects");
    Plotter::createAndSaveSimplePlot("c_delta_phi_cand", "delta_phi_cand.png", hists.h_delta_phi_cand, "#Delta #phi (LS - PLS)", "#Delta #phi [rad]", "pT2 Objects");
    Plotter::createAndSaveSimplePlot("c_delta_R_cand", "delta_R_cand.png", hists.h_delta_R_cand, "#Delta R (LS - PLS)", "#Delta R", "pT2 Objects");

    // --- Extrapolation Plots ---
    Plotter::createAndSaveSimplePlot("c_extrapolation_dist_3d", "extrapolation_dist_3d.png", hists.h_extrapolation_dist_3d,"3D Distance between extrapolated PLS and LS", "Distance [cm]", "Ideal pT2 Pairs");
    Plotter::createAndSaveSimplePlot("c_extrapolation_delta_z", "extrapolation_delta_z.png", hists.h_extrapolation_delta_z,"R-Z Extrapolation #Delta r (LS - PLS)", "#Delta r [cm]", "Ideal pT2 Pairs");
    Plotter::createAndSaveSimplePlot("c_delta_r_reverse_combined", "extrapolation_delta_r_reverse_combined.png", hists.h_extrapolation_delta_r_reverse_combined, "Inclusive R-Z Extrapolation #Delta r (LS -> All PLS Hits)", "#Delta r [cm]", "Entries");

    // --- Layer/Hit Specific Plots ---
    std::cout << "Saving Layer/Hit plots..." << std::endl;
    
    // Reverse Hits
    for (size_t i = 0; i < hists.h_extrapolation_delta_r_reverse_by_hit.size(); ++i) {
        if (hists.h_extrapolation_delta_r_reverse_by_hit[i]->GetEntries() > 0) {
            std::string canvas_name = "c_delta_r_reverse_hit_" + std::to_string(i);
            std::string file_name = "extrapolation_delta_r_reverse_hit_" + std::to_string(i) + ".png";
            std::string plot_title = "R-Z Extrapolation #Delta r (LS -> PLS Hit " + std::to_string(i) + ")";
            Plotter::createAndSaveSimplePlot(canvas_name, file_name, hists.h_extrapolation_delta_r_reverse_by_hit[i], plot_title, "#Delta r [cm]", "Entries");
        }
    }

    // Delta R by Layer
    for (size_t i = 0; i < hists.h_delta_r_by_layer.size(); ++i) {
        if (hists.h_delta_r_by_layer[i]->GetEntries() > 0) {
            std::string canvas_name = "c_delta_r_layer_" + std::to_string(i);
            std::string file_name = "extrapolation_delta_r_layer_" + std::to_string(i) + ".png";
            std::string plot_title = "Z-R Extrapolation #Delta r (Layer " + std::to_string(i) + ")";
            Plotter::createAndSaveSimplePlot(canvas_name, file_name, hists.h_delta_r_by_layer[i], plot_title, "#Delta r [cm]", "Entries");
        }
    }

    // 3D Dist by Layer
    for (size_t i = 0; i < hists.h_dist_3d_by_layer.size(); ++i) {
        if (hists.h_dist_3d_by_layer[i]->GetEntries() > 0) {
            std::string canvas_name = "c_dist_3d_layer_" + std::to_string(i);
            std::string file_name = "extrapolation_dist_3d_layer_" + std::to_string(i) + ".png";
            std::string plot_title = "3D Distance (Layer " + std::to_string(i) + ")";
            Plotter::createAndSaveSimplePlot(canvas_name, file_name, hists.h_dist_3d_by_layer[i], plot_title, "Distance [cm]", "Entries");
        }
    }

    std::cout << "Done." << std::endl;
    return 0;
}
