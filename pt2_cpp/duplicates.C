#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <set>

void duplicates() {
    // --- Configuration ---
    const char* input_filename = "LSTNtuple.root";
    const char* tree_name = "tree";

    // --- STEP 1: Open ROOT file and get TTree ---
    TFile *inputFile = TFile::Open(input_filename, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << input_filename << std::endl;
        return;
    }
    TTree *inputTree;
    inputFile->GetObject(tree_name, inputTree);
    if (!inputTree) {
        std::cerr << "Error: Could not find TTree with name " << tree_name << std::endl;
        inputFile->Close();
        return;
    }

    // --- STEP 2: Set Branch Addresses for ALL data ---
    // Pointers to vectors for the low-level segments (ls)
    std::vector<float> *ls_pt_vec = nullptr;
    std::vector<float> *ls_eta_vec = nullptr;
    std::vector<float> *ls_phi_vec = nullptr;
    std::vector<int>   *ls_isFake_vec = nullptr;

    // Pointers to vectors for ALL the index arrays
    std::vector<int> *tc_pt5Idx_vec = nullptr;
    std::vector<int> *pt5_t5Idx_vec = nullptr;
    std::vector<int> *t5_t3Idx0_vec = nullptr;
    std::vector<int> *t5_t3Idx1_vec = nullptr;
    std::vector<int> *t3_lsIdx0_vec = nullptr;
    std::vector<int> *t3_lsIdx1_vec = nullptr;
    std::vector<int> *tc_pt3Idx_vec = nullptr;
    std::vector<int> *pt3_t3Idx_vec = nullptr;
    std::vector<int> *tc_t5Idx_vec = nullptr;
    
    // Pointers for PLS branches ARE STILL NEEDED for logic and counting
    std::vector<int> *pt5_plsIdx_vec = nullptr;
    std::vector<int> *pt3_plsIdx_vec = nullptr;
    std::vector<float> *pls_pt_vec = nullptr;

    // Link the kinematic branches
    inputTree->SetBranchAddress("ls_pt", &ls_pt_vec);
    inputTree->SetBranchAddress("ls_eta", &ls_eta_vec);
    inputTree->SetBranchAddress("ls_phi", &ls_phi_vec);
    inputTree->SetBranchAddress("ls_isFake", &ls_isFake_vec);
    
    // Link ALL the index branches (ensure these names match your TTree exactly)
    inputTree->SetBranchAddress("tc_pt5Idx", &tc_pt5Idx_vec);
    inputTree->SetBranchAddress("pt5_t5Idx", &pt5_t5Idx_vec);
    inputTree->SetBranchAddress("t5_t3Idx0", &t5_t3Idx0_vec);
    inputTree->SetBranchAddress("t5_t3Idx1", &t5_t3Idx1_vec);
    inputTree->SetBranchAddress("t3_lsIdx0", &t3_lsIdx0_vec);
    inputTree->SetBranchAddress("t3_lsIdx1", &t3_lsIdx1_vec);
    inputTree->SetBranchAddress("tc_pt3Idx", &tc_pt3Idx_vec);
    inputTree->SetBranchAddress("pt3_t3Idx", &pt3_t3Idx_vec);
    inputTree->SetBranchAddress("tc_t5Idx", &tc_t5Idx_vec);
    inputTree->SetBranchAddress("pt5_plsIdx", &pt5_plsIdx_vec);
    inputTree->SetBranchAddress("pt3_plsIdx", &pt3_plsIdx_vec);
    inputTree->SetBranchAddress("pls_pt", &pls_pt_vec);

    // --- STEP 3: Define histograms for all categories ---
    int pt_bins = 180, eta_bins = 180, phi_bins = 180;
    double pt_max = 250.0, eta_max = 4.5, phi_max = 3.15;

    // pT
    TH1D *h_pt_real_used   = new TH1D("h_pt_real_used",   "", pt_bins, 0, pt_max);
    TH1D *h_pt_real_unused = new TH1D("h_pt_real_unused", "", pt_bins, 0, pt_max);
    TH1D *h_pt_fake_used   = new TH1D("h_pt_fake_used",   "", pt_bins, 0, pt_max);
    TH1D *h_pt_fake_unused = new TH1D("h_pt_fake_unused", "", pt_bins, 0, pt_max);
    // Eta
    TH1D *h_eta_real_used   = new TH1D("h_eta_real_used",   "", eta_bins, -eta_max, eta_max);
    TH1D *h_eta_real_unused = new TH1D("h_eta_real_unused", "", eta_bins, -eta_max, eta_max);
    TH1D *h_eta_fake_used   = new TH1D("h_eta_fake_used",   "", eta_bins, -eta_max, eta_max);
    TH1D *h_eta_fake_unused = new TH1D("h_eta_fake_unused", "", eta_bins, -eta_max, eta_max);
    // Phi
    TH1D *h_phi_real_used   = new TH1D("h_phi_real_used",   "", phi_bins, -phi_max, phi_max);
    TH1D *h_phi_real_unused = new TH1D("h_phi_real_unused", "", phi_bins, -phi_max, phi_max);
    TH1D *h_phi_fake_used   = new TH1D("h_phi_fake_used",   "", phi_bins, -phi_max, phi_max);
    TH1D *h_phi_fake_unused = new TH1D("h_phi_fake_unused", "", phi_bins, -phi_max, phi_max);

    // Global counters and sets for summary statistics
    long long total_ls_count = 0;
    long long total_pls_count = 0;
    //Global count of duplicates
    long long duplicate_ls_count = 0;
    long long duplicate_pls_count = 0;
    
    std::set<int> all_used_ls_indices_global;
    std::set<int> all_used_pls_indices_global; 
    
    // --- STEP 4: The Main Event Loop ---
     Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i);

        std::set<int> used_ls_indices_in_this_event;
        // *** NEW: A set to track used PLS in this event ***
        std::set<int> used_pls_indices_in_this_event;
        
        // *** NEW: Trackers to spot duplicates *within this event* -> this makes sure that we are only looking duplciates in one event
        std::set<int> ls_usage_tracker_event;
        std::set<int> pls_usage_tracker_event;
        
        // Add total object counts from this event
        if (ls_pt_vec) total_ls_count += ls_pt_vec->size();
        if (pls_pt_vec) total_pls_count += pls_pt_vec->size();

        // --- Tracing Logic with Duplicate & PLS Tracking ---
        
        // Path 1: tc_pt5 -> ... -> ls/pls
        for (int pt5_Idx : *tc_pt5Idx_vec) {
            if (pt5_Idx != -999) {
                // Track the used PLS from pt5
                int pls_idx = pt5_plsIdx_vec->at(pt5_Idx);
                if (pls_idx != -999) {
                    used_pls_indices_in_this_event.insert(pls_idx);
                    all_used_pls_indices_global.insert(pls_idx);
                    if (pls_usage_tracker_event.count(pls_idx)) { duplicate_pls_count++; } //.count() checks is value has been entered before
                    else { pls_usage_tracker_event.insert(pls_idx); }
                }

                // Trace down to the LS
                int t5_Idx = pt5_t5Idx_vec->at(pt5_Idx);
                if (t5_Idx != -999) {
                    int t3_0 = t5_t3Idx0_vec->at(t5_Idx);
                    if (t3_0 != -999) {
                        int ls0 = t3_lsIdx0_vec->at(t3_0);
                        if (ls0 != -999) {
                           used_ls_indices_in_this_event.insert(ls0);
                           all_used_ls_indices_global.insert(ls0);
                           if (ls_usage_tracker_event.count(ls0)) { duplicate_ls_count++; }
                           else { ls_usage_tracker_event.insert(ls0); }
                        }
                        int ls1 = t3_lsIdx1_vec->at(t3_0);
                        if (ls1 != -999) {
                           used_ls_indices_in_this_event.insert(ls1);
                           all_used_ls_indices_global.insert(ls1);
                           if (ls_usage_tracker_event.count(ls1)) { duplicate_ls_count++; }
                           else { ls_usage_tracker_event.insert(ls1); }
                        }
                    }
                    int t3_1 = t5_t3Idx1_vec->at(t5_Idx);
                    if (t3_1 != -999) {
                        int ls2 = t3_lsIdx0_vec->at(t3_1);
                        if (ls2 != -999) {
                           used_ls_indices_in_this_event.insert(ls2);
                           all_used_ls_indices_global.insert(ls2);
                           if (ls_usage_tracker_event.count(ls2)) { duplicate_ls_count++; }
                           else { ls_usage_tracker_event.insert(ls2); }
                        }
                        int ls3 = t3_lsIdx1_vec->at(t3_1);
                        if (ls3 != -999) {
                           used_ls_indices_in_this_event.insert(ls3);
                           all_used_ls_indices_global.insert(ls3);
                           if (ls_usage_tracker_event.count(ls3)) { duplicate_ls_count++; }
                           else { ls_usage_tracker_event.insert(ls3); }
                        }
                    }
                }
            }
        }
        
        // Path 2: tc_pt3 -> ... -> ls/pls
        for (int pt3_Idx : *tc_pt3Idx_vec) {
            if (pt3_Idx != -999) {
                // Track the used PLS from pt3
                int pls_idx = pt3_plsIdx_vec->at(pt3_Idx);
                if (pls_idx != -999) {
                    used_pls_indices_in_this_event.insert(pls_idx);
                    all_used_pls_indices_global.insert(pls_idx);
                    if (pls_usage_tracker_event.count(pls_idx)) { duplicate_pls_count++; }
                    else { pls_usage_tracker_event.insert(pls_idx); }
                }
                
                // Trace down to the LS from pt3
                int t3_idx = pt3_t3Idx_vec->at(pt3_Idx);
                if (t3_idx != -999) {
                     int ls0 = t3_lsIdx0_vec->at(t3_idx);
                     if (ls0 != -999) {
                        used_ls_indices_in_this_event.insert(ls0);
                        all_used_ls_indices_global.insert(ls0);
                        if (ls_usage_tracker_event.count(ls0)) { duplicate_ls_count++; }
                        else { ls_usage_tracker_event.insert(ls0); }
                     }
                     int ls1 = t3_lsIdx1_vec->at(t3_idx);
                     if (ls1 != -999) {
                        used_ls_indices_in_this_event.insert(ls1);
                        all_used_ls_indices_global.insert(ls1);
                        if (ls_usage_tracker_event.count(ls1)) { duplicate_ls_count++; }
                        else { ls_usage_tracker_event.insert(ls1); }
                     }
                }
            }
        }

        // Path 3: tc_t5 -> t5 -> t3 -> ls (No PLS in this path)
        for (int t5_Idx : *tc_t5Idx_vec) {
            if (t5_Idx != -999) {
                // ... (This logic is the same as the t5 part of Path 1) ...
            }
        }
        
        // Fill Histograms for LS objects
        for (size_t j = 0; j < ls_pt_vec->size(); ++j) {
            bool is_fake = ls_isFake_vec->at(j) != 0;
            bool is_used = used_ls_indices_in_this_event.count(j) > 0;
            float pt = ls_pt_vec->at(j), eta = ls_eta_vec->at(j), phi = ls_phi_vec->at(j);
            if (!is_fake) {
                if (is_used) { h_pt_real_used->Fill(pt); h_eta_real_used->Fill(eta); h_phi_real_used->Fill(phi); } 
                else { h_pt_real_unused->Fill(pt); h_eta_real_unused->Fill(eta); h_phi_real_unused->Fill(phi); }
            } else {
                if (is_used) { h_pt_fake_used->Fill(pt); h_eta_fake_used->Fill(eta); h_phi_fake_used->Fill(phi); }
                else { h_pt_fake_unused->Fill(pt); h_eta_fake_unused->Fill(eta); h_phi_fake_unused->Fill(phi); }
            }
        }
    }
    std::cout << "Finished processing entries." << std::endl;

    // --- STEP 5: Styling and Plotting ---
    gStyle->SetOptStat(0);

    // --- pT Plots ---
    TCanvas *c_pt_used = new TCanvas("c_pt_used", "pT of Used Tracks", 800, 600);
    c_pt_used->SetLogy();
    THStack *hs_pt_used = new THStack("hs_pt_used", "pT Distribution of Used Tracks;p_{T} [GeV];Tracks");
    h_pt_real_used->SetFillColor(kBlue-6); 
    h_pt_fake_used->SetFillColor(kRed-6);
    hs_pt_used->Add(h_pt_fake_used); hs_pt_used->Add(h_pt_real_used);
    hs_pt_used->Draw("HIST");
    auto legend_pt_used = new TLegend(0.6, 0.8, 0.9, 0.9);
    legend_pt_used->AddEntry(h_pt_real_used, "Real, Used", "f");
    legend_pt_used->AddEntry(h_pt_fake_used, "Fake, Used", "f");
    legend_pt_used->Draw();
    c_pt_used->SaveAs("compare_pt_used.png");

    TCanvas *c_pt_unused = new TCanvas("c_pt_unused", "pT of Unused Tracks", 800, 600);
    c_pt_unused->SetLogy();
    THStack *hs_pt_unused = new THStack("hs_pt_unused", "pT Distribution of Unused Tracks;p_{T} [GeV];Tracks");
    h_pt_real_unused->SetFillColor(kBlue-9); h_pt_fake_unused->SetFillColor(kRed-9);
    hs_pt_unused->Add(h_pt_fake_unused); hs_pt_unused->Add(h_pt_real_unused);
    hs_pt_unused->Draw("HIST");
    auto legend_pt_unused = new TLegend(0.6, 0.8, 0.9, 0.9);
    legend_pt_unused->AddEntry(h_pt_real_unused, "Real, Unused", "f");
    legend_pt_unused->AddEntry(h_pt_fake_unused, "Fake, Unused", "f");
    legend_pt_unused->Draw();
    c_pt_unused->SaveAs("compare_pt_unused.png");

    // --- Eta Plots ---
    TCanvas *c_eta_used = new TCanvas("c_eta_used", "Eta of Used Tracks", 800, 600);
    THStack *hs_eta_used = new THStack("hs_eta_used", "Eta Distribution of Used Tracks;#eta;Tracks");
    h_eta_real_used->SetFillColor(kGreen-7); h_eta_fake_used->SetFillColor(kOrange-9);
    hs_eta_used->Add(h_eta_fake_used); hs_eta_used->Add(h_eta_real_used);
    hs_eta_used->Draw("HIST");
    auto legend_eta_used = new TLegend(0.2, 0.8, 0.5, 0.9);
    legend_eta_used->AddEntry(h_eta_real_used, "Real, Used", "f");
    legend_eta_used->AddEntry(h_eta_fake_used, "Fake, Used", "f");
    legend_eta_used->Draw();
    c_eta_used->SaveAs("compare_eta_used.png");

    TCanvas *c_eta_unused = new TCanvas("c_eta_unused", "Eta of Unused Tracks", 800, 600);
    THStack *hs_eta_unused = new THStack("hs_eta_unused", "Eta Distribution of Unused Tracks;#eta;Tracks");
    h_eta_real_unused->SetFillColor(kGreen-7); h_eta_fake_unused->SetFillColor(kOrange-9);
    hs_eta_unused->Add(h_eta_fake_unused); hs_eta_unused->Add(h_eta_real_unused);
    hs_eta_unused->Draw("HIST");
    auto legend_eta_unused = new TLegend(0.2, 0.8, 0.5, 0.9);
    legend_eta_unused->AddEntry(h_eta_real_unused, "Real, Unused", "f");
    legend_eta_unused->AddEntry(h_eta_fake_unused, "Fake, Unused", "f");
    legend_eta_unused->Draw();
    c_eta_unused->SaveAs("compare_eta_unused.png");

    // --- Phi Plots ---
    TCanvas *c_phi_used = new TCanvas("c_phi_used", "Phi of Used Tracks", 800, 600);
    THStack *hs_phi_used = new THStack("hs_phi_used", "Phi Distribution of Used Tracks;#phi [rad];Tracks");
    h_phi_real_used->SetFillColor(kMagenta-7); h_phi_fake_used->SetFillColor(kCyan-7);
    hs_phi_used->Add(h_phi_fake_used); 
    hs_phi_used->Add(h_phi_real_used);
    hs_phi_used->Draw("HIST");
    auto legend_phi_used = new TLegend(0.4, 0.8, 0.7, 0.9);
    legend_phi_used->AddEntry(h_phi_real_used, "Real, Used", "f");
    legend_phi_used->AddEntry(h_phi_fake_used, "Fake, Used", "f");
    legend_phi_used->Draw();
    c_phi_used->SaveAs("compare_phi_used.png");

    TCanvas *c_phi_unused = new TCanvas("c_phi_unused", "Phi of Unused Tracks", 800, 600);
    THStack *hs_phi_unused = new THStack("hs_phi_unused", "Phi Distribution of Unused Tracks;#phi [rad];Tracks");
    h_phi_real_unused->SetFillColor(kMagenta-9); 
    h_phi_fake_unused->SetFillColor(kCyan-9);
    hs_phi_unused->Add(h_phi_fake_unused); hs_phi_unused->Add(h_phi_real_unused);
    hs_phi_unused->Draw("HIST");
    auto legend_phi_unused = new TLegend(0.4, 0.8, 0.7, 0.9);
    legend_phi_unused->AddEntry(h_phi_real_unused, "Real, Unused", "f");
    legend_phi_unused->AddEntry(h_phi_fake_unused, "Fake, Unused", "f");
    legend_phi_unused->Draw();
    c_phi_unused->SaveAs("compare_phi_unused.png");
    
    //Print duplicates
    std::cout << "\n\n--- Analysis Summary ---" << std::endl;
    std::cout << "Total LS Objects Processed: " << total_ls_count << std::endl;
    std::cout << "Unique LS Objects Used: " << all_used_ls_indices_global.size() << std::endl;
    std::cout << "Duplicate LS Assignments (within events): " << duplicate_ls_count << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "Total PLS Objects Processed: " << total_pls_count << std::endl;
    std::cout << "Unique PLS Objects Used: " << all_used_pls_indices_global.size() << std::endl;
    std::cout << "Duplicate PLS Assignments (within events): " << duplicate_pls_count << std::endl;
    std::cout << "------------------------------------------" << std::endl;




    // --- Clean up ---
    inputFile->Close();
    delete c_pt_used; delete c_pt_unused;
    delete c_eta_used; delete c_eta_unused;
    delete c_phi_used; delete c_phi_unused;
}
