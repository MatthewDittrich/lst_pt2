// 07/12/25

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

void used_unused() {
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

    // --- STEP 4: The Main Event Loop ---
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i);

        std::set<int> used_ls_indices_in_this_event;

        // Path 1: tc_pt5 -> pt5 -> t5 -> t3 -> ls
        for (int pt5_Idx : *tc_pt5Idx_vec) {
            if (pt5_Idx != -999 && pt5_Idx < pt5_t5Idx_vec->size()) {
                int t5_Idx = pt5_t5Idx_vec->at(pt5_Idx);
                if (t5_Idx != -999 && t5_Idx < t5_t3Idx0_vec->size()) {
                    int t3_0 = t5_t3Idx0_vec->at(t5_Idx);
                    if (t3_0 != -999 && t3_0 < t3_lsIdx0_vec->size()) {
                        if (t3_lsIdx0_vec->at(t3_0) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx0_vec->at(t3_0));
                        if (t3_lsIdx1_vec->at(t3_0) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx1_vec->at(t3_0));
                    }
                    int t3_1 = t5_t3Idx1_vec->at(t5_Idx);
                     if (t3_1 != -999 && t3_1 < t3_lsIdx0_vec->size()) {
                        if (t3_lsIdx0_vec->at(t3_1) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx0_vec->at(t3_1));
                        if (t3_lsIdx1_vec->at(t3_1) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx1_vec->at(t3_1));
                    }
                }
            }
        }
        
        // Path 2: tc_pt3 -> pt3 -> t3 -> ls
        for (int pt3_Idx : *tc_pt3Idx_vec) {
            if (pt3_Idx != -999 && pt3_Idx < pt3_t3Idx_vec->size()) {
                 int t3_0 = pt3_t3Idx_vec->at(pt3_Idx);
                 if (t3_0 != -999 && t3_0 < t3_lsIdx0_vec->size()) {
                    if (t3_lsIdx0_vec->at(t3_0) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx0_vec->at(t3_0));
                    if (t3_lsIdx1_vec->at(t3_0) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx1_vec->at(t3_0));
                 }
                 
            }
        }

        // Path 3: tc_t5 -> t5 -> t3 -> ls
        for (int t5_Idx : *tc_t5Idx_vec) {
            if (t5_Idx != -999 && t5_Idx < t5_t3Idx0_vec->size()) {
                int t3_0 = t5_t3Idx0_vec->at(t5_Idx);
                if (t3_0 != -999 && t3_0 < t3_lsIdx0_vec->size()) {
                    if (t3_lsIdx0_vec->at(t3_0) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx0_vec->at(t3_0));
                    if (t3_lsIdx1_vec->at(t3_0) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx1_vec->at(t3_0));
                }
                int t3_1 = t5_t3Idx1_vec->at(t5_Idx);
                if (t3_1 != -999 && t3_1 < t3_lsIdx0_vec->size()) {
                    if (t3_lsIdx0_vec->at(t3_1) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx0_vec->at(t3_1));
                    if (t3_lsIdx1_vec->at(t3_1) != -999) used_ls_indices_in_this_event.insert(t3_lsIdx1_vec->at(t3_1));
                }
            }
        }
        
        // Loop through this event's ls objects and fill histograms
        for (size_t j = 0; j < ls_pt_vec->size(); ++j) {
            bool is_fake = ls_isFake_vec->at(j) != 0;
            bool is_used = used_ls_indices_in_this_event.count(j) > 0;

            float pt = ls_pt_vec->at(j), eta = ls_eta_vec->at(j), phi = ls_phi_vec->at(j);

            if (!is_fake) { // REAL
                if (is_used) { h_pt_real_used->Fill(pt); h_eta_real_used->Fill(eta); h_phi_real_used->Fill(phi); } 
                else { h_pt_real_unused->Fill(pt); h_eta_real_unused->Fill(eta); h_phi_real_unused->Fill(phi); }
            } else { // FAKE
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
    h_pt_real_used->SetFillColor(kBlue-7); h_pt_fake_used->SetFillColor(kRed-7);
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
    h_eta_real_used->SetFillColor(kGreen+2); h_eta_fake_used->SetFillColor(kOrange-3);
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
    hs_phi_used->Add(h_phi_fake_used); hs_phi_used->Add(h_phi_real_used);
    hs_phi_used->Draw("HIST");
    auto legend_phi_used = new TLegend(0.4, 0.8, 0.7, 0.9);
    legend_phi_used->AddEntry(h_phi_real_used, "Real, Used", "f");
    legend_phi_used->AddEntry(h_phi_fake_used, "Fake, Used", "f");
    legend_phi_used->Draw();
    c_phi_used->SaveAs("compare_phi_used.png");

    TCanvas *c_phi_unused = new TCanvas("c_phi_unused", "Phi of Unused Tracks", 800, 600);
    THStack *hs_phi_unused = new THStack("hs_phi_unused", "Phi Distribution of Unused Tracks;#phi [rad];Tracks");
    h_phi_real_unused->SetFillColor(kMagenta-9); h_phi_fake_unused->SetFillColor(kCyan-9);
    hs_phi_unused->Add(h_phi_fake_unused); hs_phi_unused->Add(h_phi_real_unused);
    hs_phi_unused->Draw("HIST");
    auto legend_phi_unused = new TLegend(0.4, 0.8, 0.7, 0.9);
    legend_phi_unused->AddEntry(h_phi_real_unused, "Real, Unused", "f");
    legend_phi_unused->AddEntry(h_phi_fake_unused, "Fake, Unused", "f");
    legend_phi_unused->Draw();
    c_phi_unused->SaveAs("compare_phi_unused.png");

    // --- Clean up ---
    inputFile->Close();
    delete c_pt_used; delete c_pt_unused;
    delete c_eta_used; delete c_eta_unused;
    delete c_phi_used; delete c_phi_unused;
}
