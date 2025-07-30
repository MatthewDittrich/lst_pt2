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
#include <map>

// Include the file with our helper functions
#include "modules.C"

void main_matching() {
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
    std::vector<float> *ls_pt_vec = nullptr, *ls_eta_vec = nullptr, *ls_phi_vec = nullptr;
    std::vector<int>   *ls_isFake_vec = nullptr, *ls_simIdx_vec = nullptr;
    std::vector<int>   *tc_pt5Idx_vec = nullptr, *pt5_t5Idx_vec = nullptr;
    std::vector<int>   *t5_t3Idx0_vec = nullptr, *t5_t3Idx1_vec = nullptr;
    std::vector<int>   *t3_lsIdx0_vec = nullptr, *t3_lsIdx1_vec = nullptr;
    std::vector<int>   *tc_pt3Idx_vec = nullptr, *pt3_t3Idx_vec = nullptr;
    std::vector<int>   *tc_t5Idx_vec = nullptr, *pt5_plsIdx_vec = nullptr, *pt3_plsIdx_vec = nullptr;
    std::vector<int>   *pls_isFake_vec = nullptr, *pls_simIdx_vec = nullptr;

    inputTree->SetBranchAddress("ls_pt", &ls_pt_vec);
    inputTree->SetBranchAddress("ls_eta", &ls_eta_vec);
    inputTree->SetBranchAddress("ls_phi", &ls_phi_vec);
    inputTree->SetBranchAddress("ls_isFake", &ls_isFake_vec);
    inputTree->SetBranchAddress("ls_simIdx", &ls_simIdx_vec);
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
    inputTree->SetBranchAddress("pls_isFake", &pls_isFake_vec);
    inputTree->SetBranchAddress("pls_simIdx", &pls_simIdx_vec);

    // --- STEP 3: Define histograms for all categories ---
    int pt_bins = 180, eta_bins = 180, phi_bins = 180;
    double pt_max = 250.0, eta_max = 4.5, phi_max = 3.15;

    // Histograms for Fake LS
    TH1D *h_pt_fake_used = new TH1D("h_pt_fake_used", "", pt_bins, 0, pt_max), *h_pt_fake_unused = new TH1D("h_pt_fake_unused", "", pt_bins, 0, pt_max);
    TH1D *h_eta_fake_used = new TH1D("h_eta_fake_used", "", eta_bins, -eta_max, eta_max), *h_eta_fake_unused = new TH1D("h_eta_fake_unused", "", eta_bins, -eta_max, eta_max);
    TH1D *h_phi_fake_used = new TH1D("h_phi_fake_used", "", phi_bins, -phi_max, phi_max), *h_phi_fake_unused = new TH1D("h_phi_fake_unused", "", phi_bins, -phi_max, phi_max);

    // Histograms for Real LS (USED)
    TH1D *h_pt_real_used_matchedToUsedOnly = new TH1D("h_pt_real_used_matchedToUsedOnly", "", pt_bins, 0, pt_max), *h_pt_real_used_matchedToUnusedOnly = new TH1D("h_pt_real_used_matchedToUnusedOnly", "", pt_bins, 0, pt_max), *h_pt_real_used_matchedToBoth = new TH1D("h_pt_real_used_matchedToBoth", "", pt_bins, 0, pt_max), *h_pt_real_used_unmatched = new TH1D("h_pt_real_used_unmatched", "", pt_bins, 0, pt_max);
    
    TH1D *h_eta_real_used_matchedToUsedOnly = new TH1D("h_eta_real_used_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_used_matchedToUnusedOnly = new TH1D("h_eta_real_used_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_used_matchedToBoth = new TH1D("h_eta_real_used_matchedToBoth", "", eta_bins, -eta_max, eta_max), *h_eta_real_used_unmatched = new TH1D("h_eta_real_used_unmatched", "", eta_bins, -eta_max, eta_max);
    
    TH1D *h_phi_real_used_matchedToUsedOnly = new TH1D("h_phi_real_used_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_used_matchedToUnusedOnly = new TH1D("h_phi_real_used_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_used_matchedToBoth = new TH1D("h_phi_real_used_matchedToBoth", "", phi_bins, -phi_max, phi_max), *h_phi_real_used_unmatched = new TH1D("h_phi_real_used_unmatched", "", phi_bins, -phi_max, phi_max);

    // Histograms for Real LS (UNUSED)
    TH1D *h_pt_real_unused_matchedToUsedOnly = new TH1D("h_pt_real_unused_matchedToUsedOnly", "", pt_bins, 0, pt_max), *h_pt_real_unused_matchedToUnusedOnly = new TH1D("h_pt_real_unused_matchedToUnusedOnly", "", pt_bins, 0, pt_max), *h_pt_real_unused_matchedToBoth = new TH1D("h_pt_real_unused_matchedToBoth", "", pt_bins, 0, pt_max), *h_pt_real_unused_unmatched = new TH1D("h_pt_real_unused_unmatched", "", pt_bins, 0, pt_max);
    
    TH1D *h_eta_real_unused_matchedToUsedOnly = new TH1D("h_eta_real_unused_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_unused_matchedToUnusedOnly = new TH1D("h_eta_real_unused_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_unused_matchedToBoth = new TH1D("h_eta_real_unused_matchedToBoth", "", eta_bins, -eta_max, eta_max), *h_eta_real_unused_unmatched = new TH1D("h_eta_real_unused_unmatched", "", eta_bins, -eta_max, eta_max);
    
    TH1D *h_phi_real_unused_matchedToUsedOnly = new TH1D("h_phi_real_unused_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_unused_matchedToUnusedOnly = new TH1D("h_phi_real_unused_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_unused_matchedToBoth = new TH1D("h_phi_real_unused_matchedToBoth", "", phi_bins, -phi_max, phi_max), *h_phi_real_unused_unmatched = new TH1D("h_phi_real_unused_unmatched", "", phi_bins, -phi_max, phi_max);

    // --- STEP 4: The Main Event Loop ---
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i);

        std::set<int> used_ls_indices, used_pls_indices;
        findUsedIndices(used_ls_indices, used_pls_indices, *tc_pt5Idx_vec, *pt5_t5Idx_vec, *t5_t3Idx0_vec, *t5_t3Idx1_vec, *t3_lsIdx0_vec, *t3_lsIdx1_vec, *tc_pt3Idx_vec, *pt3_t3Idx_vec, *tc_t5Idx_vec, *pt5_plsIdx_vec, *pt3_plsIdx_vec);

        std::map<int, PlsMatchStatus> simIdx_to_pls_status;
        getPlsMatchStatus(simIdx_to_pls_status, *pls_isFake_vec, *pls_simIdx_vec, used_pls_indices);

        for (size_t j = 0; j < ls_pt_vec->size(); ++j) {
            fillHistograms(
                (used_ls_indices.count(j) > 0), (ls_isFake_vec->at(j) != 0),
                ls_pt_vec->at(j), ls_eta_vec->at(j), ls_phi_vec->at(j), ls_simIdx_vec->at(j),
                simIdx_to_pls_status,
                h_pt_fake_used, h_pt_fake_unused, h_eta_fake_used, h_eta_fake_unused, h_phi_fake_used, h_phi_fake_unused,
                h_pt_real_used_unmatched, h_pt_real_used_matchedToBoth, h_pt_real_used_matchedToUsedOnly, h_pt_real_used_matchedToUnusedOnly,
                h_eta_real_used_unmatched, h_eta_real_used_matchedToBoth, h_eta_real_used_matchedToUsedOnly, h_eta_real_used_matchedToUnusedOnly,
                h_phi_real_used_unmatched, h_phi_real_used_matchedToBoth, h_phi_real_used_matchedToUsedOnly, h_phi_real_used_matchedToUnusedOnly,
                h_pt_real_unused_unmatched, h_pt_real_unused_matchedToBoth, h_pt_real_unused_matchedToUsedOnly, h_pt_real_unused_matchedToUnusedOnly,
                h_eta_real_unused_unmatched, h_eta_real_unused_matchedToBoth, h_eta_real_unused_matchedToUsedOnly, h_eta_real_unused_matchedToUnusedOnly,
                h_phi_real_unused_unmatched, h_phi_real_unused_matchedToBoth, h_phi_real_unused_matchedToUsedOnly, h_phi_real_unused_matchedToUnusedOnly
            );
        }
    }
    std::cout << "Finished processing entries." << std::endl;

    // --- STEP 5: Styling and Plotting ---
    gStyle->SetOptStat(0);

    // pT Plots
    createAndSavePlot("c_pt_used", "pT of Used LS;p_{T} [GeV];LS Objects", "compare_pt_used_detailed.png",
                      "Used LS Categories", 0.45, 0.58, 0.88, 0.88, true,
                      h_pt_fake_used, h_pt_real_used_unmatched, h_pt_real_used_matchedToUnusedOnly, h_pt_real_used_matchedToBoth, h_pt_real_used_matchedToUsedOnly);
    createAndSavePlot("c_pt_unused", "pT of Unused LS;p_{T} [GeV];LS Objects", "compare_pt_unused_detailed.png",
                      "Unused LS Categories", 0.45, 0.58, 0.88, 0.88, true,
                      h_pt_fake_unused, h_pt_real_unused_unmatched, h_pt_real_unused_matchedToUnusedOnly, h_pt_real_unused_matchedToBoth, h_pt_real_unused_matchedToUsedOnly);

    // Eta Plots
    createAndSavePlot("c_eta_used", "Eta of Used LS;#eta;LS Objects", "compare_eta_used_detailed.png",
                      "Used LS Categories", 0.15, 0.65, 0.55, 0.88, false,
                      h_eta_fake_used, h_eta_real_used_unmatched, h_eta_real_used_matchedToUnusedOnly, h_eta_real_used_matchedToBoth, h_eta_real_used_matchedToUsedOnly);
    createAndSavePlot("c_eta_unused", "Eta of Unused LS;#eta;LS Objects", "compare_eta_unused_detailed.png",
                      "Unused LS Categories", 0.15, 0.65, 0.55, 0.88, false,
                      h_eta_fake_unused, h_eta_real_unused_unmatched, h_eta_real_unused_matchedToUnusedOnly, h_eta_real_unused_matchedToBoth, h_eta_real_unused_matchedToUsedOnly);

    // Phi Plots
    createAndSavePlot("c_phi_used", "Phi of Used LS;#phi [rad];LS Objects", "compare_phi_used_detailed.png",
                      "Used LS Categories", 0.15, 0.65, 0.55, 0.88, false,
                      h_phi_fake_used, h_phi_real_used_unmatched, h_phi_real_used_matchedToUnusedOnly, h_phi_real_used_matchedToBoth, h_phi_real_used_matchedToUsedOnly);
    createAndSavePlot("c_phi_unused", "Phi of Unused LS;#phi [rad];LS Objects", "compare_phi_unused_detailed.png",
                      "Unused LS Categories", 0.15, 0.65, 0.55, 0.88, false,
                      h_phi_fake_unused, h_phi_real_unused_unmatched, h_phi_real_unused_matchedToUnusedOnly, h_phi_real_unused_matchedToBoth, h_phi_real_unused_matchedToUsedOnly);


    // --- Clean up ---
    inputFile->Close();
    // Histograms are owned by the current directory (gDirectory) and canvases were deleted in the plot function, so minimal cleanup is needed here.
}
