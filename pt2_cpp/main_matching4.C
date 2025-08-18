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
#include "draw_maker.C"
// Include the file with our helper functions
#include "modules4.C"

// Forward declaration of the gator function if it's not in the modules file
void print_coding_gator();

void main_matching4() {
    // --- Configuration & File Opening ---
    const char* input_filename = "LSTNtuple.root";
    const char* tree_name = "tree";
    TFile *inputFile = TFile::Open(input_filename, "READ");
    if (!inputFile || inputFile->IsZombie()) { std::cerr << "Error: Could not open file " << input_filename << std::endl; return; }
    TTree *inputTree;
    inputFile->GetObject(tree_name, inputTree);
    if (!inputTree) { std::cerr << "Error: Could not find TTree with name " << tree_name << std::endl; inputFile->Close(); return; }

    // --- Set Branch Addresses ---
    std::vector<float> *ls_pt_vec = nullptr, *ls_eta_vec = nullptr, *ls_phi_vec = nullptr;
    std::vector<int>   *ls_isFake_vec = nullptr, *ls_simIdx_vec = nullptr, *tc_pt5Idx_vec = nullptr, *pt5_t5Idx_vec = nullptr, *t5_t3Idx0_vec = nullptr, *t5_t3Idx1_vec = nullptr, *t3_lsIdx0_vec = nullptr, *t3_lsIdx1_vec = nullptr;
    std::vector<int>   *tc_pt3Idx_vec = nullptr, *pt3_t3Idx_vec = nullptr, *tc_t5Idx_vec = nullptr, *pt5_plsIdx_vec = nullptr, *pt3_plsIdx_vec = nullptr;
    std::vector<float> *pls_pt_vec = nullptr, *pls_eta_vec = nullptr, *pls_phi_vec = nullptr;
    std::vector<int>   *pls_isFake_vec = nullptr, *pls_simIdx_vec = nullptr;
    std::vector<int>   *ls_mdIdx0_vec = nullptr, *ls_mdIdx1_vec = nullptr;
    std::vector<float> *md_anchor_x_vec = nullptr, *md_anchor_y_vec = nullptr, *md_anchor_z_vec = nullptr;
    std::vector<float> *md_other_x_vec = nullptr, *md_other_y_vec = nullptr, *md_other_z_vec = nullptr;
    std::vector<float> *pls_hit0_x_vec = nullptr, *pls_hit0_y_vec = nullptr, *pls_hit0_z_vec = nullptr;
    std::vector<float> *pls_hit1_x_vec = nullptr, *pls_hit1_y_vec = nullptr, *pls_hit1_z_vec = nullptr;
    std::vector<float> *pls_hit2_x_vec = nullptr, *pls_hit2_y_vec = nullptr, *pls_hit2_z_vec = nullptr;
    std::vector<float> *pls_hit3_x_vec = nullptr, *pls_hit3_y_vec = nullptr, *pls_hit3_z_vec = nullptr;
    
    inputTree->SetBranchAddress("ls_pt", &ls_pt_vec); inputTree->SetBranchAddress("ls_eta", &ls_eta_vec); inputTree->SetBranchAddress("ls_phi", &ls_phi_vec);
    inputTree->SetBranchAddress("ls_isFake", &ls_isFake_vec); inputTree->SetBranchAddress("ls_simIdx", &ls_simIdx_vec);
    inputTree->SetBranchAddress("tc_pt5Idx", &tc_pt5Idx_vec); inputTree->SetBranchAddress("pt5_t5Idx", &pt5_t5Idx_vec);
    inputTree->SetBranchAddress("t5_t3Idx0", &t5_t3Idx0_vec); inputTree->SetBranchAddress("t5_t3Idx1", &t5_t3Idx1_vec);
    inputTree->SetBranchAddress("t3_lsIdx0", &t3_lsIdx0_vec); inputTree->SetBranchAddress("t3_lsIdx1", &t3_lsIdx1_vec);
    inputTree->SetBranchAddress("tc_pt3Idx", &tc_pt3Idx_vec); inputTree->SetBranchAddress("pt3_t3Idx", &pt3_t3Idx_vec);
    inputTree->SetBranchAddress("tc_t5Idx", &tc_t5Idx_vec); inputTree->SetBranchAddress("pt5_plsIdx", &pt5_plsIdx_vec);
    inputTree->SetBranchAddress("pt3_plsIdx", &pt3_plsIdx_vec);
    inputTree->SetBranchAddress("pls_pt", &pls_pt_vec); inputTree->SetBranchAddress("pls_eta", &pls_eta_vec); inputTree->SetBranchAddress("pls_phi", &pls_phi_vec);
    inputTree->SetBranchAddress("pls_isFake", &pls_isFake_vec); inputTree->SetBranchAddress("pls_simIdx", &pls_simIdx_vec);
    inputTree->SetBranchAddress("ls_mdIdx0", &ls_mdIdx0_vec); inputTree->SetBranchAddress("ls_mdIdx1", &ls_mdIdx1_vec);
    inputTree->SetBranchAddress("md_anchor_x", &md_anchor_x_vec); inputTree->SetBranchAddress("md_anchor_y", &md_anchor_y_vec);
    inputTree->SetBranchAddress("md_anchor_z", &md_anchor_z_vec); inputTree->SetBranchAddress("md_other_x", &md_other_x_vec);
    inputTree->SetBranchAddress("md_other_y", &md_other_y_vec); inputTree->SetBranchAddress("md_other_z", &md_other_z_vec);
    inputTree->SetBranchAddress("pls_hit0_x", &pls_hit0_x_vec); inputTree->SetBranchAddress("pls_hit0_y", &pls_hit0_y_vec);
    inputTree->SetBranchAddress("pls_hit0_z", &pls_hit0_z_vec);
    inputTree->SetBranchAddress("pls_hit1_x", &pls_hit1_x_vec); inputTree->SetBranchAddress("pls_hit1_y", &pls_hit1_y_vec);
    inputTree->SetBranchAddress("pls_hit1_z", &pls_hit1_z_vec);
    inputTree->SetBranchAddress("pls_hit2_x", &pls_hit2_x_vec); inputTree->SetBranchAddress("pls_hit2_y", &pls_hit2_y_vec);
    inputTree->SetBranchAddress("pls_hit2_z", &pls_hit2_z_vec);
    inputTree->SetBranchAddress("pls_hit3_x", &pls_hit3_x_vec); inputTree->SetBranchAddress("pls_hit3_y", &pls_hit3_y_vec);
    inputTree->SetBranchAddress("pls_hit3_z", &pls_hit3_z_vec);

    // --- Define Histograms ---
    int pt_bins = 180, eta_bins = 180, phi_bins = 180;
    double pt_max = 250.0, eta_max = 4.5, phi_max = 3.15;

    TH1D *h_pt_real_ls_used = new TH1D("h_pt_real_ls_used", "", pt_bins, 0, pt_max), *h_pt_fake_ls_used = new TH1D("h_pt_fake_ls_used", "", pt_bins, 0, pt_max), *h_pt_real_ls_unused = new TH1D("h_pt_real_ls_unused", "", pt_bins, 0, pt_max), *h_pt_fake_ls_unused = new TH1D("h_pt_fake_ls_unused", "", pt_bins, 0, pt_max), *h_eta_real_ls_used = new TH1D("h_eta_real_ls_used", "", eta_bins, -eta_max, eta_max), *h_eta_fake_ls_used = new TH1D("h_eta_fake_ls_used", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_unused = new TH1D("h_eta_real_ls_unused", "", eta_bins, -eta_max, eta_max), *h_eta_fake_ls_unused = new TH1D("h_eta_fake_ls_unused", "", eta_bins, -eta_max, eta_max), *h_phi_real_ls_used = new TH1D("h_phi_real_ls_used", "", phi_bins, -phi_max, phi_max), *h_phi_fake_ls_used = new TH1D("h_phi_fake_ls_used", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_unused = new TH1D("h_phi_real_ls_unused", "", phi_bins, -phi_max, phi_max), *h_phi_fake_ls_unused = new TH1D("h_phi_fake_ls_unused", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_real_pls_used = new TH1D("h_pt_real_pls_used", "", pt_bins, 0, pt_max), *h_pt_fake_pls_used = new TH1D("h_pt_fake_pls_used", "", pt_bins, 0, pt_max), *h_pt_real_pls_unused = new TH1D("h_pt_real_pls_unused", "", pt_bins, 0, pt_max), *h_pt_fake_pls_unused = new TH1D("h_pt_fake_pls_unused", "", pt_bins, 0, pt_max), *h_eta_real_pls_used = new TH1D("h_eta_real_pls_used", "", eta_bins, -eta_max, eta_max), *h_eta_fake_pls_used = new TH1D("h_eta_fake_pls_used", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_unused = new TH1D("h_eta_real_pls_unused", "", eta_bins, -eta_max, eta_max), *h_eta_fake_pls_unused = new TH1D("h_eta_fake_pls_unused", "", eta_bins, -eta_max, eta_max), *h_phi_real_pls_used = new TH1D("h_phi_real_pls_used", "", phi_bins, -phi_max, phi_max), *h_phi_fake_pls_used = new TH1D("h_phi_fake_pls_used", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_unused = new TH1D("h_phi_real_pls_unused", "", phi_bins, -phi_max, phi_max), *h_phi_fake_pls_unused = new TH1D("h_phi_fake_pls_unused", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_fake_ls_used_detailed = new TH1D("h_pt_fake_ls_used_detailed", "", pt_bins, 0, pt_max), *h_pt_fake_ls_unused_detailed = new TH1D("h_pt_fake_ls_unused_detailed", "", pt_bins, 0, pt_max), *h_eta_fake_ls_used_detailed = new TH1D("h_eta_fake_ls_used_detailed", "", eta_bins, -eta_max, eta_max), *h_eta_fake_ls_unused_detailed = new TH1D("h_eta_fake_ls_unused_detailed", "", eta_bins, -eta_max, eta_max), *h_phi_fake_ls_used_detailed = new TH1D("h_phi_fake_ls_used_detailed", "", phi_bins, -phi_max, phi_max), *h_phi_fake_ls_unused_detailed = new TH1D("h_phi_fake_ls_unused_detailed", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_real_ls_used_unmatched = new TH1D("h_pt_real_ls_used_unmatched", "", pt_bins, 0, pt_max), *h_pt_real_ls_used_matchedToBoth = new TH1D("h_pt_real_ls_used_matchedToBoth", "", pt_bins, 0, pt_max), *h_pt_real_ls_used_matchedToUsedOnly = new TH1D("h_pt_real_ls_used_matchedToUsedOnly", "", pt_bins, 0, pt_max), *h_pt_real_ls_used_matchedToUnusedOnly = new TH1D("h_pt_real_ls_used_matchedToUnusedOnly", "", pt_bins, 0, pt_max), *h_eta_real_ls_used_unmatched = new TH1D("h_eta_real_ls_used_unmatched", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_used_matchedToBoth = new TH1D("h_eta_real_ls_used_matchedToBoth", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_used_matchedToUsedOnly = new TH1D("h_eta_real_ls_used_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_used_matchedToUnusedOnly = new TH1D("h_eta_real_ls_used_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max), *h_phi_real_ls_used_unmatched = new TH1D("h_phi_real_ls_used_unmatched", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_used_matchedToBoth = new TH1D("h_phi_real_ls_used_matchedToBoth", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_used_matchedToUsedOnly = new TH1D("h_phi_real_ls_used_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_used_matchedToUnusedOnly = new TH1D("h_phi_real_ls_used_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_real_ls_unused_unmatched = new TH1D("h_pt_real_ls_unused_unmatched", "", pt_bins, 0, pt_max), *h_pt_real_ls_unused_matchedToBoth = new TH1D("h_pt_real_ls_unused_matchedToBoth", "", pt_bins, 0, pt_max), *h_pt_real_ls_unused_matchedToUsedOnly = new TH1D("h_pt_real_ls_unused_matchedToUsedOnly", "", pt_bins, 0, pt_max), *h_pt_real_ls_unused_matchedToUnusedOnly = new TH1D("h_pt_real_ls_unused_matchedToUnusedOnly", "", pt_bins, 0, pt_max), *h_eta_real_ls_unused_unmatched = new TH1D("h_eta_real_ls_unused_unmatched", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_unused_matchedToBoth = new TH1D("h_eta_real_ls_unused_matchedToBoth", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_unused_matchedToUsedOnly = new TH1D("h_eta_real_ls_unused_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_ls_unused_matchedToUnusedOnly = new TH1D("h_eta_real_ls_unused_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max), *h_phi_real_ls_unused_unmatched = new TH1D("h_phi_real_ls_unused_unmatched", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_unused_matchedToBoth = new TH1D("h_phi_real_ls_unused_matchedToBoth", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_unused_matchedToUsedOnly = new TH1D("h_phi_real_ls_unused_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_ls_unused_matchedToUnusedOnly = new TH1D("h_phi_real_ls_unused_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_fake_pls_used_detailed = new TH1D("h_pt_fake_pls_used_detailed", "", pt_bins, 0, pt_max), *h_pt_fake_pls_unused_detailed = new TH1D("h_pt_fake_pls_unused_detailed", "", pt_bins, 0, pt_max), *h_eta_fake_pls_used_detailed = new TH1D("h_eta_fake_pls_used_detailed", "", eta_bins, -eta_max, eta_max), *h_eta_fake_pls_unused_detailed = new TH1D("h_eta_fake_pls_unused_detailed", "", eta_bins, -eta_max, eta_max), *h_phi_fake_pls_used_detailed = new TH1D("h_phi_fake_pls_used_detailed", "", phi_bins, -phi_max, phi_max), *h_phi_fake_pls_unused_detailed = new TH1D("h_phi_fake_pls_unused_detailed", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_real_pls_used_unmatched = new TH1D("h_pt_real_pls_used_unmatched", "", pt_bins, 0, pt_max), *h_pt_real_pls_used_matchedToBoth = new TH1D("h_pt_real_pls_used_matchedToBoth", "", pt_bins, 0, pt_max), *h_pt_real_pls_used_matchedToUsedOnly = new TH1D("h_pt_real_pls_used_matchedToUsedOnly", "", pt_bins, 0, pt_max), *h_pt_real_pls_used_matchedToUnusedOnly = new TH1D("h_pt_real_pls_used_matchedToUnusedOnly", "", pt_bins, 0, pt_max), *h_eta_real_pls_used_unmatched = new TH1D("h_eta_real_pls_used_unmatched", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_used_matchedToBoth = new TH1D("h_eta_real_pls_used_matchedToBoth", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_used_matchedToUsedOnly = new TH1D("h_eta_real_pls_used_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_used_matchedToUnusedOnly = new TH1D("h_eta_real_pls_used_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max), *h_phi_real_pls_used_unmatched = new TH1D("h_phi_real_pls_used_unmatched", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_used_matchedToBoth = new TH1D("h_phi_real_pls_used_matchedToBoth", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_used_matchedToUsedOnly = new TH1D("h_phi_real_pls_used_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_used_matchedToUnusedOnly = new TH1D("h_phi_real_pls_used_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);
    TH1D *h_pt_real_pls_unused_unmatched = new TH1D("h_pt_real_pls_unused_unmatched", "", pt_bins, 0, pt_max), *h_pt_real_pls_unused_matchedToBoth = new TH1D("h_pt_real_pls_unused_matchedToBoth", "", pt_bins, 0, pt_max), *h_pt_real_pls_unused_matchedToUsedOnly = new TH1D("h_pt_real_pls_unused_matchedToUsedOnly", "", pt_bins, 0, pt_max), *h_pt_real_pls_unused_matchedToUnusedOnly = new TH1D("h_pt_real_pls_unused_matchedToUnusedOnly", "", pt_bins, 0, pt_max), *h_eta_real_pls_unused_unmatched = new TH1D("h_eta_real_pls_unused_unmatched", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_unused_matchedToBoth = new TH1D("h_eta_real_pls_unused_matchedToBoth", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_unused_matchedToUsedOnly = new TH1D("h_eta_real_pls_unused_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max), *h_eta_real_pls_unused_matchedToUnusedOnly = new TH1D("h_eta_real_pls_unused_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max), *h_phi_real_pls_unused_unmatched = new TH1D("h_phi_real_pls_unused_unmatched", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_unused_matchedToBoth = new TH1D("h_phi_real_pls_unused_matchedToBoth", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_unused_matchedToUsedOnly = new TH1D("h_phi_real_pls_unused_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max), *h_phi_real_pls_unused_matchedToUnusedOnly = new TH1D("h_phi_real_pls_unused_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);
    TH1D* h_delta_pt = new TH1D("h_delta_pt", "#Delta p_{T} (LS - PLS);#Delta p_{T} [GeV];pT2 Objects", 180, -5, 5);
    TH1D* h_delta_eta = new TH1D("h_delta_eta", "#Delta #eta (LS - PLS);#Delta #eta;pT2 Objects", 180, -0.7, 0.7);
    TH1D* h_delta_phi = new TH1D("h_delta_phi", "#Delta #phi (LS - PLS);#Delta #phi [rad];pT2 Objects", 180, -1.25, 1.25);
    TH1D* h_delta_R = new TH1D("h_delta_R", "#Delta R (LS - PLS);#Delta R;pT2 Objects", 180, 0, 1);
    TH1D* h_extrapolation_dist_3d = new TH1D("h_extrapolation_dist_3d","3D Distance between extrapolated PLS and LS;Distance [cm];Ideal pT2 Pairs", 180, 0, 25.0);
    TH1D* h_extrapolation_delta_z = new TH1D("h_extrapolation_delta_z", "R-Z Extrapolation #Delta r (LS - PLS);#Delta r [cm];Ideal pT2 Pairs", 180, -15.0, 15.0);

    // --- The Main Event Loop ---
    Long64_t nEntries = inputTree->GetEntries();
    print_coding_gator();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i);
        std::set<int> used_ls_indices, used_pls_indices;
        findUsedIndices(used_ls_indices, used_pls_indices, *tc_pt5Idx_vec, *pt5_t5Idx_vec, *t5_t3Idx0_vec, *t5_t3Idx1_vec, *t3_lsIdx0_vec, *t3_lsIdx1_vec, *tc_pt3Idx_vec, *pt3_t3Idx_vec, *tc_t5Idx_vec, *pt5_plsIdx_vec, *pt3_plsIdx_vec);
        std::map<int, PlsMatchStatus> simIdx_to_pls_status;
        getPlsMatchStatus(simIdx_to_pls_status, *pls_isFake_vec, *pls_simIdx_vec, used_pls_indices);
        std::map<int, LsMatchStatus> simIdx_to_ls_status;
        getLsMatchStatus(simIdx_to_ls_status, *ls_isFake_vec, *ls_simIdx_vec, used_ls_indices);
        for (size_t j = 0; j < ls_pt_vec->size(); ++j) {
            bool is_used = (used_ls_indices.count(j) > 0);
            bool is_fake = (ls_isFake_vec->at(j) != 0);
            fillSimpleHistograms(is_used, is_fake, ls_pt_vec->at(j), ls_eta_vec->at(j), ls_phi_vec->at(j), h_pt_real_ls_used, h_pt_fake_ls_used, h_pt_real_ls_unused, h_pt_fake_ls_unused, h_eta_real_ls_used, h_eta_fake_ls_used, h_eta_real_ls_unused, h_eta_fake_ls_unused, h_phi_real_ls_used, h_phi_fake_ls_used, h_phi_real_ls_unused, h_phi_fake_ls_unused);
            fillDetailedLsHistograms(is_used, is_fake, ls_pt_vec->at(j), ls_eta_vec->at(j), ls_phi_vec->at(j), ls_simIdx_vec->at(j), simIdx_to_pls_status, h_pt_fake_ls_used_detailed, h_pt_fake_ls_unused_detailed, h_pt_real_ls_used_unmatched, h_pt_real_ls_used_matchedToBoth, h_pt_real_ls_used_matchedToUsedOnly, h_pt_real_ls_used_matchedToUnusedOnly, h_pt_real_ls_unused_unmatched, h_pt_real_ls_unused_matchedToBoth, h_pt_real_ls_unused_matchedToUsedOnly, h_pt_real_ls_unused_matchedToUnusedOnly, h_eta_fake_ls_used_detailed, h_eta_fake_ls_unused_detailed, h_eta_real_ls_used_unmatched, h_eta_real_ls_used_matchedToBoth, h_eta_real_ls_used_matchedToUsedOnly, h_eta_real_ls_used_matchedToUnusedOnly, h_eta_real_ls_unused_unmatched, h_eta_real_ls_unused_matchedToBoth, h_eta_real_ls_unused_matchedToUsedOnly, h_eta_real_ls_unused_matchedToUnusedOnly, h_phi_fake_ls_used_detailed, h_phi_fake_ls_unused_detailed, h_phi_real_ls_used_unmatched, h_phi_real_ls_used_matchedToBoth, h_phi_real_ls_used_matchedToUsedOnly, h_phi_real_ls_used_matchedToUnusedOnly, h_phi_real_ls_unused_unmatched, h_phi_real_ls_unused_matchedToBoth, h_phi_real_ls_unused_matchedToUsedOnly, h_phi_real_ls_unused_matchedToUnusedOnly);
        }
        for (size_t j = 0; j < pls_pt_vec->size(); ++j) {
            bool is_used = (used_pls_indices.count(j) > 0);
            bool is_fake = (pls_isFake_vec->at(j) != 0);
            fillSimpleHistograms(is_used, is_fake, pls_pt_vec->at(j), pls_eta_vec->at(j), pls_phi_vec->at(j), h_pt_real_pls_used, h_pt_fake_pls_used, h_pt_real_pls_unused, h_pt_fake_pls_unused, h_eta_real_pls_used, h_eta_fake_pls_used, h_eta_real_pls_unused, h_eta_fake_pls_unused, h_phi_real_pls_used, h_phi_fake_pls_used, h_phi_real_pls_unused, h_phi_fake_pls_unused);
            fillDetailedPlsHistograms(is_used, is_fake, pls_pt_vec->at(j), pls_eta_vec->at(j), pls_phi_vec->at(j), pls_simIdx_vec->at(j), simIdx_to_ls_status, h_pt_fake_pls_used_detailed, h_pt_fake_pls_unused_detailed, h_pt_real_pls_used_unmatched, h_pt_real_pls_used_matchedToBoth, h_pt_real_pls_used_matchedToUsedOnly, h_pt_real_pls_used_matchedToUnusedOnly, h_pt_real_pls_unused_unmatched, h_pt_real_pls_unused_matchedToBoth, h_pt_real_pls_unused_matchedToUsedOnly, h_pt_real_pls_unused_matchedToUnusedOnly, h_eta_fake_pls_used_detailed, h_eta_fake_pls_unused_detailed, h_eta_real_pls_used_unmatched, h_eta_real_pls_used_matchedToBoth, h_eta_real_pls_used_matchedToUsedOnly, h_eta_real_pls_used_matchedToUnusedOnly, h_eta_real_pls_unused_unmatched, h_eta_real_pls_unused_matchedToBoth, h_eta_real_pls_unused_matchedToUsedOnly, h_eta_real_pls_unused_matchedToUnusedOnly, h_phi_fake_pls_used_detailed, h_phi_fake_pls_unused_detailed, h_phi_real_pls_used_unmatched, h_phi_real_pls_used_matchedToBoth, h_phi_real_pls_used_matchedToUsedOnly, h_phi_real_pls_used_matchedToUnusedOnly, h_phi_real_pls_unused_unmatched, h_phi_real_pls_unused_matchedToBoth, h_phi_real_pls_unused_matchedToUsedOnly, h_phi_real_pls_unused_matchedToUnusedOnly);
        }
        processIdealPT2s(*pls_pt_vec, *pls_eta_vec, *pls_phi_vec, *pls_isFake_vec, *pls_simIdx_vec, *ls_pt_vec, *ls_eta_vec, *ls_phi_vec, *ls_isFake_vec, *ls_simIdx_vec, used_pls_indices, used_ls_indices, h_delta_pt, h_delta_eta, h_delta_phi, h_delta_R, h_extrapolation_dist_3d, h_extrapolation_delta_z, *pls_hit0_x_vec, *pls_hit0_y_vec, *pls_hit0_z_vec, *pls_hit1_x_vec, *pls_hit1_y_vec, *pls_hit1_z_vec, *pls_hit2_x_vec, *pls_hit2_y_vec, *pls_hit2_z_vec, *pls_hit3_x_vec, *pls_hit3_y_vec, *pls_hit3_z_vec, *ls_mdIdx0_vec, *ls_mdIdx1_vec, *md_anchor_x_vec, *md_anchor_y_vec, *md_anchor_z_vec, *md_other_x_vec, *md_other_y_vec, *md_other_z_vec);
    }
    std::cout << "Finished processing entries." << std::endl;
    // --- Styling and Plotting ---
    gStyle->SetOptStat(0);
    // Simple LS Plots
    createAndSaveSimpleStack("c_pt_ls_used", "pT of Used LS;p_{T} [GeV];LS Objects", "compare_pt_ls_used.png", "Used LS", h_pt_real_ls_used, h_pt_fake_ls_used, "Real, Used", "Fake, Used", kGreen+2, kOrange-3, true);
    createAndSaveSimpleStack("c_eta_ls_used", "Eta of Used LS;#eta;LS Objects", "compare_eta_ls_used.png", "Used LS", h_eta_real_ls_used, h_eta_fake_ls_used, "Real, Used", "Fake, Used", kGreen+2, kOrange-3);
    createAndSaveSimpleStack("c_phi_ls_used", "Phi of Used LS;#phi;LS Objects", "compare_phi_ls_used.png", "Used LS", h_phi_real_ls_used, h_phi_fake_ls_used, "Real, Used", "Fake, Used", kGreen+2, kOrange-3);
    createAndSaveSimpleStack("c_pt_ls_unused", "pT of Unused LS;p_{T} [GeV];LS Objects", "compare_pt_ls_unused.png", "Unused LS", h_pt_real_ls_unused, h_pt_fake_ls_unused, "Real, Unused", "Fake, Unused", kGreen-7, kOrange-9, true);
    createAndSaveSimpleStack("c_eta_ls_unused", "Eta of Unused LS;#eta;LS Objects", "compare_eta_ls_unused.png", "Unused LS", h_eta_real_ls_unused, h_eta_fake_ls_unused, "Real, Unused", "Fake, Unused", kGreen-7, kOrange-9);
    createAndSaveSimpleStack("c_phi_ls_unused", "Phi of Unused LS;#phi;LS Objects", "compare_phi_ls_unused.png", "Unused LS", h_phi_real_ls_unused, h_phi_fake_ls_unused, "Real, Unused", "Fake, Unused", kGreen-7, kOrange-9);
    // Simple PLS Plots
    createAndSaveSimpleStack("c_pt_pls_used", "pT of Used PLS;p_{T} [GeV];PLS Objects", "compare_pt_pls_used.png", "Used PLS", h_pt_real_pls_used, h_pt_fake_pls_used, "Real, Used", "Fake, Used", kAzure+1, kRed-7, true);
    createAndSaveSimpleStack("c_eta_pls_used", "Eta of Used PLS;#eta;PLS Objects", "compare_eta_pls_used.png", "Used PLS", h_eta_real_pls_used, h_eta_fake_pls_used, "Real, Used", "Fake, Used", kAzure+1, kRed-7);
    createAndSaveSimpleStack("c_phi_pls_used", "Phi of Used PLS;#phi;PLS Objects", "compare_phi_pls_used.png", "Used PLS", h_phi_real_pls_used, h_phi_fake_pls_used, "Real, Used", "Fake, Used", kAzure+1, kRed-7);
    createAndSaveSimpleStack("c_pt_pls_unused", "pT of Unused PLS;p_{T} [GeV];PLS Objects", "compare_pt_pls_unused.png", "Unused PLS", h_pt_real_pls_unused, h_pt_fake_pls_unused, "Real, Unused", "Fake, Unused", kAzure-9, kRed-9, true);
    createAndSaveSimpleStack("c_eta_pls_unused", "Eta of Unused PLS;#eta;PLS Objects", "compare_eta_pls_unused.png", "Unused PLS", h_eta_real_pls_unused, h_eta_fake_pls_unused, "Real, Unused", "Fake, Unused", kAzure-9, kRed-9);
    createAndSaveSimpleStack("c_phi_pls_unused", "Phi of Unused PLS;#phi;PLS Objects", "compare_phi_pls_unused.png", "Unused PLS", h_phi_real_pls_unused, h_phi_fake_pls_unused, "Real, Unused", "Fake, Unused", kAzure-9, kRed-9);
    // Detailed LS Plots
    createAndSaveDetailedStack("c_pt_ls_used_detailed", "pT of Used LS (Detailed);p_{T} [GeV];LS Objects", "compare_pt_ls_used_detailed.png", "Used Real LS Matched To PLS:", true, h_pt_fake_ls_used_detailed, h_pt_real_ls_used_unmatched, h_pt_real_ls_used_matchedToUnusedOnly, h_pt_real_ls_used_matchedToBoth, h_pt_real_ls_used_matchedToUsedOnly);
    createAndSaveDetailedStack("c_eta_ls_used_detailed", "Eta of Used LS (Detailed);#eta;LS Objects", "compare_eta_ls_used_detailed.png", "Used Real LS Matched To PLS:", false, h_eta_fake_ls_used_detailed, h_eta_real_ls_used_unmatched, h_eta_real_ls_used_matchedToUnusedOnly, h_eta_real_ls_used_matchedToBoth, h_eta_real_ls_used_matchedToUsedOnly);
    createAndSaveDetailedStack("c_phi_ls_used_detailed", "Phi of Used LS (Detailed);#phi;LS Objects", "compare_phi_ls_used_detailed.png", "Used Real LS Matched To PLS:", false, h_phi_fake_ls_used_detailed, h_phi_real_ls_used_unmatched, h_phi_real_ls_used_matchedToUnusedOnly, h_phi_real_ls_used_matchedToBoth, h_phi_real_ls_used_matchedToUsedOnly);
    createAndSaveDetailedStack("c_pt_ls_unused_detailed", "pT of Unused LS (Detailed);p_{T} [GeV];LS Objects", "compare_pt_ls_unused_detailed.png", "Unused Real LS Matched To PLS:", true, h_pt_fake_ls_unused_detailed, h_pt_real_ls_unused_unmatched, h_pt_real_ls_unused_matchedToUnusedOnly, h_pt_real_ls_unused_matchedToBoth, h_pt_real_ls_unused_matchedToUsedOnly);
    createAndSaveDetailedStack("c_eta_ls_unused_detailed", "Eta of Unused LS (Detailed);#eta;LS Objects", "compare_eta_ls_unused_detailed.png", "Unused Real LS Matched To PLS:", false, h_eta_fake_ls_unused_detailed, h_eta_real_ls_unused_unmatched, h_eta_real_ls_unused_matchedToUnusedOnly, h_eta_real_ls_unused_matchedToBoth, h_eta_real_ls_unused_matchedToUsedOnly);
    createAndSaveDetailedStack("c_phi_ls_unused_detailed", "Phi of Unused LS (Detailed);#phi;LS Objects", "compare_phi_ls_unused_detailed.png", "Unused Real LS Matched To PLS:", false, h_phi_fake_ls_unused_detailed, h_phi_real_ls_unused_unmatched, h_phi_real_ls_unused_matchedToUnusedOnly, h_phi_real_ls_unused_matchedToBoth, h_phi_real_ls_unused_matchedToUsedOnly);
    
    // Detailed PLS Plots
    createAndSaveDetailedStack("c_pt_pls_used_detailed", "pT of Used PLS (Detailed);p_{T} [GeV];PLS Objects", "compare_pt_pls_used_detailed.png", "Used Real PLS Matched To LS:", true, h_pt_fake_pls_used_detailed, h_pt_real_pls_used_unmatched, h_pt_real_pls_used_matchedToUnusedOnly, h_pt_real_pls_used_matchedToBoth, h_pt_real_pls_used_matchedToUsedOnly);
    createAndSaveDetailedStack("c_eta_pls_used_detailed", "Eta of Used PLS (Detailed);#eta;PLS Objects", "compare_eta_pls_used_detailed.png", "Used Real PLS Matched To LS:", false, h_eta_fake_pls_used_detailed, h_eta_real_pls_used_unmatched, h_eta_real_pls_used_matchedToUnusedOnly, h_eta_real_pls_used_matchedToBoth, h_eta_real_pls_used_matchedToUsedOnly);
    createAndSaveDetailedStack("c_phi_pls_used_detailed", "Phi of Used PLS (Detailed);#phi;PLS Objects", "compare_phi_pls_used_detailed.png", "Used Real PLS Matched To LS:", false, h_phi_fake_pls_used_detailed, h_phi_real_pls_used_unmatched, h_phi_real_pls_used_matchedToUnusedOnly, h_phi_real_pls_used_matchedToBoth, h_phi_real_pls_used_matchedToUsedOnly);
    createAndSaveDetailedStack("c_pt_pls_unused_detailed", "pT of Unused PLS (Detailed);p_{T} [GeV];PLS Objects", "compare_pt_pls_unused_detailed.png", "Unused Real PLS Matched To LS:", true, h_pt_fake_pls_unused_detailed, h_pt_real_pls_unused_unmatched, h_pt_real_pls_unused_matchedToUnusedOnly, h_pt_real_pls_unused_matchedToBoth, h_pt_real_pls_unused_matchedToUsedOnly);
    createAndSaveDetailedStack("c_eta_pls_unused_detailed", "Eta of Unused PLS (Detailed);#eta;PLS Objects", "compare_eta_pls_unused_detailed.png", "Unused Real PLS Matched To LS:", false, h_eta_fake_pls_unused_detailed, h_eta_real_pls_unused_unmatched, h_eta_real_pls_unused_matchedToUnusedOnly, h_eta_real_pls_unused_matchedToBoth, h_eta_real_pls_unused_matchedToUsedOnly);
    createAndSaveDetailedStack("c_phi_pls_unused_detailed", "Phi of Unused PLS (Detailed);#phi;PLS Objects", "compare_phi_pls_unused_detailed.png", "Unused Real PLS Matched To LS:", false, h_phi_fake_pls_unused_detailed, h_phi_real_pls_unused_unmatched, h_phi_real_pls_unused_matchedToUnusedOnly, h_phi_real_pls_unused_matchedToBoth, h_phi_real_pls_unused_matchedToUsedOnly);
    // Delta Plots
    createAndSaveSimplePlot("c_delta_pt", "delta_pt.png", h_delta_pt, "#Delta p_{T} (LS - PLS)", "#Delta p_{T} [GeV]", "pT2 Objects");
    createAndSaveSimplePlot("c_delta_eta", "delta_eta.png", h_delta_eta, "#Delta #eta (LS - PLS)", "#Delta #eta", "pT2 Objects");
    createAndSaveSimplePlot("c_delta_phi", "delta_phi.png", h_delta_phi, "#Delta #phi (LS - PLS)", "#Delta #phi [rad]", "pT2 Objects");
    createAndSaveSimplePlot("c_delta_R", "delta_R.png", h_delta_R, "#Delta R (LS - PLS)", "#Delta R", "pT2 Objects");
    createAndSaveSimplePlot("c_extrapolation_dist_3d", "extrapolation_dist_3d.png", h_extrapolation_dist_3d,"3D Distance between extrapolated PLS and LS", "Distance [cm]", "Ideal pT2 Pairs");

    createAndSaveSimplePlot("c_extrapolation_delta_z", "extrapolation_delta_z.png", h_extrapolation_delta_z,"R-Z Extrapolation #Delta r (LS - PLS)", "#Delta r [cm]", "Ideal pT2 Pairs");
    std::cout << "Finished creating plots!" << std::endl;
    // --- Clean up ---
    inputFile->Close();
}
