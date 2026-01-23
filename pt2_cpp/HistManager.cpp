#include "HistManager.h"
#include <TString.h>
#include <TVector2.h>
#include <cmath>

HistManager::HistManager() {
    // Constructor
}

HistManager::~HistManager() {
    // Clean up all allocated histograms
    for (TH1D* h : all_hists) {
        if (h) delete h;
    }
    all_hists.clear();
}

void HistManager::Init() {
    // --- Configuration ---
    int pt_bins = 180; double pt_max = 250.0;
    int eta_bins = 180; double eta_max = 4.5;
    int phi_bins = 180; double phi_max = 3.15;

    // Helper lambda to reduce code repetition
    auto AddHist = [&](const char* name, const char* title, int bins, double min, double max) -> TH1D* {
        TH1D* h = new TH1D(name, title, bins, min, max);
        all_hists.push_back(h);
        return h;
    };

    // =============================================================
    // 1. Simple LS Histograms
    // =============================================================
    h_pt_real_ls_used = AddHist("h_pt_real_ls_used", "", pt_bins, 0, pt_max);
    h_pt_fake_ls_used = AddHist("h_pt_fake_ls_used", "", pt_bins, 0, pt_max);
    h_pt_real_ls_unused = AddHist("h_pt_real_ls_unused", "", pt_bins, 0, pt_max);
    h_pt_fake_ls_unused = AddHist("h_pt_fake_ls_unused", "", pt_bins, 0, pt_max);

    h_eta_real_ls_used = AddHist("h_eta_real_ls_used", "", eta_bins, -eta_max, eta_max);
    h_eta_fake_ls_used = AddHist("h_eta_fake_ls_used", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_unused = AddHist("h_eta_real_ls_unused", "", eta_bins, -eta_max, eta_max);
    h_eta_fake_ls_unused = AddHist("h_eta_fake_ls_unused", "", eta_bins, -eta_max, eta_max);

    h_phi_real_ls_used = AddHist("h_phi_real_ls_used", "", phi_bins, -phi_max, phi_max);
    h_phi_fake_ls_used = AddHist("h_phi_fake_ls_used", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_unused = AddHist("h_phi_real_ls_unused", "", phi_bins, -phi_max, phi_max);
    h_phi_fake_ls_unused = AddHist("h_phi_fake_ls_unused", "", phi_bins, -phi_max, phi_max);

    // =============================================================
    // 2. Simple PLS Histograms
    // =============================================================
    h_pt_real_pls_used = AddHist("h_pt_real_pls_used", "", pt_bins, 0, pt_max);
    h_pt_fake_pls_used = AddHist("h_pt_fake_pls_used", "", pt_bins, 0, pt_max);
    h_pt_real_pls_unused = AddHist("h_pt_real_pls_unused", "", pt_bins, 0, pt_max);
    h_pt_fake_pls_unused = AddHist("h_pt_fake_pls_unused", "", pt_bins, 0, pt_max);

    h_eta_real_pls_used = AddHist("h_eta_real_pls_used", "", eta_bins, -eta_max, eta_max);
    h_eta_fake_pls_used = AddHist("h_eta_fake_pls_used", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_unused = AddHist("h_eta_real_pls_unused", "", eta_bins, -eta_max, eta_max);
    h_eta_fake_pls_unused = AddHist("h_eta_fake_pls_unused", "", eta_bins, -eta_max, eta_max);

    h_phi_real_pls_used = AddHist("h_phi_real_pls_used", "", phi_bins, -phi_max, phi_max);
    h_phi_fake_pls_used = AddHist("h_phi_fake_pls_used", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_unused = AddHist("h_phi_real_pls_unused", "", phi_bins, -phi_max, phi_max);
    h_phi_fake_pls_unused = AddHist("h_phi_fake_pls_unused", "", phi_bins, -phi_max, phi_max);

    // =============================================================
    // 3. Detailed LS Histograms
    // =============================================================
    h_pt_fake_ls_used_detailed = AddHist("h_pt_fake_ls_used_detailed", "", pt_bins, 0, pt_max);
    h_pt_fake_ls_unused_detailed = AddHist("h_pt_fake_ls_unused_detailed", "", pt_bins, 0, pt_max);
    h_eta_fake_ls_used_detailed = AddHist("h_eta_fake_ls_used_detailed", "", eta_bins, -eta_max, eta_max);
    h_eta_fake_ls_unused_detailed = AddHist("h_eta_fake_ls_unused_detailed", "", eta_bins, -eta_max, eta_max);
    h_phi_fake_ls_used_detailed = AddHist("h_phi_fake_ls_used_detailed", "", phi_bins, -phi_max, phi_max);
    h_phi_fake_ls_unused_detailed = AddHist("h_phi_fake_ls_unused_detailed", "", phi_bins, -phi_max, phi_max);

    // Real LS Match Categories
    h_pt_real_ls_used_unmatched = AddHist("h_pt_real_ls_used_unmatched", "", pt_bins, 0, pt_max);
    h_pt_real_ls_used_matchedToBoth = AddHist("h_pt_real_ls_used_matchedToBoth", "", pt_bins, 0, pt_max);
    h_pt_real_ls_used_matchedToUsedOnly = AddHist("h_pt_real_ls_used_matchedToUsedOnly", "", pt_bins, 0, pt_max);
    h_pt_real_ls_used_matchedToUnusedOnly = AddHist("h_pt_real_ls_used_matchedToUnusedOnly", "", pt_bins, 0, pt_max);

    h_eta_real_ls_used_unmatched = AddHist("h_eta_real_ls_used_unmatched", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_used_matchedToBoth = AddHist("h_eta_real_ls_used_matchedToBoth", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_used_matchedToUsedOnly = AddHist("h_eta_real_ls_used_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_used_matchedToUnusedOnly = AddHist("h_eta_real_ls_used_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max);

    h_phi_real_ls_used_unmatched = AddHist("h_phi_real_ls_used_unmatched", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_used_matchedToBoth = AddHist("h_phi_real_ls_used_matchedToBoth", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_used_matchedToUsedOnly = AddHist("h_phi_real_ls_used_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_used_matchedToUnusedOnly = AddHist("h_phi_real_ls_used_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);

    h_pt_real_ls_unused_unmatched = AddHist("h_pt_real_ls_unused_unmatched", "", pt_bins, 0, pt_max);
    h_pt_real_ls_unused_matchedToBoth = AddHist("h_pt_real_ls_unused_matchedToBoth", "", pt_bins, 0, pt_max);
    h_pt_real_ls_unused_matchedToUsedOnly = AddHist("h_pt_real_ls_unused_matchedToUsedOnly", "", pt_bins, 0, pt_max);
    h_pt_real_ls_unused_matchedToUnusedOnly = AddHist("h_pt_real_ls_unused_matchedToUnusedOnly", "", pt_bins, 0, pt_max);

    h_eta_real_ls_unused_unmatched = AddHist("h_eta_real_ls_unused_unmatched", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_unused_matchedToBoth = AddHist("h_eta_real_ls_unused_matchedToBoth", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_unused_matchedToUsedOnly = AddHist("h_eta_real_ls_unused_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max);
    h_eta_real_ls_unused_matchedToUnusedOnly = AddHist("h_eta_real_ls_unused_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max);

    h_phi_real_ls_unused_unmatched = AddHist("h_phi_real_ls_unused_unmatched", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_unused_matchedToBoth = AddHist("h_phi_real_ls_unused_matchedToBoth", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_unused_matchedToUsedOnly = AddHist("h_phi_real_ls_unused_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max);
    h_phi_real_ls_unused_matchedToUnusedOnly = AddHist("h_phi_real_ls_unused_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);

    // =============================================================
    // 4. Detailed PLS Histograms
    // =============================================================
    h_pt_fake_pls_used_detailed = AddHist("h_pt_fake_pls_used_detailed", "", pt_bins, 0, pt_max);
    h_pt_fake_pls_unused_detailed = AddHist("h_pt_fake_pls_unused_detailed", "", pt_bins, 0, pt_max);
    h_eta_fake_pls_used_detailed = AddHist("h_eta_fake_pls_used_detailed", "", eta_bins, -eta_max, eta_max);
    h_eta_fake_pls_unused_detailed = AddHist("h_eta_fake_pls_unused_detailed", "", eta_bins, -eta_max, eta_max);
    h_phi_fake_pls_used_detailed = AddHist("h_phi_fake_pls_used_detailed", "", phi_bins, -phi_max, phi_max);
    h_phi_fake_pls_unused_detailed = AddHist("h_phi_fake_pls_unused_detailed", "", phi_bins, -phi_max, phi_max);

    // Real PLS Match Categories
    h_pt_real_pls_used_unmatched = AddHist("h_pt_real_pls_used_unmatched", "", pt_bins, 0, pt_max);
    h_pt_real_pls_used_matchedToBoth = AddHist("h_pt_real_pls_used_matchedToBoth", "", pt_bins, 0, pt_max);
    h_pt_real_pls_used_matchedToUsedOnly = AddHist("h_pt_real_pls_used_matchedToUsedOnly", "", pt_bins, 0, pt_max);
    h_pt_real_pls_used_matchedToUnusedOnly = AddHist("h_pt_real_pls_used_matchedToUnusedOnly", "", pt_bins, 0, pt_max);

    h_eta_real_pls_used_unmatched = AddHist("h_eta_real_pls_used_unmatched", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_used_matchedToBoth = AddHist("h_eta_real_pls_used_matchedToBoth", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_used_matchedToUsedOnly = AddHist("h_eta_real_pls_used_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_used_matchedToUnusedOnly = AddHist("h_eta_real_pls_used_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max);

    h_phi_real_pls_used_unmatched = AddHist("h_phi_real_pls_used_unmatched", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_used_matchedToBoth = AddHist("h_phi_real_pls_used_matchedToBoth", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_used_matchedToUsedOnly = AddHist("h_phi_real_pls_used_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_used_matchedToUnusedOnly = AddHist("h_phi_real_pls_used_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);

    h_pt_real_pls_unused_unmatched = AddHist("h_pt_real_pls_unused_unmatched", "", pt_bins, 0, pt_max);
    h_pt_real_pls_unused_matchedToBoth = AddHist("h_pt_real_pls_unused_matchedToBoth", "", pt_bins, 0, pt_max);
    h_pt_real_pls_unused_matchedToUsedOnly = AddHist("h_pt_real_pls_unused_matchedToUsedOnly", "", pt_bins, 0, pt_max);
    h_pt_real_pls_unused_matchedToUnusedOnly = AddHist("h_pt_real_pls_unused_matchedToUnusedOnly", "", pt_bins, 0, pt_max);

    h_eta_real_pls_unused_unmatched = AddHist("h_eta_real_pls_unused_unmatched", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_unused_matchedToBoth = AddHist("h_eta_real_pls_unused_matchedToBoth", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_unused_matchedToUsedOnly = AddHist("h_eta_real_pls_unused_matchedToUsedOnly", "", eta_bins, -eta_max, eta_max);
    h_eta_real_pls_unused_matchedToUnusedOnly = AddHist("h_eta_real_pls_unused_matchedToUnusedOnly", "", eta_bins, -eta_max, eta_max);

    h_phi_real_pls_unused_unmatched = AddHist("h_phi_real_pls_unused_unmatched", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_unused_matchedToBoth = AddHist("h_phi_real_pls_unused_matchedToBoth", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_unused_matchedToUsedOnly = AddHist("h_phi_real_pls_unused_matchedToUsedOnly", "", phi_bins, -phi_max, phi_max);
    h_phi_real_pls_unused_matchedToUnusedOnly = AddHist("h_phi_real_pls_unused_matchedToUnusedOnly", "", phi_bins, -phi_max, phi_max);

    // =============================================================
    // 5. Delta Histograms
    // =============================================================
    h_delta_pt = AddHist("h_delta_pt", "#Delta p_{T} (LS - PLS);#Delta p_{T} [GeV];pT2 Objects", 180, -10, 10);
    h_delta_eta = AddHist("h_delta_eta", "#Delta #eta (LS - PLS);#Delta #eta;pT2 Objects", 180, -1, 1);
    h_delta_phi = AddHist("h_delta_phi", "#Delta #phi (LS - PLS);#Delta #phi [rad];pT2 Objects", 180, -3, 3);
    h_delta_R = AddHist("h_delta_R", "#Delta R (LS - PLS);#Delta R;pT2 Objects", 180, 0, 3);

    h_delta_pt_cand = AddHist("h_delta_pt_cand", "#Delta p_{T} (LS - PLS);#Delta p_{T} [GeV];pT2 Objects", 180, -10, 10);
    h_delta_eta_cand = AddHist("h_delta_eta_cand", "#Delta #eta (LS - PLS);#Delta #eta;pT2 Objects", 180, -1, 1);
    h_delta_phi_cand = AddHist("h_delta_phi_cand", "#Delta #phi (LS - PLS);#Delta #phi [rad];pT2 Objects", 180, -3, 3);
    h_delta_R_cand = AddHist("h_delta_R_cand", "#Delta R (LS - PLS);#Delta R;pT2 Objects", 180, 0, 3);

    h_extrapolation_dist_3d = AddHist("h_extrapolation_dist_3d","3D Distance between extrapolated PLS and LS;Distance [cm];Ideal pT2 Pairs", 180, 0, 100.0);
    h_extrapolation_delta_z = AddHist("h_extrapolation_delta_z", "R-Z Extrapolation #Delta r (LS - PLS);#Delta r [cm];Ideal pT2 Pairs", 180, -15.0, 15.0);
    h_extrapolation_delta_r_reverse_combined = AddHist("h_extrapolation_delta_r_reverse_combined", "Inclusive R-Z Extrapolation #Delta r (LS -> All PLS Hits);#Delta r [cm];Entries", 180, -10.0, 10.0);

    // =============================================================
    // 6. Layer/Hit Arrays
    // =============================================================
    int num_layers = 12;
    h_delta_r_by_layer.reserve(num_layers);
    h_dist_3d_by_layer.reserve(num_layers);

    for (int i = 0; i < num_layers; ++i) {
        TString name_dr = TString::Format("h_delta_r_layer_%d", i);
        TString title_dr = TString::Format("Z-R Extrapolation #Delta r (Layer %d);#Delta r [cm];Entries", i);
        h_delta_r_by_layer.push_back(AddHist(name_dr.Data(), title_dr.Data(), 180, -25.0, 25.0));

        TString name_d3 = TString::Format("h_dist_3d_layer_%d", i);
        TString title_d3 = TString::Format("3D Distance (Layer %d);Distance [cm];Entries", i);
        h_dist_3d_by_layer.push_back(AddHist(name_d3.Data(), title_d3.Data(), 180, 0, 100.0));
    }

    int n_pls_hits_max = 4; 
    h_extrapolation_delta_r_reverse_by_hit.reserve(n_pls_hits_max);
    for (int i = 0; i < n_pls_hits_max; ++i) {
        TString name = TString::Format("h_extrapolation_delta_r_reverse_hit%d", i);
        TString title = TString::Format("R-Z Extrapolation #Delta r (LS -> PLS Hit %d);#Delta r [cm];Ideal pT2 Pairs", i);
        h_extrapolation_delta_r_reverse_by_hit.push_back(AddHist(name.Data(), title.Data(), 180, -10.0, 10.0));
    }
}

void HistManager::FillSimpleHistograms(
    bool is_used, bool is_fake, 
    float pt, float eta, float phi, 
    TH1D* h_pt_r_u,   TH1D* h_pt_f_u,   TH1D* h_pt_r_un,   TH1D* h_pt_f_un, 
    TH1D* h_eta_r_u,  TH1D* h_eta_f_u,  TH1D* h_eta_r_un,  TH1D* h_eta_f_un, 
    TH1D* h_phi_r_u,  TH1D* h_phi_f_u,  TH1D* h_phi_r_un,  TH1D* h_phi_f_un
) {
    if (is_fake) {
        if (is_used) {
            h_pt_f_u->Fill(pt);  h_eta_f_u->Fill(eta);  h_phi_f_u->Fill(phi);
        } else {
            h_pt_f_un->Fill(pt); h_eta_f_un->Fill(eta); h_phi_f_un->Fill(phi);
        }
    } else {
        if (is_used) {
            h_pt_r_u->Fill(pt);  h_eta_r_u->Fill(eta);  h_phi_r_u->Fill(phi);
        } else {
            h_pt_r_un->Fill(pt); h_eta_r_un->Fill(eta); h_phi_r_un->Fill(phi);
        }
    }
}

void HistManager::FillDetailedLsHistograms(
    bool is_ls_used, bool is_ls_fake, 
    float pt, float eta, float phi, 
    int ls_sim_idx, 
    const std::map<int, PlsMatchStatus>& simIdx_to_pls_status
) {
    // Helper lambda to fill a set of (pt, eta, phi) histograms
    auto fill_all = [&](TH1D* h_pt, TH1D* h_eta, TH1D* h_phi) {
        h_pt->Fill(pt); h_eta->Fill(eta); h_phi->Fill(phi);
    };

    if (is_ls_fake) {
        if (is_ls_used) fill_all(h_pt_fake_ls_used_detailed, h_eta_fake_ls_used_detailed, h_phi_fake_ls_used_detailed);
        else            fill_all(h_pt_fake_ls_unused_detailed, h_eta_fake_ls_unused_detailed, h_phi_fake_ls_unused_detailed);
        return;
    }

    // Real LS Logic
    auto match_it = (ls_sim_idx == -999) ? simIdx_to_pls_status.end() : simIdx_to_pls_status.find(ls_sim_idx);

    if (match_it == simIdx_to_pls_status.end()) {
        // No Match
        if (is_ls_used) fill_all(h_pt_real_ls_used_unmatched, h_eta_real_ls_used_unmatched, h_phi_real_ls_used_unmatched);
        else            fill_all(h_pt_real_ls_unused_unmatched, h_eta_real_ls_unused_unmatched, h_phi_real_ls_unused_unmatched);
    } 
    else {
        const auto& status = match_it->second;
        if (status.used_count > 0 && status.unused_count > 0) {
            // Matched to BOTH
            if (is_ls_used) fill_all(h_pt_real_ls_used_matchedToBoth, h_eta_real_ls_used_matchedToBoth, h_phi_real_ls_used_matchedToBoth);
            else            fill_all(h_pt_real_ls_unused_matchedToBoth, h_eta_real_ls_unused_matchedToBoth, h_phi_real_ls_unused_matchedToBoth);
        }
        else if (status.used_count > 0) {
            // Matched to USED Only
            if (is_ls_used) fill_all(h_pt_real_ls_used_matchedToUsedOnly, h_eta_real_ls_used_matchedToUsedOnly, h_phi_real_ls_used_matchedToUsedOnly);
            else            fill_all(h_pt_real_ls_unused_matchedToUsedOnly, h_eta_real_ls_unused_matchedToUsedOnly, h_phi_real_ls_unused_matchedToUsedOnly);
        }
        else { 
            // Matched to UNUSED Only
            if (is_ls_used) fill_all(h_pt_real_ls_used_matchedToUnusedOnly, h_eta_real_ls_used_matchedToUnusedOnly, h_phi_real_ls_used_matchedToUnusedOnly);
            else            fill_all(h_pt_real_ls_unused_matchedToUnusedOnly, h_eta_real_ls_unused_matchedToUnusedOnly, h_phi_real_ls_unused_matchedToUnusedOnly);
        }
    }
}

void HistManager::FillDetailedPlsHistograms(
    bool is_pls_used, bool is_pls_fake, 
    float pt, float eta, float phi, 
    int pls_sim_idx, 
    const std::map<int, LsMatchStatus>& simIdx_to_ls_status
) {
    auto fill_all = [&](TH1D* h_pt, TH1D* h_eta, TH1D* h_phi) {
        h_pt->Fill(pt); h_eta->Fill(eta); h_phi->Fill(phi);
    };

    if (is_pls_fake) {
        if (is_pls_used) fill_all(h_pt_fake_pls_used_detailed, h_eta_fake_pls_used_detailed, h_phi_fake_pls_used_detailed);
        else             fill_all(h_pt_fake_pls_unused_detailed, h_eta_fake_pls_unused_detailed, h_phi_fake_pls_unused_detailed);
        return;
    }

    // Real PLS Logic
    auto match_it = (pls_sim_idx == -999) ? simIdx_to_ls_status.end() : simIdx_to_ls_status.find(pls_sim_idx);

    if (match_it == simIdx_to_ls_status.end()) {
        if (is_pls_used) fill_all(h_pt_real_pls_used_unmatched, h_eta_real_pls_used_unmatched, h_phi_real_pls_used_unmatched);
        else             fill_all(h_pt_real_pls_unused_unmatched, h_eta_real_pls_unused_unmatched, h_phi_real_pls_unused_unmatched);
    } 
    else {
        const auto& status = match_it->second;
        if (status.used_count > 0 && status.unused_count > 0) {
            if (is_pls_used) fill_all(h_pt_real_pls_used_matchedToBoth, h_eta_real_pls_used_matchedToBoth, h_phi_real_pls_used_matchedToBoth);
            else             fill_all(h_pt_real_pls_unused_matchedToBoth, h_eta_real_pls_unused_matchedToBoth, h_phi_real_pls_unused_matchedToBoth);
        }
        else if (status.used_count > 0) {
            if (is_pls_used) fill_all(h_pt_real_pls_used_matchedToUsedOnly, h_eta_real_pls_used_matchedToUsedOnly, h_phi_real_pls_used_matchedToUsedOnly);
            else             fill_all(h_pt_real_pls_unused_matchedToUsedOnly, h_eta_real_pls_unused_matchedToUsedOnly, h_phi_real_pls_unused_matchedToUsedOnly);
        }
        else { 
            if (is_pls_used) fill_all(h_pt_real_pls_used_matchedToUnusedOnly, h_eta_real_pls_used_matchedToUnusedOnly, h_phi_real_pls_used_matchedToUnusedOnly);
            else             fill_all(h_pt_real_pls_unused_matchedToUnusedOnly, h_eta_real_pls_unused_matchedToUnusedOnly, h_phi_real_pls_unused_matchedToUnusedOnly);
        }
    }
}

void HistManager::FillCandHistos(
    const std::vector<int>& pls_idx_vec, const std::vector<int>& ls_idx_vec,
    const std::vector<float>& pls_pt, const std::vector<float>& ls_pt,
    const std::vector<float>& pls_eta, const std::vector<float>& ls_eta,
    const std::vector<float>& pls_phi, const std::vector<float>& ls_phi
) {
    if (pls_idx_vec.size() != ls_idx_vec.size()) return;

    for (size_t i = 0; i < pls_idx_vec.size(); ++i) {
        int pls_idx = pls_idx_vec[i];
        int ls_idx = ls_idx_vec[i];

        // Safety check for indices
        if (pls_idx < 0 || pls_idx >= pls_pt.size() || ls_idx < 0 || ls_idx >= ls_pt.size()) continue;

        float d_pt = pls_pt[pls_idx] - ls_pt[ls_idx];
        float d_eta = pls_eta[pls_idx] - ls_eta[ls_idx];
        float d_phi = TVector2::Phi_mpi_pi(pls_phi[pls_idx] - ls_phi[ls_idx]);
        float d_R = std::sqrt(d_eta * d_eta + d_phi * d_phi);

        h_delta_pt_cand->Fill(d_pt);
        h_delta_eta_cand->Fill(d_eta);
        h_delta_phi_cand->Fill(d_phi);
        h_delta_R_cand->Fill(d_R);
    }
}
