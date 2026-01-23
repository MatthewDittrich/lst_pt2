#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <vector>
#include <map>
#include <string>
#include <iostream>

// ROOT Includes
#include <TH1D.h>

// Project Includes
#include "Structs.h" // For PlsMatchStatus and LsMatchStatus

class HistManager {
public:
    HistManager();
    ~HistManager();

    void Init();

    // =============================================================
    //                    HISTOGRAM POINTERS
    // =============================================================

    // -------------------------------------------------------------
    // 1. Simple LS Histograms (Real/Fake, Used/Unused)
    // -------------------------------------------------------------
    TH1D *h_pt_real_ls_used,   *h_pt_fake_ls_used;
    TH1D *h_pt_real_ls_unused, *h_pt_fake_ls_unused;
    
    TH1D *h_eta_real_ls_used,   *h_eta_fake_ls_used;
    TH1D *h_eta_real_ls_unused, *h_eta_fake_ls_unused;
    
    TH1D *h_phi_real_ls_used,   *h_phi_fake_ls_used;
    TH1D *h_phi_real_ls_unused, *h_phi_fake_ls_unused;

    // -------------------------------------------------------------
    // 2. Simple PLS Histograms (Real/Fake, Used/Unused)
    // -------------------------------------------------------------
    TH1D *h_pt_real_pls_used,   *h_pt_fake_pls_used;
    TH1D *h_pt_real_pls_unused, *h_pt_fake_pls_unused;
    
    TH1D *h_eta_real_pls_used,   *h_eta_fake_pls_used;
    TH1D *h_eta_real_pls_unused, *h_eta_fake_pls_unused;
    
    TH1D *h_phi_real_pls_used,   *h_phi_fake_pls_used;
    TH1D *h_phi_real_pls_unused, *h_phi_fake_pls_unused;

    // -------------------------------------------------------------
    // 3. Detailed LS Histograms (Matching Categories)
    // -------------------------------------------------------------
    // Fakes (Detailed)
    TH1D *h_pt_fake_ls_used_detailed,  *h_pt_fake_ls_unused_detailed;
    TH1D *h_eta_fake_ls_used_detailed, *h_eta_fake_ls_unused_detailed;
    TH1D *h_phi_fake_ls_used_detailed, *h_phi_fake_ls_unused_detailed;

    // Real Used (Categorized by match)
    TH1D *h_pt_real_ls_used_unmatched,  *h_pt_real_ls_used_matchedToBoth,  *h_pt_real_ls_used_matchedToUsedOnly,  *h_pt_real_ls_used_matchedToUnusedOnly;
    TH1D *h_eta_real_ls_used_unmatched, *h_eta_real_ls_used_matchedToBoth, *h_eta_real_ls_used_matchedToUsedOnly, *h_eta_real_ls_used_matchedToUnusedOnly;
    TH1D *h_phi_real_ls_used_unmatched, *h_phi_real_ls_used_matchedToBoth, *h_phi_real_ls_used_matchedToUsedOnly, *h_phi_real_ls_used_matchedToUnusedOnly;

    // Real Unused (Categorized by match)
    TH1D *h_pt_real_ls_unused_unmatched,  *h_pt_real_ls_unused_matchedToBoth,  *h_pt_real_ls_unused_matchedToUsedOnly,  *h_pt_real_ls_unused_matchedToUnusedOnly;
    TH1D *h_eta_real_ls_unused_unmatched, *h_eta_real_ls_unused_matchedToBoth, *h_eta_real_ls_unused_matchedToUsedOnly, *h_eta_real_ls_unused_matchedToUnusedOnly;
    TH1D *h_phi_real_ls_unused_unmatched, *h_phi_real_ls_unused_matchedToBoth, *h_phi_real_ls_unused_matchedToUsedOnly, *h_phi_real_ls_unused_matchedToUnusedOnly;

    // -------------------------------------------------------------
    // 4. Detailed PLS Histograms (Matching Categories)
    // -------------------------------------------------------------
    // Fakes (Detailed)
    TH1D *h_pt_fake_pls_used_detailed,  *h_pt_fake_pls_unused_detailed;
    TH1D *h_eta_fake_pls_used_detailed, *h_eta_fake_pls_unused_detailed;
    TH1D *h_phi_fake_pls_used_detailed, *h_phi_fake_pls_unused_detailed;

    // Real Used (Categorized by match)
    TH1D *h_pt_real_pls_used_unmatched,  *h_pt_real_pls_used_matchedToBoth,  *h_pt_real_pls_used_matchedToUsedOnly,  *h_pt_real_pls_used_matchedToUnusedOnly;
    TH1D *h_eta_real_pls_used_unmatched, *h_eta_real_pls_used_matchedToBoth, *h_eta_real_pls_used_matchedToUsedOnly, *h_eta_real_pls_used_matchedToUnusedOnly;
    TH1D *h_phi_real_pls_used_unmatched, *h_phi_real_pls_used_matchedToBoth, *h_phi_real_pls_used_matchedToUsedOnly, *h_phi_real_pls_used_matchedToUnusedOnly;

    // Real Unused (Categorized by match)
    TH1D *h_pt_real_pls_unused_unmatched,  *h_pt_real_pls_unused_matchedToBoth,  *h_pt_real_pls_unused_matchedToUsedOnly,  *h_pt_real_pls_unused_matchedToUnusedOnly;
    TH1D *h_eta_real_pls_unused_unmatched, *h_eta_real_pls_unused_matchedToBoth, *h_eta_real_pls_unused_matchedToUsedOnly, *h_eta_real_pls_unused_matchedToUnusedOnly;
    TH1D *h_phi_real_pls_unused_unmatched, *h_phi_real_pls_unused_matchedToBoth, *h_phi_real_pls_unused_matchedToUsedOnly, *h_phi_real_pls_unused_matchedToUnusedOnly;

    // -------------------------------------------------------------
    // 5. Delta Objects (pT2 and Candidates)
    // -------------------------------------------------------------
    // Ideal pT2
    TH1D *h_delta_pt, *h_delta_eta, *h_delta_phi, *h_delta_R;
    
    // Candidates
    TH1D *h_delta_pt_cand, *h_delta_eta_cand, *h_delta_phi_cand, *h_delta_R_cand;

    // -------------------------------------------------------------
    // 6. Extrapolation and Layer Logic
    // -------------------------------------------------------------
    TH1D *h_extrapolation_dist_3d;
    TH1D *h_extrapolation_delta_z;
    TH1D *h_extrapolation_delta_r_reverse_combined;

    // Layer-specific Vectors
    std::vector<TH1D*> h_delta_r_by_layer;
    std::vector<TH1D*> h_dist_3d_by_layer;
    std::vector<TH1D*> h_extrapolation_delta_r_reverse_by_hit;

    // =============================================================
    //                    FILLING METHODS
    // =============================================================

    /**
     * Fills the Simple Real/Fake Used/Unused histograms for either LS or PLS.
     */
    void FillSimpleHistograms(
        bool is_used, bool is_fake, 
        float pt, float eta, float phi, 
        TH1D* h_pt_r_u,   TH1D* h_pt_f_u,   TH1D* h_pt_r_un,   TH1D* h_pt_f_un, 
        TH1D* h_eta_r_u,  TH1D* h_eta_f_u,  TH1D* h_eta_r_un,  TH1D* h_eta_f_un, 
        TH1D* h_phi_r_u,  TH1D* h_phi_f_u,  TH1D* h_phi_r_un,  TH1D* h_phi_f_un
    );

    /**
     * Fills the Detailed LS histograms based on matching status.
     */
    void FillDetailedLsHistograms(
        bool is_ls_used, bool is_ls_fake, 
        float pt, float eta, float phi, 
        int ls_sim_idx, 
        const std::map<int, PlsMatchStatus>& simIdx_to_pls_status
    );

    /**
     * Fills the Detailed PLS histograms based on matching status.
     */
    void FillDetailedPlsHistograms(
        bool is_pls_used, bool is_pls_fake, 
        float pt, float eta, float phi, 
        int pls_sim_idx, 
        const std::map<int, LsMatchStatus>& simIdx_to_ls_status
    );
    
    // Helper to fill candidate histograms
    void FillCandHistos(
        const std::vector<int>& pls_idx_vec, const std::vector<int>& ls_idx_vec,
        const std::vector<float>& pls_pt, const std::vector<float>& ls_pt,
        const std::vector<float>& pls_eta, const std::vector<float>& ls_eta,
        const std::vector<float>& pls_phi, const std::vector<float>& ls_phi
    );

private:
    // Helper to keep track of all allocated histograms for easy deletion in destructor
    std::vector<TH1D*> all_hists; 
};

#endif // HISTMANAGER_H
