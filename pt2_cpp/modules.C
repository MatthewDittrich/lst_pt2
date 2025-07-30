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

// A simple struct to track matching status for a given simIdx
struct PlsMatchStatus {
    int used_count = 0;
    int unused_count = 0;
};

//=========================================================================================
// MODULE 1: Find Used Indices
// Identifies the indices of 'ls' and 'pls' objects that are used in higher-level objects.
//=========================================================================================
void findUsedIndices(
    std::set<int>& used_ls_indices,
    std::set<int>& used_pls_indices,
    const std::vector<int>& tc_pt5Idx_vec, const std::vector<int>& pt5_t5Idx_vec,
    const std::vector<int>& t5_t3Idx0_vec, const std::vector<int>& t5_t3Idx1_vec,
    const std::vector<int>& t3_lsIdx0_vec, const std::vector<int>& t3_lsIdx1_vec,
    const std::vector<int>& tc_pt3Idx_vec, const std::vector<int>& pt3_t3Idx_vec,
    const std::vector<int>& tc_t5Idx_vec, const std::vector<int>& pt5_plsIdx_vec,
    const std::vector<int>& pt3_plsIdx_vec)
{
    // Trace from pt5 track-candidates
    for (int pt5_Idx : tc_pt5Idx_vec) if (pt5_Idx != -999) {
        if (pt5_plsIdx_vec.at(pt5_Idx) != -999) used_pls_indices.insert(pt5_plsIdx_vec.at(pt5_Idx));
        int t5_Idx = pt5_t5Idx_vec.at(pt5_Idx);
        if (t5_Idx != -999) {
            int t3_0 = t5_t3Idx0_vec.at(t5_Idx);
            if (t3_0 != -999) {
                if (t3_lsIdx0_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_0));
                if (t3_lsIdx1_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_0));
            }
            int t3_1 = t5_t3Idx1_vec.at(t5_Idx);
            if (t3_1 != -999) {
                if (t3_lsIdx0_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_1));
                if (t3_lsIdx1_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_1));
            }
        }
    }
    // Trace from pt3 track-candidates
    for (int pt3_Idx : tc_pt3Idx_vec) if (pt3_Idx != -999) {
        if (pt3_plsIdx_vec.at(pt3_Idx) != -999) used_pls_indices.insert(pt3_plsIdx_vec.at(pt3_Idx));
        int t3_idx = pt3_t3Idx_vec.at(pt3_Idx);
        if (t3_idx != -999) {
            if (t3_lsIdx0_vec.at(t3_idx) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_idx));
            if (t3_lsIdx1_vec.at(t3_idx) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_idx));
        }
    }
    // Trace from t5 track-candidates
    for (int t5_Idx : tc_t5Idx_vec) if (t5_Idx != -999) {
        int t3_0 = t5_t3Idx0_vec.at(t5_Idx);
        if (t3_0 != -999) {
            if (t3_lsIdx0_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_0));
            if (t3_lsIdx1_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_0));
        }
        int t3_1 = t5_t3Idx1_vec.at(t5_Idx);
        if (t3_1 != -999) {
            if (t3_lsIdx0_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_1));
            if (t3_lsIdx1_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_1));
        }
    }
}

//=========================================================================================
// MODULE 2: Find PLS Matches
// Creates a map from a simulation index (simIdx) to its use status (used, unused, or both).
//=========================================================================================
void getPlsMatchStatus(
    std::map<int, PlsMatchStatus>& simIdx_to_pls_status,
    const std::vector<int>& pls_isFake_vec,
    const std::vector<int>& pls_simIdx_vec,
    const std::set<int>& used_pls_indices)
{
    for (size_t j = 0; j < pls_simIdx_vec.size(); ++j) {
        if (pls_isFake_vec.at(j) == 0) { // is REAL
            int sim_idx = pls_simIdx_vec.at(j);
            if (sim_idx != -999) {
                if (used_pls_indices.count(j) > 0) {
                    simIdx_to_pls_status[sim_idx].used_count++;
                } else {
                    simIdx_to_pls_status[sim_idx].unused_count++;
                }
            }
        }
    }
}

//=========================================================================================
// MODULE 3: Fill Histograms
// Categorizes a given LS object and fills the corresponding histograms.
//=========================================================================================
void fillHistograms(
    bool is_ls_used, bool is_ls_fake, float pt, float eta, float phi, int ls_sim_idx,
    const std::map<int, PlsMatchStatus>& simIdx_to_pls_status,
    // Histograms for Fake LS
    TH1D* h_pt_fake_used, TH1D* h_pt_fake_unused,
    TH1D* h_eta_fake_used, TH1D* h_eta_fake_unused,
    TH1D* h_phi_fake_used, TH1D* h_phi_fake_unused,
    // Histograms for USED Real LS
    TH1D* h_pt_real_used_unmatched, TH1D* h_pt_real_used_matchedToBoth,
    TH1D* h_pt_real_used_matchedToUsedOnly, TH1D* h_pt_real_used_matchedToUnusedOnly,
    TH1D* h_eta_real_used_unmatched, TH1D* h_eta_real_used_matchedToBoth,
    TH1D* h_eta_real_used_matchedToUsedOnly, TH1D* h_eta_real_used_matchedToUnusedOnly,
    TH1D* h_phi_real_used_unmatched, TH1D* h_phi_real_used_matchedToBoth,
    TH1D* h_phi_real_used_matchedToUsedOnly, TH1D* h_phi_real_used_matchedToUnusedOnly,
    // Histograms for UNUSED Real LS
    TH1D* h_pt_real_unused_unmatched, TH1D* h_pt_real_unused_matchedToBoth,
    TH1D* h_pt_real_unused_matchedToUsedOnly, TH1D* h_pt_real_unused_matchedToUnusedOnly,
    TH1D* h_eta_real_unused_unmatched, TH1D* h_eta_real_unused_matchedToBoth,
    TH1D* h_eta_real_unused_matchedToUsedOnly, TH1D* h_eta_real_unused_matchedToUnusedOnly,
    TH1D* h_phi_real_unused_unmatched, TH1D* h_phi_real_unused_matchedToBoth,
    TH1D* h_phi_real_unused_matchedToUsedOnly, TH1D* h_phi_real_unused_matchedToUnusedOnly)
{
    if (is_ls_fake) {
        if (is_ls_used) { h_pt_fake_used->Fill(pt); h_eta_fake_used->Fill(eta); h_phi_fake_used->Fill(phi); }
        else { h_pt_fake_unused->Fill(pt); h_eta_fake_unused->Fill(eta); h_phi_fake_unused->Fill(phi); }
    } else { // REAL LS
        auto match_it = (ls_sim_idx == -999) ? simIdx_to_pls_status.end() : simIdx_to_pls_status.find(ls_sim_idx);

        if (match_it == simIdx_to_pls_status.end()) { // No Match
            if (is_ls_used) { h_pt_real_used_unmatched->Fill(pt); h_eta_real_used_unmatched->Fill(eta); h_phi_real_used_unmatched->Fill(phi); }
            else { h_pt_real_unused_unmatched->Fill(pt); h_eta_real_unused_unmatched->Fill(eta); h_phi_real_unused_unmatched->Fill(phi); }
        } else {
            const auto& status = match_it->second;
            if (status.used_count > 0 && status.unused_count > 0) { // Matched to BOTH
                if (is_ls_used) { h_pt_real_used_matchedToBoth->Fill(pt); h_eta_real_used_matchedToBoth->Fill(eta); h_phi_real_used_matchedToBoth->Fill(phi); }
                else { h_pt_real_unused_matchedToBoth->Fill(pt); h_eta_real_unused_matchedToBoth->Fill(eta); h_phi_real_unused_matchedToBoth->Fill(phi); }
            } else if (status.used_count > 0) { // Matched to USED Only
                if (is_ls_used) { h_pt_real_used_matchedToUsedOnly->Fill(pt); h_eta_real_used_matchedToUsedOnly->Fill(eta); h_phi_real_used_matchedToUsedOnly->Fill(phi); }
                else { h_pt_real_unused_matchedToUsedOnly->Fill(pt); h_eta_real_unused_matchedToUsedOnly->Fill(eta); h_phi_real_unused_matchedToUsedOnly->Fill(phi); }
            } else { // Matched to UNUSED Only
                if (is_ls_used) { h_pt_real_used_matchedToUnusedOnly->Fill(pt); h_eta_real_used_matchedToUnusedOnly->Fill(eta); h_phi_real_used_matchedToUnusedOnly->Fill(phi); }
                else { h_pt_real_unused_matchedToUnusedOnly->Fill(pt); h_eta_real_unused_matchedToUnusedOnly->Fill(eta); h_phi_real_unused_matchedToUnusedOnly->Fill(phi); }
            }
        }
    }
}

//=========================================================================================
// MODULE 4: Create and Save a Plot
// Creates a styled THStack plot with a legend and saves it to a file.
//=========================================================================================
void createAndSavePlot(
    const std::string& canvasName, const std::string& stackTitle, const std::string& saveFileName,
    const std::string& legendHeader, double legend_x1, double legend_y1, double legend_x2, double legend_y2,
    bool setLogY,
    TH1D* h_cat1, TH1D* h_cat2, TH1D* h_cat3, TH1D* h_cat4, TH1D* h_cat5)
{
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    if (setLogY) canvas->SetLogy();

    THStack *hs = new THStack("hs", stackTitle.c_str());

    // Set colors
    h_cat1->SetFillColor(kGray);      // Fake
    h_cat2->SetFillColor(kRed-7);     // Unmatched
    h_cat3->SetFillColor(kAzure+1);   // Matched to Unused
    h_cat4->SetFillColor(kOrange+1);  // Matched to Both
    h_cat5->SetFillColor(kGreen+2);   // Matched to Used

    // Add in a deliberate order for stacking
    hs->Add(h_cat3); // Matched to Unused
    hs->Add(h_cat5); // Matched to Used
    hs->Add(h_cat2); // Unmatched
    hs->Add(h_cat4); // Matched to Both
    hs->Add(h_cat1); // Fake
    hs->Draw("HIST");

    auto legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2);
    //legend->SetHeader(legendHeader.c_str(), "C");
    legend->AddEntry(h_cat5, "Real, Matched to USED PLS Only", "f");
    legend->AddEntry(h_cat4, "Real, Matched to Both", "f");
    legend->AddEntry(h_cat3, "Real, Matched to UNUSED PLS Only", "f");
    legend->AddEntry(h_cat2, "Real, No PLS Match", "f");
    legend->AddEntry(h_cat1, "Fake", "f");
    legend->SetBorderSize(0);
    legend->Draw();

    canvas->SaveAs(saveFileName.c_str());
    delete canvas; // Clean up memory
}
