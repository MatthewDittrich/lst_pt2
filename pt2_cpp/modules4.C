#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TVector2.h>
#include <TVector3.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <TGraph.h>
#include <TF1.h>
#include <utility>

// Structs
struct PlsMatchStatus { int used_count = 0; int unused_count = 0; };
struct LsMatchStatus  { int used_count = 0; int unused_count = 0; };

class pT2Object {
public:
    float delta_pt, delta_eta, delta_phi, delta_R;
    pT2Object(int p_idx, int l_idx, const std::vector<float>& p_pt, const std::vector<float>& p_eta, const std::vector<float>& p_phi, const std::vector<float>& l_pt, const std::vector<float>& l_eta, const std::vector<float>& l_phi) {
        float pls_pt = p_pt.at(p_idx), pls_eta = p_eta.at(p_idx), pls_phi = p_phi.at(p_idx);
        float ls_pt = l_pt.at(l_idx), ls_eta = l_eta.at(l_idx), ls_phi = l_phi.at(l_idx);
        delta_pt = ls_pt - pls_pt; // Corrected order to match plots: LS - PLS
        delta_eta = ls_eta - pls_eta;
        delta_phi = TVector2::Phi_mpi_pi(ls_phi - pls_phi);
        delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
    }
};





double extrapolatePlsAndGetDistance(
    int pls_idx, int ls_idx,
    // PLS hit position branches
    const std::vector<float>& pls_hit0_x, const std::vector<float>& pls_hit0_y, const std::vector<float>& pls_hit0_z,
    const std::vector<float>& pls_hit1_x, const std::vector<float>& pls_hit1_y, const std::vector<float>& pls_hit1_z,
    const std::vector<float>& pls_hit2_x, const std::vector<float>& pls_hit2_y, const std::vector<float>& pls_hit2_z,
    const std::vector<float>& pls_hit3_x, const std::vector<float>& pls_hit3_y, const std::vector<float>& pls_hit3_z,
    // LS -> MD link branches
    const std::vector<int>& ls_mdIdx0, const std::vector<int>& ls_mdIdx1,
    // MD position branches
    const std::vector<float>& md_anchor_x, const std::vector<float>& md_anchor_y, const std::vector<float>& md_anchor_z,
    const std::vector<float>& md_other_x, const std::vector<float>& md_other_y, const std::vector<float>& md_other_z
) {
    // --- Part 1: Define the PLS trajectory as a line ---
    std::vector<TVector3> pls_hits;
    // A special value to check for invalid hits
    const float invalid_hit_val = -999.0; 

    if (pls_hit0_x.at(pls_idx) > invalid_hit_val) pls_hits.emplace_back(pls_hit0_x.at(pls_idx), pls_hit0_y.at(pls_idx), pls_hit0_z.at(pls_idx));
    if (pls_hit1_x.at(pls_idx) > invalid_hit_val) pls_hits.emplace_back(pls_hit1_x.at(pls_idx), pls_hit1_y.at(pls_idx), pls_hit1_z.at(pls_idx));
    if (pls_hit2_x.at(pls_idx) > invalid_hit_val) pls_hits.emplace_back(pls_hit2_x.at(pls_idx), pls_hit2_y.at(pls_idx), pls_hit2_z.at(pls_idx));
    if (pls_hit3_x.at(pls_idx) > invalid_hit_val) pls_hits.emplace_back(pls_hit3_x.at(pls_idx), pls_hit3_y.at(pls_idx), pls_hit3_z.at(pls_idx));

    // We need at least two valid hits to define a line
    if (pls_hits.size() < 2) return -1.0;

    // Use the first and last valid hits to define the line for maximum lever arm
    TVector3 p0 = pls_hits.front(); // Point on the line (inner-most hit)
    TVector3 v = pls_hits.back() - pls_hits.front(); // Direction vector of the line

    // --- Part 2: Calculate the position of the LS ---
    int md_idx0 = ls_mdIdx0.at(ls_idx);
    int md_idx1 = ls_mdIdx1.at(ls_idx);

    // Check if the minidoublet indices are valid
    if (md_idx0 < 0 || md_idx1 < 0) return -1.0;

    // The position of a minidoublet is the average of its anchor and other points
    TVector3 md0_anchor(md_anchor_x.at(md_idx0), md_anchor_y.at(md_idx0), md_anchor_z.at(md_idx0));
    TVector3 md0_other(md_other_x.at(md_idx0), md_other_y.at(md_idx0), md_other_z.at(md_idx0));
    TVector3 md0_pos = (md0_anchor + md0_other) * 0.5;

    TVector3 md1_anchor(md_anchor_x.at(md_idx1), md_anchor_y.at(md_idx1), md_anchor_z.at(md_idx1));
    TVector3 md1_other(md_other_x.at(md_idx1), md_other_y.at(md_idx1), md_other_z.at(md_idx1));
    TVector3 md1_pos = (md1_anchor + md1_other) * 0.5;

    // The position of the LS is the average of its two minidoublet positions
    TVector3 ls_pos = (md0_pos + md1_pos) * 0.5;

    // --- Part 3: Calculate the distance of closest approach ---
    // This is a standard vector projection problem.
    TVector3 vec_to_ls = ls_pos - p0;
    // Project vec_to_ls onto the direction vector v to find the closest point on the line
    double t = vec_to_ls.Dot(v) / v.Mag2();
    TVector3 closest_point = p0 + t * v;

    // The distance is the magnitude of the vector connecting ls_pos to the closest point
    double distance = (ls_pos - closest_point).Mag();

    return distance;
}

std::pair<double, double> extrapolatePlsInRZAndGetDeltaZ(
    int pls_idx, int ls_idx,
    // PLS hit position branches (x, y, z)
    const std::vector<float>& pls_hit0_x, const std::vector<float>& pls_hit0_y, const std::vector<float>& pls_hit0_z,
    const std::vector<float>& pls_hit1_x, const std::vector<float>& pls_hit1_y, const std::vector<float>& pls_hit1_z,
    const std::vector<float>& pls_hit2_x, const std::vector<float>& pls_hit2_y, const std::vector<float>& pls_hit2_z,
    const std::vector<float>& pls_hit3_x, const std::vector<float>& pls_hit3_y, const std::vector<float>& pls_hit3_z,
    // LS -> MD link branches
    const std::vector<int>& ls_mdIdx0, const std::vector<int>& ls_mdIdx1,
    // MD position branches (x, y, z)
    const std::vector<float>& md_anchor_x, const std::vector<float>& md_anchor_y, const std::vector<float>& md_anchor_z,
    const std::vector<float>& md_other_x, const std::vector<float>& md_other_y, const std::vector<float>& md_other_z
)  {
    const float invalid_val = -999.0;

    // --- Part 1: Convert PLS hits to (z, r) coordinates to define the line ---
    std::vector<TVector2> pls_hits_zr; // We use TVector2 to store (z, r) points
    
    // Helper lambda to convert (x,y,z) to a (z,r) TVector2
    auto to_zr = [](float x, float y, float z) {
        return TVector2(z, std::sqrt(x*x + y*y));
    };

    if (pls_hit0_x.at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(pls_hit0_x.at(pls_idx), pls_hit0_y.at(pls_idx), pls_hit0_z.at(pls_idx)));
    if (pls_hit1_x.at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(pls_hit1_x.at(pls_idx), pls_hit1_y.at(pls_idx), pls_hit1_z.at(pls_idx)));
    if (pls_hit2_x.at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(pls_hit2_x.at(pls_idx), pls_hit2_y.at(pls_idx), pls_hit2_z.at(pls_idx)));
    if (pls_hit3_x.at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(pls_hit3_x.at(pls_idx), pls_hit3_y.at(pls_idx), pls_hit3_z.at(pls_idx)));

    // We need at least two hits to define a line
    if (pls_hits_zr.size() < 2) return {invalid_val, invalid_val};

    TVector2 p1_zr = pls_hits_zr.front(); // Inner-most point (z1, r1)
    TVector2 p2_zr = pls_hits_zr.back();  // Outer-most point (z2, r2)

    // Safety check: avoid division by zero if all hits are at the same z
    if (std::abs(p2_zr.X() - p1_zr.X()) < 1e-9) return {invalid_val, invalid_val};
    
    // --- Part 2: Calculate the (z, r) position of the LS ---
    int md_idx0 = ls_mdIdx0.at(ls_idx);
    int md_idx1 = ls_mdIdx1.at(ls_idx);
    if (md_idx0 < 0 || md_idx1 < 0) return {invalid_val, invalid_val};

    TVector2 md0_pos_zr = to_zr((md_anchor_x.at(md_idx0)), (md_anchor_y.at(md_idx0)), (md_anchor_z.at(md_idx0)));
    TVector2 md1_pos_zr = to_zr((md_anchor_x.at(md_idx1)), (md_anchor_y.at(md_idx1)), (md_anchor_z.at(md_idx1)));
    //TVector2 ls_pos_zr = (md0_pos_zr + md1_pos_zr) * 0.5;

    // --- Part 3: Extrapolate r = f(z) and calculate the difference in r ---
    // Equation of a line: r = slope * z + intercept
    // TVector2 stores (x,y), which for us corresponds to (z,r).
    // So, point.X() is the z-coordinate and point.Y() is the r-coordinate.
    double slope = (p2_zr.Y() - p1_zr.Y()) / (p2_zr.X() - p1_zr.X()); // (r2 - r1) / (z2 - z1)
    
    // Using the point-slope form: r1 = slope * z1 + intercept  =>  intercept = r1 - slope * z1
    double intercept = p1_zr.Y() - slope * p1_zr.X();

    // Extrapolate: Predict the 'r' position at the 'z' of the LS
    double md_z_position0 = md0_pos_zr.X();
    double md_z_position1 = md1_pos_zr.X();
    double predicted_r0 = slope * md_z_position0 + intercept;
    double predicted_r1 = slope * md_z_position1 + intercept;

    // The final residual is the actual measured radius minus the predicted radius
    double delta_r0 = md0_pos_zr.Y() - predicted_r0;
    double delta_r1 = md1_pos_zr.Y() - predicted_r1;

    return {delta_r0, delta_r1};
}

// MODULE 1: Find Used Indices
void findUsedIndices(std::set<int>& used_ls_indices, std::set<int>& used_pls_indices, const std::vector<int>& tc_pt5Idx_vec, const std::vector<int>& pt5_t5Idx_vec, const std::vector<int>& t5_t3Idx0_vec, const std::vector<int>& t5_t3Idx1_vec, const std::vector<int>& t3_lsIdx0_vec, const std::vector<int>& t3_lsIdx1_vec, const std::vector<int>& tc_pt3Idx_vec, const std::vector<int>& pt3_t3Idx_vec, const std::vector<int>& tc_t5Idx_vec, const std::vector<int>& pt5_plsIdx_vec, const std::vector<int>& pt3_plsIdx_vec) {
    for (int pt5_Idx : tc_pt5Idx_vec) if (pt5_Idx != -999) { if (pt5_plsIdx_vec.at(pt5_Idx) != -999) used_pls_indices.insert(pt5_plsIdx_vec.at(pt5_Idx)); int t5_Idx = pt5_t5Idx_vec.at(pt5_Idx); if (t5_Idx != -999) { int t3_0 = t5_t3Idx0_vec.at(t5_Idx); if (t3_0 != -999) { if (t3_lsIdx0_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_0)); if (t3_lsIdx1_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_0)); } int t3_1 = t5_t3Idx1_vec.at(t5_Idx); if (t3_1 != -999) { if (t3_lsIdx0_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_1)); if (t3_lsIdx1_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_1)); } } }
    for (int pt3_Idx : tc_pt3Idx_vec) if (pt3_Idx != -999) { if (pt3_plsIdx_vec.at(pt3_Idx) != -999) used_pls_indices.insert(pt3_plsIdx_vec.at(pt3_Idx)); int t3_idx = pt3_t3Idx_vec.at(pt3_Idx); if (t3_idx != -999) { if (t3_lsIdx0_vec.at(t3_idx) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_idx)); if (t3_lsIdx1_vec.at(t3_idx) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_idx)); } }
    for (int t5_Idx : tc_t5Idx_vec) if (t5_Idx != -999) { int t3_0 = t5_t3Idx0_vec.at(t5_Idx); if (t3_0 != -999) { if (t3_lsIdx0_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_0)); if (t3_lsIdx1_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_0)); } int t3_1 = t5_t3Idx1_vec.at(t5_Idx); if (t3_1 != -999) { if (t3_lsIdx0_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_1)); if (t3_lsIdx1_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_1)); } }
}

// MODULE 2: Find PLS Match Status
void getPlsMatchStatus(std::map<int, PlsMatchStatus>& simIdx_to_pls_status, const std::vector<int>& pls_isFake_vec, const std::vector<int>& pls_simIdx_vec, const std::set<int>& used_pls_indices) {
    for (size_t j = 0; j < pls_simIdx_vec.size(); ++j) { if (pls_isFake_vec.at(j) == 0) { int sim_idx = pls_simIdx_vec.at(j); if (sim_idx != -999) { if (used_pls_indices.count(j) > 0) simIdx_to_pls_status[sim_idx].used_count++; else simIdx_to_pls_status[sim_idx].unused_count++; } } }
}

// MODULE 3: Fill Detailed LS Histograms
void fillDetailedLsHistograms(bool is_ls_used, bool is_ls_fake, float pt, float eta, float phi, int ls_sim_idx, const std::map<int, PlsMatchStatus>& simIdx_to_pls_status, TH1D* h_pt_f_u, TH1D* h_pt_f_un, TH1D* h_pt_r_u_um, TH1D* h_pt_r_u_b, TH1D* h_pt_r_u_uo, TH1D* h_pt_r_u_uno, TH1D* h_pt_r_un_um, TH1D* h_pt_r_un_b, TH1D* h_pt_r_un_uo, TH1D* h_pt_r_un_uno, TH1D* h_eta_f_u, TH1D* h_eta_f_un, TH1D* h_eta_r_u_um, TH1D* h_eta_r_u_b, TH1D* h_eta_r_u_uo, TH1D* h_eta_r_u_uno, TH1D* h_eta_r_un_um, TH1D* h_eta_r_un_b, TH1D* h_eta_r_un_uo, TH1D* h_eta_r_un_uno, TH1D* h_phi_f_u, TH1D* h_phi_f_un, TH1D* h_phi_r_u_um, TH1D* h_phi_r_u_b, TH1D* h_phi_r_u_uo, TH1D* h_phi_r_u_uno, TH1D* h_phi_r_un_um, TH1D* h_phi_r_un_b, TH1D* h_phi_r_un_uo, TH1D* h_phi_r_un_uno) {
    auto fill_all = [&](TH1D* h_pt, TH1D* h_eta, TH1D* h_phi){ h_pt->Fill(pt); h_eta->Fill(eta); h_phi->Fill(phi); };
    if (is_ls_fake) { if (is_ls_used) fill_all(h_pt_f_u, h_eta_f_u, h_phi_f_u); else fill_all(h_pt_f_un, h_eta_f_un, h_phi_f_un); return; }
    auto match_it = (ls_sim_idx == -999) ? simIdx_to_pls_status.end() : simIdx_to_pls_status.find(ls_sim_idx);
    if (match_it == simIdx_to_pls_status.end()) { if (is_ls_used) fill_all(h_pt_r_u_um, h_eta_r_u_um, h_phi_r_u_um); else fill_all(h_pt_r_un_um, h_eta_r_un_um, h_phi_r_un_um); }
    else { const auto& s = match_it->second; if (s.used_count > 0 && s.unused_count > 0) { if (is_ls_used) fill_all(h_pt_r_u_b, h_eta_r_u_b, h_phi_r_u_b); else fill_all(h_pt_r_un_b, h_eta_r_un_b, h_phi_r_un_b); }
    else if (s.used_count > 0) { if (is_ls_used) fill_all(h_pt_r_u_uo, h_eta_r_u_uo, h_phi_r_u_uo); else fill_all(h_pt_r_un_uo, h_eta_r_un_uo, h_phi_r_un_uo); }
    else { if (is_ls_used) fill_all(h_pt_r_u_uno, h_eta_r_u_uno, h_phi_r_u_uno); else fill_all(h_pt_r_un_uno, h_eta_r_un_uno, h_phi_r_un_uno); } }
}

// MODULE 4: Process Ideal pT2 objects (Now takes ls_isFake_vec)
void processIdealPT2s(const std::vector<float>& pls_pt, const std::vector<float>& pls_eta, const std::vector<float>& pls_phi, const std::vector<int>& pls_isFake, const std::vector<int>& pls_simIdx, const std::vector<float>& ls_pt, const std::vector<float>& ls_eta, const std::vector<float>& ls_phi, const std::vector<int>& ls_isFake, const std::vector<int>& ls_simIdx, const std::set<int>& used_pls_indices, const std::set<int>& used_ls_indices, TH1D* h_delta_pt, TH1D* h_delta_eta, TH1D* h_delta_phi, TH1D* h_delta_R, TH1D* h_extrapolation_dist_3d, TH1D* h_extrapolation_delta_z, const std::vector<float>& pls_hit0_x, const std::vector<float>& pls_hit0_y, const std::vector<float>& pls_hit0_z, const std::vector<float>& pls_hit1_x, const std::vector<float>& pls_hit1_y, const std::vector<float>& pls_hit1_z, const std::vector<float>& pls_hit2_x, const std::vector<float>& pls_hit2_y, const std::vector<float>& pls_hit2_z, const std::vector<float>& pls_hit3_x, const std::vector<float>& pls_hit3_y, const std::vector<float>& pls_hit3_z, const std::vector<int>& ls_mdIdx0, const std::vector<int>& ls_mdIdx1, const std::vector<float>& md_anchor_x, const std::vector<float>& md_anchor_y, const std::vector<float>& md_anchor_z, const std::vector<float>& md_other_x, const std::vector<float>& md_other_y, const std::vector<float>& md_other_z) 
{
    std::map<int, std::vector<int>> ls_simIdx_map;
    for (size_t k = 0; k < ls_pt.size(); ++k) { if (ls_isFake.at(k) == 0 && ls_simIdx.at(k) != -999) { ls_simIdx_map[ls_simIdx.at(k)].push_back(k); } }
    for (size_t i = 0; i < pls_pt.size(); ++i) { if (pls_isFake.at(i) == 0 && used_pls_indices.count(i) == 0 && pls_simIdx.at(i) != -999) { auto it = ls_simIdx_map.find(pls_simIdx.at(i)); 
    if (it != ls_simIdx_map.end()) { for (int ls_idx : it->second) { if (used_ls_indices.count(ls_idx) == 0) { pT2Object ideal_pt2(i, ls_idx, pls_pt, pls_eta, pls_phi, ls_pt, ls_eta, ls_phi); 
    h_delta_pt->Fill(ideal_pt2.delta_pt); 
    h_delta_eta->Fill(ideal_pt2.delta_eta); 
    h_delta_phi->Fill(ideal_pt2.delta_phi); 
    h_delta_R->Fill(ideal_pt2.delta_R); 
    double distance_3d = extrapolatePlsAndGetDistance(i, ls_idx, pls_hit0_x, pls_hit0_y, pls_hit0_z, pls_hit1_x, pls_hit1_y, pls_hit1_z, pls_hit2_x, pls_hit2_y, pls_hit2_z, pls_hit3_x, pls_hit3_y, pls_hit3_z, ls_mdIdx0, ls_mdIdx1, md_anchor_x, md_anchor_y, md_anchor_z, md_other_x, md_other_y, md_other_z); 
    if (distance_3d >= 0) { h_extrapolation_dist_3d->Fill(distance_3d); 
    } const float invalid_val = -999.0; 
    std::pair<double, double> delta_z = extrapolatePlsInRZAndGetDeltaZ(i, ls_idx, pls_hit0_x, pls_hit0_y, pls_hit0_z, pls_hit1_x, pls_hit1_y, pls_hit1_z, pls_hit2_x, pls_hit2_y, pls_hit2_z, pls_hit3_x, pls_hit3_y, pls_hit3_z, ls_mdIdx0, ls_mdIdx1, md_anchor_x, md_anchor_y, md_anchor_z, md_other_x, md_other_y, md_other_z); 
    if (delta_z.first > invalid_val) { h_extrapolation_delta_z->Fill(delta_z.first); } 
    if (delta_z.second > invalid_val) { h_extrapolation_delta_z->Fill(delta_z.second); }} } } } }
}

// MODULE 5: Create and Save a Detailed Stacked Plot
void createAndSaveDetailedStack(const std::string& canvasName, const std::string& stackTitle, const std::string& saveFileName, const std::string& legendHeader, bool setLogY, TH1D* h_fake, TH1D* h_real_unmatched, TH1D* h_real_matchedToUnused, TH1D* h_real_matchedToBoth, TH1D* h_real_matchedToUsed) {
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600); if (setLogY) canvas->SetLogy(); THStack *hs = new THStack("hs", stackTitle.c_str());
    h_fake->SetFillColor(kGray); h_real_unmatched->SetFillColor(kRed-7); h_real_matchedToUnused->SetFillColor(kAzure+1); h_real_matchedToBoth->SetFillColor(kOrange+1); h_real_matchedToUsed->SetFillColor(kGreen+2);
    hs->Add(h_real_matchedToUsed); hs->Add(h_real_matchedToBoth); hs->Add(h_real_matchedToUnused); hs->Add(h_real_unmatched); hs->Add(h_fake); hs->Draw("HIST");
    auto legend = new TLegend(0.40, 0.58, 0.88, 0.88); legend->SetHeader(legendHeader.c_str(), "C");
    legend->AddEntry(h_real_matchedToUsed, "Real, Matched to USED Only", "f"); legend->AddEntry(h_real_matchedToBoth, "Real, Matched to Both", "f"); legend->AddEntry(h_real_matchedToUnused, "Real, Matched to UNUSED Only", "f"); legend->AddEntry(h_real_unmatched, "Real, No Match", "f"); legend->AddEntry(h_fake, "Fake", "f");
    legend->SetBorderSize(0); legend->Draw(); canvas->SaveAs(saveFileName.c_str()); delete canvas;
}

// MODULE 6: Create and Save a Simple 1D Plot
void createAndSaveSimplePlot(const std::string& canvasName, const std::string& saveFileName, TH1D* histo, const std::string& plotTitle, const std::string& xLabel, const std::string& yLabel) {
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    histo->SetTitle(plotTitle.c_str()); histo->GetXaxis()->SetTitle(xLabel.c_str()); histo->GetYaxis()->SetTitle(yLabel.c_str());
    histo->SetLineColor(kAzure+2); histo->SetLineWidth(2); histo->Draw("HIST");
    canvas->SaveAs(saveFileName.c_str()); delete canvas;
}

// MODULE 7: Find LS Match Status
void getLsMatchStatus(std::map<int, LsMatchStatus>& simIdx_to_ls_status, const std::vector<int>& ls_isFake_vec, const std::vector<int>& ls_simIdx_vec, const std::set<int>& used_ls_indices) {
    for (size_t j = 0; j < ls_simIdx_vec.size(); ++j) { if (ls_isFake_vec.at(j) == 0) { int sim_idx = ls_simIdx_vec.at(j); if (sim_idx != -999) { if (used_ls_indices.count(j) > 0) simIdx_to_ls_status[sim_idx].used_count++; else simIdx_to_ls_status[sim_idx].unused_count++; } } }
}

// MODULE 8: Fill Detailed PLS Histograms
void fillDetailedPlsHistograms(bool is_pls_used, bool is_pls_fake, float pt, float eta, float phi, int pls_sim_idx, const std::map<int, LsMatchStatus>& simIdx_to_ls_status, TH1D* h_pt_f_u, TH1D* h_pt_f_un, TH1D* h_pt_r_u_um, TH1D* h_pt_r_u_b, TH1D* h_pt_r_u_uo, TH1D* h_pt_r_u_uno, TH1D* h_pt_r_un_um, TH1D* h_pt_r_un_b, TH1D* h_pt_r_un_uo, TH1D* h_pt_r_un_uno, TH1D* h_eta_f_u, TH1D* h_eta_f_un, TH1D* h_eta_r_u_um, TH1D* h_eta_r_u_b, TH1D* h_eta_r_u_uo, TH1D* h_eta_r_u_uno, TH1D* h_eta_r_un_um, TH1D* h_eta_r_un_b, TH1D* h_eta_r_un_uo, TH1D* h_eta_r_un_uno, TH1D* h_phi_f_u, TH1D* h_phi_f_un, TH1D* h_phi_r_u_um, TH1D* h_phi_r_u_b, TH1D* h_phi_r_u_uo, TH1D* h_phi_r_u_uno, TH1D* h_phi_r_un_um, TH1D* h_phi_r_un_b, TH1D* h_phi_r_un_uo, TH1D* h_phi_r_un_uno) {
    auto fill_all = [&](TH1D* h_pt, TH1D* h_eta, TH1D* h_phi){ h_pt->Fill(pt); h_eta->Fill(eta); h_phi->Fill(phi); };
    if (is_pls_fake) { if (is_pls_used) fill_all(h_pt_f_u, h_eta_f_u, h_phi_f_u); else fill_all(h_pt_f_un, h_eta_f_un, h_phi_f_un); return; }
    auto match_it = (pls_sim_idx == -999) ? simIdx_to_ls_status.end() : simIdx_to_ls_status.find(pls_sim_idx);
    if (match_it == simIdx_to_ls_status.end()) { if (is_pls_used) fill_all(h_pt_r_u_um, h_eta_r_u_um, h_phi_r_u_um); else fill_all(h_pt_r_un_um, h_eta_r_un_um, h_phi_r_un_um); }
    else { const auto& s = match_it->second; if (s.used_count > 0 && s.unused_count > 0) { if (is_pls_used) fill_all(h_pt_r_u_b, h_eta_r_u_b, h_phi_r_u_b); else fill_all(h_pt_r_un_b, h_eta_r_un_b, h_phi_r_un_b); }
    else if (s.used_count > 0) { if (is_pls_used) fill_all(h_pt_r_u_uo, h_eta_r_u_uo, h_phi_r_u_uo); else fill_all(h_pt_r_un_uo, h_eta_r_un_uo, h_phi_r_un_uo); }
    else { if (is_pls_used) fill_all(h_pt_r_u_uno, h_eta_r_u_uno, h_phi_r_u_uno); else fill_all(h_pt_r_un_uno, h_eta_r_un_uno, h_phi_r_un_uno); } }
}

// MODULE 9: Fill Simple Histograms
void fillSimpleHistograms(bool is_used, bool is_fake, float pt, float eta, float phi, TH1D* h_pt_r_u, TH1D* h_pt_f_u, TH1D* h_pt_r_un, TH1D* h_pt_f_un, TH1D* h_eta_r_u, TH1D* h_eta_f_u, TH1D* h_eta_r_un, TH1D* h_eta_f_un, TH1D* h_phi_r_u, TH1D* h_phi_f_u, TH1D* h_phi_r_un, TH1D* h_phi_f_un) {
    if (is_fake) { if (is_used) { h_pt_f_u->Fill(pt); h_eta_f_u->Fill(eta); h_phi_f_u->Fill(phi); } else { h_pt_f_un->Fill(pt); h_eta_f_un->Fill(eta); h_phi_f_un->Fill(phi); } }
    else { if (is_used) { h_pt_r_u->Fill(pt); h_eta_r_u->Fill(eta); h_phi_r_u->Fill(phi); } else { h_pt_r_un->Fill(pt); h_eta_r_un->Fill(eta); h_phi_r_un->Fill(phi); } }
}

// MODULE 10: Create and Save a Simple Stacked Plot
void createAndSaveSimpleStack(const std::string& canvasName, const std::string& stackTitle, const std::string& saveFileName, const std::string& legendHeader, TH1D* h_real, TH1D* h_fake, const std::string& real_label, const std::string& fake_label, Color_t real_color, Color_t fake_color, bool setLogY = false) {
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600); if (setLogY) canvas->SetLogy(); THStack *hs = new THStack("hs", stackTitle.c_str());
    h_real->SetFillColor(real_color); h_fake->SetFillColor(fake_color); hs->Add(h_fake); hs->Add(h_real); hs->Draw("HIST");
    auto legend = new TLegend(0.4, 0.8, 0.9, 0.9); legend->SetHeader(legendHeader.c_str(), "C");
    legend->AddEntry(h_real, real_label.c_str(), "f"); legend->AddEntry(h_fake, fake_label.c_str(), "f");
    legend->SetBorderSize(0); legend->Draw(); canvas->SaveAs(saveFileName.c_str()); delete canvas;
}
