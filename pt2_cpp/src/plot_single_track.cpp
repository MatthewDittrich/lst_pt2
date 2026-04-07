#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <tuple>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TArc.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TF1.h"
#include "TMath.h"
#include "TArrow.h"

#include "rootReader.h"
#include "tools.h"
#include "pt2.h"
#include "extrapolation.h"
#include "extra_cuts.h"

#include <Math/Vector2D.h>

// =========================================================================
// CIRCLE FITTING LOGIC (User Provided)
// =========================================================================
ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double R) {
    if (hits.size() < 2) return {0, 0};

    ROOT::Math::XYVector p1 = hits.front();
    ROOT::Math::XYVector p2 = hits.back();

    double dx = p2.X() - p1.X();
    double dy = p2.Y() - p1.Y();
    double d2 = dx*dx + dy*dy;
    double d = std::sqrt(d2);

    if (d > 2.0 * R || d == 0) return {0, 0};

    double midX = (p1.X() + p2.X()) / 2.0;
    double midY = (p1.Y() + p2.Y()) / 2.0;

    double h = std::sqrt(std::max(0.0, R*R - d2/4.0));

    double c1x = midX + h * (-dy / d);
    double c1y = midY + h * (dx / d);
    double c2x = midX - h * (-dy / d);
    double c2y = midY - h * (dx / d);

    if (hits.size() > 2) {
        ROOT::Math::XYVector pMid = hits[hits.size()/2];
        auto distSq = [](double x1, double y1, double x2, double y2) {
            return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
        };
        double diff1 = std::abs(std::sqrt(distSq(pMid.X(), pMid.Y(), c1x, c1y)) - R);
        double diff2 = std::abs(std::sqrt(distSq(pMid.X(), pMid.Y(), c2x, c2y)) - R);
        return (diff1 < diff2) ? ROOT::Math::XYVector(c1x, c1y) : ROOT::Math::XYVector(c2x, c2y);
    }
    return {c1x, c1y};
}

// Helper to make Z vs Angle a straight line
void unwrapAngle(std::vector<double>& angles) {
    for (size_t i = 1; i < angles.size(); ++i) {
        while (angles[i] - angles[i-1] > M_PI) angles[i] -= 2.0 * M_PI;
        while (angles[i] - angles[i-1] < -M_PI) angles[i] += 2.0 * M_PI;
    }
}

int main() {
    // =========================================================================
    // CONFIGURATION
    // =========================================================================
    int TARGET_COMBO = 8;      
    double MIN_DXY = 1.2;      
    double MAX_DXY = 1.5;      
    bool lowPT = true;

    std::string inputFile = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/ROOT_FILES/LSTNtuple_LowPT.root";
    std::string pixelMapFileDir = "/cmsuf/data/store/user/t2/users/matthew.dittrich/PT2_DATA/PIXEL_MAPS/Pixel_Maps_0p6GeV/";

    SuperbinToDetIdMap sb_POS, sb_NEG, sb_NON;
    loadSuperbinDetIdMap(pixelMapFileDir, sb_POS, sb_NEG, sb_NON);

    rootReader reader;
    if (!reader.Init(inputFile, "tree")) return 1;

    DetIdToLSMap detidToLS;
    pT2Collection pt2s;

    std::cout << "Searching for Real track: Combo " << TARGET_COMBO << ", dXY [" << MIN_DXY << "," << MAX_DXY << "]\n";

    for (Long64_t ievt = 0; ievt < reader.GetEntries(); ++ievt) {
        reader.GetEntry(ievt);
        reader.pls_origin_z.clear(); reader.pls_superbin.clear();
        detidToLS.clear(); pt2s.clear();

        for (size_t k = 0; k < reader.ls_pt->size(); ++k) {
            std::vector<int> dIds = getDetIdsForLS(reader, k);
            for (int d : dIds) detidToLS[d].push_back(k);
        }
        for (size_t j = 0; j < reader.pls_pt->size(); ++j) {
            reader.pls_origin_z.push_back(CalculatePlsZ(reader, j));
            reader.pls_superbin.push_back(CalculateSuperbin(reader, j, lowPT));
            buildPt2sForPLS(j, reader, sb_POS, sb_NEG, sb_NON, detidToLS, pt2s);
        }

        for (auto& pt2 : pt2s) {
            size_t plsIdx = pt2.pls_idx;
            size_t lsIdx  = pt2.ls_idx;

            float reco_pt = reader.pls_pt->at(plsIdx);
            if (reco_pt >= 0.8) continue;

            if (!pt2TruthFinder(reader, plsIdx, lsIdx)) continue;

            int md0_idx = reader.ls_mdIdx0->at(lsIdx);
            int md1_idx = reader.ls_mdIdx1->at(lsIdx);
            int cat0 = extra_cuts::getCategoryFromDetId(reader.md_detId->at(md0_idx));
            int cat1 = extra_cuts::getCategoryFromDetId(reader.md_detId->at(md1_idx));
            int comboIdx = extra_cuts::getConnectionIndex(cat0, cat1);

            if (comboIdx != TARGET_COMBO) continue;

            std::vector<double> heli = extrapolation::extrapolatePlsHelicallyAndGetDistance(plsIdx, lsIdx, reader);
            double dxy0 = heli[0];
            double dz0  = heli[1];
            double dxy1 = heli[2];
            double dz1  = heli[3];

            if (dxy0 >= MIN_DXY && dxy0 <= MAX_DXY) {
                std::cout << ">>> MATCH FOUND! Event: " << ievt << " Track: " << plsIdx << " dXY0: " << dxy0 << "\n";

                int simIdx = reader.pls_simIdx->at(plsIdx);
                bool hasSim = (simIdx >= 0 && static_cast<size_t>(simIdx) < reader.sim_simHitX->size());

                // --- 1. COLLECT RECO HITS ---
                std::vector<double> px, py, pz;
                std::vector<ROOT::Math::XYVector> hits_xy;
                auto add_p = [&](float x, float y, float z) {
                    if (x > -900) { px.push_back(x); py.push_back(y); pz.push_back(z); hits_xy.emplace_back(x, y); }
                };
                add_p(reader.pls_hit0_x->at(plsIdx), reader.pls_hit0_y->at(plsIdx), reader.pls_hit0_z->at(plsIdx));
                add_p(reader.pls_hit1_x->at(plsIdx), reader.pls_hit1_y->at(plsIdx), reader.pls_hit1_z->at(plsIdx));
                add_p(reader.pls_hit2_x->at(plsIdx), reader.pls_hit2_y->at(plsIdx), reader.pls_hit2_z->at(plsIdx));
                if (reader.pls_nhit->at(plsIdx) == 4)
                    add_p(reader.pls_hit3_x->at(plsIdx), reader.pls_hit3_y->at(plsIdx), reader.pls_hit3_z->at(plsIdx));

                std::vector<double> lx = {reader.md_anchor_x->at(md0_idx), reader.md_anchor_x->at(md1_idx)};
                std::vector<double> ly = {reader.md_anchor_y->at(md0_idx), reader.md_anchor_y->at(md1_idx)};
                std::vector<double> lz = {reader.md_anchor_z->at(md0_idx), reader.md_anchor_z->at(md1_idx)};

                // --- 2. CALCULATE CIRCLE & 3D OPTIMIZED PREDICTED POINTS ---
                double cr = (reco_pt / (0.3 * 3.8)) * 100.0;
                ROOT::Math::XYVector center = fitCircleWithFixedRadius(hits_xy, cr);
                double cx = center.X(), cy = center.Y();

                // Replicate Z-slope 'a' and intercept 'b' exactly like extrapolation.h
                double first_angle = std::atan2(py.front() - cy, px.front() - cx);
                double last_angle  = std::atan2(py.back() - cy, px.back() - cx);
                while (last_angle - first_angle > M_PI) last_angle -= 2*M_PI;
                while (last_angle - first_angle < -M_PI) last_angle += 2*M_PI;
                
                double a = (pz.back() - pz.front()) / (last_angle - first_angle);
                double b = pz.front() - a * first_angle;

                // 3D Newton Solver lambda to find the TRUE closest point on the helix
                auto get_3D_extrap = [&](double tx, double ty, double tz) -> std::vector<double> {
                    double phi0 = std::atan2(ty - cy, tx - cx);
                    while (phi0 - last_angle >  M_PI) phi0 -= 2*M_PI;
                    while (phi0 - last_angle < -M_PI) phi0 += 2*M_PI;

                    double p_fine = phi0;
                    for (int iter = 0; iter < 15; ++iter) {
                        double c = std::cos(p_fine), s = std::sin(p_fine);
                        double f = (tx - cx - cr*c)*(cr*s) + (ty - cy - cr*s)*(-cr*c) + a*(a*p_fine + b - tz);
                        double fp = cr*cr + a*a;
                        double next_p = p_fine - (f / fp);
                        if (std::abs(next_p - p_fine) < 1e-6) break;
                        p_fine = next_p;
                    }
                    return {cx + cr*std::cos(p_fine), cy + cr*std::sin(p_fine), a*p_fine + b, p_fine};
                };

                // Run solver for MD0 and MD1
                std::vector<double> p0 = get_3D_extrap(lx[0], ly[0], lz[0]);
                double extrap_x_0 = p0[0], extrap_y_0 = p0[1], extrap_z_0 = p0[2], extrap_a_0_raw = p0[3];

                std::vector<double> p1 = get_3D_extrap(lx[1], ly[1], lz[1]);
                double extrap_x_1 = p1[0], extrap_y_1 = p1[1], extrap_z_1 = p1[2], extrap_a_1_raw = p1[3];

                double extrap_x_arr[2] = {extrap_x_0, extrap_x_1};
                double extrap_y_arr[2] = {extrap_y_0, extrap_y_1};

                // --- 3. COLLECT SIM HITS ---
                std::vector<double> sx, sy, sz;
                if (hasSim) {
                    sx = {reader.sim_simHitX->at(simIdx).begin(), reader.sim_simHitX->at(simIdx).end()};
                    sy = {reader.sim_simHitY->at(simIdx).begin(), reader.sim_simHitY->at(simIdx).end()};
                    sz = {reader.sim_simHitZ->at(simIdx).begin(), reader.sim_simHitZ->at(simIdx).end()};
                }

                // --- 4. PLOTTING ---
                TCanvas* c1 = new TCanvas("c1", "", 1400, 650);
                c1->Divide(2, 1);

                // --- LEFT: XY VIEW ---
                c1->cd(1); gPad->SetGrid();
                TMultiGraph* mgXY = new TMultiGraph();
                
                TGraph* grP_xy = new TGraph(px.size(), px.data(), py.data());
                grP_xy->SetMarkerStyle(kFullCircle); grP_xy->SetMarkerColor(kBlue);
                
                TGraph* grL_xy = new TGraph(lx.size(), lx.data(), ly.data());
                grL_xy->SetMarkerStyle(kFullTriangleUp); grL_xy->SetMarkerColor(kGreen+2);
                grL_xy->SetMarkerSize(1.5);
                
                TGraph* grE_xy = new TGraph(2, extrap_x_arr, extrap_y_arr);
                grE_xy->SetMarkerStyle(kFullStar);
                grE_xy->SetMarkerColor(kRed);
                grE_xy->SetMarkerSize(1.5);

                mgXY->Add(grP_xy, "P"); 
                mgXY->Add(grL_xy, "P");
                mgXY->Add(grE_xy, "P");

                TGraph* grS_xy = nullptr;
                if (hasSim) {
                    grS_xy = new TGraph(sx.size(), sx.data(), sy.data());
                    grS_xy->SetMarkerStyle(kFullStar); 
                    grS_xy->SetMarkerColor(kViolet+2); 
                    mgXY->Add(grS_xy, "P");
                }
                
                mgXY->Draw("A");
                mgXY->SetTitle(Form("XY View (Entry %lld, Track %zu); X [cm]; Y [cm]", ievt, plsIdx));

                TGraph* arcGraph = new TGraph();
                for (int i = 0; i <= 1000; ++i) {
                    double theta = i * 2.0 * M_PI / 1000.0;
                    arcGraph->SetPoint(i, cx + cr * cos(theta), cy + cr * sin(theta));
                }
                arcGraph->SetLineColor(kBlue+2);
                arcGraph->SetLineWidth(2); 
                arcGraph->Draw("L SAME");

                TArrow* arrXY_0 = new TArrow(extrap_x_0, extrap_y_0, lx[0], ly[0], 0.008, "|>");
                arrXY_0->SetLineColor(kBlack); arrXY_0->SetFillColor(kBlack); arrXY_0->SetLineWidth(2);
                arrXY_0->Draw();

                TArrow* arrXY_1 = new TArrow(extrap_x_1, extrap_y_1, lx[1], ly[1], 0.008, "|>");
                arrXY_1->SetLineColor(kBlack); arrXY_1->SetFillColor(kBlack); arrXY_1->SetLineWidth(2);
                arrXY_1->Draw();

                TPaveText* ptInfo = new TPaveText(0.40, 0.15, 0.90, 0.38, "NDC");
                ptInfo->SetFillColor(0); 
                ptInfo->SetTextAlign(12);
                ptInfo->SetTextFont(42);
                ptInfo->SetTextSize(0.035);
                ptInfo->AddText("Reconstructed Track Info:");
                ptInfo->AddText(Form("PLS p_{T} = %.2f GeV", reco_pt));
                ptInfo->AddText(Form("MD0 Category: %d", cat0)); 
                ptInfo->AddText(Form("MD1 Category: %d", cat1));
                ptInfo->AddText(Form("dXY (Pred #rightarrow MD0) = %.3f cm", dxy0));
                ptInfo->AddText(Form("dXY (Pred #rightarrow MD1) = %.3f cm", dxy1));
                ptInfo->Draw();

                // --- RIGHT: Z VS ANGLE ---
                c1->cd(2); gPad->SetGrid();
                gPad->SetLeftMargin(0.15);
                
                std::vector<double> all_angles;
                for(size_t i=0; i<px.size(); ++i) all_angles.push_back(atan2(py[i]-cy, px[i]-cx));
                for(size_t i=0; i<lx.size(); ++i) all_angles.push_back(atan2(ly[i]-cy, lx[i]-cx));
                unwrapAngle(all_angles);

                std::vector<double> pa(all_angles.begin(), all_angles.begin() + px.size());
                std::vector<double> la(all_angles.begin() + px.size(), all_angles.end());

                if (hasSim) {
                    std::vector<double> sa;
                    for(size_t i=0; i<sx.size(); ++i) sa.push_back(atan2(sy[i]-cy, sx[i]-cx));
                    unwrapAngle(sa); 
                }

                // Make sure our 3D optimized angles are unwrapped perfectly alongside the MD hit angles
                double extrap_a_0 = extrap_a_0_raw;
                while (extrap_a_0 - la[0] > M_PI) extrap_a_0 -= 2.0 * M_PI;
                while (extrap_a_0 - la[0] < -M_PI) extrap_a_0 += 2.0 * M_PI;

                double extrap_a_1 = extrap_a_1_raw;
                while (extrap_a_1 - la[1] > M_PI) extrap_a_1 -= 2.0 * M_PI;
                while (extrap_a_1 - la[1] < -M_PI) extrap_a_1 += 2.0 * M_PI;

                double extrap_a_arr[2] = {extrap_a_0, extrap_a_1};
                double extrap_z_arr[2] = {extrap_z_0, extrap_z_1};

                TMultiGraph* mgZ = new TMultiGraph();
                
                TGraph* grP_z = new TGraph(pz.size(), pa.data(), pz.data());
                grP_z->SetMarkerStyle(kFullCircle); grP_z->SetMarkerColor(kBlue);
                
                TF1* fLin = new TF1("fLin", "pol1", pa.front(), pa.back());
                fLin->SetLineColor(kBlue+2); 
                fLin->SetLineWidth(2);
                grP_z->Fit(fLin, "Q");

                TGraph* grL_z = new TGraph(lz.size(), la.data(), lz.data());
                grL_z->SetMarkerStyle(kFullTriangleUp); grL_z->SetMarkerColor(kGreen+2);
                grL_z->SetMarkerSize(1.5);

                TGraph* grE_z = new TGraph(2, extrap_a_arr, extrap_z_arr);
                grE_z->SetMarkerStyle(kFullStar);
                grE_z->SetMarkerColor(kRed);
                grE_z->SetMarkerSize(1.5);
                
                mgZ->Add(grP_z, "P"); 
                mgZ->Add(grL_z, "P");
                mgZ->Add(grE_z, "P");

                mgZ->Draw("A");
                mgZ->SetTitle(Form("Z vs. Angle (Entry %lld, Track %zu); Unwrapped Angle [rad]; Z [cm]", ievt, plsIdx));
                mgZ->GetYaxis()->SetTitleOffset(1.5);

                double min_ang = std::min({pa.front(), pa.back(), la[0], la[1]});
                double max_ang = std::max({pa.front(), pa.back(), la[0], la[1]});
                fLin->SetRange(min_ang - 0.1, max_ang + 0.1);
                fLin->Draw("SAME");

                TArrow* arrZ_0 = new TArrow(extrap_a_0, extrap_z_0, la[0], lz[0], 0.008, "|>");
                arrZ_0->SetLineColor(kBlack); arrZ_0->SetFillColor(kBlack); arrZ_0->SetLineWidth(2);
                arrZ_0->Draw();

                TArrow* arrZ_1 = new TArrow(extrap_a_1, extrap_z_1, la[1], lz[1], 0.008, "|>");
                arrZ_1->SetLineColor(kBlack); arrZ_1->SetFillColor(kBlack); arrZ_1->SetLineWidth(2);
                arrZ_1->Draw();

                TPaveText* ptZInfo = new TPaveText(0.15, 0.76, 0.58, 0.88, "NDC");
                ptZInfo->SetFillColor(0); 
                ptZInfo->SetTextAlign(12);
                ptZInfo->SetTextFont(42);
                ptZInfo->SetTextSize(0.035);
                ptZInfo->AddText("Z Extrapolation Info:");
                ptZInfo->AddText(Form("dZ (Pred #rightarrow MD0) = %.3f cm", dz0));
                ptZInfo->AddText(Form("dZ (Pred #rightarrow MD1) = %.3f cm", dz1));
                ptZInfo->Draw();

                // --- LEGEND ---
                TLegend* leg = new TLegend(0.45, 0.15, 0.75, 0.40);
                leg->AddEntry(grP_xy, "PLS Hits", "P");
                leg->AddEntry(grL_xy, "MD Hits (Real)", "P");
                leg->AddEntry(grE_xy, "Predicted MD Hit", "P");
                if (hasSim) leg->AddEntry("grS_xy", "Sim Hits (True)", "P"); 
                leg->AddEntry(arcGraph, "Fitted Circle/Line", "L");
                leg->AddEntry(arrXY_0, "Distance Measured", "L");
                leg->Draw();

                c1->SaveAs(Form("Track_%lld_%zu_Combo%d.png", ievt, plsIdx, TARGET_COMBO));
                return 0; 
            }
        }
    }
    std::cout << "Finished. No matches found.\n";
    return 0;
}
