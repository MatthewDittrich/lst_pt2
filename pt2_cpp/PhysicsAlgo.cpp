#include "PhysicsAlgo.h"

// ROOT Includes
#include <TVector3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

#include <cmath>
#include <iostream>

namespace PhysicsAlgo {

    // =========================================================================
    //                            HELPER FUNCTIONS
    // =========================================================================

    ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double fixed_radius) {
        if (hits.size() < 2) return {0, 0};

        auto chi_squared_func = [&](const double *p) {
            double h = p[0];
            double k = p[1];
            double chi2 = 0.0;
            for (const auto& hit : hits) {
                double dx = hit.X() - h;
                double dy = hit.Y() - k;
                double r_squared_from_center = dx*dx + dy*dy;
                double residual = std::sqrt(r_squared_from_center) - fixed_radius;
                chi2 += residual * residual;
            }
            return chi2;
        };

        ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kMigrad);
        min.SetMaxFunctionCalls(10000);
        min.SetMaxIterations(1000);
        min.SetTolerance(0.001);
        min.SetPrintLevel(-1);

        ROOT::Math::Functor f(chi_squared_func, 2);
        min.SetFunction(f);

        double avg_x = 0.0, avg_y = 0.0;
        for (const auto& hit : hits) {
            avg_x += hit.X();
            avg_y += hit.Y();
        }
        avg_x /= hits.size();
        avg_y /= hits.size();
        
        min.SetVariable(0, "h", avg_x, 0.1);
        min.SetVariable(1, "k", avg_y, 0.1);

        min.Minimize();

        if (min.Status() == 0) {
            const double *result = min.X();
            return {result[0], result[1]};
        }
        return {0, 0}; 
    }

    // =========================================================================
    //                            MAIN ALGORITHMS
    // =========================================================================

    void findUsedIndices(std::vector<bool>& ls_is_used, std::vector<bool>& pls_is_used, const LSTReader& data) {
        // Resize if needed
        if (data.ls_pt && ls_is_used.size() < data.ls_pt->size()) ls_is_used.resize(data.ls_pt->size(), false);
        if (data.pls_pt && pls_is_used.size() < data.pls_pt->size()) pls_is_used.resize(data.pls_pt->size(), false);

        // Reset
        std::fill(ls_is_used.begin(), ls_is_used.end(), false);
        std::fill(pls_is_used.begin(), pls_is_used.end(), false);

        // Loop 1: pT5
        for (int pt5_Idx : *data.tc_pt5Idx) {
            if (pt5_Idx == -999) continue;
            int pls_idx = data.pt5_plsIdx->at(pt5_Idx);
            if (pls_idx != -999) pls_is_used[pls_idx] = true;
            
            int t5_Idx = data.pt5_t5Idx->at(pt5_Idx);
            if (t5_Idx != -999) {
                int t3_0 = data.t5_t3Idx0->at(t5_Idx);
                if (t3_0 != -999) {
                    int ls0 = data.t3_lsIdx0->at(t3_0);
                    int ls1 = data.t3_lsIdx1->at(t3_0);
                    if (ls0 != -999) ls_is_used[ls0] = true;
                    if (ls1 != -999) ls_is_used[ls1] = true;
                }
                int t3_1 = data.t5_t3Idx1->at(t5_Idx);
                if (t3_1 != -999) {
                    int ls0 = data.t3_lsIdx0->at(t3_1);
                    int ls1 = data.t3_lsIdx1->at(t3_1);
                    if (ls0 != -999) ls_is_used[ls0] = true;
                    if (ls1 != -999) ls_is_used[ls1] = true;
                }
            }
        }

        // Loop 2: pT3
        for (int pt3_Idx : *data.tc_pt3Idx) {
            if (pt3_Idx == -999) continue;
            int pls_idx = data.pt3_plsIdx->at(pt3_Idx);
            if (pls_idx != -999) pls_is_used[pls_idx] = true;
            
            int t3_idx = data.pt3_t3Idx->at(pt3_Idx);
            if (t3_idx != -999) {
                int ls0 = data.t3_lsIdx0->at(t3_idx);
                int ls1 = data.t3_lsIdx1->at(t3_idx);
                if (ls0 != -999) ls_is_used[ls0] = true;
                if (ls1 != -999) ls_is_used[ls1] = true;
            }
        }

        // Loop 3: T5
        for (int t5_Idx : *data.tc_t5Idx) {
            if (t5_Idx == -999) continue;
            int t3_0 = data.t5_t3Idx0->at(t5_Idx);
            if (t3_0 != -999) {
                int ls0 = data.t3_lsIdx0->at(t3_0);
                int ls1 = data.t3_lsIdx1->at(t3_0);
                if (ls0 != -999) ls_is_used[ls0] = true;
                if (ls1 != -999) ls_is_used[ls1] = true;
            }
            int t3_1 = data.t5_t3Idx1->at(t5_Idx);
            if (t3_1 != -999) {
                int ls0 = data.t3_lsIdx0->at(t3_1);
                int ls1 = data.t3_lsIdx1->at(t3_1);
                if (ls0 != -999) ls_is_used[ls0] = true;
                if (ls1 != -999) ls_is_used[ls1] = true;
            }
        }
    }

   /* void getPlsMatchStatus(std::map<int, PlsMatchStatus>& simIdx_to_pls_status, const std::vector<int>& pls_isFake, const std::vector<int>& pls_simIdx, const std::vector<bool>& pls_is_used) {
        for (size_t j = 0; j < pls_simIdx.size(); ++j) {
            if (pls_isFake.at(j) == 0) {
                int sim_idx = pls_simIdx.at(j);
                if (sim_idx != -999) {
                    if (pls_is_used[j]) simIdx_to_pls_status[sim_idx].used_count++;
                    else simIdx_to_pls_status[sim_idx].unused_count++;
                }
            }
        }
    }

    void getLsMatchStatus(std::map<int, LsMatchStatus>& simIdx_to_ls_status, const std::vector<int>& ls_isFake, const std::vector<int>& ls_simIdx, const std::vector<bool>& ls_is_used) {
        for (size_t j = 0; j < ls_simIdx.size(); ++j) {
            if (ls_isFake.at(j) == 0) {
                int sim_idx = ls_simIdx.at(j);
                if (sim_idx != -999) {
                    if (ls_is_used[j]) simIdx_to_ls_status[sim_idx].used_count++;
                    else simIdx_to_ls_status[sim_idx].unused_count++;
                }
            }
        }
    }

    void processIdealPT2s(const LSTReader& data, const std::vector<bool>& pls_is_used, const std::vector<bool>& ls_is_used, const std::map<int, std::vector<int>>& ls_simIdx_map, HistManager& hists, long long& ideal_pt2_count) {
        // We do NOT rebuild the map here anymore. It is passed as 'ls_simIdx_map'.
        
        // Loop over PLS
        for (size_t i = 0; i < data.pls_pt->size(); ++i) {
            // Optimization: Vector<bool> check
            if (data.pls_isFake->at(i) == 0 && !pls_is_used[i] && data.pls_simIdx->at(i) != -999) {
                auto it = ls_simIdx_map.find(data.pls_simIdx->at(i));
                if (it != ls_simIdx_map.end()) {
                    for (int ls_idx : it->second) {
                        // Optimization: Vector<bool> check
                        if (!ls_is_used[ls_idx]) {
                            // Found Ideal pT2 Pair
                            ideal_pt2_count++;
                            pT2Object ideal_pt2(i, ls_idx, *data.pls_pt, *data.pls_eta, *data.pls_phi, *data.ls_pt, *data.ls_eta, *data.ls_phi);
                            
                            hists.h_delta_pt->Fill(ideal_pt2.delta_pt);
                            hists.h_delta_eta->Fill(ideal_pt2.delta_eta);
                            hists.h_delta_phi->Fill(ideal_pt2.delta_phi);
                            hists.h_delta_R->Fill(ideal_pt2.delta_R);

                            float reco_pt_pls = data.pls_pt->at(i);
                            if (reco_pt_pls <= 2.0) {
                                int sim_idx = data.pls_simIdx->at(i);
                                if (sim_idx >= 0 && sim_idx < data.sim_pt->size()) {
                                    float true_pt = data.sim_pt->at(sim_idx);
                                    if (true_pt <= 2.0) {
                                        // 1. Helical Extrapolation
                                        auto distance_3d = extrapolatePlsHelicallyAndGetDistance(i, ls_idx, data);
                                        
                                        if (distance_3d.first >= 0) {
                                            hists.h_extrapolation_dist_3d->Fill(distance_3d.first);
                                            int md_idx0 = data.ls_mdIdx0->at(ls_idx);
                                            if (md_idx0 >= 0 && md_idx0 < data.md_layer->size()) {
                                                int layer0 = data.md_layer->at(md_idx0);
                                                if (layer0 < hists.h_dist_3d_by_layer.size()) hists.h_dist_3d_by_layer[layer0]->Fill(distance_3d.first);
                                            }
                                        }
                                        if (distance_3d.second >= 0) {
                                            hists.h_extrapolation_dist_3d->Fill(distance_3d.second);
                                            int md_idx1 = data.ls_mdIdx1->at(ls_idx);
                                            if (md_idx1 >= 0 && md_idx1 < data.md_layer->size()) {
                                                int layer1 = data.md_layer->at(md_idx1);
                                                if (layer1 < hists.h_dist_3d_by_layer.size()) hists.h_dist_3d_by_layer[layer1]->Fill(distance_3d.second);
                                            }
                                        }
                                    }
                                }
                            }

                            // 2. LS -> PLS Reverse Extrapolation (Delta R)
                            std::vector<double> delta_r_results = extrapolateLsInRZAndGetDeltaR_AllHits(i, ls_idx, data);
                            for (size_t hit_num = 0; hit_num < delta_r_results.size(); ++hit_num) {
                                if (hit_num < hists.h_extrapolation_delta_r_reverse_by_hit.size()) {
                                    hists.h_extrapolation_delta_r_reverse_by_hit[hit_num]->Fill(delta_r_results[hit_num]);
                                    hists.h_extrapolation_delta_r_reverse_combined->Fill(delta_r_results[hit_num]);
                                }
                            }

                            // 3. PLS -> LS Extrapolation (Delta Z/R)
                            std::pair<double, double> delta_r_pair = extrapolatePlsInRZAndGetDeltaZ(i, ls_idx, data);
                            const float invalid_val = -999.0;
                            int md_idx0 = data.ls_mdIdx0->at(ls_idx);
                            int md_idx1 = data.ls_mdIdx1->at(ls_idx);

                            if (delta_r_pair.first > invalid_val && md_idx0 >= 0) {
                                hists.h_extrapolation_delta_z->Fill(delta_r_pair.first);
                                int layer0 = data.md_layer->at(md_idx0);
                                if (layer0 < hists.h_delta_r_by_layer.size()) hists.h_delta_r_by_layer[layer0]->Fill(delta_r_pair.first);
                            }
                            if (delta_r_pair.second > invalid_val && md_idx1 >= 0) {
                                hists.h_extrapolation_delta_z->Fill(delta_r_pair.second);
                                int layer1 = data.md_layer->at(md_idx1);
                                if (layer1 < hists.h_delta_r_by_layer.size()) hists.h_delta_r_by_layer[layer1]->Fill(delta_r_pair.second);
                            }

                        } // End Used LS check
                    } // End LS loop
                } // End Map find
            } // End PLS Check
        } // End PLS Loop
    }
    */
    
    void runPT2Physics(int pls_idx, 
                   const LSTReader& data, 
                   const std::vector<bool>& ls_is_used, 
                   const std::map<int, std::vector<int>>& ls_simIdx_map, 
                   HistManager& hists, 
                   long long& ideal_pt2_count) 
    {
    // 1. Basic Validation: PLS must be Real and have a Sim Index
    // (Note: The caller usually checks if it is 'Used', but we check validity here)
    if (data.pls_isFake->at(pls_idx) == 0 && data.pls_simIdx->at(pls_idx) != -999) {
        
        // 2. Look for matching Line Segments in the pre-built map
        int sim_idx = data.pls_simIdx->at(pls_idx);
        auto it = ls_simIdx_map.find(sim_idx);

        if (it != ls_simIdx_map.end()) {
            
            // 3. Loop over all LSs that share this Sim Index
            for (int ls_idx : it->second) {
                
                // 4. Critical Check: Is the LS "Unused"?
                // We typically only want to form pT2s from leftover objects
                if (!ls_is_used[ls_idx]) {
                    
                    // --- FOUND AN IDEAL PT2 PAIR ---
                    ideal_pt2_count++;

                    // A. Calculate and Fill Deltas
                    pT2Object ideal_pt2(pls_idx, ls_idx, *data.pls_pt, *data.pls_eta, *data.pls_phi, *data.ls_pt, *data.ls_eta, *data.ls_phi);
                    hists.h_delta_pt->Fill(ideal_pt2.delta_pt);
                    hists.h_delta_eta->Fill(ideal_pt2.delta_eta);
                    hists.h_delta_phi->Fill(ideal_pt2.delta_phi);
                    hists.h_delta_R->Fill(ideal_pt2.delta_R);

                    // B. Helical Extrapolation (Truth Matched Low pT)
                    float reco_pt_pls = data.pls_pt->at(pls_idx);
                    if (reco_pt_pls <= 2.0) {
                        if (sim_idx >= 0 && sim_idx < data.sim_pt->size()) {
                            float true_pt = data.sim_pt->at(sim_idx);
                            
                            // Only perform expensive extrapolation for low pT truth particles
                            if (true_pt <= 2.0) {
                                auto distance_3d = extrapolatePlsHelicallyAndGetDistance(pls_idx, ls_idx, data);
                                
                                // Fill Hit 0 Distance
                                if (distance_3d.first >= 0) {
                                    hists.h_extrapolation_dist_3d->Fill(distance_3d.first);
                                    int md_idx0 = data.ls_mdIdx0->at(ls_idx);
                                    if (md_idx0 >= 0 && md_idx0 < data.md_layer->size()) {
                                        int layer0 = data.md_layer->at(md_idx0);
                                        if (layer0 < hists.h_dist_3d_by_layer.size()) {
                                            hists.h_dist_3d_by_layer[layer0]->Fill(distance_3d.first);
                                        }
                                    }
                                }
                                
                                // Fill Hit 1 Distance
                                if (distance_3d.second >= 0) {
                                    hists.h_extrapolation_dist_3d->Fill(distance_3d.second);
                                    int md_idx1 = data.ls_mdIdx1->at(ls_idx);
                                    if (md_idx1 >= 0 && md_idx1 < data.md_layer->size()) {
                                        int layer1 = data.md_layer->at(md_idx1);
                                        if (layer1 < hists.h_dist_3d_by_layer.size()) {
                                            hists.h_dist_3d_by_layer[layer1]->Fill(distance_3d.second);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // C. Reverse Extrapolation: LS -> PLS (Linear R-Z)
                    std::vector<double> delta_r_results = extrapolateLsInRZAndGetDeltaR_AllHits(pls_idx, ls_idx, data);
                    for (size_t hit_num = 0; hit_num < delta_r_results.size(); ++hit_num) {
                        // Fill per-hit histogram
                        if (hit_num < hists.h_extrapolation_delta_r_reverse_by_hit.size()) {
                            hists.h_extrapolation_delta_r_reverse_by_hit[hit_num]->Fill(delta_r_results[hit_num]);
                        }
                        // Fill combined histogram
                        hists.h_extrapolation_delta_r_reverse_combined->Fill(delta_r_results[hit_num]);
                    }

                    // D. Forward Extrapolation: PLS -> LS (Linear R-Z)
                    std::pair<double, double> delta_r_pair = extrapolatePlsInRZAndGetDeltaZ(pls_idx, ls_idx, data);
                    const float invalid_val = -999.0;
                    int md_idx0 = data.ls_mdIdx0->at(ls_idx);
                    int md_idx1 = data.ls_mdIdx1->at(ls_idx);

                    // Fill MD0
                    if (delta_r_pair.first > invalid_val && md_idx0 >= 0) {
                        hists.h_extrapolation_delta_z->Fill(delta_r_pair.first);
                        int layer0 = data.md_layer->at(md_idx0);
                        if (layer0 < hists.h_delta_r_by_layer.size()) {
                            hists.h_delta_r_by_layer[layer0]->Fill(delta_r_pair.first);
                        }
                    }
                    
                    // Fill MD1
                    if (delta_r_pair.second > invalid_val && md_idx1 >= 0) {
                        hists.h_extrapolation_delta_z->Fill(delta_r_pair.second);
                        int layer1 = data.md_layer->at(md_idx1);
                        if (layer1 < hists.h_delta_r_by_layer.size()) {
                            hists.h_delta_r_by_layer[layer1]->Fill(delta_r_pair.second);
                        }
                    }

                } // End if(!ls_is_used)
            } // End Loop over matched LSs
        } // End if(it != end)
    } // End if(valid PLS)
    }
        
    // =========================================================================
    //                            EXTRAPOLATION IMPL (Unchanged)
    // =========================================================================

    std::pair<double, double> extrapolatePlsHelicallyAndGetDistance(int pls_idx, int ls_idx, const LSTReader& data) {
        const float invalid_hit_val = -999.0f;
        const std::pair<double, double> invalid_return = {-1.0, -1.0};

        std::vector<ROOT::Math::XYVector> pls_hits_xy;
        std::vector<TVector3> pls_hits_3d;
        auto add_hit = [&](float x, float y, float z) {
            if (x > invalid_hit_val) {
                pls_hits_xy.emplace_back(x, y);
                pls_hits_3d.emplace_back(x, y, z);
            }
        };

        add_hit(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
        add_hit(data.pls_hit1_x->at(pls_idx), data.pls_hit1_y->at(pls_idx), data.pls_hit1_z->at(pls_idx));
        add_hit(data.pls_hit2_x->at(pls_idx), data.pls_hit2_y->at(pls_idx), data.pls_hit2_z->at(pls_idx));
        add_hit(data.pls_hit3_x->at(pls_idx), data.pls_hit3_y->at(pls_idx), data.pls_hit3_z->at(pls_idx));

        if (pls_hits_xy.size() < 3) return invalid_return;

        float reco_pt = data.pls_pt->at(pls_idx);
        const double B_FIELD_TESLA = 3.8; 
        double R_pred_cm = (reco_pt / (0.3 * B_FIELD_TESLA)) * 100.0;
        
        ROOT::Math::XYVector circle_center = fitCircleWithFixedRadius(pls_hits_xy, R_pred_cm);

        if (circle_center.X() == 0 && circle_center.Y() == 0) return invalid_return;
        double circle_radius = R_pred_cm;

        TGraph z_vs_angle_graph;
        double last_angle = 0;
        bool first_point = true;
        for (const auto& hit_xy : pls_hits_xy) {
            double angle = TMath::ATan2(hit_xy.Y() - circle_center.Y(), hit_xy.X() - circle_center.X());
            if (!first_point) {
                while (angle - last_angle > TMath::Pi()) angle -= 2 * TMath::Pi();
                while (angle - last_angle < -TMath::Pi()) angle += 2 * TMath::Pi();
            }
            for (const auto& hit_3d : pls_hits_3d) {
                if (std::abs(hit_3d.X() - hit_xy.X()) < 1e-6 && std::abs(hit_3d.Y() - hit_xy.Y()) < 1e-6) {
                    z_vs_angle_graph.SetPoint(z_vs_angle_graph.GetN(), angle, hit_3d.Z());
                    break;
                }
            }
            last_angle = angle;
            first_point = false;
        }

        if (z_vs_angle_graph.GetN() < 2) return invalid_return;
        z_vs_angle_graph.Fit("pol1", "Q");
        TF1* fit_func = z_vs_angle_graph.GetFunction("pol1");
        if (!fit_func) return invalid_return;

        double b = fit_func->GetParameter(0); 
        double a = fit_func->GetParameter(1); 

        int md_idx0 = data.ls_mdIdx0->at(ls_idx);
        int md_idx1 = data.ls_mdIdx1->at(ls_idx);
        if (md_idx0 < 0 || md_idx1 < 0) return invalid_return;

        TVector3 md0_anchor(data.md_anchor_x->at(md_idx0), data.md_anchor_y->at(md_idx0), data.md_anchor_z->at(md_idx0));
        TVector3 md1_anchor(data.md_anchor_x->at(md_idx1), data.md_anchor_y->at(md_idx1), data.md_anchor_z->at(md_idx1));

        auto unwrap_to_last = [&](double phi_target) {
            const int N = z_vs_angle_graph.GetN();
            const double phi_last = z_vs_angle_graph.GetX()[N-1];
            while (phi_target - phi_last >  TMath::Pi())  phi_target -= 2*TMath::Pi();
            while (phi_target - phi_last < -TMath::Pi())  phi_target += 2*TMath::Pi();
            return phi_target;
        };

        auto helixClosestPhiToPoint = [&](double xa, double ya, double za, double xc, double yc, double R, double a_param, double b_param, double phi_init) {
            double phi = phi_init;
            for (int iter = 0; iter < 12; ++iter) {
                double c = std::cos(phi);
                double s = std::sin(phi);
                double dx = xa - xc - R * c; 
                double dy = ya - yc - R * s; 
                double dz = a_param * phi + b_param - za; 

                double f = dx * (R * s) + dy * (-R * c) + a_param * dz;
                double fp = R*R + R*(dx * c + dy * s) + a_param * a_param;

                if (std::abs(fp) < 1e-12) break;
                double step = f / fp;
                phi -= step;
                if (std::abs(step) > 1.0) phi += (step > 0) ? 1.0 : -1.0; 
                if (std::abs(step) < 1e-10) break; 
            }
            return phi;
        };

        auto initial_phi = [&](double xa, double ya) {
            double phi0 = std::atan2(ya - circle_center.Y(), xa - circle_center.X());
            return unwrap_to_last(phi0);
        };

        double phiA0 = initial_phi(md0_anchor.X(), md0_anchor.Y());
        double phiA = helixClosestPhiToPoint(md0_anchor.X(), md0_anchor.Y(), md0_anchor.Z(), circle_center.X(), circle_center.Y(), circle_radius, a, b, phiA0);
        double xA = circle_center.X() + circle_radius * std::cos(phiA);
        double yA = circle_center.Y() + circle_radius * std::sin(phiA);
        double zA = a * phiA + b;
        TVector3 helixA(xA, yA, zA);
        double dist0 = (md0_anchor - helixA).Mag();

        double phiB0 = initial_phi(md1_anchor.X(), md1_anchor.Y());
        double phiB = helixClosestPhiToPoint(md1_anchor.X(), md1_anchor.Y(), md1_anchor.Z(), circle_center.X(), circle_center.Y(), circle_radius, a, b, phiB0);
        double xB = circle_center.X() + circle_radius * std::cos(phiB);
        double yB = circle_center.Y() + circle_radius * std::sin(phiB);
        double zB = a * phiB + b;
        TVector3 helixB(xB, yB, zB);
        double dist1 = (md1_anchor - helixB).Mag();

        return {dist0, dist1};
    }

    std::vector<double> extrapolateLsInRZAndGetDeltaR_AllHits(int pls_idx, int ls_idx, const LSTReader& data) {
        const float invalid_val = -999.0;
        std::vector<double> results;

        int md_idx0 = data.ls_mdIdx0->at(ls_idx);
        int md_idx1 = data.ls_mdIdx1->at(ls_idx);
        if (md_idx0 < 0 || md_idx1 < 0) return results;

        auto to_zr = [](float x, float y, float z) { return TVector2(z, std::sqrt(x*x + y*y)); };
        TVector2 p1_zr = to_zr(data.md_anchor_x->at(md_idx0), data.md_anchor_y->at(md_idx0), data.md_anchor_z->at(md_idx0));
        TVector2 p2_zr = to_zr(data.md_anchor_x->at(md_idx1), data.md_anchor_y->at(md_idx1), data.md_anchor_z->at(md_idx1));
        if (std::abs(p2_zr.X() - p1_zr.X()) < 1e-9) return results;

        double slope = (p2_zr.Y() - p1_zr.Y()) / (p2_zr.X() - p1_zr.X());
        double intercept = p1_zr.Y() - slope * p1_zr.X();

        std::vector<TVector2> target_pls_hits_zr;
        if (data.pls_hit0_x->at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx)));
        if (data.pls_hit1_x->at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(data.pls_hit1_x->at(pls_idx), data.pls_hit1_y->at(pls_idx), data.pls_hit1_z->at(pls_idx)));
        if (data.pls_hit2_x->at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(data.pls_hit2_x->at(pls_idx), data.pls_hit2_y->at(pls_idx), data.pls_hit2_z->at(pls_idx)));
        if (data.pls_hit3_x->at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(data.pls_hit3_x->at(pls_idx), data.pls_hit3_y->at(pls_idx), data.pls_hit3_z->at(pls_idx)));

        for (const auto& target_hit : target_pls_hits_zr) {
            double predicted_r = slope * target_hit.X() + intercept;
            results.push_back(target_hit.Y() - predicted_r);
        }
        return results;
    }

    std::pair<double, double> extrapolatePlsInRZAndGetDeltaZ(int pls_idx, int ls_idx, const LSTReader& data) {
        const float invalid_val = -999.0;
        std::vector<TVector2> pls_hits_zr;
        auto to_zr = [](float x, float y, float z) { return TVector2(z, std::sqrt(x*x + y*y)); };

        if (data.pls_hit0_x->at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx)));
        if (data.pls_hit1_x->at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(data.pls_hit1_x->at(pls_idx), data.pls_hit1_y->at(pls_idx), data.pls_hit1_z->at(pls_idx)));
        if (data.pls_hit2_x->at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(data.pls_hit2_x->at(pls_idx), data.pls_hit2_y->at(pls_idx), data.pls_hit2_z->at(pls_idx)));
        if (data.pls_hit3_x->at(pls_idx) > invalid_val) pls_hits_zr.push_back(to_zr(data.pls_hit3_x->at(pls_idx), data.pls_hit3_y->at(pls_idx), data.pls_hit3_z->at(pls_idx)));

        if (pls_hits_zr.size() < 2) return {invalid_val, invalid_val};

        TVector2 p1_zr = pls_hits_zr.front();
        TVector2 p2_zr = pls_hits_zr.back();

        if (std::abs(p2_zr.X() - p1_zr.X()) < 1e-9) return {invalid_val, invalid_val};
        
        int md_idx0 = data.ls_mdIdx0->at(ls_idx);
        int md_idx1 = data.ls_mdIdx1->at(ls_idx);
        if (md_idx0 < 0 || md_idx1 < 0) return {invalid_val, invalid_val};

        TVector2 md0_pos_zr = to_zr(data.md_anchor_x->at(md_idx0), data.md_anchor_y->at(md_idx0), data.md_anchor_z->at(md_idx0));
        TVector2 md1_pos_zr = to_zr(data.md_anchor_x->at(md_idx1), data.md_anchor_y->at(md_idx1), data.md_anchor_z->at(md_idx1));

        double slope = (p2_zr.Y() - p1_zr.Y()) / (p2_zr.X() - p1_zr.X());
        double intercept = p1_zr.Y() - slope * p1_zr.X();

        double predicted_r0 = slope * md0_pos_zr.X() + intercept;
        double predicted_r1 = slope * md1_pos_zr.X() + intercept;

        return {md0_pos_zr.Y() - predicted_r0, md1_pos_zr.Y() - predicted_r1};
    }
}
