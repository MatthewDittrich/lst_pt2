#include "extrapolation.h"
#include "rootReader.h"
// ROOT Includes
#include <TVector3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Functor.h>
#include <cmath>
#include <iostream>

namespace extrapolation {

    // =========================================================================
    //                            HELPER FUNCTIONS
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

    std::pair<double, double> extrapolatePlsHelicallyAndGetDistance(int pls_idx, int ls_idx, const rootReader& data) {
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

        TGraph z_vs_angle_graph;
        double first_angle = 0, first_z = 0;
        double last_angle = 0, last_z = 0;
    
        for (size_t i = 0; i < pls_hits_xy.size(); ++i) {
            double angle = TMath::ATan2(pls_hits_xy[i].Y() - circle_center.Y(), pls_hits_xy[i].X() - circle_center.X());
            if (i == 0) {
                first_angle = angle;
                first_z = pls_hits_3d[i].Z();
                last_angle = angle;
            } else {
                // Unwrap angle logic to ensure continuity
                while (angle - last_angle > TMath::Pi()) angle -= 2 * TMath::Pi();
                while (angle - last_angle < -TMath::Pi()) angle += 2 * TMath::Pi();
                last_angle = angle;
                last_z = pls_hits_3d[i].Z();
            }
        }

        double delta_phi = last_angle - first_angle;
        if (std::abs(delta_phi) < 1e-9) return invalid_return;

        double a = (last_z - first_z) / delta_phi; // slope dz/dphi
        double b = first_z - a * first_angle;      // intercept

         // 4. Target Mini-Doublets
        int md_idx0 = data.ls_mdIdx0->at(ls_idx);
        int md_idx1 = data.ls_mdIdx1->at(ls_idx);
        if (md_idx0 < 0 || md_idx1 < 0) return invalid_return;
    
        TVector3 md0_anchor(data.md_anchor_x->at(md_idx0), data.md_anchor_y->at(md_idx0), data.md_anchor_z->at(md_idx0));
        TVector3 md1_anchor(data.md_anchor_x->at(md_idx1), data.md_anchor_y->at(md_idx1), data.md_anchor_z->at(md_idx1));

        // 5. Helical Projection Logic
        auto helixClosestPhiToPoint = [&](double xa, double ya, double za, double xc, double yc, double R, double a_param, double b_param, double phi_init) {
            double phi = phi_init;
            for (int iter = 0; iter < 12; ++iter) {
                double c = std::cos(phi), s = std::sin(phi);
                double dx = xa - xc - R * c;
                double dy = ya - yc - R * s;
                double dz = a_param * phi + b_param - za;
                double f = dx * (R * s) + dy * (-R * c) + a_param * dz;
                double fp = R*R + R*(dx * c + dy * s) + a_param * a_param;
                if (std::abs(fp) < 1e-12) break;
                phi -= f / fp;
            }
            return phi;
        };

        auto get_dist = [&](const TVector3& target) {
        double phi0 = std::atan2(target.Y() - circle_center.Y(), target.X() - circle_center.X());
        
        // FIX: Use last_angle from the loop above instead of accessing the graph
        while (phi0 - last_angle >  TMath::Pi()) phi0 -= 2*TMath::Pi();
        while (phi0 - last_angle < -TMath::Pi()) phi0 += 2*TMath::Pi();

        double phi_fine = helixClosestPhiToPoint(target.X(), target.Y(), target.Z(), circle_center.X(), circle_center.Y(), R_pred_cm, a, b, phi0);
        TVector3 helix_pt(circle_center.X() + R_pred_cm * std::cos(phi_fine),
                          circle_center.Y() + R_pred_cm * std::sin(phi_fine),
                          a * phi_fine + b);
        return (target - helix_pt).Mag();
    };

    return {get_dist(md0_anchor), get_dist(md1_anchor)};
    }

    std::pair<double, double> extrapolatePlsInRZAndGetDeltaR(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -999.0;
        const std::pair<double, double> invalid_pair = {-999.0, -999.0};

        std::vector<double> pls_z, pls_r;
        auto add_pls_hit = [&](float x, float y, float z) {
            if (x > invalid_val) {
                pls_z.push_back(z);
                pls_r.push_back(std::sqrt(x*x + y*y));
            }
        };
        add_pls_hit(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
        add_pls_hit(data.pls_hit1_x->at(pls_idx), data.pls_hit1_y->at(pls_idx), data.pls_hit1_z->at(pls_idx));
        add_pls_hit(data.pls_hit2_x->at(pls_idx), data.pls_hit2_y->at(pls_idx), data.pls_hit2_z->at(pls_idx));
        add_pls_hit(data.pls_hit3_x->at(pls_idx), data.pls_hit3_y->at(pls_idx), data.pls_hit3_z->at(pls_idx));

        if (pls_z.size() < 2) return invalid_pair;

        double sumZ = 0, sumR = 0, sumZR = 0, sumZ2 = 0;
        for (size_t i = 0; i < pls_z.size(); ++i) {
            sumZ += pls_z[i];
            sumR += pls_r[i];
            sumZR += pls_z[i] * pls_r[i];
            sumZ2 += pls_z[i] * pls_z[i];
        }
        double n = pls_z.size();
        double denom = (n * sumZ2 - sumZ * sumZ);
        if (std::abs(denom) < 1e-12) return invalid_pair;

        double slope = (n * sumZR - sumZ * sumR) / denom;
        double intercept = (sumR - slope * sumZ) / n;

        int md0 = data.ls_mdIdx0->at(ls_idx);
        int md1 = data.ls_mdIdx1->at(ls_idx);
        if (md0 < 0 || md1 < 0) return invalid_pair;

        double z_ls0 = data.md_anchor_z->at(md0);
        double r_ls0_act = std::sqrt(std::pow(data.md_anchor_x->at(md0), 2) + std::pow(data.md_anchor_y->at(md0), 2));
        double z_ls1 = data.md_anchor_z->at(md1);
        double r_ls1_act = std::sqrt(std::pow(data.md_anchor_x->at(md1), 2) + std::pow(data.md_anchor_y->at(md1), 2));

        return {r_ls0_act - (slope * z_ls0 + intercept), r_ls1_act - (slope * z_ls1 + intercept)};
    }

    std::pair<double, double> extrapolateSimplePointingInRZ(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -999.0;
        const std::pair<double, double> invalid_pair = {-999.0, -999.0};

        double z_first = 0, r_first = 0, z_last = 0, r_last = 0;
        bool found_first = false, found_last = false;

        if (data.pls_hit0_x->at(pls_idx) > invalid_val) {
            z_first = data.pls_hit0_z->at(pls_idx);
            r_first = std::sqrt(std::pow(data.pls_hit0_x->at(pls_idx), 2) + std::pow(data.pls_hit0_y->at(pls_idx), 2));
            found_first = true;
        } else if (data.pls_hit1_x->at(pls_idx) > invalid_val) {
            z_first = data.pls_hit1_z->at(pls_idx);
            r_first = std::sqrt(std::pow(data.pls_hit1_x->at(pls_idx), 2) + std::pow(data.pls_hit1_y->at(pls_idx), 2));
            found_first = true;
        }

        if (data.pls_hit3_x->at(pls_idx) > invalid_val) {
            z_last = data.pls_hit3_z->at(pls_idx);
            r_last = std::sqrt(std::pow(data.pls_hit3_x->at(pls_idx), 2) + std::pow(data.pls_hit3_y->at(pls_idx), 2));
            found_last = true;
        } else if (data.pls_hit2_x->at(pls_idx) > invalid_val) {
            z_last = data.pls_hit2_z->at(pls_idx);
            r_last = std::sqrt(std::pow(data.pls_hit2_x->at(pls_idx), 2) + std::pow(data.pls_hit2_y->at(pls_idx), 2));
            found_last = true;
        }

        if (!found_first || !found_last || std::abs(z_last - z_first) < 1e-4) return invalid_pair;

        double m = (r_last - r_first) / (z_last - z_first);
        double b = r_first - m * z_first;

        int md0 = data.ls_mdIdx0->at(ls_idx), md1 = data.ls_mdIdx1->at(ls_idx);
        if (md0 < 0 || md1 < 0) return invalid_pair;

        double z0 = data.md_anchor_z->at(md0), r0 = std::sqrt(std::pow(data.md_anchor_x->at(md0), 2) + std::pow(data.md_anchor_y->at(md0), 2));
        double z1 = data.md_anchor_z->at(md1), r1 = std::sqrt(std::pow(data.md_anchor_x->at(md1), 2) + std::pow(data.md_anchor_y->at(md1), 2));

        return {r0 - (m * z0 + b), r1 - (m * z1 + b)};
    }

    double calculateDeltaAngle(int pls_idx, int ls_idx, const rootReader& data) {
        const float inv = -999.0;
        TVector3 vPixel;
        if (data.pls_hit0_x->at(pls_idx) > inv) {
            vPixel.SetXYZ(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
        } else {
            return -999.0;
        }
        int md1_idx = data.ls_mdIdx1->at(ls_idx);
        if (md1_idx < 0) return -999.0;
        TVector3 vStrip(data.md_anchor_x->at(md1_idx), data.md_anchor_y->at(md1_idx), data.md_anchor_z->at(md1_idx));
        return vPixel.Angle(vStrip); 
    }

} // End namespace extrapolation
