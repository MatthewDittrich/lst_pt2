#include "extrapolation.h"
#include "rootReader.h"
// ROOT Includes
#include <TVector3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

#include <cmath>
#include <iostream>

namespace extrapolation {

    // =========================================================================
    //                            HELPER FUNCTIONS
    // =========================================================================

        ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double R) {
        // Geometric solution: Speed over numerical minimization
        if (hits.size() < 2) return {0, 0};

        // Use the first and last hits (furthest points) to define the chord
        ROOT::Math::XYVector p1 = hits.front();
        ROOT::Math::XYVector p2 = hits.back();

        double dx = p2.X() - p1.X();
        double dy = p2.Y() - p1.Y();
        double d2 = dx*dx + dy*dy;
        double d = std::sqrt(d2);

        // Safety: If chord is longer than diameter, particle scattered or pT is wrong
        if (d > 2.0 * R || d == 0) return {0, 0};

        // Midpoint of the chord
        double midX = (p1.X() + p2.X()) / 2.0;
        double midY = (p1.Y() + p2.Y()) / 2.0;

        // Distance from midpoint to circle center
        double h = std::sqrt(std::max(0.0, R*R - d2/4.0));

        // Two possible centers perpendicular to the chord
        double c1x = midX + h * (-dy / d);
        double c1y = midY + h * (dx / d);
        double c2x = midX - h * (-dy / d);
        double c2y = midY - h * (dx / d);

        // Use a middle hit (if exists) to pick the correct center
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

    // 1. Collect valid PLS hits
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

    // 2. Helix parameters from PLS reco pT
        float reco_pt = data.pls_pt->at(pls_idx);
        const double B_FIELD_TESLA = 3.8;
        double R_pred_cm = (reco_pt / (0.3 * B_FIELD_TESLA)) * 100.0;

    // Assumes fitCircleWithFixedRadius is defined in your tools
        ROOT::Math::XYVector circle_center = fitCircleWithFixedRadius(pls_hits_xy, R_pred_cm);
        if (circle_center.X() == 0 && circle_center.Y() == 0) return invalid_return;

    // 3. Fit Z vs Angle (Helical Pitch)
        TGraph z_vs_angle_graph;
        double last_angle = 0;
        bool first_point = true;
        for (size_t i = 0; i < pls_hits_xy.size(); ++i) {
            double angle = TMath::ATan2(pls_hits_xy[i].Y() - circle_center.Y(), pls_hits_xy[i].X() - circle_center.X());
            if (!first_point) {
                while (angle - last_angle > TMath::Pi()) angle -= 2 * TMath::Pi();
                while (angle - last_angle < -TMath::Pi()) angle += 2 * TMath::Pi();
            }
            z_vs_angle_graph.SetPoint(z_vs_angle_graph.GetN(), angle, pls_hits_3d[i].Z());
            last_angle = angle;
            first_point = false;
        }

        if (z_vs_angle_graph.GetN() < 2) return invalid_return;
        z_vs_angle_graph.Fit("pol1", "Q");
        TF1* fit_func = z_vs_angle_graph.GetFunction("pol1");
        if (!fit_func) return invalid_return;

        double b = fit_func->GetParameter(0);
        double a = fit_func->GetParameter(1); // z = a*phi + b

    // 4. Get Target LS Positions (Mini-doublets)
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
        // Simple unwrap
            double phi_last = z_vs_angle_graph.GetX()[z_vs_angle_graph.GetN()-1];
            while (phi0 - phi_last >  TMath::Pi()) phi0 -= 2*TMath::Pi();
            while (phi0 - phi_last < -TMath::Pi()) phi0 += 2*TMath::Pi();

            double phi_fine = helixClosestPhiToPoint(target.X(), target.Y(), target.Z(), circle_center.X(), circle_center.Y(), R_pred_cm, a, b, phi0);
            TVector3 helix_pt(circle_center.X() + R_pred_cm * std::cos(phi_fine),
                              circle_center.Y() + R_pred_cm * std::sin(phi_fine),
                              a * phi_fine + b);
            return (target - helix_pt).Mag();
        };

        return {get_dist(md0_anchor), get_dist(md1_anchor)};
    }

}
