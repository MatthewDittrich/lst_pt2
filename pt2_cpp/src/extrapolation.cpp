#include "extrapolation.h"
#include "rootReader.h"
// ROOT Includes
#include <TVector3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

<<<<<<< Updated upstream
=======
<<<<<<< Updated upstream
=======
//#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>
>>>>>>> Stashed changes
>>>>>>> Stashed changes
#include <cmath>
#include <iostream>

namespace extrapolation {

    // =========================================================================
    //                            HELPER FUNCTIONS
    // =========================================================================

<<<<<<< Updated upstream
        ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double R) {
=======
<<<<<<< Updated upstream
        ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double R) {
=======
       ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double R) {
>>>>>>> Stashed changes
>>>>>>> Stashed changes
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

<<<<<<< Updated upstream
    
=======
<<<<<<< Updated upstream
    
=======
   
  /*     ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double fixed_radius) {
        if (hits.size() < 2) return {0, 0};

        // 1. Define the Chi-Squared function to minimize
        // p[0] is center_x (h), p[1] is center_y (k)
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

        // 2. Set up the Minuit2 Minimizer
        ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kMigrad);
        min.SetMaxFunctionCalls(10000);
        min.SetMaxIterations(1000);
        min.SetTolerance(0.001);
        min.SetPrintLevel(-1); // Silent mode

        ROOT::Math::Functor f(chi_squared_func, 2);
        min.SetFunction(f);

        // 3. Set starting point (average position of hits)
        double avg_x = 0.0, avg_y = 0.0;
        for (const auto& hit : hits) {
            avg_x += hit.X();
            avg_y += hit.Y();
        }
        avg_x /= hits.size();
        avg_y /= hits.size();

        min.SetVariable(0, "h", avg_x, 0.1);
        min.SetVariable(1, "k", avg_y, 0.1);

        // 4. Minimize
        min.Minimize();

        if (min.Status() == 0) {
            const double *result = min.X();
            return {result[0], result[1]};
        }

        // Return 0,0 if fit failed
        return {0, 0};
    }*/

>>>>>>> Stashed changes
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream

}
=======
<<<<<<< Updated upstream

}
=======
std::pair<double, double> extrapolatePlsInRZAndGetDeltaR(int pls_idx, int ls_idx, const rootReader& data) {
    const float invalid_val = -999.0;
    const std::pair<double, double> invalid_pair = {-999.0, -999.0};

    // 1. Collect PLS hits and convert to Z-R plane
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

    // 2. Linear Regression on PLS hits: R = slope * Z + intercept
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

    // 3. Get LS (Line Segment) Z positions
    int md0 = data.ls_mdIdx0->at(ls_idx);
    int md1 = data.ls_mdIdx1->at(ls_idx);
    if (md0 < 0 || md1 < 0) return invalid_pair;

    double z_ls0 = data.md_anchor_z->at(md0);
    double r_ls0_act = std::sqrt(std::pow(data.md_anchor_x->at(md0), 2) + std::pow(data.md_anchor_y->at(md0), 2));

    double z_ls1 = data.md_anchor_z->at(md1);
    double r_ls1_act = std::sqrt(std::pow(data.md_anchor_x->at(md1), 2) + std::pow(data.md_anchor_y->at(md1), 2));

    // 4. Predict R at LS Z positions and get residuals
    double r_ls0_pred = slope * z_ls0 + intercept;
    double r_ls1_pred = slope * z_ls1 + intercept;

    return {r_ls0_act - r_ls0_pred, r_ls1_act - r_ls1_pred};
}
    
    std::pair<double, double> extrapolateSimplePointingInRZ(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -999.0;
        const std::pair<double, double> invalid_pair = {-999.0, -999.0};

    // 1. Find "First" and "Last" valid PLS hits to define the slope
        double z_first = 0, r_first = 0;
        double z_last = 0, r_last = 0;
        bool found_first = false;
        bool found_last = false;

    // Identify First (Inner) Hit
        if (data.pls_hit0_x->at(pls_idx) > invalid_val) {
            z_first = data.pls_hit0_z->at(pls_idx);
            r_first = std::sqrt(std::pow(data.pls_hit0_x->at(pls_idx), 2) + std::pow(data.pls_hit0_y->at(pls_idx), 2));
            found_first = true;
        } else if (data.pls_hit1_x->at(pls_idx) > invalid_val) {
            z_first = data.pls_hit1_z->at(pls_idx);
            r_first = std::sqrt(std::pow(data.pls_hit1_x->at(pls_idx), 2) + std::pow(data.pls_hit1_y->at(pls_idx), 2));
            found_first = true;
        }

    // Identify Last (Outer) Hit
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

    // 2. Define the line: R = mZ + b
        double m = (r_last - r_first) / (z_last - z_first);
        double b = r_first - m * z_first;

    // 3. Get LS hits (MD0 and MD1)
        int md0 = data.ls_mdIdx0->at(ls_idx);
        int md1 = data.ls_mdIdx1->at(ls_idx);
        if (md0 < 0 || md1 < 0) return invalid_pair;

    // Residual for MD0
        double z0_act = data.md_anchor_z->at(md0);
        double r0_act = std::sqrt(std::pow(data.md_anchor_x->at(md0), 2) + std::pow(data.md_anchor_y->at(md0), 2));
        double dr0 = r0_act - (m * z0_act + b);

    // Residual for MD1
        double z1_act = data.md_anchor_z->at(md1);
        double r1_act = std::sqrt(std::pow(data.md_anchor_x->at(md1), 2) + std::pow(data.md_anchor_y->at(md1), 2));
        double dr1 = r1_act - (m * z1_act + b);

        return {dr0, dr1};
        }

       
        double calculateDeltaAngle(int pls_idx, int ls_idx, const rootReader& data) {
            const float inv = -999.0;
    
    // 1. Get First PLS Hit (Inner)
            TVector3 vPixel;
            if (data.pls_hit0_x->at(pls_idx) > inv) {
                vPixel.SetXYZ(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
            } else {
                return -999.0; // Essential hit missing
            }

    // 2. Get Last LS Hit (MD1 - Outer)
            int md1_idx = data.ls_mdIdx1->at(ls_idx);
            if (md1_idx < 0) return -999.0;

            TVector3 vStrip(data.md_anchor_x->at(md1_idx), data.md_anchor_y->at(md1_idx), data.md_anchor_z->at(md1_idx));

    // 3. Return the 3D angle between the two vectors (in radians)
            return vPixel.Angle(vStrip); 
        } 
    }
>>>>>>> Stashed changes
>>>>>>> Stashed changes
