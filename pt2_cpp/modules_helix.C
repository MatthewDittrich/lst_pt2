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
#include <TMatrixD.h>
#include <TVectorD.h>
#include <Math/Vector2D.h>
#include <cmath> 

// Structs
struct PlsMatchStatus { int used_count = 0; int unused_count = 0; };
struct LsMatchStatus  { int used_count = 0; int unused_count = 0; };

//Define how pt2 delta values are defined
class pT2Object {
public:
    float delta_pt, delta_eta, delta_phi, delta_R;
    pT2Object(int p_idx, int l_idx, const std::vector<float>& p_pt, const std::vector<float>& p_eta, const std::vector<float>& p_phi, const std::vector<float>& l_pt, const std::vector<float>& l_eta, const std::vector<float>& l_phi) {
        float pls_pt = p_pt.at(p_idx), pls_eta = p_eta.at(p_idx), pls_phi = p_phi.at(p_idx);
        float ls_pt = l_pt.at(l_idx), ls_eta = l_eta.at(l_idx), ls_phi = l_phi.at(l_idx);
        delta_pt = pls_pt - ls_pt; 
        delta_eta = pls_eta - ls_eta;
        delta_phi = TVector2::Phi_mpi_pi(pls_phi - ls_phi);
        delta_R = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
    }
};

////////////////////// EXTRAPOLATION//////////////////////

//-------------------------HELLICAL EXTRAPOLATION----------------------------------------------
//---------Fitting of circle using pls hits---------------
std::pair<ROOT::Math::XYVector, double> fitCircle(const std::vector<ROOT::Math::XYVector>& hits) {
    if (hits.size() < 3) {
        return {ROOT::Math::XYVector(0, 0), -1.0};
    }

    // Karimäki circle fit (minimizes bias for arcs in tracking)
    double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxxx = 0, Syyy = 0, Sxxy = 0, Sxyy = 0;
    int n = hits.size();

    for (const auto& h : hits) {
        double x = h.X();
        double y = h.Y();
        Sx   += x;
        Sy   += y;
        Sxx  += x*x;
        Syy  += y*y;
        Sxy  += x*y;
        Sxxx += x*x*x;
        Syyy += y*y*y;
        Sxxy += x*x*y;
        Sxyy += x*y*y;
    }

    double C = n*Sxx - Sx*Sx;
    double D = n*Sxy - Sx*Sy;
    double E = n*Syy - Sy*Sy;
    double G = 0.5*(n*(Sxxx + Sxyy) - Sx*(Sxx + Syy));
    double H = 0.5*(n*(Syyy + Sxxy) - Sy*(Sxx + Syy));

    double denom = C*E - D*D;
    if (std::abs(denom) < 1e-12) {
        // Points are nearly collinear
        return {ROOT::Math::XYVector(0, 0), -1.0};
    }

    double cx = (G*E - D*H) / denom;
    double cy = (C*H - D*G) / denom;

    // Radius: average distance to center
    double r_sum = 0.0;
    for (const auto& h : hits) {
        r_sum += std::hypot(h.X() - cx, h.Y() - cy);
    }
    double R = r_sum / n;

    return {ROOT::Math::XYVector(cx, cy), R};
}

ROOT::Math::XYVector fitCircleWithFixedRadius(
    const std::vector<ROOT::Math::XYVector>& hits, 
    double fixed_radius
) {
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


//-------This function uses the circle fitting to create the hellical trajectory---------
std::pair<double, double> extrapolatePlsHelicallyAndGetDistance(
    int pls_idx, int ls_idx, const std::vector<float>& pls_pt,
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
    const float invalid_hit_val = -999.0f;
    const std::pair<double, double> invalid_return = {-1.0, -1.0};

    // --- Gather PLS hits ---
    std::vector<ROOT::Math::XYVector> pls_hits_xy;
    std::vector<TVector3> pls_hits_3d;
    auto add_hit = [&](float x, float y, float z) {
        if (x > invalid_hit_val) {
            pls_hits_xy.emplace_back(x, y);
            pls_hits_3d.emplace_back(x, y, z);
        }
    };
    add_hit(pls_hit0_x.at(pls_idx), pls_hit0_y.at(pls_idx), pls_hit0_z.at(pls_idx));
    add_hit(pls_hit1_x.at(pls_idx), pls_hit1_y.at(pls_idx), pls_hit1_z.at(pls_idx));
    add_hit(pls_hit2_x.at(pls_idx), pls_hit2_y.at(pls_idx), pls_hit2_z.at(pls_idx));
    add_hit(pls_hit3_x.at(pls_idx), pls_hit3_y.at(pls_idx), pls_hit3_z.at(pls_idx));

    if (pls_hits_xy.size() < 3) return invalid_return;

    // --- Fit circle (Karimäki) ---
    float reco_pt = pls_pt.at(pls_idx);

    // Step B: Calculate the radius predicted by the reconstructed pT
    const double B_FIELD_TESLA = 3.8; // IMPORTANT: Set your detector's value
    double R_pred_cm = (reco_pt / (0.3 * B_FIELD_TESLA)) * 100.0;
    
    // Step C: Call the new fitter to find the center using the fixed radius
    ROOT::Math::XYVector circle_center = fitCircleWithFixedRadius(pls_hits_xy, R_pred_cm);

    // Step D: Check if the fit was successful and assign the final radius
    if (circle_center.X() == 0 && circle_center.Y() == 0) {
        return invalid_return; // The constrained fit failed
    }
    double circle_radius = R_pred_cm;

    // --- Build z vs angle graph ---
    TGraph z_vs_angle_graph;
    double last_angle = 0;
    bool first_point = true;
    for (const auto& hit_xy : pls_hits_xy) {
        double angle = TMath::ATan2(hit_xy.Y() - circle_center.Y(), hit_xy.X() - circle_center.X());
        if (!first_point) {
            while (angle - last_angle > TMath::Pi()) angle -= 2 * TMath::Pi();
            while (angle - last_angle < -TMath::Pi()) angle += 2 * TMath::Pi();
        }
        // find the matching 3D hit
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

    // Parameters for z(phi) = b + a * phi
    double b = fit_func->GetParameter(0); // intercept
    double a = fit_func->GetParameter(1); // slope

    // LS -> MD indices
    int md_idx0 = ls_mdIdx0.at(ls_idx);
    int md_idx1 = ls_mdIdx1.at(ls_idx);
    if (md_idx0 < 0 || md_idx1 < 0) return invalid_return;

    TVector3 md0_anchor(md_anchor_x.at(md_idx0), md_anchor_y.at(md_idx0), md_anchor_z.at(md_idx0));
    TVector3 md1_anchor(md_anchor_x.at(md_idx1), md_anchor_y.at(md_idx1), md_anchor_z.at(md_idx1));

    // Helper: unwrap phi to be continuous with last fitted angle
    auto unwrap_to_last = [&](double phi_target) {
        const int N = z_vs_angle_graph.GetN();
        const double phi_last = z_vs_angle_graph.GetX()[N-1];
        while (phi_target - phi_last >  TMath::Pi())  phi_target -= 2*TMath::Pi();
        while (phi_target - phi_last < -TMath::Pi())  phi_target += 2*TMath::Pi();
        return phi_target;
    };

    // Newton solver to find phi that minimizes D^2(phi) to anchor (closest approach)
    auto helixClosestPhiToPoint = [&](double xa, double ya, double za, double xc, double yc, double R, double a_param, double b_param, double phi_init) {
        double phi = phi_init;
        for (int iter = 0; iter < 12; ++iter) {
            double c = std::cos(phi);
            double s = std::sin(phi);
            double dx = xa - xc - R * c; // xa - x(phi)
            double dy = ya - yc - R * s; // ya - y(phi)
            double dz = a_param * phi + b_param - za; // z(phi) - za

            // f(phi) = dx*(R s) + dy*(-R c) + a*(a*phi + b - za)
            double f = dx * (R * s) + dy * (-R * c) + a_param * dz;

            // f'(phi) = R^2 + R*(dx*c + dy*s) + a^2
            double fp = R*R + R*(dx * c + dy * s) + a_param * a_param;

            if (std::abs(fp) < 1e-12) break;
            double step = f / fp;
            phi -= step;

            // clamp big steps to keep stability
            if (std::abs(step) > 1.0) {
                phi += (step > 0) ? 1.0 : -1.0; // back off a bit
            }
            if (std::abs(step) < 1e-10) break; // converged
        }
        return phi;
    };

    // For each anchor: initial guess = anchor azimuth unwrapped to last fit
    auto initial_phi = [&](double xa, double ya) {
        double phi0 = std::atan2(ya - circle_center.Y(), xa - circle_center.X());
        return unwrap_to_last(phi0);
    };

    // Anchor A (md0)
    double phiA0 = initial_phi(md0_anchor.X(), md0_anchor.Y());
    double phiA = helixClosestPhiToPoint(md0_anchor.X(), md0_anchor.Y(), md0_anchor.Z(),
                                         circle_center.X(), circle_center.Y(), circle_radius,
                                         a, b, phiA0);
    double xA = circle_center.X() + circle_radius * std::cos(phiA);
    double yA = circle_center.Y() + circle_radius * std::sin(phiA);
    double zA = a * phiA + b;
    TVector3 helixA(xA, yA, zA);
    double dist0 = (md0_anchor - helixA).Mag();

    // Anchor B (md1)
    double phiB0 = initial_phi(md1_anchor.X(), md1_anchor.Y());
    double phiB = helixClosestPhiToPoint(md1_anchor.X(), md1_anchor.Y(), md1_anchor.Z(),
                                         circle_center.X(), circle_center.Y(), circle_radius,
                                         a, b, phiB0);
    double xB = circle_center.X() + circle_radius * std::cos(phiB);
    double yB = circle_center.Y() + circle_radius * std::sin(phiB);
    double zB = a * phiB + b;
    TVector3 helixB(xB, yB, zB);
    double dist1 = (md1_anchor - helixB).Mag();

    return {dist0, dist1};
}

//--------------Z-R EXTRAPOLATION------------------------
//------------using ls hits--------------
//Extrapolation with coordinates r-z, uses ls hits instead of pls to do the fitting
std::vector<double> extrapolateLsInRZAndGetDeltaR_AllHits(
    int pls_idx, int ls_idx,
    // ALL PLS hit position branches are now needed
    const std::vector<float>& pls_hit0_x, const std::vector<float>& pls_hit0_y, const std::vector<float>& pls_hit0_z,
    const std::vector<float>& pls_hit1_x, const std::vector<float>& pls_hit1_y, const std::vector<float>& pls_hit1_z,
    const std::vector<float>& pls_hit2_x, const std::vector<float>& pls_hit2_y, const std::vector<float>& pls_hit2_z,
    const std::vector<float>& pls_hit3_x, const std::vector<float>& pls_hit3_y, const std::vector<float>& pls_hit3_z,
    // LS -> MD link branches
    const std::vector<int>& ls_mdIdx0, const std::vector<int>& ls_mdIdx1,
    // MD position branches
    const std::vector<float>& md_anchor_x, const std::vector<float>& md_anchor_y, const std::vector<float>& md_anchor_z
) {
    const float invalid_val = -999.0;
    std::vector<double> results; // We will return a vector of delta_r values

    // --- Part 1: Get the two (z,r) points from the LS's mini-doublets ---
    int md_idx0 = ls_mdIdx0.at(ls_idx);
    int md_idx1 = ls_mdIdx1.at(ls_idx);
    if (md_idx0 < 0 || md_idx1 < 0) return results; // Return empty vector on failure

    auto to_zr = [](float x, float y, float z) { return TVector2(z, std::sqrt(x*x + y*y)); };
    TVector2 p1_zr = to_zr(md_anchor_x.at(md_idx0), md_anchor_y.at(md_idx0), md_anchor_z.at(md_idx0));
    TVector2 p2_zr = to_zr(md_anchor_x.at(md_idx1), md_anchor_y.at(md_idx1), md_anchor_z.at(md_idx1));
    if (std::abs(p2_zr.X() - p1_zr.X()) < 1e-9) return results; // Return empty vector if vertical

    // --- Part 2: Create the R-Z extrapolation model (line equation) from the LS ---
    double slope = (p2_zr.Y() - p1_zr.Y()) / (p2_zr.X() - p1_zr.X());
    double intercept = p1_zr.Y() - slope * p1_zr.X();

    // --- Part 3: Create a vector of all valid PLS hits to test against ---
    std::vector<TVector2> target_pls_hits_zr;
    if (pls_hit0_x.at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(pls_hit0_x.at(pls_idx), pls_hit0_y.at(pls_idx), pls_hit0_z.at(pls_idx)));
    if (pls_hit1_x.at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(pls_hit1_x.at(pls_idx), pls_hit1_y.at(pls_idx), pls_hit1_z.at(pls_idx)));
    if (pls_hit2_x.at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(pls_hit2_x.at(pls_idx), pls_hit2_y.at(pls_idx), pls_hit2_z.at(pls_idx)));
    if (pls_hit3_x.at(pls_idx) > invalid_val) target_pls_hits_zr.push_back(to_zr(pls_hit3_x.at(pls_idx), pls_hit3_y.at(pls_idx), pls_hit3_z.at(pls_idx)));

    // --- Part 4: Loop through each PLS hit, calculate delta_r, and store it ---
    for (const auto& target_hit : target_pls_hits_zr) {
        double target_z = target_hit.X();
        double target_r_actual = target_hit.Y();
        double predicted_r = slope * target_z + intercept;
        double delta_r = target_r_actual - predicted_r;
        results.push_back(delta_r);
    }

    return results;
}

//------------------ This module also uses R-Z coordinates just fins delta z instead of delta r----------------

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
//////////////////////////////////////////////////////////////////////////////////////////////////////


//---------------MODULES TO FIND DIFFERENT OBJECTS: USED, UNUSED, REAL, FAKE, Pt2_objects------------------------------------
// MODULE 1: Find Used Indices
void findUsedIndices(
    // Output sets to be filled
    std::set<int>& used_ls_indices, 
    std::set<int>& used_pls_indices, 

    // Input data vectors (indices linking different objects)
    const std::vector<int>& tc_pt5Idx_vec, 
    const std::vector<int>& pt5_t5Idx_vec, 
    const std::vector<int>& t5_t3Idx0_vec, 
    const std::vector<int>& t5_t3Idx1_vec, 
    const std::vector<int>& t3_lsIdx0_vec, 
    const std::vector<int>& t3_lsIdx1_vec, 
    const std::vector<int>& tc_pt3Idx_vec, 
    const std::vector<int>& pt3_t3Idx_vec, 
    const std::vector<int>& tc_t5Idx_vec, 
    const std::vector<int>& pt5_plsIdx_vec, 
    const std::vector<int>& pt3_plsIdx_vec
) {
    // --- Loop 1: Process Track Candidates linked to pT5 objects ---
    for (int pt5_Idx : tc_pt5Idx_vec) {
        if (pt5_Idx == -999) continue;

        // Mark the associated PLS as "used"
        if (pt5_plsIdx_vec.at(pt5_Idx) != -999) {
            used_pls_indices.insert(pt5_plsIdx_vec.at(pt5_Idx));
        }

        // Follow the link from pT5 to its T5
        int t5_Idx = pt5_t5Idx_vec.at(pt5_Idx);
        if (t5_Idx != -999) {
            // Follow the link from T5 to its first T3
            int t3_0 = t5_t3Idx0_vec.at(t5_Idx);
            if (t3_0 != -999) {
                // Mark the two LSs associated with this T3 as "used"
                if (t3_lsIdx0_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_0));
                if (t3_lsIdx1_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_0));
            }

            // Follow the link from T5 to its second T3
            int t3_1 = t5_t3Idx1_vec.at(t5_Idx);
            if (t3_1 != -999) {
                // Mark the two LSs associated with this T3 as "used"
                if (t3_lsIdx0_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_1));
                if (t3_lsIdx1_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_1));
            }
        }
    }

    // --- Loop 2: Process Track Candidates linked to pT3 objects ---
    for (int pt3_Idx : tc_pt3Idx_vec) {
        if (pt3_Idx == -999) continue;
        
        // Mark the associated PLS as "used"
        if (pt3_plsIdx_vec.at(pt3_Idx) != -999) {
            used_pls_indices.insert(pt3_plsIdx_vec.at(pt3_Idx));
        }
        
        // Follow the link from pT3 to its T3
        int t3_idx = pt3_t3Idx_vec.at(pt3_Idx);
        if (t3_idx != -999) {
            // Mark the two LSs associated with this T3 as "used"
            if (t3_lsIdx0_vec.at(t3_idx) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_idx));
            if (t3_lsIdx1_vec.at(t3_idx) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_idx));
        }
    }

    // --- Loop 3: Process Track Candidates linked directly to T5 objects ---
    for (int t5_Idx : tc_t5Idx_vec) {
        if (t5_Idx == -999) continue;
        
        // Follow the link from T5 to its first T3
        int t3_0 = t5_t3Idx0_vec.at(t5_Idx);
        if (t3_0 != -999) {
            // Mark the two LSs associated with this T3 as "used"
            if (t3_lsIdx0_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_0));
            if (t3_lsIdx1_vec.at(t3_0) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_0));
        }
        
        // Follow the link from T5 to its second T3
        int t3_1 = t5_t3Idx1_vec.at(t5_Idx);
        if (t3_1 != -999) {
            // Mark the two LSs associated with this T3 as "used"
            if (t3_lsIdx0_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx0_vec.at(t3_1));
            if (t3_lsIdx1_vec.at(t3_1) != -999) used_ls_indices.insert(t3_lsIdx1_vec.at(t3_1));
        }
    }
}

// MODULE 2: Find PLS Match Status
void getPlsMatchStatus(
    std::map<int, PlsMatchStatus>& simIdx_to_pls_status, 
    const std::vector<int>& pls_isFake_vec, 
    const std::vector<int>& pls_simIdx_vec, 
    const std::set<int>& used_pls_indices
) {
    // Loop over every PLS object in the event
    for (size_t j = 0; j < pls_simIdx_vec.size(); ++j) {
        
        // Step 1: Process only "real" (not fake) PLS objects
        if (pls_isFake_vec.at(j) == 0) {
            
            int sim_idx = pls_simIdx_vec.at(j);

            // Step 2: Ensure the PLS has a valid link to a simulated particle
            if (sim_idx != -999) {
                
                // Step 3: Check if this PLS was used in a track candidate
                if (used_pls_indices.count(j) > 0) {
                    // If yes, increment the 'used' counter for this sim particle
                    simIdx_to_pls_status[sim_idx].used_count++;
                } else {
                    // If no, increment the 'unused' counter for this sim particle
                    simIdx_to_pls_status[sim_idx].unused_count++;
                }
            }
        }
    }
}

// MODULE 3: Find LS Match Status
void getLsMatchStatus(
    std::map<int, LsMatchStatus>& simIdx_to_ls_status, 
    const std::vector<int>& ls_isFake_vec, 
    const std::vector<int>& ls_simIdx_vec, 
    const std::set<int>& used_ls_indices
) {
    // Loop over every Line Segment (LS) object in the event
    for (size_t j = 0; j < ls_simIdx_vec.size(); ++j) {
        
        // Step 1: Process only "real" (not fake) LS objects
        if (ls_isFake_vec.at(j) == 0) {
            
            int sim_idx = ls_simIdx_vec.at(j);

            // Step 2: Ensure the LS has a valid link to a simulated particle
            if (sim_idx != -999) {
                
                // Step 3: Check if this LS was used in a track candidate
                if (used_ls_indices.count(j) > 0) {
                    // If yes, increment the 'used' counter for this sim particle
                    simIdx_to_ls_status[sim_idx].used_count++;
                } else {
                    // If no, increment the 'unused' counter for this sim particle
                    simIdx_to_ls_status[sim_idx].unused_count++;
                }
            }
        }
    }
}



// MODULE 4: Process Ideal pT2 objects 
void processIdealPT2s(const std::vector<float>& pls_pt, const std::vector<float>& pls_eta, const std::vector<float>& pls_phi, 
     const std::vector<int>& pls_isFake, const std::vector<int>& pls_simIdx, 
     const std::vector<float>& ls_pt, const std::vector<float>& ls_eta, const std::vector<float>& ls_phi, 
     const std::vector<int>& ls_isFake, const std::vector<int>& ls_simIdx, const std::set<int>& used_pls_indices, const std::set<int>& used_ls_indices, 
     TH1D* h_delta_pt, TH1D* h_delta_eta, TH1D* h_delta_phi, TH1D* h_delta_R, TH1D* h_extrapolation_dist_3d, TH1D* h_extrapolation_delta_z,
     std::vector<TH1D*>& h_delta_r_by_layer, std::vector<TH1D*>& h_extrapolation_delta_r_reverse_by_hit, TH1D* h_extrapolation_delta_r_reverse_combined,
     const std::vector<float>& pls_hit0_x, const std::vector<float>& pls_hit0_y, const std::vector<float>& pls_hit0_z, 
     const std::vector<float>& pls_hit1_x, const std::vector<float>& pls_hit1_y, const std::vector<float>& pls_hit1_z, 
     const std::vector<float>& pls_hit2_x, const std::vector<float>& pls_hit2_y, const std::vector<float>& pls_hit2_z, 
     const std::vector<float>& pls_hit3_x, const std::vector<float>& pls_hit3_y, const std::vector<float>& pls_hit3_z, 
     const std::vector<int>& ls_mdIdx0, const std::vector<int>& ls_mdIdx1, 
     const std::vector<float>& md_anchor_x, const std::vector<float>& md_anchor_y, const std::vector<float>& md_anchor_z, 
     const std::vector<float>& md_other_x, const std::vector<float>& md_other_y, const std::vector<float>& md_other_z,
     const std::vector<int>& md_layer, long long& ideal_pt2_count, const std::vector<float>& sim_pt, std::vector<TH1D*>& h_dist_3d_by_layer) 
{
    std::map<int, std::vector<int>> ls_simIdx_map;
    for (size_t k = 0; k < ls_pt.size(); ++k) { if (ls_isFake.at(k) == 0 && ls_simIdx.at(k) != -999) { ls_simIdx_map[ls_simIdx.at(k)].push_back(k); } }
    
    for (size_t i = 0; i < pls_pt.size(); ++i) { 
    	if (pls_isFake.at(i) == 0 && used_pls_indices.count(i) == 0 && pls_simIdx.at(i) != -999) { 
    	   auto it = ls_simIdx_map.find(pls_simIdx.at(i)); 
     	   if (it != ls_simIdx_map.end()) { 
     	       for (int ls_idx : it->second) { 
     		   if (used_ls_indices.count(ls_idx) == 0) { 
     		    ideal_pt2_count++;
     		    pT2Object ideal_pt2(i, ls_idx, pls_pt, pls_eta, pls_phi, ls_pt, ls_eta, ls_phi); 
    		    h_delta_pt->Fill(ideal_pt2.delta_pt); 
    		    h_delta_eta->Fill(ideal_pt2.delta_eta); 
    		    h_delta_phi->Fill(ideal_pt2.delta_phi); 
    		    h_delta_R->Fill(ideal_pt2.delta_R); 
    		
    		float reco_pt_pls = pls_pt.at(i);
                        if (reco_pt_pls <= 2.0) {	
                            
                            // --- New Truth Matching Filter ---
                            int sim_idx = pls_simIdx.at(i);
                            // Safety check for valid index
                            if (sim_idx >= 0 && sim_idx < sim_pt.size()) {
                                float true_pt = sim_pt.at(sim_idx);
                                const float pT_tolerance = 0.1; // 100 MeV

                                // Check if the pT matches the truth
                              //  if (std::abs(true_pt - reco_pt_pls) <= pT_tolerance) {
                                  if (true_pt<= 2.0) {  
                                    // If it matches, perform the extrapolation
                                    auto distance_3d = extrapolatePlsHelicallyAndGetDistance(
                                        i, ls_idx, pls_pt,
                                        pls_hit0_x, pls_hit0_y, pls_hit0_z, pls_hit1_x, pls_hit1_y, pls_hit1_z, 
                                        pls_hit2_x, pls_hit2_y, pls_hit2_z, pls_hit3_x, pls_hit3_y, pls_hit3_z, 
                                        ls_mdIdx0, ls_mdIdx1, 
                                        md_anchor_x, md_anchor_y, md_anchor_z, 
                                        md_other_x, md_other_y, md_other_z
                                    ); 
                                    
                                    // And fill the distance histogram
                                    if (distance_3d.first >= 0) {
                                        h_extrapolation_dist_3d->Fill(distance_3d.first);
                                        int md_idx0 = ls_mdIdx0.at(ls_idx);
                                        if (md_idx0 >= 0) { // Safety check
                                           int layer0 = md_layer.at(md_idx0);
                                           if (layer0 < h_dist_3d_by_layer.size()) {
                			        h_dist_3d_by_layer[layer0]->Fill(distance_3d.first);
                                    }}}
                                    if (distance_3d.second >= 0) {
                                        h_extrapolation_dist_3d->Fill(distance_3d.second);
                                        int md_idx1 = ls_mdIdx1.at(ls_idx);
        				  if (md_idx1 >= 0) { // Safety check
            				      int layer1 = md_layer.at(md_idx1);
            				      if (layer1 < h_dist_3d_by_layer.size()) {
                				  h_dist_3d_by_layer[layer1]->Fill(distance_3d.second);
            }
        }
                                        
                                    }
                                } // end of truth match if
                            } // end of safety check if
                        }
    		//USES LS and adds each hit
    		std::vector<double> delta_r_results = extrapolateLsInRZAndGetDeltaR_AllHits(
    		i, ls_idx,
    		pls_hit0_x, pls_hit0_y, pls_hit0_z,
    		pls_hit1_x, pls_hit1_y, pls_hit1_z,
    		pls_hit2_x, pls_hit2_y, pls_hit2_z,
    		pls_hit3_x, pls_hit3_y, pls_hit3_z,
    		ls_mdIdx0, ls_mdIdx1,
    		md_anchor_x, md_anchor_y, md_anchor_z);
    		
    		for (size_t hit_num = 0; hit_num < delta_r_results.size(); ++hit_num) {
    		     if (hit_num < h_extrapolation_delta_r_reverse_by_hit.size()) {
    		     h_extrapolation_delta_r_reverse_by_hit[hit_num]->Fill(delta_r_results[hit_num]);
    		     h_extrapolation_delta_r_reverse_combined->Fill(delta_r_results[hit_num]);
    			}
		}
    		///uses ls hits////
    		
    		const float invalid_val = -999.0; 
    		
    		//if (delta_r_reverse > invalid_val) {
        //h_extrapolation_delta_r_reverse->Fill(delta_r_reverse);
    //}
    		///uses pls hits///
    		std::pair<double, double> delta_r_pair = extrapolatePlsInRZAndGetDeltaZ(i, ls_idx, pls_hit0_x, pls_hit0_y, pls_hit0_z, pls_hit1_x, pls_hit1_y, pls_hit1_z, pls_hit2_x, pls_hit2_y, pls_hit2_z, pls_hit3_x, pls_hit3_y, pls_hit3_z, ls_mdIdx0, ls_mdIdx1, md_anchor_x, md_anchor_y, md_anchor_z, md_other_x, md_other_y, md_other_z);

    		int md_idx0 = ls_mdIdx0.at(ls_idx);
		int md_idx1 = ls_mdIdx1.at(ls_idx);
 
    		if (delta_r_pair.first > invalid_val && md_idx0 >= 0) {
                            h_extrapolation_delta_z->Fill(delta_r_pair.first);
                            int layer0 = md_layer.at(md_idx0);
                            if (layer0 < h_delta_r_by_layer.size()) {
                                h_delta_r_by_layer[layer0]->Fill(delta_r_pair.first);
                            }
                        }
                if (delta_r_pair.second > invalid_val && md_idx1 >= 0) {
                	     h_extrapolation_delta_z->Fill(delta_r_pair.second);
                            int layer1 = md_layer.at(md_idx1);
                            if (layer1 < h_delta_r_by_layer.size()) {
                               h_delta_r_by_layer[layer1]->Fill(delta_r_pair.second); }} } } } }
}}


//---------------CREATING HISTOGRAMS------------------

// MODULE 1: Fill Detailed LS Histograms
void fillDetailedLsHistograms(
    // Properties of the current LS
    bool is_ls_used, 
    bool is_ls_fake, 
    float pt, 
    float eta, 
    float phi, 
    int ls_sim_idx, 

    // Pre-computed map of PLS matching statuses for the event
    const std::map<int, PlsMatchStatus>& simIdx_to_pls_status, 

    // Pointers to all the histograms to be filled
    // (f=fake, r=real, u=used, un=unused, um=unmatched, b=both, uo=used_only, uno=unused_only)
    TH1D* h_pt_f_u,      TH1D* h_pt_f_un,
    TH1D* h_pt_r_u_um,   TH1D* h_pt_r_u_b,   TH1D* h_pt_r_u_uo,   TH1D* h_pt_r_u_uno,
    TH1D* h_pt_r_un_um,  TH1D* h_pt_r_un_b,  TH1D* h_pt_r_un_uo,  TH1D* h_pt_r_un_uno,
    TH1D* h_eta_f_u,     TH1D* h_eta_f_un,
    TH1D* h_eta_r_u_um,  TH1D* h_eta_r_u_b,  TH1D* h_eta_r_u_uo,  TH1D* h_eta_r_u_uno,
    TH1D* h_eta_r_un_um, TH1D* h_eta_r_un_b, TH1D* h_eta_r_un_uo, TH1D* h_eta_r_un_uno,
    TH1D* h_phi_f_u,     TH1D* h_phi_f_un,
    TH1D* h_phi_r_u_um,  TH1D* h_phi_r_u_b,  TH1D* h_phi_r_u_uo,  TH1D* h_phi_r_u_uno,
    TH1D* h_phi_r_un_um, TH1D* h_phi_r_un_b, TH1D* h_phi_r_un_uo, TH1D* h_phi_r_un_uno
) {
    // Helper lambda to fill a set of (pt, eta, phi) histograms at once
    auto fill_all = [&](TH1D* h_pt, TH1D* h_eta, TH1D* h_phi) {
        h_pt->Fill(pt);
        h_eta->Fill(eta);
        h_phi->Fill(phi);
    };

    // --- CASE 1: The LS is FAKE ---
    if (is_ls_fake) {
        if (is_ls_used) {
            fill_all(h_pt_f_u, h_eta_f_u, h_phi_f_u);
        } else {
            fill_all(h_pt_f_un, h_eta_f_un, h_phi_f_un);
        }
        return; // Processing for this LS is done
    }

    // --- From here on, we know the LS is REAL ---

    // Find the matching PLS status for this LS's sim_idx.
    // If the LS has no valid sim_idx, the iterator will be .end()
    auto match_it = (ls_sim_idx == -999) 
                    ? simIdx_to_pls_status.end() 
                    : simIdx_to_pls_status.find(ls_sim_idx);

    // --- CASE 2: The REAL LS has NO MATCHING PLS object ---
    if (match_it == simIdx_to_pls_status.end()) {
        if (is_ls_used) {
            fill_all(h_pt_r_u_um, h_eta_r_u_um, h_phi_r_u_um);
        } else {
            fill_all(h_pt_r_un_um, h_eta_r_un_um, h_phi_r_un_um);
        }
    } 
    // --- CASE 3: The REAL LS HAS a MATCHING PLS object ---
    else {
        const auto& status = match_it->second;

        // Check the status of the matched PLS(s)
        if (status.used_count > 0 && status.unused_count > 0) {
            // Matched to BOTH used and unused PLS objects
            if (is_ls_used) fill_all(h_pt_r_u_b, h_eta_r_u_b, h_phi_r_u_b);
            else            fill_all(h_pt_r_un_b, h_eta_r_un_b, h_phi_r_un_b);
        }
        else if (status.used_count > 0) {
            // Matched to USED PLS objects ONLY
            if (is_ls_used) fill_all(h_pt_r_u_uo, h_eta_r_u_uo, h_phi_r_u_uo);
            else            fill_all(h_pt_r_un_uo, h_eta_r_un_uo, h_phi_r_un_uo);
        }
        else { // This implies status.unused_count > 0
            // Matched to UNUSED PLS objects ONLY
            if (is_ls_used) fill_all(h_pt_r_u_uno, h_eta_r_u_uno, h_phi_r_u_uno);
            else            fill_all(h_pt_r_un_uno, h_eta_r_un_uno, h_phi_r_un_uno);
        }
    }
}


// MODULE 2: Create and Save a Detailed Stacked Plot
void createAndSaveDetailedStack(
    const std::string& canvasName, 
    const std::string& stackTitle, 
    const std::string& saveFileName, 
    const std::string& legendHeader, 
    bool setLogY, 
    TH1D* h_fake, 
    TH1D* h_real_unmatched, 
    TH1D* h_real_matchedToUnused, 
    TH1D* h_real_matchedToBoth, 
    TH1D* h_real_matchedToUsed
) {
    // --- 1. Canvas and THStack Setup ---
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    if (setLogY) {
        canvas->SetLogy();
    }
    THStack *hs = new THStack("hs", stackTitle.c_str());

    // --- 2. Histogram Styling (Fill Colors) ---
    h_fake->SetFillColor(kGray);
    h_real_unmatched->SetFillColor(kRed - 7);
    h_real_matchedToUnused->SetFillColor(kAzure + 1);
    h_real_matchedToBoth->SetFillColor(kOrange + 1);
    h_real_matchedToUsed->SetFillColor(kGreen + 2);

    // --- 3. Stacking Histograms ---
    // The order matters: Add from top to bottom of the stack.
    hs->Add(h_real_matchedToUsed);
    hs->Add(h_real_matchedToBoth);
    hs->Add(h_real_matchedToUnused);
    hs->Add(h_real_unmatched);
    hs->Add(h_fake);

    // --- 4. Draw the Stack ---
    hs->Draw("HIST"); 

    // --- 5. Create and Draw the Legend ---
    auto legend = new TLegend(0.40, 0.58, 0.88, 0.88);
    legend->SetHeader(legendHeader.c_str(), "C"); // "C" centers the header text
    legend->AddEntry(h_real_matchedToUsed,   "Real, Matched to USED Only",    "f");
    legend->AddEntry(h_real_matchedToBoth,   "Real, Matched to Both",         "f");
    legend->AddEntry(h_real_matchedToUnused, "Real, Matched to UNUSED Only",  "f");
    legend->AddEntry(h_real_unmatched,       "Real, No Match",                "f");
    legend->AddEntry(h_fake,                 "Fake",                          "f");
    legend->SetBorderSize(0); // Removes the box outline around the legend
    legend->Draw();

    // --- 6. Save and Clean up ---
    canvas->SaveAs(saveFileName.c_str());
    delete canvas;
}


// MODULE 3: Create and Save a Simple 1D Plot
void createAndSaveSimplePlot(const std::string& canvasName, const std::string& saveFileName, TH1D* histo, const std::string& plotTitle, const std::string& xLabel, const std::string& yLabel) {
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 700);
    histo->SetTitle(plotTitle.c_str()); 
    histo->GetXaxis()->SetTitle(xLabel.c_str()); 
    histo->GetYaxis()->SetTitle(yLabel.c_str());
    histo->SetLineColor(kAzure+2); 
    histo->SetLineWidth(2); 
    histo->Draw("HIST");
    
     double integral = histo->Integral();

    // Create a text box (TPaveText) to display the info
    // The coordinates are in normalized screen space (0 to 1)
    TPaveText *pt = new TPaveText(0.65, 0.82, 0.88, 0.88, "NDC"); // NDC = Normalized Device Coords
    pt->SetFillColor(0); // Transparent background
    pt->SetBorderSize(0);
    pt->SetTextAlign(12); // Left-aligned, vertically centered
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);

    // Format the string and add it to the box
    TString text = TString::Format("Entries: %.0f", integral);
    pt->AddText(text);
    pt->Draw("SAME");
    
    canvas->SaveAs(saveFileName.c_str()); 
    delete canvas;
}


// MODULE 4: Fill Detailed PLS Histograms
void fillDetailedPlsHistograms(
    // Properties of the current PLS
    bool is_pls_used, 
    bool is_pls_fake, 
    float pt, 
    float eta, 
    float phi, 
    int pls_sim_idx, 

    // Pre-computed map of LS matching statuses for the event
    const std::map<int, LsMatchStatus>& simIdx_to_ls_status, 

    // Pointers to all the histograms to be filled
    // (f=fake, r=real, u=used, un=unused, um=unmatched, b=both, uo=used_only, uno=unused_only)
    TH1D* h_pt_f_u,      TH1D* h_pt_f_un,
    TH1D* h_pt_r_u_um,   TH1D* h_pt_r_u_b,   TH1D* h_pt_r_u_uo,   TH1D* h_pt_r_u_uno,
    TH1D* h_pt_r_un_um,  TH1D* h_pt_r_un_b,  TH1D* h_pt_r_un_uo,  TH1D* h_pt_r_un_uno,
    TH1D* h_eta_f_u,     TH1D* h_eta_f_un,
    TH1D* h_eta_r_u_um,  TH1D* h_eta_r_u_b,  TH1D* h_eta_r_u_uo,  TH1D* h_eta_r_u_uno,
    TH1D* h_eta_r_un_um, TH1D* h_eta_r_un_b, TH1D* h_eta_r_un_uo, TH1D* h_eta_r_un_uno,
    TH1D* h_phi_f_u,     TH1D* h_phi_f_un,
    TH1D* h_phi_r_u_um,  TH1D* h_phi_r_u_b,  TH1D* h_phi_r_u_uo,  TH1D* h_phi_r_u_uno,
    TH1D* h_phi_r_un_um, TH1D* h_phi_r_un_b, TH1D* h_phi_r_un_uo, TH1D* h_phi_r_un_uno
) {
    // Helper lambda to fill a set of (pt, eta, phi) histograms at once
    auto fill_all = [&](TH1D* h_pt, TH1D* h_eta, TH1D* h_phi) {
        h_pt->Fill(pt);
        h_eta->Fill(eta);
        h_phi->Fill(phi);
    };

    // --- CASE 1: The PLS is FAKE ---
    if (is_pls_fake) {
        if (is_pls_used) {
            fill_all(h_pt_f_u, h_eta_f_u, h_phi_f_u);
        } else {
            fill_all(h_pt_f_un, h_eta_f_un, h_phi_f_un);
        }
        return; // Processing for this PLS is done
    }

    // --- From here on, we know the PLS is REAL ---

    // Find the matching LS status for this PLS's sim_idx.
    // If the PLS has no valid sim_idx, the iterator will be .end()
    auto match_it = (pls_sim_idx == -999) 
                    ? simIdx_to_ls_status.end() 
                    : simIdx_to_ls_status.find(pls_sim_idx);

    // --- CASE 2: The REAL PLS has NO MATCHING LS object ---
    if (match_it == simIdx_to_ls_status.end()) {
        if (is_pls_used) {
            fill_all(h_pt_r_u_um, h_eta_r_u_um, h_phi_r_u_um);
        } else {
            fill_all(h_pt_r_un_um, h_eta_r_un_um, h_phi_r_un_um);
        }
    } 
    // --- CASE 3: The REAL PLS HAS a MATCHING LS object ---
    else {
        const auto& status = match_it->second;

        // Check the status of the matched LS(s)
        if (status.used_count > 0 && status.unused_count > 0) {
            // Matched to BOTH used and unused LS objects
            if (is_pls_used) fill_all(h_pt_r_u_b, h_eta_r_u_b, h_phi_r_u_b);
            else             fill_all(h_pt_r_un_b, h_eta_r_un_b, h_phi_r_un_b);
        }
        else if (status.used_count > 0) {
            // Matched to USED LS objects ONLY
            if (is_pls_used) fill_all(h_pt_r_u_uo, h_eta_r_u_uo, h_phi_r_u_uo);
            else             fill_all(h_pt_r_un_uo, h_eta_r_un_uo, h_phi_r_un_uo);
        }
        else { // This implies status.unused_count > 0
            // Matched to UNUSED LS objects ONLY
            if (is_pls_used) fill_all(h_pt_r_u_uno, h_eta_r_u_uno, h_phi_r_u_uno);
            else             fill_all(h_pt_r_un_uno, h_eta_r_un_uno, h_phi_r_un_uno);
        }
    }
}
// MODULE 5: Fill Simple Histograms
void fillSimpleHistograms(
    bool is_used, 
    bool is_fake, 
    float pt, 
    float eta, 
    float phi, 
    TH1D* h_pt_r_u,   TH1D* h_pt_f_u,   TH1D* h_pt_r_un,   TH1D* h_pt_f_un, 
    TH1D* h_eta_r_u,  TH1D* h_eta_f_u,  TH1D* h_eta_r_un,  TH1D* h_eta_f_un, 
    TH1D* h_phi_r_u,  TH1D* h_phi_f_u,  TH1D* h_phi_r_un,  TH1D* h_phi_f_un
) {
    // The logic is a simple 2x2 matrix of possibilities:
    // - Fake vs. Real
    // - Used vs. Unused

    if (is_fake) {
        // --- CASE 1: The object is FAKE ---
        if (is_used) {
            // Category: Fake and Used
            h_pt_f_u->Fill(pt);
            h_eta_f_u->Fill(eta);
            h_phi_f_u->Fill(phi);
        } else {
            // Category: Fake and Unused
            h_pt_f_un->Fill(pt);
            h_eta_f_un->Fill(eta);
            h_phi_f_un->Fill(phi);
        }
    } else {
        // --- CASE 2: The object is REAL ---
        if (is_used) {
            // Category: Real and Used
            h_pt_r_u->Fill(pt);
            h_eta_r_u->Fill(eta);
            h_phi_r_u->Fill(phi);
        } else {
            // Category: Real and Unused
            h_pt_r_un->Fill(pt);
            h_eta_r_un->Fill(eta);
            h_phi_r_un->Fill(phi);
        }
    }
}

// MODULE 6: Create and Save a Simple Stacked Plot
void createAndSaveSimpleStack(
    const std::string& canvasName, 
    const std::string& stackTitle, 
    const std::string& saveFileName, 
    const std::string& legendHeader, 
    TH1D* h_real, 
    TH1D* h_fake, 
    const std::string& real_label, 
    const std::string& fake_label, 
    Color_t real_color, 
    Color_t fake_color, 
    bool setLogY = false
) {
    // --- 1. Canvas and THStack Setup ---
    TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    if (setLogY) {
        canvas->SetLogy();
    }
    THStack *hs = new THStack("hs", stackTitle.c_str());

    // --- 2. Histogram Styling ---
    h_real->SetFillColor(real_color);
    h_fake->SetFillColor(fake_color);

    // --- 3. Stacking Histograms ---
    hs->Add(h_fake);
    hs->Add(h_real); 
    
    // --- 4. Draw the Stack ---
    hs->Draw("HIST"); 

    // --- 5. Create and Draw the Legend ---
    auto legend = new TLegend(0.4, 0.8, 0.9, 0.9); 
    legend->SetHeader(legendHeader.c_str(), "C"); 
    legend->AddEntry(h_real, real_label.c_str(), "f");
    legend->AddEntry(h_fake, fake_label.c_str(), "f");
    legend->SetBorderSize(0); 
    legend->Draw(); 
    
    // --- 6. Save and Clean up ---
    canvas->SaveAs(saveFileName.c_str()); 
    delete canvas;
}
