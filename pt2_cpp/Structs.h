#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>
#include <cmath>
#include <TVector2.h>

struct PlsMatchStatus { int used_count = 0; int unused_count = 0; };
struct LsMatchStatus  { int used_count = 0; int unused_count = 0; };

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

#endif
