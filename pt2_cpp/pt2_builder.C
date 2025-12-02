#include <vector>
#include <iostream>
#include <stdexcept>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>


// Hardcoded constants
const int nEta = 25;
const int nPhi = 72;
const int nZ   = 25;
const float etaMin = -2.6;
const float etaMax = 2.6;
const float zMin = -30.0;
const float zMax = 30.0;
const float ptBins[2] = {0.8, 2.0};



// Calculate the superbin for a given pls
int getSuperbin(float pt, float z, float eta, float phi) {

    // --- Determine pt bin ---
    int ipt = (pt < 2.0) ? 0 : 1;

    // --- Determine eta bin ---
    int ieta;
    if(eta < etaMin) {
        //std::cerr << "Warning: eta (" << eta << ") below etaMin (" << etaMin << "), using closest bin.\n";
        ieta = 0;
    } 
    else if(eta > etaMax) {
        //std::cerr << "Warning: eta (" << eta << ") above etaMax (" << etaMax << "), using closest bin.\n";
        ieta = nEta - 1;
    } 
    else {
        ieta = static_cast<int>((eta - etaMin) / ((etaMax - etaMin)/nEta));
        if(ieta < 0) ieta = 0;
        if(ieta >= nEta) ieta = nEta - 1;
    }

    // --- Determine phi bin ---
    while(phi < -M_PI) phi += 2*M_PI;
    while(phi >  M_PI) phi -= 2*M_PI;
    int iphi = static_cast<int>((phi + M_PI) / (2*M_PI / nPhi));
    if(iphi < 0) iphi = 0;
    if(iphi >= nPhi) iphi = nPhi - 1;

    // --- Determine z bin ---
    int iz;
    if(z < zMin) {
        //std::cerr << "Warning: z (" << z << ") below zMin (" << zMin << "), using closest bin.\n";
        iz = 0;
    } else if(z > zMax) {
        //std::cerr << "Warning: z (" << z << ") above zMax (" << zMax << "), using closest bin.\n";
        iz = nZ - 1;
    } else {
        iz = static_cast<int>((z - zMin) / ((zMax - zMin)/nZ));
        if(iz < 0) iz = 0;
        if(iz >= nZ) iz = nZ - 1;
    }

    // --- Compute the superbin index ---
    int superbin = (ipt * nPhi * nEta * nZ) + (ieta * nPhi * nZ) + (iphi * nZ) + iz;

    return superbin;
}



//// --- Calculate the correct Z value for the pLS ---
//std::vector<float> maxHitZ(const std::vector<float>& hit3, const std::vector<float>& hit2)
//{
//    if (hit3.size() != hit2.size()) {
//        throw std::runtime_error("Error: pLS hit3 and hit2 vectors have different sizes!");
//    }
//    std::vector<float> combined(hit3.size(), -999.0f);
//    for (size_t i = 0; i < hit3.size(); ++i) {
//        if (hit2[i] == -999) {
//            throw std::runtime_error("Error: pLS hit2 is missing!");
//        }
//        combined[i] = (hit3[i] == -999) ? hit2[i] : hit3[i];
//    }
//    return combined;
//}


// Make map: superbin --> detIDs
std::unordered_map<int, std::vector<int>> loadPixelMap(const std::string& filename) {
    std::unordered_map<int, std::vector<int>> superbin_map;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Error: Pixel Map could not be loaded!");
    }
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int superbin, n_detids;
        iss >> superbin >> n_detids;
        std::vector<int> detIDs;
        int det;
        for (int i = 0; i < n_detids; ++i) {
            iss >> det;
            detIDs.push_back(det);
        }
        superbin_map[superbin] = detIDs;
    }
    infile.close();
    return superbin_map;
}



// Make map: detID --> lsIdxs
std::unordered_map<int, std::vector<int>> makeDetIDLSMap(
    const std::vector<int>& ls_mdIdx0_vec,
    const std::vector<int>& ls_mdIdx1_vec,
    const std::vector<int>& md_detId_vec)
{
    std::unordered_map<int, std::vector<int>> detID_ls_map;
    const size_t nLS = ls_mdIdx0_vec.size();
    for (size_t i = 0; i < nLS; ++i) {
        int md0 = ls_mdIdx0_vec[i];
        int md1 = ls_mdIdx1_vec[i];
        // Check index bounds
        if (md0 < 0 || md0 >= (int)md_detId_vec.size() ||
            md1 < 0 || md1 >= (int)md_detId_vec.size()) 
        {
            throw std::runtime_error("md index out of bounds in makeDetIDLSMap()");
        }
        // Get the detIDs
        int det0 = md_detId_vec[md0];
        int det1 = md_detId_vec[md1];
        // Append LS index to both detIDs
        detID_ls_map[det0].push_back(i);
        detID_ls_map[det1].push_back(i);
    }
    return detID_ls_map;
}



//std::vector<std::vector<int>> getPLSDetIDs(
//    const std::vector<float>& pls_pt_vec,
//    const std::vector<float>& pls_max_z_vec,
//    const std::vector<float>& pls_eta_vec,
//    const std::vector<float>& pls_phi_vec,
//    const std::map<int, std::vector<int>>& superbin_map
//) {
//    size_t N_pls = pls_pt_vec.size();
//    std::vector<std::vector<int>> pls_detID_vec(N_pls);
//    for (size_t i = 0; i < N_pls; ++i) {
//        int superbin = calculateSuperbin(pls_pt_vec[i], pls_max_z_vec[i], pls_eta_vec[i], pls_phi_vec[i]);
//        if (superbin_map.find(superbin) != superbin_map.end()) {
//            pls_detID_vec[i] = superbin_map.at(superbin);
//        } else {
//            throw std::runtime_error("Error: superbin not found for a given pLS!");
//        }
//    }
//    return pls_detID_vec;
//}



// --- Main pT2 Builder ---
void connectPLSandLS(
    const std::vector<float>& pls_pt_vec,
    const std::vector<float>& pls_hit3_z_vec,
    const std::vector<float>& pls_hit2_z_vec,
    const std::vector<float>& pls_eta_vec,
    const std::vector<float>& pls_phi_vec,
    const std::vector<int>& ls_mdIdx0_vec,
    const std::vector<int>& ls_mdIdx1_vec,
    const std::vector<int>& md_detId_vec,
    std::vector<int>& pls_pt2_cand_idx_vec,
    std::vector<int>& ls_pt2_cand_idx_vec
) {

    size_t N_pls = pls_pt_vec.size();
    size_t N_ls = ls_mdIdx0_vec.size();

    // Check pls vectors
    if(pls_eta_vec.size() != N_pls || pls_phi_vec.size() != N_pls || pls_hit3_z_vec.size() != N_pls || pls_hit2_z_vec.size() != N_pls){
        throw std::runtime_error("Error: pLS kinematic vectors are not same size!");
    }
    // Check ls vectors
    if(ls_mdIdx1_vec.size() != N_ls){
        throw std::runtime_error("Error: LS kinematic vectors are not same size!");
    }

    // Load pixel map once per event (TODO: Load Before Event Loop?)
    // pls superbin --> vector of detIds
    static std::unordered_map<int, std::vector<int>> superbin_map = loadPixelMap("/blue/avery/matthew.dittrich/PT2_Studies/Geometry/Pixel_Maps/pLS_map_ElCheapo.txt");

    // Create detID --> vector of lsIdx
    static std::unordered_map<int, std::vector<int>> detID_ls_map = makeDetIDLSMap(ls_mdIdx0_vec, ls_mdIdx1_vec, md_detId_vec);

    pls_pt2_cand_idx_vec.clear();
    ls_pt2_cand_idx_vec.clear();

    // Main loop over pls
    for(size_t pls_i = 0; pls_i < N_pls; ++pls_i){

        // Right now we are using the max Z value. This is wrong, but we can adjust later.
        float pls_z = (pls_hit3_z_vec[pls_i] != -999.0f) ? pls_hit3_z_vec[pls_i] : pls_hit2_z_vec[pls_i]; 

        // Compute the superbin for this pls
        int superbin = getSuperbin(pls_pt_vec[pls_i], pls_z, pls_eta_vec[pls_i], pls_phi_vec[pls_i]);

        // Get the vector of detIDs corresponding to this superbin
        const auto& pls_detID_vec = superbin_map[superbin];

        // Record the ls we have added for this particular pls
        std::unordered_set<int> added_ls;

        for (int detID : pls_detID_vec) {

            // See if this detID is used in an LS
            auto it = detID_ls_map.find(detID);
            if (it == detID_ls_map.end()) continue;

            const auto& ls_indices = it->second;

            for (int ls_idx : ls_indices) {
                if (added_ls.insert(ls_idx).second) {
                    pls_pt2_cand_idx_vec.push_back(pls_i);
                    ls_pt2_cand_idx_vec.push_back(ls_idx);
                }
            }
        }
    }
}



// Fill the basic delta histograms
void fillCandHistos(
    const std::vector<int>& pls_idx_vec,
    const std::vector<int>& ls_idx_vec,
    const std::vector<float>& pls_pt_vec,
    const std::vector<float>& ls_pt_vec,
    const std::vector<float>& pls_eta_vec,
    const std::vector<float>& ls_eta_vec,
    const std::vector<float>& pls_phi_vec,
    const std::vector<float>& ls_phi_vec,
    TH1D* h_delta_pt,
    TH1D* h_delta_eta,
    TH1D* h_delta_phi,
    TH1D* h_delta_R
) {
    for (size_t i = 0; i < pls_idx_vec.size(); ++i) {
        int pls_i = pls_idx_vec[i];
        int ls_i  = ls_idx_vec[i];

        // Compute differences
        float delta_pt  = ls_pt_vec[ls_i]  - pls_pt_vec[pls_i];
        float delta_eta = ls_eta_vec[ls_i] - pls_eta_vec[pls_i];
        float delta_phi = TVector2::Phi_mpi_pi(pls_phi_vec[pls_i] - ls_phi_vec[ls_i]); // same sign as your snippet
        float delta_R = std::sqrt(delta_eta*delta_eta + delta_phi*delta_phi);

        // Fill histograms
        h_delta_pt->Fill(delta_pt);
        h_delta_eta->Fill(delta_eta);
        h_delta_phi->Fill(delta_phi);
        h_delta_R->Fill(delta_R);
    }
}

