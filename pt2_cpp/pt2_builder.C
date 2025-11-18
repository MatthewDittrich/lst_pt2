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
const float ptBins[2] = {0.8, 2.0}; // only first bin threshold for example



// --- Calculate the superbin ---
// TODO Is this right?
int calculateSuperbin(float pt, float z, float eta, float phi) {

    // --- Determine pt bin ---
    int ipt = (pt < 2.0) ? 0 : 1;

    // --- Determine eta bin ---
    int ieta;
    if(eta < etaMin) {
//        std::cerr << "Warning: eta (" << eta << ") below etaMin (" << etaMin << "), using closest bin.\n";
        ieta = 0;
    } 
    else if(eta > etaMax) {
//        std::cerr << "Warning: eta (" << eta << ") above etaMax (" << etaMax << "), using closest bin.\n";
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
//        std::cerr << "Warning: z (" << z << ") below zMin (" << zMin << "), using closest bin.\n";
        iz = 0;
    } else if(z > zMax) {
//        std::cerr << "Warning: z (" << z << ") above zMax (" << zMax << "), using closest bin.\n";
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



// --- Calculate the correct Z value for the pLS ---
// Right now we are using the outermost hit z value. THIS MIGHT BE WRONG
std::vector<float> maxHitZ(const std::vector<float>& hit3, const std::vector<float>& hit2)
{
    if (hit3.size() != hit2.size()) {
        throw std::runtime_error("Error: pLS hit3 and hit2 vectors have different sizes!");
    }
    std::vector<float> combined(hit3.size(), -999.0f);
    for (size_t i = 0; i < hit3.size(); ++i) {
        if (hit2[i] == -999) {
            throw std::runtime_error("Error: pLS hit2 is missing!");
        }
        combined[i] = (hit3[i] == -999) ? hit2[i] : hit3[i];
    }
    return combined;
}



//// --- Read the pixelmaps ---
//std::vector<std::vector<int>> getPLSDetIDs(
//    const std::vector<float>& pls_pt_vec,
//    const std::vector<float>& pls_max_z_vec,
//    const std::vector<float>& pls_eta_vec,
//    const std::vector<float>& pls_phi_vec
//    ) {
//
//    // The outer length of this will be the number of pls
//    size_t N_pls = pls_pt_vec.size();
//    std::vector<std::vector<int>> pls_detID_vec(N_pls);
//
//    // --- Load pixel map file ---
//    std::map<int, std::vector<int>> superbin_map;
//    std::ifstream infile("/blue/avery/matthew.dittrich/PT2_Studies/Geometry/Pixel_Maps/pLS_map_ElCheapo.txt");
//    if (!infile.is_open()) {
//        throw std::runtime_error("Error: Pixel Map could not be loaded!");
//    }
//    std::string line;
//    while (std::getline(infile, line)) {
//        std::istringstream iss(line);
//        int superbin, n_detids;
//        iss >> superbin >> n_detids;
//        std::vector<int> detIDs;
//        int det;
//        for (int i = 0; i < n_detids; ++i) {
//            iss >> det;
//            detIDs.push_back(det);
//        }
//        superbin_map[superbin] = detIDs;
//    }
//    infile.close();
//
//    // Loop over pLS and fill detID vectors
//    for (size_t i = 0; i < N_pls; ++i) {
//        // Calculate superbin for this pLS
//        int superbin = calculateSuperbin(pls_pt_vec[i], pls_max_z_vec[i], pls_eta_vec[i], pls_phi_vec[i]);
//        // Lookup the detIDs from the map
//        if (superbin_map.find(superbin) != superbin_map.end()) {
//            pls_detID_vec[i] = superbin_map[superbin];
//        } 
//        else {
//            throw std::runtime_error("Error: superbin not found for a given pLS!");
//        }
//    }
//    return pls_detID_vec;
//}



//// --- Basic pT2 building ---
//void connectPLSandLS(
//    const std::vector<float>& pls_pt_vec,
//    const std::vector<float>& pls_max_z_vec,
//    const std::vector<float>& pls_eta_vec,
//    const std::vector<float>& pls_phi_vec,
//    const std::vector<int>& ls_mdIdx0_vec,
//    const std::vector<int>& ls_mdIdx1_vec,
//    const std::vector<int>& md_detId_vec,
//    std::vector<int>& pls_pt2_idx_vec,
//    std::vector<int>& ls_pt2_idx_vec
//    ) {
//
//    // --- Sanity check for vector sizes ---
//    size_t N_pls = pls_pt_vec.size();
//    if(pls_max_z_vec.size() != N_pls || pls_eta_vec.size() != N_pls || pls_phi_vec.size() != N_pls){
//        throw std::runtime_error("Error: pLS kinematic vectors are not same size!");
//    }
//
//    // --- Create nested vector where the outer length is the N_pls and the inner is the relevant detIds for the respective pLS
//    std::vector<std::vector<int>> pls_detID_vec(N_pls);
//    for(size_t i = 0; i < N_pls; ++i){
//        pls_detID_vec[i] = getPLSDetIDs(pls_pt_vec[i], pls_max_z_vec[i], pls_eta_vec[i], pls_phi_vec[i]);
//    }
//
//    // --- Build the nested ls md vector ---
//    // TODO We might be able to avoid this by just using one of the md (upper vs lower)
//    size_t N_ls = ls_mdIdx0_vec.size();
//    if(ls_mdIdx1_vec.size() != N_ls){
//        throw std::runtime_error("Error: LS md vectors are not same size!");
//    }
//    std::vector<std::vector<int>> ls_mddetID_vec(N_ls, std::vector<int>(2));
//    for(size_t i = 0; i < N_ls; ++i){
//        ls_mddetID_vec[i][0] = md_detId_vec[ls_mdIdx0_vec[i]];
//        ls_mddetID_vec[i][1] = md_detId_vec[ls_mdIdx1_vec[i]];
//    }
//
//    // --- Fill the output lists ---
//    pls_pt2_idx_vec.clear();
//    ls_pt2_idx_vec.clear();
//
//    for(size_t pls_i = 0; pls_i < N_pls; ++pls_i){
//        for(size_t ls_i = 0; ls_i < N_ls; ++ls_i){
//            for(int detID : pls_detID_vec[pls_i]){
//                if(detID == ls_mddetID_vec[ls_i][0] || detID == ls_mddetID_vec[ls_i][1]){
//                    pls_pt2_idx_vec.push_back(pls_i);
//                    ls_pt2_idx_vec.push_back(ls_i);
//                    break;
//                }
//            }
//        }
//    }
//}



// --- Read the pixelmaps once ---
std::map<int, std::vector<int>> loadPixelMap(const std::string& filename) {
    std::map<int, std::vector<int>> superbin_map;
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

std::vector<std::vector<int>> getPLSDetIDs(
    const std::vector<float>& pls_pt_vec,
    const std::vector<float>& pls_max_z_vec,
    const std::vector<float>& pls_eta_vec,
    const std::vector<float>& pls_phi_vec,
    const std::map<int, std::vector<int>>& superbin_map
) {
    size_t N_pls = pls_pt_vec.size();
    std::vector<std::vector<int>> pls_detID_vec(N_pls);
    for (size_t i = 0; i < N_pls; ++i) {
        int superbin = calculateSuperbin(pls_pt_vec[i], pls_max_z_vec[i], pls_eta_vec[i], pls_phi_vec[i]);
        if (superbin_map.find(superbin) != superbin_map.end()) {
            pls_detID_vec[i] = superbin_map.at(superbin);
        } else {
            throw std::runtime_error("Error: superbin not found for a given pLS!");
        }
    }
    return pls_detID_vec;
}

// --- Basic pT2 building ---
void connectPLSandLS(
    const std::vector<float>& pls_pt_vec,
    const std::vector<float>& pls_max_z_vec,
    const std::vector<float>& pls_eta_vec,
    const std::vector<float>& pls_phi_vec,
    const std::vector<int>& ls_mdIdx0_vec,
    const std::vector<int>& ls_mdIdx1_vec,
    const std::vector<int>& md_detId_vec,
    std::vector<int>& pls_pt2_idx_vec,
    std::vector<int>& ls_pt2_idx_vec
) {
    // --- Sanity check for vector sizes ---
    size_t N_pls = pls_pt_vec.size();
    if(pls_max_z_vec.size() != N_pls || pls_eta_vec.size() != N_pls || pls_phi_vec.size() != N_pls){
        throw std::runtime_error("Error: pLS kinematic vectors are not same size!");
    }

    // --- Load pixel map once ---
    static std::map<int, std::vector<int>> superbin_map = loadPixelMap("/blue/avery/matthew.dittrich/PT2_Studies/Geometry/Pixel_Maps/pLS_map_ElCheapo.txt");

    // --- Create nested vector where the outer length is the N_pls and the inner is the relevant detIds for the respective pLS
    std::vector<std::vector<int>> pls_detID_vec = getPLSDetIDs(pls_pt_vec, pls_max_z_vec, pls_eta_vec, pls_phi_vec, superbin_map);

    // --- Build the nested ls md vector ---
    size_t N_ls = ls_mdIdx0_vec.size();
    if(ls_mdIdx1_vec.size() != N_ls){
        throw std::runtime_error("Error: LS md vectors are not same size!");
    }
    std::vector<std::vector<int>> ls_mddetID_vec(N_ls, std::vector<int>(2));
    for(size_t i = 0; i < N_ls; ++i){
        ls_mddetID_vec[i][0] = md_detId_vec[ls_mdIdx0_vec[i]];
        ls_mddetID_vec[i][1] = md_detId_vec[ls_mdIdx1_vec[i]];
    }

    // --- Fill the output lists ---
    pls_pt2_idx_vec.clear();
    ls_pt2_idx_vec.clear();

    for(size_t pls_i = 0; pls_i < N_pls; ++pls_i){
        for(size_t ls_i = 0; ls_i < N_ls; ++ls_i){
            for(int detID : pls_detID_vec[pls_i]){
                if(detID == ls_mddetID_vec[ls_i][0] || detID == ls_mddetID_vec[ls_i][1]){
                    pls_pt2_idx_vec.push_back(pls_i);
                    ls_pt2_idx_vec.push_back(ls_i);
                    break;
                }
            }
        }
    }
}
