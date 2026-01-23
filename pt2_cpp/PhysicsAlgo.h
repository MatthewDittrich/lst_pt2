#ifndef PHYSICSALGO_H
#define PHYSICSALGO_H

#include <set>
#include <map>
#include <vector>
#include <utility>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include "Structs.h"
#include "LSTReader.h"
#include "HistManager.h"

namespace PhysicsAlgo {
    
    void findUsedIndices(std::vector<bool>& used_ls_indices, std::vector<bool>& used_pls_indices, const LSTReader& data);
    
    //void getPlsMatchStatus(std::map<int, PlsMatchStatus>& simIdx_to_pls_status, const std::vector<int>& pls_isFake, const std::vector<int>& pls_simIdx, const std::vector<bool>& used_pls_indices);
    
    //void getLsMatchStatus(std::map<int, LsMatchStatus>& simIdx_to_ls_status, const std::vector<int>& ls_isFake, const std::vector<int>& ls_simIdx, const std::vector<bool>& used_ls_indices);

    // Main extrapolation processing logic
    //void processIdealPT2s(const LSTReader& data, const std::vector<bool>& used_pls_indices, const std::vector<bool>& used_ls_indices, const std::map<int, std::vector<int>>& ls_simIdx_map, HistManager& hists, long long& ideal_pt2_count);
    void runPT2Physics(int pls_idx, 
                       const LSTReader& data, 
                       const std::vector<bool>& ls_is_used, 
                       const std::map<int, std::vector<int>>& ls_simIdx_map, 
                       HistManager& hists, 
                       long long& ideal_pt2_count);
    // Math Helpers
    std::pair<double, double> extrapolatePlsHelicallyAndGetDistance(int pls_idx, int ls_idx, const LSTReader& data);
    std::vector<double> extrapolateLsInRZAndGetDeltaR_AllHits(int pls_idx, int ls_idx, const LSTReader& data);
    std::pair<double, double> extrapolatePlsInRZAndGetDeltaZ(int pls_idx, int ls_idx, const LSTReader& data);
    
    // Internal helper
    ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double fixed_radius);
}

#endif
