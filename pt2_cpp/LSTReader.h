#ifndef LSTREADER_H
#define LSTREADER_H

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <iostream>

class LSTReader {
public:
    LSTReader();
    ~LSTReader();
    bool Init(const std::string& filename, const std::string& treename);
    void GetEntry(Long64_t i);
    Long64_t GetEntries() const;

    // Data pointers (Public for easy access by PhysicsAlgo)
    std::vector<float> *ls_pt = nullptr, *ls_eta = nullptr, *ls_phi = nullptr;
    std::vector<int>   *ls_isFake = nullptr, *ls_simIdx = nullptr;
    std::vector<int>   *tc_pt5Idx = nullptr, *pt5_t5Idx = nullptr, *t5_t3Idx0 = nullptr, *t5_t3Idx1 = nullptr, *t3_lsIdx0 = nullptr, *t3_lsIdx1 = nullptr;
    std::vector<int>   *tc_pt3Idx = nullptr, *pt3_t3Idx = nullptr, *tc_t5Idx = nullptr, *pt5_plsIdx = nullptr, *pt3_plsIdx = nullptr;
    std::vector<float> *pls_pt = nullptr, *pls_eta = nullptr, *pls_phi = nullptr;
    std::vector<int>   *pls_isFake = nullptr, *pls_simIdx = nullptr;
    std::vector<int>   *ls_mdIdx0 = nullptr, *ls_mdIdx1 = nullptr;
    std::vector<float> *md_anchor_x = nullptr, *md_anchor_y = nullptr, *md_anchor_z = nullptr;
    std::vector<float> *md_other_x = nullptr, *md_other_y = nullptr, *md_other_z = nullptr;
    std::vector<float> *pls_hit0_x = nullptr, *pls_hit0_y = nullptr, *pls_hit0_z = nullptr;
    std::vector<float> *pls_hit1_x = nullptr, *pls_hit1_y = nullptr, *pls_hit1_z = nullptr;
    std::vector<float> *pls_hit2_x = nullptr, *pls_hit2_y = nullptr, *pls_hit2_z = nullptr;
    std::vector<float> *pls_hit3_x = nullptr, *pls_hit3_y = nullptr, *pls_hit3_z = nullptr;
    std::vector<int>   *md_layer = nullptr, *md_detId = nullptr;
    std::vector<float> *sim_pt = nullptr;

private:
    TFile* inputFile;
    TTree* inputTree;
};
#endif
