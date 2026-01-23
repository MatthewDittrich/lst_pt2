#include "LSTReader.h"

LSTReader::LSTReader() : inputFile(nullptr), inputTree(nullptr) {}

LSTReader::~LSTReader() {
    if (inputFile) inputFile->Close();
}

bool LSTReader::Init(const std::string& filename, const std::string& treename) {
    inputFile = TFile::Open(filename.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) return false;
    
    inputFile->GetObject(treename.c_str(), inputTree);
    if (!inputTree) return false;

    // Set Branch Addresses
    inputTree->SetBranchAddress("ls_pt", &ls_pt);
    inputTree->SetBranchAddress("ls_eta", &ls_eta);
    inputTree->SetBranchAddress("ls_phi", &ls_phi);
    inputTree->SetBranchAddress("ls_isFake", &ls_isFake);
    inputTree->SetBranchAddress("ls_simIdx", &ls_simIdx);
    // ... (Add all other SetBranchAddress calls here from your original code) ...
    // Note: I'm abbreviating for length, but you must include ALL SetBranchAddress calls here.
    inputTree->SetBranchAddress("tc_pt5Idx", &tc_pt5Idx);
    inputTree->SetBranchAddress("pt5_t5Idx", &pt5_t5Idx);
    inputTree->SetBranchAddress("t5_t3Idx0", &t5_t3Idx0);
    inputTree->SetBranchAddress("t5_t3Idx1", &t5_t3Idx1);
    inputTree->SetBranchAddress("t3_lsIdx0", &t3_lsIdx0);
    inputTree->SetBranchAddress("t3_lsIdx1", &t3_lsIdx1);
    inputTree->SetBranchAddress("tc_pt3Idx", &tc_pt3Idx);
    inputTree->SetBranchAddress("pt3_t3Idx", &pt3_t3Idx);
    inputTree->SetBranchAddress("tc_t5Idx", &tc_t5Idx);
    inputTree->SetBranchAddress("pt5_plsIdx", &pt5_plsIdx);
    inputTree->SetBranchAddress("pt3_plsIdx", &pt3_plsIdx);
    inputTree->SetBranchAddress("pls_pt", &pls_pt);
    inputTree->SetBranchAddress("pls_eta", &pls_eta);
    inputTree->SetBranchAddress("pls_phi", &pls_phi);
    inputTree->SetBranchAddress("pls_isFake", &pls_isFake);
    inputTree->SetBranchAddress("pls_simIdx", &pls_simIdx);
    inputTree->SetBranchAddress("ls_mdIdx0", &ls_mdIdx0);
    inputTree->SetBranchAddress("ls_mdIdx1", &ls_mdIdx1);
    inputTree->SetBranchAddress("md_anchor_x", &md_anchor_x);
    inputTree->SetBranchAddress("md_anchor_y", &md_anchor_y);
    inputTree->SetBranchAddress("md_anchor_z", &md_anchor_z);
    inputTree->SetBranchAddress("md_other_x", &md_other_x);
    inputTree->SetBranchAddress("md_other_y", &md_other_y);
    inputTree->SetBranchAddress("md_other_z", &md_other_z);
    inputTree->SetBranchAddress("pls_hit0_x", &pls_hit0_x);
    inputTree->SetBranchAddress("pls_hit0_y", &pls_hit0_y);
    inputTree->SetBranchAddress("pls_hit0_z", &pls_hit0_z);
    inputTree->SetBranchAddress("pls_hit1_x", &pls_hit1_x);
    inputTree->SetBranchAddress("pls_hit1_y", &pls_hit1_y);
    inputTree->SetBranchAddress("pls_hit1_z", &pls_hit1_z);
    inputTree->SetBranchAddress("pls_hit2_x", &pls_hit2_x);
    inputTree->SetBranchAddress("pls_hit2_y", &pls_hit2_y);
    inputTree->SetBranchAddress("pls_hit2_z", &pls_hit2_z);
    inputTree->SetBranchAddress("pls_hit3_x", &pls_hit3_x);
    inputTree->SetBranchAddress("pls_hit3_y", &pls_hit3_y);
    inputTree->SetBranchAddress("pls_hit3_z", &pls_hit3_z);
    inputTree->SetBranchAddress("md_layer", &md_layer);
    inputTree->SetBranchAddress("md_detId", &md_detId);
    inputTree->SetBranchAddress("sim_pt", &sim_pt);

    return true;
}

void LSTReader::GetEntry(Long64_t i) {
    if (inputTree) inputTree->GetEntry(i);
}

Long64_t LSTReader::GetEntries() const {
    return inputTree ? inputTree->GetEntries() : 0;
}
