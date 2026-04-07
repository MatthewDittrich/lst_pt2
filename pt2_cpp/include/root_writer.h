#pragma once

#include <TFile.h>
#include <TTree.h>
#include <string>

class Pt2Writer {
public:
    Pt2Writer(const std::string& outDir_);
    ~Pt2Writer();

    void fill(
        bool pt2_isReal,

        float ls_pt,
        float ls_eta,
        float ls_phi,

        float md_0_pt,
        float md_0_eta,
        float md_0_phi,
        float md_0_anchor_x,
        float md_0_anchor_y,
        float md_0_anchor_z,
        float md_0_other_x,
        float md_0_other_y,
        float md_0_other_z,
        int md_0_layer,

        float md_1_pt,
        float md_1_eta,
        float md_1_phi,
        float md_1_anchor_x,
        float md_1_anchor_y,
        float md_1_anchor_z,
        float md_1_other_x,
        float md_1_other_y,
        float md_1_other_z,
        int md_1_layer,

        float pls_pt,
        float pls_phi,
        float pls_eta,
        float pls_hit0_x,
        float pls_hit0_y,
        float pls_hit0_z,
        float pls_hit1_x,
        float pls_hit1_y,
        float pls_hit1_z,
        float pls_hit2_x,
        float pls_hit2_y,
        float pls_hit2_z,
        float pls_hit3_x,
        float pls_hit3_y,
        float pls_hit3_z,
        int pls_charge,
        int pls_nhit,

        float pt2_delta_pt, 
        float pt2_delta_eta, 
        float pt2_delta_phi, 
        float pt2_delta_R,
        float pt2_md0_dxy, 
        float pt2_md0_dz, 
        float pt2_md1_dxy, 
        float pt2_md1_dz, 
        float pt2_md0_rz, 
        float pt2_md1_rz
    );

    // Close current ROOT file
    void close();

    // Create a new ROOT file with given filename
    void createNewFile(const std::string& filename);

private:
    std::string outDir;

    TFile* file = nullptr;
    TTree* tree = nullptr;

    size_t entryCount = 0;

    // Branch variables
    bool b_isReal;

    float b_ls_pt, b_ls_eta, b_ls_phi;

    float b_md_0_pt, b_md_0_eta, b_md_0_phi;
    float b_md_0_anchor_x, b_md_0_anchor_y, b_md_0_anchor_z;
    float b_md_0_other_x, b_md_0_other_y, b_md_0_other_z;
    int   b_md_0_layer;

    float b_md_1_pt, b_md_1_eta, b_md_1_phi;
    float b_md_1_anchor_x, b_md_1_anchor_y, b_md_1_anchor_z;
    float b_md_1_other_x, b_md_1_other_y, b_md_1_other_z;
    int   b_md_1_layer;

    float b_pls_pt, b_pls_phi, b_pls_eta;
    float b_pls_hit0_x, b_pls_hit0_y, b_pls_hit0_z;
    float b_pls_hit1_x, b_pls_hit1_y, b_pls_hit1_z;
    float b_pls_hit2_x, b_pls_hit2_y, b_pls_hit2_z;
    float b_pls_hit3_x, b_pls_hit3_y, b_pls_hit3_z;
    int   b_pls_charge;
    int   b_pls_nhit;

    float b_pt2_delta_pt, b_pt2_delta_eta, b_pt2_delta_phi, b_pt2_delta_R;
    float b_pt2_md0_dxy, b_pt2_md0_dz, b_pt2_md1_dxy, b_pt2_md1_dz, b_pt2_md0_rz, b_pt2_md1_rz;
};
