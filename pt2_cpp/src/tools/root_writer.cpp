#include "root_writer.h"
#include <filesystem>
#include <stdexcept>
#include <sstream>

namespace fs = std::filesystem;

Pt2Writer::Pt2Writer(const std::string& outDir_)
    : outDir(outDir_) 
{
    fs::create_directories(outDir);
}

void Pt2Writer::createNewFile(const std::string& filename) {
    if (file) {
        file->Write();
        file->Close();
        delete file;
        file = nullptr;
        tree = nullptr;
        entryCount = 0;
    }

    if (fs::exists(filename)) {
        throw std::runtime_error("File exists: " + filename);
    }

    file = new TFile(filename.c_str(), "RECREATE");
    tree = new TTree("tree", "pT2 dataset");

    // Branches
    tree->Branch("is_real", &b_isReal);

    tree->Branch("ls_pt", &b_ls_pt);
    tree->Branch("ls_eta", &b_ls_eta);
    tree->Branch("ls_phi", &b_ls_phi);

    tree->Branch("md_0_pt", &b_md_0_pt);
    tree->Branch("md_0_eta", &b_md_0_eta);
    tree->Branch("md_0_phi", &b_md_0_phi);
    tree->Branch("md_0_anchor_x", &b_md_0_anchor_x);
    tree->Branch("md_0_anchor_y", &b_md_0_anchor_y);
    tree->Branch("md_0_anchor_z", &b_md_0_anchor_z);
    tree->Branch("md_0_other_x", &b_md_0_other_x);
    tree->Branch("md_0_other_y", &b_md_0_other_y);
    tree->Branch("md_0_other_z", &b_md_0_other_z);
    tree->Branch("md_0_layer", &b_md_0_layer);
    tree->Branch("md_1_pt", &b_md_1_pt);
    tree->Branch("md_1_eta", &b_md_1_eta);
    tree->Branch("md_1_phi", &b_md_1_phi);
    tree->Branch("md_1_anchor_x", &b_md_1_anchor_x);
    tree->Branch("md_1_anchor_y", &b_md_1_anchor_y);
    tree->Branch("md_1_anchor_z", &b_md_1_anchor_z);
    tree->Branch("md_1_other_x", &b_md_1_other_x);
    tree->Branch("md_1_other_y", &b_md_1_other_y);
    tree->Branch("md_1_other_z", &b_md_1_other_z);
    tree->Branch("md_1_layer", &b_md_1_layer);

    tree->Branch("pls_pt", &b_pls_pt);
    tree->Branch("pls_phi", &b_pls_phi);
    tree->Branch("pls_eta", &b_pls_eta);
    tree->Branch("pls_hit0_x", &b_pls_hit0_x);
    tree->Branch("pls_hit0_y", &b_pls_hit0_y);
    tree->Branch("pls_hit0_z", &b_pls_hit0_z);
    tree->Branch("pls_hit1_x", &b_pls_hit1_x);
    tree->Branch("pls_hit1_y", &b_pls_hit1_y);
    tree->Branch("pls_hit1_z", &b_pls_hit1_z);
    tree->Branch("pls_hit2_x", &b_pls_hit2_x);
    tree->Branch("pls_hit2_y", &b_pls_hit2_y);
    tree->Branch("pls_hit2_z", &b_pls_hit2_z);
    tree->Branch("pls_hit3_x", &b_pls_hit3_x);
    tree->Branch("pls_hit3_y", &b_pls_hit3_y);
    tree->Branch("pls_hit3_z", &b_pls_hit3_z);
    tree->Branch("pls_nhit", &b_pls_nhit);
    tree->Branch("pls_charge", &b_pls_charge);

    tree->Branch("pt2_delta_pt", &b_pt2_delta_pt);
    tree->Branch("pt2_delta_eta", &b_pt2_delta_eta);
    tree->Branch("pt2_delta_phi", &b_pt2_delta_phi);
    tree->Branch("pt2_delta_R", &b_pt2_delta_R);
    tree->Branch("pt2_md0_dxy", &b_pt2_md0_dxy);
    tree->Branch("pt2_md0_dz", &b_pt2_md0_dz);
    tree->Branch("pt2_md1_dxy", &b_pt2_md1_dxy);
    tree->Branch("pt2_md1_dz", &b_pt2_md1_dz);
    tree->Branch("pt2_md0_rz", &b_pt2_md0_rz);
    tree->Branch("pt2_md1_rz", &b_pt2_md1_rz);

    file->SetCompressionLevel(1);

    entryCount = 0;
}

void Pt2Writer::fill(
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
) {
    b_isReal = pt2_isReal;

    b_ls_pt = ls_pt;
    b_ls_eta = ls_eta;
    b_ls_phi = ls_phi;

    b_md_0_pt = md_0_pt;
    b_md_0_eta = md_0_eta;
    b_md_0_phi = md_0_phi;
    b_md_0_anchor_x = md_0_anchor_x;
    b_md_0_anchor_y = md_0_anchor_y;
    b_md_0_anchor_z = md_0_anchor_z;
    b_md_0_other_x = md_0_other_x;
    b_md_0_other_y = md_0_other_y;
    b_md_0_other_z = md_0_other_z;
    b_md_0_layer = md_0_layer;
    b_md_1_pt = md_1_pt;
    b_md_1_eta = md_1_eta;
    b_md_1_phi = md_1_phi;
    b_md_1_anchor_x = md_1_anchor_x;
    b_md_1_anchor_y = md_1_anchor_y;
    b_md_1_anchor_z = md_1_anchor_z;
    b_md_1_other_x = md_1_other_x;
    b_md_1_other_y = md_1_other_y;
    b_md_1_other_z = md_1_other_z;
    b_md_1_layer = md_1_layer;

    b_pls_pt = pls_pt;
    b_pls_phi = pls_phi;
    b_pls_eta = pls_eta;
    b_pls_hit0_x = pls_hit0_x;
    b_pls_hit0_y = pls_hit0_y;
    b_pls_hit0_z = pls_hit0_z;
    b_pls_hit1_x = pls_hit1_x;
    b_pls_hit1_y = pls_hit1_y;
    b_pls_hit1_z = pls_hit1_z;
    b_pls_hit2_x = pls_hit2_x;
    b_pls_hit2_y = pls_hit2_y;
    b_pls_hit2_z = pls_hit2_z;
    b_pls_hit3_x = pls_hit3_x;
    b_pls_hit3_y = pls_hit3_y;
    b_pls_hit3_z = pls_hit3_z;
    b_pls_charge = pls_charge;
    b_pls_nhit = pls_nhit;


    b_pt2_delta_pt = pt2_delta_pt;
    b_pt2_delta_eta = pt2_delta_eta;
    b_pt2_delta_phi = pt2_delta_phi;
    b_pt2_delta_R = pt2_delta_R;
    b_pt2_md0_dxy = pt2_md0_dxy;
    b_pt2_md0_dz = pt2_md0_dz;
    b_pt2_md1_dxy = pt2_md1_dxy;
    b_pt2_md1_dz = pt2_md1_dz;
    b_pt2_md0_rz = pt2_md0_rz;
    b_pt2_md1_rz = pt2_md1_rz;

    tree->Fill();
    entryCount++;
}

void Pt2Writer::close() {
    if (file) {
        file->Write();
        file->Close();
        delete file;
        file = nullptr;
        tree = nullptr;
        entryCount =0;
    }
}

Pt2Writer::~Pt2Writer() {
    close();
}
