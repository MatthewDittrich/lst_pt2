// 06/25/25

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>
#include <TMath.h> // Include for TMath::Pi()
#include <iomanip>


void graph2() {
    // root file name
    const char* input_filename = "LSTNtuple.root";
    const char* tree_name = "tree"; 

    // branches names
    const char* pt_branch_name = "ls_pt";
    const char* eta_branch_name = "ls_eta";
    const char* phi_branch_name = "ls_phi";

    // --- Open the root file ---
    TFile *inputFile = TFile::Open(input_filename, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cout << "Error: Could not open file " << input_filename << std::endl;
        return;
    }

    // --- Get TTree  ---
    TTree *inputTree;
    inputFile->GetObject(tree_name, inputTree);
    if (!inputTree) {
        std::cout << "Error: Could not find TTree with name " << tree_name << std::endl;
        inputFile->Close();
        return;
    }
    std::cout << "Successfully accessed TTree: " << tree_name << std::endl;

    // --- Link each branch to its own C++ vector pointer ---
    std::vector<float> *pt_vec = nullptr;
    std::vector<float> *eta_vec = nullptr;
    std::vector<float> *phi_vec = nullptr;

    inputTree->SetBranchAddress(pt_branch_name, &pt_vec);
    inputTree->SetBranchAddress(eta_branch_name, &eta_vec);
    inputTree->SetBranchAddress(phi_branch_name, &phi_vec);

    // --- Define the histograms ---
    TH1D *h_pt  = new TH1D("h_pt",  "p_{T}; p_{T} [GeV];Events", 180, 0.0, 250.0);
    TH1D *h_eta = new TH1D("h_eta", "#eta; #eta;Events", 180, -4.5, 4.5);
    TH1D *h_phi = new TH1D("h_phi", "#phi; #phi [rad];Events", 180, -3.15, 3.15);

    // --- The Event Loop: Loop over events, then over tracks in each event ---
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i); // load branch data
	
        // loop over all the tracks 
        for (size_t j = 0; j < pt_vec->size(); ++j) {
            //float current_phi = phi_vec->at(j);
            h_pt->Fill(pt_vec->at(j));
            //std::cout << "phi " << current_phi << std::endl;
            h_eta->Fill(eta_vec->at(j));
            h_phi->Fill(phi_vec->at(j));
        }
    }
    std::cout << "Finished processing entries." << std::endl;
// --- Print the final yields ---
    std::cout << "========================================" << std::endl;
    std::cout << "           Final Yields" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "Tracks:" << std::endl;
    std::cout << "  pT Histogram Entries:  " << h_pt->GetEntries() << std::endl;
    std::cout << "  Eta Histogram Entries: " << h_eta->GetEntries() << std::endl;
    std::cout << "  Phi Histogram Entries: " << h_phi->GetEntries() << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    
    // --- Draw, Style, and Save the Plots ---
    
    // Plot 1: Pt
    TCanvas *c_pt = new TCanvas("c_pt", "Pt Canvas", 800, 600);
    gPad->SetLogy();
    h_pt->SetLineColor(kRed - 2);
    h_pt->SetFillColor(kRed - 7);
    h_pt->SetLineWidth(2);
    h_pt->GetXaxis()->SetLabelOffset(0.025); 
    h_pt->GetXaxis()->SetTickLength(-0.03);
    h_pt->GetXaxis()->SetTitleOffset(1.4);
    h_pt->GetYaxis()->SetLabelOffset(0.025); 
    h_pt->GetYaxis()->SetTickLength(-0.03);
    h_pt->GetYaxis()->SetTitleOffset(1.6);
    //h_pt->SetMinimum(0);
    h_pt->Draw();
    c_pt->SaveAs("lst_pt_distribution.png");
    c_pt->SaveAs("lst_pt_distribution.pdf");

    // Plot 2: Eta
    TCanvas *c_eta = new TCanvas("c_eta", "Eta Canvas", 800, 600);
    h_eta->SetFillColor(kGreen - 7);
    h_eta->SetLineColor(kGreen - 2);
    h_eta->SetLineWidth(2);
    h_eta->GetXaxis()->SetLabelOffset(0.025); 
    h_eta->GetXaxis()->SetTickLength(-0.03);
    h_eta->GetXaxis()->SetTitleOffset(1.4);
    h_eta->GetYaxis()->SetLabelOffset(0.025); 
    h_eta->GetYaxis()->SetTickLength(-0.03);
    h_eta->GetYaxis()->SetTitleOffset(1.6);
    //h_eta->SetMinimum(0);
    h_eta->Draw();
    c_eta->SaveAs("lst_eta_distribution.png");
    c_eta->SaveAs("lst_eta_distribution.pdf");
    // Plot 3: Phi
    TCanvas *c_phi = new TCanvas("c_phi", "Phi Canvas", 800, 600);
    h_phi->SetFillColor(kBlue - 7);
    h_phi->SetLineColor(kBlue - 2);
    h_phi->SetLineWidth(2);
    h_phi->GetXaxis()->SetLabelOffset(0.025); 
    h_phi->GetXaxis()->SetTickLength(-0.03);
    h_phi->GetXaxis()->SetTitleOffset(1.4);
    h_phi->GetYaxis()->SetLabelOffset(0.025); 
    h_phi->GetYaxis()->SetTickLength(-0.03);
    h_phi->GetYaxis()->SetTitleOffset(1.6);
    h_phi->SetMinimum(0);
    h_phi->Draw();
    c_phi->SaveAs("lst_phi_distribution.png");
    c_phi->SaveAs("lst_phi_distribution.pdf");    
    std::cout << "Plots saved as ls_pt_distribution.png, ls_eta_distribution.png, and ls_phi_distribution.png" << std::endl;

    inputFile->Close();
}
