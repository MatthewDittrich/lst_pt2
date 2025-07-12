// 06/25/25

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TStyle.h>
#include <TMath.h> // Include for TMath::Pi()
#include <TLegend.h>

// The main function name must match the filename
void graph_falseID() {
    // --- Configuration ---
    const char* input_filename = "LSTNtuple.root";
    const char* tree_name = "tree"; // Make sure this is the correct TTree name too!

    // The names of the branches you want to plot (UPDATED WITH YOUR NAMES)
    const char* pt_branch_name = "ls_pt";
    const char* eta_branch_name = "ls_eta";
    const char* phi_branch_name = "ls_phi";
    const char* fake_branch_name = "ls_isFake";

    // --- Open the .root file ---
    TFile *inputFile = TFile::Open(input_filename, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cout << "Error: Could not open file " << input_filename << std::endl;
        return;
    }

    // --- Get the TTree from the file ---
    TTree *inputTree;
    inputFile->GetObject(tree_name, inputTree);
    if (!inputTree) {
        std::cout << "Error: Could not find TTree with name " << tree_name << std::endl;
        inputFile->Close();
        return;
    }
    std::cout << "Successfully accessed TTree: " << tree_name << std::endl;

    // --- Link each branch to its own C++ variable ---
    // --- Link each branch to its own C++ vector pointer ---
    std::vector<float> *pt_vec = nullptr;
    std::vector<float> *eta_vec = nullptr;
    std::vector<float> *phi_vec = nullptr;
    std::vector<int> *is_fake_vec = nullptr;

    inputTree->SetBranchAddress(pt_branch_name, &pt_vec);
    inputTree->SetBranchAddress(eta_branch_name, &eta_vec);
    inputTree->SetBranchAddress(phi_branch_name, &phi_vec);
    inputTree->SetBranchAddress(fake_branch_name, &is_fake_vec);
    
    //STACK histograms
    auto hs_pt = new THStack("hs_pt", "Stacked pT Distribution;p_{T} [GeV];Tracks");
    auto hs_eta = new THStack("hs_pt", "Stacked #eta Distribution;#eta ;Tracks");
    auto hs_phi = new THStack("hs_pt", "Stacked #phi Distribution;#phi ;Tracks");

    // --- Define three separate histograms ---
    // REAL
    TH1D *h_pt_real  = new TH1D("h_pt_real",  "p_{T}; p_{T} [GeV];Events", 180, 0.0, 250.0);
    TH1D *h_eta_real = new TH1D("h_eta_real", "#eta; #eta;Events", 180, -4.5, 4.5);
    TH1D *h_phi_real = new TH1D("h_phi_real", "#phi; #phi [rad];Events", 180, -3.5, 3.5);
    
    // FALSE
    TH1D *h_pt_fake  = new TH1D("h_pt_fake",  "p_{T}; p_{T} [GeV];Events", 180, 0.0, 250.0);
    TH1D *h_eta_fake = new TH1D("h_eta_fake", "#eta; #eta;Events", 180, -4.5, 4.5);
    TH1D *h_phi_fake = new TH1D("h_phi_fake", "#phi; #phi [rad];Events", 180, -3.5, 3.5);

    // --- The Event Loop: Loop over events, then over tracks in each event ---
    Long64_t nEntries = inputTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    for (Long64_t i = 0; i < nEntries; i++) {
        inputTree->GetEntry(i); // This fills the vectors (pt_vec, eta_vec, phi_vec) for this event.

        
       // Loop over all the tracks found in this single event
        for (size_t j = 0; j < pt_vec->size(); ++j) {
            // Check the fake flag FOR EACH TRACK
            if (is_fake_vec->at(j) == 0) { // 0 is REAL
                h_pt_real->Fill(pt_vec->at(j));
                h_eta_real->Fill(eta_vec->at(j));
                h_phi_real->Fill(phi_vec->at(j));
            } else { // 1 is FAKE
                h_pt_fake->Fill(pt_vec->at(j));
                h_eta_fake->Fill(eta_vec->at(j));
                h_phi_fake->Fill(phi_vec->at(j));
            }
        }
    }
    std::cout << "Finished processing entries." << std::endl;

    // --- Draw, Style, and Save the Plots ---
    gStyle->SetOptStat(0);
    // --- pT Plot ---
    TCanvas *c_pt = new TCanvas("c_pt", "pT Comparison", 800, 600);
    c_pt->SetLogy(); // Use log scale for pT
    c_pt->SetTopMargin(0.1);
    c_pt->SetBottomMargin(0.15);
    c_pt->SetLeftMargin(0.15);
    c_pt->SetRightMargin(0.05);
    h_pt_fake->SetLineColor(kRed - 2);
    h_pt_fake->SetFillColor(kRed - 6);
    h_pt_fake->SetLineWidth(2);
    h_pt_fake->GetXaxis()->SetLabelOffset(0.025); 
    h_pt_fake->GetXaxis()->SetTickLength(-0.03);
    h_pt_fake->GetXaxis()->SetTitleOffset(1.4);
    h_pt_fake->GetYaxis()->SetLabelOffset(0.025); 
    h_pt_fake->GetYaxis()->SetTickLength(-0.03);
    h_pt_fake->GetYaxis()->SetTitleOffset(1.6);
    //h_pt_fake->Draw("hist"); // Draw fake histogram first
    h_pt_real->SetFillColor(kBlue-7);
    h_pt_real->SetLineWidth(2);
    
    //stacking histograms
    hs_pt->Add(h_pt_real);
    hs_pt->Add(h_pt_fake);
    hs_pt->Draw("hist");
    //h_pt_real->Draw("hist same"); // Draw real on the same canvas
    auto legend_pt = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend_pt->AddEntry(h_pt_real, "Real Tracks", "l");
    legend_pt->AddEntry(h_pt_fake, "Fake Tracks", "l");
    legend_pt->Draw();
    c_pt->SaveAs("compare_pt.png");
    c_pt->SaveAs("compare_pt.pdf");

    // --- Eta Plot ---
    TCanvas *c_eta = new TCanvas("c_eta", "Eta Comparison", 800, 600);
    c_eta->SetTopMargin(0.1);
    c_eta->SetBottomMargin(0.15);
    c_eta->SetLeftMargin(0.15);
    c_eta->SetRightMargin(0.05);
    h_eta_fake->SetLineColor(kGreen - 2);
    h_eta_fake->SetFillColor(kGreen - 6);
    h_eta_fake->SetLineWidth(2);
    h_eta_fake->GetXaxis()->SetLabelOffset(0.025); 
    h_eta_fake->GetXaxis()->SetTickLength(-0.03);
    h_eta_fake->GetXaxis()->SetTitleOffset(1.4);
    h_eta_fake->GetYaxis()->SetLabelOffset(0.025); 
    h_eta_fake->GetYaxis()->SetTickLength(-0.03);
    h_eta_fake->GetYaxis()->SetTitleOffset(1.6);
    //h_eta_fake->Draw("hist");
    h_eta_real->SetLineWidth(2);
    h_eta_real->SetFillColor(kBlue - 7);
    //    h_eta_real->SetLineColor(kBlue);
    //stacking histograms
    hs_eta->Add(h_eta_real);
    hs_eta->Add(h_eta_fake);
    hs_eta->Draw("hist");
    //h_eta_real->Draw("hist same");
    auto legend_eta = new TLegend(0.2, 0.8, 0.4, 0.9);
    legend_eta->AddEntry(h_eta_real, "Real Tracks", "l");
    legend_eta->AddEntry(h_eta_fake, "Fake Tracks", "l");
    legend_eta->Draw();
    c_eta->SaveAs("compare_eta.png");
    c_eta->SaveAs("compare_eta.pdf");

    // --- Phi Plot ---
    TCanvas *c_phi = new TCanvas("c_phi", "Phi Comparison", 800, 600);
    c_phi->SetTopMargin(0.1);
    c_phi->SetBottomMargin(0.15);
    c_phi->SetLeftMargin(0.15);
    c_phi->SetRightMargin(0.05);
    h_phi_fake->SetLineColor(kMagenta - 2);
    h_phi_fake->SetFillColor(kMagenta - 6);
    h_phi_fake->SetLineWidth(2);
    h_phi_fake->GetXaxis()->SetLabelOffset(0.025); 
    h_phi_fake->GetXaxis()->SetTickLength(-0.03);
    h_phi_fake->GetXaxis()->SetTitleOffset(1.4);
    h_phi_fake->GetYaxis()->SetLabelOffset(0.025); 
    h_phi_fake->GetYaxis()->SetTickLength(-0.03);
    h_phi_fake->GetYaxis()->SetTitleOffset(1.6);
    //h_phi_fake->Draw("hist");
    h_phi_real->SetFillColor(kBlue - 7);
    h_phi_real->SetLineWidth(2);
    //stacking histograms
    hs_phi->Add(h_phi_real);
    hs_phi->Add(h_phi_fake);
    hs_phi->Draw("hist");
    //h_phi_real->Draw("hist same");
    auto legend_phi = new TLegend(0.6, 0.85, 0.8, 0.95);
    legend_phi->AddEntry(h_phi_real, "Real Tracks", "l");
    legend_phi->AddEntry(h_phi_fake, "Fake Tracks", "l");
    legend_phi->Draw();
    c_phi->SaveAs("compare_phi.png");
    c_phi->SaveAs("compare_phi.pdf");
    std::cout << "Plots saved as ls_pt_distribution.png, ls_eta_distribution.png, and ls_phi_distribution.png" << std::endl;

    // --- Clean up ---
    inputFile->Close();
}
