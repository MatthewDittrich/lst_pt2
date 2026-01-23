#include "Plotter.h"

#include <iostream>

namespace Plotter {

    void createAndSaveSimpleStack(
        const std::string& canvasName, 
        const std::string& stackTitle, 
        const std::string& saveFileName, 
        const std::string& legendHeader, 
        TH1D* h_real, 
        TH1D* h_fake, 
        const std::string& real_label, 
        const std::string& fake_label, 
        Color_t real_color, 
        Color_t fake_color, 
        bool setLogY
    ) {
        // 1. Setup Canvas
        TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
        if (setLogY) {
            canvas->SetLogy();
        }

        // 2. Setup Stack
        THStack *hs = new THStack("hs", stackTitle.c_str());

        // 3. Styling
        h_real->SetFillColor(real_color);
        h_fake->SetFillColor(fake_color);
        h_real->SetLineColor(kBlack); // Optional: cleaner look
        h_fake->SetLineColor(kBlack);

        // 4. Add to Stack (Order determines visual stacking)
        // Original logic: Fake on bottom, Real on top
        hs->Add(h_fake);
        hs->Add(h_real);

        // 5. Draw
        hs->Draw("HIST");

        // 6. Legend
        // Coordinates: x1, y1, x2, y2
        auto legend = new TLegend(0.4, 0.8, 0.9, 0.9);
        legend->SetHeader(legendHeader.c_str(), "C"); // "C" centers the header
        legend->AddEntry(h_real, real_label.c_str(), "f");
        legend->AddEntry(h_fake, fake_label.c_str(), "f");
        legend->SetBorderSize(0);
        legend->Draw();

        // 7. Save and Clean
        canvas->SaveAs(saveFileName.c_str());
        
        // Cleanup: Deleting the canvas deletes primitives attached to it (like the stack drawn on it)
        delete canvas;
        // Note: We do NOT delete h_real or h_fake here, they are owned by HistManager.
    }

    void createAndSaveDetailedStack(
        const std::string& canvasName, 
        const std::string& stackTitle, 
        const std::string& saveFileName, 
        const std::string& legendHeader, 
        bool setLogY, 
        TH1D* h_fake, 
        TH1D* h_real_unmatched, 
        TH1D* h_real_matchedToUnused, 
        TH1D* h_real_matchedToBoth, 
        TH1D* h_real_matchedToUsed
    ) {
        // 1. Setup Canvas
        TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
        if (setLogY) {
            canvas->SetLogy();
        }

        // 2. Setup Stack
        THStack *hs = new THStack("hs", stackTitle.c_str());

        // 3. Styling (Colors from original macro)
        h_fake->SetFillColor(kGray);
        h_real_unmatched->SetFillColor(kRed - 7);
        h_real_matchedToUnused->SetFillColor(kAzure + 1);
        h_real_matchedToBoth->SetFillColor(kOrange + 1);
        h_real_matchedToUsed->SetFillColor(kGreen + 2);

        // 4. Add to Stack 
        // Order: MatchedToUsed (Bottom) -> ... -> Fake (Top)
        // Note: The original code added in this order, making MatchedToUsed the bottom-most layer visually?
        // Let's stick to the order in your provided macro:
        hs->Add(h_real_matchedToUsed);
        hs->Add(h_real_matchedToBoth);
        hs->Add(h_real_matchedToUnused);
        hs->Add(h_real_unmatched);
        hs->Add(h_fake);

        // 5. Draw
        hs->Draw("HIST");

        // 6. Legend
        auto legend = new TLegend(0.40, 0.58, 0.88, 0.88);
        legend->SetHeader(legendHeader.c_str(), "C");
        legend->AddEntry(h_real_matchedToUsed,   "Real, Matched to USED Only",    "f");
        legend->AddEntry(h_real_matchedToBoth,   "Real, Matched to Both",         "f");
        legend->AddEntry(h_real_matchedToUnused, "Real, Matched to UNUSED Only",  "f");
        legend->AddEntry(h_real_unmatched,       "Real, No Match",                "f");
        legend->AddEntry(h_fake,                 "Fake",                          "f");
        legend->SetBorderSize(0);
        legend->Draw();

        // 7. Save and Clean
        canvas->SaveAs(saveFileName.c_str());
        delete canvas;
    }

    void createAndSaveSimplePlot(
        const std::string& canvasName, 
        const std::string& saveFileName, 
        TH1D* histo, 
        const std::string& plotTitle, 
        const std::string& xLabel, 
        const std::string& yLabel
    ) {
        // 1. Setup Canvas
        TCanvas *canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 700);

        // 2. Styling
        histo->SetTitle(plotTitle.c_str());
        histo->GetXaxis()->SetTitle(xLabel.c_str());
        histo->GetYaxis()->SetTitle(yLabel.c_str());
        histo->SetLineColor(kAzure + 2);
        histo->SetLineWidth(2);

        // 3. Draw
        histo->Draw("HIST");

        // 4. Add Entries Box
        double integral = histo->GetEntries(); // Or Integral(), usually Entries is preferred for general plots

        TPaveText *pt = new TPaveText(0.65, 0.82, 0.88, 0.88, "NDC");
        pt->SetFillColor(0); // Transparent
        pt->SetBorderSize(0);
        pt->SetTextAlign(12); // Left-aligned, vertically centered
        pt->SetTextFont(42);
        pt->SetTextSize(0.04);

        TString text = TString::Format("Entries: %.0f", integral);
        pt->AddText(text);
        pt->Draw("SAME");

        // 5. Save and Clean
        canvas->SaveAs(saveFileName.c_str());
        delete canvas;
    }
}
