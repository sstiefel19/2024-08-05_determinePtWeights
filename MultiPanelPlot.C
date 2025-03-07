#include "TCanvas.h"
#include "TH1F.h"
#include "TPad.h"

void MultiPanelPlot() {
    // Create a canvas
    TCanvas *c = new TCanvas("c", "Multi-panel Plot", 800, 800);

    // Divide the canvas into 3x3 sub-pads with zero spacing between them
    c->Divide(3, 3, 0.0, 0.0); // No space between sub-pads

    // Create histograms for demonstration
    TH1F *histograms[9];
    for (int i = 0; i < 9; i++) {
        histograms[i] = new TH1F(Form("h%d", i), Form("Histogram %d", i), 100, -5, 5);
        histograms[i]->FillRandom("gaus", 1000); // Fill with random Gaussian data
    }

    // Loop over pads to customize margins and draw histograms
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            int padIndex = row * 3 + col + 1;
            c->cd(padIndex);

            TPad *pad = (TPad *)gPad;
            pad->SetLeftMargin(0.0);   // No left margin
            pad->SetBottomMargin(0.0); // No bottom margin
            pad->SetRightMargin(0.0);  // No right margin
            pad->SetTopMargin(0.0);    // No top margin

            if (col == 0) pad->SetLeftMargin(0.15);   // Add left margin for the first column
            if (row == 2) pad->SetBottomMargin(0.15); // Add bottom margin for the last row

            // Ensure ticks are drawn on all four sides of the subpad
            pad->SetTicks(1, 1); // Enable ticks on both x and y axes (all sides)

            histograms[padIndex - 1]->GetXaxis()->SetTickLength(0.03);
            histograms[padIndex - 1]->GetYaxis()->SetTickLength(0.03);

            // // Remove axis labels from inner subplots
            // if (row != 2) { // Remove x-axis labels unless it's the last row
            //     histograms[padIndex - 1]->GetXaxis()->SetLabelSize(0);
            // }
            // if (col != 0) { // Remove y-axis labels unless it's the first column
            //     histograms[padIndex - 1]->GetYaxis()->SetLabelSize(0);
            // }

            histograms[padIndex - 1]->Draw();
        }
    }

    c->Update();
}
