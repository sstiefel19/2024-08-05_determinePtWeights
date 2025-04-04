/**********************************************************************************
* Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                *   
* All rights reserved.                                                            *                 
*                                                                                 *            
* Author: stephan.friedrich.stiefelmaier@cern.ch                                  *           
* Version: 1.0                                                                    *                    
*                                                                                 *                
* Redistribution and use in source and binary forms, with or without              *
* modification, are permitted provided that the following conditions are met:     *
*     * Redistributions of source code must retain the above copyright            *
*       notice, this list of conditions and the following disclaimer.             *
*     * Redistributions in binary form must reproduce the above copyright         *
*       notice, this list of conditions and the following disclaimer in the       *  
*       documentation and/or other materials provided with the distribution.      *
*     * Neither the name of the <organization> nor the                            *
*       names of its contributors may be used to endorse or promote products      *
*       derived from this software without specific prior written permission.     *
*                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
* DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    * 
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
***********************************************************************************/

#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/GCo.h"           
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_computational.h" 
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_files_strings.h"
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_fits.h"
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_plotting.h"
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_utils.h"
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"           
#include "/Users/stephanstiefelmaier/analysisSoftware/utils_sstiefel_2024/include/utils_TH1.h"

#include "/Users/stephanstiefelmaier/analysisSoftware/CommonHeaders/FittingGammaConversion.h"


#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include <iostream>
#include <string.h>
#include <vector>

Style_t markerStyleData    = 20;
  Size_t  markerSizeSpectrum = 1.0;
  Style_t lineStyleDataFit   = 1;
  Double_t lineWidthFit      = 1.0;
  
  Double_t markerSizeMC      = 1.0;
  Double_t markerSizeMBMC    = 1.0;
  Int_t markerStyleMCWW      = 20; // full circle
  Int_t markerStyleMCWOW     = 24; // hollow circle
  
  Color_t colorData           = kBlack;     //GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);
  Color_t colorMCMB           = kRed-4;
  Color_t colorMCASl          = kGreen+2;    //GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);
  Color_t colorMCASh          = kGreen-2;    //GetColorDefaultColor( fOptEnergy, "", cent[i], kFALSE);
  Color_t colorFit            = kOrange+7;
  Color_t colorOldFit         = kBlue+1;

// use binning of h2, Interpolate of h1
TH1D* makeRatioDiffBinnings(TH1D* h1, TH1D* h2, const char* name, const char* title="", float ymin=-999., float ymax=-999.){

    TH1D* h1overh2 = (TH1D*)h2->Clone(name);
    h1overh2->SetTitle(title);
    
    for (int iBin=1; iBin<=h2->GetNbinsX(); ++iBin){
        
        double binContentH2 = h2->GetBinContent(iBin);
        double relErrorH2 = binContentH2 ? h2->GetBinError(iBin) / binContentH2 : 0.;
        
        h1overh2->SetBinContent(iBin, binContentH2 ? h1->Interpolate(h2->GetBinCenter(iBin)) / binContentH2 : 0.);
        h1overh2->SetBinError(iBin, relErrorH2 * h1overh2->GetBinContent(iBin));
    }
    return h1overh2;
}

/* // DrawGammaLines will draw the lines in the histogram for you
    * startX - starting point of drawing in x
    * endX - end point of drawing in x
    * startY -starting point of drawing in y
    * endY - end point of drawing in y
    * linew - line width
    */
    void DrawGammaLines(Float_t startX, Float_t endX,
                    Float_t startY, Float_t endY,
                    Float_t linew, Float_t lineColor = 4, Style_t lineStyle = 1, Float_t opacity = 1.){
        TLine * l1 = new TLine (startX,startY,endX,endY);
        l1->SetLineColor(lineColor);
        l1->SetLineWidth(linew);
        l1->SetLineStyle(lineStyle);
        if (opacity != 1.)
            l1->SetLineColorAlpha(lineColor,opacity);

        l1->Draw("same");
    }
    
    //__________________________________________________________________________________________________________
    void DrawGammaSetMarker(    TH1* histo1,
                                Style_t markerStyle,
                                Size_t markerSize,
                                Color_t markerColor,
                                Color_t lineColor ) {
        histo1->SetMarkerStyle(markerStyle);
        histo1->SetMarkerSize(markerSize);
        histo1->SetMarkerColor(markerColor);
        histo1->SetLineColor(lineColor);
        histo1->GetYaxis()->SetLabelFont(42);
        histo1->GetXaxis()->SetLabelFont(42);
        histo1->GetYaxis()->SetTitleFont(62);
        histo1->GetXaxis()->SetTitleFont(62);
    }
    //__________________________________________________________________________________________________________
    void DrawGammaSetMarker(    TH1* histo1,
                                TString xtitle = "",
                                TString ytitle = "",
                                Style_t markerStyle = 20,
                                Size_t markerSize = 1,
                                Color_t markerColor = kBlack,
                                Color_t lineColor = kBlack,
                                double textsize = 0.045,
                                double labelsize = 0.045,
                                double xoffset = 1.,
                                double yoffset = 1. ) {
        histo1->SetTitle("");
        histo1->SetStats(0);
        histo1->SetMarkerStyle(markerStyle);
        histo1->SetMarkerSize(markerSize);
        histo1->SetMarkerColor(markerColor);
        histo1->SetLineColor(lineColor);
        histo1->GetYaxis()->SetLabelFont(42);
        histo1->GetXaxis()->SetLabelFont(42);
        histo1->GetYaxis()->SetLabelSize(labelsize);
        histo1->GetXaxis()->SetLabelSize(labelsize);
        histo1->GetYaxis()->SetTitleFont(62);
        histo1->GetXaxis()->SetTitleFont(62);
        histo1->GetYaxis()->SetTitleSize(textsize);
        histo1->GetXaxis()->SetTitleSize(textsize);
        if(!xtitle.EqualTo("")) histo1->GetXaxis()->SetTitle(xtitle);
        if(!ytitle.EqualTo("")) histo1->GetYaxis()->SetTitle(ytitle);
        histo1->GetXaxis()->SetTitleOffset(xoffset);
        histo1->GetYaxis()->SetTitleOffset(yoffset);
    }

    TString ReturnFullCollisionsSystem( TString fEnergyFlagOpt, TString dummyWUP = ""){
      dummyWUP.Length();  
      if(fEnergyFlagOpt.CompareTo("7TeV") == 0){
            return  "pp, #sqrt{#it{s}} = 7 TeV";
        } else if( fEnergyFlagOpt.CompareTo("7TeVSys") == 0) {
            return  "pp, #sqrt{#it{s}} = 7 TeV";
        } else if( fEnergyFlagOpt.BeginsWith("8TeV")) {
            return  "pp, #sqrt{#it{s}} = 8 TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeV") == 0) {
            return  "pp #sqrt{#it{s}} = 13 TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeVSys") == 0) {
            return  "pp #sqrt{#it{s}} = 13 TeV";
        } else if( fEnergyFlagOpt.Contains("13TeVMult")) {
          return  "pp #sqrt{#it{s}} = 13 TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeVNoB") == 0) {
          return  "pp #sqrt{#it{s}} = 13 TeV";
        } else if( fEnergyFlagOpt.CompareTo("13TeVLowB") == 0) {
            return  "pp, #sqrt{#it{s}} = 13 TeV (low B)";
        } else if( fEnergyFlagOpt.CompareTo("13TeVRBins") == 0) {
            return  "pp, #sqrt{#it{s}} = 13 TeV (RBins)";
        } else if( fEnergyFlagOpt.CompareTo("13TeVRBinsLowB") == 0) {
            return  "pp, #sqrt{#it{s}} = 13 TeV (RBins,LowB)";
        } else if( fEnergyFlagOpt.BeginsWith("5TeV") ) {
            return  "pp, #sqrt{#it{s}} = 5.02 TeV";
        } else if( fEnergyFlagOpt.CompareTo("900GeV") == 0) {
            return  "pp, #sqrt{#it{s}} = 900 GeV";
        } else if( fEnergyFlagOpt.CompareTo("2.76TeV") == 0) {
            return  "pp, #sqrt{#it{s}} = 2.76 TeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_2.76TeV") == 0) || (fEnergyFlagOpt.CompareTo("HI") == 0) ) {
            return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV";
        } else if( (fEnergyFlagOpt.CompareTo("PbPb_5.02TeV") == 0) ) {
            return "Pb-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
        } else if( (fEnergyFlagOpt.CompareTo("XeXe_5.44TeV") == 0) ) {
            return "Xe-Xe, #sqrt{#it{s}_{_{NN}}} = 5.44 TeV";
        } else if( fEnergyFlagOpt.CompareTo("pPb_5.023TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVCent") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.02TeV") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun2") == 0 || fEnergyFlagOpt.CompareTo("pPb_5.023TeVRun1") == 0 ) {
            return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 5.02 TeV";
        } else if( fEnergyFlagOpt.Contains("pPb_8TeV")) {
            return "p-Pb, #sqrt{#it{s}_{_{NN}}} = 8.16 TeV";
        } else {
            cout << "No correct collision system specification, has been given" << endl;
            return "";
        }
    }


void SetStyleHistoTH1ForGraphs( TH1* histo,
                                TString XTitle,
                                TString YTitle,
                                Size_t xLableSize,
                                Size_t xTitleSize,
                                Size_t yLableSize,
                                Size_t yTitleSize,
                                Float_t xTitleOffset    = 1,
                                Float_t yTitleOffset    = 1,
                                Int_t xNDivisions       = 510,
                                Int_t yNDivisions       = 510,
                                Font_t textFontLabel    = 42,
                                Font_t textFontTitle    = 62
                                ){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetYaxis()->SetLabelFont(textFontLabel);
    histo->GetXaxis()->SetLabelFont(textFontLabel);
    histo->GetYaxis()->SetTitleFont(textFontTitle);
    histo->GetXaxis()->SetTitleFont(textFontTitle);

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

void SetStyleHistoTH2ForGraphs( TH2* histo,
                                TString XTitle,
                                TString YTitle,
                                Size_t xLableSize,
                                Size_t xTitleSize,
                                Size_t yLableSize,
                                Size_t yTitleSize,
                                Float_t xTitleOffset    = 1,
                                Float_t yTitleOffset    = 1,
                                Int_t xNDivisions       = 510,
                                Int_t yNDivisions       = 510,
                                Font_t textFontLabel    = 42,
                                Font_t textFontTitle    = 62
                                ){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetYaxis()->SetLabelFont(textFontLabel);
    histo->GetXaxis()->SetLabelFont(textFontLabel);
    histo->GetYaxis()->SetTitleFont(textFontTitle);
    histo->GetXaxis()->SetTitleFont(textFontTitle);

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

void SetPlotStyle_() {
// 	const Int_t nRGBs = 7;
    const Int_t nRGBs = 5;
    const Int_t nCont = 255;

    Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
    gStyle->SetNumberContours(nCont);
}

void StyleSettingsThesis_( TString format = ""){
    //gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);   //show day and time
    gStyle->SetOptStat(0);  //show statistic
    gStyle->SetPalette(1,0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTextSize(0.5);
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetLabelOffset(0.006,"xyz");
    gStyle->SetTitleFontSize(0.04);
  
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(0.7,"x");
        
    gStyle->SetCanvasColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(1);

    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadLeftMargin(0.13);


    TGaxis::SetMaxDigits(3);
    gErrorIgnoreLevel=kError;

    if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
}
