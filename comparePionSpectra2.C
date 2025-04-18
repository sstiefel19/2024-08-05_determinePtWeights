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

#include "comparePionSpectra2.h"
#include "TSystem.h"

void doMultiRound(std::map<int, std::string> const &theMapBaseDirs,
                  std::string theMeson, 
                  std::string theCent,
                  std::map<std::string, tVPars const&> const &theMap,
                  std::vector<int> const &theRounds,
                  double theLeftMargin/*=0.25*/,
                  double theRightMargin/*=0.05*/,
                  std::string const &theDir/*=""*/){
    /* 
    the left column must have a different width to accomodate the y axis captions and 
    still have equally long x axis as the others*/
    int nCols = theRounds.size();
    double denom = 1. + (nCols-1)*(1.-theLeftMargin);
    double w1 = 1./denom;
    double wi = (1.-theLeftMargin)/denom;

    size_t lColumWidth = 200; // pixels

    // the master canvas
    std::string lNameC(Form("lCMultiRound_%s_%s", theMeson.data(), theCent.data()));
    auto &lCM = *new TCanvas(lNameC.data(), lNameC.data(), nCols*lColumWidth, 1000);
    lCM.SetMargin(0., 0., 0., 0.);

    // create n parallel columns in which the canvas from fitMesonAndWriteToFile will be drawn
    double lastBoarder = 0.;
    std::vector<TPad*> lPads;
    lPads.push_back(nullptr);
    for (int i=1; i<=nCols; ++i){
        double boarder = w1 + (i-1)*wi;
        TPad* p = new TPad(Form("%s_subpad_%d", lNameC.data(), i), Form("pad_%d", i), lastBoarder, 0., boarder, 1.);
        setMarginsToZero(*p);
        p->Draw();
        lPads.push_back(p);
        lastBoarder = boarder;
    }
    
    std::string lKey = theMeson + "_" + theCent;
    tVPars const &lVector = theMap.at(lKey);
    for (size_t iPos = 0; iPos < nCols; ++iPos){    
        bool isLast = iPos == nCols-1;
        int round = theRounds[iPos];
        TCanvas* cNx1_i = 
            fitMesonAndWriteToFile(
                theMapBaseDirs,   
                round,
                theCent,
                theMeson,
                lVector[round].fitFunc.data(),
                lVector[round].fitOpts.data(),
                lVector[round].tag.data(),
                lColumWidth,
                !iPos  ? theLeftMargin : 0 /*theLeftMargin*/,
                isLast ? theRightMargin : 0. /*theRightMargin*/,
                true /*verticallyTight*/,
                theDir);

        lPads[iPos+1]->cd();
        cNx1_i->DrawClonePad();
    }

    lCM.SaveAs(Form("%s/%s.pdf", theDir.data(), lNameC.data()));
    lCM.SaveAs(Form("~/repos_cloud/thesis_writing/figures/%s.pdf", lNameC.data()));
}

//====================================================================
TCanvas* 
    fitMesonAndWriteToFile(std::map<int, std::string> const &theMapBaseDirs,
                           int round, 
                           std::string eventCutNo, 
                           std::string meson, 
                           const char* fitFunction, 
                           const char* fitOption, 
                           const char* effiPlotLabel, 
                           size_t thePlotWidth,
                           double theLeftMargin/*=0.25*/,
                           double theRightMargin/*=0.05*/,
                           bool verticallyTight/*=true*/,
                           std::string theDir/*=""*/){

    // ============================================================
    // 1) retrieve all filenames, histos, and information
    std::string const lUniqueTag(
        Form("%d_round_%d_%s_%s", 
        TDatime().GetTime(), round, eventCutNo.data(), meson.data()));

    bool isPi0 = meson == "Pi0";
    bool isCentral = eventCutNo.substr(0, 3) == "101";

    Double_t minPtPlot  = isPi0 ?  0.3 : .9;
    Double_t maxPtPlot  = isPi0 ? 30.0 : 14.0;

    Double_t lYmin{isPi0 ? 0.7 : 0.7};
    Double_t lYmax{isPi0 ? 1.3 : 1.3};
    
    std::string lPhotMesCutNo("0d200009ab770c00amd0404000_0152101500000000");
    std::string lCutNo(Form("%s_%s", eventCutNo.data(), lPhotMesCutNo.data()));
    bool isPremerged = round==3;

    auto giveFilename = [&](int round, std::string dataMC, std::string fileSuff) {
        std::string lMCconfig(Form("MB-%s_separate", (round==2 || round==3) ? "bothASpremerged" : "AS"));
        return Form("%s/%s/mesons/%s/PbPb_5.02TeV/%s_%s_GammaConvV1Correction%s_%s.root", 
                    theMapBaseDirs.at(round).data(), lMCconfig.data(), lCutNo.data(), 
                    meson.data(), dataMC.data(), fileSuff.data(), lCutNo.data());
    };

    std::string filenameData(giveFilename(round, "data", ""));
    std::string filenameData_last(giveFilename(std::max(0, round-1), "data", ""));

    TH1D* lHistoData = (TH1D*)utils_files_strings::GetObjectFromPathInFile(
        filenameData.data(), 
        "CorrectedYieldTrueEff");
                
    // get counts directly from trainfile and transform to inv yield myself
    float lMeson2GammaBR = isPi0 ? 0.98798 : 0.3931;
    std::string eventCutNoMC(eventCutNo);
    eventCutNoMC.replace(5,2,"05");
    std::string eventCutNoMCAdd(eventCutNoMC);
    eventCutNoMCAdd.replace(6,1,"2");

    // theMBAS = one of  {"mb", "as", "as2"}
    // returns nullptr for invalid round, theMBAS combinations. 
    auto computeInvariantMCYieldFromTrainFile = 
        [&](std::string theMBAS, std::string const &theWoWeights=""){
        bool lIsMB = theMBAS=="mb";
        bool lIsAS2 = !lIsMB && (theMBAS=="as2");

        // default trainconfigs
        std::string lConfig(lIsMB 
            ? "994" 
            : isPi0 
                ? "997" 
                : "996"); 
        // special cases for first iterations
        if (!lIsMB && (round<3)){
            switch (round) {
                case 1: lConfig = isPi0 ? "995" : "996";
                case 2: lConfig = isPi0 ? "997" : "998";
                default: return static_cast<TH1*>(nullptr);     
            }
        }

        // catch invalid configs
        if ( (!round && !lIsMB) || (round<4 && lIsAS2)){
            const char* pref = "computeInvariantMCYieldFromTrainFile(): INFO:"; 
            printf("%s ((!round && !lIsMB) || (round<4 && lIsAS2)) is true, returning nullptr.\n", 
                pref);
            printf("%s (round, theMBAS.data(), lConfig.data()) : (%i, %s, %s)\n", 
                pref, round, theMBAS.data(), lConfig.data());
            return static_cast<TH1*>(nullptr);
        }

        std::string asTag(isPremerged ? "AS(l+h)" : lIsAS2 ? "ASl" : "ASh");
        std::string weightsTag(theWoWeights.size() ? "WOW" : "WW");
        std::string mcTag(Form("%s", lIsMB ? "MB" : asTag.data()));
        std::string lNewHistoName(Form("MC_%s_%s", mcTag.data(), weightsTag.data()));
        std::string lHistoTitleForLeg(Form("MC %s", mcTag.data()));
        std::string abac(isCentral ? "ab" : "ac");

        printf("SFS round, theMBAS.data(), asTag.data(), weightsTag.data(), mcTag.data(), lNewHistoName.data(), lIsMB, lIsAS2 :\n\t %d %s %s %s %s %s %d %d\n", round, theMBAS.data(), asTag.data(), weightsTag.data(), mcTag.data(), lNewHistoName.data(), lIsMB, lIsAS2);
               
        std::string lTrainSubDir(lIsMB 
            ? (round<2)                 // isMB &&
                ? "mc"                  //   round < 2
                : abac.data()           //   round >= 2

            : (round==1)                // as or as2 &&
                ? "LHC20g10"            //   round=1 
                : (round<4) 
                    ? "LHC20g10-LHC24a1-LHC24a2_merged"
                    : lIsAS2 
                        ? "LHC24a2"
                        : "LHC20g10");

        lTrainSubDir.append(lIsAS2 
            ? Form("/child%s_runlist_1", isCentral ? "1" : "2")
            : "");
        
        std::string lFname(Form("GCo_%s.root", lConfig.data()));
        if (lIsMB && round<2){
            lFname = Form("GammaConvV1_%s_%s-merged.root", lConfig.data(), abac.data());
        }
        if (!lIsMB && round==1) {
            lFname = Form("GammaConvV1_%s_%s.root", lConfig.data(), meson.data());
        }

        std::string lMainDir(Form("GammaConvV1_%s/", lConfig.data()));
        std::string lEventCutNo(lIsMB ? eventCutNoMC : eventCutNoMCAdd);
        GCo lGCo(Form("%s/trains/%s/%s",
                    theMapBaseDirs.at(round).data(), lTrainSubDir.data(),lFname.data()),
                lMainDir,
                lEventCutNo,
                lPhotMesCutNo);   
        
        float nEvents = ((TH1*)lGCo.GetFromESD("NEvents"))->GetBinContent(1); 
        TH1* lHistoMesonsInRap= (TH1*)lGCo.GetFromMC(Form("MC_%s%s_Pt", meson.data(), theWoWeights.data()));
        TH1* hInvYield_WW = lHistoMesonsInRap 
            ? utils_computational::TranformD2DPtDYYieldToInvariantYield(
                *utils_TH1::DivideTH1ByBinWidths(*lHistoMesonsInRap, "divW", nullptr, nullptr),
                "inv", lNewHistoName.data(), lHistoTitleForLeg.data(), 
                1./(nEvents*1.6*lMeson2GammaBR)) // 1.6 = deltaY
            : nullptr;
        return hInvYield_WW;
    };

    // compute all weighted and unweighted inv mc yields. The weighted ones should agree very well with last iterations' fits.
    std::vector<TH1*> vInvMCYields_ww;
    std::vector<TH1*> vInvMCYields_wow;
    std::vector<std::string> vMCs({"mb"});
    if (round){vMCs.push_back("as");}
    if (round>3){vMCs.push_back("as2");}

    for (auto const &mc : vMCs){
        auto computeAndInsert = [&](std::string theWoW){
            TH1 *h = computeInvariantMCYieldFromTrainFile(mc, theWoW);
            bool ww = !static_cast<bool>(theWoW.size());
            if (h) {
                std::vector<TH1*> &lVec = ww ? vInvMCYields_ww : vInvMCYields_wow;
                printf("Adding %s to %s\n", h->GetName(), ww ? "vInvMCYields_ww" : "vInvMCYields_wow"); 
                lVec.push_back(h); 
            }
        };
        for (auto const &iWeightsQuali : std::vector<std::string>({"", "_WOWeights"})){
            computeAndInsert(iWeightsQuali);
        }
    }
    TH1 *hInvMCYield_mc_mb_nw = vInvMCYields_wow.at(0);
    
    // get fit and histo from last iteration
    std::string fitNameInFile(Form("%s_Data_5TeV_%s0", meson.data(), eventCutNo.substr(0,5).data()));
    size_t lastIt = max(0, round-1);
    std::string suff((lastIt==1) ? "0b" : Form("%zu", lastIt));
    std::string fname_weightsFile_last(Form(
        "~/2024/2024-08-05_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it%s.root", 
        suff.data()));
    printf("opening last iterations weights file %s ..\n", fname_weightsFile_last.data());
    TF1* lastItFit = (TF1*)utils_files_strings::GetObjectFromPathInFile(
        fname_weightsFile_last, 
        fitNameInFile.data());
    
    // 2) ======================= do this iterations fit =================================
    std::string fit_data_name(Form("%s_Data_5TeV_%s_it%d", meson.data(), eventCutNo.substr(0,5).data(), round)); // need the it in the name here so we dont get many  objects with the same name when doing more than one it
    TF1* fit_data = FitObject(fitFunction, fit_data_name.data(), meson.data(), NULL, minPtPlot, maxPtPlot);            
    TGraphAsymmErrors* graphYield_data = new TGraphAsymmErrors(lHistoData);
    graphYield_data->Fit(fit_data, fitOption, "", minPtPlot, maxPtPlot);

    // 2.1) fit minimum bias MC
    // std::string fit_mc_mb_name(Form("%s_LHC20e3a_5TeV_%s_it%d", meson.data(), eventCutNo.substr(0,5).data(), round)); // need the it in the name here so we dont get many  objects with the same name when doing more than one it
    // TF1* fit_mc_mb_nw = FitObject(fitFunction, fit_mc_mb_name.data(), meson.data(), NULL, minPtPlot, maxPtPlot);            
    // TGraphAsymmErrors* graphYield_mc_mb = new TGraphAsymmErrors(hInvMCYield_mc_mb_nw);
    // graphYield_mc_mb->Fit(fit_mc_mb_nw, fitOption, "", minPtPlot, maxPtPlot);

    // 2.2 try local exponential interpolations
    TF1 &f_exp_inter_mc_mb_nw = 
        utils_TH1::GlobalPieceWiseExponentialInterpolation(
            Form("%s_exp_inter", 
                 hInvMCYield_mc_mb_nw->GetName()),
            *hInvMCYield_mc_mb_nw);
    
    // 3) ========================= plotting =============================================
    bool isFirstCol = theLeftMargin;
    size_t lNrows = 5;
    
    float lXlabelSize = 0.1;
    float lXtitleSize = 0.1;
    float lXtitleOffset = 1.2;

    float lYlabelSize = 0.1;
    float lYtitleSize = 0.08;
    float lYtitleOffset = 1.7;
    float lRatioYtitleOffsets = 1.3;

    float lLegendTextSize = 0.08;

    gStyle->SetOptTitle(0); // disables histo titles as title on subpads
    gStyle->SetOptStat(0); // 

    TCanvas &cOneIt = *new TCanvas(Form("canvas_%s_%i", fit_data->GetName(), round), 
                                  Form("canvas_%s_%i", fit_data->GetName(), round), 
                                  thePlotWidth, 1000);

    cOneIt.SetMargin(0., 0., 0., 0.);    
    cOneIt.Divide(1, lNrows, 0., verticallyTight ? 0. : 0.1);

    auto getNextTab = [&](){
        int newPadNo = gPad->GetNumber() + 1;
        cOneIt.cd(newPadNo);
        bool isTop = newPadNo==1;
        bool isBottom = newPadNo==lNrows;
        TVirtualPad &p = *gPad;
        p.SetLogx();
        if (isTop) { p.SetLogy(); }
        p.SetMargin(theLeftMargin, theRightMargin, isBottom ? 0.3 : 0., isTop ? 0.03 : 0.);
        p.SetTicks(2, 2); // tx=2 enables ticks on top and bottom; ty=2 enables ticks on left and right
        return &p;
    };

    //////////////////////////////////////////////////////////////////////////////////
    // fit data and plot both together
    cout << "============================= pad 1 ===========================================\n";
    auto &pad1 = *getNextTab();    
    TH1F &histo1DSpectra = *new TH1F("histo1DSpectra", "histo1DSpectra",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( 
        &histo1DSpectra, 
        "#it{p}_{T} (GeV/#it{c})", 
        "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
        0.,  // xLableSize
        0., // xTitleSize
        isFirstCol ? lYlabelSize : 0.,  // yLableSize
        isFirstCol ? lYtitleSize : 0., // yTitleSize
        0.,  // xTitleOffset
        lYtitleOffset);  // yTitleOffset
    
    if (isFirstCol){
        histo1DSpectra.GetXaxis()->SetLabelOffset(-0.01);
        histo1DSpectra.GetYaxis()->CenterTitle(true);
        histo1DSpectra.GetYaxis()->SetLabelOffset(0.01);
    }
    
    histo1DSpectra.GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DSpectra.GetYaxis()->SetRangeUser(1E-8, 5E3);
    histo1DSpectra.DrawCopy();

    float lLegendNDC_x1 = (round==3) ? 0.51 : 0.61;
    auto xnew = [theLeftMargin](float xold){return theLeftMargin + (1.-theLeftMargin)*xold;};
    auto legend_pad1 = new TLegend(xnew(lLegendNDC_x1), 0.51, xnew(0.98), .96);
    legend_pad1->SetBorderSize(0);
            
    // plot for all iterations
    lHistoData->SetTitle(Form("%s yields for %s", meson.data(), eventCutNo.data()));
    utils_plotting::DrawAndAdd(*lHistoData,"same", colorData, 1.0, legend_pad1, "Data", "lep", lLegendTextSize, true, markerStyleData, 1.0);
    // keep this for reference
    // utils_plotting::DrawAndAdd(*lHistoMBonly_oldBin_WOW, "same", colorMCMB, 1.0, legend_pad1, "MC MB WoW", "lep", lLegendTextSize, true, markerStyleMCWOW, 1.0);
    
    typedef std::pair<Color_t, Style_t> MPair;
    auto getColorAndMarkerForHisto = [](TH1 &h){
        TString hname(h.GetName());
        Color_t lColor = hname.Contains("MB") 
            ? colorMCMB 
            : hname.Contains("+") 
                ? kGreen
                : hname.Contains("ASl")
                    ? kBlue
                    : kOrange; 
        Style_t lStyle = hname.Contains("WW") ? markerStyleMCWW : markerStyleMCWOW;
        return MPair({lColor, lStyle});
    };
    auto plotVector = [&](std::vector<TH1*> const &theVector){
        for (TH1* ih : theVector){
            if (!ih) { continue; }
            TH1& h = *ih;
            MPair lPair = getColorAndMarkerForHisto(h);           
            bool ww = lPair.second==markerStyleMCWW;
            utils_plotting::DrawAndAdd(h, "same", lPair.first, 1.0, ww ? legend_pad1 : nullptr, h.GetTitle(), "lep", lLegendTextSize, true, lPair.second, 1.0);
        }
    };

    // plot MC inv yields without weights
    plotVector(vInvMCYields_wow);
    
    // and those with weights
    if (round){
        plotVector(vInvMCYields_ww);
    }

    // move on to fits
    if (lastItFit){
        utils_plotting::DrawAndAdd(*lastItFit, "same", kBlue, 3.0, legend_pad1, "last Fit", "l", lLegendTextSize, true);
    }
    fit_data->SetRange(minPtPlot, maxPtPlot);
    // utils_plotting::DrawAndAdd(*fit_data,     "same", colorFit, 3.0, legend_pad1, "Fit Data", "l", lLegendTextSize, true);
    // utils_plotting::DrawAndAdd(*fit_mc_mb_nw, "same", colorFit+2, 3.0, legend_pad1, "Fit MC MB NW", "l", lLegendTextSize, true);
    utils_plotting::DrawAndAdd(f_exp_inter_mc_mb_nw, "same", colorFit+2, 3.0, legend_pad1, "exp inter MC MB NW", "l", lLegendTextSize, true);
        
    auto centString = [](std::string& evtCutNo){
        std::string& e = evtCutNo;
        if (e[0]=='1') { 
            std::string s(e[1]=='0' ? "" : Form("%c", e[1]));
            return Form("%s0-%c0%%",s.data(), e[2]);
        }
        else {return Form("%s",e.data());}
    };
    
    TString collisionSystem = ReturnFullCollisionsSystem("PbPb_5.02TeV");
    std::string mesDec(isPi0 ? "#pi^{0}" : "#eta");
    mesDec.append(" #rightarrow #gamma #gamma");
    
    // add text in plot
    std::vector<std::string const> const lPaveTextLines({ 
        Form("%s  %s", collisionSystem.Data(), centString(eventCutNo)),
        mesDec.data(),
        "", 
        Form("iteration: %d", round), 
        Form("%s in %s", meson.data(), eventCutNo.data()),
        effiPlotLabel});

    TPaveText &pav = utils_plotting::SetupTPaveText(
        xnew(0.045), 0.01, xnew(0.39), 0.48, lPaveTextLines, lLegendTextSize); // x1,y1,x2,y2

    //////////////////////////////////////////////////////////////////////////
    // current spectra / last iteration fit
    cout << "==================== PAD 2 ===================================\n";
    auto &pad2 = *getNextTab();

    TH1F &histo1DRatio = *new TH1F("histo1DRatio", "histo1DRatio",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs(
        &histo1DRatio, 
        "#it{p}_{T} (GeV/#it{c})", 
        "This MC over last Fit", 
        0.,  // xLableSize
        0., // xTitleSize
        isFirstCol ? lYlabelSize : 0.,  // yLableSize
        isFirstCol ? lYtitleSize : 0., // yTitleSize
        0.,  // xTitleOffset
        lRatioYtitleOffsets);  // yTitleOffset
    
    histo1DRatio.GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio.GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio.GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio.GetYaxis()->CenterTitle(true);
    histo1DRatio.GetYaxis()->SetRangeUser(0.9, 1.1);        
    histo1DRatio.DrawCopy();

    auto leg2 = new TLegend(xnew(0.144),0.73,xnew(0.44),0.92);
    leg2->SetBorderSize(0);
    leg2->Draw();

    if (lastItFit)
    {
        // get all MC histos from pad1 and calculate ratio to last it's fit (the fit that was used for weighting the MCs)
        TList* entries = legend_pad1->GetListOfPrimitives();
        for (int i = 0; i < entries->GetEntries(); ++i) {
            TLegendEntry* entry = (TLegendEntry*)entries->At(i);
            if (entry && std::string(entry->GetLabel()).find("MC") != std::string::npos) {
                TObject* obj = entry->GetObject();
                if (!obj) {continue;}
                printf("Processing entry i: %d objectName: %s label: %s from legend_pad1\n", i, obj->GetName(), entry->GetLabel());
                if (obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("MC") == 0) {
                    TH1* histoMC = dynamic_cast<TH1*>(obj);
                    if (!histoMC) {continue;}
                    TH1 &hHistoRatioToFit = *utils_TH1::DivideTH1ByTF1(*histoMC, *lastItFit, nullptr, nullptr, kTRUE /*theIntegrateTF1*/);
                    hHistoRatioToFit.SetName(Form("%s_mover_%s", histoMC->GetName(), lastItFit->GetName()));
                    utils_plotting::DrawAndAdd(hHistoRatioToFit, "same", histoMC->GetLineColor(), 1.0, leg2, entry->GetLabel(), "lp", lLegendTextSize);
                }
            }
        }
    }    
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);

    //////////////////////////////////////////////////////////////////////////
    // compare this effi to previous effis
    cout << "==================== PAD 3 ===================================\n";
    auto &pad3 = *getNextTab();

    histo1DRatio.GetYaxis()->SetTitle("this over last efficiency");
    histo1DRatio.DrawCopy();
    
    if (round) {
        printf("SFS round = %d\n", round);
        // draw ratio of effi/effiLast
        TH1 &lHistoTrueEffi     = *(TH1*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "TrueMesonEffiPt");
        TH1 *lHistoTrueEffiLast = round 
            ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameData_last.data(), 
                                                                 "TrueMesonEffiPt")
            : static_cast<TH1*>(nullptr);

        // utils_TH1::PrintBinsWithErrors(lHistoTrueEffi);    
        std::pair<TH1&, TH1&> &lBinningsAligned = 
            *utils_TH1::AlignBinnings(lHistoTrueEffi, *lHistoTrueEffiLast);
        
        TH1* lEffiOverEffiLast = utils_TH1::DivideTH1ByTH1(lBinningsAligned.first, 
                                                            lBinningsAligned.second,
                                                            nullptr,
                                                            "EffiOverLastEffi");
        // utils_TH1::PrintBinsWithErrors(*lEffiOverEffiLast);

        DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
        lEffiOverEffiLast->Draw("same");
    }

    //////////////////////////////////////////////////////////////////////////
    // show quality of fits
    cout << "==================== PAD 4 ===================================\n";

    auto &pad4 = *getNextTab();
    histo1DRatio.GetYaxis()->SetTitle("this Data over its Fit");
    histo1DRatio.GetYaxis()->SetRangeUser(lYmin,lYmax);
    histo1DRatio.DrawCopy();

    // this fit over this data
    TH1* hHistoRatioDataToFit = utils_TH1::DivideTH1ByTF1(
        *lHistoData, *fit_data, Form("hHistoRatioDataToFit_it%d", round), nullptr, false/*integrateFunction*/);
    
    // this mc_mb_nw fit over its histo
    // TH1* hHistoRatioMCMBNWToFit = utils_TH1::DivideTH1ByTF1( // todo: check integreate function yes or no also for data fits!!
    //     *hInvMCYield_mc_mb_nw, *fit_mc_mb_nw, Form("hHistoRatio_mc_mb_nw_ToFit_it%d", round), nullptr, false/*integrateFunction*/);

    // this mc_mb_nw's exponential interpolation over mc_mb_nw
    TH1* hHistoRatioMCMBNWExpInterToHisto = utils_TH1::DivideTH1ByTF1( // todo: check integreate function yes or no also for data fits!!
        *hInvMCYield_mc_mb_nw, f_exp_inter_mc_mb_nw, Form("hHistoRatio_mc_mb_nw_expInter_it%d", round), nullptr, false/*integrateFunction*/);
        
    auto *leg4 = utils_plotting::GetLegend(xnew(0.16),0.6,xnew(0.44),0.86);
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    utils_plotting::DrawAndAdd(*hHistoRatioDataToFit, "same", colorFit, 3.0, leg4, "this data over its fit", "l", lLegendTextSize, true);
    // utils_plotting::DrawAndAdd(*hHistoRatioMCMBNWToFit, "same", colorFit+2, 3.0, leg4, "this mb mc over its fit", "l", lLegendTextSize, true);
    utils_plotting::DrawAndAdd(*hHistoRatioMCMBNWExpInterToHisto, "same", colorFit+5, 3.0, leg4, "this mb mc over its exp inter", "l", lLegendTextSize, true);

    //////////////////////////////////////////////////////////////////////////
    // ratio of this weighted MCs over this data. They differ only as much as this over last true efficiency
    cout << "==================== PAD 5 ===================================\n";
    auto &pad5 = *getNextTab();
    pad5.SetTicks(1,2);
    // // Hide labels on the top X-axis
    // histo1DRatio.GetXaxis()->SetLabelOffset(999); // Move top labels far away to effectively hide them

    histo1DRatio.GetXaxis()->SetLabelSize(lXlabelSize);
    histo1DRatio.GetXaxis()->SetTitleSize(lXtitleSize);
    histo1DRatio.GetXaxis()->SetTitleOffset(lXtitleOffset);
    
    histo1DRatio.GetYaxis()->SetRangeUser(lYmin,lYmax);        
    histo1DRatio.GetYaxis()->SetTitle("MC over Data Fit");
    histo1DRatio.DrawCopy();

    auto leg5 = new TLegend(xnew(0.144),0.73,xnew(0.44),0.92);
    leg5->SetBorderSize(0);
    leg5->Draw();

    TFile file_debug("file_debug.root", "update");
    auto &dir = *file_debug.mkdir(lUniqueTag.data());
    dir.cd();
    lHistoData->Write();
    // get all MC histos from pad1 and calculate ratio to data
    TList* entries = legend_pad1->GetListOfPrimitives();
    for (int i = 0; i < entries->GetEntries(); ++i) {
        TLegendEntry* entry = (TLegendEntry*)entries->At(i);
        if (entry && std::string(entry->GetLabel()).find("MC") != std::string::npos) {
            TObject* obj = entry->GetObject();
            if (!obj) {continue;}
            printf("Processing entry i: %d objectName: %s label: %s from legend_pad1\n", i, obj->GetName(), entry->GetLabel());
            if (obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("MC") == 0) {
                TH1* histoMC = dynamic_cast<TH1*>(obj);
                if (!histoMC) {continue;}
                // TH1& histoMC_dataBin = (histoMC->GetNbinsX() == lHistoData->GetNbinsX()) 
                //     ? *histoMC 
                //     : *utils_TH1::RebinInvariantYieldHistogram(*histoMC, *lHistoData, "rebInv"); 
                // printf("SFS dividing MC %s by data %s\n", histoMC_dataBin.GetName(), lHistoData->GetName());
                
                // divide by datafit
                TH1 &ratioMCtoFit = *utils_TH1::DivideTH1ByTF1(
                    *histoMC, *fit_data, nullptr, nullptr, true /*theIntegrateTF1*/);                

                // histoMC->Dump();
                // ratioMCtoFit->Dump();
                utils_plotting::DrawAndAdd(
                    ratioMCtoFit, "same", histoMC->GetLineColor(), 3.0, 
                    leg5, histoMC->GetName()/*theLegString*/, "p", lLegendTextSize, true /*theLegDrawAlready*/, // when theLegString is empty,  objects name is used
                    histoMC->GetMarkerStyle(), histoMC->GetMarkerSize()); 

                dir.cd();
                histoMC->Write();
                ratioMCtoFit.Write();
                // histoMC_dataBin.Write();
                // ratioHist->Write();
            }
        }
    }
    file_debug.Close();    
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
        
    cout << "================= saving to file and pdfs ========================\n";
    std::string lSinglesDir(theDir + "/singles");
    gSystem->mkdir(lSinglesDir.data(), true /*recursive*/);
    cOneIt.SaveAs(Form("%s/%s_%s_it%d.pdf", lSinglesDir.data(), eventCutNo.data(), meson.data(), round));
    
    // save to file
    TFile* hfile = new TFile(Form("%s/MCSpectraInputPbPb_Stephan_it%d.root", theDir.data(), round),"UPDATE");    
    fit_data->Write();
    cOneIt.Write();
    hfile->Write();
    hfile->Close();
    hfile->Delete();    

    return &cOneIt;
}

//====================================================================   
void setMarginsToZero(TVirtualPad& vpad){
    vpad.SetLeftMargin(0.0);
    vpad.SetRightMargin(0.0);
    vpad.SetTopMargin(0.0);
    vpad.SetBottomMargin(0.0);
}

