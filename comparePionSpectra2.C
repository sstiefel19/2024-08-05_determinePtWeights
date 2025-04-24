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
        int theRound = theRounds[iPos];

        bool isPi0 = theMeson=="Pi0";        
        Double_t lXminFit = isPi0 ?  0.6 : 0.; // exclude 0.4-0.6GeV/c bin
        Double_t lXmaxFit = isPi0 ? 25 : 0.;

        TCanvas* cNx1_i = 
            fitMesonAndWriteToFile(
                theMapBaseDirs,   
                theRound,
                theCent,
                theMeson,
                lVector[theRound].fitFunc.data(),
                lVector[theRound].fitOpts.data(),
                lVector[theRound].tag.data(),
                lColumWidth,
                !iPos  ? theLeftMargin : 0 /*theLeftMargin*/,
                isLast ? theRightMargin : 0. /*theRightMargin*/,
                true /*verticallyTight*/,
                theDir,
                lXminFit ? &lXminFit : nullptr,
                lXmaxFit ? &lXmaxFit : nullptr);

        lPads[iPos+1]->cd();
        cNx1_i->DrawClonePad();
    }

    lCM.SaveAs(Form("%s/%s.pdf", theDir.data(), lNameC.data()));
    lCM.SaveAs(Form("~/repos_cloud/thesis_writing/figures/%s.pdf", lNameC.data()));
}

//====================================================================
TCanvas* 
    fitMesonAndWriteToFile(std::map<int, std::string> const &theMapBaseDirs,
                           int         theRound, 
                           std::string thEventCutNo, 
                           std::string theMeson, 
                           std::string theFitFunction, 
                           std::string theFitOption, 
                           std::string theEffiPlotLabel, 
                           size_t thePlotWidth,
                           double theLeftMargin  /* = 0.25*/,
                           double theRightMargin /* = 0.05*/,
                           bool verticallyTight  /* = true*/,
                           std::string theDir    /* = ""*/,
                           Double_t *theMinPtFit /* = nullptr*/,
                           Double_t *theMaxPtFit /* = nullptr*/,
                           bool      theSaveDNDPT /* = false*/){

    // ============================================================
    // 1) retrieve all filenames, histos, and information
    std::string const lUniqueTag(
        Form("%d_round_%d_%s_%s", 
        TDatime().GetTime(), theRound, thEventCutNo.data(), theMeson.data()));

    bool isPi0 = theMeson == "Pi0";
    bool isCentral = thEventCutNo.substr(0, 3) == "101";

    Double_t lMinPtPlot  = isPi0 ?  0.3 : .9;
    Double_t lMaxPtPlot  = isPi0 ? 30.0 : 14.0;

    Double_t lMinPtFit = theMinPtFit ? *theMinPtFit : lMinPtPlot; 
    Double_t lMaxPtFit = theMaxPtFit ? *theMaxPtFit : lMaxPtPlot; 

    Double_t lYmin_ratio{isPi0 ? 0.8 : 0.7};
    Double_t lYmax_ratio{isPi0 ? 1.2 : 1.3};

    Double_t lYmin_ratio_pad2{isPi0 ? 0.9 : 0.9};
    Double_t lYmax_ratio_pad2{isPi0 ? 1.1 : 1.1};
    
    std::string lPhotMesCutNo("0d200009ab770c00amd0404000_0152101500000000");
    std::string lCutNo(Form("%s_%s", thEventCutNo.data(), lPhotMesCutNo.data()));
    bool isPremerged = (theRound==2 || theRound==3);

    auto giveFilename = [&](int theRound, std::string dataMC, std::string fileSuff) {
        std::string lMCconfig(Form("MB-%s_separate", isPremerged ? "bothASpremerged" : "AS"));
        return Form("%s/%s/mesons/%s/PbPb_5.02TeV/%s_%s_GammaConvV1Correction%s_%s.root", 
                    theMapBaseDirs.at(theRound).data(), lMCconfig.data(), lCutNo.data(), 
                    theMeson.data(), dataMC.data(), fileSuff.data(), lCutNo.data());
    };

    std::string filenameData(giveFilename(theRound, "data", ""));
    std::string filenameData_last(giveFilename(std::max(0, theRound-1), "data", ""));

    TH1D* lHistoData_inv = (TH1D*)utils_files_strings::GetObjectFromPathInFile(
        filenameData, 
        "CorrectedYieldTrueEff",
        "clone",
        "CorrectedYieldTrueEff",
        "CorrectedYieldTrueEff");
                    
    // get counts directly from trainfile and transform to inv yield myself
    float lMeson2GammaBR = isPi0 ? 0.98798 : 0.3931;
    std::string eventCutNoMC(thEventCutNo);
    eventCutNoMC.replace(5,2,"05");
    std::string eventCutNoMCAdd(eventCutNoMC);
    eventCutNoMCAdd.replace(6,1,"2");

    // theMBAS = one of  {"mb", "as", "as2"}
    // returns nullptr for invalid theRound, theMBAS combinations. 
    auto computeInvariantMCYieldFromTrainFile = 
        [&](std::string theMBAS, std::string const &theWoWeights=""){
        bool lIsMB = theMBAS=="mb";
        bool lIsAS2 = !lIsMB && (theMBAS=="as2");

        printf("line 163: computeInvariantMCYieldFromTrainFile(): Info: called with theRound = %d, theMBAS.data() = %s, theWoWeights.data() = %s, lIsMB = %d, lIsAS2 = %d\n",
            theRound, theMBAS.data(), theWoWeights.data(), lIsMB, lIsAS2);

        // default trainconfigs
        // mb 994, pi0 997 eta 995
        std::string lConfig(lIsMB 
            // ? isCentral  // had those for testing purposes: these configs were without exp inter
            //     ? "5100"
            //     : "5130" 
            ? "994" 
            : isPi0 
                ? "997" 
                : "996"); 
        // special cases for first iterations and starting at 8
        printf("SFS line 173: theRound = %d\n", theRound);
        if (!lIsMB && ((theRound < 3) || (theRound >7))){
            printf("SFS line 175: theRound = %d\n", theRound);
            switch (theRound) {
                case 1: lConfig = isPi0 ? "995" : "996"; break;
                // case 2: lConfig = isPi0 ? "997" : "998"; break;
                case 8: lConfig = isPi0 ? "997" : "995"; break;
                // case 8: lConfig = isPi0 ? "997" : "995"; break;
                
                default: { 
printf("line 179\n");
                    return static_cast<TH1*>(nullptr);     
                }
            }
        }

        // catch invalid configs
        if ( (!theRound && !lIsMB) || (theRound<4 && lIsAS2)){
            const char* pref = "computeInvariantMCYieldFromTrainFile(): INFO:"; 
            printf("%s ((!theRound && !lIsMB) || (theRound<4 && lIsAS2)) is true, returning nullptr.\n", 
                pref);
            printf("%s (theRound, theMBAS.data(), lConfig.data()) : (%i, %s, %s)\n", 
                pref, theRound, theMBAS.data(), lConfig.data());
            return static_cast<TH1*>(nullptr);
        }

        std::string asTag(isPremerged ? "AS(l+h)" : lIsAS2 ? "ASl" : "ASh");
        std::string weightsTag(theWoWeights.size() ? "WOW" : "WW");
        std::string mcTag(Form("%s", lIsMB ? "MB" : asTag.data()));
        std::string lNewHistoName(Form("MC_%s_%s", mcTag.data(), weightsTag.data()));
        std::string lHistoTitleOnly(Form("MC %s", mcTag.data()));
        std::string lHistoTitle_wAxisCaptions(lHistoTitleOnly 
                    + ";pT (GeV/c);#frac{1}{2#pi N_ev pT} #frac{d^2 N}{dpT dy}"); 
            
        std::string abac((theRound>7) ?
              "MBabc"
            : isCentral ? 
                "ab" 
              : "ac");

        printf("SFS 208: theRound, theMBAS.data(), asTag.data(), weightsTag.data(), mcTag.data(), lNewHistoName.data(), lIsMB, lIsAS2 :\n\t %d %s %s %s %s %s %d %d\n", 
            theRound, theMBAS.data(), asTag.data(), weightsTag.data(), mcTag.data(), lNewHistoName.data(), lIsMB, lIsAS2);
               
        std::string lTrainSubDir(lIsMB 
            ? (theRound<2)              // isMB &&
                ? "mc"                  //          theRound < 2
                : abac.data()           //          theRound >= 2

            : (theRound==1)             // !MB [<=> as or as]  &&
                ? "LHC20g10"            //                         theRound=1 
                : (theRound<4)          // (!MB && !theRound==1) &&                        
                    ? "LHC20g10-LHC24a1-LHC24a2_merged" //           theRound < 4
                    : lIsAS2            // (!MB && theRound >= 4) &&               
                        ? "LHC24a2"     //                           lIsAS2
                        : "LHC20g10");  //                           !IsAS2  

        lTrainSubDir.append(lIsAS2 
            ? Form("/child%s_runlist_1", isCentral ? "1" : "2")
            : "");
        
        std::string lFname(Form("GCo_%s.root", lConfig.data()));
        if (lIsMB && theRound<2){
            lFname = Form("GammaConvV1_%s_%s-merged.root", lConfig.data(), abac.data());
        }
        if (!lIsMB && theRound==1) {
            lFname = Form("GammaConvV1_%s_%s.root", lConfig.data(), theMeson.data());
        }

        std::string lMainDir(Form("GammaConvV1_%s/", lConfig.data()));
        std::string lEventCutNo(lIsMB ? eventCutNoMC : eventCutNoMCAdd);
        GCo lGCo(Form("%s/trains/%s/%s",
                      theMapBaseDirs.at(theRound).data(), lTrainSubDir.data(),lFname.data()),
                lMainDir,
                lEventCutNo,
                lPhotMesCutNo,
                true /*_keepFileOpen*/);   
        
        float nEvents = ((TH1*)lGCo.GetFromESD("NEvents"))->GetBinContent(1); 
        TH1* lHistoMesonInRap= (TH1*)lGCo.GetFromMC(Form("MC_%s%s_Pt", theMeson.data(), theWoWeights.data()));
        // TH1* lHistoMesonInRap_daughtersInAcc= (TH1*)lGCo.GetFromMC(Form("MC_%sInAcc%s_Pt", theMeson.data(), theWoWeights.data()));
        TH1* hInvYield_WW = lHistoMesonInRap 
            ? utils_computational::TranformD2DPtDYYieldToInvariantYield(
                *utils_TH1::DivideTH1ByBinWidths(*lHistoMesonInRap, "divW", nullptr, nullptr),
                nullptr, 
                lNewHistoName.data(), 
                lHistoTitle_wAxisCaptions.data(), 
                1./(nEvents*1.6*lMeson2GammaBR)) // 1.6 = deltaY
            : nullptr;
        return hInvYield_WW;
    };

    // compute all weighted and unweighted inv mc yields. The weighted ones should agree very well with last iterations' fits.
    std::vector<TH1*> vInvMCYields_ww;
    std::vector<TH1*> vInvMCYields_wow;
    std::vector<std::string> vMCs({"mb"});
    if (theRound){vMCs.push_back("as");}
    if (theRound>3){vMCs.push_back("as2");}

    for (auto const &mc : vMCs){
        auto computeInvariantMCYieldFromTrainFile_andInsert = [&](std::string theWoW){
            TH1 *h = computeInvariantMCYieldFromTrainFile(mc, theWoW);
            bool ww = !static_cast<bool>(theWoW.size());
            if (!h){
                printf("WARNING: computeInvariantMCYieldFromTrainFile_andInsert(): computeInvariantMCYieldFromTrainFile() returned nullptr for mc = %s. Skipping iteration.\n", 
                    mc.data());
            }
            std::vector<TH1*> &lVec = ww ? vInvMCYields_ww : vInvMCYields_wow;
            printf("Adding %s to %s\n", 
                h ? h->GetName() : "nullptr", 
                ww ? "vInvMCYields_ww" : "vInvMCYields_wow"); 
            lVec.push_back(h); 
        };
    
        for (auto const &iWeightsQuali : std::vector<std::string>({"", "_WOWeights"})){
            printf("265\n");
            computeInvariantMCYieldFromTrainFile_andInsert(iWeightsQuali);
        }
    }
    
    TH1 *hInvMCYield_mc_mb_nw = vInvMCYields_wow.size() 
        ?  vInvMCYields_wow.at(0) 
        : static_cast<TH1*>(nullptr);

    TH1 *hVarMCYield_mc_mb_nw = hInvMCYield_mc_mb_nw 
        ?  utils_computational::TranformInvariantYieldToD2DPtDY(*hInvMCYield_mc_mb_nw)
        :  static_cast<TH1*>(nullptr);
    
    // get fit and histo from last iteration
    std::string fitNameInFile(Form("%s_Data_5TeV_%s0", theMeson.data(), thEventCutNo.substr(0,5).data()));
    size_t lastIt = max(0, theRound-1);
    std::string suff((lastIt==1) ? "0b" : Form("%zu", lastIt));
    std::string fname_weightsFile_last(Form(
        "~/2024/2024-08-05_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it%s.root", 
        suff.data()));
    printf("opening last iterations weights file %s ..\n", fname_weightsFile_last.data());
    TF1* lastItFit = (TF1*)utils_files_strings::GetObjectFromPathInFile(
        fname_weightsFile_last, 
        fitNameInFile,
        Form("%s_lastItFit", fitNameInFile.data()),
        Form("%s_lastItFit", fitNameInFile.data()));

    // 2) ======================= do this iterations fit =================================
    
    // first the 'tradional' way
    std::string fit_data_inv_name(Form("%s_Data_5TeV_%s0_it%d", theMeson.data(), thEventCutNo.substr(0,5).data(), theRound)); // need the it in the name here so we dont get many  objects with the same name when doing more than one it
    TF1 *lFit_data_inv = FitObject(theFitFunction.data(), fit_data_inv_name.data(), theMeson.data(), NULL, lMinPtPlot, lMaxPtPlot); 
    lHistoData_inv->Fit(lFit_data_inv, theFitOption.data(), "", lMinPtFit, lMaxPtFit);

    // then new
    TH1 &lHistoData_var = *utils_computational::TranformInvariantYieldToD2DPtDY(*lHistoData_inv);

    TF1 &lFit_data_var = utils_TF1::TransformInvYieldToDNDPT(Form("%s_multX", fit_data_inv_name.data()), *lFit_data_inv);
    lFit_data_var.SetRange(lMinPtPlot, lMaxPtPlot);
    lHistoData_var.Fit(&lFit_data_var, Form("I%s", theFitOption.data()), "", lMinPtFit, lMaxPtFit);
    lFit_data_var.SetRange(lMinPtPlot, lMaxPtPlot);

    // transform back for plotting -  should be close to lFit_data_inv
    TF1 &lFit_data_var_trans_inv = utils_TF1::TransformDNDPTToInvYield(""/*triggers autom. suffix*/, lFit_data_var);

    // compare both in invariant form
    TF1 &lRatioFits_int_over_old = utils_TF1::TF1Division("lRatioFits_int_over_old", lFit_data_var_trans_inv, *lFit_data_inv, false /*theCheckRanges*/);
 
    // 2.2 test local exponential interpolations for MB MC inv
    auto getExpInter = [](TH1 *theH){
        return theH
            ?  &utils_TH1::GlobalPieceWiseExponentialInterpolation(
                    Form("%s_exp_inter", theH->GetName()), *theH)
            : static_cast<TF1*>(nullptr);
    };
    TF1 *f_hInvMCYield_mb_nw_exp_inter = getExpInter(hInvMCYield_mc_mb_nw);
    TF1 *f_hVarMCYield_mb_nw_exp_inter = getExpInter(hVarMCYield_mc_mb_nw);

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

    TCanvas &cOneIt = *new TCanvas(Form("canvas_%s_%i", lFit_data_inv->GetName(), theRound), 
                                   Form("canvas_%s_%i", lFit_data_inv->GetName(), theRound), 
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
    TH1F &histo1DSpectra = *new TH1F("histo1DSpectra", "histo1DSpectra",1000, lMinPtPlot, lMaxPtPlot);
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
    
    histo1DSpectra.GetXaxis()->SetRangeUser(lMinPtPlot, lMaxPtPlot);
    histo1DSpectra.GetYaxis()->SetRangeUser(1E-8, 5E3);
    histo1DSpectra.DrawCopy();

    float lLegendNDC_x1 = (theRound==3) ? 0.51 : 0.61;
    auto xnew = [theLeftMargin](float xold){return theLeftMargin + (1.-theLeftMargin)*xold;};
    auto legend_pad1 = new TLegend(xnew(lLegendNDC_x1), 0.51, xnew(0.98), .96);
    legend_pad1->SetBorderSize(0);
            
    // plot for all iterations
    lHistoData_inv->SetTitle(Form("%s yields for %s", theMeson.data(), thEventCutNo.data()));
    utils_plotting::DrawAndAdd(*lHistoData_inv,"same", colorData, 1.0, legend_pad1, "Data", "lep", lLegendTextSize, true, markerStyleData, 1.0);
    
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
        printf("plotVector(): theVector.size() = %zu\n", theVector.size());
        for (TH1* ih : theVector){
            if (!ih) { 
                printf("plotVector(): WARNING: theVector contains nullptrs. Ignoring. \n");
                continue; 
            }
            printf("about to call DrawAndAdd for %s\n", ih->GetName());
            TH1& h = *ih;
            MPair lPair = getColorAndMarkerForHisto(h);           
            bool ww = lPair.second==markerStyleMCWW;
            utils_plotting::DrawAndAdd(h, "same", lPair.first, 1.0, 
                ww ? legend_pad1 : nullptr, h.GetTitle(), "lep", lLegendTextSize, true, lPair.second, 1.0);
        }
    };

    // plot MC inv yields without weights
    plotVector(vInvMCYields_wow);
    
    // and those with weights
    if (theRound){
        plotVector(vInvMCYields_ww);
    }

    // move on to fits
    if (lastItFit){
        utils_plotting::DrawAndAdd(*lastItFit, "same", kBlue, 3.0, legend_pad1, "last Fit", "l", lLegendTextSize, true);
    }

    utils_plotting::DrawAndAdd(*f_hInvMCYield_mb_nw_exp_inter, "same", colorFit+2, 2.0, legend_pad1, "MC MB NW exp inter ", "l", lLegendTextSize, true);    
    utils_plotting::DrawAndAdd(*lFit_data_inv, "same", colorFit, 3.0, legend_pad1, "Fit Data inv", "l", lLegendTextSize, true);
    utils_plotting::DrawAndAdd(lHistoData_var, "same", colorData+2, 1., legend_pad1, "lHistoData_var", "lep", lLegendTextSize, true, kCircle, 1.);
    utils_plotting::DrawAndAdd(lFit_data_var, "same", colorFit+2, 3.0, legend_pad1, "lFit_data_var", "l", lLegendTextSize, true);
    utils_plotting::DrawAndAdd(lFit_data_var_trans_inv, "same", colorFit+4, 3.0, legend_pad1, "lFit_data_var_trans_inv", "l", lLegendTextSize, true);

    // redraw data
    utils_plotting::DrawAndAdd(*lHistoData_inv, "same", colorData, 1.);
    utils_plotting::DrawAndAdd(lFit_data_var, "same", colorFit+2, 3.0);

    // done with drawing histos and functions
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
    std::vector<std::string> const lPaveTextLines({ 
        Form("%s  %s", collisionSystem.Data(), centString(thEventCutNo)),
        mesDec.data(),
        //"", // had this empty line here for optical reasons before I included the fitOptions
        Form("fitOptions: %s", theFitOption.data()), 
        Form("iteration: %d", theRound), 
        Form("%s in %s", theMeson.data(), thEventCutNo.data()),
        theEffiPlotLabel.data()});

    TPaveText &pav = utils_plotting::SetupTPaveText(
        xnew(0.045), 0.01, xnew(0.39), 0.48, lPaveTextLines, lLegendTextSize); // x1,y1,x2,y2

    //////////////////////////////////////////////////////////////////////////
    // current spectra / last iteration fit
    cout << "==================== PAD 2 ===================================\n";
    auto &pad2 = *getNextTab();

    TH1F &histo1DRatio = *new TH1F("histo1DRatio", "histo1DRatio",1000, lMinPtPlot, lMaxPtPlot);
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
    histo1DRatio.GetXaxis()->SetRangeUser(lMinPtPlot, lMaxPtPlot);
    histo1DRatio.GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio.GetYaxis()->CenterTitle(true);
    histo1DRatio.GetYaxis()->SetRangeUser(lYmin_ratio_pad2, lYmax_ratio_pad2);        
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
            if (! (entry && std::string(entry->GetLabel()).find("MC") != std::string::npos)) {
                printf("INFO: fitMesonAndWriteToFile(): skipping for i = %d\n", i);
                continue;
            }
            TObject* obj = entry->GetObject();
            if (!obj) {continue;}
            printf("Processing entry i: %d objectName: %s label: %s from legend_pad1\n", i, obj->GetName(), entry->GetLabel());

            if (!(obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("MC") == 0)) {
                printf("INFO: fitMesonAndWriteToFile(): skipping for i = %d, obj = %s\n", 
                    i, obj->GetName());
                continue;
            }
            TH1* histoMC = dynamic_cast<TH1*>(obj);
            if (!histoMC) {continue;}
            TH1 &hRatioToFit = *utils_TH1::DivideTH1ByTF1(*histoMC, *lastItFit, nullptr, nullptr, kTRUE /*theIntegrateTF1*/);
            hRatioToFit.SetName(Form("%s_mover_%s", histoMC->GetName(), lastItFit->GetName()));
            utils_plotting::DrawAndAdd(hRatioToFit, "same", histoMC->GetLineColor(), 1.0, leg2, entry->GetLabel(), "lp", lLegendTextSize);
        }
    }    
    DrawGammaLines(lMinPtPlot, lMaxPtPlot ,1., 1., 1, kBlack, 2);

    //////////////////////////////////////////////////////////////////////////
    // compare this effi to previous effis
    cout << "==================== PAD 3 ===================================\n";
    auto &pad3 = *getNextTab();

    histo1DRatio.GetYaxis()->SetTitle("this over last efficiency");
    histo1DRatio.DrawCopy();
    
    if (theRound) {
        printf("SFS theRound = %d\n", theRound);
        // draw ratio of effi/effiLast
        TH1 &lHistoTrueEffi     = *(TH1*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "TrueMesonEffiPt");
        TH1 *lHistoTrueEffiLast = theRound 
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

        DrawGammaLines(lMinPtPlot, lMaxPtPlot ,1., 1., 1, kBlack, 2);
        lEffiOverEffiLast->Draw("same");
    }

    //////////////////////////////////////////////////////////////////////////
    // show quality of fits
    cout << "==================== PAD 4 ===================================\n";

    auto &pad4 = *getNextTab();
    histo1DRatio.GetYaxis()->SetTitle("this Data over its Fit");
    histo1DRatio.GetYaxis()->SetRangeUser(lYmin_ratio,lYmax_ratio);
    histo1DRatio.DrawCopy();

    // this fit over this data
    TH1* hRatioDataToFit = utils_TH1::DivideTH1ByTF1(
        *lHistoData_inv, *lFit_data_inv, Form("hRatioDataToFit_it%d", theRound), nullptr, false/*integrateFunction*/);
    
    TH1* hRatioDataToFit_int = utils_TH1::DivideTH1ByTF1(
        lHistoData_var, lFit_data_var, Form("hRatioDataToFit_int%d", theRound), nullptr, true/*integrateFunction*/);
    
    
    auto *leg4 = utils_plotting::GetLegend(xnew(0.16),0.6,xnew(0.44),0.86);
    DrawGammaLines(lMinPtPlot, lMaxPtPlot ,1., 1., 1, kBlack, 2);

    utils_plotting::DrawAndAdd(*hRatioDataToFit, "same", colorFit, 3.0, leg4, "this data over its fit", "l", lLegendTextSize, true);
    
    // the int version
    utils_plotting::DrawAndAdd(*hRatioDataToFit_int, "same", colorFit+2, 3.0, leg4, "data_var_int", "l", lLegendTextSize, true);
    
    // the ratio int over old
    utils_plotting::DrawAndAdd(lRatioFits_int_over_old, "same", colorData, 3.0, leg4, "int_over_old", "l", lLegendTextSize, true);

    // QA PLOTS 1:
    // show quality of exp inter for MC MB nw inv yield
    // this mc_mb_nw_invYields's exponential interpolation over the histo itself
    TH1 *hQualityExpInter_mc_mb_nw_inv = (f_hInvMCYield_mb_nw_exp_inter && hInvMCYield_mc_mb_nw)
        ?   utils_computational::DivideTF1ByTH1( 
                *f_hInvMCYield_mb_nw_exp_inter, 
                *hInvMCYield_mc_mb_nw, 
                 Form("hQualityExpInter_%s_%d", 
                      f_hInvMCYield_mb_nw_exp_inter->GetName(), theRound), 
                 "hQualityExpInter_mb_nw_inv", true/*integrateFunction*/)
        :   nullptr;
    
    if (hQualityExpInter_mc_mb_nw_inv){
        utils_plotting::DrawAndAdd(*hQualityExpInter_mc_mb_nw_inv, "same", colorFit+5, 1.0, 
            leg4, "this MB_NW_inv_expInter over its histo", "lep", lLegendTextSize, true);
    }

    // QA PLOTS 2:
    // show quality of exp inter for MC MB nw var yield
    // this mc_mb_nw_varYields's exponential interpolation over the histo itself
    TH1 *hQualityExpInter_mc_mb_nw_var = (f_hVarMCYield_mb_nw_exp_inter && hVarMCYield_mc_mb_nw)
        ?   utils_computational::DivideTF1ByTH1( 
                *f_hVarMCYield_mb_nw_exp_inter, 
                *hVarMCYield_mc_mb_nw, 
                 Form("hQualityExpInter_%s_%d", 
                      f_hVarMCYield_mb_nw_exp_inter->GetName(), theRound), 
                 "hQualityExpInter_mb_nw_var", true/*integrateFunction*/)
        :   nullptr;
    
    if (hQualityExpInter_mc_mb_nw_var){
        utils_plotting::DrawAndAdd(*hQualityExpInter_mc_mb_nw_var, "same", colorFit+7, 1.0, 
            leg4, "this MB_NW_var_expInter over its histo", "lep", lLegendTextSize, true);
    }

    //////////////////////////////////////////////////////////////////////////
    // ratio of this weighted MCs over this data. They differ only as much as this over last true efficiency
    cout << "==================== PAD 5 ===================================\n";
    auto &pad5 = *getNextTab();
    pad5.SetTicks(1,2);

    histo1DRatio.GetXaxis()->SetLabelSize(lXlabelSize);
    histo1DRatio.GetXaxis()->SetTitleSize(lXtitleSize);
    histo1DRatio.GetXaxis()->SetTitleOffset(lXtitleOffset);
    
    histo1DRatio.GetYaxis()->SetRangeUser(lYmin_ratio,lYmax_ratio);        
    histo1DRatio.GetYaxis()->SetTitle("MC over Data Fit");
    histo1DRatio.DrawCopy();

    auto leg5 = new TLegend(xnew(0.144),0.73,xnew(0.44),0.92);
    leg5->SetBorderSize(0);
    leg5->Draw();

    TFile file_debug("file_debug.root", "update");
    auto &dir = *file_debug.mkdir(lUniqueTag.data());
    dir.cd();
    lHistoData_inv->Write();
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
                
                // divide by datafit
                TH1 &ratioMCtoFit = *utils_TH1::DivideTH1ByTF1(
                    *histoMC, *lFit_data_inv, nullptr, nullptr, true /*theIntegrateTF1*/);                

                utils_plotting::DrawAndAdd(
                    ratioMCtoFit, "same", histoMC->GetLineColor(), 3.0, 
                    leg5, histoMC->GetName()/*theLegString*/, "p", lLegendTextSize, true /*theLegDrawAlready*/, // when theLegString is empty,  objects name is used
                    histoMC->GetMarkerStyle(), histoMC->GetMarkerSize()); 

                dir.cd();
                histoMC->Write();
                ratioMCtoFit.Write();

            }
        }
    }
    file_debug.Close();    
    DrawGammaLines(lMinPtPlot, lMaxPtPlot ,1., 1., 1, kBlack, 2);
        
    cout << "================= saving to file and pdfs ========================\n";
    std::string lSinglesDir(theDir + "/singles");
    gSystem->mkdir(lSinglesDir.data(), true /*recursive*/);
    cOneIt.SaveAs(Form("%s/%s_%s_it%d.pdf", lSinglesDir.data(), thEventCutNo.data(), theMeson.data(), theRound));
    
    // save to file
    auto saveAllToWeightsFile = [&](bool theSaveDNDPT){
        // saves the datafit Pi0_Data_5TeV_101300_it8 without the _it8 (for instance)
        auto saveDataFit_withProperName = [](TObject const &theO){
            std::string lNewName(theO.GetName());
            size_t lPosI = lNewName.rfind("_it");
            if (lPosI != std::string::npos){
                lNewName.replace(lPosI, lNewName.size(), "");
            }
            theO.Write(lNewName.data());
            printf("Info: saveDataFit_withProperName(): saved %s as %s\n",
                theO.GetName(), lNewName.data());
        };
    
        auto saveMC_withProperName = [&](TObject const &theO){
            std::string s(theO.GetName());
            s = s.substr(3, 3);
            std::string mcId((s=="MB_")
                    ? "LHC20e3a"
                    : (s=="ASh")
                        ? "LHC20g10"
                        : "LHC24a1");
            std::string mcCutNo_suff((s=="MB_") ? "053" : "023");
            std::string lTechnicalName(
                Form("%s_%s_5TeV_%s", 
                    theMeson.data(), 
                    mcId.data(), 
                    thEventCutNo.replace(5, 3, mcCutNo_suff.data()).data()));
            theO.Write(lTechnicalName.data());
            printf("Info: saveMC_withProperName(): saved %s as %s\n", theO.GetName(), lTechnicalName.data());
            
            // LHC20e3a always needs to be saved in addition to e3b or e3c
            size_t pos = lTechnicalName.find("e3a");
            if (pos != std::string::npos){
                lTechnicalName.replace(pos, 3, isCentral ? "e3b" : "e3c");
                theO.Write(lTechnicalName.data());
                printf("Info: saveMC_withProperName(): saved %s as %s\n", theO.GetName(), lTechnicalName.data());
            }
        }; // end helper lambdas for saveAllToWeightsFile

        // saveAllToWeightsFile starts here
        TFile* hfile = new TFile(Form("%s/MCSpectraInputPbPb_Stephan_it%d%s.root", 
                                      theDir.data(), theRound, theSaveDNDPT ? "_var" : ""), 
                                 "UPDATE");
        // save data first
        TF1 &lFit = theSaveDNDPT ? lFit_data_var : lFit_data_var_trans_inv;
        lFit.SetName(lFit_data_inv->GetName()); 
        saveDataFit_withProperName(lFit);  
        
        // the MC histos
        for (TH1 const *hInv : vInvMCYields_wow){
            if (!hInv){
                printf("WARNING: fitMesonAndWriteToFile(): nullptr in vInvMCYields_wow.\n");
                continue;
            }
            TH1 const &hMC = theSaveDNDPT 
                ? *utils_computational::TranformInvariantYieldToD2DPtDY(
                    *hInv,
                    nullptr,
                    Form("%s_trans_var", 
                         hInv->GetName()),
                    Form("%s;pT (GeV/c);#frac{1}{N_ev} #frac{d^2 N}{dpT dy}", 
                         hInv->GetTitle()))
                : *hInv; 
            saveMC_withProperName(hMC);
        }

        cOneIt.Write();
        hfile->Write();
        hfile->Close();
        hfile->Delete();    
    }; // end saveAllToWeightsFile

    saveAllToWeightsFile(false /*theSaveDNDPT*/);
    saveAllToWeightsFile(true /*theSaveDNDPT*/);
    return &cOneIt;
}

//====================================================================   
void setMarginsToZero(TVirtualPad& vpad){
    vpad.SetLeftMargin(0.0);
    vpad.SetRightMargin(0.0);
    vpad.SetTopMargin(0.0);
    vpad.SetBottomMargin(0.0);
}
