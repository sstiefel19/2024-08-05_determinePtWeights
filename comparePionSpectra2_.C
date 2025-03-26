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

#include "comparePionSpectra2_.h"
#include "TSystem.h"

void doMultiRound(std::map<int, std::string> const &theMapBaseDirs,
                  std::string theMeson, 
                  std::string theCent,
                  std::map<std::string, tVPars const&> const &theMap,
                  std::vector<int> const &theRounds,
                  double theLeftMargin=0.25)
{
    std::string lNameC(Form("lCMultiRound_%s_%s", theMeson.data(), theCent.data()));

    int nCols = theRounds.size();
    double denom = 1. + (nCols-1)*(1.-theLeftMargin);
    double w1 = 1./denom;
    double wi = (1.-theLeftMargin)/denom;

    // the master canvas
    auto &lCM = *new TCanvas(lNameC.data(), lNameC.data(), 1000, 1000);
    lCM.SetLeftMargin(0.0);
    lCM.SetRightMargin(0.0);
    lCM.SetTopMargin(0.0);
    lCM.SetBottomMargin(0.0);

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
    
    // lCM.Divide(theRounds.size(), 1, 0., 0.);

    // create meaningfull id and directory to write into
    std::string lID("");
    for (int i : theRounds){
        lID += Form("%d_", i);
    }
    std::string lDir(Form("pdf_root/%s", lID.data()));
    gSystem->mkdir(lDir.data(), true /*recursive*/);

    std::string lKey = theMeson + "_" + theCent;
    tVPars const &lVector = theMap.at(lKey);
    for (size_t iPos = 0; iPos < theRounds.size(); ++iPos){
        
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
                3.,
                iPos ? 0. : theLeftMargin /*theLeftMargin*/,
                true /*verticallyTight*/,
                lDir);

        // prepare for verticallyTight layout
        squeezeAndPrepare_nSubPads(*cNx1_i, 4);
        // lCM.cd(iPos+1);
        lPads[iPos+1]->cd();
        cNx1_i->DrawClonePad();
        
    }
    

    lCM.SaveAs(Form("%s/%s.pdf", lDir.data(), lNameC.data()));
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
                           const char* mcTag, 
                           float yMaxRatio/*=3.*/,
                           double theLeftMargin/*=0.*/,
                           bool verticallyTight/*=true*/,
                           std::string theDir/*=""*/){

    // ============================================================
    // 1) retrieve all filenames, histos, and information

    bool isPi0 = meson == "Pi0";
    cout << meson << " " << isPi0 << endl;

    Double_t minPtPlot  = isPi0 ?  0.3 : .9;
    Double_t maxPtPlot  = isPi0 ? 30.0 : 14.0;

    std::string lPhotMesCutNo("0d200009ab770c00amd0404000_0152101500000000");
    std::string lCutNo(Form("%s_%s", eventCutNo.data(), lPhotMesCutNo.data()));

    auto giveFilename = [&](int round, std::string dataMC, std::string fileSuff) {
        std::string lMCconfig(Form("MB-%s_separate", (round==2 || round==3) ? "bothASpremerged" : "AS"));
        return Form("%s/%s/mesons/%s/PbPb_5.02TeV/%s_%s_GammaConvV1Correction%s_%s.root", 
                    theMapBaseDirs.at(round).data(), lMCconfig.data(), lCutNo.data(), meson.data(), dataMC.data(), fileSuff.data(), lCutNo.data());
    };


    std::string filenameData(giveFilename(round, "data", ""));
    std::string filenameData_last(giveFilename(std::max(0, round-1), "data", ""));
    std::string filenameAS(giveFilename(round,   "MC", "HistosAddSig"));
    std::string filenameAS2(giveFilename(round,  "MC", "HistosAddSig2"));

    
    TH1D* lHistoData               = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "CorrectedYieldTrueEff");
    TH1D* lHistoMBonly_dataBin_WW  = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson");
    TH1D* lHistoMBonly_oldBin_WW   = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson_oldBin");
    TH1D* lHistoMBonly_oldBin_WOW  = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson_oldBinWOWeights");
    
    // addedSig (high pt or merged)
    TH1D* lHistoMC_ASh_oldBin_WW_fd  = round ? (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson_oldBin_AddedSig") : nullptr;
    TH1D* lHistoMC_ASh_oldBin_WOW_fd  = round ? (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson_oldBinWOWeights_AddedSig") : nullptr;
            
    // alternatively 2: get counts directly from trainfile and transform to inv yield myself
    float lMeson2GammaBR = isPi0 ? 0.98798 : 0.3931;
    std::string eventCutNoMC(eventCutNo);
    eventCutNoMC.replace(5,2,"05");
    std::string eventCutNoMCAdd(eventCutNoMC);
    eventCutNoMCAdd.replace(6,1,"2");
    bool isCentral = eventCutNo.substr(0, 3) == "101";

    // theMBAS = one of  {"mb", "as", "as2"}
    // returns nullptr for invalid round, theMBAS combinations. 
    auto computeInvariantMCYieldFromTrainFile = [&](std::string theMBAS){
        bool lIsMB = theMBAS=="mb";
        bool lIsAS2 = theMBAS=="as2";

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

        bool isPremerged = (round>1) && (round<4);
        std::string asTag(isPremerged ? "asl+ash preMerged" : lIsAS2 ? "ash" : "asl");
        std::string lNewHistoName(Form("MC_%s_WW", lIsMB ? "mb" : asTag.data()));
        std::string abac(isCentral ? "ab" : "ac");
               
        
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
            lFname = Form("GammaConvV1_%s_%s.root", lConfig.data(), isPi0 ? "pi0" : "eta");
        }

        std::string lMainDir(Form("GammaConvV1_%s/", lConfig.data()));
        std::string lEventCutNo(lIsMB ? eventCutNoMC : eventCutNoMCAdd);

        GCo lGCo(Form("%s/trains/%s/%s",
                    theMapBaseDirs.at(round).data(), lTrainSubDir.data(),lFname.data()),
                 lMainDir,
                 lEventCutNo,
                 lPhotMesCutNo);   
        
        float nEvents = ((TH1*)lGCo.GetFromESD("NEvents"))->GetBinContent(1); 
        TH1* lHistoMesonsInRap= (TH1*)lGCo.GetFromMC(Form("MC_%s_Pt", meson.data()));
        TH1* hInvYield_WW = lHistoMesonsInRap 
        ? utils_computational::TranformD2DPtDYYieldToInvariantYield(
            *utils_TH1::DivideTH1ByBinWidths(*lHistoMesonsInRap, "divW", nullptr, nullptr),
            "inv", lNewHistoName.data(), nullptr, 1./(nEvents*1.6*lMeson2GammaBR)) // 1.6 = deltaY
        : nullptr;
        
        return hInvYield_WW;
    };

    // compute all weighted inv mc yields. They should agree very well with last iterations' fits.
    std::vector<TH1*> vInvMCYields_WW;
    std::vector<std::string> vMCs({"mb"});
    if (round){vMCs.push_back("as");}
    if (round>3){vMCs.push_back("as2");}
    for (auto const &mc : vMCs){
        TH1 *hInvYield = computeInvariantMCYieldFromTrainFile(mc);
        if (hInvYield) {
            printf("Adding %s to vInvMCYields_WW\n", hInvYield->GetName()); 
            vInvMCYields_WW.push_back(hInvYield); 
        }
    }
    //=======done computing all mc inv yields from train files =========
    
    // alternatively get from AS file
    //*
    std::map<std::string, std::string> mMyname_ABname{
        {"ldNdPtMC_MC_WW","MC_Meson_genPt"},
        {"ldNdPtMC_MC_WOW","MC_Meson_genPt_properBinning_WOWeights"},
        {"ldNdPtMC_MC_WW_oldBin","MC_Meson_genPt_oldBin"}, // find this in ab
        {"ldNdPtMC_MC_WOW_oldBin","MC_Meson_genPt_WOWeights"}
    };

    // AS first
    // need number of events
    TH1* lHistoNEventsASh = round 
        ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameAS.data(), "NEvents") 
        : nullptr;
    float lNEventsASh = lHistoNEventsASh ? lHistoNEventsASh->GetBinContent(1) : 0.;
    printf("SFS lNEventsASh = %f\n", lNEventsASh);
    std::map<std::string, TH1*> mAS_histos_inv; // no h because for round2 and round3 both mcs are merged
    if (lNEventsASh){
        for (auto const &p : mMyname_ABname){
            TH1D* lHisto = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameAS.data(), p.second.data());
            mAS_histos_inv[p.first] = lHisto 
                ? utils_computational::TranformD2DPtDYYieldToInvariantYield(
                    *lHisto, "inv", nullptr, nullptr, 1./(lNEventsASh*1.6*lMeson2GammaBR)) 
                : nullptr;
        }
    }
    TH1* lHistoMC_ASh_WW =  mAS_histos_inv["ldNdPtMC_MC_WW"];
    TH1* lHistoMC_ASh_WOW =  mAS_histos_inv["ldNdPtMC_MC_WOW"];
    
    TH1* lHistoMC_ASh_WW_oldBin =  mAS_histos_inv["ldNdPtMC_MC_WW_oldBin"];
    TH1* lHistoMC_ASh_WOW_oldBin =  mAS_histos_inv["ldNdPtMC_MC_WOW_oldBin"];

    
    // AS2 second
    bool AS2 = (round >= 400);
    TH1* lHistoNEventsASl = AS2 ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameAS2.data(), "NEvents") : nullptr;
    float lNEventsASl = lHistoNEventsASl ? lHistoNEventsASl->GetBinContent(1) : 0.;
    std::map<std::string, TH1*> mASl_histos_inv;
    if (AS2){
        for (auto const &p : mMyname_ABname){
            TH1D* lHisto = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameAS2.data(), p.second.data());
            mASl_histos_inv[p.first] = lHisto ? utils_computational::TranformD2DPtDYYieldToInvariantYield(*lHisto, "inv", nullptr, nullptr, 1./(lNEventsASl*1.6)) : nullptr;
        }
    }
    TH1* lHistoMC_ASl_WW = AS2 ? mASl_histos_inv["ldNdPtMC_MC_WW"] : nullptr;
    TH1* lHistoMC_ASl_WOW = AS2 ? mASl_histos_inv["ldNdPtMC_MC_WOW"] : nullptr; 
    TH1* lHistoMC_ASl_WW_oldBin = AS2 ? mASl_histos_inv["ldNdPtMC_MC_WW_oldBin"] : nullptr;
    TH1* lHistoMC_ASl_WOW_oldBin = AS2 ? mASl_histos_inv["ldNdPtMC_MC_WOW_oldBin"] : nullptr;

    
    // get fit and histo from last iteration
    std::string fitNameInFile(Form("%s_Data_5TeV_%s0", meson.data(), eventCutNo.substr(0,5).data()));
    std::string fname_weightsFile_last(Form(
        "~/2024/2024-08-05_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it%d.root", 
        max(0, round-1)));
    printf("opening last iterations weights file %s ..\n", fname_weightsFile_last.data());
    TF1* lastItFit = (TF1*)utils_files_strings::GetObjectFromPathInFile(
        fname_weightsFile_last, 
        fitNameInFile.data());
    
    // 2) ======================= do this iterations fit =================================
    std::string fitName(Form("%s_Data_5TeV_%s_it%d", meson.data(), eventCutNo.substr(0,5).data(), round)); // need the it in the name here so we dont get many  objects with the same name when doing more than one it
    TF1* fitDataYield = FitObject(fitFunction, fitName.data(), meson.data(), NULL, minPtPlot, maxPtPlot);            
    TGraphAsymmErrors* graphYieldData = new TGraphAsymmErrors(lHistoData);
    graphYieldData->Fit(fitDataYield, fitOption, "", minPtPlot, maxPtPlot);
    // graphYieldData->Fit(fitDataYield, fitOption, "", 0.6, maxPtPlot);
    
    // 3) ========================= plotting =============================================
    TCanvas* cOneIt = new TCanvas(Form("canvas_%s_%i", fitDataYield->GetName(), round), 
                              Form("canvas_%s_%i", fitDataYield->GetName(), round), 
                              300, 1000);
    
    gStyle->SetOptTitle(0); // disables histo titles as title on subpads
    gStyle->SetOptStat(0); // 
    cOneIt->SetBottomMargin(0.);
    cOneIt->SetRightMargin(0.);
    cOneIt->SetLeftMargin(0.);
    cOneIt->SetTopMargin(0.);

    cOneIt->Divide(1, 5, 0., verticallyTight ? 0. : 0.1);

    cOneIt->cd(1);
    auto* pad1 = (TPad*)gPad;
    pad1->SetLeftMargin(theLeftMargin);
    pad1->SetLogx();
    pad1->SetLogy();

    // fit data and plot both together
    // ============================= pad 1 ===================================================
    float ratio_yTitleOffsets = 1.3;
    float lYlableSize = 0.1;
    TH1F * histo1DSpectra;
    histo1DSpectra          = new TH1F("histo1DSpectra", "histo1DSpectra",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( 
        histo1DSpectra, 
        "#it{p}_{T} (GeV/#it{c})", 
        "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
        0.03,  // xLableSize
        0.04, // xTitleSize
        lYlableSize,  // yLableSize
        0.080, // yTitleSize
        0.83,  // xTitleOffset
        1.7);  // yTitleOffset
    
    histo1DSpectra->GetYaxis()->SetRangeUser(1E-8, 5E3);
    histo1DSpectra->GetYaxis()->CenterTitle(true);
    histo1DSpectra->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DSpectra->GetXaxis()->SetLabelOffset(-0.01);
    histo1DSpectra->GetYaxis()->SetLabelOffset(0.01);
    histo1DSpectra->DrawCopy();

    auto xnew = [theLeftMargin](float xold){return theLeftMargin + (1.-theLeftMargin)*xold;};
    auto legend_pad1 = new TLegend(xnew(0.6), 0.59, xnew(0.9), 0.86);
    legend_pad1->SetBorderSize(0);
            
    lHistoData->SetTitle(Form("%s yields for %s", meson.data(), eventCutNo.data()));
    fitDataYield->SetRange(minPtPlot, maxPtPlot);
    
    // plot always
    utils_plotting::DrawAndAdd(*lHistoData,              "same", colorData, 1.0, legend_pad1, "Data", "lep", .055, true, markerStyleData, 1.0);
    utils_plotting::DrawAndAdd(*lHistoMBonly_oldBin_WOW, "same", colorMCMB, 1.0, legend_pad1, "MC MB WoW", "lep", .055, true, markerStyleMCWOW, 1.0);
    
    // not in 0th iteration
    if (round){
        // MBMC WW only from 2nd iteration onwards
        if (round>1) {
            utils_plotting::DrawAndAdd(*lHistoMBonly_oldBin_WW, "same", colorMCMB, 1.0, legend_pad1, "MC MB WW", "lep", .055, true, markerStyleMCWW, 1.0);
        }

        // plot inv mc yields ww
        for (TH1* ih : vInvMCYields_WW){
            if (!ih) { continue; }
            TH1& h = *ih;
            TString hname(h.GetName());
            Color_t lColor = hname.Contains("mb") 
                ? colorMCMB 
                : hname.Contains("+") 
                    ? kPink
                    : hname.Contains("asl")
                        ? kBlue
                        : kOrange; 
            utils_plotting::DrawAndAdd(h, "same", lColor, 1.0, legend_pad1, h.GetName(), "lep", 0.055, true, markerStyleMCWW, 1.0);
        }

        // it1 only ASh
        if (round==1) {
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WOW_fd, "same", colorMCASh, 1.0, legend_pad1, "ASh MC WOW file data", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WW_fd,  "same", colorMCASh, 1.0, legend_pad1,  "ASh MC WW file data", "lep", .055, true, markerStyleMCWW, 1.0);
            // utils_plotting::DrawAndAdd(*lHistoMC_ASh_WOW,        "same", colorOldFit, 1.0, legend_pad1, "ASh MC WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            // utils_plotting::DrawAndAdd(*lHistoMC_ASh_WW,         "same", colorOldFit, 1.0, legend_pad1, "ASh MC WW", "lep", .055, true, markerStyleMCWW, 1.0);
            // utils_plotting::DrawAndAdd(*lHistoMC_ASh_WOW_oldBin, "same", colorOldFit, 1.0, legend_pad1, "ASh MC WOW fileAS", "lep", .055, true, markerStyleMCWOW, 1.0);
            // utils_plotting::DrawAndAdd(*lHistoMC_ASh_WW_oldBin,  "same", colorOldFit, 1.0, legend_pad1, "ASh MC WW fileAS", "lep", .055, true, markerStyleMCWW, 1.0);
            
            // utils_plotting::DrawAndAdd(*hInvYieldMB_WW, "same", kBlue, 1.0, legend_pad1, "MB MC WW own calc", "lep", .055, true, markerStyleMCWW, 1.0);
            // utils_plotting::DrawAndAdd(*hInvYieldAS_WW, "same", kOrange, 1.0, legend_pad1, "ASh MC WW own calc", "lep", .055, true, markerStyleMCWW, 1.0);
            // utils_plotting::DrawAndAdd(*hInvYieldAS2_WW, "same", kMagenta, 1.0, legend_pad1, "ASl MC WW own calc", "lep", .055, true, markerStyleMCWW, 1.0);
            
        }
        // it2 and it3 ASl&ASh premerged
        if (round>1 && round<4){
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WOW_fd, "same", colorMCASh, 1.0, legend_pad1, "MC (ASl&ASh) WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WW_fd,  "same", colorMCASh, 1.0, legend_pad1, "MC (ASl&ASh) WW", "lep", .055, true, markerStyleMCWW, 1.0);        
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_WOW,        "same", colorOldFit, 1.0, legend_pad1, "MC (ASl&ASh) WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_WW,         "same", colorOldFit, 1.0, legend_pad1, "MC (ASl&ASh) WW", "lep", .055, true, markerStyleMCWW, 1.0);
        }

        if (AS2){
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_WOW,        "same", colorOldFit, 1.0, legend_pad1, "ASh MC WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_WW,         "same", colorOldFit, 1.0, legend_pad1, "ASh MC WW", "lep", .055, true, markerStyleMCWW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASl_WOW_oldBin, "same", colorMCASl, 1.0, legend_pad1, "MC ASl WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASl_WW_oldBin,  "same", colorMCASl, 1.0, legend_pad1, "MC ASl WW", "lep", .055, true, markerStyleMCWW, 1.0);

        }
        
    }
    if (lastItFit){
        utils_plotting::DrawAndAdd(*lastItFit, "same", kBlue, 3.0, legend_pad1, "last Fit data", "l", .055, true);
    }
    utils_plotting::DrawAndAdd(*fitDataYield, "same", colorFit, 3.0, legend_pad1, "Fit data", "l", .055, true);
    
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
    TPaveText* pav = new TPaveText(xnew(0.045), 0.10, xnew(0.39), 0.48, "NDC"); // x1,y1,x2,y2
    pav->AddText(Form("%s  %s", collisionSystem.Data(), centString(eventCutNo)));
    pav->AddText(mesDec.data());
    pav->AddText("");
    
    pav->AddText(Form("iteration: %d", round));
    pav->AddText(Form("%s in %s", meson.data(), eventCutNo.data()));
    pav->AddText(Form("fit function: %s", fitFunction));
    pav->AddText(mcTag);
    //~ ((TText*)pav->GetListOfLines()->Last())->SetTextAlign(22); // centered hor. and vert.
    pav->SetTextSize(0.055);
    pav->SetBorderSize(0);
    pav->SetFillStyle(1001);
    pav->SetFillColor(kWhite);
    pav->SetTextAlign(11);
    pav->Draw();
        
    // show quality of fit
    // ==================== PAD 2 ===================================
    cOneIt->cd(2);
    auto* pad2 = (TPad*)gPad;
    pad2->SetLeftMargin(theLeftMargin);
    pad2->SetLogx();
    
    TH1F * histo1DRatio;
    histo1DRatio          = new TH1F("histo1DRatio", "histo1DRatio",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio, "#it{p}_{T} (GeV/#it{c})", "Ratios to Fits", 
        0.1,  // xLableSize
        0.075, // xTitleSize
        lYlableSize,  // yLableSize
        0.075, // yTitleSize
        1.4,  // xTitleOffset
        ratio_yTitleOffsets);  // yTitleOffset
    
    histo1DRatio->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio->GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio->GetYaxis()->CenterTitle(true);
    histo1DRatio->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio->GetYaxis()->SetRangeUser(gPad->GetUymin(), gPad->GetUymax());
    
    histo1DRatio->GetYaxis()->SetRangeUser(0.,2.5);
    histo1DRatio->DrawCopy();

    // this fit over this data
    TH1* hHistoRatioDataToFit = CalculateHistoRatioToFit(lHistoData, fitDataYield, kFALSE);
    hHistoRatioDataToFit->SetName(Form("hHistoRatioDataToFit_it%d", round));
    DrawGammaSetMarker(hHistoRatioDataToFit, 20, 1.0, colorFit, colorFit);        // marker style, size, color, line color
    hHistoRatioDataToFit->Draw("SAME");

    auto blegend = new TLegend(xnew(0.16),0.6,xnew(0.44),0.86);
    blegend->SetTextSize(0.05);
    blegend->AddEntry(hHistoRatioDataToFit,"this data over this fit","lep");
    blegend->SetBorderSize(0);
    blegend->Draw();
        
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    
    // compare this effi to previous effis
    // ==================== PAD 3 ===================================
    cOneIt->cd(3);
    auto* pad3 = (TPad*)gPad;   
    pad3->SetLeftMargin(theLeftMargin);
    pad3->SetLogx();
    
    TH2F *hP3 = new TH2F("hP3", "hP3", 1, minPtPlot, maxPtPlot, 1., 0.7, 1.3);
    SetStyleHistoTH2ForGraphs( 
        hP3, 
        "", 
        "this over last efficiency", 
        0.03,  // xLableSize
        0.035, // xTitleSize
        lYlableSize,  // yLableSize
        0.075, // yTitleSize
        0.9,  // xTitleOffset
        ratio_yTitleOffsets);  // yTitleOffset
    hP3->GetYaxis()->CenterTitle(true);

    hP3->Draw();

    if (round) {
        // draw ratio of effi/effiLast
        TH1 &lHistoTrueEffi     = *(TH1*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "TrueMesonEffiPt");
        TH1 *lHistoTrueEffiLast = round ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameData_last.data(), "TrueMesonEffiPt")
                                        : static_cast<TH1*>(nullptr);

        utils_TH1::PrintBinsWithErrors(lHistoTrueEffi);
        
        std::pair<TH1&, TH1&> &lBinningsAligned = 
            *utils_TH1::AlignBinnings(lHistoTrueEffi, *lHistoTrueEffiLast);
        
        TH1* lEffiOverEffiLast = utils_TH1::DivideTH1ByTH1(lBinningsAligned.first, 
                                                            lBinningsAligned.second,
                                                            nullptr,
                                                            "EffiOverLastEffi");
        
        utils_TH1::PrintBinsWithErrors(*lEffiOverEffiLast);

        
        DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
        lEffiOverEffiLast->Draw("same");
    }

    // ==================== PAD 4 ===================================
    // ratio of this weighted MCs over this data. They differ only as much as this over last true efficiency
    cOneIt->cd(4);
    auto* pad4 = (TPad*)gPad;
    pad4->SetLeftMargin(theLeftMargin);
    pad4->SetLogx();
    pad4->SetBottomMargin(0.25);
    
    histo1DRatio->GetYaxis()->SetRangeUser(0., 2.5);        
    histo1DRatio->GetYaxis()->SetTitle("MC over Data");
    histo1DRatio->DrawCopy();

    auto leg4 = new TLegend(xnew(0.144),0.73,xnew(0.44),0.92);
    leg4->SetBorderSize(0);
    leg4->Draw();

    // get all MC histos from pad1 and calculate ratio to data
    TList* entries = legend_pad1->GetListOfPrimitives();
    for (int i = 0; i < entries->GetEntries(); ++i) {
        TLegendEntry* entry = (TLegendEntry*)entries->At(i);
        // entry->Dump();
        if (entry && std::string(entry->GetLabel()).find("MC") != std::string::npos) {
            TObject* obj = entry->GetObject();
            if (!obj) {continue;}
            printf("Processing entry i: %d objectName: %s label: %s from legend_pad1\n", i, obj->GetName(), entry->GetLabel());
            if (obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("MC") == 0) {
                TH1* histoMC = dynamic_cast<TH1*>(obj);
                if (!histoMC) {continue;}
                TH1& histoMC_dataBin = (histoMC->GetNbinsX() == lHistoData->GetNbinsX()) 
                    ? *histoMC 
                    : *utils_TH1::RebinDensityHistogram(*histoMC, *lHistoData, "reb"); 
                // if (histoMC->GetNbinsX() != lHistoData->GetNbinsX()) {
                //     printf("Rebinning %s for ratio with data\n", histoMC->GetName());
                //     histoMC_dataBin = *utils_TH1::RebinDensityHistogram(*histoMC, *lHistoData, "reb");
                // }
                TH1* ratioHist = utils_TH1::DivideTH1ByTH1(histoMC_dataBin, *lHistoData, "", Form("over_%s", lHistoData->GetName()));
                utils_plotting::DrawAndAdd(*ratioHist, "same", histoMC->GetLineColor(), 1.0, leg4, entry->GetLabel(), "lp", .055);
            }
        }
    }    
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    

    // ==================== PAD 5 ===================================
    cout << "==================== PAD 5 ===================================\n";
    // this weighted MC / last iteration fit
    cOneIt->cd(5);
    auto* pad5 = (TPad*)gPad;
    pad5->SetLeftMargin(theLeftMargin);
    pad5->SetLogx();
    pad5->SetBottomMargin(0.25);
    
    histo1DRatio->GetYaxis()->SetRangeUser(0.9, 1.1);        
    histo1DRatio->GetYaxis()->SetTitle("This MC over last fit");
    histo1DRatio->DrawCopy();

    auto leg5 = new TLegend(xnew(0.144),0.73,xnew(0.44),0.92);
    leg5->SetBorderSize(0);
    leg5->Draw();

    if (lastItFit)
    {
        // get all MC histos from pad1 and calculate ratio to last it's fit (the fit that was used for weighting the MCs)
        TList* entries = legend_pad1->GetListOfPrimitives();
        for (int i = 0; i < entries->GetEntries(); ++i) {
            TLegendEntry* entry = (TLegendEntry*)entries->At(i);
            // entry->Dump();
            if (entry && std::string(entry->GetLabel()).find("MC") != std::string::npos) {
                TObject* obj = entry->GetObject();
                if (!obj) {continue;}
                printf("Processing entry i: %d objectName: %s label: %s from legend_pad1\n", i, obj->GetName(), entry->GetLabel());
                if (obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("MC") == 0) {
                    TH1* histoMC = dynamic_cast<TH1*>(obj);
                    if (!histoMC) {continue;}
                    
                    TH1* hHistoRatioToFit = CalculateHistoRatioToFit(histoMC, lastItFit, kTRUE);
                    hHistoRatioToFit->SetName(Form("%s_mover_%s", histoMC->GetName(), lastItFit->GetName()));

                    utils_plotting::DrawAndAdd(*hHistoRatioToFit, "same", histoMC->GetLineColor(), 1.0, leg5, entry->GetLabel(), "lp", .055);
                }
            }
        }
    }    
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);

    
    // ================= saving to file and pdfs ========================
    std::string lSinglesDir(theDir + "/singles");
    gSystem->mkdir(lSinglesDir.data(), true /*recursive*/);
    cOneIt->SaveAs(Form("%s/%s_%s_it%d.pdf", lSinglesDir.data(), eventCutNo.data(), meson.data(), round));
    
    // save to file
    TFile* hfile = new TFile(Form("%s/MCSpectraInputPbPb_Stephan_it%d.root", theDir.data(), round),"UPDATE");

    // MB
    lHistoMBonly_oldBin_WOW->Write(Form("%s_LHC20e3a_5TeV_%s", meson.data(), eventCutNoMC.data()));
    if (eventCutNoMC=="10130053"){
        lHistoMBonly_oldBin_WOW->Write(Form("%s_LHC20e3b_5TeV_%s", meson.data(), eventCutNoMC.data()));
    }
    if (eventCutNoMC=="13530053"){
        lHistoMBonly_oldBin_WOW->Write(Form("%s_LHC20e3c_5TeV_%s", meson.data(), eventCutNoMC.data()));
    }
    
    fitDataYield->Write();
    cOneIt->Write();
    hfile->Write();
    hfile->Close();
    hfile->Delete();

    return cOneIt;

}

//============================================================================
TH1D* getWeightedMCHistogramFromLastFitAndLastMCWOW(TH1* thisMCWOW, TF1* lastFit, TH1* lastMCWOW){
    
    TH1D* thisMCWWcalculated = (TH1D*)thisMCWOW->Clone("thisMCWWcalculated"); 
    thisMCWWcalculated->Reset("ICES");
    for (int iBin=1; iBin<=thisMCWWcalculated->GetNbinsX(); ++iBin){
        float pt = thisMCWWcalculated->GetBinCenter(iBin);
        float ptEdge = thisMCWWcalculated->GetBinLowEdge(iBin);
        float yieldlastMC = lastMCWOW->Interpolate(pt);
        
        //~ cout << ptEdge << ": " << lastFit->Eval(ptEdge)/lastMCWOW->Interpolate(ptEdge) << endl;
        
        float weight = yieldlastMC ? lastFit->Eval(pt)/yieldlastMC : 1.;
        thisMCWWcalculated->SetBinContent(iBin, thisMCWOW->GetBinContent(iBin) * weight );
        thisMCWWcalculated->SetBinError(iBin, thisMCWOW->GetBinError(iBin) * weight );
        
        //~ cout << pt << ": " << weight << endl;
        
    }
    return thisMCWWcalculated;
    
}

//====================================================================   
void setMarginsToZero(TVirtualPad& vpad){
    vpad.SetLeftMargin(0.0);
    vpad.SetRightMargin(0.0);
    vpad.SetTopMargin(0.0);
    vpad.SetBottomMargin(0.0);
}

//====================================================================
void squeezeAndPrepare_nSubPads(TVirtualPad& vpad, int n){
    
    for (int i=1; i<=n; ++i){
        vpad.cd(i);
        auto* p = (TPad*)gPad;        
        p->SetTicks(1,1); // enable ticks on both x and y axes (all sides)
        p->SetRightMargin(0.0);
        p->SetTopMargin(i==1 ? 0.1 : 0.0);
        p->SetBottomMargin(i==4 ? 0.25 : 0.0);    
    }
    vpad.Update();
}

