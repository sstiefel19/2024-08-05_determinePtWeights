#include "myStyle.h"




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

//fit
TCanvas* 
    fitMesonAndWriteToFile(int round, 
                           std::string eventCutNo, 
                           std::string meson, 
                           const char* fitFunction, 
                           const char* fitOption, 
                           const char* mcTag, 
                           float yMaxRatio=3.){

    // ============================================================
    // 1) retrieve all filenames, histos, and information

    bool isPi0 = meson == "Pi0";
    cout << meson << " " << isPi0 << endl;

    Double_t minPtPlot  = isPi0 ?  0.3 : .9;
    Double_t maxPtPlot  = isPi0 ? 30.0 : 14.0;

    std::map<int, std::string> lMapBaseDirs{
        {0, "/2023-_analysis/afterburner/2023-10-24_newCutNewSplinesWPCWFlexCock"},
        {1, "/2023-_analysis/afterburner/2023-11-05_MBwoPtWAddedSigWptw_newDataTrain_mixedMesonAmp"},
        {2, "/2023-_analysis/afterburner/2024-08-05_allASMC_ptw0b"},
        {3, "/2023-_analysis/afterburner/2024-08-14_allASMC_ptw2"},
        {4, "/2023-_analysis/afterburner/2024-08-18_allASMC_ptw3"}
    };


    auto giveFilename = [&](int round, std::string theAfterBurnerConfig, std::string file) {
        return (round >= 0) ? Form("%s/%s/%s_0d200009ab770c00amd0404000_0152101500000000/PbPb_5.02TeV/%s_%s_%s_0d200009ab770c00amd0404000_0152101500000000.root", 
                                 lMapBaseDirs.at(round).data(), theAfterBurnerConfig.data(), eventCutNo.data(), meson.data(), file.data(), eventCutNo.data())
                            : Form("");
    };
    
    std::string filename(giveFilename(round,       "MB-AS_separate/mesons", "data_GammaConvV1Correction"));
    std::string filenameLast(giveFilename(round-1, "MB-AS_separate/mesons", "data_GammaConvV1Correction"));
    
    TH1D* lHistoData               = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "CorrectedYieldTrueEff");
    TH1D* lHistoMBonly_dataBin_WW  = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "MCYield_Meson");
    TH1D* lHistoMBonly_oldBin_WW   = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "MCYield_Meson_oldBin");
    TH1D* lHistoMBonly_oldBin_WOW  = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "MCYield_Meson_oldBinWOWeights");
    
    // addedSig 
    TH1D* lHistoMC_oldBin_WW_addedSig  = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "MCYield_Meson_oldBin_AddedSig");
    TH1D* lHistoMC_oldBin_WOW_addedSig  = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "MCYield_Meson_oldBinWOWeights_AddedSig");
    
    // get fit and histo from last iteration
    std::string fitName(Form("%s_Data_5TeV_%s0", meson.data(), eventCutNo.substr(0,5).data()));
    //~ TF1* lastItFit = (TF1*)utils_files_strings::GetObjectFromPathInFile(filenameLastWeightsFile.data(), fitName.data());
    //~ lastItFit->SetLineColor(kBlue);
    //~ auto* c0 = (TCanvas*)utils_files_strings::GetObjectFromPathInFile(filenameLastWeightsFile.data(), Form("canvas_%s", fitName.data()));
    //~ TH1D* lHistoDataLastIt = (TH1D*)c0->GetPad(1)->GetPrimitive("CorrectedYieldTrueEff");

    // 2) ======================= do this iterations fit =================================
    TF1* fitDataYield = FitObject(fitFunction, fitName.data(), meson.data(), NULL, minPtPlot, maxPtPlot);            
    TGraphAsymmErrors* graphYieldData = new TGraphAsymmErrors(lHistoData);
    graphYieldData->Fit(fitDataYield, fitOption, "", minPtPlot, maxPtPlot);
    
    // 3) ========================= plotting ===============================================
    TCanvas* c5 = new TCanvas(Form("canvas_%s_%i", fitDataYield->GetName(), round), 
                              Form("canvas_%s_%i", fitDataYield->GetName(), round), 
                              500, 1500);
    
    gStyle->SetOptTitle(0); // disables histo titles as title on subpads
    gStyle->SetOptStat(0); // 
    float lLeftMargin = 0.14;
    c5->Divide(1,4);
    
    c5->cd(1);
    gPad->SetLeftMargin(lLeftMargin);
    gPad->SetRightMargin(0.0);
    gPad->SetBottomMargin(0.001);
    gPad->SetLogx();
    gPad->SetLogy();

    // fit data and plot both together
    // ============================= pad 1 ===================================================
    float lYlableSize = 0.08;
    TH1F * histo1DSpectra;
    histo1DSpectra          = new TH1F("histo1DSpectra", "histo1DSpectra",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( 
        histo1DSpectra, 
        "#it{p}_{T} (GeV/#it{c})", 
        "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 
        0.03,  // xLableSize
        0.035, // xTitleSize
        lYlableSize,  // yLableSize
        0.055, // yTitleSize
        0.83,  // xTitleOffset
        1.3);  // yTitleOffset
    
    histo1DSpectra->GetYaxis()->SetRangeUser(1E-8, 5E3);
    histo1DSpectra->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DSpectra->GetXaxis()->SetLabelOffset(-0.01);
    histo1DSpectra->GetYaxis()->SetLabelOffset(0.01);
    histo1DSpectra->DrawCopy();
        
    lHistoData->SetTitle(Form("%s yields for %s", meson.data(), eventCutNo.data()));
    DrawGammaSetMarker(lHistoData, markerStyleData, 1.0, colorData, colorData);        // marker style, size, color, line color
    lHistoData->Draw("SAME");
        
    fitDataYield->SetLineColor(colorFit);
    fitDataYield->SetRange(minPtPlot, maxPtPlot);
    fitDataYield->Draw("SAME");  
    
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

    auto alegend = new TLegend(0.6,0.6,0.9,0.9);
    alegend->AddEntry(lHistoData, "yield from data", "lep");
    alegend->AddEntry(fitDataYield, "fit data", "l");
    //~ alegend->AddEntry(lastItFit, " last fit data", "l");
    alegend->SetBorderSize(0);
    alegend->Draw();
    
    // add text in plot
    TPaveText* pav = new TPaveText(0.18, 0.10, 0.53, 0.48, "NDC"); // x1,y1,x2,y2
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
    c5->cd(2);
    gPad->SetLeftMargin(lLeftMargin);
    gPad->SetRightMargin(0.0);
    gPad->SetBottomMargin(0.001);
    gPad->SetTopMargin(0.001);
    gPad->SetLogx();
    
    TH1F * histo1DRatio;
    histo1DRatio          = new TH1F("histo1DRatio", "histo1DRatio",1000, minPtPlot, maxPtPlot);
    SetStyleHistoTH1ForGraphs( histo1DRatio, "#it{p}_{T} (GeV/#it{c})", "Ratios to Fits", 
        0.1,  // xLableSize
        0.075, // xTitleSize
        lYlableSize,  // yLableSize
        0.055, // yTitleSize
        0.83,  // xTitleOffset
        1.);  // yTitleOffset
    
    histo1DRatio->GetXaxis()->SetLabelOffset(-0.01);
    histo1DRatio->GetYaxis()->SetLabelOffset(0.01);
    histo1DRatio->GetXaxis()->SetRangeUser(minPtPlot, maxPtPlot);
    histo1DRatio->GetYaxis()->SetRangeUser(gPad->GetUymin(), gPad->GetUymax());
    
    histo1DRatio->GetYaxis()->SetRangeUser(0.,2.5);
    histo1DRatio->DrawCopy();

    // this fit over this data
    TH1* hHistoRatioDataToFit = CalculateHistoRatioToFit(lHistoData, fitDataYield, kFALSE);
    hHistoRatioDataToFit->SetName(Form("hHistoRatioDataToFit_it%d", round));
    DrawGammaSetMarker(hHistoRatioDataToFit, 20, 1.0, colorFit, colorFit);        // marker style, size, color, line color
    hHistoRatioDataToFit->Draw("SAME");

    auto blegend = new TLegend(0.16,0.6,0.44,0.86);
    blegend->SetTextSize(0.05);
    blegend->AddEntry(hHistoRatioDataToFit,"this data over this fit","lep");
    blegend->SetBorderSize(0);
    blegend->Draw();
        
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    
    // compare this effi to previous effis
    // ==================== PAD 3 ===================================
    c5->cd(3);
    gPad->SetLeftMargin(lLeftMargin);
    gPad->SetRightMargin(0.0);
    gPad->SetBottomMargin(0.01);
    gPad->SetTopMargin(0.001);
    gPad->SetLogx();
    
    TH2F *hP3 = new TH2F("hP3", "hP3", 1, minPtPlot, maxPtPlot, 1., 0.7, 1.3);
    SetStyleHistoTH2ForGraphs( 
        hP3, 
        "", 
        "this over last efficiency", 
        0.03,  // xLableSize
        0.035, // xTitleSize
        lYlableSize,  // yLableSize
        0.055, // yTitleSize
        0.83,  // xTitleOffset
        1.);  // yTitleOffset
    hP3->Draw();

    if (round) {
        // draw ratio of effi/effiLast
        TH1 &lHistoTrueEffi     = *(TH1*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "TrueMesonEffiPt");
        TH1 *lHistoTrueEffiLast = round ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameLast.data(), "TrueMesonEffiPt")
                                        : static_cast<TH1*>(nullptr);
        
        std::pair<TH1&, TH1&> &lBinningsAligned = 
            *utils_TH1::AlignBinnings(lHistoTrueEffi, *lHistoTrueEffiLast);
        
        TH1* lEffiOverEffiLast = utils_TH1::DivideTH1ByTH1(lBinningsAligned.first, 
                                                            lBinningsAligned.second,
                                                            nullptr,
                                                            "EffiOverLastEffi");
        
        DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
        lEffiOverEffiLast->Draw("same");
    }

    // ==================== PAD 4 ===================================
    // ratio of this weighted MCs over this data. They differ only as much as this over last true efficiency
    c5->cd(4);
    gPad->SetLeftMargin(lLeftMargin);
    gPad->SetRightMargin(0.0);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.01);
    gPad->SetLogx();
    
    histo1DRatio->GetYaxis()->SetRangeUser(0., 2.5);        
    histo1DRatio->GetYaxis()->SetTitle("Ratio Weighted MCs over Data");
    
    histo1DRatio->DrawCopy();

    auto leg4 = new TLegend(0.144,0.2,0.44,0.32);
    leg4->SetBorderSize(0);
    leg4->Draw();

    // DrawGammaSetMarker(hHistoRatioAddedMCWWToFit, markerStyleMC, 1.0, colorMC, colorMC);    


    TH1* hMBoverData = utils_TH1::DivideTH1ByTH1(*lHistoMBonly_dataBin_WW, *lHistoData, 
                                                 "", Form("hMBoverData_%s_%s", eventCutNo.data(), meson.data()));
    utils_plotting::DrawAndAdd(*hMBoverData, "same", colorMBMC, 1.0, leg4, "MB", "lp", .055);
    

    // sanity check for rebin inv function
    TH1* lHistoMBonly_oldBin_WW_rebin = 
        utils_TH1::RebinDensityHistogram(*lHistoMBonly_oldBin_WW,
                                         *lHistoMBonly_dataBin_WW,
                                         "reb");
    TH1* lRatio = utils_TH1::DivideTH1ByTH1(*lHistoMBonly_oldBin_WW_rebin, *lHistoMBonly_dataBin_WW, "ratio");
    utils_TH1::PrintBinsWithErrors(*lRatio);

    if (round){
        TH1* lHistoASonly_dataBin_WW_fromReb = 
            utils_TH1::RebinDensityHistogram(*lHistoMC_oldBin_WW_addedSig,
                                            *lHistoMBonly_dataBin_WW,
                                            "reb");

        TH1* hASoverData = utils_TH1::DivideTH1ByTH1(*lHistoASonly_dataBin_WW_fromReb, *lHistoData, 
                                                    "", Form("hASoverData%s_%s", eventCutNo.data(), meson.data()));
    
        utils_plotting::DrawAndAdd(*hASoverData, "same", colorMC, 1.0, leg4, "AS", "lp", .055);
    }


    // ratio mbmcww over data
    // TH1D* hRatioMBMCWWoverData = (TH1D*)makeRatioDiffBinnings(lHistoMBonly_oldBin_WW,lHistoData, "hRatioMBMCWWoverData", "hRatioMBMCWWoverData");
    // hRatioMBMCWWoverData->Draw("SAME");
    
    DrawGammaLines(minPtPlot, maxPtPlot ,1., 1., 1, kBlack, 2);
    
    
    // ================= saving to file and pdfs ========================
    std::string eventCutNoMC(eventCutNo);
    eventCutNoMC.replace(5,2,"05");
    std::string eventCutNoMCAdd(eventCutNoMC);
    eventCutNoMCAdd.replace(6,1,"2");
    
    c5->SaveAs(Form("%s_%s.pdf", eventCutNo.data(), meson.data()));
    
    // save to file
    TFile* hfile = new TFile(Form("MCSpectraInputPbPb_Stephan_it%d.root", round),"UPDATE");

    // MB
    lHistoMBonly_oldBin_WOW->Write(Form("%s_LHC20e3a_5TeV_%s", meson.data(), eventCutNoMC.data()));
    if (eventCutNoMC=="10130053"){
        lHistoMBonly_oldBin_WOW->Write(Form("%s_LHC20e3b_5TeV_%s", meson.data(), eventCutNoMC.data()));
    }
    if (eventCutNoMC=="13530053"){
        lHistoMBonly_oldBin_WOW->Write(Form("%s_LHC20e3c_5TeV_%s", meson.data(), eventCutNoMC.data()));
    }
    
    fitDataYield->Write();
    c5->Write();
    hfile->Write();
    hfile->Close();
    hfile->Delete();

    return c5;

}

/*
    TString cent[nCentClasses]              = {"0-10%","0-20%","60-80%","0-5%","5-10%","10-20%","40-60%","20-40%","20-50%"};
    TString fitFunctionsPi0[nCentClasses] = {"rad",      "rad",  "doHag", "doHag", "doHag", "rad", "oHag", "rad",   "doHag"};
    TString fitFunctionsEta[nCentClasses] = {"modkfunc", "oHag", "oHag",  "doHag", "doHag", "rad", "oHag", "doHag", "oHag"};
    */
// fitMesonAndWriteToFile(
    //     round,
    //     "13530e03",
    //     "Pi0",
    //     "oHag",
    //     "EX0FM",
    //     "efficiency from MB + AS");

    struct sPars {
        std::string fitFunc;
        std::string fitOpts;
        std::string tag;
    };

    typedef std::vector<sPars> tVPars;

void doMultiRound(std::string theMeson, 
                  std::string theCent,
                  std::map<std::string, tVPars const&> const &theMap)
{
    std::string lNameC(Form("lCMultiRound_%s_%s", theMeson.data(), theCent.data()));

    auto &lCM = *new TCanvas(lNameC.data(), lNameC.data(), 2000, 1500);
    lCM.Divide(5,1);
    
    std::string lKey = theMeson + "_" + theCent;
    tVPars const &lVector = theMap.at(lKey);
    for (int round=0; round < 5; ++round){

        auto ci = 
            fitMesonAndWriteToFile(
                round,
                theCent,
                theMeson,
                lVector[round].fitFunc.data(),
                lVector[round].fitOpts.data(),
                lVector[round].tag.data());

        lCM.cd(round+1);
        ci->DrawClonePad();
    }
    lCM.SaveAs(Form("%s.pdf", lNameC.data()));
}

    // ===================================================
void comparePionSpectra2(){

    tVPars vPi0_101 = {{"oHag", "FM", "efficiency from MB"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"}
                      };
    
    tVPars vPi0_135 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"}
                      };

    tVPars vEta_101 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"}
                      };

    tVPars vEta_135 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"}
                      };
    
    std::map<std::string, tVPars const&> theMap;
    theMap.insert({"Pi0_10130e03", vPi0_101});
    theMap.insert({"Pi0_13530e03", vPi0_135});
    theMap.insert({"Eta_10130e03", vEta_101});
    theMap.insert({"Eta_13530e03", vEta_135});
    
    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    // doMultiRound("Pi0", "10130e03", theMap);

    for (auto meson : std::vector<std::string>{"Pi0", "Eta"}){
        for (auto evtcut : std::vector<std::string>{"10130e03", "13530e03"}){
            printf("%s %s\n", meson.data(), evtcut.data());
            doMultiRound(meson, evtcut, theMap);
        }
    }

    // int round = 4;


    // fitMesonAndWriteToFile(
    //     round,
    //     "10130e03",
    //     "Pi0",
    //     "tcmDoublePow",
    //     "FM",
    //     "efficiency from MB + AS");

    
    // fitMesonAndWriteToFile(
    //     round,
    //     "13530e03",
    //     "Pi0",
    //     "tcmDoublePow",
    //     "FM",
    //     "efficiency from MB + AS");
        
    // fitMesonAndWriteToFile(
    //     round,
    //     "10130e03",
    //     "Eta",
    //     "oHag",
    //     "EX0FM",
    //     "efficiency from MB + AS");    

    // fitMesonAndWriteToFile(
    //     round,
    //     "13530e03",
    //     "Eta",
    //     "tcmDoublePow",
    //     "FM",
    //     "efficiency from MB + AS");
}

