#include "comparePionSpectra2.h"


// ===============================================================
void comparePionSpectra2(){

    tVPars vPi0_101 = {{"oHag", "FM", "efficiency from MB"},                    // 0
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh"},      // 1
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 2
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 3
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 4
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 5
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 6
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"} // 7
                      };
    
    tVPars vPi0_135 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"}
                      };

    tVPars vEta_101 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"}
                      };

    tVPars vEta_135 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
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

    std::vector<int> lRounds{0, 1, 3, 5, 6};

    // doMultiRound("Pi0", "10130e03", theMap, lRounds, 0.3/*theLeftMargin*/);
    // return;
    for (auto meson : std::vector<std::string>{"Pi0", "Eta"}){
        for (auto evtcut : std::vector<std::string>{"10130e03", "13530e03"}){
            printf("%s %s\n", meson.data(), evtcut.data());
            doMultiRound(meson, evtcut, theMap, lRounds, 0.3/*theLeftMargin*/);
        }
    }
    return;


    int round = 1;

    fitMesonAndWriteToFile(
        round,
        "10130e03",
        "Pi0",
        "tcmDoublePow",
        "FM",
        "efficiency from MB + AS");
    return;
    fitMesonAndWriteToFile(
        round,
        "13530e03",
        "Pi0",
        "tcmDoublePow",
        "FM",
        "efficiency from MB + AS");
        
    fitMesonAndWriteToFile(
        round,
        "10130e03",
        "Eta",
        "oHag",
        "EX0FM",
        "efficiency from MB + AS");    

    fitMesonAndWriteToFile(
        round,
        "13530e03",
        "Eta",
        "tcmDoublePow",
        "FM",
        "efficiency from MB + AS");

    // mod
    // fitMesonAndWriteToFile(
    //     round,
    //     "13530e03",
    //     "Eta",
    //     "oHag",
    //     "EX0FM",
    //     "efficiency from MB + AS");
}

//====================================================================

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

//    
void setMarginsToZero(TVirtualPad& vpad){
    vpad.SetLeftMargin(0.0);
    vpad.SetRightMargin(0.0);
    vpad.SetTopMargin(0.0);
    vpad.SetBottomMargin(0.0);
}

void doMultiRound(std::string theMeson, 
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
    
    std::string lKey = theMeson + "_" + theCent;
    tVPars const &lVector = theMap.at(lKey);
    for (size_t iPos = 0; iPos < theRounds.size(); ++iPos){
        
        int round = theRounds[iPos];
        TCanvas* cNx1_i = 
            fitMesonAndWriteToFile(
                round,
                theCent,
                theMeson,
                lVector[round].fitFunc.data(),
                lVector[round].fitOpts.data(),
                lVector[round].tag.data(),
                3.,
                iPos ? 0. : theLeftMargin /*theLeftMargin*/);

        // prepare for verticallyTight layout
        squeezeAndPrepare_nSubPads(*cNx1_i, 4);
        // lCM.cd(iPos+1);
        lPads[iPos+1]->cd();
        cNx1_i->DrawClonePad();
        
    }

    lCM.SaveAs(Form("%s.pdf", lNameC.data()));
    lCM.SaveAs(Form("/repos_cloud/thesis_writing/figures/%s.pdf", lNameC.data()));
    
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

//====================================================================
TCanvas* 
    fitMesonAndWriteToFile(int round, 
                           std::string eventCutNo, 
                           std::string meson, 
                           const char* fitFunction, 
                           const char* fitOption, 
                           const char* mcTag, 
                           float yMaxRatio,
                           double theLeftMargin,
                           bool verticallyTight){

    // ============================================================
    // 1) retrieve all filenames, histos, and information

    bool isPi0 = meson == "Pi0";
    cout << meson << " " << isPi0 << endl;

    Double_t minPtPlot  = isPi0 ?  0.3 : .9;
    Double_t maxPtPlot  = isPi0 ? 30.0 : 14.0;

    // all iterations
    std::map<int, std::string> lMapBaseDirs{
        {0, "/2023-_analysis/afterburner/2023-10-24_newCutNewSplinesWPCWFlexCock"},
        {1, "/2023-_analysis/afterburner/2023-11-05_MBwoPtWAddedSigWptw_newDataTrain_mixedMesonAmp"},
        {2, "/2023-_analysis/afterburner/2024-08-05_allASMC_ptw0b"},
        {3, "/2023-_analysis/afterburner/2024-08-14_allASMC_ptw2"}, 
        {4, "/2023-_analysis/afterburner/2024-10-25_allASMC_ptw3_multiEffiMerge_limPt"},
        {5, "/2023-_analysis/afterburner/2024-10-31_allASMC_ptw4"},
        {6, "/2023-_analysis/afterburner/2024-11-04_allASMC_ptw5"},
        {7, "/2023-_analysis/afterburner/2024-11-07_allASMC_ptw6"}
    };

    // leave out /2023-11-05_MBwoPtWAddedSigWptw_newDataTrain_mixedMesonAmp
    // std::map<int, std::string> lMapBaseDirs{
    //     {0, "/2023-_analysis/afterburner/2023-10-24_newCutNewSplinesWPCWFlexCock"},
    //     {1, "/2023-_analysis/afterburner/2024-08-05_allASMC_ptw0b"},
    //     {2, "/2023-_analysis/afterburner/2024-08-14_allASMC_ptw2"},
    //     {3, "/2023-_analysis/afterburner/2024-10-25_allASMC_ptw3_multiEffiMerge_limPt"}
    // };

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
    std::string fitName(Form("%s_Data_5TeV_%s_it%d", meson.data(), eventCutNo.substr(0,5).data(), round));
    //~ TF1* lastItFit = (TF1*)utils_files_strings::GetObjectFromPathInFile(filenameLastWeightsFile.data(), fitName.data());
    //~ lastItFit->SetLineColor(kBlue);
    //~ auto* c0 = (TCanvas*)utils_files_strings::GetObjectFromPathInFile(filenameLastWeightsFile.data(), Form("canvas_%s", fitName.data()));
    //~ TH1D* lHistoDataLastIt = (TH1D*)c0->GetPad(1)->GetPrimitive("CorrectedYieldTrueEff");

    // 2) ======================= do this iterations fit =================================
    TF1* fitDataYield = FitObject(fitFunction, fitName.data(), meson.data(), NULL, minPtPlot, maxPtPlot);            
    TGraphAsymmErrors* graphYieldData = new TGraphAsymmErrors(lHistoData);
    graphYieldData->Fit(fitDataYield, fitOption, "", minPtPlot, maxPtPlot);
    // graphYieldData->Fit(fitDataYield, fitOption, "", 0.6, maxPtPlot);
    
    // 3) ========================= plotting ===============================================
    TCanvas* cOneIt = new TCanvas(Form("canvas_%s_%i", fitDataYield->GetName(), round), 
                              Form("canvas_%s_%i", fitDataYield->GetName(), round), 
                              300, 1000);
    
    gStyle->SetOptTitle(0); // disables histo titles as title on subpads
    gStyle->SetOptStat(0); // 
    cOneIt->SetBottomMargin(0.);
    cOneIt->SetRightMargin(0.);
    cOneIt->SetLeftMargin(0.);
    cOneIt->SetTopMargin(0.);

    cOneIt->Divide(1, 4, 0., verticallyTight ? 0. : 0.1);

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
    auto alegend = new TLegend(xnew(0.6), 0.59, xnew(0.9), 0.86);
    alegend->SetBorderSize(0);
            
    lHistoData->SetTitle(Form("%s yields for %s", meson.data(), eventCutNo.data()));
    fitDataYield->SetRange(minPtPlot, maxPtPlot);

    utils_plotting::DrawAndAdd(*lHistoData,              "same", colorData, 1.0, alegend, "Data", "lep", .055, true, markerStyleData, 1.0);
    utils_plotting::DrawAndAdd(*lHistoMBonly_oldBin_WOW, "same", colorMCMB, 1.0, alegend, "MB MC WoW", "lep", .055, true, markerStyleMCWOW, 1.0);
    if (round){
        utils_plotting::DrawAndAdd(*lHistoMC_oldBin_WOW_addedSig, "same", colorMCASh, 1.0, alegend, "ASh MC WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
        utils_plotting::DrawAndAdd(*lHistoMC_oldBin_WW_addedSig, "same", colorMCASh, 1.0, alegend,  "ASh MC WW", "lep", .055, true, markerStyleMCWW, 1.0);
    }
    utils_plotting::DrawAndAdd(*fitDataYield,            "same", colorFit, 3.0, alegend, "Fit data", "l", .055, true);
   
    
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
        TH1 &lHistoTrueEffi     = *(TH1*)utils_files_strings::GetObjectFromPathInFile(filename.data(), "TrueMesonEffiPt");
        TH1 *lHistoTrueEffiLast = round ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameLast.data(), "TrueMesonEffiPt")
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

    // DrawGammaSetMarker(hHistoRatioAddedMCWWToFit, markerStyleMC, 1.0, colorMCASh, colorMCASh);    


    TH1* hMBoverData = utils_TH1::DivideTH1ByTH1(*lHistoMBonly_dataBin_WW, *lHistoData, 
                                                 "", Form("hMBoverData_%s_%s", eventCutNo.data(), meson.data()));
    utils_plotting::DrawAndAdd(*hMBoverData, "same", colorMCMB, 1.0, leg4, round ? "MB WW" : "MB MC WoW", "lp", .055, true, markerStyleMCWOW);
    

    // sanity check for rebin inv function
    TH1* lHistoMBonly_oldBin_WW_rebin = 
        utils_TH1::RebinDensityHistogram(*lHistoMBonly_oldBin_WW,
                                         *lHistoMBonly_dataBin_WW,
                                         "reb");
    TH1* lRatio = utils_TH1::DivideTH1ByTH1(*lHistoMBonly_oldBin_WW_rebin, *lHistoMBonly_dataBin_WW, "ratio");
    utils_TH1::PrintBinsWithErrors(*lRatio);

    if (round){
        // calc ratio AS ww over data
        TH1* lHistoAShonly_dataBin_WW_fromReb = 
            utils_TH1::RebinDensityHistogram(*lHistoMC_oldBin_WW_addedSig,
                                            *lHistoMBonly_dataBin_WW,
                                            "reb");

        TH1* hAShWWoverData = utils_TH1::DivideTH1ByTH1(*lHistoAShonly_dataBin_WW_fromReb, *lHistoData, 
                                                    "", Form("hAShWWoverData%s_%s", eventCutNo.data(), meson.data()));
        
        // calc ratio AS wow over data
        TH1* lHistoAShonly_dataBin_WOW_fromReb = 
            utils_TH1::RebinDensityHistogram(*lHistoMC_oldBin_WOW_addedSig,
                                            *lHistoMBonly_dataBin_WW,
                                            "reb");
        TH1* hAShWOWoverData = utils_TH1::DivideTH1ByTH1(*lHistoAShonly_dataBin_WOW_fromReb, *lHistoData, 
                                                    "", Form("hAShWOWoverData%s_%s", eventCutNo.data(), meson.data()));

        // plot
        utils_plotting::DrawAndAdd(*hAShWOWoverData, "same", colorMCASh, 1.0, leg4, "ASh WoW", "lp", .055);
        utils_plotting::DrawAndAdd(*hAShWWoverData, "same", colorMCASh, 1.0, leg4, "ASh WW", "lp", .055);
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
    
    cOneIt->SaveAs(Form("%s_%s_it%d.pdf", eventCutNo.data(), meson.data(), round));
    
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

