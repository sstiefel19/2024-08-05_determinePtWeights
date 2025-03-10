#include "comparePionSpectra2.h"





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
TCanvas* 
    fitMesonAndWriteToFile(std::map<int, std::string> const &theMapBaseDirs,
                           int round, 
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

    auto giveFilename = [&](int round, std::string dataMC, std::string fileSuff) {
        std::string lCutNo(Form("%s_0d200009ab770c00amd0404000_0152101500000000", eventCutNo.data()));
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
    TH1D* lHistoMC_ASh_oldBin_WW  = round ? (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson_oldBin_AddedSig") : nullptr;
    TH1D* lHistoMC_ASh_oldBin_WOW  = round ? (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameData.data(), "MCYield_Meson_oldBinWOWeights_AddedSig") : nullptr;
            


    // alternativ get from AS file
    std::map<std::string, std::string> mMyname_ABname{
        {"ldNdPtMC_MC_WW","MC_Meson_genPt"},
        {"ldNdPtMC_MC_WOW","MC_Meson_genPt_properBinning_WOWeights"},
        {"ldNdPtMC_MC_WW_oldBin","MC_Meson_genPt_oldBin"},
        {"ldNdPtMC_MC_WOW_oldBin","MC_Meson_genPt_WOWeights"}
    };

    // AS first
    // need number of events
    TH1* lHistoNEventsASh = round ? (TH1*)utils_files_strings::GetObjectFromPathInFile(filenameAS.data(), "NEvents") : nullptr;
    float lNEventsASh = lHistoNEventsASh ? lHistoNEventsASh->GetBinContent(1) : 0.;
    std::map<std::string, TH1*> mAS_histos_inv; // no h because for round2 and round3 both mcs are merged
    if (round){
        for (auto const &p : mMyname_ABname){
            TH1D* lHisto = (TH1D*)utils_files_strings::GetObjectFromPathInFile(filenameAS.data(), p.second.data());
            mAS_histos_inv[p.first] = lHisto ? utils_computational::TranformD2DPtDYYieldToInvariantYield(*lHisto, "inv", nullptr, nullptr, 1./(lNEventsASh*1.6)) : nullptr;
        }
    }
    TH1* lHistoMC_ASh_WW = round ? mAS_histos_inv["ldNdPtMC_MC_WW"] : nullptr;
    TH1* lHistoMC_ASh_WOW = round ? mAS_histos_inv["ldNdPtMC_MC_WOW"] : nullptr;
    TH1* lHistoMC_ASh_WW_oldBin = round ? mAS_histos_inv["ldNdPtMC_MC_WW_oldBin"] : nullptr;
    TH1* lHistoMC_ASh_WOW_oldBin = round ? mAS_histos_inv["ldNdPtMC_MC_WOW_oldBin"] : nullptr;


    // AS2 second
    bool AS2 = (round >= 4);
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
        // it1 only ASh
        if (round==1) {
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WOW, "same", colorMCASh, 1.0, legend_pad1, "ASh MC WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WW,  "same", colorMCASh, 1.0, legend_pad1,  "ASh MC WW", "lep", .055, true, markerStyleMCWW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_WOW,        "same", colorOldFit, 1.0, legend_pad1, "ASh MC WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_WW,         "same", colorOldFit, 1.0, legend_pad1, "ASh MC WW", "lep", .055, true, markerStyleMCWW, 1.0);
        }
        // it2 and it3 ASl&ASh premerged
        if (round>1 && round<4){
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WOW, "same", colorMCASh, 1.0, legend_pad1, "MC (ASl&ASh) WOW", "lep", .055, true, markerStyleMCWOW, 1.0);
            utils_plotting::DrawAndAdd(*lHistoMC_ASh_oldBin_WW,  "same", colorMCASh, 1.0, legend_pad1, "MC (ASl&ASh) WW", "lep", .055, true, markerStyleMCWW, 1.0);        
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

    // TH1* hMBoverData = utils_TH1::DivideTH1ByTH1(*lHistoMBonly_dataBin_WW, *lHistoData, 
    //                                              "", Form("hMBoverData_%s_%s", eventCutNo.data(), meson.data()));
    // utils_plotting::DrawAndAdd(*hMBoverData, "same", colorMCMB, 1.0, leg4, round ? "MB WW" : "MC MB WoW", "lp", .055, true, markerStyleMCWOW);
    

    // sanity check for rebin inv function
    // TH1* lHistoMBonly_oldBin_WW_rebin = 
    //     utils_TH1::RebinDensityHistogram(*lHistoMBonly_oldBin_WW,
    //                                      *lHistoMBonly_dataBin_WW,
    //                                      "reb");
    // TH1* lRatio = utils_TH1::DivideTH1ByTH1(*lHistoMBonly_oldBin_WW_rebin, *lHistoMBonly_dataBin_WW, "ratio");
    // utils_TH1::PrintBinsWithErrors(*lRatio);

    // if (round){
    //     // calc ratio AS ww over data
    //     TH1* lHistoAShonly_dataBin_WW_fromReb = 
    //         utils_TH1::RebinDensityHistogram(*lHistoMC_ASh_oldBin_WW,
    //                                         *lHistoMBonly_dataBin_WW,
    //                                         "reb");

    //     TH1* hAShWWoverData = utils_TH1::DivideTH1ByTH1(*lHistoAShonly_dataBin_WW_fromReb, *lHistoData, 
    //                                                 "", Form("hAShWWoverData%s_%s", eventCutNo.data(), meson.data()));
        
    //     // calc ratio AS wow over data
    //     TH1* lHistoAShonly_dataBin_WOW_fromReb = 
    //         utils_TH1::RebinDensityHistogram(*lHistoMC_ASh_oldBin_WOW,
    //                                         *lHistoMBonly_dataBin_WW,
    //                                         "reb");
    //     TH1* hAShWOWoverData = utils_TH1::DivideTH1ByTH1(*lHistoAShonly_dataBin_WOW_fromReb, *lHistoData, 
    //                                                 "", Form("hAShWOWoverData%s_%s", eventCutNo.data(), meson.data()));

    //     // plot
    //     utils_plotting::DrawAndAdd(*hAShWOWoverData, "same", colorMCASh, 1.0, leg4, "ASh WoW", "lp", .055);
    //     utils_plotting::DrawAndAdd(*hAShWWoverData, "same", colorMCASh, 1.0, leg4, "ASh WW", "lp", .055);
    // }


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

