#include "myStyle.h"

struct sPars {
        std::string fitFunc;
        std::string fitOpts;
        std::string tag;
    };

    typedef std::vector<sPars> tVPars;

void doMultiRound(std::map<int, std::string> const &theMapBaseDirs,
                  std::string theMeson, 
                  std::string theCent,
                  std::map<std::string, tVPars const&> const &theMap,
                  std::vector<int> const &theRounds,
                  double theLeftMargin=0.25);

TCanvas* 
    fitMesonAndWriteToFile(std::map<int, std::string> const &theMapBaseDirs,
                           int round, 
                           std::string eventCutNo, 
                           std::string meson, 
                           const char* fitFunction, 
                           const char* fitOption, 
                           const char* mcTag, 
                           float yMaxRatio=3.,
                           double theLeftMargin=0.0,
                           bool verticallyTight=true);

TH1D* getWeightedMCHistogramFromLastFitAndLastMCWOW(TH1* thisMCWOW, TF1* lastFit, TH1* lastMCWOW);
void setMarginsToZero(TVirtualPad& vpad);
void squeezeAndPrepare_nSubPads(TVirtualPad& c, int n);

// ===============================================================
/*
    TString cent[nCentClasses]              = {"0-10%","0-20%","60-80%","0-5%","5-10%","10-20%","40-60%","20-40%","20-50%"};
    TString fitFunctionsPi0[nCentClasses] = {"rad",      "rad",  "doHag", "doHag", "doHag", "rad", "oHag", "rad",   "doHag"};
    TString fitFunctionsEta[nCentClasses] = {"modkfunc", "oHag", "oHag",  "doHag", "doHag", "rad", "oHag", "doHag", "oHag"};
    */

    /*
    next:
    todo:
        Check histo MC_Pi0InAcc_Pt in file ...ConvV1WithoutCorrectionAddSig to see if I can inv yield
        from this histo. Compare then to the one in the normal data output file.*/
void comparePionSpectra2(){

    tVPars vPi0_101 = {{"oHag", "FM", "efficiency from MB"},                    // 0
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh"},      // 1
                       {"tcmDoublePow", "FM", "efficiency from MB + (ASl&ASh)"},// 2
                       {"tcmDoublePow", "FM", "efficiency from MB + (ASl&ASh)"},// 3
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 4
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 5
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},// 6
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"} // 7
                      };
    
    tVPars vPi0_135 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + (ASl&ASh)"},
                       {"tcmDoublePow", "FM", "efficiency from MB + (ASl&ASh)"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh + ASl"}
                      };

    tVPars vEta_101 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh"},
                       {"oHag", "EX0FM", "efficiency from MB + (ASl&ASh)"},
                       {"oHag", "EX0FM", "efficiency from MB + (ASl&ASh)"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"},
                       {"oHag", "EX0FM", "efficiency from MB + ASh + ASl"}
                      };

    tVPars vEta_135 = {{"oHag", "EX0FM", "efficiency from MB"},
                       {"tcmDoublePow", "FM", "efficiency from MB + ASh"},
                       {"tcmDoublePow", "FM", "efficiency from MB + (ASl&ASh)"},
                       {"tcmDoublePow", "FM", "efficiency from MB + (ASl&ASh)"},
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
    /* std::map<int, std::string> lMapBaseDirs{
        {0, "/2023-_analysis/afterburner/2023-10-24_newCutNewSplinesWPCWFlexCock"},
        {1, "/2023-_analysis/afterburner/2024-08-05_allASMC_ptw0b"},
        {2, "/2023-_analysis/afterburner/2024-08-14_allASMC_ptw2"},
        {3, "/2023-_analysis/afterburner/2024-10-25_allASMC_ptw3_multiEffiMerge_limPt"}
    }; */

    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    // std::vector<int> lRounds{0, 1, 2, 3};
    std::vector<int> lRounds{0, 4, 5, 6, 7};

    // doMultiRound(lMapBaseDirs, "Pi0", "10130e03", theMap, lRounds, 0.3/*theLeftMargin*/);
    // return;

    for (auto meson : std::vector<std::string>{"Pi0", "Eta"}){
        for (auto evtcut : std::vector<std::string>{"10130e03", "13530e03"}){
            printf("%s %s\n", meson.data(), evtcut.data());
            doMultiRound(lMapBaseDirs, meson, evtcut, theMap, lRounds, 0.3/*theLeftMargin*/);
        }
    }
    return;


    int round = 1;

    fitMesonAndWriteToFile(
        lMapBaseDirs,
        round,
        "10130e03",
        "Pi0",
        "tcmDoublePow",
        "FM",
        "efficiency from MB + AS");
    return;
    fitMesonAndWriteToFile(
        lMapBaseDirs,
        round,
        "13530e03",
        "Pi0",
        "tcmDoublePow",
        "FM",
        "efficiency from MB + AS");
        
    fitMesonAndWriteToFile(
        lMapBaseDirs,
        round,
        "10130e03",
        "Eta",
        "oHag",
        "EX0FM",
        "efficiency from MB + AS");    

    fitMesonAndWriteToFile(
        lMapBaseDirs,
        round,
        "13530e03",
        "Eta",
        "tcmDoublePow",
        "FM",
        "efficiency from MB + AS");

    // mod
    // fitMesonAndWriteToFile(
    //    lMapBaseDirs,
    //     round,
    //     "13530e03",
    //     "Eta",
    //     "oHag",
    //     "EX0FM",
    //     "efficiency from MB + AS");
}



