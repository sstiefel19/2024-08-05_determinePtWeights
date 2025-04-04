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
                  std::map<std::string, tVPars const&> const &theMap_mesonCent_params,
                  std::vector<int> const &theRounds,
                  double theLeftMargin=0.25,
                  double theRightMargin=0.05,
                  std::string const &theDir="");

TCanvas* 
    fitMesonAndWriteToFile(std::map<int, std::string> const &theMapBaseDirs,
                           int round, 
                           std::string eventCutNo, 
                           std::string meson, 
                           const char* fitFunction, 
                           const char* fitOption, 
                           const char* effiPlotLabel, 
                           size_t thePlotWidth, // in pixels
                           double theLeftMargin=0.25,
                           double theRightMargin=0.05,
                           bool verticallyTight=true,
                           std::string theDir="");

void setMarginsToZero(TVirtualPad& vpad);

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
    
    std::map<std::string, tVPars const&> lMap_mesonCent_params;
    lMap_mesonCent_params.insert({"Pi0_10130e03", vPi0_101});
    lMap_mesonCent_params.insert({"Pi0_13530e03", vPi0_135});
    lMap_mesonCent_params.insert({"Eta_10130e03", vEta_101});
    lMap_mesonCent_params.insert({"Eta_13530e03", vEta_135});

        // all iterations
    std::map<int, std::string> lMapBaseDirs{
        {0, "~/2023-_analysis/afterburner/2023-10-24_newCutNewSplinesWPCWFlexCock"},
        {1, "~/2023-_analysis/afterburner/2023-11-05_MBwoPtWAddedSigWptw_newDataTrain_mixedMesonAmp"},
        {2, "~/2023-_analysis/afterburner/2024-08-05_allASMC_ptw0b"},
        {3, "~/2023-_analysis/afterburner/2024-08-14_allASMC_ptw2"}, 
        {4, "~/2023-_analysis/afterburner/2024-10-25_allASMC_ptw3_multiEffiMerge_limPt"},
        {5, "~/2023-_analysis/afterburner/2024-10-31_allASMC_ptw4"},
        {6, "~/2023-_analysis/afterburner/2024-11-04_allASMC_ptw5"},
        {7, "~/2023-_analysis/afterburner/2024-11-07_allASMC_ptw6"}
        // {7, "~/2025/2025-03-18_AliPhysics_testPtWeights/afterburner_ptw6"}
    };

    gROOT->Reset();   
    gROOT->SetStyle("Plain");

    std::vector<int> lRounds{0, 2, 3, 4, 5, 6, 7};
    // std::vector<int> lRounds{0, 2, 3, 4, 5, 6};
    // std::vector<int> lRounds{6, 7};

    // create meaningfull id and directory to write into
    // std::string lID(Form("%d_%d_", TDatime().GetDate(), TDatime().GetTime()));
    std::string lID(Form("%d_", TDatime().GetDate()));
    for (int i : lRounds){
        lID += Form("%d_", i);
    }
    std::string lDir(Form("pdf_root/%s", lID.data()));
    gSystem->mkdir(lDir.data(), true /*recursive*/);

    float lLeftMargin = 0.3;
    float lRightMargin = 0.02;

    // run one config
    // doMultiRound(lMapBaseDirs, "Pi0", "10130e03", lMap_mesonCent_params, lRounds,
    //              lLeftMargin/*theLeftMargin*/, lRightMargin/*theRightMargin*/, lDir);
    // doMultiRound(lMapBaseDirs, "Eta", "10130e03", lMap_mesonCent_params, lRounds,
    //              lLeftMargin/*theLeftMargin*/, lRightMargin/*theRightMargin*/, lDir);
    
    //  return;

    for (auto meson : std::vector<std::string>{"Pi0", "Eta"}){
        for (auto evtcut : std::vector<std::string>{"10130e03", "13530e03"}){
            printf("%s %s\n", meson.data(), evtcut.data());
            doMultiRound(lMapBaseDirs, 
                         meson, 
                         evtcut, 
                         lMap_mesonCent_params, 
                         lRounds, 
                         lLeftMargin/*theLeftMargin*/,
                         lRightMargin/*theRightMargin*/,
                         lDir);
        }
    }
    return;

/*
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
*/
}



