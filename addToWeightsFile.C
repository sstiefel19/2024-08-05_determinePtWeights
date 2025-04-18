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

#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"


#include <iostream>
#include <string.h>
#include <vector>



void copyObjectFromFileToFile(std::string theFnameSource="", std::string theObjectNameSource="", std::string theFnameTarget="", std::string theObjectNameTarget=""){
    
    
    TFile lFileSrc(theFnameSource.data());
    TObject* lObjSrc = lFileSrc.Get(theObjectNameSource.data());
    if (!lObjSrc) std::cout << theObjectNameSource << " was not found\n";
    
    TFile lFileTrg(theFnameTarget.data(), "UPDATE");
    lObjSrc->Write(theObjectNameTarget != "" ? theObjectNameTarget.data() : theObjectNameSource.data());
    
    lFileSrc.Close();
    lFileTrg.Close();
}

void addToWeightsFile_LHC20g10(){
        std::string motherDir("/afterburner/2024-01-05_round1_MBptw0b_ASptw0/");

        std::string targetFileName("newUploadedFiles/MCSpectraInputPbPb_Stephan_it2.root");

        // 0-10% Pi0
        {   
            std::string cutNo("10130e03_0d200009ab770c00amd0404000_0152101500000000");
            std::string meson("Pi0");
            std::string fnameRegular(motherDir + cutNo + "/PbPb_5.02TeV/" + meson +"_data_GammaConvV1Correction_" + cutNo +".root");
            
            copyObjectFromFileToFile(
                fnameRegular,
                "MCYield_Meson_oldBinWOWeights_AddedSig",
                targetFileName,
                meson + "_LHC20g10_5TeV_10130023");
        }
        
        // 0-10% Eta
        {
            std::string cutNo("10130e03_0d200009ab770c00amd0404000_0152101500000000");
            std::string meson("Eta");
            std::string fnameRegular(motherDir + cutNo + "/PbPb_5.02TeV/" + meson +"_data_GammaConvV1Correction_" + cutNo +".root");
            
            copyObjectFromFileToFile(
                fnameRegular,
                "MCYield_Meson_oldBinWOWeights_AddedSig",
                targetFileName,
                meson + "_LHC20g10_5TeV_10130023");
        }
        
        // 30-50
        //  Pi0
        {   
            std::string cutNo("13530e03_0d200009ab770c00amd0404000_0152101500000000");
            std::string meson("Pi0");
            std::string fnameRegular(motherDir + cutNo + "/PbPb_5.02TeV/" + meson +"_data_GammaConvV1Correction_" + cutNo +".root");

            copyObjectFromFileToFile(
                fnameRegular,
                "MCYield_Meson_oldBinWOWeights_AddedSig",
                targetFileName,
                meson + "_LHC20g10_5TeV_13530023");
        }
        
        // 30-50% Eta
        {
            std::string cutNo("13530e03_0d200009ab770c00amd0404000_0152101500000000");
            std::string meson("Eta");
            std::string fnameRegular(motherDir + cutNo + "/PbPb_5.02TeV/" + meson +"_data_GammaConvV1Correction_" + cutNo +".root");
            
            copyObjectFromFileToFile(
                fnameRegular,
                "MCYield_Meson_oldBinWOWeights_AddedSig",
                targetFileName,
                meson + "_LHC20g10_5TeV_13530023");
        }
}


void addToWeightsFile_generic(std::string theSourceFileName,
                              std::string theMCname,
                              std::string theEventCutStr,
                              std::string theMeson,
                              std::string theTargetFileName){

        std::string theFullHistoName(theMeson + "_" + theMCname + "_" + "5TeV" + "_" + theEventCutStr);

        copyObjectFromFileToFile(
            theSourceFileName,
            theFullHistoName,
            theTargetFileName,
            theFullHistoName);
        
}


void addToWeightsFile(){

    std::string lSourceFileName("it2/MCSpectraInputPbPb_Stephan_it2.root");
    std::string lTargetFileName("MCSpectraInputPbPb_Stephan_it6.root");

    for (auto const &MCname : std::vector<std::string>({"LHC24a1", "LHC20g10"})){
        for (auto const &meson : std::vector<std::string>({"Pi0", "Eta"})){
            for (auto const &evtCut : std::vector<std::string>({"10130023", "13530023"})){
                addToWeightsFile_generic(lSourceFileName,
                                          MCname,
                                          evtCut,
                                          meson,
                                          lTargetFileName);
            }    
        }
    }
}