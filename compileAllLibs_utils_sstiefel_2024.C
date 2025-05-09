#include "TROOT.h"

void compileAllLibs_utils_sstiefel_2024()
{
    bool u1 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_utils.cpp+");
    bool u2 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_computational.cpp+");
    bool u3 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_files_strings.cpp+");
    bool u4 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/GCo.cpp+");
    bool u5 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_TF1.cpp+");
   
    bool u7 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/TH1_ExponentialInterpolation.cpp+");
    bool u8 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_TH1.cpp+");
    bool u6 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_plotting.cpp+");
    bool u9 = gROOT->ProcessLine(".L ~/analysisSoftware/utils_sstiefel_2024/src/utils_fits.cpp+");
}