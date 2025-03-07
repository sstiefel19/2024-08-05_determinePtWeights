#include "myStyle.h"

struct sPars {
        std::string fitFunc;
        std::string fitOpts;
        std::string tag;
    };

    typedef std::vector<sPars> tVPars;



void comparePionSpectra2();
void doMultiRound(std::string theMeson, 
                  std::string theCent,
                  std::map<std::string, tVPars const&> const &theMap,
                  std::vector<int> const &theRounds,
                  double theLeftMargin=0.25);

TCanvas* 
    fitMesonAndWriteToFile(int round, 
                           std::string eventCutNo, 
                           std::string meson, 
                           const char* fitFunction, 
                           const char* fitOption, 
                           const char* mcTag, 
                           float yMaxRatio=3.,
                           double theLeftMargin=0.0,
                           bool verticallyTight=true);

TH1D* getWeightedMCHistogramFromLastFitAndLastMCWOW(TH1* thisMCWOW, TF1* lastFit, TH1* lastMCWOW);
void squeezeAndPrepare_nSubPads(TVirtualPad& c, int n);

