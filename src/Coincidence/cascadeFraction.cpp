#include <iostream>
#include <string>
#include <ausa/util/FileUtil.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Default.h>
#include <Math/Vector3D.h>
#include <TROOT.h>
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <fstream>
#include "TFitResult.h"
#include <regex>
#include <Math/Vector3D.h>
#include <TROOT.h>
#include <TVector3.h>
#include "include/runner.h"
#include <vector>
#include <TLorentzVector.h>
#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include <TCanvas.h>
#include <TLine.h>
#include <cmath>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace ROOT;

void cascadeFraction(string in){
    int GV = stoi(regex_replace(in, regex(R"([\D])"), ""));
    int energy = stoi(regex_replace(in, regex(R"([\D])"), "")) * 1.169;
    //skab en pointer til root-filen der er blevet lavet af analyse.
    string analyzed = "analyzed/N"+regex_replace(in, regex(R"([\D])"), "") + "gv.root";
    TFile *myFile = TFile::Open(analyzed.c_str());
    //Hent træet
    TTree *t = (TTree*)myFile->Get("a");

    double_t E[1000];
    double_t scatterAngle[1000];
    Short_t BI[1000];
    Short_t FI[1000];
    Short_t id[1000];
    double_t solid[1000];
    Short_t canBeAlpha[1000];
    Short_t canBeCascade[1000];
    double_t recoilE[1000];
    UInt_t mul;
    t->SetBranchAddress("solidAngle",solid);
    t->SetBranchAddress("id",id);
    t->SetBranchAddress("BI",BI);
    t->SetBranchAddress("FI",FI);
    t->SetBranchAddress("cmE",E); // OOOBBBBSSS SKAL VÆRE CME
    t->SetBranchAddress("scatterAngle",scatterAngle);
    t->SetBranchAddress("mul",&mul);
    t->SetBranchAddress("canBeAlpha",canBeAlpha);
    t->SetBranchAddress("canBeCascade",canBeCascade);
    t->SetBranchAddress("recoilE",recoilE);
    auto entries = t->GetEntries();

    //skriv peak positioner til en .txt fil

    //Lav et array af histogrammer og vinkler på tilsvarende indekser
    TH1D *histograms[10000] = {};
    double_t solidAngles[10000] = {};
    double_t angles[10000] = {};
    double_t CMangles[10000] = {};
    UInt_t pixelInfo[10000][3] = {};

    auto mHe = Ion("He4").getMass();
    auto mN15 = Ion("N15").getMass();
    auto mC12 = Ion("C12").getMass();

    auto beamVector = constructBeamVector(Ion("H1"),Ion("N15"),energy);
    auto boostVector = TMath::Sqrt((energy+PROTON_MASS)*(energy+PROTON_MASS)-PROTON_MASS*PROTON_MASS)/(energy+PROTON_MASS + mN15) * TVector3(0,0,1);
    beamVector.Boost(-boostVector);
    auto cmEnergy = beamVector[3];
    double expectedE0 = (pow(cmEnergy,2)+pow(mHe,2) - pow(mC12,2))/(2*cmEnergy) - mHe;
    double expectedE1 = (pow(cmEnergy,2)+pow(mHe,2) - pow(mC12+4439,2))/(2*cmEnergy) - mHe;

    long cascadeCount = 0;
    long alpha0Count = 0;

    int lastPrinted = 0;
    for (Int_t i = 0; i < entries; i++) {
        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            //hvis vi ikke har set denne vinkel før skal vi lave et nyt histogram for denne vinkel.
            short currentAlpha = 0;
            currentAlpha += canBeAlpha[j];
            short currentCascade = 0;
            currentCascade += canBeCascade[j];
            //vil kun se Alphaer
            if(currentAlpha == 0) continue;
            double currentAngle = 0;
            currentAngle += scatterAngle[j];
            short currentFI = 0;
            short currentBI = 0;
            short currentid = 0;
            double currentE = 0;
            double currentSolid = 0;
            currentFI += FI[j];
            currentBI += BI[j];
            currentid += id[j];
            currentE += E[j];
            currentSolid += solid[j];

            //skab to histogrammer:
            if((currentE > lround(expectedE0) - 150) and currentE < (lround(expectedE0) + 149)) {
                alpha0Count++;
            }
            if(currentCascade == 1){
                cascadeCount++;
            }
            //gider ikke events, der sker i dårlige strips
        }
    }

    //for hvert vinkelhistogram: fit og gem. Læg resultater ind i en .txt fil:
    string saveto = "AngCross/cascadeFraction" + to_string(GV) +".txt";
    ofstream mytxt (saveto);
    mytxt << "CascadeCount\tAlpha0counts\n";

    mytxt << to_string(cascadeCount) + "\t" + to_string(alpha0Count) + "\n";

    mytxt.close();
}