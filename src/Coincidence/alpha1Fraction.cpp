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

double gauss2(double *x, double *par){
    double_t result = par[2]*TMath::Gaus(x[0],par[0],par[1], "kTRUE") + par[2] + x[0]*par[3];// + par[5]*TMath::Gaus(x[0],par[3],par[4], "kTRUE") ;
    return result;
}

void angularAlpha1Cross(string in){
    auto delta = findCurrent(in)[0];
    int GV = stoi(regex_replace(in, regex(R"([\D])"), ""));
    int energy = stoi(regex_replace(in, regex(R"([\D])"), "")) * 1.167;
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

    string name = GV + "alpha0";
    auto alpha0hist = new TH1D(name.c_str(), name.c_str(), 300, int(expectedE0+0.5)-150, int(expectedE0+0.5)+149);

    string name1 = GV + "alpha1";
    auto alpha1hist = new TH1D(name1.c_str(), name.c_str(), 200, int(expectedE1+0.5)-350, int(expectedE1+0.5)+149);

    int lastPrinted = 0;
    for (Int_t i = 0; i < entries; i++) {
        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            //hvis vi ikke har set denne vinkel før skal vi lave et nyt histogram for denne vinkel.
            short currentAlpha = 0;
            currentAlpha += canBeAlpha[j];
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
            if((currentE > lround(expectedE0) - 150) and currentE < (lround(expectedE0) + 149)){
                alpha0hist -> Fill(currentE);
            }
            if((currentE > lround(expectedE1) - 350) and currentE < (lround(expectedE1) + 149)){
                alpha1hist -> Fill(currentE);
            }
            //gider ikke events, der sker i dårlige strips
        }
    }

    //for hvert vinkelhistogram: fit og gem. Læg resultater ind i en .txt fil:
    string saveto = "AngCross/alpha1FractionwithLin" + to_string(GV) +".txt";
    ofstream mytxt (saveto);
    mytxt << "Alpha0s\tCAlpha0serror\tAlpha1counts\n";

    string histRoot =
            "alphaAccCalib/Alpha1fraction/" + to_string(GV) + ".root";
    TFile output(histRoot.c_str(), "RECREATE");
    output.cd();

    TCanvas *c1 = new TCanvas;
    TLine *l = new TLine(expectedE0, 0, expectedE0, alpha0hist->GetMaximum());
    l->SetLineColor(kBlack);
    TF1 *fit = new TF1("fit", gauss2, expectedE0 - 150, expectedE0 + 150, 3);
    fit->SetParameters(expectedE0, 10, alpha0hist->GetMaximum(),expectedE0-150, 10, alpha0hist->GetMaximum());
    fit->SetParLimits(0,expectedE0-50,expectedE0+50);
    TFitResultPtr fp = alpha0hist->Fit("fit", "s && Q && S", "", expectedE0 - 150, expectedE0 + 150);
    alpha0hist->Draw();
    l->Draw();
    c1->Write();

    TCanvas *c2 = new TCanvas;
    TLine *l2 = new TLine(expectedE1, 0, expectedE1, alpha1hist->GetMaximum());
    l2->SetLineColor(kBlack);
    TF1 *fit2 = new TF1("fit2", gauss2, expectedE1 - 50, expectedE1 + 40, 5);
    fit2->SetParameters(expectedE1, 10, alpha1hist->GetMaximum(),1,1);
    fit2->SetParLimits(0,expectedE1-50,expectedE1+50);
    TFitResultPtr fp2 = alpha1hist->Fit("fit2", "s && Q && S", "", expectedE1 - 50, expectedE1 + 40);
    alpha1hist->Draw();
    l2->Draw();
    c2->Write();

    mytxt << to_string(fp2 -> Parameter(2)) + "\t" + to_string(fp2 -> Error(2)) + "\t" + to_string(alpha0hist -> GetEntries()) + "\n";

    mytxt.close();
}