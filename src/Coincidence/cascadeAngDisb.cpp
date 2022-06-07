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

void angularCascade(string in){
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
    Short_t canBeCascade[1000];
    Short_t canBeAlpha[1000];
    UInt_t mul;
    t->SetBranchAddress("solidAngle",solid);
    t->SetBranchAddress("id",id);
    t->SetBranchAddress("BI",BI);
    t->SetBranchAddress("FI",FI);
    t->SetBranchAddress("cmE",E); // OOOBBBBSSS SKAL VÆRE CME
    t->SetBranchAddress("scatterAngle",scatterAngle);
    t->SetBranchAddress("mul",&mul);
    t->SetBranchAddress("canBeCascade",canBeCascade);
    t->SetBranchAddress("canBeAlpha",canBeAlpha);
    auto entries = t->GetEntries();

    //skriv peak positioner til en .txt fil

    //Lav et array af histogrammer og vinkler på tilsvarende indekser
    TH1D *histograms[10000] = {};
    TH1D *histogramsa0[10000] = {};
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
    double expectedE = (pow(cmEnergy,2)+pow(mHe,2) - pow(mC12,2))/(2*cmEnergy) - mHe;

    double maxAngle = 0;
    double minAngle = 180;

    int lastPrinted = 0;
    for (Int_t i = 0; i < entries; i++) {
        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            //hvis vi ikke har set denne vinkel før skal vi lave et nyt histogram for denne vinkel.
            short currentCascade = canBeCascade[j];
            short currentAlpha = canBeAlpha[j];
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
            auto boolAndIndex = findPixel(pixelInfo, currentFI, currentBI, currentid, lastPrinted+1);
            //case for hvis der endnu ikke findes et histogram for denne pixel.
            if (!get<0>(boolAndIndex)) {
                //skab nyt histogram til at indeholde events ved denne vinkel
                string name = "ID: " + to_string(currentid) + " FI: " + to_string(currentFI) + " BI: " +
                              to_string(currentBI) + " angle: " + to_string(currentAngle);
                histograms[lastPrinted] = new TH1D(name.c_str(), name.c_str(), 2000, 1000,6000);//int(expectedE+0.5)-1000, int(expectedE+0.5)+999);
                string namea0 = "IDa0: " + to_string(currentid) + " FI: " + to_string(currentFI) + " BI: " +
                              to_string(currentBI) + " angle: " + to_string(currentAngle);
                histogramsa0[lastPrinted] = new TH1D(namea0.c_str(), name.c_str(), 2000, 1000,6000);
                pixelInfo[lastPrinted][0] = currentFI;
                pixelInfo[lastPrinted][1] = currentBI;
                pixelInfo[lastPrinted][2] = currentid;
                solidAngles[lastPrinted] = currentSolid;
                angles[lastPrinted] = currentAngle;
                //fyld energien ind i det skabte histogram
                if(currentCascade == 1) {
                    histograms[lastPrinted]->Fill(currentE);
                }
                if(currentAlpha == 1 and currentE > expectedE - 300 and currentE < expectedE + 300){
                    histogramsa0[lastPrinted]->Fill(currentE);
                }
                //læg en til lastPrinted, så vi er klar til næste gang der er en ny vinkel
                lastPrinted++;
                if(currentAngle > maxAngle){
                    maxAngle = currentAngle;
                }
                if(currentAngle < minAngle){
                    minAngle = currentAngle;
                }
            }
            else {
                //der fandtes allerede et histogram for dette pixel
                if(currentAlpha == 1 and currentE > expectedE - 300 and currentE < expectedE + 300){
                    histogramsa0[get<1>(boolAndIndex)]->Fill(currentE);
                }
                if(currentCascade == 1){
                    histograms[get<1>(boolAndIndex)]->Fill(currentE);
                }
            }
        }
    }
    vector<int> uniqueAngles = {};
    vector<double> uniqueSolids = {};
    vector<double> uniqueCMAngles = {};
    TH1D *uniqueHistograms[10000];
    TH1D *uniquea0Histograms[10000];
    int k = 0;

    maxAngle = 180;
    minAngle = 0;

    for(int i = int(minAngle); i < maxAngle; i=i+5){
        uniqueAngles.push_back(i);
        string name = to_string(energy) + "angle" + to_string(int(minAngle)+i*5)+"+-2.5";
        string namea0 = to_string(energy) + "a0angle" + to_string(int(minAngle)+i*5)+"+-2.5";
        auto newhist = new TH1D(name.c_str(), name.c_str(), 2000,1000,6000);//, int(expectedE+0.5)-1000, int(expectedE+0.5)+999);
        auto newhista0 = new TH1D(namea0.c_str(), name.c_str(), 2000,1000,6000);
        uniqueSolids.push_back(0);
        uniqueCMAngles.push_back(0);
        uniqueHistograms[k] = newhist;
        uniquea0Histograms[k] = newhista0;
        k++;
    }

    //lav en liste af de vinkler, vi har ramt med
    for(int i = 0; i < lastPrinted; i++){
        double angleInt = angles[i];
        int k = 0;
        //hvis den er indenfor +-1.5 af de vinkler vi allerede har fundet, så
        for(auto angle : uniqueAngles){
            //cout << angleInt << " og " << (angle-2.5) << "og" << (angle+2.5) << endl;
            if(angleInt > (angle-2.5) and angleInt < (angle+2.5)){
                uniqueHistograms[k]->Add(histograms[i]);
                uniquea0Histograms[k]->Add(histogramsa0[i]);
                uniqueSolids[k] += solidAngles[i];
                uniqueCMAngles[k] = CMangles[i];
                break;
            }
            k++;
        }
        //case: der findes et histogram for denne vinkel. Læg pixelens solidangle til solidangle for denne vinkel.
        //læg histogrammerne sammen.
    }

    //for hvert vinkelhistogram: fit og gem. Læg resultater ind i en .txt fil:
    string saveto = "AngCross/angCascadeCross" + to_string(GV) +".txt";
    ofstream mytxt (saveto);
    mytxt << "Angle\tDirectCounts\tSolidAngle\tAlpha0\tVCharge\n";

    for(int i = 0; i < uniqueAngles.size(); i++) {
        auto cmEHist = uniqueHistograms[i];
        mytxt << to_string(uniqueAngles[i]) + "\t" +
                 "\t"
                 + to_string(uniqueHistograms[i]->GetEntries()) + "\t" + to_string(uniqueSolids[i]) + "\t"
                 + to_string(uniquea0Histograms[i]->GetEntries())+ "\t" + to_string(delta)+ "\n";
    }
    mytxt.close();
}