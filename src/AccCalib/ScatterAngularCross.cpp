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
    double_t result = par[2]*TMath::Gaus(x[0],par[0],par[1], "kTRUE") + par[5]*TMath::Gaus(x[0],par[3],par[4], "kTRUE");
    return result;
}

void angularScatterCross(string in){
    auto delta = findCurrent(in)[0];
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
    UInt_t mul;
    t->SetBranchAddress("solidAngle",solid);
    t->SetBranchAddress("id",id);
    t->SetBranchAddress("BI",BI);
    t->SetBranchAddress("FI",FI);
    t->SetBranchAddress("cmE",E); // OOOBBBBSSS SKAL VÆRE CME
    t->SetBranchAddress("scatterAngle",scatterAngle);
    t->SetBranchAddress("mul",&mul);
    auto entries = t->GetEntries();

    //skriv peak positioner til en .txt fil

    //Lav et array af histogrammer og vinkler på tilsvarende indekser
    TH1D *histograms[10000] = {};
    double_t solidAngles[10000] = {};
    double_t angles[10000] = {};
    UInt_t pixelInfo[10000][3] = {};

    auto mHe = Ion("He4").getMass();
    auto mN15 = Ion("N15").getMass();
    auto mC12 = Ion("C12").getMass();

    auto beamVector = constructBeamVector(Ion("H1"),Ion("N15"),energy);
    auto boostVector = TMath::Sqrt((energy+PROTON_MASS)*(energy+PROTON_MASS)-PROTON_MASS*PROTON_MASS)/(energy+PROTON_MASS + mN15) * TVector3(0,0,1);
    beamVector.Boost(-boostVector);
    auto cmEnergy = beamVector[3];
    double expectedE = (pow(cmEnergy,2)+pow(PROTON_MASS,2) - pow(mN15,2))/(2*cmEnergy) - PROTON_MASS;

    double maxAngle = 0;
    double minAngle = 180;
    int lastPrinted = 0;
    for (Int_t i = 0; i < entries; i++) {
        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            //hvis vi ikke har set denne vinkel før skal vi lave et nyt histogram for denne vinkel.
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
            //gider ikke events der er for lave til at være alpha0
            if(currentid == 1){
                continue;
            }
            if(!(currentE > lround(expectedE) -200) and currentE < (lround(expectedE)+199)){continue;}
            //gider ikke events, der sker i detektorkanter
            /*
            if(currentid == 0 || currentid == 1){
                if(currentFI >= 15 || currentFI <= 2 || currentBI <= 2 || currentBI >= 15){
                    continue;
                }
            }
            if(currentid ==  2) {
                if (currentBI == 1 || currentBI == 24) {
                    continue;
                }
            }*/
            auto boolAndIndex = findPixel(pixelInfo, currentFI, currentBI, currentid, lastPrinted+1);
            //case for hvis der endnu ikke findes et histogram for denne pixel.
            if (!get<0>(boolAndIndex)) {
                //skab nyt histogram til at indeholde events ved denne vinkel
                string name = "ID: " + to_string(currentid) + " FI: " + to_string(currentFI) + " BI: " +
                              to_string(currentBI) + " angle: " + to_string(currentAngle);
                histograms[lastPrinted] = new TH1D(name.c_str(), name.c_str(), 2000, 0,2000);//int(expectedE+0.5)-1000, int(expectedE+0.5)+999);
                pixelInfo[lastPrinted][0] = currentFI;
                pixelInfo[lastPrinted][1] = currentBI;
                pixelInfo[lastPrinted][2] = currentid;
                solidAngles[lastPrinted] = currentSolid;
                angles[lastPrinted] = currentAngle;

                //fyld energien ind i det skabte histogram
                histograms[lastPrinted]->Fill(currentE);
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
                histograms[get<1>(boolAndIndex)]->Fill(currentE);
            }
        }
    }
    vector<int> uniqueAngles = {};
    vector<double> uniqueSolids = {};
    TH1D *uniqueHistograms[10000];
    int k = 0;

    for(int i = int(minAngle); i < maxAngle; i=i+5){
        uniqueAngles.push_back(i);
        string name = to_string(energy) + "angle" + to_string(int(minAngle)+i*5)+"+-2.5";
        auto newhist = new TH1D(name.c_str(), name.c_str(), 2000,0,2000);//, int(expectedE+0.5)-1000, int(expectedE+0.5)+999);
        uniqueSolids.push_back(0);
        uniqueHistograms[k] = newhist;
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
                uniqueSolids[k] += solidAngles[i];
                break;
            }
            k++;
        }
        //case: der findes et histogram for denne vinkel. Læg pixelens solidangle til solidangle for denne vinkel.
        //læg histogrammerne sammen.
    }

    //for hvert vinkelhistogram: fit og gem. Læg resultater ind i en .txt fil:
    string saveto = "angScatterCross" + to_string(energy) +".txt";
    ofstream mytxt (saveto);
    mytxt << "Angle\tCounts\tCountErr\tDirectCounts\tSolidAngle\tVCharge\tCMangle\n";

    for(int i = 0; i < uniqueAngles.size(); i++) {
        string histRoot =
                "alphaAccCalib/uniqueAngleHists/" + to_string(energy) + "angle" + to_string(uniqueAngles[i]) + ".root";
        TFile output(histRoot.c_str(), "RECREATE");
        output.cd();

        auto cmEHist = uniqueHistograms[i];
        TCanvas *c1 = new TCanvas;
        TLine *l = new TLine(expectedE, 0, expectedE, cmEHist->GetMaximum());
        l->SetLineColor(kBlack);
        if (cmEHist->GetEntries() > 0) {
            TF1 *fit = new TF1("fit", gauss2, expectedE - 300, expectedE + 300, 6);
            fit->SetParameters(expectedE, 10, cmEHist->GetMaximum(),expectedE-100,10,cmEHist->GetMaximum());
            TFitResultPtr fp = cmEHist->Fit("fit", "s && Q && S", "", expectedE - 300, expectedE + 300);
            cmEHist->Draw();
            l->Draw();
            c1->Write();
            mytxt << to_string(uniqueAngles[i]) + "\t" + to_string(fp->Parameter(2)) + "\t" + to_string(fp->Error(2)) +
                     "\t"
                     + to_string(uniqueHistograms[i]->GetEntries()) + "\t" + to_string(uniqueSolids[i]) + "\t" +
                     to_string(delta)  + "\n";
        }
    }
    mytxt.close();
}