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

void txtCreator(int GV){
    string analyzed = "analyzed/N"+ to_string(GV)+"gv.root";
    TFile *myFile = TFile::Open(analyzed.c_str());
    //Hent træet
    TTree *t = (TTree*)myFile->Get("a");

    double_t E[1000];
    double_t scatterAngle[1000];
    Short_t BI[1000];
    Short_t FI[1000];
    double_t BT[1000];
    double_t FT[1000];
    Short_t canBeAlpha[1000];
    Short_t canBeCascade[1000];
    Short_t id[1000];
    double_t solid[1000];
    double_t recoilE[1000];
    double_t angDiff[1000];
    double_t pCM[1000];
    double_t timeDiff[1000];
    UInt_t mul;
    t->SetBranchAddress("FT",FT);
    t->SetBranchAddress("BT",BT);
    t->SetBranchAddress("solidAngle",solid);
    t->SetBranchAddress("angDiff",angDiff);
    t->SetBranchAddress("timeDiff",timeDiff);
    t->SetBranchAddress("pCM",pCM);
    t->SetBranchAddress("canBeAlpha",canBeAlpha);
    t->SetBranchAddress("canBeCascade",canBeCascade);
    t->SetBranchAddress("recoilE",recoilE);
    t->SetBranchAddress("id",id);
    t->SetBranchAddress("BI",BI);
    t->SetBranchAddress("FI",FI);
    t->SetBranchAddress("cmE",E); // OOOBBBBSSS SKAL VÆRE CME
    t->SetBranchAddress("scatterAngle",scatterAngle);
    t->SetBranchAddress("mul",&mul);
    auto entries = t->GetEntries();
    cout << entries << endl;

    long t0count = 0;
    long tAllcount = 0;
    //skriv events til en .txt fil
    string saveto = "txts/noFT0events"+ to_string(GV)+".txt";
    ofstream mytxt (saveto);
    for (Int_t i = 0; i < entries-10; i++) {
        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            //hvis vi ikke har set denne vinkel før skal vi lave et nyt histogram for denne vinkel.
            double currentAngle = 0;
            currentAngle += scatterAngle[j];
            short currentFI = 0;
            short currentBI = BI[j];
            short currentID = id[j];
            double currentE = 0;
            currentFI += FI[j];
            currentE += E[j];
            double currentRecoil = 0;
            currentRecoil += recoilE[j];
            double currentAngDiff = angDiff[j];
            double currentpCM = pCM[j];
            double currentTimeDiff = timeDiff[j];
            short currentCanBeAlpha = canBeAlpha[j];
            double currentFT = FT[j];
            double currentBT = BT[j];
            short currentCanBeCascade = canBeCascade[j];
            if(currentFT == 0 || currentBT == 0){
                t0count++;
            }
            tAllcount ++;
        }
    }
    mytxt << t0count << "\t" << tAllcount;
    mytxt.close();
}

int main(int argc, char *argv[]) {
    /*
    vector<short> GVs = {813,771,879,835,860,880,898,920,950,975,1000,1034,1045,1055};
    for(auto gv : GVs){
        txtCreator(gv);
    }*/
    txtCreator(879);
}