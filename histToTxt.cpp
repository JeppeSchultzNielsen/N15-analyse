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
#include <include/runner.h>
#include <vector>
#include <TLorentzVector.h>
#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include <TCanvas.h>
#include <TLine.h>
#include <cmath>
#include <TH2D.h>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace ROOT;

int main(int argc, char *argv[]){
    string in = argv[1];
    string name = argv[2];

    TFile *myFile = TFile::Open(in.c_str());
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
    t->SetBranchAddress("cmE",E);
    t->SetBranchAddress("scatterAngle",scatterAngle);
    t->SetBranchAddress("mul",&mul);

    double lower_x = 70;
    double upper_x = 80;
    double lower_y = 4000;
    double upper_y = 5000;

    int steps_x = lround(upper_x - lower_x);
    int steps_y = lround(upper_y - lower_y);

    auto entries = t->GetEntries();
    auto hist =  new TH2D(name.c_str(), name.c_str(), steps_x, lower_x, upper_x, steps_y, lower_y, upper_y);
    for (Int_t i = 0; i < entries; i++) {
        //hent entry
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            auto x_data = scatterAngle[j];
            auto y_data = E[j];
            if(x_data > lower_x and x_data < upper_x and y_data > lower_y and y_data < upper_y){
                hist -> Fill(x_data,y_data);
            }
        }
    }
    string histRoot = "txtHists/"+name + ".root";
    TFile output(histRoot.c_str(), "RECREATE");
    output.cd();
    hist -> Write();
}