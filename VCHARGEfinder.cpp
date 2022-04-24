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
#include "include/runner.h"
#include <regex>
#include <array>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace ROOT;

vector<double> findCurrent(string in) {
    //energien denne fil blev optaget ved er givet i dens titel
    string energyString = regex_replace(in, regex(R"([\D])"), "");
    //skab en pointer til root-filen der er blevet lavet af analyse.
    string matched = "match/N" + energyString + "gvm.root";
    cout << matched << endl;
    TFile *myFile = TFile::Open(matched.c_str());
    TTree *t = (TTree *) myFile->Get("a101");

    UInt_t vcharge;
    t ->SetBranchAddress("VCHARGE",&vcharge);
    UInt_t entries = t -> GetEntries();
    t -> GetEntry(0);

    UInt_t lowCharge = vcharge;
    UInt_t highCharge = vcharge;
    for(int i = 0; i < entries; i++){
        t->GetEntry(i);
        if(vcharge < lowCharge){
            lowCharge = vcharge;
        }
        if(vcharge > highCharge){
            highCharge = vcharge;
        }
    }
    double delta = highCharge-lowCharge;
    double high = highCharge;
    double low = lowCharge;
    return{delta,high,low};
}