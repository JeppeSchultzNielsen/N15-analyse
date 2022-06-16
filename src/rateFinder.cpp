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
#include <vector>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace ROOT;

void findRate(int GV) {
    string analyzed = "match/N" + to_string(GV) + "gvm.root";
    TFile *myFile = TFile::Open(analyzed.c_str());
    //Hent træet
    TTree *t = (TTree *) myFile->Get("a101");
    auto entries = t -> GetEntries();
    UInt_t mul;
    UInt_t CLOCK;
    t->SetBranchAddress("CLOCK", &CLOCK);
    t->SetBranchAddress("mul",&mul);
    UInt_t N;
    t->SetBranchAddress("___N___", &N);
    t->GetEntry(0);
    long N1 = N;
    long CLOCK1 = CLOCK;
    t->GetEntry(entries-1);
    long N2 = N;
    long CLOCK2 = CLOCK;

    string saveto = "txts/rate"+ to_string(GV)+".txt";
    ofstream mytxt (saveto);
    mytxt << N2-N1 << "\t" << CLOCK2-CLOCK1 << "\n";
    mytxt.close();
}

int main(int argc, char *argv[]){
    vector<short> GVs = {813,771,879,835,860,880,898,920,950,975,1000,1034,1045,1055};
    for(auto gv : GVs){
        findRate(gv);
    }
}