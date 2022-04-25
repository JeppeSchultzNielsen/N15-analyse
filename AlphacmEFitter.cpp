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

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace ROOT;

double gauss(double *x, double *par){
    double_t result = par[2]*TMath::Gaus(x[0],par[0],par[1], "kTRUE");
    return result;
}

//løber igennem listen af pixels og returnerer true hvis pixelen er i listen
tuple<bool,int> findPixel(UInt_t toSearch[10000][3], UInt_t FI, UInt_t BI, UInt_t id, int loopuntil){
    for (int i = 0; i<loopuntil; i++){
        if(toSearch[i][0] == FI && toSearch[i][1] == BI && toSearch[i][2] == id){
            return make_tuple(true,i);
        }
    }
    return make_tuple(false,-1);
}

TLorentzVector constructBeamVector(const Ion& beam,
                                   const Ion& targetIon,
                                   double beamEnergy) {
    TLorentzVector plbeam( TVector3(0,0,sqrt(2*beamEnergy*beam.getMass())), beamEnergy+beam.getMass() );
    TLorentzVector pltarget( TVector3(0,0,0), targetIon.getMass() );
    return plbeam + pltarget;
}

vector<double> AlphacmEfitter(string in, double factor){
    //energien denne fil blev optaget ved er givet i dens titel
    string energyString = regex_replace(in, regex(R"([\D])"), "");
    double energy = stoi(energyString) * factor;
    //skab en pointer til root-filen der er blevet lavet af analyse.
    string analyzed = "analyzed/N" + energyString + "gv.root";
    TFile *myFile = TFile::Open(analyzed.c_str());
    //Hent træet
    TTree *t = (TTree*)myFile->Get("a");
    double_t cmE[100];
    double_t scatterAngle[100];
    UInt_t mul;
    Short_t id[1000];
    Short_t BI[1000];
    Short_t FI[1000];
    double_t solid[1000];
    t->SetBranchAddress("id",id);
    t->SetBranchAddress("BI",BI);
    t->SetBranchAddress("FI",FI);
    t->SetBranchAddress("cmE",cmE);
    t->SetBranchAddress("mul",&mul);
    t->SetBranchAddress("solidAngle",solid);

    auto entries = t->GetEntries();

    auto mHe = Ion("He4").getMass();
    auto mN15 = Ion("N15").getMass();
    auto mC12 = Ion("C12").getMass();

    auto beamVector = constructBeamVector(Ion("H1"),Ion("N15"),energy);
    auto boostVector = TMath::Sqrt((energy+PROTON_MASS)*(energy+PROTON_MASS)-PROTON_MASS*PROTON_MASS)/(energy+PROTON_MASS + mN15) * TVector3(0,0,1);
    beamVector.Boost(-boostVector);
    auto cmEnergy = beamVector[3];
    double expectedE = (pow(cmEnergy,2)+pow(mHe,2) - pow(mC12,2))/(2*cmEnergy) - mHe;
    auto cmEHist = new TH1D(energyString.c_str(), energyString.c_str(), 2000, int(expectedE+0.5)-1000, int(expectedE+0.5)+999);

    UInt_t pixelInfo[10000][3] = {};

    int lastPrinted = 0;

    double_t totalSolid = 0;

    for (Int_t i = 0; i < entries; i++) {
        t->GetEntry(i);
        //loop over alle hits i denne entry
        for (Int_t j = 0; j < mul; j++) {
            double_t currentEnergy = cmE[j];
            if(currentEnergy < int(expectedE+0.5)+999 && currentEnergy > (int(expectedE+0.5)-1000)) {
                cmEHist->Fill(currentEnergy);
            }
            short currentFI = 0;
            short currentBI = 0;
            short currentid = 0;
            double currentSolid = 0;
            currentFI += FI[j];
            currentBI += BI[j];
            currentid += id[j];
            currentSolid += solid[j];
            //se om denne pixel er en ny pixel
            auto boolAndIndex = findPixel(pixelInfo, currentFI, currentBI, currentid, lastPrinted+1);
            if(!get<0>(boolAndIndex)){
                pixelInfo[lastPrinted][0] = currentFI;
                pixelInfo[lastPrinted][1] = currentBI;
                pixelInfo[lastPrinted][2] = currentid;
                totalSolid += currentSolid;
                //læg en til lastPrinted, så vi er klar til næste gang der er en ny vinkel
                lastPrinted++;
            }
        }
    }
    string histRoot = "alphaAccCalib/cmEhists/"+energyString +".root";
    TFile output(histRoot.c_str(), "RECREATE");

    output.cd();
    TCanvas *c1= new TCanvas;
    TLine *l=new TLine(expectedE,0,expectedE,cmEHist -> GetMaximum());
    l->SetLineColor(kBlack);

    TF1 *fit = new TF1("fit", gauss, expectedE-20, expectedE + 100, 3);
    fit->SetParameters(expectedE,10,cmEHist->GetMaximum());
    TFitResultPtr fp = cmEHist->Fit("fit","s && Q","",expectedE-20,expectedE + 200);
    cmEHist -> Draw();
    l-> Draw();
    c1 -> Write();
    return {energy,expectedE,fp ->Parameter(0), fp -> Error(0), fp ->Parameter(2), fp -> Error(2), totalSolid, cmEHist -> GetEntries()};
}