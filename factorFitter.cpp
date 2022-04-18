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
#include <TGraphErrors.h>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace ROOT;

double cmE(double *x, double *par){
    double gev = x[0];
    double a = par[0];
    double b = par[3];
    double mb = par[1];
    double mt = par[2];

    double accE = a*gev+b;
    double beta = TMath::Sqrt((pow((accE+mb),2)-pow(mb,2)))/(accE+mb+mt);
    double pz = TMath::Sqrt(pow((accE+mb),2)-pow(mb,2));
    double en = accE + mb + mt;
    double gamma = 1/TMath::Sqrt(1-pow(beta,2));
    double enCM = gamma*(en - beta*pz);
    return (pow(enCM,2)+pow(mb,2)-pow(mt,2))/(2*enCM) - mb;
}

vector<double> factorFitter(vector<double> gevs, vector<double> measureds, vector<double> errors, Ion target, Ion beam, string factor){
    int n = gevs.size();
    auto gevGraph = new TGraphErrors();
    for(int i = 0; i < n; i++){
        gevGraph -> SetPoint(i, gevs[i],measureds[i]);
        gevGraph ->SetPointError(i,errors[i],0);
    }

    TF1 *func = new TF1("fit",cmE,*std::min_element(measureds.begin(), measureds.end())-100,*std::max_element(measureds.begin(), measureds.end())+100,4);
    func -> SetParameters(1.17, beam.getMass(), target.getMass(), 20);
    func -> FixParameter(1,beam.getMass());
    func -> FixParameter(2,target.getMass());
    TFitResultPtr fp = gevGraph -> Fit(func, "s && Q && F");

    auto *canvas = new TCanvas;
    gevGraph -> Draw("*");
    string graphRoot = "aFits/"+factor+".root";
    TFile output(graphRoot.c_str(), "RECREATE");
    output.cd();
    gevGraph -> Draw();
    canvas -> Write();
    cout << "b: " << fp ->Parameter(3) << endl;
    return{fp ->Parameter(0), fp -> Error(0)};
}