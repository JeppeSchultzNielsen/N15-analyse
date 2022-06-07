//
// Created by jeppe on 3/10/22.
//
#include <string>
#include <vector>
#include <ausa/eloss/Ion.h>
#include <TLorentzVector.h>
using namespace std;
using namespace AUSA::EnergyLoss;
using namespace ROOT;

#ifndef ALUSCATTERING_RUNNER_H
#define ALUSCATTERING_RUNNER_H

void createFile(string in, double energyGV, double factor, Ion targetIon);
std::vector<double> thickness(string in, int angle);
double gaussSum(double *x, double *par);
std::vector<double> cmEfitter(string in, string factor, Ion target);
vector<double> factorFitter(vector<double> gevs, vector<double> measureds, vector<double> errors, Ion target, Ion beam, Ion ms, Ion ml, string factor);
void createFileN15(string in, double energyGV, double factor, Ion targetIon);
void createFileCutOff(string in, double energyGV, double factor, Ion targetIon);
vector<double> AlphacmEfitter(string in, double factor);
vector<double> findCurrent(string in);
TLorentzVector constructBeamVector(const Ion& beam, const Ion& targetIon,double beamEnergy);
tuple<bool,int> findPixel(UInt_t toSearch[10000][3], UInt_t FI, UInt_t BI, UInt_t id, int loopuntil);
double gauss(double *x, double *par);
void angularCross(string in);
void createFileCoin(string in, double energyGV, double factor, Ion targetIon, Ion recoilIon);
void angularScatterCross(string in);
vector<double> AlphacmEfitterCoin(string in, double factor, double gammaLoss);
void fixFile(string gvString, UInt_t startClock);
void angularAlpha1Cross(string in);
void cascadeFraction(string in);
void angularCascade(string in);
#endif //ALUSCATTERING_RUNNER_H
