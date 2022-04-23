//
// Created by jeppe on 3/10/22.
//
#include <string>
#include <vector>
#include <ausa/eloss/Ion.h>
using namespace std;
using namespace AUSA::EnergyLoss;

#ifndef ALUSCATTERING_RUNNER_H
#define ALUSCATTERING_RUNNER_H

void createFile(string in, double energyGV, double factor, Ion targetIon);
std::vector<double> thickness(string in, int angle);
double gaussSum(double *x, double *par);
std::vector<double> findCurrent(string in);
std::vector<double> cmEfitter(string in, string factor, Ion target);
vector<double> factorFitter(vector<double> gevs, vector<double> measureds, vector<double> errors, Ion target, Ion beam, Ion ms, Ion ml, string factor);
void createFileN15(string in, double energyGV, double factor, Ion targetIon);
void createFileCutOff(string in, double energyGV, double factor, Ion targetIon);
vector<double> AlphacmEfitter(string in, double factor);
#endif //ALUSCATTERING_RUNNER_H
