#include "include/runner.h"
#include "include/runner2.h"
#include <filesystem>
#include <string>
#include <iostream>
#include "dirent.h"
#include <tuple>
#include <regex>
#include <vector>
#include <fstream>
using namespace std;


vector<vector<double>> findFactors(int runs, Ion target, string adresses[], int noAdresses){
    double factor = 1.17;
    vector<double> peakPositions = {};
    vector<double> gvs = {};
    vector<double> factors = {};
    vector<double> errors = {};
    vector<double> factorerrors = {};

    string saveto = "AlphaCMs.txt";
    ofstream mytxt (saveto);
    mytxt << "GV\tAlphaCM\t\AlphaCMerr\n";

    string ionString = "N";
    if(target == Ion("Al27")){
        ionString = "A";
    }

    for(int k = 0; k < runs; k++) {
        for (int j = 0; j < noAdresses; j++) {
            string targetStr = adresses[j].substr(0, 1);
            double energyGV = stoi(regex_replace(adresses[j], regex(R"([\D])"), ""));
            string energyString = to_string(energyGV);
            if (targetStr == ionString) {
                createFileN15(adresses[j], energyGV, factor, target);
                auto result = AlphacmEfitter(adresses[j], factor);
                peakPositions.push_back(result[2]);
                gvs.push_back(energyGV);
                errors.push_back(result[3]);
                mytxt << to_string(energyGV) + "\t" + to_string(result[2]) + "\t" + to_string(result[3]) + "\n";
            }
        }
        auto result = factorFitter(gvs, peakPositions, errors, Ion(target), Ion("p"), Ion("He4"), Ion("C12"), to_string(factor));
        factor = result[0];
        factors.push_back(factor);
        factorerrors.push_back(result[1]);
    }
    mytxt.close();
    return {factors,factorerrors};
}

int main(int argc, char *argv[]) {
    //det her stykke kode er copy pastet herfra
    //https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    //det finder navnene på alle filer i match mappen.
    string adresses[100];
    int i = 0;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir("match")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            string name = ent->d_name;
            if (!(name == ".") && !(name == "..")) {
                adresses[i] = name;
                i++;
            }
        }
        closedir(dir);
    } else {
        /* could not open directory */
        perror("");
        return EXIT_FAILURE;
    }
    string saveto2 = "accFitAlpha.txt";
    ofstream mytxt2 (saveto2);
    mytxt2 << "RunNo.\tFactor\tError\n";
    auto aFactors = findFactors(2,Ion("N15"), adresses, i);
    for(int k = 0; k < aFactors[0].size(); k++){
        cout << aFactors.at(0).at(k) << "+-" << aFactors.at(1).at(k) << endl;
        mytxt2 << "alpha"+to_string(k)+"\t"+to_string(aFactors.at(0).at(k)) + "\t" + to_string(aFactors.at(1).at(k)) + "\n";
    }
    mytxt2.close();
}