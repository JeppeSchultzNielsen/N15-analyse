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

    string ionString = "N";
    if(target == Ion("Al27")){
        ionString = "A";
    }

    for(int k = 0; k < runs; k++) {
        for (int j = 0; j < noAdresses; j++) {
            string targetStr = adresses[j].substr(0, 1);
            double energyGV = stoi(regex_replace(adresses[j], regex(R"([\D])"), ""));
            string energyString = to_string(energyGV);
            if (targetStr == ionString && (energyGV == 1710 || energyGV == 2137 || energyGV == 2564) ) {
                createFile(adresses[j], energyGV, factor, target);
                auto result = cmEfitter(adresses[j], to_string(factor),target);
                peakPositions.push_back(result[2]);
                gvs.push_back(energyGV);
                errors.push_back(result[3]);
            }
        }
        auto result = factorFitter(gvs, peakPositions, errors, Ion(target), Ion("p"), to_string(factor));
        factor = result[0];
        factors.push_back(factor);
        factorerrors.push_back(result[1]);
    }
    return {factors,factorerrors};
}

int main(int argc, char *argv[]){
    //det her stykke kode er copy pastet herfra
    //https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    //det finder navnene på alle filer i match mappen.
    string adresses[100];
    int i = 0;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir ("match")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            string name = ent -> d_name;
            if(!(name == ".") && !(name == "..")){
                adresses[i] = name;
                i++;
            }
        }
        closedir (dir);
    } else {
        /* could not open directory */
        perror ("");
        return EXIT_FAILURE;
    }

    string saveto2 = "accFit.txt";
    ofstream mytxt2 (saveto2);
    mytxt2 << "RunNo.\tFactor\tError\n";

    auto nFactors = findFactors(5,Ion("N15"),adresses,i);
    auto alFactors = findFactors(5,Ion("Al27"),adresses,i);
    auto cFactors = findFactors(5,Ion("C12"),adresses,i);

    for(int k = 0; k < nFactors[0].size(); k++){
        cout << nFactors.at(0).at(k) << "+-" << nFactors.at(1).at(k) << endl;
        mytxt2 << "N"+to_string(k)+"\t"+to_string(nFactors.at(0).at(k)) + "\t" + to_string(nFactors.at(1).at(k)) + "\n";
    }
    for(int k = 0; k < alFactors[0].size(); k++){
        cout << alFactors.at(0).at(k) << "+-" << alFactors.at(1).at(k) << endl;
        mytxt2 << "Al"+to_string(k)+"\t"+to_string(alFactors.at(0).at(k)) + "\t" + to_string(alFactors.at(1).at(k)) + "\n";
    }
    for(int k = 0; k < nFactors[0].size(); k++){
        cout << cFactors.at(0).at(k) << "+-" << cFactors.at(1).at(k) << endl;
        mytxt2 << "C"+to_string(k)+"\t"+to_string(cFactors.at(0).at(k)) + "\t" + to_string(cFactors.at(1).at(k)) + "\n";
    }
    mytxt2.close();
}