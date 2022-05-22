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

    vector<double> counts = {};
    vector<double> countErrors = {};
    vector<double> GV = {};
    vector<double> Vcharge = {};
    vector<double> solidAngle = {};
    vector<double> totalCount = {};
    /*
    for(int k = 0; k < i; k++) {
        double energyGV = stoi(regex_replace(adresses[k], regex(R"([\D])"), ""));
        createFileCoin(adresses[k], energyGV, 1.167, Ion("N15"), Ion("C12"));
        auto result = AlphacmEfitterCoin(adresses[k], 1.167, 0);
        auto current = findCurrent(adresses[k]);
        GV.push_back(energyGV);
        counts.push_back(result[4]);
        countErrors.push_back(result[5]);
        Vcharge.push_back(current[0]);
        solidAngle.push_back(result[6]);
        totalCount.push_back(result[7]);
    }*/
    createFileCoin("N879gvm.root", 879, 1.167, Ion("N15"), Ion("C12"));
    createFileCoin("N771gvm.root", 771, 1.167, Ion("N15"), Ion("C12"));
}