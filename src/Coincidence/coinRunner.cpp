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
        if(energyGV == 771 || energyGV == 835 || energyGV == 950 || energyGV == 1710 || energyGV == 1045 || energyGV == 920 || energyGV == 860 || energyGV == 880 ||
                energyGV == 898 ||energyGV == 2564 ||energyGV == 1000 ||energyGV == 1055 || energyGV == 975 ||energyGV == 1034 || energyGV == 2460 || energyGV == 1400) continue;
        createFileCoin(adresses[k], energyGV, 1.169, Ion("N15"), Ion("C12"));
        angularCascade(adresses[k]);
    }*/
    createFileCoin("N771gvm.root", 771, 1.169, Ion("N15"), Ion("C12"));
}