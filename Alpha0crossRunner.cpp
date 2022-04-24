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
    vector<double> Vcharge {};

    for(int k = 0; k < i; k++){
        string targetStr = adresses[k].substr(0, 1);
        double energyGV = stoi(regex_replace(adresses[k], regex(R"([\D])"), ""));
        string energyString = to_string(energyGV);

        double factor = 1.16;
        if(targetStr == "N"){
            //createFileN15(adresses[j], energyGV, factor, target);
            auto result = AlphacmEfitter(adresses[k],factor);
            auto current = findCurrent(adresses[k]);
            GV.push_back(energyGV);
            counts.push_back(result[5]);
            countErrors.push_back(result[6]);
            Vcharge.push_back(current[0]);
        }
    }
}