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

    for(int k = 0; k < i; k++){
        string targetStr = adresses[k].substr(0, 1);
        double energyGV = stoi(regex_replace(adresses[k], regex(R"([\D])"), ""));
        string energyString = to_string(energyGV);

        double factor = 1.169;
        if(targetStr == "N"){
            if(energyGV == 771 || energyGV == 879) continue;
                //createFileN15(adresses[k], energyGV, factor, Ion("N15"));
                //angularCross(adresses[k]);
                /*
                auto result = AlphacmEfitter(adresses[k], factor);
                auto current = findCurrent(adresses[k]);
                GV.push_back(energyGV);
                counts.push_back(result[4]);
                countErrors.push_back(result[5]);
                Vcharge.push_back(current[0]);
                solidAngle.push_back(result[6]);
                totalCount.push_back(result[7]);
                angularCross(adresses[k]);
                 */

        }
    }

    string saveto = "Alpha0Cross.txt";
    ofstream mytxt (saveto);
    mytxt << "GV\tCounts\tCountErr\tVCharge\tSolidAngle\n";
    //skriv til en txt
    for(int k = 0; k < i; k++){
        string targetStr = adresses[k].substr(0, 1);
        int energyGV = stoi(regex_replace(adresses[k], regex(R"([\D])"), ""));
        if(energyGV == 771 || energyGV == 879) continue;
        //createFileN15(adresses[k], energyGV, 1.169, Ion("N15"));
        /*
        if(targetStr == "N"){
            mytxt << to_string(GV[k]) + "\t" + to_string(counts[k]) + "\t" + to_string(countErrors[k]) + "\t" +
                     to_string(Vcharge[k]) + "\t" + to_string(solidAngle[k]) + "\t" + to_string(totalCount[k]) + "\n";
            if(energyGV == 880 || energyGV == 1034){
                angularCross(adresses[k]);
            }
        }*/
    }
    createFileN15("N879gvm.root", 879, 1.169, Ion("N15"));
    mytxt.close();
    //createFileN15("N879gvm.root", 879, 1.169, Ion("N15"));
    //createFileN15("N880gvm.root", 880, 1.169, Ion("N15"));
    //angularCross("N880gvm.root");
}