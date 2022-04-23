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

    for (int j = 0; j < 1; j++) {
        string targetStr = adresses[j].substr(0, 1);
        double energyGV = stoi(regex_replace(adresses[j], regex(R"([\D])"), ""));
        string energyString = to_string(energyGV);
        double factor = 1.164;
        cout << "hej" << endl;
        if (targetStr == "N") {
            createFileN15(adresses[j], energyGV, factor, Ion("N15"));

        }
    }
}