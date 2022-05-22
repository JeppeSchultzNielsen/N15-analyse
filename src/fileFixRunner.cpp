#include <include/runner.h>
#include <vector>
#include <string>

int main(){
    vector<short> gvStrings = {813,835,860,879,898,920,950,975,1055,1400,1710,2137,2564};
    vector<UInt_t> startClocks = {3563650000, 3564900000, 3307100000, 2040000000, 3565970000, 3308500000,
                                  3567200000,3309800000,3312100000,3571100000,3822875000,3823700000,3824730000};
    for(int i = 0; i < gvStrings.size(); i++){
        fixFile(to_string(gvStrings[i]),startClocks[i]);
    }
}