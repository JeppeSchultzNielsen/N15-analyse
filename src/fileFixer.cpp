#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TSeqCollection.h>
#include <include/runner.h>
#include <string>

void fixFile(string gvString, UInt_t startClock){
    TString openstring = "match/N"+gvString+"gvm.root";
    TString savestring = "fixed/N"+gvString+"gvm.root";
    std::unique_ptr<TFile>  oldfile(TFile::Open(openstring));
    TTree* oldtree = oldfile->Get<TTree>("a101");


    UInt_t CLOCK;
    oldtree->SetBranchAddress("CLOCK", &CLOCK);

    TFile newfile(savestring, "recreate");
    auto newtree = oldtree->CloneTree(0);
    cout << startClock << endl;
    for (int i = 0; oldtree->LoadTree(i) >= 0; i++) {
        oldtree->GetEntry(i);
        if (CLOCK > startClock) {
            newtree->Fill();
        }
    }
    newtree->Print();
    newfile.Write();
}
