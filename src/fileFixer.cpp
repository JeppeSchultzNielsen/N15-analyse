#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TSeqCollection.h>
#include <include/runner.h>
#include <string>
#include <TObjString.h>
#include <TObjArray.h>
#include <TDirectoryFile.h>

void fixFile(string gvString, UInt_t startClock){
    TString openstring = "match/N"+gvString+"gvm.root";
    TString savestring = "fixed/N"+gvString+"gvm.root";
    std::unique_ptr<TFile>  oldfile(TFile::Open(openstring));
    TTree* oldtree = oldfile->Get<TTree>("a101");
    TObjString* hash;
    oldfile->GetObject("AUSALIB_HASH", hash);
    TObjString* branch = oldfile->Get<TObjString>("AUSALIB_BRANCH");
    TObjString* guid = oldfile->Get<TObjString>("GUID");
    TObjArray* detectors = oldfile->Get<TObjArray>("detectors");
    TDirectoryFile* sortQC = oldfile->Get<TDirectoryFile>("sortQC");
    TObjString* matcher_config = oldfile->Get<TObjString>("MATCHER_CONFIG");

    UInt_t CLOCK;
    oldtree->SetBranchAddress("CLOCK", &CLOCK);
    UInt_t N;
    oldtree->SetBranchAddress("___N___", &N);

    TFile newfile(savestring, "recreate");
    auto newtree = oldtree->CloneTree(0);
    cout << startClock << endl;
    long n = 0;
    for (int i = 0; oldtree->LoadTree(i) >= 0; i++) {
        oldtree->GetEntry(i);
        if (CLOCK > startClock) {
            N = n;
            newtree->Fill();
            n++;
        }
    }
    newtree->Print();
    newfile.Write();
    newfile.WriteObject(hash, "AUSALIB_HASH");
    newfile.WriteObject(branch, "AUSALIB_BRANCH");
    newfile.WriteObject(guid, "GUID");
    newfile.WriteObject(detectors, "detectors");
    newfile.WriteObject(sortQC, "sortQC");
    newfile.WriteObject(matcher_config, "MATCHER_CONFIG");
}
