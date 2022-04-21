#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/setup/SingleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Default.h>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include <Math/Vector3D.h>
#include <TROOT.h>
#include <ctime>
#include <libconfig.h++>
#include "include/Hit.h"
#include "include/runner.h"
#include <regex>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;
using namespace libconfig;

class MyAnalysisCutOff : public AbstractSortedAnalyzer {
public:
    //Initialiser variable
    int NUM;
    TTree *t;
    unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
    unique_ptr<DynamicBranchVector<double>> v_E, v_BE, v_FE, v_theta, v_dE, v_solang, v_cmE, v_cmE2, v_angDiff, v_pCM;
    unique_ptr<DynamicBranchVector<short>> v_i;
    unique_ptr<DynamicBranchVector<short>> v_F, v_B;
    unique_ptr<DynamicBranchVector<double>> v_ang, v_SAng;
    unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;

    UInt_t mul{}, TPATTERN{}, TPROTONS{}, EGPS{};
    vector<Hit> hits;

    unique_ptr<EnergyLossRangeInverter> SiCalc;
    vector<unique_ptr<EnergyLossRangeInverter>> targetCalcs;
    unique_ptr<EnergyLossRangeInverter> SiCalcC12;
    vector<unique_ptr<EnergyLossRangeInverter>> targetCalcsC12;
    vector<double> deadlayerF, deadlayerB, deadlayerP;
    Target &target;
    bool simulation = false;
    double accEnergy;
    TVector3 beta;
    TVector3 beta2;

    //Constructor for analyseklassen. Vi initialiserer det TTree, som vi vil ende med at gemme alle vores ting i.
    //Der laves også nogle energitabsberegninger. Gad vide mon hvad de skal bruges til.
    MyAnalysisCutOff(Target &target, TFile *output, double in, double factor, Ion targetIon) : target(target) {
        NUM = 0;
        accEnergy = in*factor;

        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);

        v_dir = make_unique<DynamicBranchVector<TVector3>>(*t, "dir");
        v_pos = make_unique<DynamicBranchVector<TVector3>>(*t, "pos");

        v_angDiff = make_unique<DynamicBranchVector<double>>(*t, "angDiff", "mul");
        v_pCM = make_unique<DynamicBranchVector<double>>(*t, "pCM", "mul");

        v_theta = make_unique<DynamicBranchVector<double>>(*t, "theta", "mul");
        v_ang = make_unique<DynamicBranchVector<double>>(*t, "angle", "mul");
        v_SAng = make_unique<DynamicBranchVector<double>>(*t, "scatterAngle", "mul");
        v_solang = make_unique<DynamicBranchVector<double>>(*t, "solidAngle", "mul");

        v_E = make_unique<DynamicBranchVector<double>>(*t, "E", "mul");
        v_BE = make_unique<DynamicBranchVector<double>>(*t, "BE", "mul");
        v_FE = make_unique<DynamicBranchVector<double>>(*t, "FE", "mul");

        v_FT = make_unique<DynamicBranchVector<double>>(*t, "FT", "mul");
        v_BT = make_unique<DynamicBranchVector<double>>(*t, "BT", "mul");

        v_dE = make_unique<DynamicBranchVector<double>>(*t, "dE", "mul");

        v_i = make_unique<DynamicBranchVector<short>>(*t, "id", "mul");

        v_F = make_unique<DynamicBranchVector<short>>(*t, "FI", "mul");
        v_B = make_unique<DynamicBranchVector<short>>(*t, "BI", "mul");

        v_cmE = make_unique<DynamicBranchVector<double>>(*t, "cmE", "mul");
        v_cmE2 = make_unique<DynamicBranchVector<double>>(*t, "cmE2", "mul");

        beta = TMath::Sqrt((accEnergy+PROTON_MASS)*(accEnergy+PROTON_MASS)-PROTON_MASS*PROTON_MASS)/(accEnergy+PROTON_MASS + targetIon.getMass()) * TVector3(0,0,1);

        cout << beta.Z() << endl;
        beta2 = constructBeamVector(Ion("H1"),targetIon,accEnergy).BoostVector();
        cout << beta2.Z() << endl;

        //t->Branch("TPATTERN", &TPATTERN);
        //t->Branch("TPROTONS", &TPROTONS);
        //t->Branch("EGPS", &EGPS);

        SiCalc = defaultRangeInverter("He4", "Silicon");
        SiCalcC12 = defaultRangeInverter("C12", "Silicon");
        for (auto &layer: target.getLayers()) {
            targetCalcsC12.push_back(defaultRangeInverter(Ion::predefined("C12"), layer.getMaterial()));
        }
    }

    //Setup bliver kørt som det første under run(). Først kaldes setup funktionen fra den abstrakte analyse klasse.
    //men i den abstrakte klasse er setup ikke initialiseret, så hvad mon det betyder? Efter det loopes der over
    //alle detektorene, og så finder vi længden af dead layers for hver detektor. Det gemmer vi så i "deadlayer"
    //-vektorene.
    void setup(const SortedSetupOutput &output) override {
        AbstractSortedAnalyzer::setup(output);
        for (size_t i = 0; i < output.dssdCount(); ++i) {
            auto dl = getFrontDeadLayer(output.getDssdOutput(i).detector());
            auto dlB = getBackDeadLayer(output.getDssdOutput(i).detector());
            //auto dlP = getFrontDeadLayer(output.getSingleOutput(i).detector());
            deadlayerF.push_back(dl);
            deadlayerB.push_back(dlB);
            //deadlayerP.push_back(dlP);
        }
    }

    //Hjælpefunktion til analyze (analyze bliver kørt for hvert event i output)
    void findHits() {
        //først finder jeg multipli

        //Vi looper over alle detektorene
        for (size_t i = 0; i < output.dssdCount(); i++) {
            //Hent outputs fra hver detektor. Find også multipliciteten i detektoren for det pågældende event?
            auto &o = output.getDssdOutput(i);
            //auto &p = output.getSingleOutput(i);
            auto &d = o.detector();
            auto m = AUSA::mul(o);

            //Der loopes over alle "events" indeni eventet, det vil sige multipliciteten i det givne event.
            for (UInt_t j = 0; j < m; j++) {
                //initialiser et "hit".
                Hit hit;

                //giv hittet dets dE, forskellen mellem frontenergi og backenergy for hvert matchede event.
                auto dE = fEnergy(o, j) - bEnergy(o, j);
                hit.dE = dE;

                //for Dssd detektorene regnes energien ud for eventet (gennemsnittet af den afsatte energi i
                //front og back) og front og back energier assignes til nogle variable.
                auto eFDssd = fEnergy(o, j);
                auto eBDssd = bEnergy(o, j);
                auto eDssd = energy(o,j);
                //  auto ePad = p.energy(0);

                // vi assigner de detektorsegments, eventet er sket i. Det gemmes også i hit.
                auto BI = bSeg(o, j);
                auto FI = fSeg(o, j);
                hit.fseg = short(FI);
                hit.bseg = short(BI);

                // vi finder koordinatet i 3d-koordsystemet som svarer til pixelpositionen, der blev ramt her.
                // vi kan så finde retningen (en vektor) for dette event. Dette tillader også, at vi kan finde
                // vinkler.
                auto position = o.detector().getUniformPixelPosition(FI, BI);
                auto origin = target.getCenter();
                hit.position = position;
                auto direction = (position - origin).Unit();
                hit.direction = direction;
                hit.theta = hit.direction.Theta();
                hit.solidAngle = o.detector().getPixelSolidAngle(FI,BI);

                //tilføj en vektor, der beskriver beamens retning. Beregn vinklen mellem denne og direction, og assign
                //det til hittet.
                auto beamDirection = TVector3(0,0,1);
                auto scatterAngle = hit.direction.Angle(beamDirection);
                hit.scatterAngle = scatterAngle;

                //assign de tider, hittet er sket i.
                if (!simulation) {
                    hit.TF = fTime(o, j);
                    hit.TB = bTime(o, j);
                    //    hit.TPad = p.time(0);
                } else {
                    hit.TF = 42;
                    hit.TB = 42;
                    //  hit.TPad = 42;
                }

                // angle er vinklen mellem detektorens normalvektor og retningen, partiklen kommer ind i. Denne
                // vinkel kan så bruges til at udregne hvor lang en rejse, partiklen har gennem dødlaget på detektoren
                auto angle = hit.direction.Angle(-d.getNormal());
                hit.angle = angle;

                auto tF = deadlayerF[i] / abs(cos(angle));
                auto tB = deadlayerB[i] / abs(cos(angle));
                //auto tP = deadlayerP[i] / abs(cos(angle));

                //initialiser nogle variable (lidt mærkeligt det her? De kunne bare være initialiseret som de
                //originale energier, eDssd osv - de bliver tillagt senere
                double E = 0.0;
                double FE = 0.0;
                double BE = 0.0;

                //assigner energier til hittet efter der er blevet lavet energikorrektioner følgende det døde
                //lag i detektoren.
                E += eDssd;
                E += SiCalc->getTotalEnergyCorrection(E, tF);
                FE += eFDssd;
                FE += SiCalc->getTotalEnergyCorrection(FE, tF);
                BE += eBDssd;
                BE += SiCalc->getTotalEnergyCorrection(BE, tF);

                hit.E = E;
                hit.BE = BE;
                hit.FE = FE;

                //jeg ved ikke hvad der sker her, det er noget mere energikorrektion.
                auto &from = position;
                for (auto &intersection: target.getIntersections(from, target.getCenter())) {
                    auto &calc = targetCalcs[intersection.index];
                    hit.E += calc->getTotalEnergyCorrection(hit.E, intersection.transversed);
                    hit.FE += calc->getTotalEnergyCorrection(hit.FE, intersection.transversed);
                    hit.BE += calc->getTotalEnergyCorrection(hit.BE, intersection.transversed);
                }

                //assign detektorindekset til hittet. Der laves også en Lorentz-vektor, ligner godt nok, at vi
                //antager, at det er alpha-partikler - skal nok ændres hvis det skal bruges.
                hit.index = i;
                hit.lVector = {TMath::Sqrt((hit.E + ALPHA_MASS)*(hit.E + ALPHA_MASS) - ALPHA_MASS*ALPHA_MASS) * hit.direction, hit.E + ALPHA_MASS};
                hit.lVector.Boost(-1*beta);
                hit.cmEnergy = hit.lVector[3] - ALPHA_MASS;

                hit.lVector2 = {sqrt(2 * hit.E * ALPHA_MASS) * hit.direction, hit.E + ALPHA_MASS};
                hit.lVector2.Boost(-1*beta2);
                hit.cmEnergy2 = hit.lVector2[3] - ALPHA_MASS;

                //vedhæft dette hit til vores liste af hits.
                hits.emplace_back(move(hit));
            }
        }
    }

    vector<Hit> findCoincidences(vector<Hit> hits){
        if( hits.size() != 2){
            return {};
        }
        return hits;
    }

    //hjælpefunktion til Analyze. Bliver kørt efter findhits.
    void doAnalysis() {
        //hvis der ikke er nogen hits i dette event stopper vi med det samme. Men hvordan ville der nogensinde
        //være et event hvor der ikke var nogle hits?
        auto coinHits = findCoincidences(hits);
        if (coinHits.empty()) return;
        //multipliciteten af eventet er antallet af hits.
        mul = coinHits.size();
        int i = 0;
        //for hvert hit lægger vi hittets information ind i vores dynamicbranches
        for (auto &hit: coinHits) {
            if(i == 0){
                v_angDiff -> add(hit.lVector.Vect().Angle(coinHits[1].lVector.Vect())*TMath::RadToDeg());
                v_pCM -> add((hit.lVector.Vect() + coinHits[1].lVector.Vect()).Mag());
            }
            if(i == 1){
                v_angDiff -> add(hit.lVector.Vect().Angle(coinHits[0].lVector.Vect())*TMath::RadToDeg());
                v_pCM -> add((hit.lVector.Vect() - coinHits[0].lVector.Vect()).Mag());
            }
            i++;
            v_pos->add(hit.position);
            v_dir->add(hit.direction);
            v_theta->add(hit.theta * TMath::RadToDeg());
            v_SAng->add(hit.scatterAngle * TMath::RadToDeg());

            v_E->add(hit.E);
            v_BE->add(hit.BE);
            v_FE->add(hit.FE);

            v_dE->add(hit.dE);
            v_ang->add(hit.angle * TMath::RadToDeg());
            v_solang -> add(hit.solidAngle);

            v_i->add(static_cast<short>(hit.index));
            v_F->add(hit.fseg);
            v_B->add(hit.bseg);
            v_FT->add(hit.TF);
            v_BT->add(hit.TB);
            v_cmE->add(hit.cmEnergy);
            v_cmE2->add(hit.cmEnergy2);
        }
    }

    //denne kode bliver kørt for hvert event. Først clearer vi hvad helt præcist? Alle dynamic branch vectors.
    //Så assigner vi til TPATTERN, TPROTONS og EGPS. Jeg tror ikke det betyder noget. Så kører vi FindHits på det
    //pågældende event, og så kører vi analysen på den pågældende event. Alt informationen, der ligger i vores
    // dynamic branch vectors bliver fyldt ind i t som er et TTree, så clearer vi hits og lægger en til NUM, som
    // nok holder styr på, hvor mange events der er.
    void analyze() override {
        clear();
        //TPATTERN = output.getScalerOutput("TPATTERN").getValue();
        //TPROTONS = output.getScalerOutput("TPROTONS").getValue();
        //EGPS = output.getScalerOutput("EGPS").getValue();
        findHits();
        doAnalysis();
        if (mul > 0) { t->Fill(); }
        hits.clear();
        NUM++;
    }

    static TLorentzVector constructBeamVector(const Ion& beam,
                                              const Ion& targetIon,
                                              double beamEnergy) {
        TLorentzVector plbeam( TVector3(0,0,sqrt(2*beamEnergy*beam.getMass())), beamEnergy+beam.getMass() );
        TLorentzVector pltarget( TVector3(0,0,0), targetIon.getMass() );
        return plbeam + pltarget;
    }


    void clear() {
        mul = 0;
        AUSA::clear(
                *v_E, *v_theta,
                *v_i, *v_FE, *v_BE,
                *v_F, *v_B, *v_SAng,
                *v_ang, *v_pos, *v_dir,
                *v_dE, *v_FT, *v_BT, *v_cmE, *v_cmE2, *v_angDiff, *v_pCM
        );
    }

    //bliver kaldt efter sidste event. Den abstrakte terminate kaldes (hvorfor?) og jeg tror gDirectoyry WriteTObject
    //gemmer vores træ ved at skrive den ind i hukommelsen?
    void terminate() override {
        AbstractSortedAnalyzer::terminate();
        gDirectory->WriteTObject(t);
    }
};

template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

//essentielt main
void createFileCutOff(string in, double energyGV, double factor, Ion targetIon){
    //læs setup og target filer.
    string targetStr = in.substr(0,1);

    auto setup = JSON::readSetupFromJSON("setup/setup.json");
    auto target = JSON::readTargetFromJSON("setup/"+targetStr +"target.json");

    SortedReader reader{*setup};
    reader.add("match/" + in);
    reader.setVerbose(true);

    string analyzed = "analyzed/"+targetStr+to_string_with_precision(energyGV,0) +"gv.root";
    TString outfile = analyzed;

    cout << "Reading from: " << in << endl;
    cout << "Printing to:  " << outfile << endl;

    TFile output(outfile, "RECREATE");
    auto analysis = make_shared<MyAnalysisCutOff>(target, &output, energyGV, factor, targetIon);
    reader.attach(analysis);
    reader.run();
}