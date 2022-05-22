#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Default.h>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include <Math/Vector3D.h>
#include <TROOT.h>
#include <ctime>
#include "include/Hit.h"
#include "include/runner.h"

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;

class MyAnalysisCutOff : public AbstractSortedAnalyzer {
public:
    //Initialiser variable
    int NUM;
    TTree *t;
    unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
    unique_ptr<DynamicBranchVector<double>> v_E, v_BE, v_FE, v_theta, v_dE, v_solang, v_cmE, v_cmE2, v_angDiff, v_pCM, v_timeDiff, v_recoilE;
    unique_ptr<DynamicBranchVector<short>> v_i;
    unique_ptr<DynamicBranchVector<short>> v_F, v_B;
    unique_ptr<DynamicBranchVector<double>> v_ang, v_SAng;
    unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;
    unique_ptr<DynamicBranchVector<short>> v_alpha;

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
    int GV;
    TVector3 beta;
    TVector3 beta2;
    double c12_mass;

    //Constructor for analyseklassen. Vi initialiserer det TTree, som vi vil ende med at gemme alle vores ting i.
    //Der laves også nogle energitabsberegninger. Gad vide mon hvad de skal bruges til.
    MyAnalysisCutOff(Target &target, TFile *output, double in, double factor, Ion targetIon, Ion recoilIon) : target(target) {
        NUM = 0;
        GV = in;
        accEnergy = in*factor;
        c12_mass = recoilIon.getMass();

        t = new TTree("a", "a");
        t->Branch("mul", &mul);
        t->Branch("num", &NUM);

        v_dir = make_unique<DynamicBranchVector<TVector3>>(*t, "dir");
        v_pos = make_unique<DynamicBranchVector<TVector3>>(*t, "pos");

        v_angDiff = make_unique<DynamicBranchVector<double>>(*t, "angDiff", "mul");
        v_timeDiff = make_unique<DynamicBranchVector<double>>(*t, "timeDiff", "mul");
        v_pCM = make_unique<DynamicBranchVector<double>>(*t, "pCM", "mul");
        v_recoilE = make_unique<DynamicBranchVector<double>>(*t, "recoilE", "mul");
        v_alpha = make_unique<DynamicBranchVector<short>>(*t, "canBeAlpha", "mul");

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
        SiCalcC12 = defaultRangeInverter(recoilIon.getName(), "Silicon");
        for (auto &layer: target.getLayers()) {
            targetCalcs.push_back(defaultRangeInverter(Ion::predefined("He4"), layer.getMaterial()));
            targetCalcsC12.push_back(defaultRangeInverter(Ion::predefined(recoilIon.getName()), layer.getMaterial()));
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
                if(GV == 879 || GV == 771){
                    eDssd = eBDssd;
                }
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
                hit.thickness = tF;

                //initialiser nogle variable (lidt mærkeligt det her? De kunne bare være initialiseret som de
                //originale energier, eDssd osv - de bliver tillagt senere
                double E = 0.0;
                double FE = 0.0;
                double BE = 0.0;
                E += eDssd;
                FE += eFDssd;
                BE += eBDssd;
                //assigner energier til hittet efter der er blevet lavet energikorrektioner følgende det døde
                //lag i detektoren.
                hit.E = E;
                hit.BE = BE;
                hit.FE = FE;
                hit.index = i;
                //vedhæft dette hit til vores liste af hits.
                hits.emplace_back(move(hit));
            }
        }
    }

    Hit correctEnergy(Hit hit){
        auto A = hit.canBeAlpha;
        auto tF = hit.thickness;
        auto from = hit.position;
        if(A){
            hit.E += SiCalc->getTotalEnergyCorrection(hit.E, tF);
            hit.FE += SiCalc->getTotalEnergyCorrection(hit.FE, tF);
            hit.BE += SiCalc->getTotalEnergyCorrection(hit.BE, tF);
            for (auto &intersection: target.getIntersections(from, target.getCenter())) {
                auto &calc = targetCalcs[intersection.index];
                hit.E += calc->getTotalEnergyCorrection(hit.E, intersection.transversed);
                hit.FE += calc->getTotalEnergyCorrection(hit.FE, intersection.transversed);
                hit.BE += calc->getTotalEnergyCorrection(hit.BE, intersection.transversed);
            }
        }
        else{
            hit.E += SiCalcC12->getTotalEnergyCorrection(hit.E, tF);
            hit.FE += SiCalcC12->getTotalEnergyCorrection(hit.FE, tF);
            hit.BE += SiCalcC12->getTotalEnergyCorrection(hit.BE, tF);
            for (auto &intersection: target.getIntersections(from, target.getCenter())) {
                auto &calc = targetCalcsC12[intersection.index];
                hit.E += calc->getTotalEnergyCorrection(hit.E, intersection.transversed);
                hit.FE += calc->getTotalEnergyCorrection(hit.FE, intersection.transversed);
                hit.BE += calc->getTotalEnergyCorrection(hit.BE, intersection.transversed);
            }
        }
        return hit;
    }

    //hjælpefunktion til Analyze. Bliver kørt efter findhits.
    void doAnalysis() {
        mul = 0;
        if(hits.size() < 2){
            return;
        };
        //for hvert hit lægger vi hittets information ind i vores dynamicbranches
        for (int i = 0; i < hits.size(); i++) {
            for (int j = i+1; j < hits.size(); j++) {
                //her kan betingelser indsættes
                if(hits.at(i).index == hits.at(j).index) continue;
                Hit alpha;
                Hit c12;
                //hittet med størst energi er alphaen.
                if(hits[i].E > hits[j].E){
                    alpha = hits[i];
                    c12 = hits[j];
                    alpha.canBeAlpha = true;
                    c12.canBeAlpha = false;
                }
                else{
                    alpha = hits[j];
                    c12 = hits[i];
                    alpha.canBeAlpha = true;
                    c12.canBeAlpha = false;
                }

                //lav energikorrektion. Funktionen tager højde for om der er tale om en alpha eller c12.
                //alpha = correctEnergy(alpha);
                //c12 = correctEnergy(c12);

                //skab lorentzvektorer
                alpha.lVector = {TMath::Sqrt((alpha.E + ALPHA_MASS)*(alpha.E + ALPHA_MASS) - ALPHA_MASS*ALPHA_MASS) * alpha.direction, alpha.E + ALPHA_MASS};
                c12.lVector = {TMath::Sqrt((c12.E + c12_mass)*(c12.E + c12_mass) - c12_mass*c12_mass) * c12.direction, c12.E + c12_mass};

                //boost dem til CM
                alpha.lVector.Boost(-1*beta);
                c12.lVector.Boost(-1*beta);

                //jeg kan nu finde deres CM energier:
                alpha.cmEnergy = alpha.lVector[3] - ALPHA_MASS;
                c12.cmEnergy = c12.lVector[3] - c12_mass;

                //jeg kan finde deres fælles tidslige afstand, impuls i CM og vinkeladskillelse i CM
                auto timeDiff = abs(alpha.TF - c12.TF);
                auto cmP = (alpha.lVector + c12.lVector).Vect().Mag();
                auto angDiff = alpha.lVector.Vect().Angle(c12.lVector.Vect())*TMath::RadToDeg();

                if(angDiff < 160 || cmP > 50000) continue;

                if(hits.at(i).index == hits.at(j).index) continue;

                v_alpha->add(1); //det første hit er alpha - tilføj 1 "true"
                v_alpha->add(0); //det andet hit er rekylet

                //nu kan jeg lægge værdierne ind i træet.
                v_pos->add(alpha.position);
                v_pos->add(c12.position);

                v_dir->add(alpha.direction);
                v_theta->add(alpha.theta * TMath::RadToDeg());
                v_SAng->add(alpha.scatterAngle * TMath::RadToDeg());

                v_dir->add(c12.direction);
                v_theta->add(c12.theta * TMath::RadToDeg());
                v_SAng->add(c12.scatterAngle * TMath::RadToDeg());

                v_E->add(alpha.E);
                v_BE->add(alpha.BE);
                v_FE->add(alpha.FE);

                v_dE->add(alpha.dE);
                v_ang->add(alpha.angle * TMath::RadToDeg());
                v_solang -> add(alpha.solidAngle);

                v_E->add(c12.E);
                v_BE->add(c12.BE);
                v_FE->add(c12.FE);

                v_dE->add(c12.dE);
                v_ang->add(c12.angle * TMath::RadToDeg());
                v_solang -> add(c12.solidAngle);

                v_i->add(static_cast<short>(alpha.index));
                v_F->add(alpha.fseg);
                v_B->add(alpha.bseg);
                v_FT->add(alpha.TF);
                v_BT->add(alpha.TB);
                v_cmE->add(alpha.cmEnergy);
                v_cmE2->add(alpha.cmEnergy2);

                v_i->add(static_cast<short>(c12.index));
                v_F->add(c12.fseg);
                v_B->add(c12.bseg);
                v_FT->add(c12.TF);
                v_BT->add(c12.TB);
                v_cmE->add(c12.cmEnergy);
                v_cmE2->add(c12.cmEnergy2);

                //læg hver ind to gange - jeg vil have at tidsforskellen står for hver partikel.
                v_timeDiff->add(timeDiff);
                v_timeDiff->add(timeDiff);

                v_angDiff->add(angDiff);
                v_angDiff->add(angDiff);

                v_pCM->add(cmP);
                v_pCM->add(cmP);

                v_recoilE->add(c12.cmEnergy);
                v_recoilE->add(alpha.cmEnergy);

                //jeg har lagt 2 elementer ind i træet; mul skal øges to gange.
                mul++;
                mul++;
            }
        }
    }

    //denne kode bliver kørt for hvert event. Først clearer vi hvad helt præcist? Alle dynamic branch vectors.
    //Så assigner vi til TPATTERN, TPROTONS og EGPS. Jeg tror ikke det betyder noget. Så kører vi FindHits på det
    //pågældende event, og så kører vi analysen på den pågældende event. Alt informationen, der ligger i vores
    // dynamic branch vectors bliver fyldt ind i t som er et TTree, så clearer vi hits og lægger en til NUM, som
    // nok holder styr på, hvor mange events der er.
    void analyze() override {
        clear();
        findHits();
        doAnalysis();
        if (mul > 0) {t->Fill(); }
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
                *v_dE, *v_FT, *v_BT, *v_cmE, *v_cmE2, *v_angDiff, *v_pCM, *v_alpha, *v_timeDiff,
                *v_recoilE
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
void createFileCoin(string in, double energyGV, double factor, Ion targetIon, Ion recoilIon){
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
    auto analysis = make_shared<MyAnalysisCutOff>(target, &output, energyGV, factor, targetIon, recoilIon);
    reader.attach(analysis);
    reader.run();
}