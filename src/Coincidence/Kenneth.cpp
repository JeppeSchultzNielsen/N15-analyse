#include <iostream>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Ion.h>
#include <ausa/eloss/Default.h>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include <TLorentzVector.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include "include/Hit.h"


using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;

namespace {
    const auto ALPHA = new ParticleType{"He4"};
    const auto BE8 = new ParticleType{"Be8"};
    const auto C12 = new ParticleType{"C12"};
    const auto O16 = new ParticleType{"O16"};
}

class O16Analysis : public AbstractSortedAnalyzer {
public:
    O16Analysis(double beamEnergy, Target& target) : target(target) {
        intialFourVector = constructBeamVector(Ion{"H1"}, Ion{"N15"}, beamEnergy);
        cmBoost = -intialFourVector.BoostVector();

        t = new TTree("a", "a");
        t->Branch("mul", &mul);

        v_pC = std::make_unique<DynamicBranchVector<TVector3>>(*t, "pC");
        v_pA = std::make_unique<DynamicBranchVector<TVector3>>(*t, "pA");
        v_pCM = std::make_unique<DynamicBranchVector<TVector3>>(*t, "pCM");
        v_dP = std::make_unique<DynamicBranchVector<TVector3>>(*t, "dP");
        v_V = std::make_unique<DynamicBranchVector<double>>(*t, "V", "mul");
        v_vC = std::make_unique<DynamicBranchVector<double>>(*t, "VC", "mul");
        v_vA = std::make_unique<DynamicBranchVector<double>>(*t, "VA", "mul");
        v_vCcm = std::make_unique<DynamicBranchVector<double>>(*t, "VCcm", "mul");
        v_vAcm = std::make_unique<DynamicBranchVector<double>>(*t, "VAcm", "mul");

        v_DEC = std::make_unique<DynamicBranchVector<double>>(*t, "DEC", "mul");
        v_EC = std::make_unique<DynamicBranchVector<double>>(*t, "EC", "mul");
        v_ECcm = std::make_unique<DynamicBranchVector<double>>(*t, "ECcm", "mul");
        v_DEA = std::make_unique<DynamicBranchVector<double>>(*t, "DEA", "mul");
        v_EA = std::make_unique<DynamicBranchVector<double>>(*t, "EA", "mul");
        v_EAcm = std::make_unique<DynamicBranchVector<double>>(*t, "EAcm", "mul");
        v_ex = std::make_unique<DynamicBranchVector<double>>(*t, "ex", "mul");

        v_TC = std::make_unique<DynamicBranchVector<double>>(*t, "TC", "mul");
        v_TA = std::make_unique<DynamicBranchVector<double>>(*t, "TA", "mul");

        v_iC = std::make_unique<DynamicBranchVector<short>>(*t, "iC", "mul");
        v_iA = std::make_unique<DynamicBranchVector<short>>(*t, "iA", "mul");



        ACalc = defaultRangeInverter("He4", "Silicon");
        C12Calc = defaultRangeInverter("C12", "Silicon");
        for (auto& layer : target.getLayers()) {
            AtargetCalcs.push_back(defaultRangeInverter(Ion::predefined("He4"), layer.getMaterial()));
            C12targetCalcs.push_back(defaultRangeInverter(Ion::predefined("C12"), layer.getMaterial()));
        }
    }

    virtual void analyze() override {
        clear();

        findHits();

        doAnalysis();
    }

    double correctEnergy(ParticleType *const pType, Hit *hit,
                         EnergyLossRangeInverter &detCalc,
                         vector<unique_ptr<EnergyLossRangeInverter>> &targetCalcs) {

        auto E  = hit->deposited + detCalc.getTotalEnergyCorrection(hit->deposited, hit->thickness);

        auto& from = hit->position;
        TVector3 to{0,0,0};

        for (auto& intersection : target.getIntersections(from, to)) {
            // energy loss calculator for target layer
            auto& calc = targetCalcs[intersection.index];
            E += calc->getTotalEnergyCorrection(E, intersection.transversed);
        }
        return E;
    }

    TLorentzVector calc4Vector(double mass, double E, const TVector3& dir) {
        return TLorentzVector{sqrt(2 * E * mass) * dir, E + mass};
    }

    void doAnalysis() {
        if (hits.size() < 2) return;

        for (size_t i = 0; i < hits.size(); i++) {
            auto& h0 = hits[i];

            for (size_t j = i+1; j < hits.size(); j++) {
                auto& h1 = hits[j];
                if (h0.index == h1.index // No self coincidence
                    //|| abs(h0.T - h1.T) > 100E3 // < 100ns difference
                        ) continue;
                if(h0.index == 2 || h0.index == 3 || h1.index == 2 || h1.index == 3) continue;
                if(h0.index == 0 && h0.bseg == 13) continue;
                if(h1.index == 0 && h1.bseg == 13) continue;

                // Particle ID
                // Alphas has highest energy
                Hit *hC, *hA;

                if (h0.deposited < h1.deposited) {
                    hC = &h0;
                    hA = &h1;
                }
                else {
                    hC = &h1;
                    hA = &h0;
                }

                auto eC = correctEnergy(C12, hC, *C12Calc, C12targetCalcs);
                auto eA = correctEnergy(ALPHA, hA, *ACalc, AtargetCalcs);

                auto pC = calc4Vector(C12->mass, eC, hC->direction);
                auto pA = calc4Vector(ALPHA->mass, eA, hA->direction);

                v_dP->add((pC + pA - intialFourVector).Vect());

                pC.Boost(cmBoost);
                auto eC_cm = pC.E() - C12->mass;

                pA.Boost(cmBoost);
                auto eA_cm = pA.E() - ALPHA->mass;

                auto p = pC + pA;
//                p.Boost(cmBoost);

                auto m = p.M();
                auto eX = m - O16->mass;

                auto p3C = pC.Vect();
                auto p3A = pA.Vect();
                v_pC->add(p3C);
                v_pA->add(p3A);
                v_pCM->add(p.Vect());

                v_V->add(pA.Angle(p3C) * TMath::RadToDeg());
                v_vC->add(hC->theta * TMath::RadToDeg());
                v_vA->add(hA->theta * TMath::RadToDeg());
                v_vCcm->add(p3C.Theta() * TMath::RadToDeg());
                v_vAcm->add(p3A.Theta() * TMath::RadToDeg());

                v_DEC->add(hC->deposited);
                v_EC->add(eC);
                v_ECcm->add(eC_cm);

                v_DEA->add(hA->deposited);
                v_EA->add(eA);
                v_EAcm->add(eA_cm);

                v_ex->add(eX);

                v_TC->add(hC->T);
                v_TA->add(hA->T);

                v_iC->add(static_cast<short>(hC->index));
                v_iA->add(static_cast<short>(hA->index));

                mul++;
            }
        }
        t->Fill();
    }

    virtual void setup(const SortedSetupOutput& output) override {
        AbstractSortedAnalyzer::setup(output);

        for (size_t i = 0; i < output.dssdCount(); ++i) {
            deadLayer.push_back(getFrontDeadLayer(output.getDssdOutput(i).detector()));
        }
    }


    void findHits() {
        for (size_t i = 0; i < output.dssdCount(); i++) {
            auto& o = output.getDssdOutput(i);
            auto& d = o.detector();
            auto m = AUSA::mul(o);

            for (UInt_t j = 0; j < m; j++) {
                Hit hit;
                if(i == 1){
                    hit.deposited = bEnergy(o, j);
                }
                else hit.deposited = energy(o, j);

//                if (hit.deposited < 200) continue;

                const auto front = fSeg(o, j);
                const auto back = bSeg(o, j);

                hit.position = d.getUniformPixelPosition(front, back);
                hit.direction = hit.position.Unit();
                hit.theta = hit.direction.Theta();

                // CFD on back side of S3 - front side of W1
                hit.T = (i <= 1) ? bTime(o, j): fTime(o,j);


                auto angle = hit.direction.Angle(d.getNormal());
                hit.thickness = deadLayer[i] / abs(cos(angle));
                hit.index = i;

                hits.emplace_back(move(hit));
            }
        }
    }

    virtual void terminate() override {
        AbstractSortedAnalyzer::terminate();

        gDirectory->WriteTObject(t);
    }


    TLorentzVector constructBeamVector(const Ion& beam,
                                       const Ion& target,
                                       double beamEnergy) {
        TLorentzVector plbeam( TVector3(0,0,sqrt(2*beamEnergy*beam.getMass())), beamEnergy+beam.getMass() );
        TLorentzVector pltarget( TVector3(0,0,0), target.getMass() );
        return plbeam + pltarget;
    }

    void clear() {
        mul = 0;
        AUSA::clear(
                *v_pC, *v_pA, *v_pCM, *v_dP,
                *v_EC, *v_EA, *v_DEC, *v_DEA, *v_ex, *v_ECcm, *v_EAcm, *v_TC, *v_TA, *v_V, *v_vC, *v_vA, *v_vCcm, *v_vAcm,
                *v_iC, *v_iA
        );

        hits.clear();
    }

    TLorentzVector intialFourVector;
    TVector3 cmBoost;

    TTree* t;

    std::unique_ptr<DynamicBranchVector<TVector3>> v_pC, v_pA, v_pCM, v_dP;
    std::unique_ptr<DynamicBranchVector<double>> v_EC, v_EA, v_DEC, v_DEA,
            v_ex,
            v_ECcm, v_EAcm,
            v_TC, v_TA,
            v_V, v_vC, v_vA, v_vCcm, v_vAcm;
    std::unique_ptr<DynamicBranchVector<short>> v_iC, v_iA;
    UInt_t mul;

    vector<Hit> hits;

    // Silicon eloss calc. Use for detector DL
    std::unique_ptr<EnergyLossRangeInverter> ACalc, C12Calc;

    // Target layer eloss calc.
    vector<std::unique_ptr<EnergyLossRangeInverter>> AtargetCalcs, C12targetCalcs;
    vector<double> deadLayer;
    Target& target;
};


int main(int argc, char* argv[]) {

    auto setup = JSON::readSetupFromJSON("setup/setup.json");
    auto target = JSON::readTargetFromJSON("setup/Ntarget.json");

    vector<string> input;
    for (int i = 2; i < argc; i++) input.push_back(argv[i]);

    SortedReader reader{*setup};
    reader.add("match/N771gvm.root");
    reader.setVerbose(true);

    auto base = stripFileExtension(extractFileName(findLongestCommonSubstring(input)));

    TFile output("analyzed/kenneth771_O16.root", "RECREATE");

    reader.attach(std::make_shared<O16Analysis>(771*1.169, target));

    reader.run();
}
