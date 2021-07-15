
#include "HadronGasPipiBremsstrahlung.h"
#include <cmath>

HadronGasPipiBremsstrahlung::HadronGasPipiBremsstrahlung(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


void HadronGasPipiBremsstrahlung::analyticRates(
            double T, double muB, std::vector<double> &Eq,
            std::vector<double> &eqrate_ptr) {
    // parameterization taken from e-Print: 1411.7012 [hep-ph]
    double T2 = T*T;
    double T3 = T2*T;

    double abT = -16.28 + 62.45*T - 93.4*T2 - 7.5*T3;
    double bbT = -35.54 + 414.8*T - 2054*T2 + 3718.8*T3;
    double ybT = 0.7364 - 10.72*T + 56.32*T2 - 103.5*T3;
    double dbT = -2.51 + 58.152*T - 318.24*T2 + 610.7*T3;
    for (unsigned int i = 0; i < Eq.size(); i++) {
        double Eq_local = Eq[i];
        double logB = (abT + bbT*Eq_local + ybT*Eq_local*Eq_local
                       + dbT/(Eq_local + 0.2));
        eqrate_ptr[i] = exp(logB);
    }
}
