#include "HadronGasRhoSpectralFunction.h"

#include <cmath>

HadronGasRhoSpectralFunction::HadronGasRhoSpectralFunction(
    std::shared_ptr<ParameterReader> paraRdr_in, std::string emissionProcess)
    : ThermalPhoton {paraRdr_in, emissionProcess} {}

void HadronGasRhoSpectralFunction::analyticRates(
    double T, std::vector<double> &Eq, std::vector<double> &eqrate_ptr) {
    // parameterization taken from e-Print: 1411.7012 [hep-ph]
    double T2 = T * T;
    double T3 = T2 * T;

    double aT = -31.21 + 353.61 * T - 1739.4 * T2 + 3105 * T3;
    double bT = -5.513 - 42.2 * T + 333 * T2 + -570 * T3;
    double cT = -6.153 + 57 * T - 134.61 * T2 + 8.31 * T3;

    for (unsigned int i = 0; i < Eq.size(); i++) {
        double Eq_local = Eq[i];
        double logR = aT * Eq_local + bT + cT / (Eq_local + 0.2);  // log(R_0)
        eqrate_ptr[i] = exp(logR);
    }
}

void HadronGasRhoSpectralFunction::NetBaryonCorrection(
    double T, double muB, std::vector<double> &Eq,
    std::vector<double> &eqrate_ptr) {
    muB = std::abs(muB);
    if (muB < 1e-8) return;

    // set the upper limit of muB to 0.4 GeV for this parameterization
    muB = std::min(0.4, muB);

    double T2 = T * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;
    double muB2 = muB * muB;
    double muB3 = muB2 * muB;

    double nT = -0.04 + 2.3 * T - 12.8 * T2;
    double pT = 23.66 - 354 * T + 1175 * T2;
    double rT = -54.3 + 742.6 * T - 2350 * T2;
    double sT =
        (-22.11 + 808.7 * T - 11604.4 * T2 + 81700 * T3 - 282480 * T4
         + 384116 * T5);
    double vT = -1.6 - 121.7 * T + 1775 * T2 - 5516 * T3;
    double wT = -9.874 + 469 * T - 4371.5 * T2 + 11000 * T3;
    double alphaT =
        (84.784 - 3028.6 * T + 42434 * T2 - 291390 * T3 + 981000 * T4
         - 1295400 * T5);
    double betaT = 59.64 - 726.46 * T + 1093.4 * T2 + 4256 * T3;
    double etaT = -73.9 + 458.66 * T + 2450 * T2 - 12348 * T3;
    double dmT = nT * muB + pT * muB2 + rT * muB3;
    double kmT = sT * muB + vT * muB2 + wT * muB3;
    double mmT = alphaT * muB + betaT * muB2 + etaT * muB3;

    for (unsigned int i = 0; i < Eq.size(); i++) {
        double Eq_local = Eq[i];
        double logFrho = dmT - kmT / (Eq_local * Eq_local) - mmT / Eq_local;
        eqrate_ptr[i] *= exp(logFrho);
    }
}
