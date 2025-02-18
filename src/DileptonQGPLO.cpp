#include "DileptonQGPLO.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_fermi_dirac.h>

#include <cmath>

#include "data_struct.h"

using PhysConsts::alphaEM;
using PhysConsts::hbarC;

void error_handler_lo(
    const char *reason, const char *file, int line, int gsl_errno) {
    // fprintf(stderr, "GSL ERROR: %s:%d: %s (error code: %d)\n", file, line,
    // reason, gsl_errno);
}
DileptonQGPLO::DileptonQGPLO(
    std::shared_ptr<ParameterReader> paraRdr_in, std::string emissionProcess)
    : ThermalDilepton {paraRdr_in, emissionProcess} {}

double DileptonQGPLO::nB(double x) {
    double e = exp(-x);
    return e / (1. - e);
};

double DileptonQGPLO::mllFactor(double x) {
    // lepton pair kinematic factor
    if (4. * x > 1.)
        return 0.;
    else
        return (1. + 2. * x) * sqrt(1. - 4. * x);
}

// double DileptonQGPLO::l1f(double x) { return +gsl_sf_fermi_dirac_0(-x); }
//
// double DileptonQGPLO::l2f(double x) { return -gsl_sf_fermi_dirac_1(-x); }
//
// double DileptonQGPLO::l3f(double x) { return -gsl_sf_fermi_dirac_2(-x); }

double DileptonQGPLO::l1f(double x) {
    gsl_set_error_handler(&error_handler_lo);
    gsl_sf_result result;
    int status = gsl_sf_fermi_dirac_0_e(-x, &result);

    if (status != GSL_SUCCESS) {
        // printf("An error occurred: %f\n", result.val);
        // not printf, set result.val=0.0
    }
    return result.val;
}
double DileptonQGPLO::l2f(double x) {
    gsl_set_error_handler(&error_handler_lo);
    gsl_sf_result result;
    int status = gsl_sf_fermi_dirac_1_e(-x, &result);

    if (status != GSL_SUCCESS) {
        // printf("An error occurred: %f\n", result.val);
        // not printf, set result.val=0.0
    }
    return -result.val;
}
double DileptonQGPLO::l3f(double x) {
    gsl_set_error_handler(&error_handler_lo);
    gsl_sf_result result;
    int status = -gsl_sf_fermi_dirac_2_e(-x, &result);
    if (status != GSL_SUCCESS) {
        // printf("An error occurred: %f\n", result.val);
        // not printf, set result.val=0.0
    }

    return -result.val;
}

void DileptonQGPLO::rho_LO(
    double o, double k, double mu, double &rT, double &rL) {
    // leading order result, see (2.4) of 1910.09567
    // (mu dependence not actually published yet...)
    double rV, r00;
    double K2 = o * o - k * k;  // note: omega, k are in units of T!
    double kp = .5 * (o + k), km = .5 * std::abs(o - k);
    double somk = 1;
    if ((o - k) < 0) somk = -1;

    rV = l1f(kp + mu) + l1f(kp - mu) - l1f(km + mu) - l1f(km - mu);
    rV = rV / k + .5 * (somk + 1.);
    rV *= K2 * 3. * OOFP;

    r00 = l3f(kp + mu) + l3f(kp - mu) - l3f(km + mu) - l3f(km - mu);
    r00 = r00 * 2. / k + l2f(kp + mu) + l2f(kp - mu)
          + somk * (l2f(km + mu) + l2f(km - mu));
    r00 = r00 * 6. + .5 * k * k * (somk + 1.);
    r00 *= -OOFP;

    rL = -K2 * r00 / (k * k);
    rT = .5 * (rV - rL);
}

// PRC. 93, 044902, 2016
void DileptonQGPLO::analyticRates(
    const double E, const double k, const double muB, const double T,
    const double m_l, double &rateTot, double &rateT, double &rateL) {
    double M2 = E * E - k * k;

    double prefactor = mllFactor(m_l * m_l / M2) * (2. / 3.) * pow(hbarC, -4.)
                       * pow(alphaEM, 2.);

    double rhoT_app, rhoL_app;
    rho_LO(E / T, k / T, muB / T, rhoT_app, rhoL_app);

    double prefac =
        prefactor * nB(E / T) * pow(T, 2.) / (3. * pow(M_PI, 3.) * M2);

    rateT = prefac * rhoT_app;
    rateL = prefac * rhoL_app;

    rateTot = 2. * rateT + rateL;
}
