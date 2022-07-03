
#include "QGP2to2Total.h"
#include "data_struct.h"
#include <cmath>

using PhysConsts::hbarC;

QGP2to2Total::QGP2to2Total(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


void QGP2to2Total::analyticRates(
            double T, std::vector<double> &Eq,
            std::vector<double> &eqrate_ptr) {
    const double T2  = T*T;
    const double aem = 1./137.;
    const double gs  = 2.;

    const double prefac = 4./(3.*pow(2*M_PI, 3))*aem*gs*gs*T2/pow(hbarC, 4);
    const double logfac = log(sqrt(3.)/gs);

    for (unsigned int i = 0; i < Eq.size(); i++) {
        double x = Eq[i]/T;
        double c2= 0.041/x - 0.3615 + 1.01*exp(-1.35*x);
        double kkT = prefac*(1./(exp(x) + 1.))*(logfac + log(2.*x)/2. + c2);
        eqrate_ptr[i] = kkT;
    }
}

void QGP2to2Total::NetBaryonCorrection(
        double T, double muB, std::vector<double> &Eq,
        std::vector<double> &eqrate_ptr) {

    muB = std::abs(muB);
    if (muB < 1e-8) return;

    // set the upper limit to 0.4 GeV for this parameterization
    muB = std::min(0.4, muB);

    const double muB_over_T = muB/T;
    const double kmuB = 1. + 1./(M_PI*M_PI)*muB_over_T*muB_over_T;

    for (unsigned int i = 0; i < Eq.size(); i++) {
        eqrate_ptr[i] *= kmuB;
    }
}
