
#include "QGP2to2Total.h"
#include <cmath>

QGP2to2Total::QGP2to2Total(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


void QGP2to2Total::analyticRates(
            double T, std::vector<double> &Eq,
            std::vector<double> &eqrate_ptr) {
    const double T2 =T*T;
    const double aem=1/137;
    const double gs =2;

    for (unsigned int i = 0; i < Eq.size(); i++) {
        double k = Eq[i];
        double x = k/T;
        double c2= (0.041/x)-0.3615+1.01*exp(-1.35*x);
        double kkT = 4/(3*pow((2*M_PI),3))*aem*gs*gs*T2
                    *(1/(exp(x)+1))*(log(sqrt(3)/gs)+log(2*x)/2+c2);
        eqrate_ptr[i] = kkT;
    }
}

void QGP2to2Total::NetBaryonCorrection(
        double T, double muB, std::vector<double> &Eq,
        std::vector<double> &eqrate_ptr) {

    if (std::abs(muB) < 1e-8)
        return;

    double kmuB = (1+(1/(M_PI*M_PI))*(muB/T)*(muB/T));

    for (unsigned int i = 0; i < Eq.size(); i++) {
        eqrate_ptr[i] *= kmuB;
    }
}
