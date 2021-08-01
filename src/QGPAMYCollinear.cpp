
#include "QGPAMYCollinear.h"
#include <cmath>

QGPAMYCollinear::QGPAMYCollinear(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


void QGPAMYCollinear::analyticRates(
            double T, std::vector<double> &Eq,
            std::vector<double> &eqrate_ptr) {
    const double T2 =T*T;
    const double aem=1/137;
    const double gs =2;
    
    for (unsigned int i = 0; i < Eq.size(); i++) {
        double k = Eq[i];
        double x = k/T;
        double kak = (4/3)*aem*gs*gs*T2*(1/(exp(x)+1));
        double cba = (sqrt(6)/2)*(((0.548*log(12.28+(1/x)))/pow(x,(3/2)))
                     +((0.133*x)/(sqrt(1+(x/16.27)))));
        double kamy = (1/pow((2*M_PI),3))*kak*cba;
        eqrate_ptr[i] = kamy;
    }
}

void QGPAMYCollinear::NetBaryonCorrection(
        double T, double muB, std::vector<double> &Eq,
        std::vector<double> &eqrate_ptr) {

    if (std::abs(muB) < 1e-8)
        return;

    for (unsigned int i = 0; i < Eq.size(); i++) {
        double x = Eq[i]/T;
        double famy = 1+((-0.008673644436239206+(0.0009992379594915468*x))*(muB/T))
                      /(1-(0.05587549314290918*x)+(0.20423253055206605*x*x))
                      +(0.24476805976626895+(0.04706777450080506*x)*(muB/T)*(muB/T))
                      /(1-(0.001394126856568418*x)+(0.20423253055206605*x*x));
        eqrate_ptr[i] *= famy;
    }
}
