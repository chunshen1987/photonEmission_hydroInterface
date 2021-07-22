
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
    const double aem = 1/137;
    const double muB = 0.0;
    for (unsigned int i = 0; i < Eq.size(); i++) {
        double k = Eq[i];
        double x = k/T;

        double kak = (16/3)*aem*T2*(1/(exp(x)+1));
        double cba = (sqrt(6)/2)*(((0.548*log(12.28+(1/x)))/pow(x,(3/2)))
                     +((0.133*x)/(sqrt(1+(x/16.27)))));
        double kamy = (1/pow((2*M_PI),3))*kak*cba;
        
        double famy = 1+((-0.008673644436239206+(0.0009992379594915468*x))*(muB/T))
                      /(1-(0.05587549314290918*x)+(0.20423253055206605*x*x))
                      +(0.24476805976626895+(0.04706777450080506*x)*(muB/T)*(muB/T))
                      /(1-(0.001394126856568418*x)+(0.20423253055206605*x*x));
        double c2= (0.041/x)-0.3615+1.01*exp(-1.35*x);
        double kkT = 16/(3*pow((2*M_PI),3))*aem*T2*(1/(exp(x)+1))*(log(sqrt(3)/2)+log(2*x)/2+c2);
        double kmuB = kkT*famy;
        eqrate_ptr[i] = kmuB;
    }
}
