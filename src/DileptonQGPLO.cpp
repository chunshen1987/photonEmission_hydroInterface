#include "DileptonQGPLO.h"
#include "data_struct.h"
#include <cmath>

using PhysConsts::hbarC;

DileptonQGPLO::DileptonQGPLO(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalDilepton{paraRdr_in, emissionProcess} {}


// PRC. 93, 044902, 2016
void DileptonQGPLO::analyticRates(
    const double T, const double MInv, std::vector<double> &Eq,
    std::vector<double> &eqrate_ptr) {

    const double aem = 1./137.;
    const double Qu  = 2./3.;
    const double Qd  = -1./3.;
    const double Qs  = -1./3.;

    double prefac = (Qu*Qu+Qd*Qd+Qs*Qs)*aem*aem/(2.*pow(M_PI, 4))/pow(hbarC, 4);
    for (unsigned int i = 0; i < Eq.size(); i++) {
        double p = sqrt(Eq[i]*Eq[i] - MInv*MInv);
        double x = Eq[i]/T;
        double y = p/T;
        double fq = 1./(exp(x) - 1.);
        double log_r = log(cosh(0.25*(x+y))/cosh(0.25*(x-y)));

        eqrate_ptr[i] = prefac/y*fq*log_r;
    }
}

