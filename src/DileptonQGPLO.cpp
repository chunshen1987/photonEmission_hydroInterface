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
    const double T, const double MInv, const double Eq,
    double &eqrate) {

    //const double aem = 1./137.;
    //const double Qu  = 2./3.;
    //const double Qd  = -1./3.;
    //const double Qs  = -1./3.;

    //double prefac = (Qu*Qu+Qd*Qd+Qs*Qs)*aem*aem/(2.*pow(M_PI, 4))/pow(hbarC, 4);
    double prefac = 0.00012024462067159784;

    double p = sqrt(Eq*Eq - MInv*MInv);
    double x = Eq/T;
    double y = p/T;
    double fq = 1./(exp(x) - 1.);
    double log_r = log(cosh(0.25*(x+y))/cosh(0.25*(x-y)));

    eqrate = prefac/y*fq*log_r;
}

