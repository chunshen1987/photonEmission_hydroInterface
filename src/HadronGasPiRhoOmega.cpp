
#include "HadronGasPiRhoOmega.h"
#include <cmath>

HadronGasPiRhoOmega::HadronGasPiRhoOmega(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess
        ) : ThermalPhoton{paraRdr_in, emissionProcess} {}


void HadronGasPiRhoOmega::analyticRates(
            double T, double muB, std::vector<double> &Eq,
            std::vector<double> &eqrate_ptr) {
    // parameterization taken from e-Print: 1506.09205 [hep-ph]
    const double T2=T*T;
    const double T3=T2*T;

    const double a1=-35.8991+460.425*T-2592.04*T2+5342.32*T3;
    const double a2=-41.9725+601.952*T-3587.8 *T2+7604.97*T3;
    const double a3=0.740436-16.7159*T+133.526*T2-347.589*T3;
    const double a4= 2.00611-3.79343*T+29.3101*T2-72.8725*T3;
    const double a5=-8.33046+121.091*T-801.676*T2+1712.16*T3;
    const double a6=17.9029 -  388.5*T+2779.03*T2- 6448.4*T3;
    const double a7=-15.622 +340.651*T-2483.18*T2+5870.61*T3;

    const double b1=-29.6866+331.769*T-1618.66*T2+2918.53*T3;
    const double b2=-15.3332+90.2225*T-300.185*T2+428.386*T3;
    const double b3=-7.35061+109.288*T-630.396*T2+1227.69*T3;
    const double b4=-10.6044+  109.1*T-500.718*T2+872.951*T3;

    const double d1=-29.4663 +291.356*T-1301.27*T2+2102.12*T3;
    const double d2=-45.081  +688.929*T-4150.15*T2+8890.76*T3;
    const double d3=-0.260076+8.92875*T- 60.868*T2+ 136.57*T3;
    const double d4= 2.2663  -8.30596*T+49.3342*T2-90.8501*T3;
    const double d5= 10.2955 -317.077*T+2412.15*T2- 6020.9*T3;
    const double d6= 3.12251 -47.5277*T+ 222.61*T2-  241.9*T3;
    const double d7=-3.39045 +56.5927*T- 336.97*T2+622.756*T3;

    for (unsigned int i = 0; i < Eq.size(); i++) {
        double q0 = Eq[i];
        // pi+rho->gamma+omega
        double FFpiro = exp(a1*q0 +a2+a3*pow(q0,a4)+a5*pow((q0+a6),a7));
        // rho+omega->gamma+pi
        double FFomro = exp(b1*q0+b2+b3/(q0+0.2)+b4/pow((q0+0.2),2));
        // pi+omega->gamma+rho
        double FFompi = exp(d1*q0+d2+d3*pow(q0,d4)+d5*pow((q0+d6),d7));

        eqrate_ptr[i] = FFpiro + FFompi + FFomro;
    }
}
