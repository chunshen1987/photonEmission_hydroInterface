#ifndef SRC_DILEPTONQGPLO_H
#define SRC_DILEPTONQGPLO_H

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalDilepton.h"

class DileptonQGPLO : public ThermalDilepton {
  public:
    DileptonQGPLO(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);
    ~DileptonQGPLO() {}

    double nB(double x);
    double l1f(double x);
    double l2f(double x);
    double l3f(double x);

    void rho_LO(double o, double k, double mu, double &rT, double &rL);
    double mllFactor(double x);

    void analyticRates(
        const double E, const double k, const double muB, const double T,
        const double m_l, double &rateTot, double &rateT, double &rateL);

  private:
    const double OOFP = 0.0795774715459476678844418816863;  // 1/(4*pi)
};

#endif  // SRC_DILEPTONQGPLO_H
