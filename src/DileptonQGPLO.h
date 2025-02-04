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

    void analyticRates(
        const double T, const double MInv, const double Eq, double &eqrate);
};

#endif  // SRC_DILEPTONQGPLO_H
