#ifndef SRC_DILEPTONQGPLO_H
#define SRC_DILEPTONQGPLO_H

#include <vector>
#include <memory>

#include "ThermalDilepton.h"
#include "ParameterReader.h"

class DileptonQGPLO : public ThermalDilepton {
 public:
    DileptonQGPLO(std::shared_ptr<ParameterReader> paraRdr_in,
                  std::string emissionProcess);
    ~DileptonQGPLO() {}

    void analyticRates(
        double T, std::vector<double> &Eq, std::vector<double> &eqrate_ptr);
};

#endif  // SRC_DILEPTONQGPLO_H
