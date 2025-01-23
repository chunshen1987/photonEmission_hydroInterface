#ifndef SRC_HADRONGASRHOSPECTRALFUNCTION
#define SRC_HADRONGASRHOSPECTRALFUNCTION

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalPhoton.h"

class HadronGasRhoSpectralFunction : public ThermalPhoton {
  public:
    HadronGasRhoSpectralFunction(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);
    ~HadronGasRhoSpectralFunction() {}
    void analyticRates(
        double T, std::vector<double> &Eq, std::vector<double> &eqrate_ptr);
    void NetBaryonCorrection(
        double T, double muB, std::vector<double> &Eq,
        std::vector<double> &eqrate_ptr);
};

#endif
