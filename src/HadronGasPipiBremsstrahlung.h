#ifndef SRC_HADRONGASPIPIBREMSSTRAHLUNG
#define SRC_HADRONGASPIPIBREMSSTRAHLUNG

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalPhoton.h"

class HadronGasPipiBremsstrahlung : public ThermalPhoton {
  public:
    HadronGasPipiBremsstrahlung(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);
    ~HadronGasPipiBremsstrahlung() {}
    void analyticRates(
        double T, std::vector<double> &Eq, std::vector<double> &eqrate_ptr);
};

#endif
