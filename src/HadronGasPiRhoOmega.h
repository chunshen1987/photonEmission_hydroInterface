#ifndef SRC_HADRONGASPIRHOOMEGA
#define SRC_HADRONGASPIRHOOMEGA

#include <memory>
#include <vector>

#include "ParameterReader.h"
#include "ThermalPhoton.h"

class HadronGasPiRhoOmega : public ThermalPhoton {
  public:
    HadronGasPiRhoOmega(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);
    ~HadronGasPiRhoOmega() {}
    void analyticRates(
        double T, std::vector<double> &Eq, std::vector<double> &eqrate_ptr);
};

#endif
