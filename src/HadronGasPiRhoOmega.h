#ifndef SRC_HADRONGASPIRHOOMEGA
#define SRC_HADRONGASPIRHOOMEGA

#include <vector>
#include <memory>

#include "ThermalPhoton.h"
#include "ParameterReader.h"

class HadronGasPiRhoOmega : public ThermalPhoton {
 public:
    HadronGasPiRhoOmega(std::shared_ptr<ParameterReader> paraRdr_in,
                        std::string emissionProcess);
    ~HadronGasPiRhoOmega() {}
    void analyticRates(double T, std::vector<double> &Eq,
                       std::vector<double> &eqrate_ptr);
};

#endif
