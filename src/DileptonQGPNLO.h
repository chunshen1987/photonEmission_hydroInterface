#ifndef SRC_DILEPTONQGPNLO_H
#define SRC_DILEPTONQGPNLO_H

#include <memory>
#include <string>
#include <vector>

#include "ParameterReader.h"
#include "ThermalDilepton.h"

class DileptonQGPNLO : public ThermalDilepton {
  public:
    DileptonQGPNLO(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);
    ~DileptonQGPNLO();

    double nF(double x);
    double nB(double x);

    double l1f(double x);
    double l2f(double x);
    double l3f(double x);

    void rho_LO(double o, double k, double mu, double &rT, double &rL);
    double rhoV_AMY(double alpha, double mu, double k);
    void rho_OPE(
        double o, double k, double alpha, double mu, double &rT, double &rL);

    int getIdx(double xval, std::vector<double> &xTable);
    void interp(
        const double alpha, const double muOverT, const double MOverT,
        const double kOverT, double &resRhoT, double &resRhoL);
    int approx_rho(double *input, double &rT, double &rL);
    double mllFactor(double x);

    void readInEmissionTables(std::string emissionProcess);
    void getRateFromTable(
        const double E, const double k, const double alpha_s, const double muB,
        const double T, const double m_l, double &rateTot, double &rateT,
        double &rateL);

  private:
    const double OOFP = 0.0795774715459476678844418816863;  // 1/(4*pi)
    const int Nc = 3;

    std::string ratePath_;

    double ****rateRhoT;
    double ****rateRhoL;

    std::vector<double> alphaS_list;
    std::vector<double> muoverT_list;
    std::vector<double> MoverT_list;
    std::vector<double> koverT_list;
};

#endif  // SRC_DILEPTONQGPNLO_H
