#ifndef SRC_THERMALDILEPTON_H_
#define SRC_THERMALDILEPTON_H_

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Arsenal.h"
#include "ParameterReader.h"
#include "Table2D.h"

class ThermalDilepton {
  private:
    std::shared_ptr<ParameterReader> paraRdr;

    int np, nphi, nrapidity, nMInv_;
    int norder_;
    int neta;
    std::string rate_path_;

    double mlsq_;

    double dy, dMInv_;

    bool bRateTable_;
    bool bShearVisCorr_;
    bool bBulkVisCorr_;

    double alphaS_;

    // photon spectra parameters
    std::string emissionProcess_name;
    double *p, *p_weight;
    double *phi, *phi_weight;
    std::vector<double> y;
    std::vector<double> theta;
    std::vector<double> MInv_;

    double ****dNpTdpTdphidydM_eq;
    double ****dNpTdpTdphidydM_eqT;
    double ****dNpTdpTdphidydM_eqL;

    double **vnMInv_cos_eq, **vnMInv_sin_eq;
    double **vnMInv_cos_eqT, **vnMInv_sin_eqT;
    double **vnMInv_cos_eqL, **vnMInv_sin_eqL;
    double ***vnMInvpT_cos_eq, ***vnMInvpT_sin_eq;
    double ***vnMInvpT_cos_eqT, ***vnMInvpT_sin_eqT;
    double ***vnMInvpT_cos_eqL, ***vnMInvpT_sin_eqL;

  public:
    ThermalDilepton(
        std::shared_ptr<ParameterReader> paraRdr_in,
        std::string emissionProcess);

    virtual ~ThermalDilepton();

    void setRateTableFlag(bool bRateTable) { bRateTable_ = bRateTable; }
    bool getRateTableFlag() const { return (bRateTable_); }
    double get_dy() const { return (dy); }
    double get_dMInv() const { return (dMInv_); }

    double getPhotonp(int i) { return (p[i]); }
    double getPhoton_pweight(int i) { return (p_weight[i]); }
    double getPhotonphi(int i) { return (phi[i]); }
    double getPhoton_phiweight(int i) { return (phi_weight[i]); }
    double getPhotontheta(int i) { return (theta[i]); }
    double getPhotonrapidity(int i) { return (y[i]); }

    double getDileptonMinv(int i) { return (MInv_[i]); }

    virtual void analyticRates(
        const double E, const double k, const double muB, const double T,
        const double m_l, double &rateTot, double &rateT, double &rateL) {
        rateTot = 0.0;
        rateT = 0.0;
        rateL = 0.0;
    };

    virtual void NetBaryonCorrection(
        double T, double muB, std::vector<double> &Eq,
        std::vector<double> &eqrate_ptr) {}

    virtual void getRateFromTable(
        const double E, const double k, const double MInv, const double alpha_s,
        const double muB, const double T, const double m_l, double &rateTot,
        double &rateT, double &rateL) {
        rateTot = 0.0;
        rateT = 0.0;
        rateL = 0.0;
    };

    // void checkAnalyticRates();

    void getEmissionRate(
        std::vector<double> &Eq, const double T, const double muB,
        std::vector<double> &eqrate_ptr, std::vector<double> &eqrateT_ptr,
        std::vector<double> &eqrateL_ptr);
    void calThermalDileptonemission(
        std::vector<double> &Eq, int Tb_length, double T,
        std::vector<double> &volume, double fraction);
    void calThermalDileptonemission_3d(
        std::vector<double> &Eq, double T, double muB, double volume,
        double fraction);

    void calPhoton_SpvnpT(
        double ****dNpTdpTdphidydM, double ***vnMInvpT_cos,
        double ***vnMInvpT_sin, double **vnMInv_cos, double **vnMInv_sin);
    void calPhoton_SpvnpT_shell();
    void outputPhoton_SpvnpT(
        std::string path, std::string type_str, double ****dNpTdpTdphidydM,
        double ***vnMInvpT_cos, double ***vnMInvpT_sin, double **vnMInv_cos,
        double **vnMInv_sin);
    void outputPhoton_SpvnpT_shell(std::string path);
};

#endif  // SRC_THERMALDILEPTON_H
