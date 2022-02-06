#ifndef SRC_THERMALPHOTON_H_
#define SRC_THERMALPHOTON_H_

#include <string>
#include <fstream>
#include <memory>
#include <vector>

#include "Arsenal.h"
#include "Table2D.h"
#include "ParameterReader.h"

class ThermalPhoton {
 private:
    std::shared_ptr<ParameterReader> paraRdr;

    int np, nphi, nrapidity;
    int norder;
    int neta;
    std::string rate_path_;

    double dy;

    bool bRateTable_;
    bool bShearVisCorr_;
    bool bBulkVisCorr_;

    // photon emission rate
    std::unique_ptr<Table2D> Photonemission_eqrateTable_ptr;
    std::unique_ptr<Table2D> Photonemission_viscous_rateTable_ptr;
    std::unique_ptr<Table2D> Photonemission_bulkvis_rateTable_ptr;

    double** Emission_eqrateTb_ptr;
    double** Emission_viscous_rateTb_ptr;
    double** Emission_bulkvis_rateTb_ptr;
    std::vector<double> EmissionrateTb_Yidxptr;
    double EmissionrateTb_Xmin;
    double EmissionrateTb_Ymin;
    int EmissionrateTb_sizeX;
    int EmissionrateTb_sizeY;
    double EmissionrateTb_dX;
    double EmissionrateTb_dY;

    // photon spectra parameters
    std::string emissionProcess_name;
    double *p, *p_weight;
    double *phi, *phi_weight;
    std::vector<double> y;
    std::vector<double> theta;

    double ***dNd2pTdphidy_eq, ***dNd2pTdphidy_vis, ***dNd2pTdphidy_tot;
    double ***dNd2pTdphidy_bulkvis;
    double ***dNd2pTdphidy_vis_deltaf_restricted;
    double ***dNd2pTdphidy_bulkvis_deltaf_restricted;

    double **vnpT_cos_eq, **vnpT_sin_eq;
    double ***vnypT_cos_eq, ***vnypT_sin_eq;
    double **vnpT_cos_vis, **vnpT_sin_vis;
    double ***vnypT_cos_vis, ***vnypT_sin_vis;
    double **vnpT_cos_vis_deltaf_restricted;
    double **vnpT_sin_vis_deltaf_restricted;
    double ***vnypT_cos_vis_deltaf_restricted;
    double ***vnypT_sin_vis_deltaf_restricted;
    double **vnpT_cos_bulkvis, **vnpT_sin_bulkvis;
    double ***vnypT_cos_bulkvis, ***vnypT_sin_bulkvis;
    double **vnpT_cos_bulkvis_deltaf_restricted;
    double **vnpT_sin_bulkvis_deltaf_restricted;
    double ***vnypT_cos_bulkvis_deltaf_restricted;
    double ***vnypT_sin_bulkvis_deltaf_restricted;
    double **vnpT_cos_tot, **vnpT_sin_tot;
    double ***vnypT_cos_tot, ***vnypT_sin_tot;

    std::vector<double> vn_cos_eq;
    std::vector<double> vn_sin_eq;
    std::vector<double> vn_cos_vis;
    std::vector<double> vn_sin_vis;
    std::vector<double> vn_cos_vis_deltaf_restricted;
    std::vector<double> vn_sin_vis_deltaf_restricted;
    std::vector<double> vn_cos_bulkvis;
    std::vector<double> vn_sin_bulkvis;
    std::vector<double> vn_cos_bulkvis_deltaf_restricted;
    std::vector<double> vn_sin_bulkvis_deltaf_restricted;
    std::vector<double> vn_cos_tot;
    std::vector<double> vn_sin_tot;

    // matrix for cuts on temperature and proper time
    int nTcut, nTaucut;
    double Tcut_high, Tcut_low;
    double Taucut_high, Taucut_low;
    double *****dNd2pTdphidydTdtau_eq, *****dNd2pTdphidydTdtau_tot;
    double *****dNd2pTdphidydTdtau_vis, *****dNd2pTdphidydTdtau_bulkvis;
    double **dNdydTdtau_eq, **dNdydTdtau_tot;
    double **dNdydTdtau_vis, **dNdydTdtau_bulkvis;
    double ***vndTdtau_cos_eq, ***vndTdtau_sin_eq;
    double ***vndTdtau_cos_vis, ***vndTdtau_sin_vis;
    double ***vndTdtau_cos_bulkvis, ***vndTdtau_sin_bulkvis;
    double ***vndTdtau_cos_tot, ***vndTdtau_sin_tot;

    // matrix for cuts on x and proper time
    int n_xperp_cut, n_tau_cut_xtau;
    double xperp_high, xperp_low;
    double tau_cut_high, tau_cut_low;
    double *****dNd2pTdphidydxperpdtau_eq, *****dNd2pTdphidydxperpdtau_tot;
    double *****dNd2pTdphidydxperpdtau_vis;
    double *****dNd2pTdphidydxperpdtau_bulkvis;
    double **dNdydxperpdtau_eq, **dNdydxperpdtau_tot;
    double **dNdydxperpdtau_vis, **dNdydxperpdtau_bulkvis;
    double ***vndxperpdtau_cos_eq, ***vndxperpdtau_sin_eq;
    double ***vndxperpdtau_cos_vis, ***vndxperpdtau_sin_vis;
    double ***vndxperpdtau_cos_bulkvis, ***vndxperpdtau_sin_bulkvis;
    double ***vndxperpdtau_cos_tot, ***vndxperpdtau_sin_tot;

 public:
    ThermalPhoton(std::shared_ptr<ParameterReader> paraRdr_in,
                  std::string emissionProcess);

    virtual ~ThermalPhoton();

    void setupEmissionrateFromFile(
            double Xmin, double dX, double Ymin, double dY,
            bool bShearVisCorr, bool bBulkVisCorr);
    void readEmissionrate(std::string);

    void setupEmissionrateFromParametrization(
        double Xmin, double dX, int nX,
        double Ymin, double dY, int nY);

    double get_dy() {return(dy);}

    double getPhotonp(int i) {return(p[i]);}
    double getPhoton_pweight(int i) {return(p_weight[i]);}
    double getPhotonphi(int i) {return(phi[i]);}
    double getPhoton_phiweight(int i) {return(phi_weight[i]);}
    double getPhotontheta(int i) {return(theta[i]);}
    double getPhotonrapidity(int i) {return(y[i]);}
    double getPhotonSpMatrix_eq(int i, int j, int k) {
        return(dNd2pTdphidy_eq[i][j][k]);
    }
    double getPhotonSpMatrix_tot(int i, int j, int k) {
        return(dNd2pTdphidy_tot[i][j][k]);
    }

    virtual void analyticRates(double T, std::vector<double> &Eq,
                               std::vector<double> &eqrate_ptr);
    virtual void NetBaryonCorrection(double T, double muB,
                                     std::vector<double> &Eq,
                                     std::vector<double> &eqrate_ptr) {}

    virtual void analyticRatesShearVis(double T, std::vector<double> &Eq,
                                       std::vector<double> &eqrate_ptr);
    virtual void analyticRatesBulkVis(double T, std::vector<double> &Eq,
                                      std::vector<double> &eqrate_ptr);

    void getPhotonemissionRate(std::vector<double> &Eq,
                               std::vector<double> &pi_zz,
                               std::vector<double> &bulkPi,
                               const double T, const double muB,
                               std::vector<double> &eqrate_ptr,
                               std::vector<double> &visrate_ptr,
                               std::vector<double> &bulkvis_ptr);
    void calThermalPhotonemission(
        std::vector<double> &Eq, std::vector<double> &pi_zz,
        std::vector<double> &bulkPi, int Tb_length, double T,
        std::vector<double> &volume, double fraction);
    void calThermalPhotonemission_3d(
        std::vector<double> &Eq, std::vector<double> &pi_zz,
        std::vector<double> &bulkPi, double T, double muB, double volume,
        double fraction);
    void calThermalPhotonemissiondTdtau(
        std::vector<double> &Eq, std::vector<double> &pi_zz,
        std::vector<double> &bulkPi, int Tb_length, double T, double tau,
        std::vector<double> &volume, double fraction);
    void calThermalPhotonemissiondTdtau_3d(
        std::vector<double> &Eq, std::vector<double> &pi_zz,
        std::vector<double> &bulkPi,
        double T, double muB, double tau, double volume, double fraction);
    void calThermalPhotonemissiondxperpdtau(
        std::vector<double> &Eq, std::vector<double> &pi_zz,
        std::vector<double> &bulkPi, int Tb_length,
        double T, double x_local, double tau,
        std::vector<double> &volume, double fraction);
    void calThermalPhotonemissiondxperpdtau_3d(
        std::vector<double> &Eq, std::vector<double> &pi_zz,
        std::vector<double> &bulkPi, double T, double muB,
        double x_local, double tau, double volume, double fraction);

    void calPhoton_SpvnpT(double ***dNd2pTdphipy,
                          double ***vnypT_cos, double *** vnypT_sin,
                          double **vnpT_cos, double **vnpT_sin,
                          std::vector<double> &vn_cos,
                          std::vector<double> &vn_sin);
    void calPhoton_SpvnpT_shell();
    void calPhoton_SpvnpT_dTdtau();
    void calPhoton_SpvnpT_dxperpdtau();
    void outputPhoton_SpvnpT(std::string path, std::string type_str,
                             double ***dNd2pTdphidy,
                             double ***vnypT_cos, double ***vnypT_sin,
                             double **vnpT_cos, double **vnpT_sin,
                             std::vector<double> &vn_cos,
                             std::vector<double> &vn_sin);
    void outputPhoton_SpvnpT_shell(std::string path);
    void outputPhoton_SpvnpTdTdtau(std::string path);
    void outputPhoton_SpvnpTdxperpdtau(std::string path);
    void output_photon_spectra_dTdtau(std::string path);
    void interpolation2D_bilinear(double varX, std::vector<double> &varY,
                                  double** Table2D_ptr,
                                  std::vector<double> &results);

    void update_rates_with_polyakov_suppression();
    double get_polyakov_suppression_factor(double T_in_GeV);
};
#endif  // SRC_THERMALPHOTON_H_
