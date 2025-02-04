/////////////////////////////////////////////////////////////////////////
//  To do in the future:
//      change the integration routines into gaussian integration
/////////////////////////////////////////////////////////////////////////

#include "ThermalDilepton.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Arsenal.h"
#include "ParameterReader.h"
#include "Table2D.h"
#include "gauss_quadrature.h"

using namespace std;
using ARSENAL::createA2DMatrix;
using ARSENAL::createA3DMatrix;
using ARSENAL::createA4DMatrix;
using ARSENAL::deleteA2DMatrix;
using ARSENAL::deleteA3DMatrix;
using ARSENAL::deleteA4DMatrix;

ThermalDilepton::ThermalDilepton(
    std::shared_ptr<ParameterReader> paraRdr_in, std::string emissionProcess) {
    paraRdr = paraRdr_in;
    emissionProcess_name = emissionProcess;

    neta = paraRdr->getVal("neta");
    np = paraRdr->getVal("np");
    nphi = paraRdr->getVal("nphi");
    nrapidity = paraRdr->getVal("nrapidity");
    nMInv_ = paraRdr->getVal("nMInv");
    norder_ = paraRdr->getVal("norder");
    rate_path_ = "ph_rates/";

    bRateTable_ = false;
    bShearVisCorr_ = false;
    bBulkVisCorr_ = false;

    // initial variables for photon spectra
    double p_i = paraRdr->getVal("photon_q_i");
    double p_f = paraRdr->getVal("photon_q_f");
    double phi_i = paraRdr->getVal("photon_phi_q_i");
    double phi_f = paraRdr->getVal("photon_phi_q_f");
    double y_i = paraRdr->getVal("photon_y_i");
    double y_f = paraRdr->getVal("photon_y_f");
    if (nrapidity > 1) {
        dy = (y_f - y_i) / (nrapidity - 1 + 1e-100);
    } else {
        dy = 1.0;
    }
    double m_i = paraRdr->getVal("dilepton_mass_i");
    double m_f = paraRdr->getVal("dilepton_mass_f");

    p = new double[np];
    p_weight = new double[np];
    phi = new double[nphi];
    phi_weight = new double[nphi];

    gauss_quadrature(np, 1, 0.0, 0.0, p_i, p_f, p, p_weight);
    gauss_quadrature(nphi, 1, 0.0, 0.0, phi_i, phi_f, phi, phi_weight);

    y.resize(nrapidity, 0);
    theta.resize(nrapidity, 0);
    for (int i = 0; i < nrapidity; i++) {
        y[i] = y_i + i * dy;
        theta[i] = acos(tanh(y[i]));  // rapidity's corresponding polar angle
    }

    if (nMInv_ > 1) {
        dMInv_ = (m_f - m_i) / (nMInv_ - 1 + 1e-100);
    } else {
        dMInv_ = 1.0;
    }
    MInv_.resize(nMInv_, 0);
    for (int i = 0; i < nMInv_; i++) {
        MInv_[i] = m_i + i * dMInv_;
    }

    dNpTdpTdphidydM_eq = createA4DMatrix(nMinv_, np, nphi, nrapidity, 0.);

    vnMInvpT_cos_eq = createA3DMatrix(norder_, nMInv_, np, 0.);
    vnMInvpT_sin_eq = createA3DMatrix(norder_, nMInv_, np, 0.);

    vnMInv_cos_eq = createA2DMatrix(norder_, nMInv_, 0.);
    vnMInv_sin_eq = createA2DMatrix(norder_, nMInv_, 0.);
}

ThermalDilepton::~ThermalDilepton() {
    if (bRateTable_) {
        int TbsizeX = EmissionrateTb_sizeX;
        deleteA3DMatrix(Emission_eqrateTb_ptr, TbsizeX);
    }

    delete[] p;
    delete[] p_weight;
    delete[] phi;
    delete[] phi_weight;

    deleteA4DMatrix(dNpTdpTdphidydM_eq, nMInv_, np, nphi);
    deleteA3DMatrix(vnMInvpT_cos_eq, norder_, nMInv_);
    deleteA3DMatrix(vnMInvpT_sin_eq, norder_, nMInv_);
    deleteA2DMatrix(vnMInv_cos_eq, norder_);
    deleteA2DMatrix(vnMInv_sin_eq, norder_);
}


void ThermalDilepton::analyticRates(
    const double T, const double MInv, std::vector<double> &Eq,
    std::vector<double> &eqrate_ptr) {
    for (unsigned int i = 0; i < Eq.size(); i++) {
        eqrate_ptr[i] = 1e-16;
    }
}

void ThermalDilepton::analyticRatesShearVis(
    double T, vector<double> &Eq, std::vector<double> &visrate_ptr) {
    for (unsigned int i = 0; i < visrate_ptr.size(); i++) {
        visrate_ptr[i] = 0.;
    }
}

void ThermalDilepton::analyticRatesBulkVis(
    double T, vector<double> &Eq, std::vector<double> &bulkvis_ptr) {
    for (unsigned int i = 0; i < bulkvis_ptr.size(); i++) {
        bulkvis_ptr[i] = 0.;
    }
}

void ThermalDilepton::checkAnalyticRates() {
    ofstream checkRates("checkPhotonRates.dat");
    double Emin = 0.05;
    double Tmin = 0.1;
    double dE = 0.05;
    double dT = 0.002;
    int nE = 80;
    int nT = 351;
    vector<double> Eq(nE, 0);
    for (int iE = 0; iE < nE; iE++) {
        Eq[iE] = Emin + iE * dE;
    }
    vector<double> eqrate(nE, 0);
    for (int iT = 0; iT < nT; iT++) {
        double T_local = Tmin + iT * dT;
        analyticRates(T_local, 0, Eq, eqrate);
        for (const auto rate_i : eqrate) {
            checkRates << std::scientific << std::setprecision(6)
                       << std::setw(10) << rate_i << "  ";
        }
        checkRates << std::endl;
    }
    checkRates.close();
}

void ThermalDilepton::getEmissionRate(
    vector<double> &Eq, const double Minv, const double T, const double muB,
    vector<double> &eqrate_ptr) {
    if (bRateTable_) {
        analyticRates(T, Minv, Eq, eqrate_ptr);
    } else {
        analyticRates(T, Minv, Eq, eqrate_ptr);
    }
    NetBaryonCorrection(T, muB, Eq, eqrate_ptr);
}

void ThermalDilepton::calThermalDileptonemission_3d(
    vector<double> &Eq, vector<double> &pi_zz, vector<double> &bulkPi,
    double T, double muB, double volume, double fraction) {
    const int Tb_length = Eq.size();
    if (Tb_length != nrapidity * nphi * np) {
        std::cout << "The length of Eq array is not right! Please check!"
                  << std::endl;
        exit(1);
    }
    const double volfrac = volume * fraction;
    // photon emission equilibrium rate at local rest cell
    vector<double> em_eqrate(Tb_length, 0);
    // photon emission viscous correction at local rest cell
    vector<double> em_visrate(Tb_length, 0);
    // photon emission bulk viscous correction at local rest cell
    vector<double> em_bulkvis(Tb_length, 0);

    getEmissionRate(Eq, T, muB, em_eqrate);

    int idx = 0;
    for (int k = 0; k < nrapidity; k++) {
        for (int m = 0; m < nphi; m++) {
            for (int l = 0; l < np; l++) {
                double local_eq = em_eqrate[idx];

                double temp_eq_sum = local_eq * volfrac;

                dNpTdpTdphidydM_eq[l][m][k] += temp_eq_sum;
                idx++;
            }
        }
    }
}

void ThermalDilepton::calThermalDileptonemission(
    vector<double> &Eq, vector<double> &pi_zz, vector<double> &bulkPi,
    int Tb_length, double T, vector<double> &volume, double fraction) {
    // photon emission equilibrium rate at local rest cell
    vector<double> em_eqrate(Tb_length, 0);

    getEmissionRate(Eq, T, 0., em_eqrate);

    int n_pt_point = nrapidity * np * nphi;

    double temp_eq_sum;
    int idx = 0;
    for (int k = 0; k < nrapidity; k++) {
        for (int m = 0; m < nphi; m++) {
            for (int l = 0; l < np; l++) {
                temp_eq_sum = 0.0;
                for (int i = 0; i < neta; i++) {
                    double local_eq = em_eqrate[idx + i * n_pt_point];
                    temp_eq_sum += local_eq * volume[i] * fraction;
                }
                dNpTdpTdphidydM_eq[l][m][k] += temp_eq_sum;
                idx++;
            }
        }
    }
}


void ThermalDilepton::calPhoton_SpvnpT(
    double ***dNd2pTdphidy, double ***vnypT_cos, double ***vnypT_sin,
    double **vnpT_cos, double **vnpT_sin, vector<double> &vn_cos,
    vector<double> &vn_sin) {
    // calculate the photon spectra and differential vn
    const double eps = 1e-15;
    for (int i = 0; i < np; i++) {
        for (int k = 0; k < nrapidity; k++) {
            for (int j = 0; j < nphi; j++) {
                double weight = phi_weight[j];
                for (int order = 0; order < norder; order++) {
                    double cos_tmp =
                        (dNd2pTdphidy[i][j][k] * cos(order * phi[j]) * weight);
                    double sin_tmp =
                        (dNd2pTdphidy[i][j][k] * sin(order * phi[j]) * weight);
                    vnypT_cos[order][i][k] += cos_tmp;
                    vnypT_sin[order][i][k] += sin_tmp;
                    if (std::abs(y[k]) < 0.5) {
                        // only integrate mid-rapidity
                        vnpT_cos[order][i] += cos_tmp * dy;
                        vnpT_sin[order][i] += sin_tmp * dy;
                    }
                }
            }
        }
        double p_weight_factor = p[i] * p_weight[i];
        for (int order = 0; order < norder; order++) {
            vn_cos[order] += vnpT_cos[order][i] * p_weight_factor;
            vn_sin[order] += vnpT_sin[order][i] * p_weight_factor;

            // vn(pT)
            if (order > 0) {
                vnpT_cos[order][i] /= (vnpT_cos[0][i] + eps);
                vnpT_sin[order][i] /= (vnpT_cos[0][i] + eps);
                for (int k = 0; k < nrapidity; k++) {
                    vnypT_cos[order][i][k] /= (vnypT_cos[0][i][k] + eps);
                    vnypT_sin[order][i][k] /= (vnypT_cos[0][i][k] + eps);
                }
            }
        }
        vnpT_cos[0][i] /= (2 * M_PI);  // dN/(2pi dy pT dpT)
    }
    for (int order = 1; order < norder; order++) {
        // vn
        vn_cos[order] /= (vn_cos[0] + eps);
        vn_sin[order] /= (vn_cos[0] + eps);
    }
}

void ThermalDilepton::calPhoton_SpvnpT_shell() {
    calPhoton_SpvnpT(
        dNd2pTdphidy_eq, vnypT_cos_eq, vnypT_sin_eq, vnpT_cos_eq, vnpT_sin_eq,
        vn_cos_eq, vn_sin_eq);
}


void ThermalDilepton::outputPhoton_SpvnpT(
    string path, string type_str, double ***dNd2pTdphidy, double ***vnypT_cos,
    double ***vnypT_sin, double **vnpT_cos, double **vnpT_sin,
    vector<double> &vn_cos, vector<double> &vn_sin) {
    ostringstream filename_stream_SpMatrix;
    ostringstream filename_stream_Spvn;
    ostringstream filename_stream_inte_Spvn;

    filename_stream_SpMatrix << path << emissionProcess_name << "_Spvn_"
                             << type_str << "_ypTdiff.dat";
    filename_stream_Spvn << path << emissionProcess_name << "_Spvn_" << type_str
                         << ".dat";
    filename_stream_inte_Spvn << path << emissionProcess_name << "_Spvn_"
                              << type_str << "_inte.dat";

    ofstream fphotonSpMatrix(filename_stream_SpMatrix.str().c_str());
    ofstream fphotonSpvn(filename_stream_Spvn.str().c_str());
    ofstream fphotoninteSpvn(filename_stream_inte_Spvn.str().c_str());

    for (int k = 0; k < nrapidity; k++) {
        for (int j = 0; j < np; j++) {
            fphotonSpMatrix << scientific << setprecision(6) << setw(16) << y[k]
                            << "  " << p[j] << "  " << vnypT_cos[0][j][k]
                            << "  ";
            for (int order = 1; order < norder; order++) {
                fphotonSpMatrix << scientific << setprecision(6) << setw(16)
                                << vnypT_cos[order][j][k] << "  "
                                << vnypT_sin[order][j][k] << "  ";
            }
            fphotonSpMatrix << endl;
        }
    }

    for (int i = 0; i < np; i++) {
        fphotonSpvn << scientific << setprecision(6) << setw(16) << p[i] << "  "
                    << vnpT_cos[0][i] << "  ";
        for (int order = 1; order < norder; order++) {
            fphotonSpvn << scientific << setprecision(6) << setw(16)
                        << vnpT_cos[order][i] << "  " << vnpT_sin[order][i]
                        << "  ";
        }
        fphotonSpvn << endl;
    }

    for (int order = 0; order < norder; order++) {
        fphotoninteSpvn << scientific << setprecision(6) << setw(16) << order
                        << "   " << vn_cos[order] << "   " << vn_sin[order]
                        << endl;
    }

    fphotonSpMatrix.close();
    fphotonSpvn.close();
    fphotoninteSpvn.close();
}

void ThermalDilepton::outputPhoton_SpvnpT_shell(string path) {
    outputPhoton_SpvnpT(
        path, "eq", dNd2pTdphidy_eq, vnypT_cos_eq, vnypT_sin_eq, vnpT_cos_eq,
        vnpT_sin_eq, vn_cos_eq, vn_sin_eq);
}
