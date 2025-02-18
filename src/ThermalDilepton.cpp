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
#include "data_struct.h"
#include "gauss_quadrature.h"

using namespace std;
using ARSENAL::createA2DMatrix;
using ARSENAL::createA3DMatrix;
using ARSENAL::createA4DMatrix;
using ARSENAL::deleteA2DMatrix;
using ARSENAL::deleteA3DMatrix;
using ARSENAL::deleteA4DMatrix;

using PhysConsts::MElectron;

ThermalDilepton::ThermalDilepton(
    std::shared_ptr<ParameterReader> paraRdr_in, std::string emissionProcess) {
    paraRdr = paraRdr_in;
    emissionProcess_name = emissionProcess;

    mlsq_ = MElectron * MElectron;

    neta = paraRdr->getVal("neta");
    np = paraRdr->getVal("np");
    nphi = paraRdr->getVal("nphi");
    nrapidity = paraRdr->getVal("nrapidity");
    nMInv_ = paraRdr->getVal("nMInv");
    norder_ = paraRdr->getVal("norder");
    rate_path_ = "ph_rates/";

    alphaS_ = paraRdr->getVal("alpha_s");

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

    dNpTdpTdphidydM_eq = createA4DMatrix(nMInv_, np, nphi, nrapidity, 0.);
    dNpTdpTdphidydM_eqT = createA4DMatrix(nMInv_, np, nphi, nrapidity, 0.);
    dNpTdpTdphidydM_eqL = createA4DMatrix(nMInv_, np, nphi, nrapidity, 0.);

    vnMInvpT_cos_eq = createA3DMatrix(norder_, nMInv_, np, 0.);
    vnMInvpT_sin_eq = createA3DMatrix(norder_, nMInv_, np, 0.);
    vnMInvpT_cos_eqT = createA3DMatrix(norder_, nMInv_, np, 0.);
    vnMInvpT_sin_eqT = createA3DMatrix(norder_, nMInv_, np, 0.);
    vnMInvpT_cos_eqL = createA3DMatrix(norder_, nMInv_, np, 0.);
    vnMInvpT_sin_eqL = createA3DMatrix(norder_, nMInv_, np, 0.);

    vnMInv_cos_eq = createA2DMatrix(norder_, nMInv_, 0.);
    vnMInv_sin_eq = createA2DMatrix(norder_, nMInv_, 0.);
    vnMInv_cos_eqT = createA2DMatrix(norder_, nMInv_, 0.);
    vnMInv_sin_eqT = createA2DMatrix(norder_, nMInv_, 0.);
    vnMInv_cos_eqL = createA2DMatrix(norder_, nMInv_, 0.);
    vnMInv_sin_eqL = createA2DMatrix(norder_, nMInv_, 0.);
}

ThermalDilepton::~ThermalDilepton() {
    delete[] p;
    delete[] p_weight;
    delete[] phi;
    delete[] phi_weight;

    deleteA4DMatrix(dNpTdpTdphidydM_eq, nMInv_, np, nphi);
    deleteA4DMatrix(dNpTdpTdphidydM_eqT, nMInv_, np, nphi);
    deleteA4DMatrix(dNpTdpTdphidydM_eqL, nMInv_, np, nphi);
    deleteA3DMatrix(vnMInvpT_cos_eq, norder_, nMInv_);
    deleteA3DMatrix(vnMInvpT_sin_eq, norder_, nMInv_);
    deleteA3DMatrix(vnMInvpT_cos_eqT, norder_, nMInv_);
    deleteA3DMatrix(vnMInvpT_sin_eqT, norder_, nMInv_);
    deleteA3DMatrix(vnMInvpT_cos_eqL, norder_, nMInv_);
    deleteA3DMatrix(vnMInvpT_sin_eqL, norder_, nMInv_);
    deleteA2DMatrix(vnMInv_cos_eq, norder_);
    deleteA2DMatrix(vnMInv_sin_eq, norder_);
    deleteA2DMatrix(vnMInv_cos_eqT, norder_);
    deleteA2DMatrix(vnMInv_sin_eqT, norder_);
    deleteA2DMatrix(vnMInv_cos_eqL, norder_);
    deleteA2DMatrix(vnMInv_sin_eqL, norder_);
}

// void ThermalDilepton::checkAnalyticRates() {
//     ofstream checkRates("checkPhotonRates.dat");
//     double Emin = 0.05;
//     double Tmin = 0.1;
//     double dE = 0.05;
//     double dT = 0.002;
//     int nE = 80;
//     int nT = 351;
//     vector<double> Eq(nE, 0);
//     for (int iE = 0; iE < nE; iE++) {
//         Eq[iE] = Emin + iE * dE;
//     }
//     vector<double> eqrate(nE, 0);
//     for (int iT = 0; iT < nT; iT++) {
//         double T_local = Tmin + iT * dT;
//         for (int iE = 0; iE < nE; iE++) {
//             int iM = static_cast<int>(iE / npoints) % nMInv_;
//             double MInv = MInv_[iM];
//             analyticRates(T_local, 1.0, Eq[iE], eqrate[iE]);
//         }
//         for (const auto rate_i : eqrate) {
//             checkRates << std::scientific << std::setprecision(6)
//                        << std::setw(10) << rate_i << "  ";
//         }
//         checkRates << std::endl;
//     }
//     checkRates.close();
// }

void ThermalDilepton::getEmissionRate(
    vector<double> &Eq, const double T, const double muB,
    vector<double> &eqrate_ptr, vector<double> &eqrateT_ptr,
    vector<double> &eqrateL_ptr) {
    int npoints = np * nphi * nrapidity;
    if (!bRateTable_) {
        for (int i = 0; i < Eq.size(); i++) {
            int iM = static_cast<int>(i / npoints) % nMInv_;
            double MInv = MInv_[iM];
            double k = sqrt(Eq[i] * Eq[i] - MInv * MInv);
            double rateTot, rateT, rateL;
            analyticRates(Eq[i], k, muB, T, MElectron, rateTot, rateT, rateL);
            eqrate_ptr[i] = rateTot;
            eqrateT_ptr[i] = rateT;
            eqrateL_ptr[i] = rateL;
        }
        NetBaryonCorrection(T, muB, Eq, eqrate_ptr);
    } else {
        for (int i = 0; i < Eq.size(); i++) {
            int iM = static_cast<int>(i / npoints) % nMInv_;
            double MInv = MInv_[iM];
            double k = sqrt(Eq[i] * Eq[i] - MInv * MInv);
            double rateTot, rateT, rateL;
            getRateFromTable(
                Eq[i], k, MInv, alphaS_, muB, T, mlsq_, rateTot, rateT, rateL);
            eqrate_ptr[i] = rateTot;
            eqrateT_ptr[i] = rateT;
            eqrateL_ptr[i] = rateL;
        }
    }
}

void ThermalDilepton::calThermalDileptonemission_3d(
    vector<double> &Eq, double T, double muB, double volume, double fraction) {
    const int Tb_length = Eq.size();
    if (Tb_length != nMInv_ * nrapidity * nphi * np) {
        std::cout << "The length of Eq array is not right! Please check!"
                  << std::endl;
        exit(1);
    }
    const double volfrac = volume * fraction;
    // photon emission equilibrium rate at local rest cell
    vector<double> em_eqrate(Tb_length, 0);
    vector<double> em_eqrateT(Tb_length, 0);
    vector<double> em_eqrateL(Tb_length, 0);

    getEmissionRate(Eq, T, muB, em_eqrate, em_eqrateT, em_eqrateL);

    int idx = 0;
    for (int im = 0; im < nMInv_; im++) {
        for (int k = 0; k < nrapidity; k++) {
            for (int m = 0; m < nphi; m++) {
                for (int l = 0; l < np; l++) {
                    dNpTdpTdphidydM_eq[im][l][m][k] += em_eqrate[idx] * volfrac;
                    dNpTdpTdphidydM_eqT[im][l][m][k] +=
                        em_eqrateT[idx] * volfrac;
                    dNpTdpTdphidydM_eqL[im][l][m][k] +=
                        em_eqrateL[idx] * volfrac;
                    idx++;
                }
            }
        }
    }
}

void ThermalDilepton::calThermalDileptonemission(
    vector<double> &Eq, int Tb_length, double T, vector<double> &volume,
    double fraction) {
    // photon emission equilibrium rate at local rest cell
    vector<double> em_eqrate(Tb_length, 0);
    vector<double> em_eqrateT(Tb_length, 0);
    vector<double> em_eqrateL(Tb_length, 0);

    getEmissionRate(Eq, T, 0., em_eqrate, em_eqrateT, em_eqrateL);

    int n_pt_point = nMInv_ * nrapidity * np * nphi;

    double temp_eq_sum, temp_eq_sumT, temp_eq_sumL;
    int idx = 0;
    for (int im = 0; im < nMInv_; im++) {
        for (int k = 0; k < nrapidity; k++) {
            for (int m = 0; m < nphi; m++) {
                for (int l = 0; l < np; l++) {
                    temp_eq_sum = 0.0;
                    temp_eq_sumT = 0.0;
                    temp_eq_sumL = 0.0;
                    for (int i = 0; i < neta; i++) {
                        // integrate over eta
                        double local_eq = em_eqrate[idx + i * n_pt_point];
                        temp_eq_sum += local_eq * volume[i] * fraction;

                        double local_eqT = em_eqrateT[idx + i * n_pt_point];
                        temp_eq_sumT += local_eqT * volume[i] * fraction;

                        double local_eqL = em_eqrateL[idx + i * n_pt_point];
                        temp_eq_sumL += local_eqL * volume[i] * fraction;
                    }
                    dNpTdpTdphidydM_eq[im][l][m][k] += temp_eq_sum;
                    dNpTdpTdphidydM_eqT[im][l][m][k] += temp_eq_sumT;
                    dNpTdpTdphidydM_eqL[im][l][m][k] += temp_eq_sumL;
                    idx++;
                }
            }
        }
    }
}

void ThermalDilepton::calPhoton_SpvnpT(
    double ****dNpTdpTdphidydM, double ***vnMInvpT_cos, double ***vnMInvpT_sin,
    double **vnMInv_cos, double **vnMInv_sin) {
    // calculate the photon spectra and differential vn
    const double eps = 1e-15;
    for (int im = 0; im < nMInv_; im++) {
        for (int i = 0; i < np; i++) {
            for (int k = 0; k < nrapidity; k++) {
                for (int j = 0; j < nphi; j++) {
                    double weight = phi_weight[j];
                    for (int order = 0; order < norder_; order++) {
                        double cos_tmp =
                            (dNpTdpTdphidydM[im][i][j][k] * cos(order * phi[j])
                             * weight);
                        double sin_tmp =
                            (dNpTdpTdphidydM[im][i][j][k] * sin(order * phi[j])
                             * weight);
                        if (std::abs(y[k]) < 0.5) {
                            vnMInvpT_cos[order][im][i] += cos_tmp * dy;
                            vnMInvpT_sin[order][im][i] += sin_tmp * dy;
                        }
                    }
                }
            }
            double p_weight_factor = p[i] * p_weight[i];
            for (int order = 0; order < norder_; order++) {
                vnMInv_cos[order][im] +=
                    vnMInvpT_cos[order][im][i] * p_weight_factor;
                vnMInv_sin[order][im] +=
                    vnMInvpT_sin[order][im][i] * p_weight_factor;

                // vn(pT)
                if (order > 0) {
                    vnMInvpT_cos[order][im][i] /=
                        (vnMInvpT_cos[0][im][i] + eps);
                    vnMInvpT_sin[order][im][i] /=
                        (vnMInvpT_cos[0][im][i] + eps);
                }
            }
            vnMInvpT_cos[0][im][i] /= (2 * M_PI);  // dN/(2pi dy pT dpT)
        }
        for (int order = 1; order < norder_; order++) {
            // vn
            vnMInv_cos[order][im] /= (vnMInv_cos[0][im] + eps);
            vnMInv_sin[order][im] /= (vnMInv_cos[0][im] + eps);
        }
    }
}

void ThermalDilepton::calPhoton_SpvnpT_shell() {
    calPhoton_SpvnpT(
        dNpTdpTdphidydM_eq, vnMInvpT_cos_eq, vnMInvpT_sin_eq, vnMInv_cos_eq,
        vnMInv_sin_eq);
    calPhoton_SpvnpT(
        dNpTdpTdphidydM_eqT, vnMInvpT_cos_eqT, vnMInvpT_sin_eqT, vnMInv_cos_eqT,
        vnMInv_sin_eqT);
    calPhoton_SpvnpT(
        dNpTdpTdphidydM_eqL, vnMInvpT_cos_eqL, vnMInvpT_sin_eqL, vnMInv_cos_eqL,
        vnMInv_sin_eqL);
}

void ThermalDilepton::outputPhoton_SpvnpT(
    string path, string type_str, double ****dNpTdpTdphidydM,
    double ***vnMInvpT_cos, double ***vnMInvpT_sin, double **vnMInv_cos,
    double **vnMInv_sin) {
    ostringstream filename_stream_SpMatrix;
    ostringstream filename_stream_Spvn;
    ostringstream filename_stream_inte_Spvn;

    filename_stream_SpMatrix << path << emissionProcess_name << "_Spvn_"
                             << type_str << "_MInvpTdiff.dat";
    filename_stream_Spvn << path << emissionProcess_name << "_Spvn_" << type_str
                         << "_MInv.dat";

    ofstream fphotonSpMatrix(filename_stream_SpMatrix.str().c_str());
    ofstream fphotonSpvn(filename_stream_Spvn.str().c_str());

    for (int im = 0; im < nMInv_; im++) {
        for (int j = 0; j < np; j++) {
            fphotonSpMatrix << scientific << setprecision(6) << setw(16)
                            << MInv_[im] << "  " << p[j] << "  "
                            << vnMInvpT_cos[0][im][j];
            for (int order = 1; order < norder_; order++) {
                fphotonSpMatrix << scientific << setprecision(6) << setw(16)
                                << "  " << vnMInvpT_cos[order][im][j] << "  "
                                << vnMInvpT_sin[order][im][j];
            }
            fphotonSpMatrix << endl;
        }
    }

    for (int im = 0; im < nMInv_; im++) {
        fphotonSpvn << scientific << setprecision(6) << setw(16) << MInv_[im]
                    << "  " << vnMInv_cos[0][im];
        for (int order = 1; order < norder_; order++) {
            fphotonSpvn << scientific << setprecision(6) << setw(16) << "  "
                        << vnMInv_cos[order][im] << "  "
                        << vnMInv_sin[order][im];
        }
        fphotonSpvn << endl;
    }

    fphotonSpMatrix.close();
    fphotonSpvn.close();
}

void ThermalDilepton::outputPhoton_SpvnpT_shell(string path) {
    outputPhoton_SpvnpT(
        path, "eq", dNpTdpTdphidydM_eq, vnMInvpT_cos_eq, vnMInvpT_sin_eq,
        vnMInv_cos_eq, vnMInv_sin_eq);
    outputPhoton_SpvnpT(
        path, "eqT", dNpTdpTdphidydM_eqT, vnMInvpT_cos_eqT, vnMInvpT_sin_eqT,
        vnMInv_cos_eqT, vnMInv_sin_eqT);
    outputPhoton_SpvnpT(
        path, "eqL", dNpTdpTdphidydM_eqL, vnMInvpT_cos_eqL, vnMInvpT_sin_eqL,
        vnMInv_cos_eqL, vnMInv_sin_eqL);
}
