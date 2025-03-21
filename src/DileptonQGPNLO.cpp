#include "DileptonQGPNLO.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_fermi_dirac.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Arsenal.h"
#include "data_struct.h"

using namespace std;
using PhysConsts::alphaEM;
using PhysConsts::hbarC;

using ARSENAL::createA4DMatrix;
using ARSENAL::deleteA4DMatrix;

void error_handler_nlo(
    const char *reason, const char *file, int line, int gsl_errno) {
    // fprintf(stderr, "GSL ERROR: %s:%d: %s (error code: %d)\n", file, line,
    // reason, gsl_errno);
}
DileptonQGPNLO::DileptonQGPNLO(
    std::shared_ptr<ParameterReader> paraRdr_in, std::string emissionProcess)
    : ThermalDilepton {paraRdr_in, emissionProcess} {
    ratePath_ = "ph_rates/";
    readInEmissionTables(emissionProcess);
    prefac1_ = (2. / 3.) * pow(hbarC, -4.) * pow(alphaEM, 2.);
    prefac2_ = 3. * pow(M_PI, 3.);
}

DileptonQGPNLO::~DileptonQGPNLO() {
    if (getRateTableFlag()) {
        deleteA4DMatrix(
            rateRhoT, alphaS_list.size(), muoverT_list.size(),
            MoverT_list.size());
        deleteA4DMatrix(
            rateRhoL, alphaS_list.size(), muoverT_list.size(),
            MoverT_list.size());
    }
}

/*--------------------------------------------------------------------*/
// some pQCD results

double DileptonQGPNLO::nF(double x) {
    double e = exp(-x);
    return e / (1. + e);
}

double DileptonQGPNLO::nB(double x) {
    double e = exp(-x);
    return e / (1. - e);
};

// double DileptonQGPNLO::l1f(double x) { return +gsl_sf_fermi_dirac_0(-x); }

// double DileptonQGPNLO::l2f(double x) { return -gsl_sf_fermi_dirac_1(-x); }

// double DileptonQGPNLO::l3f(double x) { return -gsl_sf_fermi_dirac_2(-x); }

double DileptonQGPNLO::l1f(double x) {
    gsl_set_error_handler(&error_handler_nlo);
    gsl_sf_result result;
    int status = gsl_sf_fermi_dirac_0_e(-x, &result);

    if (status != GSL_SUCCESS) {
        // printf("An error occurred: %f\n", result.val);
        // not printf, set result.val=0.0
    }
    return result.val;
}
double DileptonQGPNLO::l2f(double x) {
    gsl_set_error_handler(&error_handler_nlo);
    gsl_sf_result result;
    int status = gsl_sf_fermi_dirac_1_e(-x, &result);

    if (status != GSL_SUCCESS) {
        // printf("An error occurred: %f\n", result.val);
        // not printf, set result.val=0.0
    }
    return -result.val;
}
double DileptonQGPNLO::l3f(double x) {
    gsl_set_error_handler(&error_handler_nlo);
    gsl_sf_result result;
    int status = -gsl_sf_fermi_dirac_2_e(-x, &result);
    if (status != GSL_SUCCESS) {
        // printf("An error occurred: %f\n", result.val);
        // not printf, set result.val=0.0
    }

    return -result.val;
}

void DileptonQGPNLO::rho_LO(
    double o, double k, double mu, double &rT, double &rL) {
    // leading order result, see (2.4) of 1910.09567
    // (mu dependence not actually published yet...)
    double rV, r00;
    double K2 = o * o - k * k;  // note: omega, k are in units of T!
    double kp = .5 * (o + k), km = .5 * std::abs(o - k);
    double somk = 1;
    if ((o - k) < 0) somk = -1;

    rV = l1f(kp + mu) + l1f(kp - mu) - l1f(km + mu) - l1f(km - mu);
    rV = rV / k + .5 * (somk + 1.);
    rV *= K2 * 3. * OOFP;

    r00 = l3f(kp + mu) + l3f(kp - mu) - l3f(km + mu) - l3f(km - mu);
    r00 = r00 * 2. / k + l2f(kp + mu) + l2f(kp - mu)
          + somk * (l2f(km + mu) + l2f(km - mu));
    r00 = r00 * 6. + .5 * k * k * (somk + 1.);
    r00 *= -OOFP;

    rL = -K2 * r00 / (k * k);
    rT = .5 * (rV - rL);
}

double DileptonQGPNLO::rhoV_AMY(
    double alpha, double mu, double k) {  // (1.9) and (1.10) of hep-ph/0111107
    double nf = 3.;
    double Cf = (Nc * Nc - 1) / (2. * Nc);
    double minf2, fcol, f2t2;
    minf2 = 4. * M_PI * alpha * (Cf / 4.) * (1. + 1. * mu * mu / (M_PI * M_PI));
    f2t2 = .041 / k - .3615 + 1.01 * exp(-1.35 * k);
    fcol = sqrt(1. + nf / 6.)
           * (.548 * log(12.28 + 1. / k) * pow(k, -1.5)
              + .133 * k / sqrt(1. + k / 16.27));

    return .5 * (Nc / M_PI) * minf2 * (1. - nF(k - mu) - nF(k + mu))
           * (-.5 * log(minf2) + .5 * log(2. * k) + f2t2 + fcol);
}

void DileptonQGPNLO::rho_OPE(
    double o, double k, double alpha, double mu, double &rT, double &rL) {
    // operator product expansion, see appendix D of 1910.07552
    // (mu dependence not actually published yet...)
    double K2 = o * o - k * k;  // note: omega, k are in units of T!
    double coeff = 1. + 96. * mu * mu * OOFP * OOFP + 768. * pow(mu * OOFP, 4);
    // = T^4 + (6/pi^2).T^2.mu^2 + (3/pi^4).mu^4 )
    rT = OOFP * K2;
    rL = OOFP * K2;
    rT += alpha * (4. * OOFP * OOFP * K2 + coeff / (OOFP * OOFP * 27. * K2));
    rL += alpha
          * (4. * OOFP * OOFP * K2
             + coeff * (o * o + k * k) / (OOFP * OOFP * 27. * K2 * K2));
}

int DileptonQGPNLO::getIdx(double xval, std::vector<double> &xTable) {
    // binary search
    if (xval > xTable[xTable.size() - 1]) return xTable.size() - 2;
    if (xval < xTable[0]) return 0;
    int iU = xTable.size() - 1;
    int iL = 0;
    int iM = (iU + iL) / 2;
    while (iU - iL > 1) {
        if (xval >= xTable[iM])
            iL = iM;
        else
            iU = iM;
        iM = (iU + iL) / 2;
    }
    return iL;
}

void DileptonQGPNLO::interp(
    const double alpha, const double muOverT, const double MOverT,
    const double kOverT, double &resRhoT, double &resRhoL) {
    const int i = getIdx(alpha, alphaS_list);
    const int j = getIdx(muOverT, muoverT_list);
    const int k = getIdx(MOverT, MoverT_list);
    const int l = getIdx(kOverT, koverT_list);

    double a = (alpha - alphaS_list[i]) / (alphaS_list[i + 1] - alphaS_list[i]);
    double b =
        (muOverT - muoverT_list[j]) / (muoverT_list[j + 1] - muoverT_list[j]);
    double c =
        (MOverT - MoverT_list[k]) / (MoverT_list[k + 1] - MoverT_list[k]);
    double d =
        (kOverT - koverT_list[l]) / (koverT_list[l + 1] - koverT_list[l]);

    // avoid overflows
    a = std::max(0., std::min(1., a));
    b = std::max(0., std::min(1., b));
    c = std::max(0., std::min(1., c));
    d = std::max(0., std::min(1., d));

    resRhoT = 0;
    resRhoL = 0;
    for (int ialpha = 0; ialpha < 2; ialpha++) {
        for (int imu = 0; imu < 2; imu++) {
            for (int im = 0; im < 2; im++) {
                for (int ik = 0; ik < 2; ik++) {
                    double weight = std::abs(
                        ((1 - ialpha) - a) * ((1 - imu) - b) * ((1 - im) - c)
                        * ((1 - ik) - d));
                    resRhoT +=
                        weight * rateRhoT[i + ialpha][j + imu][k + im][l + ik];
                    resRhoL +=
                        weight * rateRhoL[i + ialpha][j + imu][k + im][l + ik];
                }
            }
        }
    }
}

void DileptonQGPNLO::readInEmissionTables(std::string emissionProcess) {
    std::ostringstream eqrate_filename_stream;
    eqrate_filename_stream << ratePath_ << "rate_" << emissionProcess
                           << "_eqrate.dat";
    std::cout << "reading in file: [" << eqrate_filename_stream.str() << "]"
              << std::endl;

    std::ifstream fin;
    fin.open(eqrate_filename_stream.str().c_str());
    if (fin.is_open() == false) {
        std::cout << "DileptonQGPNLO::readInEmissionTables error: "
                  << "the data file cannot be opened." << std::endl;
        exit(-1);
    }

    const int nAlphaS = 11;
    const int nMuB = 7;
    const int nM = 42;
    const int nK = 39;
    setRateTableFlag(true);
    rateRhoT = createA4DMatrix(nAlphaS, nMuB, nM, nK, 0.);
    rateRhoL = createA4DMatrix(nAlphaS, nMuB, nM, nK, 0.);

    fin.ignore(256, '\n');

    double alpha, muB, M, k;
    for (int ialpha = 0; ialpha < nAlphaS; ialpha++) {
        for (int imuB = 0; imuB < nMuB; imuB++) {
            for (int im = 0; im < nM; im++) {
                for (int ik = 0; ik < nK; ik++) {
                    fin >> alpha >> muB >> M >> k
                        >> rateRhoT[ialpha][imuB][im][ik]
                        >> rateRhoL[ialpha][imuB][im][ik];
                    if (ialpha == 0 && imuB == 0 && im == 0)
                        koverT_list.push_back(k);
                    if (ialpha == 0 && imuB == 0 && ik == 0)
                        MoverT_list.push_back(M);
                    if (ialpha == 0 && im == 0 && ik == 0)
                        muoverT_list.push_back(muB);
                    if (imuB == 0 && im == 0 && ik == 0)
                        alphaS_list.push_back(alpha);
                }
            }
        }
    }
    std::cout << " ... done!" << std::endl;
}

/*--------------------------------------------------------------------*/
// interpolation & extrapolation
int DileptonQGPNLO::approx_rho(double *input, double &rT, double &rL) {
    // int function returns: -1 if below LC (not needed for dileptons)
    //                       0 if interpolated (i.e. input values fall within
    //                       grid) 1 if extrapolated (*)
    //
    // (*) extrapolation will NOT work properly for large alpha_s, large muB, or
    // if k falls outside the grid! In these cases, just the boundary value is
    // used.

    double EoverT = input[0];
    double koverT = input[1];
    double MoverT = input[2];
    double alphaS = input[3];
    double muoverT = input[4];

    double al_min = alphaS_list[0];
    double M_min = MoverT_list[0];
    double M_max = MoverT_list[MoverT_list.size() - 1];

    if (alphaS < al_min - 1e-5) {
        // use LO result to extrapolate if alpha < alpha_min
        double LO_rT, LO_rL;
        rho_LO(EoverT, koverT, muoverT, LO_rT, LO_rL);
        double min_rT, min_rL;
        interp(al_min, muoverT, MoverT, koverT, min_rT, min_rL);

        double weight = alphaS / al_min;
        rT = weight * min_rT + (1. - weight) * LO_rT;
        rL = weight * min_rL + (1. - weight) * LO_rL;
        return 1;
    }
    if (MoverT > M_max + 1e-5) {
        // use OPE expansion to extrapolate if M > M_max
        double OPE_rT, OPE_rL;
        rho_OPE(EoverT, koverT, alphaS, muoverT, OPE_rT, OPE_rL);
        rT = OPE_rT;
        rL = OPE_rL;
        return 1;
    } else if (MoverT < M_min - 1e-5) {
        // use AMY approx. formula if M < M_min
        double min_rT, min_rL;
        interp(alphaS, muoverT, M_min, koverT, min_rT, min_rL);

        double weight = MoverT / M_min;
        rT = weight * min_rT
             + (1. - weight) * rhoV_AMY(alphaS, muoverT, koverT) / 2.;
        rL = weight * min_rL;
        return 1;
    }

    interp(alphaS, muoverT, MoverT, koverT, rT, rL);
    return 0;
}

double DileptonQGPNLO::mllFactor(double x) {
    // lepton pair kinematic factor
    double component = 1. - 4. * x;
    if (component < 0.)
        return 0.;
    else
        return (1. + 2. * x) * sqrt(component);
}

void DileptonQGPNLO::getRateFromTable(
    const double E, const double k, const double MInv, const double alpha_s,
    const double muB, const double T, const double mlsq, double &rateTot,
    double &rateT, double &rateL) {
    // NB: quantities are defined in the local rest frame, i.e. u_\mu=(1,0,0,0)

    if (E < k) {
        // skip 'DIS region'
        rateT = 0.;
        rateL = 0.;
        rateTot = 0.;
        return;
    }

    double M2 = MInv * MInv;

    // double prefactor = mllFactor(mlsq / M2) * (2. / 3.) * pow(hbarC, -4.)
    //                    * pow(alphaEM, 2.);
    //  units: GeV^-4.fm^-4
    double prefactor = mllFactor(mlsq / M2) * prefac1_;

    // m_l, T, muB, o and k must all be given in units of GeV!
    double in[5] = {E / T, k / T, MInv / T, alpha_s, muB / T / 3.};
    double rhoT_app, rhoL_app;
    approx_rho(in, rhoT_app, rhoL_app);

    double prefac = prefactor * nB(E / T) * T * T / (prefac2_ * M2);

    rateT = prefac * rhoT_app;
    rateL = prefac * rhoL_app;

    rateTot = 2. * rateT + rateL;
}
