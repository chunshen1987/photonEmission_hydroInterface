// Hydroinfo_MUSIC.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions
// that return interpolated data at a given space-time point

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

#include "./Hydroinfo_h5.h"
#include "./Hydroinfo_MUSIC.h"

using namespace std;

Hydroinfo_MUSIC::Hydroinfo_MUSIC() {
    lattice_2D = new vector<fluidCell_2D>;
    lattice_3D = new vector<fluidCell_3D>;
}

Hydroinfo_MUSIC::~Hydroinfo_MUSIC() {
    if (boost_invariant) {
        lattice_2D->clear();
    } else {
        lattice_3D->clear();
    }
    delete lattice_2D;
    delete lattice_3D;
}

void Hydroinfo_MUSIC::readHydroData(
    double tau0, double taumax, double dtau,
    double xmax, double zmax, double dx, double dz,
    int nskip_tau, int nskip_x, int nskip_z, int whichHydro) {
// all hydro data is stored in tau steps (not t) - the t and z in the MARTINI
// evolution is converted to tau when accessing the hydro data
    lattice_2D->clear();
    lattice_3D->clear();

    // get hydro grid information
    hydroTau0 = tau0;
    hydroTauMax = taumax;
    hydroDtau = dtau*nskip_tau;
    hydroXmax = xmax;
    hydroZmax = zmax;
    hydroDx = dx*nskip_x;
    hydroDz = dz*nskip_z;

    hydroWhichHydro = whichHydro;
    use_tau_eta_coordinate = 1;

    if (use_tau_eta_coordinate == 0) {
        cout << "Hydroinfo_MUSIC:: Warning hydro grid is set to "
             << "cartesian coordinates, please make sure this is correct!"
             << endl;
    }

    if (whichHydro != 6 && whichHydro != 8) {
        cout << "Hydroinfo_MUSIC:: This option is obsolete! whichHydro = "
             << whichHydro << endl;
        exit(1);
    } else if (whichHydro == 6) {
        // 3+1D MUSIC hydro (Schenke, Jeon, Gale)
        cout << "Using 3+1D Jeon Schenke hydro reading data ..." << endl;
        boost_invariant = false;

        itaumax = static_cast<int>((taumax-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*xmax/dx+0.001);
        ietamax = static_cast<int>(2*zmax/dz+0.001);

        // read in temperature, QGP fraction , flow velocity
        // The name of the evolution file: evolution_name
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
            "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";
        cout << "Evolution file name = " << evolution_name << endl;
        ifstream fin;
        fin.open(evolution_name.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }
        ifstream fin1;
        fin1.open(evolution_name_Wmunu.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Wmunu << endl;
            exit(1);
        }
        ifstream fin2;
        fin2.open(evolution_name_Pi.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Pi << endl;
            exit(1);
        }

        double T, vx, vy, vz, QGPfrac;
        fluidCell_3D newCell;
        int ik = 0;
        while (!fin.eof()) {
            ik++;
            fin >> T;
            fin >> QGPfrac;
            fin >> vx;
            fin >> vy;
            fin >> vz;
            newCell.temperature = T;
            newCell.vx = vx;
            newCell.vy = vy;
            newCell.vz = vz;

            lattice_3D->push_back(newCell);
            if (ik%50000 == 0)
                cout << "o" << flush;
        }
        cout << ik << endl;
        fin.close();
        fin1.close();
        fin2.close();
    } else if (whichHydro == 8) {
        // event-by-event (2+1)-d MUSIC hydro from JF
        // there are two slices in medium in eta_s
        // one at eta_s = -15. and the other at eta_s = 0.0
        // only the medium at middle rapidity will be kept in the memory
        boost_invariant = true;
        cout << "Reading event-by-event hydro evolution data from JF ..."
             << endl;
        ixmax = static_cast<int>(2*xmax/dx+0.001);
        ietamax = 1;

        int n_eta = 2;  // there are two slices in eta_s
        // number of fluid cell in the transverse plane
        int num_fluid_cell_trans = ixmax*ixmax;

        // read in hydro evolution
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
                "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "resutls/evolution_bulk_pressure_xyeta.dat";

        std::FILE *fin;
        string evolution_file_name = evolution_name;
        cout << "Evolution file name = " << evolution_file_name << endl;
        fin = std::fopen(evolution_file_name.c_str(), "rb");

        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_file_name << endl;
            exit(1);
        }

        std::FILE *fin1;
        fin1 = std::fopen(evolution_name_Wmunu.c_str(), "rb");
        if (fin1 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Wmunu << endl;
            exit(1);
        }

        std::FILE *fin2;
        fin2 = std::fopen(evolution_name_Pi.c_str(), "rb");
        if (fin1 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Pi << endl;
            exit(1);
        }

        int ik = 0;
        fluidCell_2D newCell;
        float T, QGPfrac, vx, vy, vz;
        int size = sizeof(float);
        while (true) {
            int status = 0;
            status = std::fread(&T, size, 1, fin);
            status += std::fread(&QGPfrac, size, 1, fin);
            status += std::fread(&vx, size, 1, fin);
            status += std::fread(&vy, size, 1, fin);
            status += std::fread(&vz, size, 1, fin);

            if (status != 5) {  // this is the end of file
                break;
            }

            int ieta_idx = static_cast<int>(ik/num_fluid_cell_trans) % n_eta;
            int itau_idx = static_cast<int>(ik/(num_fluid_cell_trans*n_eta));
            ik++;
            if (itau_idx%nskip_tau != 0)  // skip in tau
                continue;

            // print out tau information
            double tau_local = hydroTau0 + itau_idx*hydroDtau/nskip_tau;
            if ((ik-1)%(num_fluid_cell_trans*n_eta) == 0) {
                cout << "read in tau frame: " << itau_idx
                     << " tau_local = " << setprecision(3) << tau_local
                     << " fm ..."<< endl;
            }

            if (ieta_idx == (n_eta-1)) {
                // store the hydro medium at eta_s = 0.0
                // vz = 0 at eta_s = 0
                double gamma_L = 1./sqrt(1. - vz*vz);
                if (gamma_L > 1.01) {
                    cout << "gamma_L :" << gamma_L << endl;
                }
                newCell.temperature = T;
                // convert vx and vy to longitudinal co-moving frame
                newCell.vx = vx*gamma_L;
                newCell.vy = vy*gamma_L;

                // pi^\mu\nu tensor
                newCell.pi00 = 0.0;
                newCell.pi01 = 0.0;
                newCell.pi02 = 0.0;
                newCell.pi11 = 0.0;
                newCell.pi12 = 0.0;
                newCell.pi22 = 0.0;
                newCell.pi33 = 0.0;

                // bulk pressure
                newCell.bulkPi = 0.0;

                lattice_2D->push_back(newCell);
            }
        }
        std::fclose(fin);
        std::fclose(fin1);
        std::fclose(fin2);
        cout << endl;
        cout << "number of fluid cells: " << lattice_2D->size() << endl;
    }

    // One final step for easy automation of MARTINI:
    // hydroTauMax is reset for the case where writing to evolution.dat
    // ended early (due to all cells freezing out):
    if (whichHydro == 6) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*static_cast<int>(
                        static_cast<double>(lattice_3D->size())
                        /((2.*hydroXmax/hydroDx+1.)*(2.*hydroXmax/hydroDx+1.)
                        *2.*(hydroZmax/hydroDz))));
        itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);
    }
    if (whichHydro == 8) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*static_cast<int>(
                        static_cast<double>(lattice_2D->size())
                        /((2.*hydroXmax/hydroDx)*(2.*hydroXmax/hydroDx)) - 1));
        itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);
    }

    cout << "hydroTau0 = " << hydroTau0 << endl;
    cout << "hydroTauMax = " << hydroTauMax << endl;
    cout << "hydroDtau = " << hydroDtau << endl;
    cout << "hydroXmax = " << hydroXmax << endl;
    cout << "hydroDx = " << hydroDx << endl;
    cout << "hydroZmax = " << hydroZmax << endl;
    cout << "hydroDz = " << hydroDz << endl;
}


void Hydroinfo_MUSIC::getHydroValues(double x, double y,
                                     double z, double t, fluidCell* info) {
// For interpolation of evolution files in tau-eta coordinates. Only the
// reading of MUSIC's evolution_xyeta.dat file is implemented here.
// For simplicity, hydroZmax refers to MUSIC's eta_size, and similarly for
// hydroDz; however, x, y, z, and t are as usual to stay compatible with
// MARTINI.
    double tau, eta;
    if (use_tau_eta_coordinate == 1) {
        if (t*t > z*z) {
            tau = sqrt(t*t-z*z);
            eta = 0.5*log((t+z)/(t-z));
        } else {
            tau = 0.;
            eta = 0.;
        }
    } else {
        // if the medium is given in cartesian coordinates
        // set tau and eta to t and z
        tau = t;
        eta = z;
    }

    int ieta = floor((hydroZmax+eta)/hydroDz + 0.0001);
    if (hydroWhichHydro == 8)
        ieta = 0;

    int itau = floor((tau-hydroTau0)/hydroDtau + 0.0001);
    int ix = floor((hydroXmax+x)/hydroDx + 0.0001);
    int iy = floor((hydroXmax+y)/hydroDx + 0.0001);

    double xfrac = (x - (static_cast<double>(ix)*hydroDx - hydroXmax))/hydroDx;
    double yfrac = (y - (static_cast<double>(iy)*hydroDx - hydroXmax))/hydroDx;
    double etafrac = (eta/hydroDz - static_cast<double>(ieta)
                      + 0.5*static_cast<double>(ietamax));
    double taufrac = (tau - hydroTau0)/hydroDtau - static_cast<double>(itau);

    if (ix < 0 || ix >= ixmax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - x out of range x=" << x
             << ", ix=" << ix << ", ixmax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (iy < 0 || iy >= ixmax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy
             << ", iymax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (itau < 0 || itau >= itaumax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << itaumax
             << endl;
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: tau= " << tau
             << ", hydroTauMax = " << hydroTauMax
             << ", hydroDtau = " << hydroDtau << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (ieta < 0 || ieta >= ietamax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "eta out of range, ieta=" << ieta << ", ietamax=" << ietamax
             << endl;
        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }

  // The array of positions on the 4-dimensional rectangle:
  int position[2][2][2][2];
  for (int ipx = 0; ipx < 2; ipx++) {
        int px;
        if (ipx == 0 || ix == ixmax-1)
            px = ix;
        else
            px = ix + 1;
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == ixmax-1)
                py = iy;
            else
                py = iy + 1;
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1)
                    peta = ieta;
                else
                    peta = ieta + 1;
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax-1)
                        ptau = itau;
                    else
                        ptau = itau + 1;
                    position[ipx][ipy][ipeta][iptau] = (
                                px + ixmax*(py + ixmax*(peta + ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    double T = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double pi00 = 0.0;
    double pi01 = 0.0;
    double pi02 = 0.0;
    double pi03 = 0.0;
    double pi11 = 0.0;
    double pi12 = 0.0;
    double pi13 = 0.0;
    double pi22 = 0.0;
    double pi23 = 0.0;
    double pi33 = 0.0;
    double bulkPi = 0.0;

    fluidCell_2D *HydroCell_2D_ptr1, *HydroCell_2D_ptr2;
    fluidCell_3D *HydroCell_3D_ptr1, *HydroCell_3D_ptr2;
    for (int iptau = 0; iptau < 2; iptau++) {
        double taufactor;
        if (iptau == 0)
            taufactor = 1. - taufrac;
        else
            taufactor = taufrac;
        for (int ipeta = 0; ipeta < 2; ipeta++) {
            double etafactor;
            if (ipeta == 0)
                etafactor = 1. - etafrac;
            else
                etafactor = etafrac;
            for (int ipy = 0; ipy < 2; ipy++) {
                double yfactor;
                if (ipy == 0)
                    yfactor = 1. - yfrac;
                else
                    yfactor = yfrac;

                double prefrac = yfactor*etafactor*taufactor;

                if (boost_invariant) {
                    HydroCell_2D_ptr1 = (
                            &(*lattice_2D)[position[0][ipy][ipeta][iptau]]);
                    HydroCell_2D_ptr2 = (
                            &(*lattice_2D)[position[1][ipy][ipeta][iptau]]);
                    T += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->temperature
                                  + xfrac*HydroCell_2D_ptr2->temperature);
                    vx += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->vx
                                    + xfrac*HydroCell_2D_ptr2->vx);
                    vy += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->vy
                                    + xfrac*HydroCell_2D_ptr2->vy);
                    pi00 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi00
                                    + xfrac*HydroCell_2D_ptr2->pi00);
                    pi01 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi01
                                    + xfrac*HydroCell_2D_ptr2->pi01);
                    pi02 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi02
                                    + xfrac*HydroCell_2D_ptr2->pi02);
                    pi11 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi11
                                    + xfrac*HydroCell_2D_ptr2->pi11);
                    pi12 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi12
                                    + xfrac*HydroCell_2D_ptr2->pi12);
                    pi22 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi22
                                    + xfrac*HydroCell_2D_ptr2->pi22);
                    pi33 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi33
                                    + xfrac*HydroCell_2D_ptr2->pi33);
                    bulkPi += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->bulkPi
                                    + xfrac*HydroCell_2D_ptr2->bulkPi);
                } else {
                    HydroCell_3D_ptr1 = (
                            &(*lattice_3D)[position[0][ipy][ipeta][iptau]]);
                    HydroCell_3D_ptr2 = (
                            &(*lattice_3D)[position[1][ipy][ipeta][iptau]]);
                    T += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->temperature
                                  + xfrac*HydroCell_3D_ptr2->temperature);
                    vx += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->vx
                                    + xfrac*HydroCell_3D_ptr2->vx);
                    vy += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->vy
                                    + xfrac*HydroCell_3D_ptr2->vy);
                    vz += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->vz
                                    + xfrac*HydroCell_3D_ptr2->vz);
                    pi00 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi00
                                    + xfrac*HydroCell_3D_ptr2->pi00);
                    pi01 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi01
                                    + xfrac*HydroCell_3D_ptr2->pi01);
                    pi02 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi02
                                    + xfrac*HydroCell_3D_ptr2->pi02);
                    pi03 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi03
                                    + xfrac*HydroCell_3D_ptr2->pi03);
                    pi11 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi11
                                    + xfrac*HydroCell_3D_ptr2->pi11);
                    pi12 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi12
                                    + xfrac*HydroCell_3D_ptr2->pi12);
                    pi13 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi13
                                    + xfrac*HydroCell_3D_ptr2->pi13);
                    pi22 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi22
                                    + xfrac*HydroCell_3D_ptr2->pi22);
                    pi23 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi23
                                    + xfrac*HydroCell_3D_ptr2->pi23);
                    pi33 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi33
                                    + xfrac*HydroCell_3D_ptr2->pi33);
                    bulkPi += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->bulkPi
                                    + xfrac*HydroCell_3D_ptr2->bulkPi);
                }
            }
        }
    }

    if (boost_invariant) {      // for boost invariant medium
        vz = z/t;               // Bjorken flow in the lab frame
        double gamma_L_inv = sqrt(1. - vz*vz);
        vx = vx*gamma_L_inv;    // convert vx and vy to lab frame
        vy = vy*gamma_L_inv;
    }

    info->temperature = T;
    info->vx = vx;
    info->vy = vy;
    info->vz = vz;

    info->ed = 1.0;
    info->sd = 0.0;
    info->pressure = 0.0;

    info->pi[0][0] = pi00;
    info->pi[0][1] = pi01;
    info->pi[0][2] = pi02;
    info->pi[0][3] = pi03;
    info->pi[1][0] = pi01;
    info->pi[1][1] = pi11;
    info->pi[1][2] = pi12;
    info->pi[1][3] = pi13;
    info->pi[2][0] = pi02;
    info->pi[2][1] = pi12;
    info->pi[2][2] = pi22;
    info->pi[2][3] = pi23;
    info->pi[3][0] = pi03;
    info->pi[3][1] = pi13;
    info->pi[3][2] = pi23;
    info->pi[3][3] = pi33;

    info->bulkPi = bulkPi;
    return;
}

void Hydroinfo_MUSIC::output_temperature_evolution(string filename_base) {
    fluidCell *hydroInfo = new fluidCell;
    for (int i = 0; i < itaumax; i++) {
        double tau = hydroTau0 + i*hydroDtau;
        ostringstream filename;
        filename << filename_base << "_tau_" << tau << ".dat";
        ofstream temp_evo(filename.str().c_str());
        for (int ix = 0; ix < ixmax; ix++) {
            double x_local = -hydroXmax + ix*hydroDx;
            for (int iy = 0; iy < ixmax; iy++) {
                double y_local = -hydroXmax + iy*hydroDx;
                getHydroValues(x_local, y_local, 0.0, tau, hydroInfo);
                double temp_local = hydroInfo->temperature;
                temp_evo << scientific << setw(16) << setprecision(8)
                         << temp_local << "   ";
            }
            temp_evo << endl;
        }
        temp_evo.close();
    }
    delete hydroInfo;
}

void Hydroinfo_MUSIC::update_grid_info(
    double tau0, double tau_max, double dtau,
    double x_max, double dx, double z_max, double dz) {
    hydroTau0 = tau0;
    hydroTauMax = tau_max;
    hydroDtau = dtau;
    hydroXmax = x_max;
    hydroDx = dx;
    hydroZmax = z_max;
    hydroDz = dz;
    if (hydroWhichHydro == 8) {
        itaumax = static_cast<int>((tau_max-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*x_max/dx+0.001);
        ietamax = static_cast<int>(2*z_max/dz+0.001);
    }
    if (hydroWhichHydro == 6) {
        itaumax = static_cast<int>((tau_max-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*x_max/dx+0.001);
        ietamax = static_cast<int>(2*z_max/dz+0.001);
    }
}
