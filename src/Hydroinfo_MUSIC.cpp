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
    lattice = new vector<fluidCell>;
}

Hydroinfo_MUSIC::~Hydroinfo_MUSIC() {
    lattice->clear();
    delete lattice;
}

void Hydroinfo_MUSIC::readHydroData(
    double tau0, double taumax, double dtau,
    double xmax, double zmax, double dx, double dz,
    int nskip_tau, int nskip_x, int nskip_z,
    int whichHydro, double Tfinal, int taueta_coord, string evolution_name) {
// all hydro data is stored in tau steps (not t) - the t and z in the MARTINI
// evolution is converted to tau when accessing the hydro data
    lattice->clear();

    // get hydro grid information
    hydroTau0 = tau0;
    hydroTauMax = taumax;
    hydroDtau = dtau;
    hydroXmax = xmax;
    hydroZmax = zmax;
    hydroDx = dx;
    hydroDz = dz;

    hydroWhichHydro = whichHydro;
    use_tau_eta_coordinate = taueta_coord;

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
        cout << "Evolution file name = " << evolution_name << endl;
        ifstream fin;
        fin.open(evolution_name.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }

        double T, vx, vy, vz, QGPfrac;
        fluidCell newCell;
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

            lattice->push_back(newCell);
            if (ik%50000 == 0)
                cout << "o" << flush;
        }
        cout << ik << endl;
        fin.close();
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

        // read in temperature, QGP fraction , flow velocity
        std::FILE *fin;
        string evolution_file_name = evolution_name;
        cout << "Evolution file name = " << evolution_file_name << endl;
        fin = std::fopen(evolution_file_name.c_str(), "rb");

        float T, QGPfrac, vx, vy, vz;
        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_file_name << endl;
            exit(1);
        }

        int ik = 0;
        fluidCell newCell;
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
                newCell.temperature = T;
                // convert vx and vy to longitudinal co-moving frame
                newCell.vx = vx*gamma_L;
                newCell.vy = vy*gamma_L;
                // vz is Bjorken flow no need to assign values
                newCell.vz = 0.0;

                lattice->push_back(newCell);
            }
        }
        std::fclose(fin);
        cout << endl;
        cout << "number of fluid cells: " << lattice->size() << endl;
    }

    // One final step for easy automation of MARTINI:
    // hydroTauMax is reset for the case where writing to evolution.dat
    // ended early (due to all cells freezing out):
    if (whichHydro == 6) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*(int)((double)lattice->size()
                        /((2.*hydroXmax/hydroDx+1.)*(2.*hydroXmax/hydroDx+1.)
                        *2.*(hydroZmax/hydroDz))));
        itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);
    }
    if (whichHydro == 8) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*(int)((double)lattice->size()
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

    int ieta = floor((hydroZmax+eta)/hydroDz+0.0001);
    if (hydroWhichHydro == 8)
        ieta = 0;

    int itau = floor((tau-hydroTau0)/hydroDtau+0.0001);
    int ix = floor((hydroXmax+x)/hydroDx+0.0001);
    int iy = floor((hydroXmax+y)/hydroDx+0.0001);

    double xfrac = (x - ((double)ix*hydroDx-hydroXmax))/hydroDx;
    double yfrac = (y - ((double)iy*hydroDx-hydroXmax))/hydroDx;
    double etafrac = eta/hydroDz-(double)ieta+0.5*(double)ietamax;
    double taufrac = (tau-hydroTau0)/hydroDtau-(double)itau;

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
            px = ix+1;
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == ixmax-1)
                py = iy;
            else
                py = iy+1;
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1)
                    peta = ieta;
                else
                    peta = ieta+1;
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax-1)
                        ptau = itau;
                    else
                        ptau = itau+1;
                    position[ipx][ipy][ipeta][iptau] = (
                                    px+ixmax*(py+ixmax*(peta+ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    double T = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;

    fluidCell HydroCell_temp1, HydroCell_temp2;
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

                HydroCell_temp1 = (*lattice)[position[0][ipy][ipeta][iptau]];
                HydroCell_temp2 = (*lattice)[position[1][ipy][ipeta][iptau]];

                T += prefrac*((1. - xfrac)*HydroCell_temp1.temperature
                              + xfrac*HydroCell_temp2.temperature);
                vx += prefrac*((1. - xfrac)*HydroCell_temp1.vx
                               + xfrac*HydroCell_temp2.vx);
                vy += prefrac*((1. - xfrac)*HydroCell_temp1.vy
                               + xfrac*HydroCell_temp2.vy);
                if (hydroWhichHydro != 8) {
                    vz += prefrac*((1. - xfrac)*HydroCell_temp1.vz
                                   + xfrac*HydroCell_temp2.vz);
                }
            }
        }
    }

    if (hydroWhichHydro == 8) {  // for boost invariant medium
        vz = z/t;               // Bjorken flow in the lab frame
        double gamma_L_inv = sqrt(1. - vz*vz);
        vx = vx*gamma_L_inv;    // convert vx and vy to lab frame
        vy = vy*gamma_L_inv;
    }

    info->temperature = T;
    info->vx = vx;
    info->vy = vy;
    info->vz = vz;
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
