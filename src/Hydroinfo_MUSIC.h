// Hydroinfo_MUSIC.h is a part of the MARTINI event generator.
// Copyright (C) 2009 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions
// that return interpolated data at a given space-time point

#ifndef SRC_HYDROINFO_MUSIC_H_
#define SRC_HYDROINFO_MUSIC_H_

#include <vector>
#include <string>
#include "./Hydroinfo_h5.h"

class Hydroinfo_MUSIC {
 private:
    double hydroTau0;       // tau_0 in the hydro data files
    double hydroTauMax;     // tau_max in the hydro data files
    double hydroDtau;       // step dtau in fm/c in the hydro data files
    double hydroXmax;       // maximum x in fm in the hydro data files
                            // [-xmax, +xmax] for both x and y
    double hydroZmax;       // maximum z in fm in the hydro data files
                            // [-zmax, +zmax] for 3D hydro
    double hydroDx;         // step dx in fm in the hydro data files
    double hydroDz;         // step dz in fm in the hydro data files in
                            // the z-direction for 3D hydro

    int hydroWhichHydro;    // choose a hydro evolution model to use
    int use_tau_eta_coordinate;

    bool boost_invariant;

    int itaumax, ixmax, ietamax;

    std::vector<fluidCell> *lattice;     // array to store hydro information

 public:
    Hydroinfo_MUSIC();       // constructor
    ~Hydroinfo_MUSIC();      // destructor

    double get_hydro_tau_max() {return(hydroTauMax);}

    void readHydroData(double tau0, double taumax, double dtau,
                       double xmax, double zmax, double dx, double dz,
                       int nskip_tau, int nskip_x, int nskip_z,
                       int whichHydro, std::string evolution_name);

    void getHydroValues(double x, double y, double z, double t, 
                        fluidCell *info);
    void output_temperature_evolution(std::string filename_base);
    void update_grid_info(double tau0, double tau_max, double dtau,
                          double x_max, double dx, double z_max, double dz);
};

#endif  // SRC_HYDROINFO_MUSIC_H_

