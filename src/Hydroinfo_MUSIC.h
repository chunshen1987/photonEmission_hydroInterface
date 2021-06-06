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
#include "data_struct.h"


class Hydroinfo_MUSIC {
 private:
    float hydroTau0;       // tau_0 in the hydro data files
    float hydroTauMax;     // tau_max in the hydro data files
    float hydroDtau;       // step dtau in fm/c in the hydro data files
    float hydroXmax;       // maximum x in fm in the hydro data files
                            // [-xmax, +xmax] for both x and y
    float hydro_eta_max;       // maximum z in fm in the hydro data files
                            // [-zmax, +zmax] for 3D hydro
    float hydroDx;         // step dx in fm in the hydro data files
    float hydroDeta;         // step dz in fm in the hydro data files in
                            // the z-direction for 3D hydro

    int nskip_tau, nskip_x, nskip_eta;

    int hydroWhichHydro;    // choose a hydro evolution model to use
    int use_tau_eta_coordinate;

    bool boost_invariant_;

    int itaumax, ixmax, ietamax;
    int turn_on_shear;
    int turn_on_bulk;
    int turn_on_rhob;
    int turn_on_diff;

    // array to store hydro information
    std::vector<fluidCell_2D> lattice_2D;
    std::vector<fluidCell_3D> lattice_3D;
    std::vector<fluidCell_3D_new> lattice_new_;
    std::vector<fluidCell_3D_ideal> lattice_3D_ideal;
    std::vector<int> idx_map_;

 public:
    Hydroinfo_MUSIC();       // constructor
    ~Hydroinfo_MUSIC();      // destructor

    float get_hydro_tau_max() {return(hydroTauMax);}
    float get_hydro_tau0() {return(hydroTau0);}
    float get_hydro_dtau() {return(hydroDtau);}
    float get_hydro_dx() {return(hydroDx);}
    float get_hydro_deta() {return(hydroDeta);}
    float get_hydro_eta_max() {return(hydro_eta_max);}
    float get_hydro_x_max() {return(hydroXmax);}
    int get_hydro_Nskip_tau() {return(nskip_tau);}
    int get_hydro_Nskip_x() {return(nskip_x);}
    int get_hydro_Nskip_eta() {return(nskip_eta);}
    int get_number_of_fluid_cells_3d() {return(lattice_new_.size());}

    bool isBoostInvariant() const {return(boost_invariant_);}

    void readHydroData(int whichHydro, int nskip_tau_in);

    void getHydroValues(float x, float y, float z, float t,
                        fluidCell *info);
    void get_hydro_cell_info_3d(int cell_id, fluidCell_3D_new &info);
    void output_temperature_evolution(std::string filename_base);
    void update_grid_info(float tau0, float tau_max, float dtau,
                          float x_max, float dx, float z_max, float dz);
};

#endif  // SRC_HYDROINFO_MUSIC_H_

