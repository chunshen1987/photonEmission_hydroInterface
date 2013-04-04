#ifndef Hydroinfo_h5_H
#define Hydroinfo_h5_H

#include<fstream>
#include<sstream>
#include<string>

#include "ParameterReader.h"
#include "hdf5.h"

using namespace std;

struct fluidCell {
   double ed, sd, vx, vy, temperature, pressure;
   double pi[4][4];
   double bulkPi;
};

class HydroinfoH5
{
   private:
      ParameterReader* paraRdr;
      int Visflag;  // flag to determine whether to read evolutions for viscous variables

      int Buffersize;
      hid_t H5file_id, H5groupEventid;

      int grid_XL, grid_XH, grid_YL, grid_YH;
      int grid_Framenum;
      double grid_X0, grid_Y0;
      double grid_Xmax, grid_Ymax;
      double grid_Tau0, grid_dTau, grid_dx, grid_dy;
      double grid_Taumax;

      int dimensionX, dimensionY;
      double ***ed, ***sd, ***vx, ***vy, ***Temperature, ***Pressure;
      double ***pi00, ***pi01, ***pi02, ***pi03, ***pi11, ***pi12, ***pi13, ***pi22, ***pi23, ***pi33;
      double ***BulkPi;

   public:
      HydroinfoH5(string file, ParameterReader* paraRdr_in);
      ~HydroinfoH5();
      
      void readHydrogridInfo();
      void printHydrogridInfo();
      int readH5Attribute_int(hid_t id, string attributeName);
      double readH5Attribute_double(hid_t id, string attributeName);

      void readHydroinfoBuffered_total();
      void readHydroinfoSingleframe(int frameIdx);
      void readH5Dataset_double(hid_t id, string datasetName, double** dset_data);

      int getNumberofFrames() {return((int)grid_Framenum);};
      double getHydrogridDX() {return(grid_dx);};
      double getHydrogridDY() {return(grid_dy);};
      double getHydrogridDTau() {return(grid_dTau);};
      double getHydrogridTau0() {return(grid_Tau0);};
      double getHydrogridTaumax() {return(grid_Taumax);};
      double getHydrogridNX() {return(grid_XH - grid_XL + 1);};
      double getHydrogridNY() {return(grid_YH - grid_YL + 1);};
      void getHydroinfoOnlattice(int frameIdx, int xIdx, int yIdx, fluidCell* fluidCellptr);
      void getHydroinfo(double tau, double x, double y, fluidCell* fluidCellptr);
      void setZero_fluidCell(fluidCell* fluidCellptr);

};

#endif
