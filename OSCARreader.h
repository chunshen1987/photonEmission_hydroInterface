#ifndef OSCARREADER_H
#define OSCARREADER_H

#include<fstream>

#include "parameter.h"

using namespace std;

typedef struct
{
      char filename[15];
      char hydrotype[10];     //viscous or ideal
      char datainf[15];       //specify data only contain isothemal hypersurface or whole history of the evolution
      char initial[10];       //type of initialize energy density
      char nucleartype[10];   //type of collision nucleus
      char impactb[8];        //impact parameter
      char maxentropy[10];    //maximum entropy density at first step
      char thermaltime[10];   //thermalization time
      char etas[15];         //value of shear viscosity
      char EOS[10];           //type equation of state
      char geometry[15];      //geometry of the model: (2+1)d or (3+1)d
      int timestep;           //number of timesteps in the file
      int n_x;                 //number of grids along x
      int n_y;                 //number of grids along y
      double ti;              //starting time
      double tf;              //ending time
      double xi;              //starting value of x
      double xf;              //ending value of x
      double yi;              //starting value of y
      double yf;              //ending value of y
}OSCARheader;

typedef struct
{
      double ed[nx][ny], pl[nx][ny], temp[nx][ny]; //energy density, pressure, and temperature
      double vx[nx][ny], vy[nx][ny]; //flow velocity in the transverse plane 
      double pi00[nx][ny], pi01[nx][ny], pi02[nx][ny], pi03[nx][ny];
      double pi11[nx][ny], pi12[nx][ny], pi13[nx][ny]; 
      double pi22[nx][ny], pi23[nx][ny];
      double pi33[nx][ny], PPi[nx][ny]; //shear and bulk pressure tensor
      double eta_s[nx][ny], taupi[nx][ny], tauPib[nx][ny]; // transport coefficients, \eta/s, and relaxation time
      double r_qgp[nx][ny]; //parameter: 1 for QGP phase and 0 for hadronic phase
      int itime;
}readindata;

class OSCARreader
{
   private:
      ifstream* fp;
      OSCARheader* headerptr;
   
   public:
      OSCARreader(string );
      ~OSCARreader();

      void readheader(); //read the header of OSCAR file
      void printinfo(); //print out the information of the read in file
      void printframe(readindata* dataptr, string filename);

      double getGriddt() {return((headerptr->tf-headerptr->ti)/(headerptr->timestep-1));};
      //double getGriddt() {return(0.02);}
      double getGriddx() {return((headerptr->xf-headerptr->xi)/(headerptr->n_x-1));};
      double getGriddy() {return((headerptr->yf-headerptr->yi)/(headerptr->n_y-1));};
      double getGridt0() {return(headerptr->ti);};
      int getTimestep() {return(headerptr->timestep);};

      void readframe_2Dboostinvariant(readindata*); //read one frame (one time step) of data
      void readnextframe(readindata** , readindata** ); // read in next frame
};

#endif
