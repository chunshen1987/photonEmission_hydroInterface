#ifndef ARSENAL_H
#define ARSENAL_H

#include<fstream>
#include<string>
#include<vector>
#include "Table2D.h"

using namespace std;

double Simpson_sum(double* , int, double);
//double interpolation2D_bilinear(Table2D* , double , double , int);

vector< vector<double>* >* readBlockData(istream &stream_in);
void releaseBlockData(vector< vector<double>* >* data);
vector<double> stringToDoubles(string str);

void
gauss (int npts, int job, double a, double b, double xpts[], double weights[]); //generate points and weights for gaussian quadrature, code is in gauss.cpp

double interpCubicDirect(vector<double>* x, vector<double>* y, double xx);
double interpCubicMono(vector<double>* x, vector<double>* y, double xx);
double interpLinearDirect(vector<double>* x, vector<double>* y, double xx);
double interpLinearMono(vector<double>* x, vector<double>* y, double xx);
double interpNearestDirect(vector<double>* x, vector<double>* y, double xx);
double interpNearestMono(vector<double>* x, vector<double>* y, double xx);

double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy=1e-10);

double invertTableDirect_hook(double xx);
double invertTableDirect(vector<double>* x, vector<double>* y, double y0, double x0, double relative_accuracy=1e-10);

long binarySearch(vector<double>* A, double value);

void outputFunctionerror(string function_name, string massage, double value, int level);

double cubeInterp(double x, double y, double z, double A000, double A100, double A010, double A110, double A001, double A101, double A011, double A111);

string toLower(string str);
string trim(string str);
double stringToDouble(string str);

#endif
