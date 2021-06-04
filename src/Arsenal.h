#ifndef ARSENAL_H
#define ARSENAL_H

#include <fstream>
#include <string>
#include <vector>

namespace ARSENAL {

double Simpson_sum(double* , int, double);

std::vector< std::vector<double>* >* readBlockData(std::istream &stream_in);
void releaseBlockData(std::vector< std::vector<double>* >* data);
std::vector<double> stringToDoubles(std::string str);

double interpCubicDirect(std::vector<double>* x, std::vector<double>* y, double xx);
double interpCubicMono(std::vector<double>* x, std::vector<double>* y, double xx);
double interpLinearDirect(std::vector<double>* x, std::vector<double>* y, double xx);
double interpLinearMono(std::vector<double>* x, std::vector<double>* y, double xx);
double interpNearestDirect(std::vector<double>* x, std::vector<double>* y, double xx);
double interpNearestMono(std::vector<double>* x, std::vector<double>* y, double xx);

double invertFunc(double (*func)(double), double y, double xL, double xR, double dx, double x0, double relative_accuracy=1e-10);

double invertTableDirect_hook(double xx);
double invertTableDirect(std::vector<double>* x, std::vector<double>* y, double y0, double x0, double relative_accuracy=1e-10);

long binarySearch(std::vector<double>* A, double value);

void outputFunctionerror(std::string function_name, std::string massage, double value, int level);

std::string toLower(std::string str);
std::string trim(std::string str);
double stringToDouble(std::string str);

double** createA2DMatrix(const int n1, const int n2, const double init);

void deleteA2DMatrix(double **mat, const int n1);

double*** createA3DMatrix(const int n1, const int n2, const int n3, const double init);

void deleteA3DMatrix(double ***mat, const int n1, const int n2);

double***** createA5DMatrix(const int n1, const int n2, const int n3, const int n4, const int n5, const double init);

void deleteA5DMatrix(double *****mat, const int n1, const int n2, const int n3, const int n4);

};

#endif
