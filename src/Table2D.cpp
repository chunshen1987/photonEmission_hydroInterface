#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "stdlib.h"
#include <cmath>
#include "Arsenal.h"
#include "Table2D.h"

Table2D::Table2D() {};

Table2D::Table2D(std::string filename) { loadTableFromFile(filename); };

//----------------------------------------------------------------------
void Table2D::loadTableFromFile(std::string data_filename)
// The Table data file (data_filename) is assumed to be a n-column file:
{
    std::ostringstream filename_stream;
    filename_stream << data_filename;
    std::fstream fs(filename_stream.str().c_str());
    if (fs.is_open()==false) {
        std::cout << "Table2D::loadTableFromFile error: the data file cannot be opened."
                  << std::endl;
        std::cout << "Filename: " << data_filename << std::endl;
        exit(-1);
    }
    data = ARSENAL::readBlockData(fs);
    tb_sizeY = (*data).size();
    tb_sizeX = (*(*data)[0]).size();
}

void Table2D::outputTabletoFile(std::string filename) {
    std::ostringstream filename_stream;
    filename_stream << filename << ".dat";
    std::ofstream output(filename_stream.str().c_str());
    for (int i=0; i<tb_sizeX; i++) {
        for (int j=0; j<tb_sizeY; j++)
            output << std::scientific << std::setw(16) << std::setprecision(6)
                   << (*(*data)[j])[i] << "   ";
            output << std::endl;
    }
}
