#ifndef DATA_STRUCT_H_
#define DATA_STRUCT_H_

namespace PhysConsts {
    const double hbarC = 0.19733;
}


struct fluidCell {
   float ed, sd, temperature, pressure;
   float vx, vy, vz;
   float pi[4][4];
   float bulkPi;
};


class fluidCell_2D {
 public:
    float temperature = 0.;
    float ux = 0.;
    float uy = 0.;
    float ueta = 0.;
    // the shear stress tensor are already divided by e+P
    float pi00 = 0.;
    float pi01 = 0.;
    float pi02 = 0.;
    float pi11 = 0.;
    float pi12 = 0.;
    float pi22 = 0.;
    float pi33 = 0.;
    float bulkPi = 0.;

    fluidCell_2D operator + (fluidCell_2D const &obj) {
        fluidCell_2D res;
        res.temperature = temperature + obj.temperature;
        res.ux = ux + obj.ux;
        res.uy = uy + obj.uy;
        res.ueta = ueta + obj.ueta;
        res.pi00 = pi00 + obj.pi00;
        res.pi01 = pi01 + obj.pi01;
        res.pi02 = pi02 + obj.pi02;
        res.pi11 = pi11 + obj.pi11;
        res.pi12 = pi12 + obj.pi12;
        res.pi22 = pi22 + obj.pi22;
        res.pi33 = pi33 + obj.pi33;
        res.bulkPi = bulkPi + obj.bulkPi;
        return(res);
    }

    fluidCell_2D operator * (const double a) {
        fluidCell_2D res;
        res.temperature = temperature*a;
        res.ux = ux*a;
        res.uy = uy*a;
        res.ueta = ueta*a;
        res.pi00 = pi00*a;
        res.pi01 = pi01*a;
        res.pi02 = pi02*a;
        res.pi11 = pi11*a;
        res.pi12 = pi12*a;
        res.pi22 = pi22*a;
        res.pi33 = pi33*a;
        res.bulkPi = bulkPi*a;
        return(res);
    }
};


class fluidCell_3D_ideal {
 public:
    int itau = 0;
    int ix = 0;
    int iy = 0;
    int ieta = 0;
    float temperature = 0;
    float ed = 0;
    float pressure = 0;
    float ux = 0;
    float uy = 0;
    float uz = 0;

    fluidCell_3D_ideal operator + (fluidCell_3D_ideal const &obj) {
        fluidCell_3D_ideal res;
        res.temperature = temperature + obj.temperature;
        res.ed = ed + obj.ed;
        res.pressure = pressure + obj.pressure;
        res.ux = ux + obj.ux;
        res.uy = uy + obj.uy;
        res.uz = uz + obj.uz;
        return(res);
    }

    fluidCell_3D_ideal operator * (const double a) {
        fluidCell_3D_ideal res;
        res.temperature = temperature*a;
        res.ed = ed*a;
        res.pressure = pressure*a;
        res.ux = ux*a;
        res.uy = uy*a;
        res.uz = uz*a;
        return(res);
    }
};


class fluidCell_3D {
 public:
    float temperature = 0.;
    float ux = 0.;
    float uy = 0.;
    float ueta = 0.;
    // the shear stress tensor are already divided by e+P
    float pi00 = 0;
    float pi01 = 0;
    float pi02 = 0;
    float pi03 = 0;
    float pi11 = 0;
    float pi12 = 0;
    float pi13 = 0;
    float pi22 = 0;
    float pi23 = 0;
    float pi33 = 0;
    float bulkPi = 0;

    fluidCell_3D operator + (fluidCell_3D const &obj) {
        fluidCell_3D res;
        res.temperature = temperature + obj.temperature;
        res.ux = ux + obj.ux;
        res.uy = uy + obj.uy;
        res.ueta = ueta + obj.ueta;
        res.pi00 = pi00 + obj.pi00;
        res.pi01 = pi01 + obj.pi01;
        res.pi02 = pi02 + obj.pi02;
        res.pi03 = pi03 + obj.pi03;
        res.pi11 = pi11 + obj.pi11;
        res.pi12 = pi12 + obj.pi12;
        res.pi13 = pi13 + obj.pi13;
        res.pi22 = pi22 + obj.pi22;
        res.pi23 = pi23 + obj.pi23;
        res.pi33 = pi33 + obj.pi33;
        res.bulkPi = bulkPi + obj.bulkPi;
        return(res);
    }

    fluidCell_3D operator * (const double a) {
        fluidCell_3D res;
        res.temperature = temperature*a;
        res.ux = ux*a;
        res.uy = uy*a;
        res.ueta = ueta*a;
        res.pi00 = pi00*a;
        res.pi01 = pi01*a;
        res.pi02 = pi02*a;
        res.pi03 = pi03*a;
        res.pi11 = pi11*a;
        res.pi12 = pi12*a;
        res.pi13 = pi13*a;
        res.pi22 = pi22*a;
        res.pi23 = pi23*a;
        res.pi33 = pi33*a;
        res.bulkPi = bulkPi*a;
        return(res);
    }
};


class fluidCell_3D_new {
 public:
    int itau = 0;
    int ix = 0;
    int iy = 0;
    int ieta = 0;
    float temperature = 0.;
    float ux = 0;
    float uy = 0;
    float ueta = 0;
    float cs2 = 0;
    float muB = 0;
    // the shear stress tensor are already divided by e+P
    float pi11 = 0;
    float pi12 = 0;
    float pi13 = 0;
    float pi22 = 0;
    float pi23 = 0;
    float bulkPi = 0;

    fluidCell_3D_new operator + (fluidCell_3D_new const &obj) {
        fluidCell_3D_new res;
        res.temperature = temperature + obj.temperature;
        res.ux = ux + obj.ux;
        res.uy = uy + obj.uy;
        res.ueta = ueta + obj.ueta;
        res.cs2 = cs2 + obj.cs2;
        res.muB = muB + obj.muB;
        res.pi11 = pi11 + obj.pi11;
        res.pi12 = pi12 + obj.pi12;
        res.pi13 = pi13 + obj.pi13;
        res.pi22 = pi22 + obj.pi22;
        res.pi23 = pi23 + obj.pi23;
        res.bulkPi = bulkPi + obj.bulkPi;
        return(res);
    }

    fluidCell_3D_new operator * (const double a) {
        fluidCell_3D_new res;
        res.temperature = temperature*a;
        res.ux = ux*a;
        res.uy = uy*a;
        res.ueta = ueta*a;
        res.cs2 = cs2*a;
        res.muB = muB*a;
        res.pi11 = pi11*a;
        res.pi12 = pi12*a;
        res.pi13 = pi13*a;
        res.pi22 = pi22*a;
        res.pi23 = pi23*a;
        res.bulkPi = bulkPi*a;
        return(res);
    }
};


#endif  // DATA_STRUCT_H_
