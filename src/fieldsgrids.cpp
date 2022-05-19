#include<iostream>
#include "vector3.h"
#include "fieldsgrids.h"
#include "parameters.h"

GridsPoints::GridsPoints( double px_in, double py_in, double pz_in,
                          double ex_in, double ey_in, double ez_in,
                          double dEx_in, double dEy_in, double dEz_in,
                          double bx_in, double by_in, double bz_in,
                          double dBx_in,double dBy_in,double dBz_in,
                          double vx_in, double vy_in, double vz_in,
                          double vHx_in, double vHy_in, double vHz_in,
                          double vHex_in, double vHey_in, double vHez_in,
                          double vOx_in, double vOy_in, double vOz_in,
                          double vHx_in_cumu, double vHy_in_cumu, double vHz_in_cumu,
                          double vHex_in_cumu, double vHey_in_cumu, double vHez_in_cumu,
                          double vOx_in_cumu, double vOy_in_cumu, double vOz_in_cumu,
                          double vex_in, double vey_in, double vez_in,
                          double gradBx_in, double gradBy_in, double gradBz_in,
                          double gradPex_in, double gradPey_in, double gradPez_in,
                          double density_H_in, double density_He_in, double density_O_in,
                          double density_H_in_cumu, double density_He_in_cumu, double density_O_in_cumu,
                          double temperature_in,
                          double temperature_H_in,
                          double temperature_He_in,
                          double temperature_O_in,
                          double potential_in,
                          double psd_H_in,
                          int stopSign_in)
{
//    ex = ex_in; ey = ey_in; ez = ez_in;
//    density = density_in;
//    FIJKtoXYZ( face_in, gi_in, gj_in, gk_in); // calculate px py pz
//    XYZtoB( pos3.x(), pos3.y(), pos3.z()); // calculate bx, by, bz

    pos3 = Vector3( px_in, py_in, pz_in);
    e3 =   Vector3( ex_in, ey_in, ez_in);
    dE3 =  Vector3( dEx_in, dEy_in, dEz_in);
    b3 =   Vector3( bx_in, by_in, bz_in);
    dB3=   Vector3( dBx_in, dBy_in, dBz_in);
    v3 =   Vector3( vx_in, vy_in, vz_in);
    vH3=   Vector3( vHx_in, vHy_in, vHz_in);
    vHe3=  Vector3( vHex_in, vHey_in, vHez_in);
    vO3=   Vector3( vOx_in, vOy_in, vOz_in);
    vH3_cumu=   Vector3( vHx_in, vHy_in, vHz_in);
    vHe3_cumu=  Vector3( vHex_in, vHey_in, vHez_in);
    vO3_cumu=   Vector3( vOx_in, vOy_in, vOz_in);
    ve3=   Vector3( vex_in, vey_in, vez_in);
    gradB3 = Vector3( gradBx_in, gradBy_in, gradBz_in);
    gradPe = Vector3( gradPex_in, gradPey_in, gradPez_in);

    density_H = density_H_in;
    density_He= density_He_in;
    density_O = density_O_in;
    density_H_cumu = density_H_in;
    density_He_cumu= density_He_in;
    density_O_cumu = density_O_in;
    
    temperature = temperature_in;
    temperature_H = temperature_H_in;
    temperature_He = temperature_He_in;
    temperature_O = temperature_O_in;
    
    potential = potential_in;
    psd_H = psd_H_in;    
    stopSign = stopSign_in;
//    face = 0; gi = 0; gj = 0; gk =0;
}

GridsPoints::GridsPoints(const GridsPoints& other)
{
  //  pos3.SetVector3( other.pos3);
    pos3 = Vector3( other.pos3);
    e3 = Vector3( other.e3);
    dE3 = Vector3( other.dE3);
    b3 = Vector3( other.b3);
    dB3= Vector3( other.dB3);
    v3 = Vector3( other.v3);
    vH3 = Vector3( other.vH3);
    vHe3 = Vector3( other.vHe3);
    vO3 = Vector3( other.vO3);
    vH3_cumu = Vector3( other.vH3_cumu);
    vHe3_cumu = Vector3( other.vHe3_cumu);
    vO3_cumu = Vector3( other.vO3_cumu);
    ve3= Vector3( other.ve3);
    gradB3 = Vector3( other.gradB3);
    gradPe = Vector3( other.gradPe);

    density_H = other.density_H;
    density_He = other.density_He;
    density_O = other.density_O;
    density_H_cumu = other.density_H_cumu;
    density_He_cumu = other.density_He_cumu;
    density_O_cumu = other.density_O_cumu;

    temperature = other.temperature;
    temperature = other.temperature;
    temperature = other.temperature;
    temperature = other.temperature;
    potential = other.potential;
    psd_H = other.psd_H;

    stopSign = other.stopSign;
}

GridsPoints::GridsPoints()
{
    pos3 = Vector3( 0.0, 0.0, 0.0);
    e3 =   Vector3( 0.0, 0.0, 0.0);
    dE3 =   Vector3( 0.0, 0.0, 0.0);
    b3 =   Vector3( 0.0, 0.0, 0.0);
    dB3=   Vector3( 0.0, 0.0, 0.0);
    v3 =   Vector3( 0.0, 0.0, 0.0);
    vH3 =  Vector3( 0.0, 0.0, 0.0);
    vHe3 = Vector3( 0.0, 0.0, 0.0);
    vO3 =  Vector3( 0.0, 0.0, 0.0);
    vH3_cumu =  Vector3( 0.0, 0.0, 0.0);
    vHe3_cumu = Vector3( 0.0, 0.0, 0.0);
    vO3_cumu =  Vector3( 0.0, 0.0, 0.0);
    ve3=   Vector3( 0.0, 0.0, 0.0);
    gradB3= Vector3( 0.0, 0.0, 0.0);
    gradPe= Vector3( 0.0, 0.0, 0.0);

    density_H = 0.0;
    density_He= 0.0;
    density_O = 0.0;
    density_H_cumu = 0.0;
    density_He_cumu= 0.0;
    density_O_cumu = 0.0;

    temperature = 0.0;
    temperature_H = 0.0;
    temperature_He = 0.0;
    temperature_O = 0.0;
    potential = 0.0;
    psd_H = 0.0;

    stopSign= 0;
//    face = 0; gi = 0; gj = 0; gk =0;
}