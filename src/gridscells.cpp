#include <iostream>
#include <vector>
#include "vector3.h"
#include "gridscells.h"
#include "fieldsgrids.h"
#include "parameters.h"
#include "particles.h"
#include "module_base.h"
using std::vector;

// Constructors
GridsCells::GridsCells()
{
    gradPe = Vector3(0.0, 0.0, 0.0);
    gradPotential = Vector3(0.0, 0.0, 0.0);
    b3Cell = Vector3(0.0, 0.0, 0.0);
    particles_H = NULL;
    particles_He = NULL;
    particles_O = NULL;
    velDist_H = NULL;
    velDist_He = NULL;
    velDist_O = NULL;
    par_random = NULL;
    volume = 0.0;
}
GridsCells::GridsCells(double gradPe_x, double gradPe_y, double gradPe_z,
                       double gradPo_x, double gradPo_y, double gradPo_z,
                       double b3_x, double b3_y, double b3_z,
                       vector<Particles> part_H,
                       vector<Particles> part_He,
                       vector<Particles> part_O,
                       vector<vector<double>> velDist_HH,
                       vector<vector<double>> velDist_HEHE,
                       vector<vector<double>> velDist_OO,
                       double* par_r,
                       double vol)
{
    gradPe = Vector3(gradPe_x, gradPe_y, gradPe_z);
    gradPotential = Vector3(gradPo_x, gradPo_y, gradPo_z);
    b3Cell = Vector3(b3_x, b3_y, b3_z);
    particles_H = &part_H;
    particles_He = &part_He;
    particles_O = &part_O;
    //
    velDist_H = &velDist_HH;
    velDist_H = &velDist_HEHE;
    velDist_H = &velDist_OO;
    //
    par_random = par_r;
    volume = vol;
}