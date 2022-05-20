#ifndef _MODULE_H_1_
#define _MODULE_H_1_
#include <iostream>
#include "parameters.h"
//#include "module_1.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include <cmath>
#include <limits>
#include <bitset>
#include "module_base.h"
#include "gridscells.h"
#include "vector3.h"
#include <algorithm>
#include <omp.h>
#include <array>
#include<bits/stdc++.h>

using std::max;
using std::vector;
using std::array;
using namespace std;
//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position
void ParticlesLists(vector<Particles> &listsPtr_in,
                    GridsPoints *****ptrArray_in,
                    double ***ptrVolumeCellArray_in,
                    double mi0);

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position
// Generate lists of particles for bot and top region temp
void ParticlesListsTemp(vector<Particles> &listsPtrTemp_in,
                        GridsPoints *****ptrArray_in,
                        double ***ptrVolumeCellArray_in,
                        double mi0,
                        int ionType_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Transform the ubit64 to vector of position

Vector3 Uint64ToVector3(uint_64 intPos_in);

//************************************************************************
//************************************************************************
// FUNCTION create a new particles
// in cell (face, i, j, k).
// uint_64, pos3, vel_para, vel_perp
Particles CreateParticles(GridsPoints *****ptrArray_in, Particles &particles, double mi0);

//************************************************************************
//************************************************************************
// FUNCTION
// Initialize particles in each cell
void ParticlesInCellsInitial(GridsPoints *****ptrArrayGrids,
                             GridsCells ****ptrArrayCells,
                             int type);

//************************************************************************
//************************************************************************
// FUNCTION
// Reset particles at bottom boundary layer
void ParticlesAtBottomCells(GridsPoints *****ptrArrayGrids,
                            GridsCells ****ptrArrayCells);

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random Vector3
// in cell (face, i, j, k).
// Remind that there are limited position for particles to locate, we just
// need to random the position in smaller cells
uint_64 UniDisInCell(int face_in, int i_in, int j_in, int k_in);

//************************************************************************
//************************************************************************
// FUNCTION return maxwell distrubtion random double number
// Only select 0.10-0.90
inline double MaxwellDisV(double iKT, double bulkV, double mi0_in)
{
    double sigma = sqrt(iKT / mi0_in);
    double temp = 0.0;
    while (temp < 0.005 || temp > 0.995)
    {
        temp = dRand();
    }
    return erfinv(2.0 * temp - 1.0) * sqrt(2.0) * sigma + bulkV;
}


void MaxwellDisV(GridsPoints *****ptrArray_in, uint_64 posUint, double mi0, double temperature, Vector3 &vp, double &mu);
Vector3 MaxwellDisV(GridsPoints *****ptrArray_in, uint_64 posUint, double mi0);
Vector3 MaxwellDisV_1(GridsPoints *****ptrArray_in, uint_64 posUint, double mi0, int n);

//************************************************************************
//************************************************************************
// FUNCTION return maxwell distrubtion random double number of particle
// energy for magnetic moment/ adiabatic invarient 0.1~10 eV
// the range for random number is -1 ~ 1 and then negetive it to make the
// major particles is between 0.1~1 eV. should make sigma_in < 0.2
// Only select 0.25-0.75
double MaxwellDisEnergy(GridsPoints *****ptrArray_in, uint_64 posUint);

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random double number
inline double UniformDisX(double x1, double x2)
{
    if (x1 < x2)
        return dRand() * (x2 - x1) + x1;
    else
        return dRand() * (x1 - x2) + x2;
}

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random Vector3
// return a Vector3 starting at the end of v1 and pointing to v2
inline Vector3 UniformDisVector3(Vector3 v1, Vector3 v2)
{
    Vector3 temp = v2.MinusProduct(v1).ScaleProduct(dRand());
    return temp;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Update info in the grids due to the info of particles and related
// weighting
void UpdateInfoGrids(GridsPoints *****ptrArray_in,
                     vector<Particles> ptrParticlesList_H_in,
                     vector<Particles> ptrParticlesList_He_in,
                     vector<Particles> ptrParticlesList_O_in,
                     vector<Particles> ptrParticlesListTemp_H_in,
                     vector<Particles> ptrParticlesListTemp_He_in,
                     vector<Particles> ptrParticlesListTemp_O_in,
                     double ***ptrVolumeGridArray_in,
                     int timeline_in, int updateInfoPeriod_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Go through particles vectors in temp domain
// 1. Boris method
// 2. Insert particles into main domain
//************************************************************************
//************************************************************************
void IterateParticlesTemp(GridsPoints *****ptrArray_in,
                          vector<Particles> &ptrParticlesList_in,
                          vector<Particles> &ptrParticlesListTemp_in,
                          vector<Particles> &ptrParticlesListTempBackup_in,
                          vector<int> &ptrParticlesList_out_in,
                          double ***ptrVolumeWeightGridArray,
                          Vector3 ******ptrVelWeightGridsArray,
                          double ******ptrMassWeightGridsArray,
                          double mi0_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Go through particles vectors in main domain
// 1. Update grids info
// 2. Boris method
// 3. Generate vector which shows the locations of the particles in 2
//************************************************************************
//************************************************************************
void IterateParticlesMain(GridsPoints *****ptrArray_in,
                          vector<Particles> &ptrParticlesList_in,
                          vector<int> &ptrParticlesList_out_in,
                          double ***ptrVolumeWeightGridArray,
                          Vector3 ******ptrVelWeightGridsArray,
                          double ******ptrMassWeightGridsArray,
                          double mi0_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Short for updating the weighting of particles on grids
//************************************************************************
//************************************************************************
inline void UpdateDueToWgtDown(GridsPoints *****ptrArray_in,
                               struct structg tempStr,
                               double tempNumber,
                               Vector3 tempVel,
                               double mi0_in)
{
    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);

    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);

    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);

    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
};

inline void UpdateDueToWgtUp(GridsPoints *****ptrArray_in,
                             struct structg tempStr,
                             double tempNumber,
                             Vector3 tempVel,
                             double mi0_in)
{
    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, mi0_in);

    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, mi0_in);

    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, mi0_in);

    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, mi0_in);
};

inline void UpdateDueToWgt(Vector3 ******ptrVelWeightGridsArray,
                           double ******ptrMassWeightGridsArray,
                           struct structg tempStr,
                           double number_in,
                           Vector3 vp_in,
                           double mi0_in,
                           int ipos,
                           int jpos,
                           int kpos,
                           int posP,
                           double w000,
                           double w100,
                           double w010,
                           double w110,
                           double w001,
                           double w101,
                           double w011,
                           double w111)
{
    int i = tempStr.ig - fieldsGridsSize / 2 * ipos;
    int j = tempStr.jg - fieldsGridsSize / 2 * jpos;
    int k = tempStr.kg - fieldsGridsSize / 2 * kpos;
    int ionType_in;
    if (mi0_in == mi0_H)
        ionType_in = 1;
    if (mi0_in == mi0_He)
        ionType_in = 4;
    if (mi0_in == mi0_O)
        ionType_in = 16;

    switch (ionType_in)
    {
    case 1:
    {
        ptrMassWeightGridsArray[0][tempStr.face][posP][i][j][k] += number_in * w000; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i][j][k] = ptrVelWeightGridsArray[0][tempStr.face][posP][i][j][k].PlusProduct(Vector3(number_in * vp_in.x() * w000, number_in * vp_in.y() * w000, number_in * vp_in.z() * w000));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i + 1][j][k] += number_in * w100; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j][k] = ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j][k].PlusProduct(Vector3(number_in * vp_in.x() * w100, number_in * vp_in.y() * w100, number_in * vp_in.z() * w100));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i][j + 1][k] += number_in * w010; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i][j + 1][k] = ptrVelWeightGridsArray[0][tempStr.face][posP][i][j + 1][k].PlusProduct(Vector3(number_in * vp_in.x() * w010, number_in * vp_in.y() * w010, number_in * vp_in.z() * w010));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i][j][k + 1] += number_in * w001; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i][j][k + 1] = ptrVelWeightGridsArray[0][tempStr.face][posP][i][j][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w001, number_in * vp_in.y() * w001, number_in * vp_in.z() * w001));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i + 1][j + 1][k] += number_in * w110; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j + 1][k] = ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j + 1][k].PlusProduct(Vector3(number_in * vp_in.x() * w110, number_in * vp_in.y() * w110, number_in * vp_in.z() * w110));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i + 1][j][k + 1] += number_in * w101; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j][k + 1] = ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w101, number_in * vp_in.y() * w101, number_in * vp_in.z() * w101));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i][j + 1][k + 1] += number_in * w011; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i][j + 1][k + 1] = ptrVelWeightGridsArray[0][tempStr.face][posP][i][j + 1][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w011, number_in * vp_in.y() * w011, number_in * vp_in.z() * w011));

        ptrMassWeightGridsArray[0][tempStr.face][posP][i + 1][j + 1][k + 1] += number_in * w111; // acutally is number not number density
        ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j + 1][k + 1] = ptrVelWeightGridsArray[0][tempStr.face][posP][i + 1][j + 1][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w111, number_in * vp_in.y() * w111, number_in * vp_in.z() * w111));

        break;
    }
    case 4:
    {
        ptrMassWeightGridsArray[1][tempStr.face][posP][i][j][k] += number_in * w000; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i][j][k] = ptrVelWeightGridsArray[1][tempStr.face][posP][i][j][k].PlusProduct(Vector3(number_in * vp_in.x() * w000, number_in * vp_in.y() * w000, number_in * vp_in.z() * w000));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i + 1][j][k] += number_in * w100; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j][k] = ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j][k].PlusProduct(Vector3(number_in * vp_in.x() * w100, number_in * vp_in.y() * w100, number_in * vp_in.z() * w100));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i][j + 1][k] += number_in * w010; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i][j + 1][k] = ptrVelWeightGridsArray[1][tempStr.face][posP][i][j + 1][k].PlusProduct(Vector3(number_in * vp_in.x() * w010, number_in * vp_in.y() * w010, number_in * vp_in.z() * w010));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i][j][k + 1] += number_in * w001; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i][j][k + 1] = ptrVelWeightGridsArray[1][tempStr.face][posP][i][j][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w001, number_in * vp_in.y() * w001, number_in * vp_in.z() * w001));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i + 1][j + 1][k] += number_in * w110; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j + 1][k] = ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j + 1][k].PlusProduct(Vector3(number_in * vp_in.x() * w110, number_in * vp_in.y() * w110, number_in * vp_in.z() * w110));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i + 1][j][k + 1] += number_in * w101; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j][k + 1] = ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w101, number_in * vp_in.y() * w101, number_in * vp_in.z() * w101));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i][j + 1][k + 1] += number_in * w011; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i][j + 1][k + 1] = ptrVelWeightGridsArray[1][tempStr.face][posP][i][j + 1][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w011, number_in * vp_in.y() * w011, number_in * vp_in.z() * w011));

        ptrMassWeightGridsArray[1][tempStr.face][posP][i + 1][j + 1][k + 1] += number_in * w111; // acutally is number not number density
        ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j + 1][k + 1] = ptrVelWeightGridsArray[1][tempStr.face][posP][i + 1][j + 1][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w111, number_in * vp_in.y() * w111, number_in * vp_in.z() * w111));

        break;
    }
    default:
    {
        ptrMassWeightGridsArray[2][tempStr.face][posP][i][j][k] += number_in * w000; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i][j][k] = ptrVelWeightGridsArray[2][tempStr.face][posP][i][j][k].PlusProduct(Vector3(number_in * vp_in.x() * w000, number_in * vp_in.y() * w000, number_in * vp_in.z() * w000));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i + 1][j][k] += number_in * w100; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j][k] = ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j][k].PlusProduct(Vector3(number_in * vp_in.x() * w100, number_in * vp_in.y() * w100, number_in * vp_in.z() * w100));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i][j + 1][k] += number_in * w010; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i][j + 1][k] = ptrVelWeightGridsArray[2][tempStr.face][posP][i][j + 1][k].PlusProduct(Vector3(number_in * vp_in.x() * w010, number_in * vp_in.y() * w010, number_in * vp_in.z() * w010));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i][j][k + 1] += number_in * w001; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i][j][k + 1] = ptrVelWeightGridsArray[2][tempStr.face][posP][i][j][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w001, number_in * vp_in.y() * w001, number_in * vp_in.z() * w001));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i + 1][j + 1][k] += number_in * w110; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j + 1][k] = ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j + 1][k].PlusProduct(Vector3(number_in * vp_in.x() * w110, number_in * vp_in.y() * w110, number_in * vp_in.z() * w110));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i + 1][j][k + 1] += number_in * w101; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j][k + 1] = ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w101, number_in * vp_in.y() * w101, number_in * vp_in.z() * w101));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i][j + 1][k + 1] += number_in * w011; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i][j + 1][k + 1] = ptrVelWeightGridsArray[2][tempStr.face][posP][i][j + 1][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w011, number_in * vp_in.y() * w011, number_in * vp_in.z() * w011));

        ptrMassWeightGridsArray[2][tempStr.face][posP][i + 1][j + 1][k + 1] += number_in * w111; // acutally is number not number density
        ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j + 1][k + 1] = ptrVelWeightGridsArray[2][tempStr.face][posP][i + 1][j + 1][k + 1].PlusProduct(Vector3(number_in * vp_in.x() * w111, number_in * vp_in.y() * w111, number_in * vp_in.z() * w111));

        break;
    }
    }
}

// calculate dW from Gaussian distribution
inline double gau(double sigma)
{
    double temp = dRand();
    double dW = erfinv(2.0 * temp - 1.0) * sqrt(2.0) * sigma;
    return dW;
}

// calculate sigma
// ***************** collision functions ************************

//***********************************************************************
// Coulomb logarithm (Gven in Nanbu and Yonemura, JCP, 145, 639, 1998)
// vv- parallel velocity (m/s);
inline double LnAab(int a, int b, double Te, double ne, double Ta, double Tb, Vector3 &vva, Vector3 &vvb)
{
    //! I calculate including transfering unit of ne into cm-3
    //const double cab[3][3] = {{0.690010011, -1.450056153, -0.226280721},
    //                          {-1.450056153, -2.082578711, -1.612575082},
    //                          {-0.226280721, -1.612575082, -0.696284350}};
    //  ! kbm=3kb/m (O+; H+; He+)
    const double kbm[3] = {1547.777412, 24764.4386, 6191.10965};
    double gav2 = kbm[a] * Ta + kbm[b] * Tb + vva.MinusProduct(vvb).norm2();
    double lnAab = cab[a][b] + 0.50 * log(Te / ne) + log(gav2);
    return lnAab;
}

//***************************************************************
// Particles collisions
inline void OperationsInCells(GridsPoints *****ptrArrayGrids,
                       GridsCells ****ptrArrayCells,
                       int timeline)
{
    //double Te, ne, TH, THe, TO;
    Vector3 vH_para, vHe_para, vO_para;
    double lnAab;
    // index for cell, i, j, k
#pragma omp parallel for collapse(4) 
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 1; k < fieldsGridsSize-1; k++)
                {   
                    // get B3 at corners of cells, for faster access
                    vector<Vector3> b3atvertex(8);
                    b3atvertex[0] = ptrArrayGrids[f][i+ 1][j + 1][k]->B3();
                    b3atvertex[1] = ptrArrayGrids[f][i+ 2][j + 1][k]->B3();
                    b3atvertex[2] = ptrArrayGrids[f][i+ 2][j + 2][k]->B3();
                    b3atvertex[3] = ptrArrayGrids[f][i+ 1][j + 2][k]->B3();
                    b3atvertex[4] = ptrArrayGrids[f][i+ 1][j + 1][k + 1]->B3();
                    b3atvertex[5] = ptrArrayGrids[f][i+ 2][j + 1][k + 1]->B3();
                    b3atvertex[6] = ptrArrayGrids[f][i+ 2][j + 2][k + 1]->B3();
                    b3atvertex[7] = ptrArrayGrids[f][i+ 1][j + 2][k + 1]->B3();

                    //************************************************************
                    // new version of collision
                    // optimize collisions
                    // 1. split or combine
                    // 2. change to absolute velocity 
                    // 3. self collision
                    // 4. mutual collision
                    //************************************************************
                    // (1) split or combine particles && (2) random particles position in array
                    ptrArrayCells[f][i][j][k].SplitOrCombine(0);
                    ptrArrayCells[f][i][j][k].SplitOrCombine(1);
                    ptrArrayCells[f][i][j][k].SplitOrCombine(2);
                    // absolute velocity: trans velocity of guiding center to abs velocity
                    ptrArrayCells[f][i][j][k].TransAbsoluteVelocity(0,
                                                                    b3atvertex);
                    ptrArrayCells[f][i][j][k].TransAbsoluteVelocity(1,
                                                                    b3atvertex);
                    ptrArrayCells[f][i][j][k].TransAbsoluteVelocity(2,
                                                                    b3atvertex);
                    if (collision_particles_control == 1 && timeline % collision_perPeriod == 0 ) // check collisions off/on and operation period
                    {
                        // assume
                        lnAab = 10.0;
                        //// LnA, need species (a,b), Te, ne, Ta, Tb, va_parallel, vb_parallel
                        //Vec3 b3Cell = ptrArrayCells[f][i][j][k].B3Cell().NormalizedVector();
                        ////
                        //Te = ptrArrayGrids[f][i+1][j+1][k]->Temperature() + ptrArrayGrids[f][i+2][j+1][k]->Temperature()
                        //    +ptrArrayGrids[f][i+1][j+2][k]->Temperature() + ptrArrayGrids[f][i+2][j+2][k]->Temperature()
                        //    +ptrArrayGrids[f][i+1][j+1][k+1]->Temperature() + ptrArrayGrids[f][i+2][j+1][k+1]->Temperature()
                        //    +ptrArrayGrids[f][i+1][j+2][k+1]->Temperature() + ptrArrayGrids[f][i+2][j+2][k+1]->Temperature();
                        //Te = Te / 8.0;
                        //ne = ptrArrayGrids[f][i+1][j+1][k]->Density() + ptrArrayGrids[f][i+2][j+1][k]->Density()
                        //    +ptrArrayGrids[f][i+1][j+2][k]->Density() + ptrArrayGrids[f][i+2][j+2][k]->Density()
                        //    +ptrArrayGrids[f][i+1][j+1][k+1]->Density() + ptrArrayGrids[f][i+2][j+1][k+1]->Density()
                        //    +ptrArrayGrids[f][i+1][j+2][k+1]->Density() + ptrArrayGrids[f][i+2][j+2][k+1]->Density();
                        //ne = ne / 8.0;
                        ////  velocity in cell and Temperature
                        //vH_para=Vec3( ptrArrayGrids[f][i+1][j+1][k]->VelH3().x() + ptrArrayGrids[f][i+2][j+1][k]->VelH3().x()
                        //                +ptrArrayGrids[f][i+1][j+2][k]->VelH3().x() + ptrArrayGrids[f][i+2][j+2][k]->VelH3().x()
                        //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelH3().x() + ptrArrayGrids[f][i+2][j+1][k+1]->VelH3().x()
                        //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelH3().x() + ptrArrayGrids[f][i+2][j+2][k+1]->VelH3().x(),
                        //                ptrArrayGrids[f][i+1][j+1][k]->VelH3().y() + ptrArrayGrids[f][i+2][j+1][k]->VelH3().y()
                        //                +ptrArrayGrids[f][i+1][j+2][k]->VelH3().y() + ptrArrayGrids[f][i+2][j+2][k]->VelH3().y()
                        //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelH3().y() + ptrArrayGrids[f][i+2][j+1][k+1]->VelH3().y()
                        //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelH3().y() + ptrArrayGrids[f][i+2][j+2][k+1]->VelH3().y(),
                        //                ptrArrayGrids[f][i+1][j+1][k]->VelH3().z() + ptrArrayGrids[f][i+2][j+1][k]->VelH3().z()
                        //                +ptrArrayGrids[f][i+1][j+2][k]->VelH3().z() + ptrArrayGrids[f][i+2][j+2][k]->VelH3().z()
                        //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelH3().z() + ptrArrayGrids[f][i+2][j+1][k+1]->VelH3().z()
                        //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelH3().z() + ptrArrayGrids[f][i+2][j+2][k+1]->VelH3().z()
                        //                );
                        //vH_para = vH_para.ScaleProduct(0.125);
                        //TH = ptrArrayCells[f][i][j][k].TempColl(0, vH_para);
                        //vHe_para=Vec3( ptrArrayGrids[f][i+1][j+1][k]->VelHe3().x() + ptrArrayGrids[f][i+2][j+1][k]->VelHe3().x()
                        //                 +ptrArrayGrids[f][i+1][j+2][k]->VelHe3().x() + ptrArrayGrids[f][i+2][j+2][k]->VelHe3().x()
                        //                 +ptrArrayGrids[f][i+1][j+1][k+1]->VelHe3().x() + ptrArrayGrids[f][i+2][j+1][k+1]->VelHe3().x()
                        //                 +ptrArrayGrids[f][i+1][j+2][k+1]->VelHe3().x() + ptrArrayGrids[f][i+2][j+2][k+1]->VelHe3().x(),
                        //                 ptrArrayGrids[f][i+1][j+1][k]->VelHe3().y() + ptrArrayGrids[f][i+2][j+1][k]->VelHe3().y()
                        //                 +ptrArrayGrids[f][i+1][j+2][k]->VelHe3().y() + ptrArrayGrids[f][i+2][j+2][k]->VelHe3().y()
                        //                 +ptrArrayGrids[f][i+1][j+1][k+1]->VelHe3().y() + ptrArrayGrids[f][i+2][j+1][k+1]->VelHe3().y()
                        //                 +ptrArrayGrids[f][i+1][j+2][k+1]->VelHe3().y() + ptrArrayGrids[f][i+2][j+2][k+1]->VelHe3().y(),
                        //                 ptrArrayGrids[f][i+1][j+1][k]->VelHe3().z() + ptrArrayGrids[f][i+2][j+1][k]->VelHe3().z()
                        //                 +ptrArrayGrids[f][i+1][j+2][k]->VelHe3().z() + ptrArrayGrids[f][i+2][j+2][k]->VelHe3().z()
                        //                 +ptrArrayGrids[f][i+1][j+1][k+1]->VelHe3().z() + ptrArrayGrids[f][i+2][j+1][k+1]->VelHe3().z()
                        //                 +ptrArrayGrids[f][i+1][j+2][k+1]->VelHe3().z() + ptrArrayGrids[f][i+2][j+2][k+1]->VelHe3().z()
                        //                 );
                        //vHe_para = vHe_para.ScaleProduct(0.125);
                        //THe = ptrArrayCells[f][i][j][k].TempColl(1, vHe_para);
                        //vO_para=Vec3( ptrArrayGrids[f][i+1][j+1][k]->VelO3().x() + ptrArrayGrids[f][i+2][j+1][k]->VelO3().x()
                        //                +ptrArrayGrids[f][i+1][j+2][k]->VelO3().x() + ptrArrayGrids[f][i+2][j+2][k]->VelO3().x()
                        //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelO3().x() + ptrArrayGrids[f][i+2][j+1][k+1]->VelO3().x()
                        //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelO3().x() + ptrArrayGrids[f][i+2][j+2][k+1]->VelO3().x(),
                        //                ptrArrayGrids[f][i+1][j+1][k]->VelO3().y() + ptrArrayGrids[f][i+2][j+1][k]->VelO3().y()
                        //                +ptrArrayGrids[f][i+1][j+2][k]->VelO3().y() + ptrArrayGrids[f][i+2][j+2][k]->VelO3().y()
                        //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelO3().y() + ptrArrayGrids[f][i+2][j+1][k+1]->VelO3().y()
                        //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelO3().y() + ptrArrayGrids[f][i+2][j+2][k+1]->VelO3().y(),
                        //                ptrArrayGrids[f][i+1][j+1][k]->VelO3().z() + ptrArrayGrids[f][i+2][j+1][k]->VelO3().z()
                        //                +ptrArrayGrids[f][i+1][j+2][k]->VelO3().z() + ptrArrayGrids[f][i+2][j+2][k]->VelO3().z()
                        //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelO3().z() + ptrArrayGrids[f][i+2][j+1][k+1]->VelO3().z()
                        //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelO3().z() + ptrArrayGrids[f][i+2][j+2][k+1]->VelO3().z()
                        //                );
                        //vO_para = vO_para.ScaleProduct(0.125);
                        //TO = ptrArrayCells[f][i][j][k].TempColl(2, vO_para);
                        ////  parallel velocity
                        //vH_para = b3Cell.ScaleProduct( vH_para.DotProduct( b3Cell));
                        //vHe_para= b3Cell.ScaleProduct( vHe_para.DotProduct( b3Cell));
                        //vO_para = b3Cell.ScaleProduct( vO_para.DotProduct( b3Cell));
                        //  H - H
                        //lnAab = LnAab(0, 0, Te, ne, TH, TH, vH_para, vH_para);
                        ptrArrayCells[f][i][j][k].self_coll(0, lnAab);
                        // He - He
                        //lnAab = LnAab(1, 1, Te, ne, THe, THe, vHe_para, vHe_para);
                        ptrArrayCells[f][i][j][k].self_coll(1, lnAab);
                        // O - O
                        //lnAab = LnAab(2, 2, Te, ne, TO, TO, vO_para, vO_para);
                        ptrArrayCells[f][i][j][k].self_coll(2, lnAab);
                        // H - He
                        //lnAab = LnAab(0, 1, Te, ne, TH, THe, vH_para, vHe_para);
                        ptrArrayCells[f][i][j][k].mutual_coll(0, 1, lnAab);
                        // H - O
                        //lnAab = LnAab(0, 2, Te, ne, TH, TO, vH_para, vO_para);
                        ptrArrayCells[f][i][j][k].mutual_coll(0, 2, lnAab);
                        // He - O
                        //lnAab = LnAab(1, 2, Te, ne, THe, TO, vHe_para, vO_para);
                        ptrArrayCells[f][i][j][k].mutual_coll(1, 2, lnAab);
                    }
                    //
                    // guiding center velocity: trans abs velocity to guiding center velocity
                    ptrArrayCells[f][i][j][k].TransGuidingCenterVelocity(   0,
                                                                            b3atvertex);
                    ptrArrayCells[f][i][j][k].TransGuidingCenterVelocity(   1,
                                                                            b3atvertex);
                    ptrArrayCells[f][i][j][k].TransGuidingCenterVelocity(   2,
                                                                            b3atvertex);
                    // wave-particles interaction
                    ptrArrayCells[f][i][j][k].WaveParticlesInteraction( 0, tstep, b3atvertex);
                    ptrArrayCells[f][i][j][k].WaveParticlesInteraction( 1, tstep, b3atvertex);
                    ptrArrayCells[f][i][j][k].WaveParticlesInteraction( 2, tstep, b3atvertex);
                }
            }
        }
    }
}

//************************************************************************
// // Collision in each cells   (old version)
// inline void Cell_coll(GridsPoints *****ptrArrayGrids,
//                       GridsCells ****ptrArrayCells,
//                       int timeline)
// {
//     //double Te, ne, TH, THe, TO;
//     Vector3 vH_para, vHe_para, vO_para;
//     double lnAab;
//     // index for cell, i, j, k
// #pragma omp parallel for collapse(4)
//     for (int f = 0; f < totalFace; f++)
//     {
//         for (int i = 0; i < fieldsGridsSize; i++)
//         {
//             for (int j = 0; j < fieldsGridsSize; j++)
//             {
//                 for (int k = 0; k < fieldsGridsSize; k++)
//                 {
//                     // speciesRandom
//                     ptrArrayCells[f][i][j][k].SplitOrCombine(0);
//                     ptrArrayCells[f][i][j][k].SplitOrCombine(1);
//                     ptrArrayCells[f][i][j][k].SplitOrCombine(2);
//                 //    ptrArrayCells[f][i][j][k].SpeciesRandom(0);
//                 //    ptrArrayCells[f][i][j][k].SpeciesRandom(1);
//                 //    ptrArrayCells[f][i][j][k].SpeciesRandom(2);
//                     if (collision_particles_control == 1 && timeline % collision_perPeriod ==0) // collision is on
//                     {
//                         // assume
//                         lnAab = 10.0;
//                         //// LnA, need species (a,b), Te, ne, Ta, Tb, va_parallel, vb_parallel
//                         //Vector3 b3Cell = ptrArrayCells[f][i][j][k].B3Cell().NormalizedVector();
//                         ////
//                         //Te = ptrArrayGrids[f][i+1][j+1][k]->Temperature() + ptrArrayGrids[f][i+2][j+1][k]->Temperature()
//                         //    +ptrArrayGrids[f][i+1][j+2][k]->Temperature() + ptrArrayGrids[f][i+2][j+2][k]->Temperature()
//                         //    +ptrArrayGrids[f][i+1][j+1][k+1]->Temperature() + ptrArrayGrids[f][i+2][j+1][k+1]->Temperature()
//                         //    +ptrArrayGrids[f][i+1][j+2][k+1]->Temperature() + ptrArrayGrids[f][i+2][j+2][k+1]->Temperature();
//                         //Te = Te / 8.0;
//                         //ne = ptrArrayGrids[f][i+1][j+1][k]->Density() + ptrArrayGrids[f][i+2][j+1][k]->Density()
//                         //    +ptrArrayGrids[f][i+1][j+2][k]->Density() + ptrArrayGrids[f][i+2][j+2][k]->Density()
//                         //    +ptrArrayGrids[f][i+1][j+1][k+1]->Density() + ptrArrayGrids[f][i+2][j+1][k+1]->Density()
//                         //    +ptrArrayGrids[f][i+1][j+2][k+1]->Density() + ptrArrayGrids[f][i+2][j+2][k+1]->Density();
//                         //ne = ne / 8.0;
//                         ////  velocity in cell and Temperature
//                         //vH_para=Vector3( ptrArrayGrids[f][i+1][j+1][k]->VelH3().x() + ptrArrayGrids[f][i+2][j+1][k]->VelH3().x()
//                         //                +ptrArrayGrids[f][i+1][j+2][k]->VelH3().x() + ptrArrayGrids[f][i+2][j+2][k]->VelH3().x()
//                         //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelH3().x() + ptrArrayGrids[f][i+2][j+1][k+1]->VelH3().x()
//                         //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelH3().x() + ptrArrayGrids[f][i+2][j+2][k+1]->VelH3().x(),
//                         //                ptrArrayGrids[f][i+1][j+1][k]->VelH3().y() + ptrArrayGrids[f][i+2][j+1][k]->VelH3().y()
//                         //                +ptrArrayGrids[f][i+1][j+2][k]->VelH3().y() + ptrArrayGrids[f][i+2][j+2][k]->VelH3().y()
//                         //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelH3().y() + ptrArrayGrids[f][i+2][j+1][k+1]->VelH3().y()
//                         //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelH3().y() + ptrArrayGrids[f][i+2][j+2][k+1]->VelH3().y(),
//                         //                ptrArrayGrids[f][i+1][j+1][k]->VelH3().z() + ptrArrayGrids[f][i+2][j+1][k]->VelH3().z()
//                         //                +ptrArrayGrids[f][i+1][j+2][k]->VelH3().z() + ptrArrayGrids[f][i+2][j+2][k]->VelH3().z()
//                         //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelH3().z() + ptrArrayGrids[f][i+2][j+1][k+1]->VelH3().z()
//                         //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelH3().z() + ptrArrayGrids[f][i+2][j+2][k+1]->VelH3().z()
//                         //                );
//                         //vH_para = vH_para.ScaleProduct(0.125);
//                         //TH = ptrArrayCells[f][i][j][k].TempColl(0, vH_para);
//                         //vHe_para=Vector3( ptrArrayGrids[f][i+1][j+1][k]->VelHe3().x() + ptrArrayGrids[f][i+2][j+1][k]->VelHe3().x()
//                         //                 +ptrArrayGrids[f][i+1][j+2][k]->VelHe3().x() + ptrArrayGrids[f][i+2][j+2][k]->VelHe3().x()
//                         //                 +ptrArrayGrids[f][i+1][j+1][k+1]->VelHe3().x() + ptrArrayGrids[f][i+2][j+1][k+1]->VelHe3().x()
//                         //                 +ptrArrayGrids[f][i+1][j+2][k+1]->VelHe3().x() + ptrArrayGrids[f][i+2][j+2][k+1]->VelHe3().x(),
//                         //                 ptrArrayGrids[f][i+1][j+1][k]->VelHe3().y() + ptrArrayGrids[f][i+2][j+1][k]->VelHe3().y()
//                         //                 +ptrArrayGrids[f][i+1][j+2][k]->VelHe3().y() + ptrArrayGrids[f][i+2][j+2][k]->VelHe3().y()
//                         //                 +ptrArrayGrids[f][i+1][j+1][k+1]->VelHe3().y() + ptrArrayGrids[f][i+2][j+1][k+1]->VelHe3().y()
//                         //                 +ptrArrayGrids[f][i+1][j+2][k+1]->VelHe3().y() + ptrArrayGrids[f][i+2][j+2][k+1]->VelHe3().y(),
//                         //                 ptrArrayGrids[f][i+1][j+1][k]->VelHe3().z() + ptrArrayGrids[f][i+2][j+1][k]->VelHe3().z()
//                         //                 +ptrArrayGrids[f][i+1][j+2][k]->VelHe3().z() + ptrArrayGrids[f][i+2][j+2][k]->VelHe3().z()
//                         //                 +ptrArrayGrids[f][i+1][j+1][k+1]->VelHe3().z() + ptrArrayGrids[f][i+2][j+1][k+1]->VelHe3().z()
//                         //                 +ptrArrayGrids[f][i+1][j+2][k+1]->VelHe3().z() + ptrArrayGrids[f][i+2][j+2][k+1]->VelHe3().z()
//                         //                 );
//                         //vHe_para = vHe_para.ScaleProduct(0.125);
//                         //THe = ptrArrayCells[f][i][j][k].TempColl(1, vHe_para);
//                         //vO_para=Vector3( ptrArrayGrids[f][i+1][j+1][k]->VelO3().x() + ptrArrayGrids[f][i+2][j+1][k]->VelO3().x()
//                         //                +ptrArrayGrids[f][i+1][j+2][k]->VelO3().x() + ptrArrayGrids[f][i+2][j+2][k]->VelO3().x()
//                         //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelO3().x() + ptrArrayGrids[f][i+2][j+1][k+1]->VelO3().x()
//                         //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelO3().x() + ptrArrayGrids[f][i+2][j+2][k+1]->VelO3().x(),
//                         //                ptrArrayGrids[f][i+1][j+1][k]->VelO3().y() + ptrArrayGrids[f][i+2][j+1][k]->VelO3().y()
//                         //                +ptrArrayGrids[f][i+1][j+2][k]->VelO3().y() + ptrArrayGrids[f][i+2][j+2][k]->VelO3().y()
//                         //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelO3().y() + ptrArrayGrids[f][i+2][j+1][k+1]->VelO3().y()
//                         //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelO3().y() + ptrArrayGrids[f][i+2][j+2][k+1]->VelO3().y(),
//                         //                ptrArrayGrids[f][i+1][j+1][k]->VelO3().z() + ptrArrayGrids[f][i+2][j+1][k]->VelO3().z()
//                         //                +ptrArrayGrids[f][i+1][j+2][k]->VelO3().z() + ptrArrayGrids[f][i+2][j+2][k]->VelO3().z()
//                         //                +ptrArrayGrids[f][i+1][j+1][k+1]->VelO3().z() + ptrArrayGrids[f][i+2][j+1][k+1]->VelO3().z()
//                         //                +ptrArrayGrids[f][i+1][j+2][k+1]->VelO3().z() + ptrArrayGrids[f][i+2][j+2][k+1]->VelO3().z()
//                         //                );
//                         //vO_para = vO_para.ScaleProduct(0.125);
//                         //TO = ptrArrayCells[f][i][j][k].TempColl(2, vO_para);
//                         ////  parallel velocity
//                         //vH_para = b3Cell.ScaleProduct( vH_para.DotProduct( b3Cell));
//                         //vHe_para= b3Cell.ScaleProduct( vHe_para.DotProduct( b3Cell));
//                         //vO_para = b3Cell.ScaleProduct( vO_para.DotProduct( b3Cell));
//                         //  H - H
//                         //lnAab = LnAab(0, 0, Te, ne, TH, TH, vH_para, vH_para);
//                         ptrArrayCells[f][i][j][k].self_coll(0, lnAab, ptrArrayGrids);
//                         // He - He
//                         //lnAab = LnAab(1, 1, Te, ne, THe, THe, vHe_para, vHe_para);
//                         ptrArrayCells[f][i][j][k].self_coll(1, lnAab, ptrArrayGrids);
//                         // O - O
//                         //lnAab = LnAab(2, 2, Te, ne, TO, TO, vO_para, vO_para);
//                         ptrArrayCells[f][i][j][k].self_coll(2, lnAab, ptrArrayGrids);
//                         // H - He
//                         //lnAab = LnAab(0, 1, Te, ne, TH, THe, vH_para, vHe_para);
//                         ptrArrayCells[f][i][j][k].mutual_coll(0, 1, lnAab, ptrArrayGrids);
//                         // H - O
//                         //lnAab = LnAab(0, 2, Te, ne, TH, TO, vH_para, vO_para);
//                         ptrArrayCells[f][i][j][k].mutual_coll(0, 2, lnAab, ptrArrayGrids);
//                         // He - O
//                         //lnAab = LnAab(1, 2, Te, ne, THe, TO, vHe_para, vO_para);
//                         ptrArrayCells[f][i][j][k].mutual_coll(1, 2, lnAab, ptrArrayGrids);
//                     }
//                     // velDist
//                     ptrArrayCells[f][i][j][k].velDistArray(0, ptrArrayGrids);
//                     ptrArrayCells[f][i][j][k].velDistArray(1, ptrArrayGrids);
//                     ptrArrayCells[f][i][j][k].velDistArray(2, ptrArrayGrids);
//                 }
//             }
//         }
//     }
// }
// *****************************************************************
 //Move particles via boris method
 // old version, not good, only used in test
 inline void IterateParticlesInCells(GridsPoints *****ptrArrayGrids,
                                     GridsCells ****ptrArrayCells,
                                     int a)
 {
     vector<Particles> ptrPartemp;
     double mi0;
     if ( a == 0)
     {mi0 = mi0_H;
     }
     else if (a ==1)
     {mi0 = mi0_He;
     }
     else if (a ==2)
     {mi0 = mi0_O;
     }
     // index for cell
 #pragma omp parallel for collapse(4) firstprivate(a)
     for (int f = 0; f < totalFace; f++)
     {
         for (int i = 0; i < fieldsGridsSize; i++)
         {
             for (int j = 0; j < fieldsGridsSize; j++)
             {
                 for (int k = 0; k < fieldsGridsSize; k++)
                 {
                     // private for openmp
                     struct structg tempStr;
                     struct structPar tempStrPar;
                     vector<Particles> *ptrParticles = NULL;
                     //
                     if ( a == 0)
                     {
                         ptrParticles = ptrArrayCells[f][i][j][k].Particles_H();
                     }
                     else if (a ==1)
                     {
                         ptrParticles = ptrArrayCells[f][i][j][k].Particles_He();
                     }
                     else if (a ==2)
                     {
                         ptrParticles = ptrArrayCells[f][i][j][k].Particles_O();
                     }
                     else
                     {
                         std::cout << " particles Type error in ->IterateParticlesInCells" << endl;
                         exit(0);
                     }
                     for (auto iter = ptrParticles->begin(); iter != ptrParticles->end(); ++iter)
                     {
                         //
                         if (!iter->AliveParticle())
                             continue;
                         // std::cout << " OppoPosUint check " << std::endl;
                         // update grids info
                         tempStr = iter->InttoStrp1();
                         double tempNumber = iter->WeightNi();
                         Vector3 tempVel = iter->VelParticles();
                         StructPar(tempStr, tempStrPar);
                         ////
                         //int ipos = tempStr.ig >> (fieldsGridsLevel - 1) | 0;
                         //int kpos = tempStr.kg >> (fieldsGridsLevel - 1) | 0;
                         //int jpos = tempStr.jg >> (fieldsGridsLevel - 1) | 0;
                         //int posP = kpos * 4 + jpos * 2 + ipos;
 #pragma omp critical
                     //    if (tempStr.kg >= tempGridsCellLevelBot && tempStr.kg <= fieldsGridsSize * grid_domain - 1 - tempGridsCellLevelTop)
                         if (tempStrPar.kg >= 0 && tempStrPar.kg <= fieldsGridsSize * grid_domain - 1 )
                         {
                         //    if (tempStr.kg > tempGridsCellLevelBot)
                         if (tempStr.kg > 0)
                             {
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 1][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w000, tempNumber, tempVel, mi0);
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 1][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w100, tempNumber, tempVel, mi0);
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 2][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w010, tempNumber, tempVel, mi0);
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 2][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w110, tempNumber, tempVel, mi0);
                             }
                         //    if (tempStr.kg < fieldsGridsSize * grid_domain - tempGridsCellLevelTop - 1)
                         if (tempStr.kg < fieldsGridsSize * grid_domain - 1)
                             {
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 1][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w001, tempNumber, tempVel, mi0);
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 1][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w101, tempNumber, tempVel, mi0);
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 2][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w011, tempNumber, tempVel, mi0);
                                 ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 2][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w111, tempNumber, tempVel, mi0);
                             }
                         }
                         // Trans into temp Par list: in domain, not in current cell
                         int check = 0;
                         check = iter->BorisMethod(tempStrPar, ptrArrayGrids, mi0, 0);
                         if (check == 0) // still in simulation domain
                         {
                             tempStr = iter->InttoStrp1();
                             if (tempStr.face == f && tempStr.ig == i && tempStr.jg == j && tempStr.kg == k)     // in the same cell
                             {
                                 continue;
                             }
                             else    // in different cell
                             {
 #pragma omp critical
                                 ptrPartemp.push_back(*iter);
                                 // deal with the old particles, recreate a new particles in bottom layer
                                 if (k < tempGridsCellLevelBot)
                                 {
                                     Particles tempBackup = *iter;
                                     *iter = CreateParticles(ptrArrayGrids, tempBackup, mi0);
                                 }
                                 else
                                 {
                                     iter->SetOutParticles();
                                 }
                             }
                         } else if (check == 1) // out the simulation domain
                         {
                             iter->SetOutParticles();
                         }                  
                     }
                     //
                     ptrArrayCells[f][i][j][k].SpeciesRandom(a);
                 }
             }
         }
     }
     // iterator tempPar back to cells
     struct structg tempStr;
     for (auto iter = ptrPartemp.begin(); iter < ptrPartemp.end(); ++iter)
     {
         if (!iter->AliveParticle())
             continue;
         tempStr = iter->InttoStrp1();
         if (mi0 == mi0_H)
         {
             ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_H()->push_back(*iter);
         }
         else if (mi0 == mi0_He)
         {
             ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_He()->push_back(*iter);
         }
         else if (mi0 == mi0_O)
         {
             ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_O()->push_back(*iter);
         }
     }
     //
     ptrPartemp.clear();
     ptrPartemp.shrink_to_fit();
}

// *****************************************************************
// Average velDist in each cells
inline void AverVelDistInCells(GridsCells ****ptrArrayCells, int average_timestep)
{
    // f, i, j, k are index of cells
#pragma omp parallel for collapse(4)
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    ptrArrayCells[f][i][j][k].AverVelDisArray(average_timestep, 0);
                    ptrArrayCells[f][i][j][k].AverVelDisArray(average_timestep, 1);
                    ptrArrayCells[f][i][j][k].AverVelDisArray(average_timestep, 2);
                }
            }
        }
    }
}
// Set zero velDist
inline void ResetVelDistInCells(GridsCells ****ptrArrayCells)
{
#pragma omp parallel for collapse(4)
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    ptrArrayCells[f][i][j][k].ResetVelDistArray(1);
                }
            }
        }
    }
}

//*********************************************************************
inline void MoveParticles(GridsPoints *****ptrArrayGrids,
                   GridsCells ****ptrArrayCells,
                   int a)
{
    vector<Particles> ptrPartemp;
    ptrPartemp.reserve(1000000);
    double mi0;
    if (a == 0)
    {
        mi0 = mi0_H;
    }
    else if (a == 1)
    {
        mi0 = mi0_He;
    }
    else if (a == 2)
    {
        mi0 = mi0_O;
    }
    else {
        std::cout << " Particle type error at the beginning of ->MoveParticles: a = " << a << endl;
    }

    // f, i, j, k are index of cells
#pragma omp parallel for collapse(4) firstprivate(a)
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    // set up a temp gridspoints to store information of particles in cells
                    //vector<GridsPoints> gridsPointsAtVertex(8);
                    //gridsPointsAtVertex[0] = *ptrArrayGrids[f][i+ 1][j + 1][k]; // w000
                    //gridsPointsAtVertex[1] = *ptrArrayGrids[f][i+ 2][j + 1][k]; // w100
                    //gridsPointsAtVertex[2] = *ptrArrayGrids[f][i+ 2][j + 2][k]; // w110
                    //gridsPointsAtVertex[3] = *ptrArrayGrids[f][i+ 1][j + 2][k]; // w010
                    //gridsPointsAtVertex[4] = *ptrArrayGrids[f][i+ 1][j + 1][k + 1]; // w001
                    //gridsPointsAtVertex[5] = *ptrArrayGrids[f][i+ 2][j + 1][k + 1]; // w101
                    //gridsPointsAtVertex[6] = *ptrArrayGrids[f][i+ 2][j + 2][k + 1]; // w111
                    //gridsPointsAtVertex[7] = *ptrArrayGrids[f][i+ 1][j + 2][k + 1]; // w011
                    //
                    struct structg tempStr = {0,0,0,0,0,0,0,0.0,0.0,0.0};
                    struct structPar tempStrPar = {0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                    vector<Particles> *ptrParticles = NULL;
                    //
                    if (a == 0)
                    {
                        ptrParticles = ptrArrayCells[f][i][j][k].Particles_H();
                    }
                    else if (a == 1)
                    {
                        ptrParticles = ptrArrayCells[f][i][j][k].Particles_He();
                    }
                    else if (a == 2)
                    {
                        ptrParticles = ptrArrayCells[f][i][j][k].Particles_O();
                    }
                    else
                    {
                        std::cout << " particles Type error in ->MoveParticles: a = " << a << endl;
                        exit(0);
                    }
                    
                    // normal version
                    for (auto iter = ptrParticles->begin(); iter != ptrParticles->end(); ++iter)
                    {
                        //
                        if (!iter->AliveParticle())
                            continue;
                        tempStr = iter->InttoStrp1();
                        StructPar(tempStr, tempStrPar);
                        // Trans into temp Par list: in domain, not in current cell
                        int check = 0;
                        check = iter->BorisMethod(tempStrPar, ptrArrayGrids, mi0, 0);
                        //   
                        if (check == 0) // still in simulation domain
                        {
                            tempStr = iter->InttoStrp1();   
                            if (tempStr.face == f && tempStr.ig == i && tempStr.jg == j && tempStr.kg == k) // in the same cell
                            {
                                continue;
                            }
                            else // in different cell
                            {
#pragma omp critical
                                ptrPartemp.push_back(*iter);
                                // deal with the old particles, recreate a new particles in bottom layer
                                if (k < tempGridsCellLevelBot)
                                {
                                    Particles tempBackup = *iter;
                                    *iter = CreateParticles(ptrArrayGrids, tempBackup, mi0);
                                }
                                else
                                {
                                    iter->SetOutParticles();
                                }
                            }
                        }
                        else if(check == 1) // out the simulation domain
                        {   
                            iter->SetOutParticles();
                        }
                    }
                    //
                    ptrArrayCells[f][i][j][k].SpeciesRearrange(a);
                    //
                    // // test version**********************************************************************
                    // for (auto iter = ptrParticles->begin(); iter != ptrParticles->end(); ++iter)
                    // {
                    //     //
                    //     if (!iter->AliveParticle())
                    //         continue;
                    //     Particles testPar_1 = *iter;
                    //     structg testStr_1, testStr_2, testStr_3;
                    //     testStr_1 = iter->InttoStrp1();
                    //     structPar testStrPar_1, testStrPar_2;
                    //     StructPar(testStr_1, testStrPar_1);
                    //     int check = 0;
                    //     check = iter->UpdateUint_64_test();
                    //     testStr_2 = iter->InttoStrp1();
                    //     StructPar(testStr_2, testStrPar_2);
                    //     Particles testPar_2 = *iter;
                    //     //
                    //     check = iter->UpdateUint_64_test();
                    //     testStr_3 = iter->InttoStrp1();
                    //     //
                    //     if(testStr_3.face != testStr_2.face || testStr_3.ig != testStr_2.ig || testStr_3.jg != testStr_2.jg || testStr_3.kg != testStr_2.kg)
                    //     {
                    //     std::cout << fixed << setprecision(16) << "test1 " 
                    //                 << iter->PosParticles().x() << " " << testPar_1.PosParticles().x() << " "
                    //                 << iter->PosParticles().y() << " " << testPar_1.PosParticles().y() << " "
                    //                 << iter->PosParticles().z() << " " << testPar_1.PosParticles().z() << " " 
                    //                 << testStr_2.face << " " << testStr_1.face << " " << k 
                    //                 << testStr_2.ig << " " << testStr_1.ig << " " << i
                    //                 << testStr_2.jg << " " << testStr_1.jg << " " << j
                    //                 << testStr_2.kg << " " << testStr_1.kg << " " << k
                    //                 << "\n"
                    //                 << std::bitset<64>(iter->PosUint()) << " \n" 
                    //                 << std::bitset<64>(testPar_2.PosUint()) << " \n" 
                    //                 << std::bitset<64>(testPar_1.PosUint()) << "\n";
                    //     std::cin.get();
                    //     }
                    //     if(testStr_1.face != testStr_2.face || testStr_1.ig != testStr_2.ig || testStr_1.jg != testStr_2.jg || testStr_1.kg != testStr_2.kg)
                    //     {
                    //     std::cout << fixed << setprecision(16) << "test2 " 
                    //                 << iter->PosParticles().x() << " " << testPar_1.PosParticles().x() << " "
                    //                 << iter->PosParticles().y() << " " << testPar_1.PosParticles().y() << " "
                    //                 << iter->PosParticles().z() << " " << testPar_1.PosParticles().z() << " " 
                    //                 << testStr_2.face << " " << testStr_1.face << " " << k << " i "
                    //                 << testStr_2.ig << " " << testStr_1.ig << " " << i << " j "
                    //                 << testStr_2.jg << " " << testStr_1.jg << " " << j << " k "
                    //                 << testStr_2.kg << " " << testStr_1.kg << " " << k
                    //                 << "\n"
                    //                 << std::bitset<64>(iter->PosUint()) << " \n" 
                    //                 << std::bitset<64>(testPar_2.PosUint()) << " \n" 
                    //                 << std::bitset<64>(testPar_1.PosUint()) << "\n";
                    //     std::cin.get();
                    //     }
                    // }
                }
            }
        }
    }
    // sequence version
    //// iterator tempPar back to cells
    //std::cout << " Num Particles moved out= " << npoints << " Species= " << a << "\n";
    //struct structg tempStr;
    //for (auto iter = ptrPartemp.begin(); iter < ptrPartemp.end(); ++iter)
    //{
    //    if (!iter->AliveParticle())
    //        continue;
    //    tempStr = iter->InttoStrp1();
    //    if (a==0)
    //    {
    //        ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_H()->push_back(*iter);
    //    }
    //    else if (a==1)
    //    {
    //        ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_He()->push_back(*iter);
    //    }
    //    else if (a==2)
    //    {
    //        ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_O()->push_back(*iter);
    //    }
    //}
    ////
    //ptrPartemp.clear();
    //ptrPartemp.shrink_to_fit();

    //  parallel version
    //  iterator tempPar back to cells
    int npoints = ptrPartemp.size();
    int thread_num, num_threads, sizeThread, startThread;
    //
    std::cout << " Num Particles moved out of original cells = " << npoints << " Species= " << a << "\n";
    //
    struct structg tempStr;
    #pragma omp parallel default(shared) private( thread_num, num_threads, sizeThread, startThread, tempStr)
    {
        thread_num = omp_get_thread_num();
        num_threads= omp_get_num_threads();
        sizeThread = npoints / num_threads;  // size for each thread
        startThread= thread_num * sizeThread;   // start point for each thread
        if( thread_num == num_threads-1)    // if last thread, it may have more
        sizeThread = npoints - startThread;
        //
        auto iter_zero = ptrPartemp.begin();
        auto iter_start = iter_zero + startThread;
        //
        for( auto iter = iter_start; iter < iter_start + sizeThread; ++iter)
        {
            if (!iter->AliveParticle())
                continue;
            tempStr = iter->InttoStrp1();
            //
            #pragma omp critical
            if (a==0)
            {
                ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_H()->push_back(*iter);
            }
            else if (a==1)
            {
                ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_He()->push_back(*iter);
            }
            else if (a==2)
            {
                ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_O()->push_back(*iter);
            }
        }
    #pragma omp barrier
    }
    //
    ptrPartemp.clear();
    ptrPartemp.shrink_to_fit();
}

//*****************************************************************************
inline void AssignParticles(GridsPoints *****ptrArrayGrids,
                     GridsCells ****ptrArrayCells,
                     int timeline,
                     int a)
{
    vector<Particles> ptrPartemp;
    /*double mi0;
    if (a == 0)
    {
        mi0 = mi0_H;
    }
    else if (a == 1)
    {
        mi0 = mi0_He;
    }
    else if (a == 2)
    {
        mi0 = mi0_O;
    }*/
    // index for cell
#pragma omp parallel for collapse(4) firstprivate(a)
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 1; k < fieldsGridsSize-1; k++)
                {
                    // set up a temp gridspoints to store information of particles in cells
                    vector<GridsPoints> gridsPointsAtVertex(8);
                    gridsPointsAtVertex[0] = *ptrArrayGrids[f][i+ 1][j + 1][k];
                    gridsPointsAtVertex[1] = *ptrArrayGrids[f][i+ 2][j + 1][k];
                    gridsPointsAtVertex[2] = *ptrArrayGrids[f][i+ 2][j + 2][k];
                    gridsPointsAtVertex[3] = *ptrArrayGrids[f][i+ 1][j + 2][k];
                    gridsPointsAtVertex[4] = *ptrArrayGrids[f][i+ 1][j + 1][k + 1];
                    gridsPointsAtVertex[5] = *ptrArrayGrids[f][i+ 2][j + 1][k + 1];
                    gridsPointsAtVertex[6] = *ptrArrayGrids[f][i+ 2][j + 2][k + 1];
                    gridsPointsAtVertex[7] = *ptrArrayGrids[f][i+ 1][j + 2][k + 1];
                    // velDist array assign
                    if (timeline % updateVelDist == 0)
                        ptrArrayCells[f][i][j][k].velDistArray(a, ptrArrayGrids);
                    // private for openmp
                    struct structg tempStr;
                    struct structPar tempStrPar; 
                    tempStrPar.face=0;
                    tempStrPar.ig=0;
                    tempStrPar.jg=0;
                    tempStrPar.kg=0;
                    vector<Particles> *ptrParticles = NULL;
                    //
                    if (a == 0)
                    {
                        ptrParticles = ptrArrayCells[f][i][j][k].Particles_H();
                    }
                    else if (a == 1)
                    {
                        ptrParticles = ptrArrayCells[f][i][j][k].Particles_He();
                    }
                    else if (a == 2)
                    {
                        ptrParticles = ptrArrayCells[f][i][j][k].Particles_O();
                    }
                    else
                    {
                        std::cout << " particles Type error in ->AssignParticles: a = " << a << endl;
                        exit(0);
                    }
                    for(auto iter = ptrParticles->begin(); iter != ptrParticles->end(); ++iter)
                    {
                        //
                        if(!iter->AliveParticle())
                            continue;
                        // std::cout << " OppoPosUint check " << std::endl;
                        // update grids info
                        tempStr = iter->InttoStrp1();
                        // use normalized_N to avoid exceed limit of double range
                        // iter->WeightNi() may be similar to 10^25 
                        double tempNumber = iter->WeightNi() / normalized_N;
                        StructPar(tempStr, tempStrPar);
                        Vector3 tempVel = iter->VelParticles();
                        //
                        #pragma omp critical
                        //  if (tempStr.kg >= tempGridsCellLevelBot && tempStr.kg <= fieldsGridsSize * grid_domain - 1 - tempGridsCellLevelTop)
                        if(tempStrPar.kg > 0 && tempStrPar.kg < fieldsGridsSize * grid_domain )
                        {
                            //    if (tempStr.kg > tempGridsCellLevelBot)
                            if(tempStr.kg > 1)
                            {
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 1][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w000, tempNumber, tempVel, mi0);
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 1][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w100, tempNumber, tempVel, mi0);
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 2][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w010, tempNumber, tempVel, mi0);
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 2][tempStrPar.kg]->UpdateDueToWgt(tempStrPar.w110, tempNumber, tempVel, mi0);
                                gridsPointsAtVertex[0].UpdateDueToWgt(tempStrPar.w000, tempNumber, tempVel, a);
                                gridsPointsAtVertex[1].UpdateDueToWgt(tempStrPar.w100, tempNumber, tempVel, a);
                                gridsPointsAtVertex[2].UpdateDueToWgt(tempStrPar.w110, tempNumber, tempVel, a);
                                gridsPointsAtVertex[3].UpdateDueToWgt(tempStrPar.w010, tempNumber, tempVel, a);
                            }
                            //    if (tempStr.kg < fieldsGridsSize * grid_domain - tempGridsCellLevelTop - 1)
                            if(tempStr.kg < fieldsGridsSize * grid_domain - 1)
                            {
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 1][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w001, tempNumber, tempVel, mi0);
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 1][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w101, tempNumber, tempVel, mi0);
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 2][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w011, tempNumber, tempVel, mi0);
                                //ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 2][tempStrPar.kg + 1]->UpdateDueToWgt(tempStrPar.w111, tempNumber, tempVel, mi0);
                                gridsPointsAtVertex[4].UpdateDueToWgt(tempStrPar.w001, tempNumber, tempVel, a);
                                gridsPointsAtVertex[5].UpdateDueToWgt(tempStrPar.w101, tempNumber, tempVel, a);
                                gridsPointsAtVertex[6].UpdateDueToWgt(tempStrPar.w111, tempNumber, tempVel, a);
                                gridsPointsAtVertex[7].UpdateDueToWgt(tempStrPar.w011, tempNumber, tempVel, a);
                            }
                        }
                    }
                    // critical add to gridsPoints
                    #pragma omp critical
                    {
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 1][tempStrPar.kg]->UpdateDueToWgt(gridsPointsAtVertex[0], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 1][tempStrPar.kg]->UpdateDueToWgt(gridsPointsAtVertex[1], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 2][tempStrPar.kg]->UpdateDueToWgt(gridsPointsAtVertex[3], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 2][tempStrPar.kg]->UpdateDueToWgt(gridsPointsAtVertex[2], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 1][tempStrPar.kg + 1]->UpdateDueToWgt(gridsPointsAtVertex[4], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 1][tempStrPar.kg + 1]->UpdateDueToWgt(gridsPointsAtVertex[5], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 1][tempStrPar.jg + 2][tempStrPar.kg + 1]->UpdateDueToWgt(gridsPointsAtVertex[7], a);
                        ptrArrayGrids[tempStrPar.face][tempStrPar.ig + 2][tempStrPar.jg + 2][tempStrPar.kg + 1]->UpdateDueToWgt(gridsPointsAtVertex[6], a);
                    }     
                }
            }
        }
    }
}

#endif