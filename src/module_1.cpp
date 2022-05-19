#include <iostream>
#include <list>
#include <memory>
#include <string>
#include "parameters.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module_1.h"
#include <cmath>
#include "H5Cpp.h"
#include <bitset>
#include <omp.h>
#include "module_base.h"
using std::cout;
using std::endl;
using std::make_shared;
using std::shared_ptr;


//************************************************************************
//************************************************************************
// FUNCTION
// Initialize particles in each cell
void ParticlesInCellsInitial(GridsPoints *****ptrArrayGrids,
                             GridsCells ****ptrArrayCells,
                             int type)
{ // F, I, J, K, index of cell
#pragma omp parallel for collapse(4)
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    if (k >= tempGridsCellLevelBot && type == 1)
                        continue;
                    // reset bottom particles
                    ptrArrayCells[f][i][j][k].Particles_H()->clear();
                    ptrArrayCells[f][i][j][k].Particles_He()->clear();
                    ptrArrayCells[f][i][j][k].Particles_O()->clear();
                    // number density
                    double N_H, N_He, N_O;
                    N_H = ptrArrayGrids[f][i + 1][j + 1][k]->Density_H() * 0.125 + ptrArrayGrids[f][i + 1][j + 1][k + 1]->Density_H() * 0.125 +
                          ptrArrayGrids[f][i + 2][j + 1][k]->Density_H() * 0.125 + ptrArrayGrids[f][i + 2][j + 1][k + 1]->Density_H() * 0.125 +
                          ptrArrayGrids[f][i + 1][j + 2][k]->Density_H() * 0.125 + ptrArrayGrids[f][i + 1][j + 2][k + 1]->Density_H() * 0.125 +
                          ptrArrayGrids[f][i + 2][j + 2][k]->Density_H() * 0.125 + ptrArrayGrids[f][i + 2][j + 2][k + 1]->Density_H() * 0.125;
                    N_He = ptrArrayGrids[f][i + 1][j + 1][k]->Density_He() * 0.125 + ptrArrayGrids[f][i + 1][j + 1][k + 1]->Density_He() * 0.125 +
                           ptrArrayGrids[f][i + 2][j + 1][k]->Density_He() * 0.125 + ptrArrayGrids[f][i + 2][j + 1][k + 1]->Density_He() * 0.125 +
                           ptrArrayGrids[f][i + 1][j + 2][k]->Density_He() * 0.125 + ptrArrayGrids[f][i + 1][j + 2][k + 1]->Density_He() * 0.125 +
                           ptrArrayGrids[f][i + 2][j + 2][k]->Density_He() * 0.125 + ptrArrayGrids[f][i + 2][j + 2][k + 1]->Density_He() * 0.125;
                    N_O = ptrArrayGrids[f][i + 1][j + 1][k]->Density_O() * 0.125 + ptrArrayGrids[f][i + 1][j + 1][k + 1]->Density_O() * 0.125 +
                          ptrArrayGrids[f][i + 2][j + 1][k]->Density_O() * 0.125 + ptrArrayGrids[f][i + 2][j + 1][k + 1]->Density_O() * 0.125 +
                          ptrArrayGrids[f][i + 1][j + 2][k]->Density_O() * 0.125 + ptrArrayGrids[f][i + 1][j + 2][k + 1]->Density_O() * 0.125 +
                          ptrArrayGrids[f][i + 2][j + 2][k]->Density_O() * 0.125 + ptrArrayGrids[f][i + 2][j + 2][k + 1]->Density_O() * 0.125;
                    // weightN of each simulation particles (count)
                    double N_H_simu = N_H / iniParticleNumberPerCell * ptrArrayCells[f][i][j][k].Volume();
                    double N_He_simu = N_He / iniParticleNumberPerCell * ptrArrayCells[f][i][j][k].Volume();
                    double N_O_simu = N_O / iniParticleNumberPerCell * ptrArrayCells[f][i][j][k].Volume();
                    //
                    uint_64 intPos;
                    Vector3 tempVector3, vVel;
                    double mu_simu;
                    Particles tempP;
                    // H
                    for (int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                        // calculate random position
                        intPos = UniDisInCell(f, i, j, k);
                        tempVector3 = Uint64ToVector3(intPos);
                        // Generate vp and mu
                        MaxwellDisV(ptrArrayGrids, intPos, mi0_H, 0.0, vVel, mu_simu);
                        // put the particles at the end of list
                        tempP = Particles(intPos, tempVector3, vVel, N_H_simu, mu_simu);
                        ptrArrayCells[f][i][j][k].Particles_H()->push_back(tempP);
                    }
                    // He
                    for (int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                        // calculate random position
                        intPos = UniDisInCell(f, i, j, k);
                        tempVector3 = Uint64ToVector3(intPos);
                        // Generate vp and mu
                        MaxwellDisV(ptrArrayGrids, intPos, mi0_He, 0.0, vVel, mu_simu);
                        // put the particles at the end of list
                        tempP = Particles(intPos, tempVector3, vVel, N_He_simu, mu_simu);
                        ptrArrayCells[f][i][j][k].Particles_He()->push_back(tempP);
                    }
                    // O
                    for (int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                        // calculate random position
                        intPos = UniDisInCell(f, i, j, k);
                        tempVector3 = Uint64ToVector3(intPos);
                        // Generate vp and mu
                        MaxwellDisV(ptrArrayGrids, intPos, mi0_O, 0.0, vVel, mu_simu);
                        // put the particles at the end of list
                        tempP = Particles(intPos, tempVector3, vVel, N_O_simu, mu_simu);
                        ptrArrayCells[f][i][j][k].Particles_O()->push_back(tempP);
                    }
                }
                //
            }
        }
    }
}

//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position
// Generate lists of particles
void ParticlesLists(GridsPoints *****ptrArray_in,
                    GridsCells ****ptrArrayCells,
                    double mi0)
{
    vector<Particles> *ptrPar;
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 1; i <= fieldsGridsSize; i++)
        {
            for (int j = 1; j <= fieldsGridsSize; j++)
            {
                for (int k = 1; k < fieldsGridsSize - 1; k++)
                {
                    //double r = ptrArray_in[f][i][j][k]->Pos3().norm() / radius;
                    double N;
                    if (mi0 == mi0_H)
                    {
                        N = ptrArray_in[f][i][j][k]->Density_H() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_H() * 0.125 +
                            ptrArray_in[f][i + 1][j][k]->Density_H() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_H() * 0.125 +
                            ptrArray_in[f][i][j + 1][k]->Density_H() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_H() * 0.125 +
                            ptrArray_in[f][i + 1][j + 1][k]->Density_H() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_H() * 0.125;
                        ptrPar = ptrArrayCells[f][i][j][k].Particles_H();
                    }
                    else if (mi0 == mi0_He)
                    {
                        N = ptrArray_in[f][i][j][k]->Density_He() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_He() * 0.125 +
                            ptrArray_in[f][i + 1][j][k]->Density_He() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_He() * 0.125 +
                            ptrArray_in[f][i][j + 1][k]->Density_He() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_He() * 0.125 +
                            ptrArray_in[f][i + 1][j + 1][k]->Density_He() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_He() * 0.125;
                        ptrPar = ptrArrayCells[f][i][j][k].Particles_He();
                    }
                    else //if (mi0 == mi0_O)
                    {
                        N = ptrArray_in[f][i][j][k]->Density_O() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_O() * 0.125 +
                            ptrArray_in[f][i + 1][j][k]->Density_O() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_O() * 0.125 +
                            ptrArray_in[f][i][j + 1][k]->Density_O() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_O() * 0.125 +
                            ptrArray_in[f][i + 1][j + 1][k]->Density_O() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_O() * 0.125;
                        ptrPar = ptrArrayCells[f][i][j][k].Particles_O();
                    }
                    // weightNi of each simulation particle
                    double Ni_simu = N / iniParticleNumberPerCell * ptrArrayCells[f][i][j][k].Volume();
                    //
                    for (int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                        // calculate random position
                        uint_64 intPos = UniDisInCell(f, i, j, k);
                        Vector3 tempVector3 = Uint64ToVector3(intPos);
                        // parallel vel && drift vel
                        Vector3 vVel = MaxwellDisV(ptrArray_in, intPos, mi0);
                        // perpendicular vel
                        double mu_simu = MaxwellDisEnergy(ptrArray_in, intPos);
                        // put the particles at the end of list
                        Particles tempP = Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
                        //
                        ptrPar->push_back(tempP);
                    }
                }
            }
        }
    }
}

//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position
// Generate lists of particles
void ParticlesLists(vector<Particles> &listsPtr_in,
                    GridsPoints *****ptrArray_in,
                    double ***ptrVolumeCellArray_in,
                    double mi0)
{
    //list<Particles> listP;
    //    auto listsPtr = make_shared<list<Particles>>();
    //    vector<Particles>* listsPtr = new vector<Particles>;
    //   shared_ptr<Particles> p1 = make_shared<Particles> (); // test
    //    double scaleHeight = ikT / mi0 / gravity; // assume the T is const for initiallization

    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 1; i <= fieldsGridsSize; i++)
        {
            for (int j = 1; j <= fieldsGridsSize; j++)
            {
                for (int k = 1; k < fieldsGridsSize - 1; k++)
                {
                    // number of real particles ( notice the unit is number density)
                    //double N = N0 * exp(-1.0 * (ptrArray_in[f][i][j][k]->Pos3().norm() - radius) / scaleHeight);
                    //double r = ptrArray_in[f][i][j][k]->Pos3().norm() / radius;
                    double N;
                    if (mi0 == mi0_H)
                    {
                        N = ptrArray_in[f][i][j][k]->Density_H() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_H() * 0.125 +
                            ptrArray_in[f][i + 1][j][k]->Density_H() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_H() * 0.125 +
                            ptrArray_in[f][i][j + 1][k]->Density_H() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_H() * 0.125 +
                            ptrArray_in[f][i + 1][j + 1][k]->Density_H() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_H() * 0.125;
                    }
                    else if (mi0 == mi0_He)
                    {
                        N = ptrArray_in[f][i][j][k]->Density_He() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_He() * 0.125 +
                            ptrArray_in[f][i + 1][j][k]->Density_He() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_He() * 0.125 +
                            ptrArray_in[f][i][j + 1][k]->Density_He() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_He() * 0.125 +
                            ptrArray_in[f][i + 1][j + 1][k]->Density_He() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_He() * 0.125;
                    }
                    else //if (mi0 == mi0_O)
                    {
                        N = ptrArray_in[f][i][j][k]->Density_O() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_O() * 0.125 +
                            ptrArray_in[f][i + 1][j][k]->Density_O() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_O() * 0.125 +
                            ptrArray_in[f][i][j + 1][k]->Density_O() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_O() * 0.125 +
                            ptrArray_in[f][i + 1][j + 1][k]->Density_O() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_O() * 0.125;
                    }
                    // weightNi of each simulation particle
                    double Ni_simu = N / iniParticleNumberPerCell * ptrVolumeCellArray_in[i][j][k];

                    for (int t = 1; t <= iniParticleNumberPerCell; t++)
                    {
                        // calculate random position
                        //test
                        /*
                    Vector3 temp1 = ptrArray_in[f][i][j][k]->Pos3();
                    Vector3 temp2 = ptrArray_in[f][i][j][k+1]->Pos3();
                    Vector3 temp = UniformDisVector3( ptrArray_in[f][i][j][k]->Pos3(),ptrArray_in[f][i][j][k+1]->Pos3());

                    Vector3 vPos = ptrArray_in[f][i][j][k]->Pos3().PlusProduct(UniformDisVector3( ptrArray_in[f][i][j][k]->Pos3(),ptrArray_in[f][i][j][k+1]->Pos3()));
                    vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[f][i][j][k]->Pos3(),ptrArray_in[f][i][j+1][k]->Pos3()));
                    vPos = vPos.PlusProduct( UniformDisVector3( ptrArray_in[f][i][j][k]->Pos3(),ptrArray_in[f][i+1][j][k]->Pos3()));
                    uint_64 intPos = vPos.Uint_64_Trans();
                    */
                        uint_64 intPos = UniDisInCell(f, i, j, k);
                        Vector3 tempVector3 = Uint64ToVector3(intPos);

                        //          std::cout << "1    " << std::bitset<64>(intPos) << std::endl;

                        // calculate random velocity
                        //           std::cout<< MaxwellDisV( ikT, ptrArray_in[f][i][j][k]->Vel3().x()) << std::endl;
                        //           std:: cout << ptrArray_in[f][i][j][k]->B3().norm()<< " B "<< std::endl;

                        // parallel vel && drift vel
                        Vector3 vVel = MaxwellDisV(ptrArray_in, intPos, mi0);

                        // perpendicular vel
                        double mu_simu = MaxwellDisEnergy(ptrArray_in, intPos);

                        // put the particles at the end of list
                        Particles tempP = Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
                        listsPtr_in.push_back(tempP);
                    }
                }
            }
        }
    }
    cout << "Particles initial size" << listsPtr_in.size() << endl;
    //    listsP.emplace_back(radius*3, radius*5, radius*1, 0.0, 0.0, 0.0);

    //    for( auto ptrL1 = listsPtr->begin(); ptrL1 != listsPtr->end(); ptrL1++)
    //   {  ptrL1->PrintP();  }
}

//************************************************************************
//************************************************************************
// FUNCTION
// Transform the ubit64 to vector of position
//************************************************************************
//************************************************************************
Vector3 Uint64ToVector3(uint_64 intPos_in)
{
    double px, py, pz;
    double temp[2];
    uint_64 posUint = intPos_in;
    uint_64 face = 0, ip = 0, jp = 0, kp = 0;

    face = posUint >> 61;
    for (int i = 0; i < particlesGridsLevel; i++)
    {
        ip = (ip << 1) + ((posUint >> (60 - i * 3)) & 1);
        jp = (jp << 1) + ((posUint >> (59 - i * 3)) & 1);
        kp = (kp << 1) + ((posUint >> (58 - i * 3)) & 1);
    }
    // 2. transfor to double x y z
    // 2.1 radial

    // cellSize1 = /particlesGridsSize * fieldsGridsSize
    // particles located at center of cell

    //    double L = LMin * pow(10, logRatio *  ( (kp +0.5)/ cellSize1 ));
    double L = 0.0;
    if (grid_domain == 1)
        L = LMin * pow(ratioK, kp + 0.5);
    //L = LMin * pow(10, logRatio * (kp + 0.5F));
    else if (grid_domain == 2)
    {
        if (kp < (uint_64)particlesGridsSize)
            L = LMin + (LMid - LMin) * sinh((kp + 0.5) / grid_N1) / const_sinh1;
        else if (kp >= (uint_64)particlesGridsSize)
            L = LMid + (LMax - LMid) * sinh((kp + 0.5) / grid_N2) / const_sinh2;
    }
    else
    {
        std::cout << " grid_domain error \n";
        std::cin.get();
    }
    //
    uint_64 testkp = static_cast<uint_64>(floor(log10(L / LMin) / logRatio));
    if (kp != testkp)
        std::cout << " checkkp " << testkp << " " << kp << "\n";
    //
    // 2.2 IgJg to ST note 0<ST<1
    temp[0] = (1.0 / particlesGridsSize) * (ip + 0.5);
    temp[1] = (1.0 / particlesGridsSize) * (jp + 0.5);
    // 2.3 ST to UV note -1<UV<1
    // tan(theta) = temp[i] = ST, note 0<ST<1
    temp[0] = tan(( temp[0] - 0.5) * pi / 2.0);
    temp[1] = tan(( temp[1] - 0.5) * pi / 2.0);
    // // // faster access version
    // for (int i = 0; i <= 1; i++)
    // {
    //     if (temp[i] >= 0.5)
    //     {
    //         temp[i] = (1.0 / 3.0) * (4.0 * temp[i] * temp[i] - 1.0);
    //     }
    //     else
    //     {
    //         temp[i] = (1.0 / 3.0) * (1.0 - 4.0 * (1.0 - temp[i]) * (1.0 - temp[i]));
    //     }
    // }
    // 2.4 UV to xyz
    double k = L * radius / sqrt(1.0 + temp[0] * temp[0] + temp[1] * temp[1]);
    switch (face)
    {
    case 0:
        px = 1.0;
        py = temp[0];
        pz = temp[1];
        break;
    case 1:
        px = -1.0 * temp[0];
        py = 1.0;
        pz = temp[1];
        break;
    case 2:
        px = -1.0 * temp[1];
        py = temp[0];
        pz = 1.0;
        break;
    case 3:
        px = -1.0;
        py = -1.0 * temp[0];
        pz = temp[1];
        break;
    case 4:
        px = temp[0];
        py = -1.0;
        pz = temp[1];
        break;
    default:
        px = temp[1];
        py = temp[0];
        pz = -1.0;
        break;
    }

    px *= k;
    py *= k;
    pz *= k;
    Vector3 tempV = Vector3(px, py, pz);
    return tempV;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Initialization the particles for velocity and position
// Generate lists of particles for bot and top region temp
//************************************************************************
//************************************************************************
void ParticlesListsTemp(vector<Particles> &listsPtrTemp_in,
                        GridsPoints *****ptrArray_in,
                        double ***ptrVolumeCellArray_in,
                        double mi0,
                        int ionType_in)
{
    double N, Ni_simu;
    int particleNumberPerCell=0;
    // for bottom temp domain
    // i, j, k are index of cells
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 1; i <= fieldsGridsSize; i++)
        {
            for (int j = 1; j <= fieldsGridsSize; j++)
            {
                if (update_type == 1)
                {
                    for (int k = 0; k < tempGridsCellLevelBot; k++)
                    {
                        // number density
                        switch (ionType_in)
                        {
                        case 1:
                            N = (ptrArray_in[f][i][j][k]->Density_H() + ptrArray_in[f][i + 1][j][k]->Density_H() + ptrArray_in[f][i][j + 1][k]->Density_H() + ptrArray_in[f][i + 1][j + 1][k]->Density_H()) * 0.25;
                            particleNumberPerCell = tempParticleNumberPerCellH;
                            break;
                        case 4:
                            N = (ptrArray_in[f][i][j][k]->Density_He() + ptrArray_in[f][i + 1][j][k]->Density_He() + ptrArray_in[f][i][j + 1][k]->Density_He() + ptrArray_in[f][i + 1][j + 1][k]->Density_He()) * 0.25;
                            particleNumberPerCell = tempParticleNumberPerCellHe;
                            break;
                        default:
                            N = (ptrArray_in[f][i][j][k]->Density_O() + ptrArray_in[f][i + 1][j][k]->Density_O() + ptrArray_in[f][i][j + 1][k]->Density_O() + ptrArray_in[f][i + 1][j + 1][k]->Density_O()) * 0.25;
                            particleNumberPerCell = tempParticleNumberPerCellO;
                            break;
                        }
                        // weight of each simulation particle
                        Ni_simu = N / particleNumberPerCell * ptrVolumeCellArray_in[i][j][k];
                        for (int t = 1; t <= particleNumberPerCell; t++)
                        {
                            // calculate random position
                            uint_64 intPos = UniDisInCell(f, i, j, k);
                            Vector3 tempVector3 = Uint64ToVector3(intPos);

                            // parallel vel && Drift vel
                            Vector3 vVel = MaxwellDisV(ptrArray_in, intPos, mi0);

                            // perpendicular vel
                            double mu_simu = MaxwellDisEnergy(ptrArray_in, intPos);

                            // put the particles at the end of list
                            Particles tempP = Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
                            listsPtrTemp_in.push_back(tempP);
                        }
                    }

                    for (int k = fieldsGridsSize * grid_domain - tempGridsCellLevelTop; k < fieldsGridsSize * grid_domain; k++)
                    {
                        // number density
                        switch (ionType_in)
                        {
                        case 1:
                            N = ptrArray_in[f][i][j][k]->Density_H() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_H() * 0.125 +
                                ptrArray_in[f][i + 1][j][k]->Density_H() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_H() * 0.125 +
                                ptrArray_in[f][i][j + 1][k]->Density_H() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_H() * 0.125 +
                                ptrArray_in[f][i + 1][j + 1][k]->Density_H() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_H() * 0.125;
                            particleNumberPerCell = tempParticleNumberPerCellH;
                            break;
                        case 4:
                            N = ptrArray_in[f][i][j][k]->Density_He() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_He() * 0.125 +
                                ptrArray_in[f][i + 1][j][k]->Density_He() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_He() * 0.125 +
                                ptrArray_in[f][i][j + 1][k]->Density_He() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_He() * 0.125 +
                                ptrArray_in[f][i + 1][j + 1][k]->Density_He() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_He() * 0.125;
                            particleNumberPerCell = tempParticleNumberPerCellHe;
                            break;
                        default:
                            N = ptrArray_in[f][i][j][k]->Density_O() * 0.125 + ptrArray_in[f][i][j][k + 1]->Density_O() * 0.125 +
                                ptrArray_in[f][i + 1][j][k]->Density_O() * 0.125 + ptrArray_in[f][i + 1][j][k + 1]->Density_O() * 0.125 +
                                ptrArray_in[f][i][j + 1][k]->Density_O() * 0.125 + ptrArray_in[f][i][j + 1][k + 1]->Density_O() * 0.125 +
                                ptrArray_in[f][i + 1][j + 1][k]->Density_O() * 0.125 + ptrArray_in[f][i + 1][j + 1][k + 1]->Density_O() * 0.125;
                            particleNumberPerCell = tempParticleNumberPerCellO;
                            break;
                        }

                        // weight of each simulation particle
                        // Ni_simu = N / tempParticleNumberPerCell *mi0 * ptrVolumeCellArray_in[i][j][k];
                        Ni_simu = N / particleNumberPerCell * ptrVolumeCellArray_in[i][j][k];

                        for (int t = 1; t <= particleNumberPerCell; t++)
                        {
                            // calculate random position
                            uint_64 intPos = UniDisInCell(f, i, j, k);
                            Vector3 tempVector3 = Uint64ToVector3(intPos);

                            // parallel vel
                            Vector3 vVel = MaxwellDisV(ptrArray_in, intPos, mi0);

                            // perpendicular vel
                            double mu_simu = MaxwellDisEnergy(ptrArray_in, intPos);
                            // put the particles at the end of list
                            Particles tempP = Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
                            listsPtrTemp_in.push_back(tempP);
                        }
                    }
                }
                else
                {
                    for (int k = 0; k < tempGridsCellLevelBot; k++)
                    {
                        // number density
                        switch (ionType_in)
                        {
                        case 1:
                            N = (ptrArray_in[f][i][j][k]->Density_H() + ptrArray_in[f][i + 1][j][k]->Density_H() + ptrArray_in[f][i][j + 1][k]->Density_H() + ptrArray_in[f][i + 1][j + 1][k]->Density_H()) * 0.25;
                            particleNumberPerCell = tempParticleNumberPerCellH;
                            break;
                        case 4:
                            N = (ptrArray_in[f][i][j][k]->Density_He() + ptrArray_in[f][i + 1][j][k]->Density_He() + ptrArray_in[f][i][j + 1][k]->Density_He() + ptrArray_in[f][i + 1][j + 1][k]->Density_He()) * 0.25;
                            particleNumberPerCell = tempParticleNumberPerCellHe;
                            break;
                        default:
                            N = (ptrArray_in[f][i][j][k]->Density_O() + ptrArray_in[f][i + 1][j][k]->Density_O() + ptrArray_in[f][i][j + 1][k]->Density_O() + ptrArray_in[f][i + 1][j + 1][k]->Density_O()) * 0.25;
                            particleNumberPerCell = tempParticleNumberPerCellO;
                            break;
                        }
                        // weight of each simulation particle
                        Ni_simu = N / particleNumberPerCell * ptrVolumeCellArray_in[i][j][k];
                        if (Ni_simu < 0.000001)
                        {
                            std::cout << " Ni_simu = 0 " << Ni_simu << std::endl;
                            std::cout << " face " << f << " " << i << " " << j << " " << k << " \n";
                        }
                        for (int t = 1; t <= particleNumberPerCell; t++)
                        {
                            // calculate random position
                            uint_64 intPos = UniDisInCell(f, i, j, k);
                            Vector3 tempVector3 = Uint64ToVector3(intPos);

                            // parallel vel && Drift vel
                            Vector3 vVel = MaxwellDisV(ptrArray_in, intPos, mi0);

                            // perpendicular vel
                            double mu_simu = MaxwellDisEnergy(ptrArray_in, intPos);

                            // put the particles at the end of list
                            Particles tempP = Particles(intPos, tempVector3, vVel, Ni_simu, mu_simu);
                            listsPtrTemp_in.push_back(tempP);
                        }
                    }
                }
            }
        }
    }

    cout << "ParticlesTemp initial size" << listsPtrTemp_in.size() << " " << particleNumberPerCell << endl;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Go through particles vectors in main domain
// 1. Update grids info
// 2. Boris method
// 3. Generate vector which shows the locations of the particles in 2
void IterateParticlesMain(GridsPoints *****ptrArray_in,
                          vector<Particles> &ptrParticlesList_in,
                          vector<int> &ptrParticlesList_out_in,
                          double ***ptrVolumeWeightGridArray,
                          Vector3 ******ptrVelWeightGridsArray,
                          double ******ptrMassWeightGridsArray,
                          double mi0_in)
{
    int npoints = ptrParticlesList_in.size();
    int thread_num, num_threads, sizeThread, startThread;

    //   omp_set_num_threads(1);
#pragma omp parallel default(shared) private(thread_num, num_threads, sizeThread, startThread)
    {
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();

        sizeThread = npoints / num_threads;    // size for each thread
        startThread = thread_num * sizeThread; // start point for each thread
        if (thread_num == num_threads - 1)     // if last thread, it may have more
            sizeThread = npoints - startThread;

        auto iter_zero = ptrParticlesList_in.begin();
        auto iter_start = iter_zero + startThread;

        omp_lock_t mylock;
        omp_init_lock(&mylock);

        for (auto iter = iter_start; iter < iter_start + sizeThread; iter++)
        {
            if (iter->PosUint() == 0)
                continue;
            // update grids info
            struct structg tempStr = iter->InttoStrp1();
            double tempNumber = iter->WeightNi();
            Vector3 tempVel = iter->VelParticles();

            // openmp method

            /*            double w000 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg+1][tempStr.kg+1];
            double w100 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg+1][tempStr.kg+1];
            double w010 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg][tempStr.kg+1];
            double w110 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg][tempStr.kg+1];
            double w001 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg+1][tempStr.kg];
            double w101 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg+1][tempStr.kg];
            double w011 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg][tempStr.kg];
            double w111 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg][tempStr.kg];


            double w000 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w100 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w010 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w110 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w001 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w101 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w011 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w111 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            
     
            double totalWeight = w000+ w100+ w010+ w110+ w001+ w101+ w011+ w111;
            w000 = w000/totalWeight;
            w100 = w100/totalWeight;
            w010 = w010/totalWeight;
            w110 = w110/totalWeight;
            w001 = w001/totalWeight;
            w101 = w101/totalWeight;
            w011 = w011/totalWeight;
            w111 = w111/totalWeight;

           */

            //            double w000 =1.0* (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw);
            //            double w100 =1.0* (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw);
            //            double w010 =1.0* (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw);
            //            double w110 =1.0* (1 +  tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw);
            //            double w001 =1.0* (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw);
            //            double w101 =1.0* (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw);
            //            double w011 =1.0* (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw);
            //            double w111 =1.0* (1 +  tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw);

            double w000 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (cellSize1 - tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
            double w100 = 1.0 * (0.5 + tempStr.iw) * (cellSize1 - tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
            double w010 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (1 + tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
            double w110 = 1.0 * (0.5 + tempStr.iw) * (0.5 + tempStr.jw) * (cellSize1 - tempStr.kw - 0.5);
            double w001 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (cellSize1 - tempStr.jw - 0.5) * (0.5 + tempStr.kw);
            double w101 = 1.0 * (0.5 + tempStr.iw) * (cellSize1 - tempStr.jw - 0.5) * (0.5 + tempStr.kw);
            double w011 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (0.5 + tempStr.jw) * (0.5 + tempStr.kw);
            double w111 = 1.0 * (0.5 + tempStr.iw) * (0.5 + tempStr.jw) * (0.5 + tempStr.kw);

            double totalWeight = w000 + w100 + w010 + w110 + w001 + w101 + w011 + w111;
            w000 = w000 / totalWeight;
            w100 = w100 / totalWeight;
            w010 = w010 / totalWeight;
            w110 = w110 / totalWeight;
            w001 = w001 / totalWeight;
            w101 = w101 / totalWeight;
            w011 = w011 / totalWeight;
            w111 = w111 / totalWeight;

            //           std::cout << w000 << " " << w100 << " " << w010 << " " << w001 << " " << w101 << " " << w110 << " " << w011 << " " <<w111 << std::endl;
            //           std::cout << tempStr.iw << " " <<  tempStr.jw << " " << tempStr.kw << std::endl;
            //           std::cin.get();

            //int ipos = tempStr.ig >> (fieldsGridsLevel - 1) | 0;
            //int kpos = tempStr.kg >> (fieldsGridsLevel - 1) | 0;
            //int jpos = tempStr.jg >> (fieldsGridsLevel - 1) | 0;
            //int posP = kpos * 4 + jpos * 2 + ipos;

#pragma omp critical
            if (tempStr.kg >= tempGridsCellLevelBot - 1 && tempStr.kg <= fieldsGridsSize * grid_domain - tempGridsCellLevelTop)
            {
                //  std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
                //  std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " "
                //  << "tempNumber "<< tempNumber << std::endl;
                /*
                if( tempStr.kg > tempGridsCellLevel - coverGridsCellLevel){
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                }
                if( tempStr.kg < fieldsGridsSize - tempGridsCellLevel + coverGridsCellLevel - 1){
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                }
            */

                if (tempStr.kg >= tempGridsCellLevelBot)
                {
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(w000, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(w100, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(w010, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(w110, tempNumber, tempVel, mi0_in);
                }
                if (tempStr.kg <= fieldsGridsSize * grid_domain - tempGridsCellLevelTop - 1)
                {
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(w001, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(w101, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(w011, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(w111, tempNumber, tempVel, mi0_in);
                }
            }

            //    #pragma omp barrier

            // boris method, update particle info
            // update vector_out_particles
            if (iter->PosUint() == 0)
                continue;
            else
            {
                int check = 0; // check whether in the main domain or not, "0" means in; "1" means out
                // update velocity // update position

                //            std::cout << "before " << iter->PosParticles().x() << " " << iter->VelParticles().x() ;

                //         check = iter->BorisMethod(&tempStr, ptrArray_in, mi0_in, 0);

                //            std::cout << " after " << iter->PosParticles().x() << " " << iter->VelParticles().x() << "\n";
                // check if still in the main domain
                if (check == 1) // out of the domain
                {
                    iter->SetOutParticles();
                    int tempint = iter - ptrParticlesList_in.begin();
#pragma omp critical
                    ptrParticlesList_out_in.push_back(tempint);
                }
            }
        }

#pragma omp barrier
    }
}

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
                          double mi0_in)
{

    int npoints = ptrParticlesListTemp_in.size();
    int thread_num, num_threads, sizeThread, startThread;

    //   omp_set_num_threads(1);
#pragma omp parallel default(shared) private(thread_num, num_threads, sizeThread, startThread)
    {
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        sizeThread = npoints / num_threads;    // size for each thread
        startThread = thread_num * sizeThread; // start point for each thread
        if (thread_num == num_threads - 1)     // if last thread, it may have more
            sizeThread = npoints - startThread;

        auto iter_zero = ptrParticlesListTemp_in.begin();
        auto iter_start = iter_zero + startThread;

        for (auto iter = iter_start; iter < iter_start + sizeThread; iter++)
        {
            if (iter->PosUint() == 0)
                continue;

            //Particles temp = *iter;
            struct structg tempStr = iter->InttoStrp1();
            double tempNumber = iter->WeightNi();
            Vector3 tempVel = iter->VelParticles();

            int loc = iter - iter_zero;
            Particles tempBackup = ptrParticlesListTempBackup_in[loc];

            // openmp method
            /*
            double w000 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg+1][tempStr.kg+1];
            double w100 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg+1][tempStr.kg+1];
            double w010 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg][tempStr.kg+1];
            double w110 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg][tempStr.kg+1];
            double w001 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg+1][tempStr.kg];
            double w101 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg+1][tempStr.kg];
            double w011 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig+1][tempStr.jg][tempStr.kg];
            double w111 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)*ptrVolumeWeightGridArray[tempStr.ig][tempStr.jg][tempStr.kg];

 
            
            int ipos = tempStr.ig >> (fieldsGridsLevel-1) | 0;
            int kpos = tempStr.kg >> (fieldsGridsLevel-1) | 0;
            int jpos = tempStr.jg >> (fieldsGridsLevel-1) | 0;
            int posP = kpos * 4 + jpos *2 + ipos;
                             

            double w000 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w100 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w010 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w110 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(cellSize1 - tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w001 = (cellSize1 - tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w101 = (1 + tempStr.iw)*(cellSize1 - tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w011 = (cellSize1 - tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
            double w111 = (1 +  tempStr.iw)*(1 + tempStr.jw)*(1 + tempStr.kw)/cellSize1/cellSize1/cellSize1;
     */

            double w000 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (cellSize1 - tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
            double w100 = 1.0 * (0.5 + tempStr.iw) * (cellSize1 - tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
            double w010 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (1 + tempStr.jw - 0.5) * (cellSize1 - tempStr.kw - 0.5);
            double w110 = 1.0 * (0.5 + tempStr.iw) * (0.5 + tempStr.jw) * (cellSize1 - tempStr.kw - 0.5);
            double w001 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (cellSize1 - tempStr.jw - 0.5) * (0.5 + tempStr.kw);
            double w101 = 1.0 * (0.5 + tempStr.iw) * (cellSize1 - tempStr.jw - 0.5) * (0.5 + tempStr.kw);
            double w011 = 1.0 * (cellSize1 - tempStr.iw - 0.5) * (0.5 + tempStr.jw) * (0.5 + tempStr.kw);
            double w111 = 1.0 * (0.5 + tempStr.iw) * (0.5 + tempStr.jw) * (0.5 + tempStr.kw);

            double totalWeight = w000 + w100 + w010 + w110 + w001 + w101 + w011 + w111;
            w000 = w000 / totalWeight;
            w100 = w100 / totalWeight;
            w010 = w010 / totalWeight;
            w110 = w110 / totalWeight;
            w001 = w001 / totalWeight;
            w101 = w101 / totalWeight;
            w011 = w011 / totalWeight;
            w111 = w111 / totalWeight;

            if (w000 < 0 || w000 < 0 || w000 < 0 || w000 < 0 || w000 < 0 || w000 < 0 || w000 < 0 || w000 < 0)
            {
                int pause;
                std::cout << "pause :";
                std::cin >> pause;
            }
#pragma omp critical
            if (tempStr.kg >= tempGridsCellLevelBot - 1 && tempStr.kg <= fieldsGridsSize * grid_domain - tempGridsCellLevelTop)
            {
                //  std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
                //  std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " "
                //  << "tempNumber "<< tempNumber << std::endl;
                if (tempStr.kg >= tempGridsCellLevelBot)
                {
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(w000, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(w100, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(w010, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(w110, tempNumber, tempVel, mi0_in);
                }
                if (tempStr.kg <= fieldsGridsSize * grid_domain - tempGridsCellLevelTop - 1)
                {
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(w001, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(w101, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(w011, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(w111, tempNumber, tempVel, mi0_in);
                }
            }

            /*
            #pragma omp critical
            if( tempStr.kg > tempGridsCellLevel - coverGridsCellLevel - 1 && tempStr.kg < fieldsGridsSize - tempGridsCellLevel + coverGridsCellLevel)
            {
            //  std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
            //  std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " " 
            //  << "tempNumber "<< tempNumber << std::endl;
                if( tempStr.kg > tempGridsCellLevel - coverGridsCellLevel){
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         cellSize1 - tempStr.kw, tempNumber, tempVel, mi0_in);
                }
                if( tempStr.kg < fieldsGridsSize - tempGridsCellLevel + coverGridsCellLevel - 1){
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         cellSize1 - tempStr.jw, tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( cellSize1 - tempStr.iw, tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                    ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->UpdateDueToWgt( tempStr.iw + 1,         tempStr.jw + 1,         tempStr.kw + 1,         tempNumber, tempVel, mi0_in);
                }
                
            }
*/
            // check whether in the main domain or not, "0" means in "1" means out
            // update velocity  // update position

            int check = 1; // in temp domain
            //Particles tempcheck = *iter;
            //    check = iter->BorisMethod(&tempStr, ptrArray_in, mi0_in, 1);

            //         if( check == 0)
            //         std::cout << iter->VelParticles().x() << " " << tempcheck.VelParticles().x() << iter->PosParticles().x() << " " << tempcheck.PosParticles().x() << " \n";

            // check if still in the main domain
            if (check == 0) // in the central domain, inject into particles array
            {
#pragma omp critical
                if (ptrParticlesList_out_in.size() > 0)
                {
                    auto temp_pos_out = ptrParticlesList_out_in.end() - 1;
                    ptrParticlesList_in[*temp_pos_out] = *iter;
                    temp_pos_out = ptrParticlesList_out_in.erase(temp_pos_out);
                }
                else
                {
                    ptrParticlesList_in.push_back(*iter);
                }
                //    iter->SetOutParticles();
                //         std::cout << "before " <<iter->VelParticles().x() << " " << tempcheck.VelParticles().x() << " " << tempBackup.VelParticles().x() << " \n";

                *iter = CreateParticles(ptrArray_in, tempBackup, mi0_in);

                //         std::cout << "after* " <<iter->VelParticles().x() << " " << tempcheck.VelParticles().x() << " " << tempBackup.VelParticles().x() << " \n";
            }
            else if (check == 2) // out of whole domain
            {

                *iter = CreateParticles(ptrArray_in, tempBackup, mi0_in);
            }
            else // random particles
            {
                double tempRandom = (double)rand() / (RAND_MAX);
                if (tempRandom < updateParticlesRate)
                {
                    Vector3 vVel = MaxwellDisV(ptrArray_in, iter->PosUint(), mi0_in);
                    double mu_simu = MaxwellDisEnergy(ptrArray_in, iter->PosUint());
                    iter->UpdateTempPar(vVel, mu_simu);
                }
            }
        }

#pragma omp barrier
    }
}

//************************************************************************
//************************************************************************
// FUNCTION
// Update info in the grids due to the info of particles and related
// weighting
// input: ptrArray_in for the grids address
//        ptrParticlesLists for the particles list address
//************************************************************************
//************************************************************************
void UpdateInfoGrids(GridsPoints *****ptrArray_in,
                     vector<Particles> ptrParticlesList_H_in,
                     vector<Particles> ptrParticlesList_He_in,
                     vector<Particles> ptrParticlesList_O_in,
                     vector<Particles> ptrParticlesListTemp_H_in,
                     vector<Particles> ptrParticlesListTemp_He_in,
                     vector<Particles> ptrParticlesListTemp_O_in,
                     double ***ptrVolumeGridArray_in,
                     int timeline_in, int updateInfoPeriod_in)
{
    if (timeline_in == 0 || (timeline_in - 1) % updateInfoPeriod_in == 0) // timeline_in should start with 1
    {
        // reset, clear the previous value: density and velocity
        for (int face = 0; face < totalFace; face++)
        {
            for (int i = 1; i < fieldsGridsSize + 2; i++)
            {
                for (int j = 1; j < fieldsGridsSize + 2; j++)
                {
                    for (int k = 1 + tempGridsCellLevelBot; k < fieldsGridsSize * grid_domain - tempGridsCellLevelTop; k++)
                    {
                        ptrArray_in[face][i][j][k]->ResetParameters();
                    }
                }
            }
        }
    }

    // For H particles in main domain
    for (auto iter = ptrParticlesList_H_in.begin(); iter != ptrParticlesList_H_in.end(); ++iter)
    {
        // get the weighting info of this particle
        struct structg tempStr = iter->InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = iter->WeightNi();

        //std::cout << temp.VelParticles().x() << " " << temp.VelParticles().y() << " " << temp.VelParticles().z() << std::endl;
        //std::cout << temp.WeightMi() << " <== " << std::endl;
        //int pause;
        //std::cin >>pause;
        //std::cout << "checkUp "<< tempStr.face << " " << tempStr.ig << " " << tempStr.jg << " " << tempStr.kg << " " << tempStr.iw << " " << tempStr.jw << " " << tempStr.kw << " " ;
        //std::cout << "str.v "<< tempStr.vx << " " << tempStr.vy << " " << tempStr.vz << " " ;
        //std::cout << "tempVel "<< tempVel.x() << " " << tempVel.y() << " " << tempVel.z() << " ";
        //std::cout << "tempNumber "<< tempNumber << std::endl;
        //int pause ;
        //std::cin >> pause;
        //std::cout << tempStr.face <<  tempStr.ig+1 << tempStr.jg+1 << tempStr.kg << " " << tempStr.iw << tempStr.jw << tempStr.kw << " mass " << tempNumber << " xxxx "; ;

        // get the info of velocity

        Vector3 tempVel = iter->VelParticles();
        // update density and velocity in the grids

        if (tempStr.kg > tempGridsCellLevelBot)
        {
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, 1);
        }
        if (tempStr.kg < fieldsGridsSize - 1 - tempGridsCellLevelTop)
        {
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, 1);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, 1);
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, 1);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, 1);
        }
        /*     
        std::cout << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg]->Density() << " nnnn "
        << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+1][tempStr.kg+1]->Density() << "  "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+1][tempStr.kg+1]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+1][tempStr.jg+2][tempStr.kg+1]->Density() << " "
        << ptrArray_in[tempStr.face][tempStr.ig+2][tempStr.jg+2][tempStr.kg+1]->Density() << std::endl;
        */
        // int pause ;
        // std::cin >> pause;
    }

    // For He particles in main domain
    for (auto iter = ptrParticlesList_He_in.begin(); iter != ptrParticlesList_He_in.end(); ++iter)
    {
        // get the weighting info of this particle
        struct structg tempStr = iter->InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = iter->WeightNi();
        // get the info of velocity
        Vector3 tempVel = iter->VelParticles();
        // update density and velocity in the grids
        if (tempStr.kg > tempGridsCellLevelBot)
        {
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, 4);
        }
        if (tempStr.kg < fieldsGridsSize - 1 - tempGridsCellLevelTop)
        {
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, 4);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, 4);
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, 4);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, 4);
        }
    }

    // For O particles in main domain
    for (auto iter = ptrParticlesList_O_in.begin(); iter != ptrParticlesList_O_in.end(); ++iter)
    {
        // get the weighting info of this particle
        struct structg tempStr = iter->InttoStrp1();
        // get the info of mass(weight of each simulation particle)
        double tempNumber = iter->WeightNi();
        // get the info of velocity
        Vector3 tempVel = iter->VelParticles();
        // update density and velocity in the grids
        if (tempStr.kg > tempGridsCellLevelBot)
        {
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, cellSize1 - tempStr.kw, tempNumber, tempVel, 16);
        }
        if (tempStr.kg < fieldsGridsSize - 1 - tempGridsCellLevelTop)
        {
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, 16);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 1][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, cellSize1 - tempStr.jw, tempStr.kw + 1, tempNumber, tempVel, 16);
            ptrArray_in[tempStr.face][tempStr.ig + 1][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(cellSize1 - tempStr.iw, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, 16);
            ptrArray_in[tempStr.face][tempStr.ig + 2][tempStr.jg + 2][tempStr.kg + 1]->UpdateDueToWgt(tempStr.iw + 1, tempStr.jw + 1, tempStr.kw + 1, tempNumber, tempVel, 16);
        }
    }

    // finish culmulating and average the density and velocity
    if (timeline_in % updateInfoPeriod_in == 0)
    {
        for (int face = 0; face < totalFace; face++)
        {
            for (int i = 1; i < fieldsGridsSize + 2; i++)
            {
                for (int j = 1; j < fieldsGridsSize + 2; j++)
                {
                    for (int k = 1 + tempGridsCellLevelBot; k < fieldsGridsSize - tempGridsCellLevelTop; k++)
                    {
                        //check stopsign
                        if (ptrArray_in[face][i][j][k]->StopSign() == 1)
                            continue;
                        // set volume
                        double volume = ptrVolumeGridArray_in[i - 1][j - 1][k]; // face of ptrArray is greater than that of ptrVolumeGridArray

                        //                    std::cout << face << i << j << k << " " ;
                        //                    std::cout << " volume " << volume << " density " << ptrArray_in[face][i][j][k]->Density() << " ==> " ;

                        ptrArray_in[face][i][j][k]->UpdateDueToWgt(ptrArray_in, volume, updateInfoPeriod_in);
                        // set stopSign
                        ptrArray_in[face][i][j][k]->SetStopSign(1);

                        //                    std::cout << ptrArray_in[face][i][j][k]->Density() << std::endl;
                    }
                }
            }
        }

#pragma omp parallel for collapse(4)
        // reset stopsign
        for (int face = 0; face < totalFace; face++)
        {
            for (int i = 1; i < fieldsGridsSize + 2; i++)
            {
                for (int j = 1; j < fieldsGridsSize + 2; j++)
                {
                    for (int k = 1; k < fieldsGridsSize; k++)
                    {
                        // set stopSign
                        ptrArray_in[face][i][j][k]->SetStopSign(0);
                    }
                }
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION create a new particles
// in cell (face, i, j, k).
// uint_64, pos3, vel_para, vel_perp
Particles CreateParticles(GridsPoints *****ptrArray_in, Particles &particles, double mi0)
{
    // 1. uint_64
    //    uint_64 parCell = static_cast<uint_64>( floor(cellSize1 *cellSize1 *cellSize1 *dRand()));
    double num = 1.0 * cellSize1 * cellSize1 * cellSize1 * dRand();

    uint_64 parCell = static_cast<uint_64>(floor(num));

    uint_64 pos64 = particles.PosUint();
    //Vector3 pVel = particles.VelParticles();

    //    std::cout << std::bitset<64>(parCell) << std::endl;
    //    std::cout << std::bitset<64>(pos64) <<  std::endl;

    //    pos64 = ((pos64 >> cellBitlength) << cellBitlength);

    //    std::cout << std::bitset<64>(pos64) <<  std::endl;

    pos64 = ((pos64 >> cellBitlength) << cellBitlength) + (parCell << blankBitlength);
    Vector3 tempVector3 = Uint64ToVector3(pos64);

    //    std::cout << tempVector3.norm() << "\n";
    // parallel vel && drift velocity
    //Vector3 vVel = MaxwellDisV( ptrArray_in, pos64, mi0);
    ////   std::cout << " velp "<< pVel.x() << " " << pVel.y() << " " << pVel.z() << std::endl;
    ////   std::cout << " velv "<< vVel.x() << " " << vVel.y() << " " << vVel.z() << std::endl;
    ////   std::cin.get();
    //// perpendicular vel
    //double mu_simu = MaxwellDisEnergy( ptrArray_in, pos64);
    //
    Vector3 vVel;
    double mu_simu;
    // weightNi of simulation particles
    double Ni_simu = particles.WeightNi();
    //
    MaxwellDisV(ptrArray_in, pos64, mi0, 0.0, vVel, mu_simu);
    Particles tempP = Particles(pos64, tempVector3, vVel, Ni_simu, mu_simu);

    return tempP;
}

//************************************************************************
//************************************************************************
// FUNCTION return uniform distrubtion random Vector3
// in cell (face, i, j, k).
// Remind that there are limited position for particles to locate, we just
// need to random the position in smaller cells
uint_64 UniDisInCell(int face_in, int i_in, int j_in, int k_in)
{
    // 1. Randomly get three ints in 3 direction
    //uint_64 ipar = static_cast<uint_64>(floor(cellSize1 * dRand()));
    //uint_64 jpar = static_cast<uint_64>(floor(cellSize1 * dRand()));
    //uint_64 kpar = static_cast<uint_64>(floor(cellSize1 * dRand()));
    //
    uint_64 ipar = static_cast<uint_64>(rand() % cellSize1);
    uint_64 jpar = static_cast<uint_64>(rand() % cellSize1);
    uint_64 kpar = static_cast<uint_64>(rand() % cellSize1);
    // 2. Transfer the base point from grids location for posUint calculation
    uint_64 face = face_in;
    uint_64 ig = i_in;
    uint_64 jg = j_in;
    uint_64 kg = k_in;

    //    std::cout << ipar << " " << jpar << " " << kpar << " " << std::endl;

    // 3. Transfer the location in int to the location in Uint_64
    uint_64 posUint = face << 61;
    for (int i = 0; i < fieldsGridsLevel; i++)
    {
        posUint += (((ig >> (fieldsGridsLevel - 1 - i)) & 1) << (60 - i * 3))
                + (((jg >> (fieldsGridsLevel - 1 - i)) & 1) << (60 - 1 - i * 3))
                + (((kg >> (fieldsGridsLevel - 1 - i) & 1)) << (60 - 2 - i * 3));
    }

    //    std::cout << "1    " << std::bitset<64>(posUint) << std::endl;
    int cellLevel = particlesGridsLevel - fieldsGridsLevel;
    for (int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        int j = i - fieldsGridsLevel;
        posUint += (((ipar >> (cellLevel - 1 - j)) & 1) << (60 - i * 3))
                + (((jpar >> (cellLevel - 1 - j)) & 1) << (60 - 1 - i * 3))
                + (((kpar >> (cellLevel - 1 - j)) & 1) << (60 - 2 - i * 3));
    }
    if (k_in >= fieldsGridsSize)
        posUint = posUint + 1;

    /*    std::cout << "2    " << std::bitset<64>(posUint) << std::endl; 

    struct structg strg = {0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 0.0};
    strg.face = posUint >> 61;
    for( int i = 0; i < fieldsGridsLevel; i++) 
    {
        strg.ig = (strg.ig << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jg = (strg.jg << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kg = (strg.kg << 1) + ((posUint >> 60-2 - i*3) & 1);
    }

    std::cout << "3    " << std::bitset<64>(posUint) << std::endl; 

    for( int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg.iw = (strg.iw << 1) + ((posUint >> 60   - i*3) & 1);
        strg.jw = (strg.jw << 1) + ((posUint >> 60-1 - i*3) & 1);
        strg.kw = (strg.kw << 1) + ((posUint >> 60-2 - i*3) & 1);
    }

    
   std::cout <<"face " << std::bitset<64>(strg.face) << std::endl;
    std::cout << "ig   " <<std::bitset<64>(strg.ig) << std::endl;
    std::cout <<"jg   " << std::bitset<64>(strg.jg) << std::endl;
    std::cout << "kg   " <<std::bitset<64>(strg.kg) << std::endl << std::endl;
    std::cout << "iw   " <<std::bitset<64>(strg.iw) << std::endl;
    std::cout <<"jw   " << std::bitset<64>(strg.jw) << std::endl;
    std::cout << "kw   " <<std::bitset<64>(strg.kw) << std::endl;

int pause ;
std::cin >> pause;
*/

    return posUint;
}

//************************************************************************
// Generate vp and mu
void MaxwellDisV(GridsPoints *****ptrArray_in, uint_64 posUint, double mi0, double temperature, Vector3 &vp, double &mu)
{
    struct structg strg_in = {0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
    strg_in.face = posUint >> 61;
    for (int i = 0; i < fieldsGridsLevel; i++)
    {
        strg_in.ig = (strg_in.ig << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jg = (strg_in.jg << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kg = (strg_in.kg << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }

    for (int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in.iw = (strg_in.iw << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jw = (strg_in.jw << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kw = (strg_in.kw << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }
    strg_in.vx = 0.0;
    strg_in.vy = 0.0;
    strg_in.vz = 0.0; // not need
                      //   strg_in.mass = 0.0; // not need
                      //************************************************************************
                      //   if( strg_in.kg >= fieldsGridsSize - tempGridsCellLevel)
                      //   strg_in.kg = strg_in.kg + tempGridsCellLevel - fieldsGridsSize;

    Vector3 tempb1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->B3();
    Vector3 tempb8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->B3();

    Vector3 tempe1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->E3();
    Vector3 tempe2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->E3();
    Vector3 tempe3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->E3();
    Vector3 tempe4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->E3();
    Vector3 tempe5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->E3();
    Vector3 tempe6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->E3();
    Vector3 tempe7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->E3();
    Vector3 tempe8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->E3();

    Vector3 tempv1;
    Vector3 tempv2;
    Vector3 tempv3;
    Vector3 tempv4;
    Vector3 tempv5;
    Vector3 tempv6;
    Vector3 tempv7;
    Vector3 tempv8;

    //    double w1 = 1.0- 1.0*(strg_in.iw +1) * (strg_in.jw +1) * (strg_in.kw +1) / cellSize1/ cellSize1/ cellSize1;
    //    double w2 = 1.0- 1.0*(cellSize1- strg_in.iw)* (strg_in.jw +1) * (strg_in.kw +1) / cellSize1/ cellSize1/ cellSize1;
    //    double w3 = 1.0- 1.0*(cellSize1- strg_in.iw) * (cellSize1- strg_in.jw) * (strg_in.kw +1) / cellSize1/ cellSize1/ cellSize1;
    //    double w4 = 1.0- 1.0*(strg_in.iw +1) * (cellSize1- strg_in.jw )* (strg_in.kw +1)/ cellSize1/ cellSize1/ cellSize1;
    //    double w5 = 1.0- 1.0*(strg_in.iw +1) * (strg_in.jw +1) * (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;
    //    double w6 = 1.0- 1.0*(cellSize1- strg_in.iw)* (strg_in.jw +1) * (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;
    //    double w7 = 1.0- 1.0*(cellSize1- strg_in.iw) * (cellSize1- strg_in.jw) * (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;
    //    double w8 = 1.0- 1.0*(strg_in.iw +1) * (cellSize1- strg_in.jw )* (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;

    double w1 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (strg_in.jw + 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w2 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (strg_in.jw + 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w3 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (cellSize1 - strg_in.jw - 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w4 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (cellSize1 - strg_in.jw - 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w5 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (strg_in.jw + 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w6 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (strg_in.jw + 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w7 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (cellSize1 - strg_in.jw - 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w8 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (cellSize1 - strg_in.jw - 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;

    Vector3 tempb;
    tempb.Setx(tempb1.x() * w1 + tempb2.x() * w2 + tempb3.x() * w3 + tempb4.x() * w4 + tempb5.x() * w5 + tempb6.x() * w6 + tempb7.x() * w7 + tempb8.x() * w8);

    tempb.Sety(tempb1.y() * w1 + tempb2.y() * w2 + tempb3.y() * w3 + tempb4.y() * w4 + tempb5.y() * w5 + tempb6.y() * w6 + tempb7.y() * w7 + tempb8.y() * w8);

    tempb.Setz(tempb1.z() * w1 + tempb2.z() * w2 + tempb3.z() * w3 + tempb4.z() * w4 + tempb5.z() * w5 + tempb6.z() * w6 + tempb7.z() * w7 + tempb8.z() * w8);

    Vector3 tempe;
    tempe.Setx(tempe1.x() * w1 + tempe2.x() * w2 + tempe3.x() * w3 + tempe4.x() * w4 + tempe5.x() * w5 + tempe6.x() * w6 + tempe7.x() * w7 + tempe8.x() * w8);

    tempe.Sety(tempe1.y() * w1 + tempe2.y() * w2 + tempe3.y() * w3 + tempe4.y() * w4 + tempe5.y() * w5 + tempe6.y() * w6 + tempe7.y() * w7 + tempe8.y() * w8);

    tempe.Setz(tempe1.z() * w1 + tempe2.z() * w2 + tempe3.z() * w3 + tempe4.z() * w4 + tempe5.z() * w5 + tempe6.z() * w6 + tempe7.z() * w7 + tempe8.z() * w8);

    Vector3 tempv;
    if (mi0 == mi0_H)
    {
        tempv1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->VelH3();
        tempv2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->VelH3();
        tempv3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->VelH3();
        tempv4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->VelH3();
        tempv5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->VelH3();
        tempv6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->VelH3();
        tempv7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->VelH3();
        tempv8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->VelH3();

        tempv.Setz(tempv1.z() * w1 + tempv2.z() * w2 + tempv3.z() * w3 + tempv4.z() * w4 + tempv5.z() * w5 + tempv6.z() * w6 + tempv7.z() * w7 + tempv8.z() * w8);
    }
    else if (mi0 == mi0_He)
    {
        tempv1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->VelHe3();
        tempv2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->VelHe3();
        tempv3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->VelHe3();
        tempv4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->VelHe3();
        tempv5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->VelHe3();
        tempv6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->VelHe3();
        tempv7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->VelHe3();
        tempv8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->VelHe3();

        tempv.Setz(tempv1.z() * w1 + tempv2.z() * w2 + tempv3.z() * w3 + tempv4.z() * w4 + tempv5.z() * w5 + tempv6.z() * w6 + tempv7.z() * w7 + tempv8.z() * w8);
    }
    else if (mi0 == mi0_O)
    {
        tempv1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->VelO3();
        tempv2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->VelO3();
        tempv3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->VelO3();
        tempv4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->VelO3();
        tempv5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->VelO3();
        tempv6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->VelO3();
        tempv7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->VelO3();
        tempv8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->VelO3();

        tempv.Setz(tempv1.z() * w1 + tempv2.z() * w2 + tempv3.z() * w3 + tempv4.z() * w4 + tempv5.z() * w5 + tempv6.z() * w6 + tempv7.z() * w7 + tempv8.z() * w8);
    }
    else
    {
        std::cout << " MaxwellDisV tempV error \n";
        std::cout << " MaxwellDisV/n";
        std::cin.get();
    }
    // tempv - bulk velocity of guiding center; parallel components- v_para; perp components- v_drift
    //Vector3 tempv_para = tempb.NormalizedVector().ScaleProduct(tempv.DotProduct(tempb.NormalizedVector()));
    //Vector3 tempv_drift= tempv.MinusProduct(tempv_para);
    // tempV - velocity of a particle in coordinates of guiding center;
    Vector3 tempV = Vector3(MaxwellDisV(ikT + boltzmann_k * temperature, 0.0, mi0),
                            MaxwellDisV(ikT + boltzmann_k * temperature, 0.0, mi0),
                            MaxwellDisV(ikT + boltzmann_k * temperature, 0.0, mi0));
    Vector3 tempV_para = tempb.NormalizedVector().ScaleProduct(tempV.DotProduct(tempb.NormalizedVector()));
    Vector3 tempV_gyro = tempV.MinusProduct(tempV_para);
    // vp has components: 1. perp: tempv_drift
    //                    2. para: tempv_para + tempV_para
    // mu has components into only one: tempV_gyro
    vp = tempv.PlusProduct(tempV_para); // vel of guiding center
    if (tempb.norm2() > 0.0)
        mu = tempV_gyro.norm2() * 0.5 * mi0 / tempb.norm();
    else
        mu = 0.0;
}

//********************************************************************************
Vector3 MaxwellDisV(GridsPoints *****ptrArray_in, uint_64 posUint, double mi0)
{
    struct structg strg_in = {0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
    strg_in.face = posUint >> 61;
    for (int i = 0; i < fieldsGridsLevel; i++)
    {
        strg_in.ig = (strg_in.ig << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jg = (strg_in.jg << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kg = (strg_in.kg << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }

    for (int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in.iw = (strg_in.iw << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jw = (strg_in.jw << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kw = (strg_in.kw << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }
    strg_in.vx = 0.0;
    strg_in.vy = 0.0;
    strg_in.vz = 0.0; // not need
                      //   strg_in.mass = 0.0; // not need
                      //************************************************************************
                      //   if( strg_in.kg >= fieldsGridsSize - tempGridsCellLevel)
                      //   strg_in.kg = strg_in.kg + tempGridsCellLevel - fieldsGridsSize;

    Vector3 tempb1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->B3();
    Vector3 tempb8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->B3();

    Vector3 tempe1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->E3();
    Vector3 tempe2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->E3();
    Vector3 tempe3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->E3();
    Vector3 tempe4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->E3();
    Vector3 tempe5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->E3();
    Vector3 tempe6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->E3();
    Vector3 tempe7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->E3();
    Vector3 tempe8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->E3();

    Vector3 tempv1;
    Vector3 tempv2;
    Vector3 tempv3;
    Vector3 tempv4;
    Vector3 tempv5;
    Vector3 tempv6;
    Vector3 tempv7;
    Vector3 tempv8;

    //    double w1 = 1.0- 1.0*(strg_in.iw +1) * (strg_in.jw +1) * (strg_in.kw +1) / cellSize1/ cellSize1/ cellSize1;
    //    double w2 = 1.0- 1.0*(cellSize1- strg_in.iw)* (strg_in.jw +1) * (strg_in.kw +1) / cellSize1/ cellSize1/ cellSize1;
    //    double w3 = 1.0- 1.0*(cellSize1- strg_in.iw) * (cellSize1- strg_in.jw) * (strg_in.kw +1) / cellSize1/ cellSize1/ cellSize1;
    //    double w4 = 1.0- 1.0*(strg_in.iw +1) * (cellSize1- strg_in.jw )* (strg_in.kw +1)/ cellSize1/ cellSize1/ cellSize1;
    //    double w5 = 1.0- 1.0*(strg_in.iw +1) * (strg_in.jw +1) * (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;
    //    double w6 = 1.0- 1.0*(cellSize1- strg_in.iw)* (strg_in.jw +1) * (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;
    //    double w7 = 1.0- 1.0*(cellSize1- strg_in.iw) * (cellSize1- strg_in.jw) * (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;
    //    double w8 = 1.0- 1.0*(strg_in.iw +1) * (cellSize1- strg_in.jw )* (cellSize1- strg_in.kw) / cellSize1/ cellSize1/ cellSize1;

    double w1 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (strg_in.jw + 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w2 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (strg_in.jw + 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w3 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (cellSize1 - strg_in.jw - 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w4 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (cellSize1 - strg_in.jw - 0.5) * (strg_in.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w5 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (strg_in.jw + 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w6 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (strg_in.jw + 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w7 = 1.0 - 1.0 * (cellSize1 - strg_in.iw - 0.5) * (cellSize1 - strg_in.jw - 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    double w8 = 1.0 - 1.0 * (strg_in.iw + 0.5) * (cellSize1 - strg_in.jw - 0.5) * (cellSize1 - strg_in.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;

    Vector3 tempb;
    tempb.Setx(tempb1.x() * w1 + tempb2.x() * w2 + tempb3.x() * w3 + tempb4.x() * w4 + tempb5.x() * w5 + tempb6.x() * w6 + tempb7.x() * w7 + tempb8.x() * w8);

    tempb.Sety(tempb1.y() * w1 + tempb2.y() * w2 + tempb3.y() * w3 + tempb4.y() * w4 + tempb5.y() * w5 + tempb6.y() * w6 + tempb7.y() * w7 + tempb8.y() * w8);

    tempb.Setz(tempb1.z() * w1 + tempb2.z() * w2 + tempb3.z() * w3 + tempb4.z() * w4 + tempb5.z() * w5 + tempb6.z() * w6 + tempb7.z() * w7 + tempb8.z() * w8);

    Vector3 tempe;
    tempe.Setx(tempe1.x() * w1 + tempe2.x() * w2 + tempe3.x() * w3 + tempe4.x() * w4 + tempe5.x() * w5 + tempe6.x() * w6 + tempe7.x() * w7 + tempe8.x() * w8);

    tempe.Sety(tempe1.y() * w1 + tempe2.y() * w2 + tempe3.y() * w3 + tempe4.y() * w4 + tempe5.y() * w5 + tempe6.y() * w6 + tempe7.y() * w7 + tempe8.y() * w8);

    tempe.Setz(tempe1.z() * w1 + tempe2.z() * w2 + tempe3.z() * w3 + tempe4.z() * w4 + tempe5.z() * w5 + tempe6.z() * w6 + tempe7.z() * w7 + tempe8.z() * w8);

    Vector3 tempv;
    if (mi0 == mi0_H)
    {
        tempv1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->VelH3();
        tempv2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->VelH3();
        tempv3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->VelH3();
        tempv4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->VelH3();
        tempv5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->VelH3();
        tempv6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->VelH3();
        tempv7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->VelH3();
        tempv8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->VelH3();

        tempv.Setz(tempv1.z() * w1 + tempv2.z() * w2 + tempv3.z() * w3 + tempv4.z() * w4 + tempv5.z() * w5 + tempv6.z() * w6 + tempv7.z() * w7 + tempv8.z() * w8);
    }
    else if (mi0 == mi0_He)
    {
        tempv1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->VelHe3();
        tempv2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->VelHe3();
        tempv3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->VelHe3();
        tempv4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->VelHe3();
        tempv5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->VelHe3();
        tempv6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->VelHe3();
        tempv7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->VelHe3();
        tempv8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->VelHe3();

        tempv.Setz(tempv1.z() * w1 + tempv2.z() * w2 + tempv3.z() * w3 + tempv4.z() * w4 + tempv5.z() * w5 + tempv6.z() * w6 + tempv7.z() * w7 + tempv8.z() * w8);
    }
    else if (mi0 == mi0_O)
    {
        tempv1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->VelO3();
        tempv2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->VelO3();
        tempv3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->VelO3();
        tempv4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->VelO3();
        tempv5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->VelO3();
        tempv6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->VelO3();
        tempv7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->VelO3();
        tempv8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->VelO3();

        tempv.Setz(tempv1.z() * w1 + tempv2.z() * w2 + tempv3.z() * w3 + tempv4.z() * w4 + tempv5.z() * w5 + tempv6.z() * w6 + tempv7.z() * w7 + tempv8.z() * w8);
    }
    else
    {
        std::cout << " MaxwellDisV tempV error \n";
    }
    // tempv - bulk velocity; tempV - velocity of a particle
    Vector3 tempV = Vector3(MaxwellDisV(ikT, tempv.x(), mi0),
                            MaxwellDisV(ikT, tempv.y(), mi0),
                            MaxwellDisV(ikT, tempv.z(), mi0));
    //    if( mi0==mi0_H) std::cout << tempV.x() << " " << tempV.y() << " " << tempV.z() << " \n";
    Vector3 tempDriftV = {0.0, 0.0, 0.0};
    if (tempb.norm2() > 0.0)
        tempDriftV = tempe.CrossProduct(tempb).ScaleProduct(1.0 / tempb.norm2());
    //    std::cout << " driftV " << tempDriftV.x() << "\n";
    return tempb.NormalizedVector().ScaleProduct(tempV.DotProduct(tempb.NormalizedVector())).PlusProduct(tempDriftV);
}

Vector3 MaxwellDisV_1(GridsPoints *****ptrArray_in, uint_64 posUint, double mi0, int n)
{
    //************************************************************************
    struct structg strg_in = {0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
    strg_in.face = posUint >> 61;
    for (int i = 0; i < fieldsGridsLevel; i++)
    {
        strg_in.ig = (strg_in.ig << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jg = (strg_in.jg << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kg = (strg_in.kg << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }

    for (int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in.iw = (strg_in.iw << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jw = (strg_in.jw << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kw = (strg_in.kw << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }
    strg_in.vx = 0.0;
    strg_in.vy = 0.0;
    strg_in.vz = 0.0; // not need
                      //   strg_in.mass = 0.0; // not need
                      //************************************************************************
                      //   if( strg_in.kg >= fieldsGridsSize - tempGridsCellLevel)
                      //   strg_in.kg = strg_in.kg + tempGridsCellLevel - fieldsGridsSize;

    Vector3 tempb1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->B3();
    Vector3 tempb8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->B3();

    double w1 = 1.0 - 1.0 * (strg_in.iw + 1) * (strg_in.jw + 1) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w2 = 1.0 - 1.0 * (cellSize1 - strg_in.iw) * (strg_in.jw + 1) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w3 = 1.0 - 1.0 * (cellSize1 - strg_in.iw) * (cellSize1 - strg_in.jw) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w4 = 1.0 - 1.0 * (strg_in.iw + 1) * (cellSize1 - strg_in.jw) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w5 = 1.0 - 1.0 * (strg_in.iw + 1) * (strg_in.jw + 1) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;
    double w6 = 1.0 - 1.0 * (cellSize1 - strg_in.iw) * (strg_in.jw + 1) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;
    double w7 = 1.0 - 1.0 * (cellSize1 - strg_in.iw) * (cellSize1 - strg_in.jw) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;
    double w8 = 1.0 - 1.0 * (strg_in.iw + 1) * (cellSize1 - strg_in.jw) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;

    Vector3 tempb;
    tempb.Setx(tempb1.x() * w1 + tempb2.x() * w2 + tempb3.x() * w3 + tempb4.x() * w4 + tempb5.x() * w5 + tempb6.x() * w6 + tempb7.x() * w7 + tempb8.x() * w8);

    tempb.Sety(tempb1.y() * w1 + tempb2.y() * w2 + tempb3.y() * w3 + tempb4.y() * w4 + tempb5.y() * w5 + tempb6.y() * w6 + tempb7.y() * w7 + tempb8.y() * w8);

    tempb.Setz(tempb1.z() * w1 + tempb2.z() * w2 + tempb3.z() * w3 + tempb4.z() * w4 + tempb5.z() * w5 + tempb6.z() * w6 + tempb7.z() * w7 + tempb8.z() * w8);

    Vector3 tempV = Vector3(MaxwellDisV(ikT, ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Vel_e3().x(), mi0),
                            MaxwellDisV(ikT, ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Vel_e3().y(), mi0),
                            MaxwellDisV(ikT, ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Vel_e3().z(), mi0));

    return tempb.NormalizedVector().ScaleProduct(tempV.DotProduct(tempb.NormalizedVector()));
}

//************************************************************************
//************************************************************************
// FUNCTION return maxwell distrubtion random double number of particle
// energy for magnetic moment/ adiabatic invarient 10~1000 eV
// the range for random number is 1 ~ 3 and then negetive it to make the
// major particles is between 10~1000 eV. should make sigma_in < 0.2
double MaxwellDisEnergy(GridsPoints *****ptrArray_in, uint_64 posUint)
{
    double temp = 0.0;
    double sigma_in = sigma_MaxwellDis;
    double mu_in = mu_MaxwellDis;
    while (temp < 0.15 || temp > 0.85)
    {
        temp = dRand();
    }
    double x_temp = energy1eV * 0.6667 * pow(10.0, erfinv(2.0 * temp - 1.0) * sqrt(2.0) * sigma_in + mu_in);

    //************************************************************************
    struct structg strg_in = {0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
    strg_in.face = posUint >> 61;
    for (int i = 0; i < fieldsGridsLevel; i++)
    {
        strg_in.ig = (strg_in.ig << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jg = (strg_in.jg << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kg = (strg_in.kg << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }

    for (int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
    {
        strg_in.iw = (strg_in.iw << 1) + ((posUint >> (60 - i * 3)) & 1);
        strg_in.jw = (strg_in.jw << 1) + ((posUint >> (60 - 1 - i * 3)) & 1);
        strg_in.kw = (strg_in.kw << 1) + ((posUint >> (60 - 2 - i * 3)) & 1);
    }
    strg_in.vx = 0.0;
    strg_in.vy = 0.0;
    strg_in.vz = 0.0; // not need
                      //   strg_in.mass = 0.0; // not need
                      //************************************************************************
    if (strg_in.kg >= fieldsGridsSize - tempGridsCellLevelBot)
        strg_in.kg = strg_in.kg + tempGridsCellLevelBot - fieldsGridsSize;

    Vector3 tempb1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 tempb3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 tempb5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 tempb7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->B3();
    Vector3 tempb8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->B3();

    double w1 = 1 - (strg_in.iw + 1) * (strg_in.jw + 1) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w2 = 1 - (cellSize1 - strg_in.iw) * (strg_in.jw + 1) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w3 = 1 - (cellSize1 - strg_in.iw) * (cellSize1 - strg_in.jw) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w4 = 1 - (strg_in.iw + 1) * (cellSize1 - strg_in.jw) * (strg_in.kw + 1) / cellSize1 / cellSize1 / cellSize1;
    double w5 = 1 - (strg_in.iw + 1) * (strg_in.jw + 1) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;
    double w6 = 1 - (cellSize1 - strg_in.iw) * (strg_in.jw + 1) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;
    double w7 = 1 - (cellSize1 - strg_in.iw) * (cellSize1 - strg_in.jw) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;
    double w8 = 1 - (strg_in.iw + 1) * (cellSize1 - strg_in.jw) * (cellSize1 - strg_in.kw) / cellSize1 / cellSize1 / cellSize1;

    Vector3 tempb;
    tempb.Setx(tempb1.x() * w1 + tempb2.x() * w2 + tempb3.x() * w3 + tempb4.x() * w4 + tempb5.x() * w5 + tempb6.x() * w6 + tempb7.x() * w7 + tempb8.x() * w8);

    tempb.Sety(tempb1.y() * w1 + tempb2.y() * w2 + tempb3.y() * w3 + tempb4.y() * w4 + tempb5.y() * w5 + tempb6.y() * w6 + tempb7.y() * w7 + tempb8.y() * w8);

    tempb.Setz(tempb1.z() * w1 + tempb2.z() * w2 + tempb3.z() * w3 + tempb4.z() * w4 + tempb5.z() * w5 + tempb6.z() * w6 + tempb7.z() * w7 + tempb8.z() * w8);

    return x_temp / tempb.norm();
}
