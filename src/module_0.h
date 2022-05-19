#ifndef _MODULE_H_0_
#define _MODULE_H_0_
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "parameters.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include "module_0.h"
#include <cmath>
#include <limits>
#include <bitset>
#include <vector>
#include "module_base.h"
#include "gridscells.h"
#include <fstream>
#include <string>
#include "geopack.h"

using std::vector;

//************************************************************************
//************************************************************************
// FUNCTION // Return a authogonal vector to a face when provied a location
// (face, i, j, k). The norm is 1.
//  We want to point out all the six face orthogonal unit vector.

inline Vector3 UnitVectorL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().OrthoUnitVector(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());
}

inline Vector3 UnitVectorT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3().OrthoUnitVector(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
}

inline Vector3 UnitVectorR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3().OrthoUnitVector(ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3());
}

inline Vector3 UnitVectorBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3().OrthoUnitVector(ptrArray_in[face_in][i_in][j_in][k_in]->Pos3());
}

inline Vector3 UnitVectorF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
    return temp1.PlusProduct(temp2).NormalizedVector();
}

inline Vector3 UnitVectorBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    return UnitVectorF(ptrArray_in, face_in, i_in, j_in, k_in).ScaleProduct(-1.0);
}

//************************************************************************
//************************************************************************
// FUNCTION // Return a average field vector on a face when provied a
// location (face, i, j, k). The norm is of the average of four vector at the
// nearest gridspoints.
// We want to point out all the six face field vector for E and B and P

inline Vector3 FaceFieldVectorL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch (field_in)
    {
    case 'E':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->E3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->E3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->E3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->E3());
        break;
    }
    case 'B':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->B3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->B3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->B3());
        break;
    }
    case 'D':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->DB3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->DB3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->DB3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->DB3());
        break;
    }
    case 'P':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3());
        break;
    }
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch (field_in)
    {
    case 'E':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->E3());
        temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->E3());
        break;
    }
    case 'B':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3());
        temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->B3());
        break;
    }
    case 'P':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
        temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3());
        break;
    }
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch (field_in)
    {
    case 'E':
    {
        temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->E3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->E3());
        break;
    }
    case 'B':
    {
        temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->B3());
        break;
    }
    case 'P':
    {
        temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3());
        break;
    }
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch (field_in)
    {
    case 'E':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in]->E3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->E3());
        break;
    }
    case 'B':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->B3());
        break;
    }
    case 'P':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
        temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3());
        break;
    }
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25); // 1/4 is not acceptable !
};

inline Vector3 FaceFieldVectorF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch (field_in)
    {
    case 'E':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->E3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->E3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->E3());
        break;
    }
    case 'B':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->B3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->B3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->B3());
        break;
    }
    case 'P':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3());
        break;
    }
    }
    //    std::cout << " PPP" << temp1.PlusProduct(temp2).ScaleProduct(0.25).x() << " " << temp1.PlusProduct(temp2).y() << " " << temp1.PlusProduct(temp2).z() << " " ;

    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

inline Vector3 FaceFieldVectorBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, char field_in)
{
    Vector3 temp1, temp2;
    switch (field_in)
    {
    case 'E':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->E3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->E3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->E3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->E3());
        break;
    }
    case 'B':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->B3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->B3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3());
        break;
    }
    case 'P':
    {
        temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());
        temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3().PlusProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
        break;
    }
    }
    return temp1.PlusProduct(temp2).ScaleProduct(0.25);
};

//************************************************************************
//************************************************************************
// FUNCTION
// Return a average density of ions( electrons) on a face when provied a
// location (face, i, j, k). The norm is of the average of four densities
// at the nearest gridspoints.
// We want to point out all the six face field vector for E and B and P

inline double FaceDensityL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density();

    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceDensityT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceDensityR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceDensityBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceDensityF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceDensityBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density();
    return (temp1 + temp2) * (1.0 / 4.0);
};

//************************************************************************
//************************************************************************
// FUNCTION
// Return a average number density of ions( electrons) on a face when provied a
// location (face, i, j, k). The norm is of the average of four densities
// at the nearest gridspoints.

inline double FaceNumberDensityL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_O();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_H() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_He() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_O() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_O();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNumberDensityT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_O();
    temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_H() +
            ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_O();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNumberDensityR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_O();

    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_H() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_He() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_O();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNumberDensityBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_O();

    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_O();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNumberDensityF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_H() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_He() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in + 1]->Density_O() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Density_O();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_H() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_He() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Density_O();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNumberDensityBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Density_O();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_H() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_H() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_He() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_He() +
            ptrArray_in[face_in][i_in + 1][j_in][k_in]->Density_O() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Density_O();
    return (temp1 + temp2) * (1.0 / 4.0);
};

//************************************************************************
//************************************************************************
// FUNCTION
// Return a average Temperature of electrons on a face when provied a
// location (face, i, j, k). The norm is of the average of four temperature
// at the nearest gridspoints.
// We want to point out all the six face field vector for E and B and P

inline double FaceTemperatureL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Temperature() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Temperature();

    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceTemperatureT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Temperature();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceTemperatureR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Temperature();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceTemperatureBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Temperature();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceTemperatureF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Temperature() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Temperature();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Temperature();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceTemperatureBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Temperature();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Temperature() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Temperature();
    return (temp1 + temp2) * (1.0 / 4.0);
};

//************************************************************************
//************************************************************************
// FUNCTION
// Return a average>B3().norm(of electrons on a face when provied a
// location (face, i, j, k). The norm is of the average of four temperature
// at the nearest gridspoints.
// We want to point out all the six face field vector for E and B and P

inline double FaceNormBL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->B3().norm() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->B3().norm();

    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNormBT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->B3().norm();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNormBR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->B3().norm();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNormBBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->B3().norm();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNormBF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->B3().norm() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->B3().norm();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FaceNormBBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->B3().norm();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->B3().norm() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->B3().norm();
    return (temp1 + temp2) * (1.0 / 4.0);
};

//************************************************************************
//************************************************************************
// FUNCTION
// Return a average Potential() on a face when provied a
// location (face, i, j, k). The norm is of the average of four potential
// at the nearest gridspoints.

inline double FacePotentialL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Potential() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Potential();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Potential() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Potential();

    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FacePotentialT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->Potential() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Potential();
    temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Potential() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Potential();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FacePotentialR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Potential() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Potential();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Potential() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Potential();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FacePotentialBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Potential() + ptrArray_in[face_in][i_in + 1][j_in][k_in]->Potential();
    temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Potential() + ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Potential();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FacePotentialF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Potential() + ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Potential();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Potential() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Potential();
    return (temp1 + temp2) * (1.0 / 4.0);
};

inline double FacePotentialBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double temp1, temp2;
    temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Potential() + ptrArray_in[face_in][i_in][j_in + 1][k_in]->Potential();
    temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Potential() + ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Potential();
    return (temp1 + temp2) * (1.0 / 4.0);
};

//************************************************************************
//************************************************************************
// FUNCTION // Return a vector athogonal to a face when provied a
// location (face, i, j, k).
// We want to point out all the six vectors whose norm is the face area.
// For the side area, use bigger tri-angle minus smaller tri-angle to represent
// the norm.
// For the front and back area, use two small tri-angle sum to represent the
// norm, and use average vector, as UnitVectorF and UnitVectorBack do, to
// represent the direction.

inline Vector3 AreaVectorL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().CrossProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3().CrossProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.50);
}

inline Vector3 AreaVectorT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3().CrossProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3().CrossProduct(ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.50);
};

inline Vector3 AreaVectorR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3().CrossProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3().CrossProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.50);
};

inline Vector3 AreaVectorBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3().CrossProduct(ptrArray_in[face_in][i_in][j_in][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3().CrossProduct(ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3());
    return temp2.MinusProduct(temp1).ScaleProduct(0.50);
};

inline Vector3 AreaVectorF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3().MinusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3().MinusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->Pos3());

    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->Pos3().MinusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->Pos3().MinusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->Pos3());
    /*    std::cout << " check " << ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in][j_in+1][k_in]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in+1][j_in][k_in]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in+1][j_in+1][k_in]->Pos3().x() << " - " <<
                      ptrArray_in[face_in][i_in][j_in][k_in+1]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in][j_in+1][k_in+1]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in+1][j_in][k_in+1]->Pos3().x() << 
                      " - " << ptrArray_in[face_in][i_in+1][j_in+1][k_in+1]->Pos3().x() << " - " << 
                      temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5).x() << " " 
                      <<temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5).y() << " " 
                      << temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5).z();
                      int pause ;
                      std::cin >> pause;
                      */
    return temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5);
};

inline Vector3 AreaVectorBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().MinusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3().MinusProduct(ptrArray_in[face_in][i_in + 1][j_in][k_in]->Pos3());

    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().MinusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->Pos3().MinusProduct(ptrArray_in[face_in][i_in][j_in + 1][k_in]->Pos3());

    return temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(-0.5);
};

////////////////////////////////////// dual ///////////////////////////
inline Vector3 AreaVectorL(Vector3 *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in][j_in][k_in]);
    Vector3 temp2 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in][j_in][k_in]);

    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]);
    Vector3 temp4 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]);

    return temp1.CrossProduct(temp2).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5);
}

inline Vector3 AreaVectorT(Vector3 *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in]);
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in]);

    Vector3 temp3 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]);
    Vector3 temp4 = ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]);

    return temp1.CrossProduct(temp2).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5);
};

inline Vector3 AreaVectorR(Vector3 *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in]);
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in]);

    Vector3 temp3 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]);
    Vector3 temp4 = ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]);

    return temp1.CrossProduct(temp2).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(-0.5);
};

inline Vector3 AreaVectorBot(Vector3 *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in][j_in][k_in]);
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->MinusProduct(*ptrArray_in[face_in][i_in][j_in][k_in]);

    Vector3 temp3 = ptrArray_in[face_in][i_in + 1][j_in][k_in]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]);
    Vector3 temp4 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]);

    return temp1.CrossProduct(temp2).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(-0.5);
};

inline Vector3 AreaVectorF(Vector3 *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]);
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in + 1]);

    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]);
    Vector3 temp4 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in + 1]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in + 1]);

    return temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(0.5);
};

inline Vector3 AreaVectorBack(Vector3 *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 temp1 = ptrArray_in[face_in][i_in][j_in][k_in]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in]);
    Vector3 temp2 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in + 1][j_in][k_in]);

    Vector3 temp3 = ptrArray_in[face_in][i_in][j_in][k_in]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in]);
    Vector3 temp4 = ptrArray_in[face_in][i_in + 1][j_in + 1][k_in]->MinusProduct(*ptrArray_in[face_in][i_in][j_in + 1][k_in]);

    return temp2.CrossProduct(temp1).PlusProduct(temp3.CrossProduct(temp4)).ScaleProduct(-0.5);
};

//************************************************************************
//************************************************************************
// FUNCTION, update E
// i, j, k are the index of main cell, range are [0-fsize+1][0-fsize+1][0-fsize]
// while the value should be in range [1-fsize+1][1-fsize+1][?]
inline Vector3 EOnLength(GridsPoints *****ptrArray, int f, int i, int j, int k, int dir)
{
    Vector3 tempE;
    if (dir == 0)
    {
        // two points at ends of length are: (f, i, j, k) and (f, i+1, j, k)
        // four related main cells are: (f, i, j, k); (f, i, j-1, k); (f, i, j-1, k-1); (f, i, j, k-1)
        // weight of gridspoints: 4 - (f,i,j,k) (f,i+1,j.k)
        //                      : 2 - (f,i,j+1,k) (f,i+1,j+1,k)
        //                      : 2 - (f,i,j-1,k) (f,i+1,j-1,k)
        //                      : 2 - (f,i,j,k+1) (f,i+1,j,k+1)
        //                      : 2 - (f,i,j,k-1) (f,i+1,j,k-1)
        // weight of girdspoints: 1 - (f,i,j+1,k+1) (f,i+1,j+1,k+1)
        //                      : 1 - (f,i,j+1,k-1) (f,i+1,j+1,k-1)
        //                      : 1 - (f,i,j-1,k+1) (f,i+1,j-1,k+1)
        //                      : 1 - (f,i,j-1,k-1) (f,i+1,j-1,k-1)
        // total 32.
        if (j > 0 && j < fieldsGridsSize + 2 && i > 0 && i < fieldsGridsSize + 1)
        {
            tempE =
                Vector3(4.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i + 1][j][k]->E3().x()) +
                            2.0 * (ptrArray[f][i][j + 1][k]->E3().x() + ptrArray[f][i + 1][j + 1][k]->E3().x() +
                                   ptrArray[f][i][j - 1][k]->E3().x() + ptrArray[f][i + 1][j - 1][k]->E3().x() +
                                   ptrArray[f][i][j][k + 1]->E3().x() + ptrArray[f][i + 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i][j][k - 1]->E3().x() + ptrArray[f][i + 1][j][k - 1]->E3().x()) +
                            1.0 * (ptrArray[f][i][j + 1][k + 1]->E3().x() + ptrArray[f][i + 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j + 1][k - 1]->E3().x() + ptrArray[f][i + 1][j + 1][k - 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k + 1]->E3().x() + ptrArray[f][i + 1][j - 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k - 1]->E3().x() + ptrArray[f][i + 1][j - 1][k - 1]->E3().x()),
                        4.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i + 1][j][k]->E3().y()) +
                            2.0 * (ptrArray[f][i][j + 1][k]->E3().y() + ptrArray[f][i + 1][j + 1][k]->E3().y() +
                                   ptrArray[f][i][j - 1][k]->E3().y() + ptrArray[f][i + 1][j - 1][k]->E3().y() +
                                   ptrArray[f][i][j][k + 1]->E3().y() + ptrArray[f][i + 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i][j][k - 1]->E3().y() + ptrArray[f][i + 1][j][k - 1]->E3().y()) +
                            1.0 * (ptrArray[f][i][j + 1][k + 1]->E3().y() + ptrArray[f][i + 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j + 1][k - 1]->E3().y() + ptrArray[f][i + 1][j + 1][k - 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k + 1]->E3().y() + ptrArray[f][i + 1][j - 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k - 1]->E3().y() + ptrArray[f][i + 1][j - 1][k - 1]->E3().y()),
                        4.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i + 1][j][k]->E3().z()) +
                            2.0 * (ptrArray[f][i][j + 1][k]->E3().z() + ptrArray[f][i + 1][j + 1][k]->E3().z() +
                                   ptrArray[f][i][j - 1][k]->E3().z() + ptrArray[f][i + 1][j - 1][k]->E3().z() +
                                   ptrArray[f][i][j][k + 1]->E3().z() + ptrArray[f][i + 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i][j][k - 1]->E3().z() + ptrArray[f][i + 1][j][k - 1]->E3().z()) +
                            1.0 * (ptrArray[f][i][j + 1][k + 1]->E3().z() + ptrArray[f][i + 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j + 1][k - 1]->E3().z() + ptrArray[f][i + 1][j + 1][k - 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k + 1]->E3().z() + ptrArray[f][i + 1][j - 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k - 1]->E3().z() + ptrArray[f][i + 1][j - 1][k - 1]->E3().z()))
                    .ScaleProduct(1.0 / 32.0);
        }
        else
        {
            std::cout << " direction 0 range out E on length\n";
        }
    }
    else if (dir == 1)
    {
        // two points at ends of length are: (f, i, j, k) and (f, i, j+1, k)
        // four releted main cells are: (f, i, j, k); (f, i-1, j, k); (f, i, j, k-1); (f, i-1, j, k-1)

        // weight of gridspoints: 4 - (f,i,j,k) (f,i,j+1,k)
        //                      : 2 - (f,i+1,j,k) (f,i+1,j+1,k)
        //                      : 2 - (f,i-1,j,k) (f,i-1,j+1,k)
        //                      : 2 - (f,i,j,k+1) (f,i,j+1,k+1)
        //                      : 2 - (f,i,j,k-1) (f,i,j+1,k-1)
        // weight of girdspoints: 1 - (f,i+1,j,k+1) (f,i+1,j+1,k+1)
        //                      : 1 - (f,i+1,j,k-1) (f,i+1,j+1,k-1)
        //                      : 1 - (f,i-1,j,k+1) (f,i-1,j+1,k+1)
        //                      : 1 - (f,i-1,j,k-1) (f,i-1,j+1,k-1)
        //
        if (i > 0 && i < fieldsGridsSize + 2 && j > 0 && j < fieldsGridsSize + 1)
        {
            tempE =
                Vector3(4.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i][j + 1][k]->E3().x()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().x() + ptrArray[f][i + 1][j + 1][k]->E3().x() +
                                   ptrArray[f][i - 1][j][k]->E3().x() + ptrArray[f][i - 1][j + 1][k]->E3().x() +
                                   ptrArray[f][i][j][k + 1]->E3().x() + ptrArray[f][i][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j][k - 1]->E3().x() + ptrArray[f][i][j + 1][k - 1]->E3().x()) +
                            1.0 * (ptrArray[f][i + 1][j][k + 1]->E3().x() + ptrArray[f][i + 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i + 1][j][k - 1]->E3().x() + ptrArray[f][i + 1][j + 1][k - 1]->E3().x() +
                                   ptrArray[f][i - 1][j][k + 1]->E3().x() + ptrArray[f][i - 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j][k - 1]->E3().x() + ptrArray[f][i - 1][j + 1][k - 1]->E3().x()),
                        4.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i][j + 1][k]->E3().y()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().y() + ptrArray[f][i + 1][j + 1][k]->E3().y() +
                                   ptrArray[f][i - 1][j][k]->E3().y() + ptrArray[f][i - 1][j + 1][k]->E3().y() +
                                   ptrArray[f][i][j][k + 1]->E3().y() + ptrArray[f][i][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j][k - 1]->E3().y() + ptrArray[f][i][j + 1][k - 1]->E3().y()) +
                            1.0 * (ptrArray[f][i + 1][j][k + 1]->E3().y() + ptrArray[f][i + 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i + 1][j][k - 1]->E3().y() + ptrArray[f][i + 1][j + 1][k - 1]->E3().y() +
                                   ptrArray[f][i - 1][j][k + 1]->E3().y() + ptrArray[f][i - 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j][k - 1]->E3().y() + ptrArray[f][i - 1][j + 1][k - 1]->E3().y()),
                        4.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i][j + 1][k]->E3().z()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().z() + ptrArray[f][i + 1][j + 1][k]->E3().z() +
                                   ptrArray[f][i - 1][j][k]->E3().z() + ptrArray[f][i - 1][j + 1][k]->E3().z() +
                                   ptrArray[f][i][j][k + 1]->E3().z() + ptrArray[f][i][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j][k - 1]->E3().z() + ptrArray[f][i][j + 1][k - 1]->E3().z()) +
                            1.0 * (ptrArray[f][i + 1][j][k + 1]->E3().z() + ptrArray[f][i + 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i + 1][j][k - 1]->E3().z() + ptrArray[f][i + 1][j + 1][k - 1]->E3().z() +
                                   ptrArray[f][i - 1][j][k + 1]->E3().z() + ptrArray[f][i - 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j][k - 1]->E3().z() + ptrArray[f][i - 1][j + 1][k - 1]->E3().z()))
                    .ScaleProduct(1.0 / 32.0);
        }
        else
        {
            std::cout << " direction 1 range out E on length\n";
        }
    }
    else if (dir == 2)
    {
        // two points at ends of length are: (f, i, j, k) and (f, i, j, k+1)

        // four releted main cells are: (f, i, j, k); (f, i-1, j, k); (f, i, j-1, k); (f, i-1, j-1, k)
        // weight of gridspoints: 4 - (f,i,j,k) (f,i,j,k+1)
        //                      : 2 - (f,i+1,j,k) (f,i+1,j,k+1)
        //                      : 2 - (f,i-1,j,k) (f,i-1,j,k+1)
        //                      : 2 - (f,i,j+1,k) (f,i,j+1,k+1)
        //                      : 2 - (f,i,j-1,k) (f,i,j-1,k+1)
        //
        //                      : 1 - (f,i+1,j+1,k) (f,i+1,j+1,k+1)
        //                      : 1 - (f,i+1,j-1,k) (f,i+1,j-1,k+1)
        //                      : 1 - (f,i-1,j+1,k) (f,i-1,j+1,k+1)
        //                      : 1 - (f,i-1,j-1,k) (f,i-1,j-1,k+1)
        if (i == 1 && j == fieldsGridsSize + 1)
        {
            tempE =
                Vector3(3.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i][j][k + 1]->E3().x()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().x() + ptrArray[f][i + 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i][j + 1][k]->E3().x() + ptrArray[f][i][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k]->E3().x() + ptrArray[f][i][j - 1][k + 1]->E3().x()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().x() + ptrArray[f][i + 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().x() + ptrArray[f][i + 1][j - 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().x() + ptrArray[f][i - 1][j - 1][k + 1]->E3().x()),
                        3.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i][j][k + 1]->E3().y()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().y() + ptrArray[f][i + 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i][j + 1][k]->E3().y() + ptrArray[f][i][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k]->E3().y() + ptrArray[f][i][j - 1][k + 1]->E3().y()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().y() + ptrArray[f][i + 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().y() + ptrArray[f][i + 1][j - 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().y() + ptrArray[f][i - 1][j - 1][k + 1]->E3().y()),
                        3.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i][j][k + 1]->E3().z()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().z() + ptrArray[f][i + 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i][j + 1][k]->E3().z() + ptrArray[f][i][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k]->E3().z() + ptrArray[f][i][j - 1][k + 1]->E3().z()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().z() + ptrArray[f][i + 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().z() + ptrArray[f][i + 1][j - 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().z() + ptrArray[f][i - 1][j - 1][k + 1]->E3().z()))
                    .ScaleProduct(1.0 / 24.0);
        }
        else if (i == 1 && j == 1)
        {
            tempE =
                Vector3(3.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i][j][k + 1]->E3().x()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().x() + ptrArray[f][i + 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i][j + 1][k]->E3().x() + ptrArray[f][i][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k]->E3().x() + ptrArray[f][i][j - 1][k + 1]->E3().x()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().x() + ptrArray[f][i + 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().x() + ptrArray[f][i + 1][j - 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().x() + ptrArray[f][i - 1][j + 1][k + 1]->E3().x()),
                        3.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i][j][k + 1]->E3().y()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().y() + ptrArray[f][i + 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i][j + 1][k]->E3().y() + ptrArray[f][i][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k]->E3().y() + ptrArray[f][i][j - 1][k + 1]->E3().y()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().y() + ptrArray[f][i + 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().y() + ptrArray[f][i + 1][j - 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().y() + ptrArray[f][i - 1][j + 1][k + 1]->E3().y()),
                        3.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i][j][k + 1]->E3().z()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().z() + ptrArray[f][i + 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i][j + 1][k]->E3().z() + ptrArray[f][i][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k]->E3().z() + ptrArray[f][i][j - 1][k + 1]->E3().z()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().z() + ptrArray[f][i + 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().z() + ptrArray[f][i + 1][j - 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().z() + ptrArray[f][i - 1][j + 1][k + 1]->E3().z()))
                    .ScaleProduct(1.0 / 24.0);
        }
        else if (i == fieldsGridsSize + 1 && j == 1)
        {
            tempE =
                Vector3(3.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i][j][k + 1]->E3().x()) +
                            2.0 * (ptrArray[f][i - 1][j][k]->E3().x() + ptrArray[f][i - 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i][j + 1][k]->E3().x() + ptrArray[f][i][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k]->E3().x() + ptrArray[f][i][j - 1][k + 1]->E3().x()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().x() + ptrArray[f][i + 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().x() + ptrArray[f][i - 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().x() + ptrArray[f][i - 1][j - 1][k + 1]->E3().x()),
                        3.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i][j][k + 1]->E3().y()) +
                            2.0 * (ptrArray[f][i - 1][j][k]->E3().y() + ptrArray[f][i - 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i][j + 1][k]->E3().y() + ptrArray[f][i][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k]->E3().y() + ptrArray[f][i][j - 1][k + 1]->E3().y()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().y() + ptrArray[f][i + 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().y() + ptrArray[f][i - 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().y() + ptrArray[f][i - 1][j - 1][k + 1]->E3().y()),
                        3.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i][j][k + 1]->E3().z()) +
                            2.0 * (ptrArray[f][i - 1][j][k]->E3().z() + ptrArray[f][i - 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i][j + 1][k]->E3().z() + ptrArray[f][i][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k]->E3().z() + ptrArray[f][i][j - 1][k + 1]->E3().z()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().z() + ptrArray[f][i + 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().z() + ptrArray[f][i - 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().z() + ptrArray[f][i - 1][j - 1][k + 1]->E3().z()))
                    .ScaleProduct(1.0 / 24.0);
        }
        else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
        {
            tempE =
                Vector3(3.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i][j][k + 1]->E3().x()) +
                            2.0 * (ptrArray[f][i - 1][j][k]->E3().x() + ptrArray[f][i - 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i][j + 1][k]->E3().x() + ptrArray[f][i][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k]->E3().x() + ptrArray[f][i][j - 1][k + 1]->E3().x()) +
                            1.0 * (ptrArray[f][i + 1][j - 1][k]->E3().x() + ptrArray[f][i + 1][j - 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().x() + ptrArray[f][i - 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().x() + ptrArray[f][i - 1][j - 1][k + 1]->E3().x()),
                        3.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i][j][k + 1]->E3().y()) +
                            2.0 * (ptrArray[f][i - 1][j][k]->E3().y() + ptrArray[f][i - 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i][j + 1][k]->E3().y() + ptrArray[f][i][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k]->E3().y() + ptrArray[f][i][j - 1][k + 1]->E3().y()) +
                            1.0 * (ptrArray[f][i + 1][j - 1][k]->E3().y() + ptrArray[f][i + 1][j - 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().y() + ptrArray[f][i - 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().y() + ptrArray[f][i - 1][j - 1][k + 1]->E3().y()),
                        3.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i][j][k + 1]->E3().z()) +
                            2.0 * (ptrArray[f][i - 1][j][k]->E3().z() + ptrArray[f][i - 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i][j + 1][k]->E3().z() + ptrArray[f][i][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k]->E3().z() + ptrArray[f][i][j - 1][k + 1]->E3().z()) +
                            1.0 * (ptrArray[f][i + 1][j - 1][k]->E3().z() + ptrArray[f][i + 1][j - 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().z() + ptrArray[f][i - 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().z() + ptrArray[f][i - 1][j - 1][k + 1]->E3().z()))
                    .ScaleProduct(1.0 / 24.0);
        }
        else
        {
            tempE =
                Vector3(4.0 * (ptrArray[f][i][j][k]->E3().x() + ptrArray[f][i][j][k + 1]->E3().x()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().x() + ptrArray[f][i + 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j][k]->E3().x() + ptrArray[f][i - 1][j][k + 1]->E3().x() +
                                   ptrArray[f][i][j + 1][k]->E3().x() + ptrArray[f][i][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i][j - 1][k]->E3().x() + ptrArray[f][i][j - 1][k + 1]->E3().x()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().x() + ptrArray[f][i + 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().x() + ptrArray[f][i + 1][j - 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().x() + ptrArray[f][i - 1][j + 1][k + 1]->E3().x() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().x() + ptrArray[f][i - 1][j - 1][k + 1]->E3().x()),
                        4.0 * (ptrArray[f][i][j][k]->E3().y() + ptrArray[f][i][j][k + 1]->E3().y()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().y() + ptrArray[f][i + 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j][k]->E3().y() + ptrArray[f][i - 1][j][k + 1]->E3().y() +
                                   ptrArray[f][i][j + 1][k]->E3().y() + ptrArray[f][i][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i][j - 1][k]->E3().y() + ptrArray[f][i][j - 1][k + 1]->E3().y()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().y() + ptrArray[f][i + 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().y() + ptrArray[f][i + 1][j - 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().y() + ptrArray[f][i - 1][j + 1][k + 1]->E3().y() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().y() + ptrArray[f][i - 1][j - 1][k + 1]->E3().y()),
                        4.0 * (ptrArray[f][i][j][k]->E3().z() + ptrArray[f][i][j][k + 1]->E3().z()) +
                            2.0 * (ptrArray[f][i + 1][j][k]->E3().z() + ptrArray[f][i + 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j][k]->E3().z() + ptrArray[f][i - 1][j][k + 1]->E3().z() +
                                   ptrArray[f][i][j + 1][k]->E3().z() + ptrArray[f][i][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i][j - 1][k]->E3().z() + ptrArray[f][i][j - 1][k + 1]->E3().z()) +
                            1.0 * (ptrArray[f][i + 1][j + 1][k]->E3().z() + ptrArray[f][i + 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i + 1][j - 1][k]->E3().z() + ptrArray[f][i + 1][j - 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j + 1][k]->E3().z() + ptrArray[f][i - 1][j + 1][k + 1]->E3().z() +
                                   ptrArray[f][i - 1][j - 1][k]->E3().z() + ptrArray[f][i - 1][j - 1][k + 1]->E3().z()))
                    .ScaleProduct(1.0 / 32.0);
        }
    }
    return tempE;
}

// Return a integration of circuit E towards outside
inline double EIntegrationL(GridsPoints *****ptrArray, Vector3 ******ptrEVectorFaceArray_dual, int f, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray[f][I][J][K + 1]->Pos3().MinusProduct(ptrArray[f][I][J][K]->Pos3());
    Vector3 temp2 = ptrArray[f][I][J + 1][K + 1]->Pos3().MinusProduct(ptrArray[f][I][J][K + 1]->Pos3());
    Vector3 temp3 = ptrArray[f][I][J + 1][K]->Pos3().MinusProduct(ptrArray[f][I][J + 1][K + 1]->Pos3());
    Vector3 temp4 = ptrArray[f][I][J][K]->Pos3().MinusProduct(ptrArray[f][I][J + 1][K]->Pos3());

    double EL = ptrEVectorFaceArray_dual[f][I][J][K][2]->DotProduct(temp1) +
                ptrEVectorFaceArray_dual[f][I][J][K + 1][1]->DotProduct(temp2) +
                ptrEVectorFaceArray_dual[f][I][J + 1][K][2]->DotProduct(temp3) +
                ptrEVectorFaceArray_dual[f][I][J][K][1]->DotProduct(temp4);
    return EL;
}

inline double EIntegrationR(GridsPoints *****ptrArray, Vector3 ******ptrEVectorFaceArray_dual, int f, int i, int j, int k)
{
    int I = i + 1;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray[f][I][J][K + 1]->Pos3().MinusProduct(ptrArray[f][I][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray[f][I][J + 1][K + 1]->Pos3().MinusProduct(ptrArray[f][I][J][K + 1]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray[f][I][J + 1][K]->Pos3().MinusProduct(ptrArray[f][I][J + 1][K + 1]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray[f][I][J][K]->Pos3().MinusProduct(ptrArray[f][I][J + 1][K]->Pos3()).ScaleProduct(-1.0);

    double ER = ptrEVectorFaceArray_dual[f][I][J][K][2]->DotProduct(temp1) +
                ptrEVectorFaceArray_dual[f][I][J][K + 1][1]->DotProduct(temp2) +
                ptrEVectorFaceArray_dual[f][I][J + 1][K][2]->DotProduct(temp3) +
                ptrEVectorFaceArray_dual[f][I][J][K][1]->DotProduct(temp4);
    return ER;
}

inline double EIntegrationBack(GridsPoints *****ptrArray, Vector3 ******ptrEVectorFaceArray_dual, int f, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray[f][I][J + 1][K]->Pos3().MinusProduct(ptrArray[f][I][J][K]->Pos3());
    Vector3 temp2 = ptrArray[f][I + 1][J + 1][K]->Pos3().MinusProduct(ptrArray[f][I][J + 1][K]->Pos3());
    Vector3 temp3 = ptrArray[f][I + 1][J][K]->Pos3().MinusProduct(ptrArray[f][I + 1][J + 1][K]->Pos3());
    Vector3 temp4 = ptrArray[f][I][J][K]->Pos3().MinusProduct(ptrArray[f][I + 1][J][K]->Pos3());

    double EBack = ptrEVectorFaceArray_dual[f][I][J][K][1]->DotProduct(temp1) +
                   ptrEVectorFaceArray_dual[f][I][J + 1][K][0]->DotProduct(temp2) +
                   ptrEVectorFaceArray_dual[f][I + 1][J][K][1]->DotProduct(temp3) +
                   ptrEVectorFaceArray_dual[f][I][J][K][0]->DotProduct(temp4);
    return EBack;
}

inline double EIntegrationF(GridsPoints *****ptrArray, Vector3 ******ptrEVectorFaceArray_dual, int f, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k + 1;
    Vector3 temp1 = ptrArray[f][I][J + 1][K]->Pos3().MinusProduct(ptrArray[f][I][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray[f][I + 1][J + 1][K]->Pos3().MinusProduct(ptrArray[f][I][J + 1][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray[f][I + 1][J][K]->Pos3().MinusProduct(ptrArray[f][I + 1][J + 1][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray[f][I][J][K]->Pos3().MinusProduct(ptrArray[f][I + 1][J][K]->Pos3()).ScaleProduct(-1.0);

    double EF = ptrEVectorFaceArray_dual[f][I][J][K][1]->DotProduct(temp1) +
                ptrEVectorFaceArray_dual[f][I][J + 1][K][0]->DotProduct(temp2) +
                ptrEVectorFaceArray_dual[f][I + 1][J][K][1]->DotProduct(temp3) +
                ptrEVectorFaceArray_dual[f][I][J][K][0]->DotProduct(temp4);
    return EF;
}

inline double EIntegrationBot(GridsPoints *****ptrArray, Vector3 ******ptrEVectorFaceArray_dual, int f, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray[f][I + 1][J][K]->Pos3().MinusProduct(ptrArray[f][I][J][K]->Pos3());
    Vector3 temp2 = ptrArray[f][I + 1][J][K + 1]->Pos3().MinusProduct(ptrArray[f][I + 1][J][K]->Pos3());
    Vector3 temp3 = ptrArray[f][I][J][K + 1]->Pos3().MinusProduct(ptrArray[f][I + 1][J][K + 1]->Pos3());
    Vector3 temp4 = ptrArray[f][I][J][K]->Pos3().MinusProduct(ptrArray[f][I][J][K + 1]->Pos3());

    double EBot = ptrEVectorFaceArray_dual[f][I][J][K][0]->DotProduct(temp1) +
                  ptrEVectorFaceArray_dual[f][I + 1][J][K][2]->DotProduct(temp2) +
                  ptrEVectorFaceArray_dual[f][I][J][K + 1][0]->DotProduct(temp3) +
                  ptrEVectorFaceArray_dual[f][I][J][K][2]->DotProduct(temp4);
    return EBot;
}

inline double EIntegrationT(GridsPoints *****ptrArray, Vector3 ******ptrEVectorFaceArray_dual, int f, int i, int j, int k)
{
    int I = i;
    int J = j + 1;
    int K = k;
    Vector3 temp1 = ptrArray[f][I + 1][J][K]->Pos3().MinusProduct(ptrArray[f][I][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray[f][I + 1][J][K + 1]->Pos3().MinusProduct(ptrArray[f][I + 1][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray[f][I][J][K + 1]->Pos3().MinusProduct(ptrArray[f][I + 1][J][K + 1]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray[f][I][J][K]->Pos3().MinusProduct(ptrArray[f][I][J][K + 1]->Pos3()).ScaleProduct(-1.0);

    double ET = ptrEVectorFaceArray_dual[f][I][J][K][0]->DotProduct(temp1) +
                ptrEVectorFaceArray_dual[f][I + 1][J][K][2]->DotProduct(temp2) +
                ptrEVectorFaceArray_dual[f][I][J][K + 1][0]->DotProduct(temp3) +
                ptrEVectorFaceArray_dual[f][I][J][K][2]->DotProduct(temp4);
    return ET;
}

// Return a integration of circuit E towards outside ( on edge)
inline double EIntegrationL(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    int I = i_in;
    int J = j_in;
    int K = k_in;
    Vector3 temp1 = ptrArray_in[face_in][I][J][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][I][J + 1][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K + 1]->Pos3());
    Vector3 temp3 = ptrArray_in[face_in][I][J + 1][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J + 1][K + 1]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][I][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J + 1][K]->Pos3());

    double EL = ptrArray_in[face_in][I][J][K]->DE3().DotProduct(temp1.PlusProduct(temp4)) + ptrArray_in[face_in][I][J][K + 1]->DE3().DotProduct(temp1.PlusProduct(temp2)) + ptrArray_in[face_in][I][J + 1][K + 1]->DE3().DotProduct(temp2.PlusProduct(temp3)) + ptrArray_in[face_in][I][J + 1][K]->DE3().DotProduct(temp3.PlusProduct(temp4));
    return EL / 2.0;
}

inline double EIntegrationR(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    int I = i_in + 1;
    int J = j_in;
    int K = k_in;
    Vector3 temp1 = ptrArray_in[face_in][I][J][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray_in[face_in][I][J + 1][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K + 1]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray_in[face_in][I][J + 1][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J + 1][K + 1]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray_in[face_in][I][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J + 1][K]->Pos3()).ScaleProduct(-1.0);

    double ER = ptrArray_in[face_in][I][J][K]->DE3().DotProduct(temp1.PlusProduct(temp4)) + ptrArray_in[face_in][I][J][K + 1]->DE3().DotProduct(temp1.PlusProduct(temp2)) + ptrArray_in[face_in][I][J + 1][K + 1]->DE3().DotProduct(temp2.PlusProduct(temp3)) + ptrArray_in[face_in][I][J + 1][K]->DE3().DotProduct(temp3.PlusProduct(temp4));
    return ER / 2.0;
}

inline double EIntegrationF(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    int I = i_in;
    int J = j_in;
    int K = k_in + 1;
    Vector3 temp1 = ptrArray_in[face_in][I][J + 1][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray_in[face_in][I + 1][J + 1][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J + 1][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray_in[face_in][I + 1][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J + 1][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray_in[face_in][I][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J][K]->Pos3()).ScaleProduct(-1.0);

    double EF = ptrArray_in[face_in][I][J][K]->DE3().DotProduct(temp1.PlusProduct(temp4)) + ptrArray_in[face_in][I][J + 1][K]->DE3().DotProduct(temp1.PlusProduct(temp2)) + ptrArray_in[face_in][I + 1][J + 1][K]->DE3().DotProduct(temp2.PlusProduct(temp3)) + ptrArray_in[face_in][I + 1][J][K]->DE3().DotProduct(temp3.PlusProduct(temp4));
    return EF / 2.0;
}

inline double EIntegrationBack(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    int I = i_in;
    int J = j_in;
    int K = k_in;
    Vector3 temp1 = ptrArray_in[face_in][I][J + 1][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][I + 1][J + 1][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J + 1][K]->Pos3());
    Vector3 temp3 = ptrArray_in[face_in][I + 1][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J + 1][K]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][I][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J][K]->Pos3());

    double EB = ptrArray_in[face_in][I][J][K]->DE3().DotProduct(temp1.PlusProduct(temp4)) + ptrArray_in[face_in][I][J + 1][K]->DE3().DotProduct(temp1.PlusProduct(temp2)) + ptrArray_in[face_in][I + 1][J + 1][K]->DE3().DotProduct(temp2.PlusProduct(temp3)) + ptrArray_in[face_in][I + 1][J][K]->DE3().DotProduct(temp3.PlusProduct(temp4));
    return EB / 2.0;
}

inline double EIntegrationT(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    int I = i_in;
    int J = j_in + 1;
    int K = k_in;
    Vector3 temp1 = ptrArray_in[face_in][I + 1][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray_in[face_in][I + 1][J][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J][K]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray_in[face_in][I][J][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J][K + 1]->Pos3()).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray_in[face_in][I][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K + 1]->Pos3()).ScaleProduct(-1.0);

    double ET = ptrArray_in[face_in][I][J][K]->DE3().DotProduct(temp1.PlusProduct(temp4)) + ptrArray_in[face_in][I + 1][J][K]->DE3().DotProduct(temp1.PlusProduct(temp2)) + ptrArray_in[face_in][I + 1][J][K + 1]->DE3().DotProduct(temp2.PlusProduct(temp3)) + ptrArray_in[face_in][I][J][K + 1]->DE3().DotProduct(temp3.PlusProduct(temp4));
    return ET / 2.0;
}

inline double EIntegrationBot(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    int I = i_in;
    int J = j_in;
    int K = k_in;
    Vector3 temp1 = ptrArray_in[face_in][I + 1][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K]->Pos3());
    Vector3 temp2 = ptrArray_in[face_in][I + 1][J][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J][K]->Pos3());
    Vector3 temp3 = ptrArray_in[face_in][I][J][K + 1]->Pos3().MinusProduct(ptrArray_in[face_in][I + 1][J][K + 1]->Pos3());
    Vector3 temp4 = ptrArray_in[face_in][I][J][K]->Pos3().MinusProduct(ptrArray_in[face_in][I][J][K + 1]->Pos3());

    double EBot = ptrArray_in[face_in][I][J][K]->DE3().DotProduct(temp1.PlusProduct(temp4)) + ptrArray_in[face_in][I + 1][J][K]->DE3().DotProduct(temp1.PlusProduct(temp2)) + ptrArray_in[face_in][I + 1][J][K + 1]->DE3().DotProduct(temp2.PlusProduct(temp3)) + ptrArray_in[face_in][I][J][K + 1]->DE3().DotProduct(temp3.PlusProduct(temp4));
    return EBot / 2.0;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Return a integration of circuit dB towards outside
// i, j, k are the index of cells in dual cell array
inline double dBIntegrationL(Vector3 ******ptrArrayB_main, Vector3 *****ptrArray_dual, int face, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray_dual[face][I][J][K + 1]->MinusProduct(*ptrArray_dual[face][I][J][K]);
    Vector3 temp2 = ptrArray_dual[face][I][J + 1][K + 1]->MinusProduct(*ptrArray_dual[face][I][J][K + 1]);
    Vector3 temp3 = ptrArray_dual[face][I][J + 1][K]->MinusProduct(*ptrArray_dual[face][I][J + 1][K + 1]);
    Vector3 temp4 = ptrArray_dual[face][I][J][K]->MinusProduct(*ptrArray_dual[face][I][J + 1][K]);

    double curldBL = ptrArrayB_main[face][I][J][K + 1][2]->DotProduct(temp1) + ptrArrayB_main[face][I][J + 1][K + 1][1]->DotProduct(temp2) + ptrArrayB_main[face][I][J + 1][K + 1][2]->DotProduct(temp3) + ptrArrayB_main[face][I][J + 1][K][1]->DotProduct(temp4);

    return curldBL;
}

inline double dBIntegrationR(Vector3 ******ptrArrayB_main, Vector3 *****ptrArray_dual, int face, int i, int j, int k)
{
    int I = i + 1;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray_dual[face][I][J][K + 1]->MinusProduct(*ptrArray_dual[face][I][J][K]).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray_dual[face][I][J + 1][K + 1]->MinusProduct(*ptrArray_dual[face][I][J][K + 1]).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray_dual[face][I][J + 1][K]->MinusProduct(*ptrArray_dual[face][I][J + 1][K + 1]).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray_dual[face][I][J][K]->MinusProduct(*ptrArray_dual[face][I][J + 1][K]).ScaleProduct(-1.0);

    double curldBR = ptrArrayB_main[face][I][J][K + 1][2]->DotProduct(temp1) + ptrArrayB_main[face][I][J + 1][K + 1][1]->DotProduct(temp2) + ptrArrayB_main[face][I][J + 1][K + 1][2]->DotProduct(temp3) + ptrArrayB_main[face][I][J + 1][K][1]->DotProduct(temp4);
    return curldBR;
}

inline double dBIntegrationBack(Vector3 ******ptrArrayB_main, Vector3 *****ptrArray_dual, int face, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray_dual[face][I][J + 1][K]->MinusProduct(*ptrArray_dual[face][I][J][K]);
    Vector3 temp2 = ptrArray_dual[face][I + 1][J + 1][K]->MinusProduct(*ptrArray_dual[face][I][J + 1][K]);
    Vector3 temp3 = ptrArray_dual[face][I + 1][J][K]->MinusProduct(*ptrArray_dual[face][I + 1][J + 1][K]);
    Vector3 temp4 = ptrArray_dual[face][I][J][K]->MinusProduct(*ptrArray_dual[face][I + 1][J][K]);

    double curldBBack = ptrArrayB_main[face][I][J + 1][K][1]->DotProduct(temp1) + ptrArrayB_main[face][I + 1][J + 1][K][0]->DotProduct(temp2) + ptrArrayB_main[face][I + 1][J + 1][K][1]->DotProduct(temp3) + ptrArrayB_main[face][I + 1][J][K][0]->DotProduct(temp4);
    return curldBBack;
}

inline double dBIntegrationF(Vector3 ******ptrArrayB_main, Vector3 *****ptrArray_dual, int face, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k + 1;
    Vector3 temp1 = ptrArray_dual[face][I][J + 1][K]->MinusProduct(*ptrArray_dual[face][I][J][K]).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray_dual[face][I + 1][J + 1][K]->MinusProduct(*ptrArray_dual[face][I][J + 1][K]).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray_dual[face][I + 1][J][K]->MinusProduct(*ptrArray_dual[face][I + 1][J + 1][K]).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray_dual[face][I][J][K]->MinusProduct(*ptrArray_dual[face][I + 1][J][K]).ScaleProduct(-1.0);

    double curldBF = ptrArrayB_main[face][I][J + 1][K][1]->DotProduct(temp1) + ptrArrayB_main[face][I + 1][J + 1][K][0]->DotProduct(temp2) + ptrArrayB_main[face][I + 1][J + 1][K][1]->DotProduct(temp3) + ptrArrayB_main[face][I + 1][J][K][0]->DotProduct(temp4);
    return curldBF;
}

inline double dBIntegrationBot(Vector3 ******ptrArrayB_main, Vector3 *****ptrArray_dual, int face, int i, int j, int k)
{
    int I = i;
    int J = j;
    int K = k;
    Vector3 temp1 = ptrArray_dual[face][I + 1][J][K]->MinusProduct(*ptrArray_dual[face][I][J][K]);
    Vector3 temp2 = ptrArray_dual[face][I + 1][J][K + 1]->MinusProduct(*ptrArray_dual[face][I + 1][J][K]);
    Vector3 temp3 = ptrArray_dual[face][I][J][K + 1]->MinusProduct(*ptrArray_dual[face][I + 1][J][K + 1]);
    Vector3 temp4 = ptrArray_dual[face][I][J][K]->MinusProduct(*ptrArray_dual[face][I][J][K + 1]);

    double curldBBot = ptrArrayB_main[face][I + 1][J][K][0]->DotProduct(temp1) + ptrArrayB_main[face][I + 1][J][K + 1][2]->DotProduct(temp2) + ptrArrayB_main[face][I + 1][J][K + 1][0]->DotProduct(temp3) + ptrArrayB_main[face][I][J][K + 1][2]->DotProduct(temp4);
    return curldBBot;
}

inline double dBIntegrationT(Vector3 ******ptrArrayB_main, Vector3 *****ptrArray_dual, int face, int i, int j, int k)
{
    int I = i;
    int J = j + 1;
    int K = k;
    Vector3 temp1 = ptrArray_dual[face][I + 1][J][K]->MinusProduct(*ptrArray_dual[face][I][J][K]).ScaleProduct(-1.0);
    Vector3 temp2 = ptrArray_dual[face][I + 1][J][K + 1]->MinusProduct(*ptrArray_dual[face][I + 1][J][K]).ScaleProduct(-1.0);
    Vector3 temp3 = ptrArray_dual[face][I][J][K + 1]->MinusProduct(*ptrArray_dual[face][I + 1][J][K + 1]).ScaleProduct(-1.0);
    Vector3 temp4 = ptrArray_dual[face][I][J][K]->MinusProduct(*ptrArray_dual[face][I][J][K + 1]).ScaleProduct(-1.0);

    double curldBT = ptrArrayB_main[face][I + 1][J][K][0]->DotProduct(temp1) + ptrArrayB_main[face][I + 1][J][K + 1][2]->DotProduct(temp2) + ptrArrayB_main[face][I + 1][J][K + 1][0]->DotProduct(temp3) + ptrArrayB_main[face][I][J][K + 1][2]->DotProduct(temp4);
    return curldBT;
}

//************************************************************************
//************************************************************************
// FUNCTION // Return a volume for a cell. Input (face, i, j, k) and return
// a double. The value is (average of back and front area) * (norm of vector
// between back and front face vector)

inline double CellVolume(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    double normV = FaceFieldVectorF(ptrArray_in, face_in, i_in, j_in, k_in, 'P').PlusProduct(FaceFieldVectorBack(ptrArray_in, face_in, i_in, j_in, k_in, 'P')).norm() / 2.0;

    double avgArea = AreaVectorBack(ptrArray_in, face_in, i_in, j_in, k_in).PlusProduct(AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).ScaleProduct(-1.0)).norm() / 2.0;
    /*    std::cout << " test" << AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).x() << " " 
                << AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).y() << " "
                << AreaVectorF(ptrArray_in, face_in, i_in, j_in, k_in).z() << " "
                << normV * avgArea;
    int pause;
    std::cin >> pause;
*/
    return normV * avgArea;
}

//************************************************************************
// FUNCTION // Return vector3
// index f, i, j, k are for ptrArrayGrids
inline Vector3 CellB3(GridsPoints *****ptrArray, int f, int i, int j, int k)
{
    double x = (ptrArray[f][i][j][k]->B3().x() + ptrArray[f][i + 1][j][k]->B3().x() +
                ptrArray[f][i][j + 1][k]->B3().x() + ptrArray[f][i + 1][j + 1][k]->B3().x() +
                ptrArray[f][i][j][k + 1]->B3().x() + ptrArray[f][i + 1][j][k + 1]->B3().x() +
                ptrArray[f][i][j + 1][k + 1]->B3().x() + ptrArray[f][i + 1][j + 1][k + 1]->B3().x()) *0.125F;
    double y = (ptrArray[f][i][j][k]->B3().y() + ptrArray[f][i + 1][j][k]->B3().y() +
                ptrArray[f][i][j + 1][k]->B3().y() + ptrArray[f][i + 1][j + 1][k]->B3().y() +
                ptrArray[f][i][j][k + 1]->B3().y() + ptrArray[f][i + 1][j][k + 1]->B3().y() +
                ptrArray[f][i][j + 1][k + 1]->B3().y() + ptrArray[f][i + 1][j + 1][k + 1]->B3().y()) *0.125F;
    double z = (ptrArray[f][i][j][k]->B3().z() + ptrArray[f][i + 1][j][k]->B3().z() +
                ptrArray[f][i][j + 1][k]->B3().z() + ptrArray[f][i + 1][j + 1][k]->B3().z() +
                ptrArray[f][i][j][k + 1]->B3().z() + ptrArray[f][i + 1][j][k + 1]->B3().z() +
                ptrArray[f][i][j + 1][k + 1]->B3().z() + ptrArray[f][i + 1][j + 1][k + 1]->B3().z()) *0.125F;
    Vector3 temp = Vector3(x, y, z);
    return temp;
}

//************************************************************************
//************************************************************************
// FUNCTION
// As in the updating curlField and gradientPe array, some variables are
// repeating calculating, it is suitable to put them in one function.
// Therefore, we need three matrix of curlB, curlE, and gradientPe.
// Assume they are curlB, curlE and gradPe, respectively.

void updateCellMatrix(Vector3 ****curlB_in, Vector3 ****curlE_in,
                      Vector3 ****gradPe_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Value gradient field of Pe.

Vector3 ***ValueGradientPe(Vector3 ***gradientArray_in, double ***ptrVolumeCellArray_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Value gradient field of Pe.
// gradientArray_in is in size of ( fsize+2 * fsize+2 * fsize) with vector3
// ptrVolumeCellArray is in size of ( fsize+2 * fsize+2 * fsize) with double
// Pe = n k T, in which n is the number density, k is the boltzmann constant, and T is the Te
Vector3 *****ValueGradient(GridsPoints *****ptrArray_in, Vector3 *****gradientArray_in, double ***ptrVolumeCellArray_in, char char_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Value the matrix field using finite volume method, put in the pointer
// of the MatrixField, value it, and return the pointer.
// Notice that the cell at corners should be absent in calculation.
Vector3 ***ValueCurlField(Vector3 ***curlArray_in, double ***ptrVolumeCellArray_in, GridsPoints *****ptrArray_in, int face_in, char field_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density
// at each grids points
// (fieldGridsSize+1)*(fieldGridsSize+1)*(fieldGridsSize+1) for "double"

double ***VolumeGridsField(double ***ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density
// at each cell
// (fieldGridsSize * fieldGridsSize * fieldGridsSize) for "double" type

double ***VolumeCellsField(GridsPoints *****ptrArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density
// at each cell
// (fieldGridsSize * fieldGridsSize * fieldGridsSize) for "double" type
double ***VolumeWeightGridsField(GridsPoints *****ptrArray_in);

inline double WeightCalculation(double x, double y, double z)
{
    double det = 1.0 / (x * x + y * y) * (x * x + z * z) * (1.0 - y * y * z * z / (x * x + y * y) / (x * x + z * z)) / (x * x + y * y + z * z) / logRatio / logRatio * radius * radius * LMin * LMin / log(10.0) / log(10.0);
    //            std::cout << x << " " << y << " " << z << " " << det<< std::endl;
    return 1.0 / pow(det, 0.5);
}
inline double WeightCalculation(const Vector3 &pos3)
{
    double x = pos3.x();
    double y = pos3.y();
    double z = pos3.z();
    double det = 1.0 / (x * x + y * y) * (x * x + z * z) * (1.0 - y * y * z * z / (x * x + y * y) / (x * x + z * z)) / (x * x + y * y + z * z) / logRatio / logRatio * radius * radius * LMin * LMin / log(10.0) / log(10.0);
    //            std::cout << x << " " << y << " " << z << " " << det<< std::endl;
    return 1.0 / pow(det, 0.5);
}
//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)

void updateGrids_nocurrent(Vector3 ***gradPe_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// UpdateVe3
void UpdateVe3(Vector3 ***curlField_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)
// Update E at grids for ve ( with current)
void UpdateE3(Vector3 ***gradPe_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// UpdateB3 vased on faraday's law

void UpdateB3(Vector3 ***curlField_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.

void updateGrids_withcurrent(Vector3 ***ptrVectorCellArray_in, GridsPoints *****ptrArray_in, int face_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the temprature

void Titheridge_Te(GridsPoints *****ptrArray_in, int);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the PSD, potential

void PP_update_(GridsPoints *****ptrArray);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the gradient of normal of B

void GradBNorm(GridsPoints *****ptrArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Sec convection velocity due to a ideal convectional cell and a dipole
// magnetic field

void SetConvectionVel(GridsPoints *****ptrArray_in, int face_in, int i_in, int k_in, int j_in);

//************************************************************************
//************************************************************************
// Function
// Set velocity due to earth rotation
void SetRotationalVel(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in);

//************************************************************************
//************************************************************************
// Function
// Set initial condition

void SetInitialCondition(GridsPoints *****ptrArray_in, Vector3 ***ptrVectorCellArray_in, double ***ptrVolumeCellArray_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the grad|B| on the gridspoints

void UpdateGrad(Vector3 *****gradBNorm_in, GridsPoints *****ptrArray_in, char char_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
void PrintOutHdf5_const(GridsPoints *****ptrArray_in);
void PrintOutHdf5(GridsPoints *****ptrArray_in, int i_in);
//************************************************************************
//************************************************************************
// FUNCTION
// Printout the particles
void PrintOutHdf5_Particles(int timeline,
                            vector<Particles> &ptrParticlesList_H,
                            vector<Particles> &ptrParticlesList_He,
                            vector<Particles> &ptrParticlesList_O);

void PrintOutHdf5_Grids(int timeline,
                        GridsPoints *****ptrArray);

void PrintOutHdf5_Particles_Grids(int timeline_in,
                                  vector<Particles> &ptrParticlesList_H,
                                  vector<Particles> &ptrParticlesList_He,
                                  vector<Particles> &ptrParticlesList_O,
                                  GridsPoints *****ptrArray);

void PrintOutHdf5_Particles_Grids(int timeline_in,
                                  GridsCells ****ptrArrayCells,
                                  GridsPoints *****ptrArray);

void PrintOutHdf5_DivB(int timeline,
                       GridsPoints *****ptrArray,
                       Vector3 *****ptrBVectorFaceArray,
                       double ***ptrVolumeCellArray,
                       Vector3 *****ptrDivBVectorCellArray,
                       Vector3 *****ptrposVectorCellArray);
//  FOr velDist plot
void PrintOutHdf5Cells(GridsCells ****ptrArrayCells, int timeline, int h5FileCheck);
//************************************************************************
//  FOr velDist plot,
void PrintOutHdf5Cells(GridsCells ****ptrArrayCells, int timeline);
//************************************************************************
//************************************************************************
// Read particles and grids data
int ReadSavedData(GridsPoints *****ptrArray,
                  vector<Particles> &ptrParticlesList_H,
                  vector<Particles> &ptrParticlesList_He,
                  vector<Particles> &ptrParticlesList_O,
                  vector<int> &ptrParticlesList_out_H,
                  vector<int> &ptrParticlesList_out_He,
                  vector<int> &ptrParticlesList_out_O);

int ReadSavedData(GridsPoints *****ptrArray,
                  GridsCells ****ptrArrayCells);

//************************************************************************
//************************************************************************
// Read dat file for initialing bottom layers
inline int CountLines(char *filename)
{
    std::ifstream ReadFile;
    int n = 0;
    std::string tmp;
    ReadFile.open(filename, std::ios::in); //ios::in, only read
    if (ReadFile.fail())                   //
    {
        return 0;
    }
    else //
    {
        while (getline(ReadFile, tmp, '\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}

/*inline void ReadDataFileForBottom(char *filename, double data[][19][2])
{

    std::ifstream datfile; //("./GITMdata//trans_t021221_040000.DAT");
    datfile.open(filename, std::ios::in);
    using namespace std;
    if (!datfile)
    {
        cout << " Unable to open data file "<<filename;
        exit(1);
    }
    else
    {
        cout << " Open data file successfully \n";
    }

    int totalLine = CountLines(filename);
    int i = 0; // line
    int k = 0;
    char buffer[1024];
    while (!datfile.eof() && i < totalLine) // the last line won't copy twice
    //while(getline(datfile,buffer))
    {
        double test[19];
        //
        datfile.getline(buffer, 1024);
        sscanf(buffer, "%*f%*f%*f %*f %*f%*f%*f%*f %*f%*f%*f%*f%*f %*f%*f%*f %*f%*f %*f",
               &test[0],
               &test[1],
               &test[2],
               &test[3],
               &test[4],
               &test[5],
               &test[6],
               &test[7],
               &test[8],
               &test[9],
               &test[10],
               &test[11],
               &test[12],
               &test[13],
               &test[14],
               &test[15],
               &test[16],
               &test[17],
               &test[18]);
        data[i][0][k] = test[0];
        data[i][1][k] = test[1];
        data[i][2][k] = test[2];
        data[i][3][k] = test[3];
        data[i][4][k] = test[4];
        data[i][5][k] = test[5];
        data[i][6][k] = test[6];
        data[i][7][k] = test[7];
        data[i][8][k] = test[8];
        data[i][9][k] = test[9];
        data[i][10][k] = test[10];
        data[i][11][k] = test[11];
        data[i][12][k] = test[12];
        data[i][13][k] = test[13];
        data[i][14][k] = test[14];
        data[i][15][k] = test[15];
        data[i][16][k] = test[16];
        data[i][17][k] = test[17];
        data[i][18][k] = test[18];
        i++;
        //Lontitude
        //Latitude
        //altitude
        //
        //'[O_4SP_!U+!N]'O(4SP)+ density
        //'[O!D2!U+!N]'O2+ density
        //'[N!D2!U+!N]'N2+ density
        //'[N!U+!N]'N+ density
        //'[NO!U+!N]'NO+ density
        //'[O(!U2!ND)!U+!N]'O(2D)+ density
        //'[O(!U2!NP)!U+!N]'O(2P)+ density
        //'[H!U+!N]'H+ density
        //'[He!U+!N]'He+ density
        //'[e-]'electron density
        //
        //'V!Di!N(east)'  eastward ion velocity
        //'V!Di!N(north)' northward ion velocity
        //'V!Di!N(up)'    upward ion velocity
        //
        //'eTemperature' electron temperature
        //'iTemperature' ion temperature
        //
        //'potential"
    }
    datfile.close();
    std::cout << " Finished loading data \n";
}*/

// filename coding
// 1. given a (x, y, z) point, trans into (long, lat) in GMLT
// 2. trans into (LONG, LAT) in GEO coord
// 3. find out related index in IRI and TIEGCM array
// 4- to achieve 3, we need a array to store IRI and TIEGCM data
// 4-1 [201][201][4] for IRI, [LONG every 1.8 degree][LAT every 0.9 degree][Ti, NO+, NHe+, NH+]
// 4-2 [73][38][4], [LONG every 5 degree][LAT every 5 degree][PHI, v_long, v_lat, v_ip(vertical)] note tiegcm data type
// 5- to achieve 3, we need a array to store long and lat for every girds because they're constant
// Therefore, the procedure is as following:
// 1. Read IRI file, store data into [201][201][4] array, unchangable
// 2. Read TIEGCM file, store data into [73][38][4] array, unchangable
// 3. Create a [totalface][fieldsize+1][fieldsize+1][fieldsize+1][long][lat] to store related/ field line
// along point(long, lat) in GM->GEO coordinates, unchangable (maybe not because the calculation is fast)
// 4. With a timestep, change approciated (long) to find related index position in IRI and TIEGCM array
// 5. Interpolate the value from the four vertex
inline void InputDataFileFromExternalData(GridsPoints *****ptrArrayGrids,
                                          int timeline,
                                          string workdir,
                                          vector<vector<vector<double>>> *iriData,
                                          vector<vector<vector<double>>> *tiegcmData)
{
    //double pi = acos(-1.0);
    string iriDataPath = workdir + "/ExternalData/INPUT/IRI.txt";
    string tiegcmDataPath = workdir + "/ExternalData/INPUT/TIEGCM.txt";

    std::ifstream iriFile, tiegcmFile;
    //
    iriFile.open(iriDataPath, std::ios::in);
    tiegcmFile.open(tiegcmDataPath, std::ios::in);
    if (!iriFile || !tiegcmFile)
    {
        std::cout << " Unable to open lower boundary data file" << iriDataPath <<" or " << tiegcmDataPath <<endl;
        exit(1);
    }
    else
    {
        std::cout << " Open lower bounday data files successfully \n";
    }
    // apply continue memory to store IRI data
    // iri: Ti, H+, He+, O+
    int rowIri = 201;
    int columnIri = 201;
    //char bufferIri[2048];
    //double tempIri[4];
    double iri[16];
    int    k;

    for (int i = 0; i < rowIri; i++)
    {
        for (int j = 0; j < columnIri; j++)
        {
            for (k=0; k< 16; k++) iriFile >> iri[k];

            (*iriData)[i][j][0] = iri[5]; // Ti
            (*iriData)[i][j][1] = iri[8]; // NH+
            (*iriData)[i][j][2] = iri[9]; // NHe+
            (*iriData)[i][j][3] = iri[7]; // NO+

            /*iriFile.getline(bufferIri, 2048);
            sscanf(bufferIri, "%*f%*f%*f%*f%*f %*f%*f%*f%*f%*f %*f%*f%*f%*f%*f %*f",
                   &tempIri[0],               // Ti
                   &tempIri[1],               //NO+
                   &tempIri[2],               // NH+
                   &tempIri[3]);              // NHe+
            (*iriData)[i][j][0] = tempIri[0]; // Ti
            (*iriData)[i][j][1] = tempIri[2]; // NH+
            (*iriData)[i][j][2] = tempIri[3]; // NHe+
            (*iriData)[i][j][3] = tempIri[1]; // NO+*/
        }
    }
    // TIEGCM
    int rowTiegcm = 38; // lat, note tiegcm data type!
    int columeTiegcm = 73;
    //char bufferTiegcm[2048];
    //double tempTiegcm[4];
    double tiegcm[47];

    for (int i = 0; i < rowTiegcm; i++)
    {
        for (int j = 0; j < columeTiegcm; j++)
        {
            for (k=0; k<47; k++) tiegcmFile >> tiegcm[k];

            (*tiegcmData)[i][j][0] = tiegcm[17];
            (*tiegcmData)[i][j][1] = tiegcm[18];
            (*tiegcmData)[i][j][2] = tiegcm[19];
            (*tiegcmData)[i][j][3] = tiegcm[20];

            /*tiegcmFile.getline(bufferTiegcm, 2048);
            sscanf(bufferTiegcm, "%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f %*f%*f%*f%*f%*f%*f%*f%*f%*f%*f %*f%*f%*f%*f%*f%*f%*f%*f%*f%*f %*f%*f%*f%*f%*f%*f%*f%*f%*f%*f %*f%*f%*f%*f%*f%*f%*f",
                   &tempTiegcm[0],  // PHI potential
                   &tempTiegcm[1],  //vi_lon
                   &tempTiegcm[2],  //vi_lat
                   &tempTiegcm[3]); // vi_vertical
            (*tiegcmData)[i][j][0] = tempTiegcm[0];
            (*tiegcmData)[i][j][1] = tempTiegcm[1];
            (*tiegcmData)[i][j][2] = tempTiegcm[2];
            (*tiegcmData)[i][j][3] = tempTiegcm[3];*/
        }
    }
}
inline void UpdateDatafileFromExternalDataForBottomBoundary(GridsPoints *****ptrArrayGrids,
                                                            int timeline,
                                                            string workdir,
                                                            vector<vector<vector<double>>> *iriData,
                                                            vector<vector<vector<double>>> *tiegcmData)
{
    //double pi = acos(-1.0);

    double xgeo, ygeo, zgeo;
    double xmag, ymag, zmag;
    double x, y, z;
    double ra, thet, phis;
    int a, b;
    double aDouble, bDouble;
    double ti, nH, nHe, nO;
    double phi, vi_lat, vi_vert; //, vi_lon
    double vxGEO, vyGEO, vzGEO, vxMAG, vyMAG, vzMAG;

    // weight for interpolation
    double i1, i2, j1, j2;
    double w1, w2, w3, w4;

    // f, i, j, k are index of grids points
//#pragma omp parallel for collapse(4)
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    //std::cout << " index " << f << " " << i << " " << j << " " << k << "\n";
                    // private for openmp
                    // (xyz) in GM -> (long,lat) in GEO
                    // get (xmag, ymag, zmag)
                    x = ptrArrayGrids[f][i][j][k]->Pos3().x();
                    y = ptrArrayGrids[f][i][j][k]->Pos3().y();
                    z = ptrArrayGrids[f][i][j][k]->Pos3().z();
                    // trans into (xmag, ymag, zmag) on earth at 500km  // why 400 km does not work
                    DipoleRelatedOnEarthPolar(x, y, z, xmag, ymag, zmag, AltitudeMin-50.0);
                    // trans into (xgeo, ygeo, zgeo)
                    GEOMAG_08(xgeo, ygeo, zgeo, xmag, ymag, zmag, -1);
                    // trans into (thet, phis) in GEO
                    SPHCAR_08(ra, thet, phis, xgeo, ygeo, zgeo, -1);
                    // rotation by earth, phi changes (decrease)
                    phis = phis - omega_earth * timeline * tstep;
                    while (phis < 0.0)
                        phis += 2 * pi;
                    // bottom layer, for only one layer profile, we set k =0 and 1 are the same density
                    if (k < tempGridsCellLevelBot)
                    {
                        // locate the index in IRI and TIEGCM array by "thet" and "phis", 
                        // and do an approximate interpolation
                        // IRI: ti, nH, nHe, nO
                        aDouble = thet * 180.0 / pi / 0.9;
                        bDouble = phis * 180.0 / pi / 1.8;
                        a = static_cast<int>(aDouble);
                        b = static_cast<int>(bDouble);
                        // weight calculation
                        i1 = aDouble - a;
                        i2 = 1.0 - i1;
                        j1 = bDouble - b;
                        j2 = 1.0 - j1;
                        w1 = i2 * j2;
                        w2 = i1 * j2;
                        w3 = i2 * j1;
                        w4 = i1 * j1;
                        //
                        ti = (*iriData)[a][b][0] * w1 + (*iriData)[a][b + 1][0] * w3 + (*iriData)[a + 1][b][0] * w2 + (*iriData)[a + 1][b + 1][0] * w4;
                        nH = (*iriData)[a][b][1] * w1 + (*iriData)[a][b + 1][1] * w3 + (*iriData)[a + 1][b][1] * w2 + (*iriData)[a + 1][b + 1][1] * w4;
                        nHe = (*iriData)[a][b][2] * w1 + (*iriData)[a][b + 1][2] * w3 + (*iriData)[a + 1][b][2] * w2 + (*iriData)[a + 1][b + 1][2] * w4;
                        nO = (*iriData)[a][b][3] * w1 + (*iriData)[a][b + 1][3] * w3 + (*iriData)[a + 1][b][3] * w2 + (*iriData)[a + 1][b + 1][3] * w4;
                        // corrected value
                        nH = nH * 1000.0;
                        nHe = nHe *1000.0;
                        nO = nO *1000.0;
                        //
                        // TIEGCM: phi, (vxMAG, vyMAG, vzMAG), + 0.5 for a is because the data type of TIEGCM
                        aDouble = thet * 180.0 / pi / 5.0 + 0.5;
                        bDouble = phis * 180.0 / pi / 5.0;
                        a = static_cast<int>(aDouble);
                        b = static_cast<int>(bDouble);
                        // weight calculation
                        i1 = aDouble - a;
                        i2 = 1.0 - i1;
                        j1 = bDouble - b;
                        j2 = 1.0 - j1;
                        w1 = i2 * j2;
                        w2 = i1 * j2;
                        w3 = i2 * j1;
                        w4 = i1 * j1;
                        phi = (*tiegcmData)[a][b][0] * w1 + (*tiegcmData)[a][b + 1][0] * w3 + (*tiegcmData)[a + 1][b][0] * w2 + (*tiegcmData)[a + 1][b + 1][0] * w4;
                        //vi_lon = (*tiegcmData)[a][b][1] * w1 + (*tiegcmData)[a][b + 1][1] * w3 + (*tiegcmData)[a + 1][b][1] * w2 + (*tiegcmData)[a + 1][b + 1][1] * w4;
                        vi_lat = (*tiegcmData)[a][b][2] * w1 + (*tiegcmData)[a][b + 1][2] * w3 + (*tiegcmData)[a + 1][b][2] * w2 + (*tiegcmData)[a + 1][b + 1][2] * w4;
                        vi_vert = (*tiegcmData)[a][b][3] * w1 + (*tiegcmData)[a][b + 1][3] * w3 + (*tiegcmData)[a + 1][b][3] * w2 + (*tiegcmData)[a + 1][b + 1][3] * w4;
                        // trans vi in (long, lat, vertical) into (vxGEO, vyGEO, vzGEO)
                        BCARSPRR_08(thet, phis, vi_vert, vi_lat, vi_lat, vxGEO, vyGEO, vzGEO);
                        // trans velocity (vxGEO, vyGEO, vzGEO) into (vxMAG, vyMAG, vzMAG)
                        GEOMAG_08(vxGEO, vyGEO, vzGEO, vxMAG, vyMAG, vzMAG, 1);
//#pragma omp critical
                        {
                            // value the related values
                            ptrArrayGrids[f][i][j][k]->SetTemperatureIons(ti);                // Ti
                            ptrArrayGrids[f][i][j][k]->Density_H(nH);                         // nH
                            ptrArrayGrids[f][i][j][k]->Density_He(nHe);                        // nHe
                            ptrArrayGrids[f][i][j][k]->Density_O(nO);                         // nO
                            ptrArrayGrids[f][i][j][k]->SetPotential(phi);                     // PHI
                            ptrArrayGrids[f][i][j][k]->SetVel3(Vector3(vxMAG, vyMAG, vzMAG)); // vi
                        }
                        //
                    }
                    else
                    {
                        // TIEGCM: phi, only!,  - 0.5 for a is because the data type of TIEGCM
                        aDouble = thet * 180.0 / pi / 5.0 + 0.5;
                        bDouble = phis * 180.0 / pi / 5.0;
                        a = static_cast<int>(aDouble);
                        b = static_cast<int>(bDouble);
                        // weight calculation
                        i1 = aDouble - a;
                        i2 = 1.0 - i1;
                        j1 = bDouble - b;
                        j2 = 1.0 - j1;
                        w1 = i2 * j2;
                        w2 = i1 * j2;
                        w3 = i2 * j1;
                        w4 = i1 * j1;
                        phi = (*tiegcmData)[a][b][0] * w1 + (*tiegcmData)[a][b + 1][0] * w3 + (*tiegcmData)[a + 1][b][0] * w2 + (*tiegcmData)[a + 1][b + 1][0] * w4;
//#pragma omp critical
{
                        ptrArrayGrids[f][i][j][k]->SetPotential(phi); // Potential       
                        if( k == tempGridsCellLevelBot)
                        {
                            // value the related values
                            ptrArrayGrids[f][i][j][k]->SetTemperatureIons(ptrArrayGrids[f][i][j][0]->temperature);                // Ti
                            ptrArrayGrids[f][i][j][k]->Density_H(ptrArrayGrids[f][i][j][0]->density_H);                         // nH
                            ptrArrayGrids[f][i][j][k]->Density_He(ptrArrayGrids[f][i][j][0]->density_He);                        // nHe
                            ptrArrayGrids[f][i][j][k]->Density_O(ptrArrayGrids[f][i][j][0]->density_O);                         // nO
                            ptrArrayGrids[f][i][j][k]->SetVel3(ptrArrayGrids[f][i][j][0]->v3); // vi
                        } 
}
                    }
                }
            }
        }
    }

    // std::string path = "./ExternalData/IRIDATA/pointdata_392696028509_20060101_010000";
    // //1
    // std::string::size_type iPos = path.find_last_of('/') + 1;
    // std::string filename = path.substr(iPos, path.length() - iPos);
    // std::cout << " 1 " << filename << std::endl;
    // //2
    // string name = filename.substr(0, filename.rfind("."));
    // cout << " 2 " << name << endl;
    // //3
    // string suffix_str = filename.substr(filename.find_last_of('.') + 1);
    // cout << " 2 " << suffix_str << endl;
    // //
    // std::ifstream ReadFile;
    // int n = 0;
    // std::string tmp;
    // ReadFile.open(filename, std::ios::in);
}
// Update potential based on external datafile 
// Assume the bottom layer/region has updated potentials
// and known the footprint (FP) points of each gridsPoints, interpolate potential
// from bottom layer/region grids to the FP of each gridsPoints.
// Assume the potentials are the same along B field lines
inline void UpdatePotentialAtGrids( GridsPoints *****_ptrArray,
                               foot3 ****_footArray)
{
    // face ,i, j, k are index of gridsPoints
    for(int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for( int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for( int k = 0; k < fieldsGridsSize+1; k++)
                {
                    // assume a "particle" at the position of FP
                    Particles unrealPar = Particles();
                    struct structg tempStr;
                    struct structPar tempStrPar;
                    unrealPar.SetPosPar(Vector3(_footArray[face][i][j][k].x(),
                                                _footArray[face][i][j][k].y(),
                                                _footArray[face][i][j][k].z()));
                    if(unrealPar.PosParticles().norm() < LMin)
                    {
                        std::cout << " can't locate footprint in simulation domain \n";
                        exit (1);
                    }
                    tempStr = unrealPar.InttoStrp1();
                    StructPar(tempStr, tempStrPar);
                    // get potential at position of unreal Particle
                    double tempPotential = unrealPar.PotentialAtUnrealPar(tempStrPar, _ptrArray);
                    // set potential
                    _ptrArray[face][i][j][k]->SetPotential(tempPotential);
                }
            }
        }
    }
}
//************************************************************************
//************************************************************************
// FUNCTION
// Create grids and cells

GridsPoints *****GridsCreation();

GridsPoints *****GridsCreation(GridsPoints *****ptrArray, int gridsSize);

GridsCells ****GridsCellsCreation(GridsCells *mem_GridsCells,
                                  GridsPoints *****ptrArray);

// Create a array of pointer to store position of dual cell grid point(corners)
// [face][fieldsGridsSize+2][fieldsGridsSize+2][fieldsGridsSize]
Vector3 *****GridsCreation_dual(GridsPoints *****ptrArray);

void InitializeTempGrids(GridsPoints *****ptrArray, GridsPoints *****ptrArray_bot, GridsPoints *****ptrArray_top, int gridsSize);
//************************************************************************
//************************************************************************
// FUNCTION // Set up a matrix to store the curl E or B for Faraday's Law
// and for Ampere's Law, or the gradient of Pe.
// The size of the matrix should be 1 smaller than
// the size of gridspoints in main doman which is a cubic, which is [fsize+2].
// Therefore, it is [fsize+2][fsize+2][fsize]
// For each face, 8 corner cell should be excluded. ?
// Notice that the curl E or B is at the center of each cell.
// The data structure is array of Vector3, which is created in heap. Return
// a pointer(may not need to be a smart pointer), and would not need to
// delete, or would be deleted as a smart pointer.
void VectorCellField(Vector3 ***&cellArray);

void VectorCellField_Vel(Vector3 ***&cellArray);

void VectorCellField_Grad(Vector3 ***&cellArray);

void DEL_VectorCellField(Vector3 ***&cellArray);

void test_VectorCellField(Vector3 **&cellArray);

//************************************************************************
//************************************************************************
// Create Cell centered field array
Vector3 *****PtrVectorCellArray(Vector3 *mem_gradPeVectorCellArray);

//************************************************************************
//************************************************************************
// FUNCTION
// finish culmulating and average the density and velocity
void CalculatingAveragedPhoVatGrids(GridsPoints *****ptrArray_in,
                                    double ***ptrVolumeGridArray_in,
                                    int updateInfoPeriod_in,
                                    Vector3 ******ptrVelWeightGridsArray,
                                    double ******ptrMassWeightGridsArray);

void CalculatingAveragedPhoVatGrids(GridsPoints *****ptrArray_in,
                                    double ***ptrVolumeGridArray_in,
                                    int updateInfoPeriod_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Set zero for pho and v at each points
void ResetPhoVatGrids(GridsPoints *****ptrArray_in);

//************************************************************************
//************************************************************************
// Function
// initial the bot boundary for the  velocity of magnetic field line
void SetRotationalVeBotBoundary(GridsPoints *****ptrArray_in,
                                int timeline_in);

void SetRotationalVeBotBoundary(GridsPoints *****ptrArray_in,
                                Vector3 *****ptrEVectorCellArray,
                                Vector3 *****ptrVeleVectorCellArray,
                                Vector3 *****ptrGradVectorCellArray,
                                int timeline_in);

//************************************************************************
//************************************************************************
// Function
// initial the top boundary for the  velocity of magnetic field line
void SetConvectionVelTopBoundary(GridsPoints *****ptrArray_in, int timeline_in);

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the B on face
// The size of this array is [direction * face * (fsize+1) * (fsize+1) * (fsize+1)]
Vector3 *****BVectorFaceArray(GridsPoints *****ptrArray_in);

// The size of this array is [face * (fsize+2) * (fsize+2) * (fsize) * direction
Vector3 ******BVectorFaceArray(Vector3 ******ptrBVectorFaceArray_main);

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the B on face

// The size of this array is [face * (fsize+1) * (fsize+1) * (fsize) * direction
Vector3 ******EVectorFaceArray(Vector3 ******ptrEVectorFaceArray);

//************************************************************************
//************************************************************************
// Prerun 1.6 // Create const array of B at center of each cells
// [totalface * fsize+2 * fsize+2 * fsize +2]
Vector3 *****BVectorCellArray(GridsPoints *****ptrArray);

//************************************************************************
//************************************************************************
//
int ******BVectorFaceArray(int ******ptrBCheckFaceArray_main);
int ******EVectorFaceArray(int ******ptrECheckFaceArray_dual);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the B on each face from BVectorFaceArray
void BVectorFaceArrayUpdate(GridsPoints *****ptrArray,
                            Vector3 ******ptrBVectorFaceArray_main,
                            Vector3 ******ptrBVectorFaceArray_main_backup,
                            int ******ptrBCheckFaceArray_main);

//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the curl B at the center of each cell
Vector3 *****CurlBCellArray(GridsPoints *****ptrArray,
                            Vector3 *****ptrVectorCellArray,
                            Vector3 *****ptrBVectorCellArray,
                            double ***ptrVolumeCellArray);

//************************************************************************
//************************************************************************
void UpdateECellArray(GridsPoints *****ptrArray_in,
                      Vector3 *****ptrEVectorCellArray,
                      Vector3 *****ptrGradVectorCellArray);
// FUNCTION
// Calculate the curl B at the center of each cell
void UpdateECellArray(GridsPoints *****ptrArray,
                      Vector3 *****ptrEVectorCellArray,
                      Vector3 *****ptrVeleVectorCellArray,
                      Vector3 *****ptrVectorCellArray,
                      Vector3 *****ptrGradVectorCellArray);

// ******************************************************
// Update the dB and E on each grids
void BVectorGridsArrayUpdate(GridsPoints *****ptrArray,
                             Vector3 ******ptrBVectorFaceArray);

void EVectorGridsArrayUpdate(GridsPoints *****ptrArray,
                             Vector3 ******ptrEFace_dual);

// ******************************************************
// Create grids array for openmp weighting calculation
inline Vector3 ******PtrVelWeightGridsArray(Vector3 *mem_VelWeightGridsArray)
{
    Vector3 ******ptrVelWeightGridsArray = new Vector3 *****[3];
    for (int ionType = 0; ionType < 3; ionType++)
    {
        ptrVelWeightGridsArray[ionType] = new Vector3 ****[totalFace];
        for (int face = 0; face < totalFace; face++)
        {
            ptrVelWeightGridsArray[ionType][face] = new Vector3 ***[8];
            for (int pos = 0; pos < 8; pos++)
            {
                ptrVelWeightGridsArray[ionType][face][pos] = new Vector3 **[fieldsGridsSize / 2 + 1];
                for (int i = 0; i < fieldsGridsSize / 2 + 1; i++)
                {
                    ptrVelWeightGridsArray[ionType][face][pos][i] = new Vector3 *[fieldsGridsSize / 2 + 1];
                    for (int j = 0; j < fieldsGridsSize / 2 + 1; j++)
                    {
                        ptrVelWeightGridsArray[ionType][face][pos][i][j] = mem_VelWeightGridsArray + ionType * totalFace * 8 * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + face * 8 * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + pos * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + i * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + j * (fieldsGridsSize / 2 + 1);

                        for (int k = 0; k < fieldsGridsSize / 2 + 1; k++)
                            ptrVelWeightGridsArray[ionType][face][pos][i][j][k] = Vector3(0.0, 0.0, 0.0);
                    }
                }
            }
        }
    }
    return ptrVelWeightGridsArray;
};

inline double ******PtrMassWeightGridsArray(double *mem_MassWeightGridsArray)
{
    double ******ptrMassWeightGridsArray = new double *****[3];
    for (int ionType = 0; ionType < 3; ionType++)
    {
        ptrMassWeightGridsArray[ionType] = new double ****[totalFace];
        for (int face = 0; face < totalFace; face++)
        {
            ptrMassWeightGridsArray[ionType][face] = new double ***[8];
            for (int pos = 0; pos < 8; pos++)
            {
                ptrMassWeightGridsArray[ionType][face][pos] = new double **[fieldsGridsSize / 2 + 1];
                for (int i = 0; i < fieldsGridsSize / 2 + 1; i++)
                {
                    ptrMassWeightGridsArray[ionType][face][pos][i] = new double *[fieldsGridsSize / 2 + 1];
                    for (int j = 0; j < fieldsGridsSize / 2 + 1; j++)
                    {
                        ptrMassWeightGridsArray[ionType][face][pos][i][j] = mem_MassWeightGridsArray + ionType * totalFace * 8 * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + face * 8 * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + pos * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + i * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) + j * (fieldsGridsSize / 2 + 1);

                        for (int k = 0; k < fieldsGridsSize / 2 + 1; k++)
                            ptrMassWeightGridsArray[ionType][face][pos][i][j][k] = 0.0;
                    }
                }
            }
        }
    }
    return ptrMassWeightGridsArray;
};

// calculate curl B on face of E dual cell ( current j)
void EVectorFaceArrayUpdate(Vector3 *****ptrArray_dual,
                            Vector3 ******ptrBVectorFaceArray_main,
                            Vector3 ******ptrEVectorFaceArray_dual,
                            Vector3 ******ptrBVectorFaceArray_dual_backup,
                            int ******ptrECheckFaceArray_dual);
// calculate curl dB
Vector3 CurldBOnDualCellFace(Vector3 *****ptrArray_dual,
                             Vector3 ******ptrBFace_main,
                             int dir,
                             int face,
                             int i,
                             int j,
                             int k);
// Calculate ve at grid points of main cell
// type: for update_type == 0, no current case, type = 0 for parallel E, type = 1 for pependicular E
void UpdateVeGridsMain(GridsPoints *****ptrArray,
                       Vector3 ******ptrEVectorFaceArray_dual,
                       int type);

//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the temprature of ions on the gridspoints
void UpdateTempIons(GridsPoints *****ptrArray, GridsCells ****ptrArrayCells);

#endif