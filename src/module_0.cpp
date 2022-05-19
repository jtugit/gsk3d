#include <iostream>
#include <list>
#include <memory>
#include <string>
#include "parameters.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module_0.h"
#include <cmath>
#include "H5Cpp.h"
#include <bitset>
#include <vector>
#include "module_base.h"
#include "gridscells.h"
#include <fstream>

using std::cout;
using std::endl;
using std::make_shared;
using std::shared_ptr;

//************************************************************************
//************************************************************************
// Generate a grids of pointers which point to the
// objects of GridsPoints in heap.
//
//************************************************************************
//************************************************************************

GridsPoints *****GridsCreation()
{
    GridsPoints *****ptrArray;
    ptrArray = new GridsPoints ****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrArray[face] = new GridsPoints ***[fieldsGridsSize + 3];
        for (int i = 0; i <= fieldsGridsSize + 2; i++)
        {
            ptrArray[face][i] = new GridsPoints **[fieldsGridsSize + 3];
            for (int j = 0; j <= fieldsGridsSize + 2; j++)
            {
                ptrArray[face][i][j] = new GridsPoints *[fieldsGridsSize * grid_domain + 1];
            }
        }
    }
    // face 0 (in the center)
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j <= fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
            {
                ptrArray[0][i][j][k] = new GridsPoints();
                ptrArray[0][i][j][k]->InttoPos3(0, i, j, k);
                ptrArray[0][i][j][k]->XYZtoB(ptrArray[0][i][j][k]->Pos3());
                ptrArray[0][i][j][k]->XYZtoVel();
                ptrArray[0][i][j][k]->XYZtoE();
                ptrArray[0][i][j][k]->XYZtoDensity( 0);
                ptrArray[0][i][j][k]->SetStopSign(0);
                ptrArray[0][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // face 1 (on the right)
    // share len with face 0
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[1][1][j][k] = ptrArray[0][fieldsGridsSize + 1][j][k];
        }
    }

    for (int i = 2; i <= fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j <= fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
            {
                ptrArray[1][i][j][k] = new GridsPoints();
                ptrArray[1][i][j][k]->InttoPos3(1, i, j, k);
                ptrArray[1][i][j][k]->XYZtoB(ptrArray[1][i][j][k]->Pos3());
                ptrArray[1][i][j][k]->XYZtoVel();
                ptrArray[1][i][j][k]->XYZtoE();
                ptrArray[1][i][j][k]->XYZtoDensity( 0);
                ptrArray[1][i][j][k]->SetStopSign(0);
                ptrArray[1][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 2 (on the top)
    // share len with face 0
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[2][i][1][k] = ptrArray[0][i][fieldsGridsSize + 1][k];
        }
    }
    // share len with face 1
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[2][fieldsGridsSize + 1][j][k] = ptrArray[1][j][fieldsGridsSize + 1][k];
        }
    }
    for (int i = 1; i <= fieldsGridsSize; i++)
    {
        for (int j = 2; j <= fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
            {
                ptrArray[2][i][j][k] = new GridsPoints();
                ptrArray[2][i][j][k]->InttoPos3(2, i, j, k);
                ptrArray[2][i][j][k]->XYZtoB(ptrArray[2][i][j][k]->Pos3());
                ptrArray[2][i][j][k]->XYZtoVel();
                ptrArray[2][i][j][k]->XYZtoE();
                ptrArray[2][i][j][k]->XYZtoDensity( 0);
                ptrArray[2][i][j][k]->SetStopSign(0);
                ptrArray[2][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // face 4 (on the left)
    // share len with face 0
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[4][fieldsGridsSize + 1][j][k] = ptrArray[0][1][j][k];
        }
    }
    // share len with face 2
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[4][i][fieldsGridsSize + 1][k] = ptrArray[2][1][fieldsGridsSize + 2 - i][k];
        }
    }

    for (int i = 1; i <= fieldsGridsSize; i++)
    {
        for (int j = 1; j <= fieldsGridsSize; j++)
        {
            for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
            {
                ptrArray[4][i][j][k] = new GridsPoints();
                ptrArray[4][i][j][k]->InttoPos3(4, i, j, k);
                ptrArray[4][i][j][k]->XYZtoB(ptrArray[4][i][j][k]->Pos3());
                ptrArray[4][i][j][k]->XYZtoVel();
                ptrArray[4][i][j][k]->XYZtoE();
                ptrArray[4][i][j][k]->XYZtoDensity( 0);
                ptrArray[4][i][j][k]->SetStopSign(0);
                ptrArray[4][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 5 (on the bottom)
    // share len with face 0
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[5][i][fieldsGridsSize + 1][k] = ptrArray[0][i][1][k];
        }
    }
    // share len with face 1
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[5][fieldsGridsSize + 1][j][k] = ptrArray[1][fieldsGridsSize + 2 - j][1][k];
        }
    }
    // share len with face 4
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[5][1][j][k] = ptrArray[4][j][1][k];
        }
    }
    for (int i = 2; i <= fieldsGridsSize; i++)
    {
        for (int j = 1; j <= fieldsGridsSize; j++)
        {
            for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
            {
                ptrArray[5][i][j][k] = new GridsPoints();
                ptrArray[5][i][j][k]->InttoPos3(5, i, j, k);
                ptrArray[5][i][j][k]->XYZtoB(ptrArray[5][i][j][k]->Pos3());
                ptrArray[5][i][j][k]->XYZtoVel();
                ptrArray[5][i][j][k]->XYZtoE();
                ptrArray[5][i][j][k]->XYZtoDensity( 0);
                ptrArray[5][i][j][k]->SetStopSign(0);
                ptrArray[5][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 3
    // share len with face 1
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[3][1][j][k] = ptrArray[1][fieldsGridsSize + 1][j][k];
        }
    }
    // share len with face 2
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[3][i][fieldsGridsSize + 1][k] = ptrArray[2][fieldsGridsSize + 2 - i][fieldsGridsSize + 1][k];
        }
    }
    // share len with face 4
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[3][fieldsGridsSize + 1][j][k] = ptrArray[4][1][j][k];
        }
    }

    // share len with face 5
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[3][i][1][k] = ptrArray[5][fieldsGridsSize + 2 - i][1][k];
        }
    }
    for (int i = 2; i <= fieldsGridsSize; i++)
    {
        for (int j = 2; j <= fieldsGridsSize; j++)
        {
            for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
            {
                ptrArray[3][i][j][k] = new GridsPoints();
                ptrArray[3][i][j][k]->InttoPos3(3, i, j, k);
                ptrArray[3][i][j][k]->XYZtoB(ptrArray[3][i][j][k]->Pos3());
                ptrArray[3][i][j][k]->XYZtoVel();
                ptrArray[3][i][j][k]->XYZtoE();
                ptrArray[3][i][j][k]->XYZtoDensity( 0);
                ptrArray[3][i][j][k]->SetStopSign(0);
                ptrArray[3][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // info lens for adjant face
    //face 0
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[0][i][0][k] = ptrArray[5][i][fieldsGridsSize][k];     // bot
            ptrArray[0][fieldsGridsSize + 2][i][k] = ptrArray[1][2][i][k]; // right
            ptrArray[0][i][fieldsGridsSize + 2][k] = ptrArray[2][i][2][k]; // top
            ptrArray[0][0][i][k] = ptrArray[4][fieldsGridsSize][i][k];     // left
        }
    }
    //face 1
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[1][i][0][k] = ptrArray[5][fieldsGridsSize][fieldsGridsSize + 2 - i][k]; // bot
            ptrArray[1][fieldsGridsSize + 2][i][k] = ptrArray[3][2][i][k];                   // right
            ptrArray[1][i][fieldsGridsSize + 2][k] = ptrArray[2][fieldsGridsSize][i][k];     // top
            ptrArray[1][0][i][k] = ptrArray[0][fieldsGridsSize][i][k];                       // left
        }
    }
    //face 2
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[2][i][0][k] = ptrArray[0][i][fieldsGridsSize][k];                                         // bot
            ptrArray[2][fieldsGridsSize + 2][i][k] = ptrArray[1][i][fieldsGridsSize][k];                       // right
            ptrArray[2][i][fieldsGridsSize + 2][k] = ptrArray[3][fieldsGridsSize + 2 - i][fieldsGridsSize][k]; // top
            ptrArray[2][0][i][k] = ptrArray[4][fieldsGridsSize + 2 - i][fieldsGridsSize][k];                   // left
        }
    }
    //face 3
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[3][i][0][k] = ptrArray[5][fieldsGridsSize + 2 - i][2][k];                                 // bot
            ptrArray[3][fieldsGridsSize + 2][i][k] = ptrArray[4][2][i][k];                                     // right
            ptrArray[3][i][fieldsGridsSize + 2][k] = ptrArray[2][fieldsGridsSize + 2 - i][fieldsGridsSize][k]; // top
            ptrArray[3][0][i][k] = ptrArray[1][fieldsGridsSize][i][k];                                         // left
        }
    }
    //face 4
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[4][i][0][k] = ptrArray[5][2][i][k];                                         // bot
            ptrArray[4][fieldsGridsSize + 2][i][k] = ptrArray[0][2][i][k];                       // right
            ptrArray[4][i][fieldsGridsSize + 2][k] = ptrArray[2][2][fieldsGridsSize + 2 - i][k]; // top
            ptrArray[4][0][i][k] = ptrArray[3][fieldsGridsSize][i][k];                           // left
        }
    }
    //face 5
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[5][i][0][k] = ptrArray[3][fieldsGridsSize + 2 - i][2][k];                   // bot
            ptrArray[5][fieldsGridsSize + 2][i][k] = ptrArray[1][fieldsGridsSize + 2 - i][2][k]; // right
            ptrArray[5][i][fieldsGridsSize + 2][k] = ptrArray[0][i][2][k];                       // top
            ptrArray[5][0][i][k] = ptrArray[4][i][2][k];                                         // left
        }
    }

    // Info for not used pointers at four corners of each face
    for (int face = 0; face < totalFace; face++)
    {
        for (int k = 0; k <= fieldsGridsSize * grid_domain; k++)
        {
            ptrArray[face][0][0][k] = new GridsPoints();
            ptrArray[face][0][fieldsGridsSize + 2][k] = new GridsPoints();
            ptrArray[face][fieldsGridsSize + 2][0][k] = new GridsPoints();
            ptrArray[face][fieldsGridsSize + 2][fieldsGridsSize + 2][k] = new GridsPoints();
        }
    }

    //////////// test of location for peak points

    cout << " fieldsGridsSize " << fieldsGridsSize << endl;
    //    cout << "sizeof " << sizeof(ptrArray[0])/sizeof(ptrArray[0][0][0][0]) << endl;

    /*    // face 0, bottom left
    cout <<"0 "<< ptrArray[0][1][1][0]<<" 5 "<< ptrArray[5][1][fieldsGridsSize+1][0] <<" 4 "<< ptrArray[4][fieldsGridsSize+1][1][0] << endl;
    // face 0, bottom right
    cout <<"0 "<< ptrArray[0][fieldsGridsSize+1][1][0]<<" 5 "<< ptrArray[5][fieldsGridsSize+1][fieldsGridsSize+1][0] <<" 1 "<< ptrArray[1][1][1][0] << endl;
    // face 0, top left
    cout <<"0 "<< ptrArray[0][1][fieldsGridsSize+1][0]<<" 2 "<< ptrArray[2][1][1][0] <<" 4 "<< ptrArray[4][fieldsGridsSize+1][fieldsGridsSize+1][0] << endl;
    // face 0, top right
    cout <<"0 "<< ptrArray[0][fieldsGridsSize+1][fieldsGridsSize+1][0]<<" 1 "<< ptrArray[1][1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][fieldsGridsSize+1][1][0] << endl;
    // face 3, bottom left
    cout <<"3 "<< ptrArray[3][1][1][0]<<" 1 "<< ptrArray[1][fieldsGridsSize+1][1][0] <<" 5 "<< ptrArray[5][fieldsGridsSize+1][1][0] << endl;
    // face 3, bottom right
    cout <<"3 "<< ptrArray[3][fieldsGridsSize+1][1][0]<<" 4 "<< ptrArray[4][1][1][0] <<" 5 "<< ptrArray[5][1][1][0] << endl;
    // face 3, top left
    cout <<"3 "<< ptrArray[3][1][fieldsGridsSize+1][0]<<" 1 "<< ptrArray[1][fieldsGridsSize+1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][fieldsGridsSize+1][fieldsGridsSize+1][0] << endl;
    // face 3, top right
    cout <<"3 "<< ptrArray[3][fieldsGridsSize+1][fieldsGridsSize+1][0]<<" 4 "<< ptrArray[4][1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][1][fieldsGridsSize+1][0] << endl;
    //face 0 
    for( int i =0; i < 6; i++)
    {   
        cout << "face " << i << endl;
        cout << "botleft " << ptrArray[i][0][1][0] << " " << ptrArray[i][1][0][0] << endl;
        cout << "botright" << ptrArray[i][fieldsGridsSize+2][1][0] << " " << ptrArray[i][fieldsGridsSize+1][0][0] << endl;
        cout << "topleft " << ptrArray[i][0][fieldsGridsSize+1][0] << " " << ptrArray[i][1][fieldsGridsSize+2][0] << endl;
        cout << "topright" << ptrArray[i][fieldsGridsSize+2][fieldsGridsSize+1][0] << " " << ptrArray[i][fieldsGridsSize+1][fieldsGridsSize+2][0] << endl << endl;        
    }
*/
    return ptrArray;
}

///////////////////////////////////////////////////////////////
// Create cells
GridsCells ****GridsCellsCreation(GridsCells *mem_GridsCells,
                                  GridsPoints *****ptrArray)
{
    // f, i, j, k, index of cell of ptrArrayCell
    GridsCells ****ptrArrayCells;
    ptrArrayCells = new GridsCells ***[totalFace];
    for (int f = 0; f < totalFace; f++)
    {
        ptrArrayCells[f] = new GridsCells **[fieldsGridsSize];
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            ptrArrayCells[f][i] = new GridsCells *[fieldsGridsSize];
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                //    ptrArrayCells[f][i][j] = new GridsCells[fieldsGridsSize];
                ptrArrayCells[f][i][j] = mem_GridsCells + f * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain +
                                         i * fieldsGridsSize * fieldsGridsSize * grid_domain +
                                         j * fieldsGridsSize * grid_domain;
                for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                {
                    ptrArrayCells[f][i][j][k] = GridsCells();
                    ptrArrayCells[f][i][j][k].InitialParticlesCells();
                    ptrArrayCells[f][i][j][k].InitialRandomIndexPar();
                    //ptrArrayCells[f][i][j][k].InitialParVelDist();
                    ptrArrayCells[f][i][j][k].InitialParVelDistArray();
                    // function is for index of ptrArrayGrids
                    ptrArrayCells[f][i][j][k].SetVolume(CellVolume(ptrArray, f, i + 1, j + 1, k));
                    //
                    ptrArrayCells[f][i][j][k].SetB3(CellB3(ptrArray, f, i + 1, j + 1, k));
                }
            }
        }
    }
    return ptrArrayCells;
}

// Create an array of pointer to store position of dual cell grid point(corners)
// [face][fieldsGridsSize+2][fieldsGridsSize+2][fieldsGridsSize]
Vector3 *****GridsCreation_dual(GridsPoints *****ptrArray)
{
    Vector3 *****ptrArray_dual;
    ptrArray_dual = new Vector3 ****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrArray_dual[face] = new Vector3 ***[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrArray_dual[face][i] = new Vector3 **[fieldsGridsSize + 2];
            for (int j = 0; j < fieldsGridsSize + 2; j++)
            {
                ptrArray_dual[face][i][j] = new Vector3 *[fieldsGridsSize * grid_domain];
            }
        }
    }

    // i, j, k are index of main cell, the same with grid points of dual cell
    // face 0 ( to us)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                ptrArray_dual[0][i][j][k] = new Vector3((ptrArray[0][i][j][k]->Pos3().x() +
                                                         ptrArray[0][i + 1][j][k]->Pos3().x() +
                                                         ptrArray[0][i][j + 1][k]->Pos3().x() +
                                                         ptrArray[0][i + 1][j + 1][k]->Pos3().x() +
                                                         ptrArray[0][i][j][k + 1]->Pos3().x() +
                                                         ptrArray[0][i + 1][j][k + 1]->Pos3().x() +
                                                         ptrArray[0][i][j + 1][k + 1]->Pos3().x() +
                                                         ptrArray[0][i + 1][j + 1][k + 1]->Pos3().x()) *
                                                            0.1250,
                                                        (ptrArray[0][i][j][k]->Pos3().y() +
                                                         ptrArray[0][i + 1][j][k]->Pos3().y() +
                                                         ptrArray[0][i][j + 1][k]->Pos3().y() +
                                                         ptrArray[0][i + 1][j + 1][k]->Pos3().y() +
                                                         ptrArray[0][i][j][k + 1]->Pos3().y() +
                                                         ptrArray[0][i + 1][j][k + 1]->Pos3().y() +
                                                         ptrArray[0][i][j + 1][k + 1]->Pos3().y() +
                                                         ptrArray[0][i + 1][j + 1][k + 1]->Pos3().y()) *
                                                            0.1250,
                                                        (ptrArray[0][i][j][k]->Pos3().z() +
                                                         ptrArray[0][i + 1][j][k]->Pos3().z() +
                                                         ptrArray[0][i][j + 1][k]->Pos3().z() +
                                                         ptrArray[0][i + 1][j + 1][k]->Pos3().z() +
                                                         ptrArray[0][i][j][k + 1]->Pos3().z() +
                                                         ptrArray[0][i + 1][j][k + 1]->Pos3().z() +
                                                         ptrArray[0][i][j + 1][k + 1]->Pos3().z() +
                                                         ptrArray[0][i + 1][j + 1][k + 1]->Pos3().z()) *
                                                            0.1250);
            }
        }
    }
    // face 1 ( on the right)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                ptrArray_dual[1][i][j][k] = new Vector3((ptrArray[1][i][j][k]->Pos3().x() +
                                                         ptrArray[1][i + 1][j][k]->Pos3().x() +
                                                         ptrArray[1][i][j + 1][k]->Pos3().x() +
                                                         ptrArray[1][i + 1][j + 1][k]->Pos3().x() +
                                                         ptrArray[1][i][j][k + 1]->Pos3().x() +
                                                         ptrArray[1][i + 1][j][k + 1]->Pos3().x() +
                                                         ptrArray[1][i][j + 1][k + 1]->Pos3().x() +
                                                         ptrArray[1][i + 1][j + 1][k + 1]->Pos3().x()) *
                                                            0.1250,
                                                        (ptrArray[1][i][j][k]->Pos3().y() +
                                                         ptrArray[1][i + 1][j][k]->Pos3().y() +
                                                         ptrArray[1][i][j + 1][k]->Pos3().y() +
                                                         ptrArray[1][i + 1][j + 1][k]->Pos3().y() +
                                                         ptrArray[1][i][j][k + 1]->Pos3().y() +
                                                         ptrArray[1][i + 1][j][k + 1]->Pos3().y() +
                                                         ptrArray[1][i][j + 1][k + 1]->Pos3().y() +
                                                         ptrArray[1][i + 1][j + 1][k + 1]->Pos3().y()) *
                                                            0.1250,
                                                        (ptrArray[1][i][j][k]->Pos3().z() +
                                                         ptrArray[1][i + 1][j][k]->Pos3().z() +
                                                         ptrArray[1][i][j + 1][k]->Pos3().z() +
                                                         ptrArray[1][i + 1][j + 1][k]->Pos3().z() +
                                                         ptrArray[1][i][j][k + 1]->Pos3().z() +
                                                         ptrArray[1][i + 1][j][k + 1]->Pos3().z() +
                                                         ptrArray[1][i][j + 1][k + 1]->Pos3().z() +
                                                         ptrArray[1][i + 1][j + 1][k + 1]->Pos3().z()) *
                                                            0.1250);
            }
        }
    }
    // face 2( on the top)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                ptrArray_dual[2][i][j][k] = new Vector3((ptrArray[2][i][j][k]->Pos3().x() +
                                                         ptrArray[2][i + 1][j][k]->Pos3().x() +
                                                         ptrArray[2][i][j + 1][k]->Pos3().x() +
                                                         ptrArray[2][i + 1][j + 1][k]->Pos3().x() +
                                                         ptrArray[2][i][j][k + 1]->Pos3().x() +
                                                         ptrArray[2][i + 1][j][k + 1]->Pos3().x() +
                                                         ptrArray[2][i][j + 1][k + 1]->Pos3().x() +
                                                         ptrArray[2][i + 1][j + 1][k + 1]->Pos3().x()) *
                                                            0.1250,
                                                        (ptrArray[2][i][j][k]->Pos3().y() +
                                                         ptrArray[2][i + 1][j][k]->Pos3().y() +
                                                         ptrArray[2][i][j + 1][k]->Pos3().y() +
                                                         ptrArray[2][i + 1][j + 1][k]->Pos3().y() +
                                                         ptrArray[2][i][j][k + 1]->Pos3().y() +
                                                         ptrArray[2][i + 1][j][k + 1]->Pos3().y() +
                                                         ptrArray[2][i][j + 1][k + 1]->Pos3().y() +
                                                         ptrArray[2][i + 1][j + 1][k + 1]->Pos3().y()) *
                                                            0.1250,
                                                        (ptrArray[2][i][j][k]->Pos3().z() +
                                                         ptrArray[2][i + 1][j][k]->Pos3().z() +
                                                         ptrArray[2][i][j + 1][k]->Pos3().z() +
                                                         ptrArray[2][i + 1][j + 1][k]->Pos3().z() +
                                                         ptrArray[2][i][j][k + 1]->Pos3().z() +
                                                         ptrArray[2][i + 1][j][k + 1]->Pos3().z() +
                                                         ptrArray[2][i][j + 1][k + 1]->Pos3().z() +
                                                         ptrArray[2][i + 1][j + 1][k + 1]->Pos3().z()) *
                                                            0.1250);
            }
        }
    }
    // face 4( on the left)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                ptrArray_dual[4][i][j][k] = new Vector3((ptrArray[4][i][j][k]->Pos3().x() +
                                                         ptrArray[4][i + 1][j][k]->Pos3().x() +
                                                         ptrArray[4][i][j + 1][k]->Pos3().x() +
                                                         ptrArray[4][i + 1][j + 1][k]->Pos3().x() +
                                                         ptrArray[4][i][j][k + 1]->Pos3().x() +
                                                         ptrArray[4][i + 1][j][k + 1]->Pos3().x() +
                                                         ptrArray[4][i][j + 1][k + 1]->Pos3().x() +
                                                         ptrArray[4][i + 1][j + 1][k + 1]->Pos3().x()) *
                                                            0.1250,
                                                        (ptrArray[4][i][j][k]->Pos3().y() +
                                                         ptrArray[4][i + 1][j][k]->Pos3().y() +
                                                         ptrArray[4][i][j + 1][k]->Pos3().y() +
                                                         ptrArray[4][i + 1][j + 1][k]->Pos3().y() +
                                                         ptrArray[4][i][j][k + 1]->Pos3().y() +
                                                         ptrArray[4][i + 1][j][k + 1]->Pos3().y() +
                                                         ptrArray[4][i][j + 1][k + 1]->Pos3().y() +
                                                         ptrArray[4][i + 1][j + 1][k + 1]->Pos3().y()) *
                                                            0.1250,
                                                        (ptrArray[4][i][j][k]->Pos3().z() +
                                                         ptrArray[4][i + 1][j][k]->Pos3().z() +
                                                         ptrArray[4][i][j + 1][k]->Pos3().z() +
                                                         ptrArray[4][i + 1][j + 1][k]->Pos3().z() +
                                                         ptrArray[4][i][j][k + 1]->Pos3().z() +
                                                         ptrArray[4][i + 1][j][k + 1]->Pos3().z() +
                                                         ptrArray[4][i][j + 1][k + 1]->Pos3().z() +
                                                         ptrArray[4][i + 1][j + 1][k + 1]->Pos3().z()) *
                                                            0.1250);
            }
        }
    }
    // face 5( on the bottom)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                ptrArray_dual[5][i][j][k] = new Vector3((ptrArray[5][i][j][k]->Pos3().x() +
                                                         ptrArray[5][i + 1][j][k]->Pos3().x() +
                                                         ptrArray[5][i][j + 1][k]->Pos3().x() +
                                                         ptrArray[5][i + 1][j + 1][k]->Pos3().x() +
                                                         ptrArray[5][i][j][k + 1]->Pos3().x() +
                                                         ptrArray[5][i + 1][j][k + 1]->Pos3().x() +
                                                         ptrArray[5][i][j + 1][k + 1]->Pos3().x() +
                                                         ptrArray[5][i + 1][j + 1][k + 1]->Pos3().x()) *
                                                            0.1250,
                                                        (ptrArray[5][i][j][k]->Pos3().y() +
                                                         ptrArray[5][i + 1][j][k]->Pos3().y() +
                                                         ptrArray[5][i][j + 1][k]->Pos3().y() +
                                                         ptrArray[5][i + 1][j + 1][k]->Pos3().y() +
                                                         ptrArray[5][i][j][k + 1]->Pos3().y() +
                                                         ptrArray[5][i + 1][j][k + 1]->Pos3().y() +
                                                         ptrArray[5][i][j + 1][k + 1]->Pos3().y() +
                                                         ptrArray[5][i + 1][j + 1][k + 1]->Pos3().y()) *
                                                            0.1250,
                                                        (ptrArray[5][i][j][k]->Pos3().z() +
                                                         ptrArray[5][i + 1][j][k]->Pos3().z() +
                                                         ptrArray[5][i][j + 1][k]->Pos3().z() +
                                                         ptrArray[5][i + 1][j + 1][k]->Pos3().z() +
                                                         ptrArray[5][i][j][k + 1]->Pos3().z() +
                                                         ptrArray[5][i + 1][j][k + 1]->Pos3().z() +
                                                         ptrArray[5][i][j + 1][k + 1]->Pos3().z() +
                                                         ptrArray[5][i + 1][j + 1][k + 1]->Pos3().z()) *
                                                            0.1250);
            }
        }
    }
    // face 3 ( on the back)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                ptrArray_dual[3][i][j][k] = new Vector3((ptrArray[3][i][j][k]->Pos3().x() +
                                                         ptrArray[3][i + 1][j][k]->Pos3().x() +
                                                         ptrArray[3][i][j + 1][k]->Pos3().x() +
                                                         ptrArray[3][i + 1][j + 1][k]->Pos3().x() +
                                                         ptrArray[3][i][j][k + 1]->Pos3().x() +
                                                         ptrArray[3][i + 1][j][k + 1]->Pos3().x() +
                                                         ptrArray[3][i][j + 1][k + 1]->Pos3().x() +
                                                         ptrArray[3][i + 1][j + 1][k + 1]->Pos3().x()) *
                                                            0.1250,
                                                        (ptrArray[3][i][j][k]->Pos3().y() +
                                                         ptrArray[3][i + 1][j][k]->Pos3().y() +
                                                         ptrArray[3][i][j + 1][k]->Pos3().y() +
                                                         ptrArray[3][i + 1][j + 1][k]->Pos3().y() +
                                                         ptrArray[3][i][j][k + 1]->Pos3().y() +
                                                         ptrArray[3][i + 1][j][k + 1]->Pos3().y() +
                                                         ptrArray[3][i][j + 1][k + 1]->Pos3().y() +
                                                         ptrArray[3][i + 1][j + 1][k + 1]->Pos3().y()) *
                                                            0.1250,
                                                        (ptrArray[3][i][j][k]->Pos3().z() +
                                                         ptrArray[3][i + 1][j][k]->Pos3().z() +
                                                         ptrArray[3][i][j + 1][k]->Pos3().z() +
                                                         ptrArray[3][i + 1][j + 1][k]->Pos3().z() +
                                                         ptrArray[3][i][j][k + 1]->Pos3().z() +
                                                         ptrArray[3][i + 1][j][k + 1]->Pos3().z() +
                                                         ptrArray[3][i][j + 1][k + 1]->Pos3().z() +
                                                         ptrArray[3][i + 1][j + 1][k + 1]->Pos3().z()) *
                                                            0.1250);
            }
        }
    }

    // info lens for adjacent face
    // face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[0][i][0][k] = ptrArray_dual[5][i][fieldsGridsSize][k];     // bot
            ptrArray_dual[0][fieldsGridsSize + 1][i][k] = ptrArray_dual[1][1][i][k]; // right
            ptrArray_dual[0][i][fieldsGridsSize + 1][k] = ptrArray_dual[2][i][1][k]; // top
            ptrArray_dual[0][0][i][k] = ptrArray_dual[4][fieldsGridsSize][i][k];     // left
        }
    }

    //face 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[1][i][0][k] = ptrArray_dual[5][fieldsGridsSize][fieldsGridsSize + 1 - i][k]; // bot
            ptrArray_dual[1][fieldsGridsSize + 1][i][k] = ptrArray_dual[3][1][i][k];                   // right
            ptrArray_dual[1][i][fieldsGridsSize + 1][k] = ptrArray_dual[2][fieldsGridsSize][i][k];     // top
            ptrArray_dual[1][0][i][k] = ptrArray_dual[0][fieldsGridsSize][i][k];                       // left
        }
    }

    //face 2
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[2][i][0][k] = ptrArray_dual[0][i][fieldsGridsSize][k];                                         // bot
            ptrArray_dual[2][fieldsGridsSize + 1][i][k] = ptrArray_dual[1][i][fieldsGridsSize][k];                       // right
            ptrArray_dual[2][i][fieldsGridsSize + 1][k] = ptrArray_dual[3][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // top
            ptrArray_dual[2][0][i][k] = ptrArray_dual[4][fieldsGridsSize + 1 - i][fieldsGridsSize][k];                   // left
        }
    }

    //face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[3][i][0][k] = ptrArray_dual[5][fieldsGridsSize + 1 - i][1][k];                                 // bot
            ptrArray_dual[3][fieldsGridsSize + 1][i][k] = ptrArray_dual[4][1][i][k];                                     // right
            ptrArray_dual[3][i][fieldsGridsSize + 1][k] = ptrArray_dual[2][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // top
            ptrArray_dual[3][0][i][k] = ptrArray_dual[1][fieldsGridsSize][i][k];                                         // left
        }
    }

    //face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[4][i][0][k] = ptrArray_dual[5][1][i][k];                                         // bot
            ptrArray_dual[4][fieldsGridsSize + 1][i][k] = ptrArray_dual[0][1][i][k];                       // right
            ptrArray_dual[4][i][fieldsGridsSize + 1][k] = ptrArray_dual[2][1][fieldsGridsSize + 1 - i][k]; // top
            ptrArray_dual[4][0][i][k] = ptrArray_dual[3][fieldsGridsSize][i][k];                           // left
        }
    }

    //face 5
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[5][i][0][k] = ptrArray_dual[3][fieldsGridsSize + 1 - i][1][k];                   // bot
            ptrArray_dual[5][fieldsGridsSize + 1][i][k] = ptrArray_dual[1][fieldsGridsSize + 1 - i][1][k]; // right
            ptrArray_dual[5][i][fieldsGridsSize + 1][k] = ptrArray_dual[0][i][1][k];                       // top
            ptrArray_dual[5][0][i][k] = ptrArray_dual[4][i][1][k];                                         // left
        }
    }
    // not used pointfor( int face = 0; face < totalFace; face++)

    for (int face = 0; face < totalFace; face++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrArray_dual[face][0][0][k] = new Vector3(0.0, 0.0, 0.0);
            ptrArray_dual[face][0][fieldsGridsSize + 1][k] = new Vector3(0.0, 0.0, 0.0);
            ptrArray_dual[face][fieldsGridsSize + 1][0][k] = new Vector3(0.0, 0.0, 0.0);
            ptrArray_dual[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k] = new Vector3(0.0, 0.0, 0.0);
        }
    }

    return ptrArray_dual;
}

GridsPoints *****GridsCreation(GridsPoints *****ptrArray, int gridsSize)
{
    //    GridsPoints **a = new GridsPoints* [totalFace * (fieldsGridsSize+1) * (fieldsGridsSize+1) * (fieldsGridsSize+1)];
    //    GridsPoints *ptrArray[totalFace][fieldsGridsSize+1][fieldsGridsSize+1][radialGridsSize+1];

    ptrArray = new GridsPoints ****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrArray[face] = new GridsPoints ***[fieldsGridsSize + 3];
        for (int i = 0; i <= fieldsGridsSize + 2; i++)
        {
            ptrArray[face][i] = new GridsPoints **[fieldsGridsSize + 3];
            for (int j = 0; j <= fieldsGridsSize + 2; j++)
            {
                ptrArray[face][i][j] = new GridsPoints *[gridsSize + 1];
            }
        }
    }
    // face 0 (to us)
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j <= fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k <= gridsSize; k++)
            {
                ptrArray[0][i][j][k] = new GridsPoints();
                ptrArray[0][i][j][k]->InttoPos3(0, i, j, k);
                ptrArray[0][i][j][k]->XYZtoB(ptrArray[0][i][j][k]->Pos3());
                ptrArray[0][i][j][k]->XYZtoVel();
                ptrArray[0][i][j][k]->XYZtoE();
                ptrArray[0][i][j][k]->XYZtoDensity( 0);
                ptrArray[0][i][j][k]->SetStopSign(0);
                ptrArray[0][i][j][k]->SetTemperature(0.0);
    /*            if( k == fieldsGridsSize) std::cout << ptrArray[0][i][j][k]->Vel3().x() << " "
                    << ptrArray[0][i][j][k]->Vel3().y() << " "
                    << ptrArray[0][i][j][k]->Vel3().z() << std::endl;
    */      }
        }
    }

    // face 1 (on the right)
    // share len with face 0
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[1][1][j][k] = ptrArray[0][fieldsGridsSize + 1][j][k];
        }
    }

    for (int i = 2; i <= fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j <= fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k <= gridsSize; k++)
            {
                ptrArray[1][i][j][k] = new GridsPoints();
                ptrArray[1][i][j][k]->InttoPos3(1, i, j, k);
                ptrArray[1][i][j][k]->XYZtoB(ptrArray[1][i][j][k]->Pos3());
                ptrArray[1][i][j][k]->XYZtoVel();
                ptrArray[1][i][j][k]->XYZtoE();
                ptrArray[1][i][j][k]->XYZtoDensity( 0);
                ptrArray[1][i][j][k]->SetStopSign(0);
                ptrArray[1][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 2 (on the top)
    // share len with face 0
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[2][i][1][k] = ptrArray[0][i][fieldsGridsSize + 1][k];
        }
    }
    // share len with face 1
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[2][fieldsGridsSize + 1][j][k] = ptrArray[1][j][fieldsGridsSize + 1][k];
        }
    }
    for (int i = 1; i <= fieldsGridsSize; i++)
    {
        for (int j = 2; j <= fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k <= gridsSize; k++)
            {
                ptrArray[2][i][j][k] = new GridsPoints();
                ptrArray[2][i][j][k]->InttoPos3(2, i, j, k);
                ptrArray[2][i][j][k]->XYZtoB(ptrArray[2][i][j][k]->Pos3());
                ptrArray[2][i][j][k]->XYZtoVel();
                ptrArray[2][i][j][k]->XYZtoE();
                ptrArray[2][i][j][k]->XYZtoDensity( 0);
                ptrArray[2][i][j][k]->SetStopSign(0);
                ptrArray[2][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // face 4 (on the left)
    // share len with face 0
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[4][fieldsGridsSize + 1][j][k] = ptrArray[0][1][j][k];
        }
    }
    // share len with face 2
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[4][i][fieldsGridsSize + 1][k] = ptrArray[2][1][fieldsGridsSize + 2 - i][k];
        }
    }

    for (int i = 1; i <= fieldsGridsSize; i++)
    {
        for (int j = 1; j <= fieldsGridsSize; j++)
        {
            for (int k = 0; k <= gridsSize; k++)
            {
                ptrArray[4][i][j][k] = new GridsPoints();
                ptrArray[4][i][j][k]->InttoPos3(4, i, j, k);
                ptrArray[4][i][j][k]->XYZtoB(ptrArray[4][i][j][k]->Pos3());
                ptrArray[4][i][j][k]->XYZtoVel();
                ptrArray[4][i][j][k]->XYZtoE();
                ptrArray[4][i][j][k]->XYZtoDensity( 0);
                ptrArray[4][i][j][k]->SetStopSign(0);
                ptrArray[4][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 5 (on the bottom)
    // share len with face 0
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[5][i][fieldsGridsSize + 1][k] = ptrArray[0][i][1][k];
        }
    }
    // share len with face 1
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[5][fieldsGridsSize + 1][j][k] = ptrArray[1][fieldsGridsSize + 2 - j][1][k];
        }
    }
    // share len with face 4
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[5][1][j][k] = ptrArray[4][j][1][k];
        }
    }
    for (int i = 2; i <= fieldsGridsSize; i++)
    {
        for (int j = 1; j <= fieldsGridsSize; j++)
        {
            for (int k = 0; k <= gridsSize; k++)
            {
                ptrArray[5][i][j][k] = new GridsPoints();
                ptrArray[5][i][j][k]->InttoPos3(5, i, j, k);
                ptrArray[5][i][j][k]->XYZtoB(ptrArray[5][i][j][k]->Pos3());
                ptrArray[5][i][j][k]->XYZtoVel();
                ptrArray[5][i][j][k]->XYZtoE();
                ptrArray[5][i][j][k]->XYZtoDensity( 0);
                ptrArray[5][i][j][k]->SetStopSign(0);
                ptrArray[5][i][j][k]->SetTemperature(0.0);
            }
        }
    }
    // face 3
    // share len with face 1
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[3][1][j][k] = ptrArray[1][fieldsGridsSize + 1][j][k];
        }
    }
    // share len with face 2
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[3][i][fieldsGridsSize + 1][k] = ptrArray[2][fieldsGridsSize + 2 - i][fieldsGridsSize + 1][k];
        }
    }
    // share len with face 4
    for (int j = 1; j <= fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[3][fieldsGridsSize + 1][j][k] = ptrArray[4][1][j][k];
        }
    }

    // share len with face 5
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[3][i][1][k] = ptrArray[5][fieldsGridsSize + 2 - i][1][k];
        }
    }
    for (int i = 2; i <= fieldsGridsSize; i++)
    {
        for (int j = 2; j <= fieldsGridsSize; j++)
        {
            for (int k = 0; k <= gridsSize; k++)
            {
                ptrArray[3][i][j][k] = new GridsPoints();
                ptrArray[3][i][j][k]->InttoPos3(3, i, j, k);
                ptrArray[3][i][j][k]->XYZtoB(ptrArray[3][i][j][k]->Pos3());
                ptrArray[3][i][j][k]->XYZtoVel();
                ptrArray[3][i][j][k]->XYZtoE();
                ptrArray[3][i][j][k]->XYZtoDensity( 0);
                ptrArray[3][i][j][k]->SetStopSign(0);
                ptrArray[3][i][j][k]->SetTemperature(0.0);
            }
        }
    }

    // info lens for adjant face
    //face 0
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[0][i][0][k] = ptrArray[5][i][fieldsGridsSize][k];     // bot
            ptrArray[0][fieldsGridsSize + 2][i][k] = ptrArray[1][2][i][k]; // right
            ptrArray[0][i][fieldsGridsSize + 2][k] = ptrArray[2][i][2][k]; // top
            ptrArray[0][0][i][k] = ptrArray[4][fieldsGridsSize][i][k];     // left
        }
    }
    //face 1
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[1][i][0][k] = ptrArray[5][fieldsGridsSize][fieldsGridsSize + 2 - i][k]; // bot
            ptrArray[1][fieldsGridsSize + 2][i][k] = ptrArray[3][2][i][k];                   // right
            ptrArray[1][i][fieldsGridsSize + 2][k] = ptrArray[2][fieldsGridsSize][i][k];     // top
            ptrArray[1][0][i][k] = ptrArray[0][fieldsGridsSize][i][k];                       // left
        }
    }
    //face 2
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[2][i][0][k] = ptrArray[0][i][fieldsGridsSize][k];                                         // bot
            ptrArray[2][fieldsGridsSize + 2][i][k] = ptrArray[1][i][fieldsGridsSize][k];                       // right
            ptrArray[2][i][fieldsGridsSize + 2][k] = ptrArray[3][fieldsGridsSize + 2 - i][fieldsGridsSize][k]; // top
            ptrArray[2][0][i][k] = ptrArray[4][fieldsGridsSize + 2 - i][fieldsGridsSize][k];                   // left
        }
    }
    //face 3
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[3][i][0][k] = ptrArray[5][fieldsGridsSize + 2 - i][2][k];                                 // bot
            ptrArray[3][fieldsGridsSize + 2][i][k] = ptrArray[4][2][i][k];                                     // right
            ptrArray[3][i][fieldsGridsSize + 2][k] = ptrArray[2][fieldsGridsSize + 2 - i][fieldsGridsSize][k]; // top
            ptrArray[3][0][i][k] = ptrArray[1][fieldsGridsSize][i][k];                                         // left
        }
    }
    //face 4
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[4][i][0][k] = ptrArray[5][2][i][k];                                         // bot
            ptrArray[4][fieldsGridsSize + 2][i][k] = ptrArray[0][2][i][k];                       // right
            ptrArray[4][i][fieldsGridsSize + 2][k] = ptrArray[2][2][fieldsGridsSize + 2 - i][k]; // top
            ptrArray[4][0][i][k] = ptrArray[3][fieldsGridsSize][i][k];                           // left
        }
    }
    //face 5
    for (int i = 1; i <= fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[5][i][0][k] = ptrArray[3][fieldsGridsSize + 2 - i][2][k];                   // bot
            ptrArray[5][fieldsGridsSize + 2][i][k] = ptrArray[1][fieldsGridsSize + 2 - i][2][k]; // right
            ptrArray[5][i][fieldsGridsSize + 2][k] = ptrArray[0][i][2][k];                       // top
            ptrArray[5][0][i][k] = ptrArray[4][i][2][k];                                         // left
        }
    }

    // Info for not used pointers at four corners of each face
    for (int face = 0; face < totalFace; face++)
    {
        for (int k = 0; k <= gridsSize; k++)
        {
            ptrArray[face][0][0][k] = new GridsPoints();
            ptrArray[face][0][fieldsGridsSize + 2][k] = new GridsPoints();
            ptrArray[face][fieldsGridsSize + 2][0][k] = new GridsPoints();
            ptrArray[face][fieldsGridsSize + 2][fieldsGridsSize + 2][k] = new GridsPoints();
        }
    }

    //////////// test of location for peak points

    cout << "fieldsGridsSize " << fieldsGridsSize << endl;
    //    cout << "sizeof " << sizeof(ptrArray[0])/sizeof(ptrArray[0][0][0][0]) << endl;

    /*    // face 0, bottom left
    cout <<"0 "<< ptrArray[0][1][1][0]<<" 5 "<< ptrArray[5][1][fieldsGridsSize+1][0] <<" 4 "<< ptrArray[4][fieldsGridsSize+1][1][0] << endl;
    // face 0, bottom right
    cout <<"0 "<< ptrArray[0][fieldsGridsSize+1][1][0]<<" 5 "<< ptrArray[5][fieldsGridsSize+1][fieldsGridsSize+1][0] <<" 1 "<< ptrArray[1][1][1][0] << endl;
    // face 0, top left
    cout <<"0 "<< ptrArray[0][1][fieldsGridsSize+1][0]<<" 2 "<< ptrArray[2][1][1][0] <<" 4 "<< ptrArray[4][fieldsGridsSize+1][fieldsGridsSize+1][0] << endl;
    // face 0, top right
    cout <<"0 "<< ptrArray[0][fieldsGridsSize+1][fieldsGridsSize+1][0]<<" 1 "<< ptrArray[1][1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][fieldsGridsSize+1][1][0] << endl;
    // face 3, bottom left
    cout <<"3 "<< ptrArray[3][1][1][0]<<" 1 "<< ptrArray[1][fieldsGridsSize+1][1][0] <<" 5 "<< ptrArray[5][fieldsGridsSize+1][1][0] << endl;
    // face 3, bottom right
    cout <<"3 "<< ptrArray[3][fieldsGridsSize+1][1][0]<<" 4 "<< ptrArray[4][1][1][0] <<" 5 "<< ptrArray[5][1][1][0] << endl;
    // face 3, top left
    cout <<"3 "<< ptrArray[3][1][fieldsGridsSize+1][0]<<" 1 "<< ptrArray[1][fieldsGridsSize+1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][fieldsGridsSize+1][fieldsGridsSize+1][0] << endl;
    // face 3, top right
    cout <<"3 "<< ptrArray[3][fieldsGridsSize+1][fieldsGridsSize+1][0]<<" 4 "<< ptrArray[4][1][fieldsGridsSize+1][0] <<" 2 "<< ptrArray[2][1][fieldsGridsSize+1][0] << endl;
    //face 0 
    for( int i =0; i < 6; i++)
    {   
        cout << "face " << i << endl;
        cout << "botleft " << ptrArray[i][0][1][0] << " " << ptrArray[i][1][0][0] << endl;
        cout << "botright" << ptrArray[i][fieldsGridsSize+2][1][0] << " " << ptrArray[i][fieldsGridsSize+1][0][0] << endl;
        cout << "topleft " << ptrArray[i][0][fieldsGridsSize+1][0] << " " << ptrArray[i][1][fieldsGridsSize+2][0] << endl;
        cout << "topright" << ptrArray[i][fieldsGridsSize+2][fieldsGridsSize+1][0] << " " << ptrArray[i][fieldsGridsSize+1][fieldsGridsSize+2][0] << endl << endl;        
    }
*/
    return ptrArray;
}

//************************************************************************
// Initialize top and bot temp grids
//************************************************************************

void InitializeTempGrids(GridsPoints *****ptrArray, GridsPoints *****ptrArray_bot, GridsPoints *****ptrArray_top, int gridsSize)
{
    int face, i, j, k;
    int k_top;
    for (face = 0; face < totalFace; face++)
    {
        for (i = 0; i < fieldsGridsSize; i++)
        {
            for (j = 0; j < fieldsGridsSize; j++)
            {
                for (k = 0; k < gridsSize + 1; k++)
                {
                    k_top = fieldsGridsSize - gridsSize + k;
                    ptrArray_bot[face][i][j][k]->CopyGridsPoints(*ptrArray[face][i][j][k]);
                    ptrArray_top[face][i][j][k]->CopyGridsPoints(*ptrArray[face][i][j][k_top]);
                }
            }
        }
    }
}

//************************************************************************
// temp value
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

void VectorCellField(Vector3 ***&cellArray)
{
    static Vector3 *mem_VectorCellField = new Vector3[(fieldsGridsSize + 2) * (fieldsGridsSize + 2) * fieldsGridsSize];
    cellArray = new Vector3 **[fieldsGridsSize + 2];
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        cellArray[i] = new Vector3 *[fieldsGridsSize + 2];
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            //            cellArray[i][j] = new Vector3[fieldsGridsSize];
            cellArray[i][j] = mem_VectorCellField + i * (fieldsGridsSize + 2) * (fieldsGridsSize) + j * fieldsGridsSize;
            for (int k = 0; k < fieldsGridsSize; k++)
            {
                cellArray[i][j][k] = Vector3(0.0, 0.0, 0.0);
            }
        }
    }
}

void VectorCellField_Vel(Vector3 ***&cellArray)
{
    static Vector3 *mem_VectorCellField_Vel = new Vector3[(fieldsGridsSize + 2) * (fieldsGridsSize + 2) * fieldsGridsSize];
    cellArray = new Vector3 **[fieldsGridsSize + 2];
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        cellArray[i] = new Vector3 *[fieldsGridsSize + 2];
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            cellArray[i][j] = mem_VectorCellField_Vel + i * (fieldsGridsSize + 2) * (fieldsGridsSize) + j * fieldsGridsSize;
            for (int k = 0; k < fieldsGridsSize; k++)
            {
                cellArray[i][j][k] = Vector3(0.0, 0.0, 0.0);
            }
        }
    }
}

void VectorCellField_Grad(Vector3 ***&cellArray)
{
    static Vector3 *mem_VectorCellField_Grad = new Vector3[(fieldsGridsSize + 2) * (fieldsGridsSize + 2) * fieldsGridsSize];
    cellArray = new Vector3 **[fieldsGridsSize + 2];
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        cellArray[i] = new Vector3 *[fieldsGridsSize + 2];
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            cellArray[i][j] = mem_VectorCellField_Grad + i * (fieldsGridsSize + 2) * (fieldsGridsSize) + j * fieldsGridsSize;
            for (int k = 0; k < fieldsGridsSize; k++)
            {
                cellArray[i][j][k] = Vector3(0.0, 0.0, 0.0);
            }
        }
    }
}

// DELET the NEW items
void DEL_VectorCellField(Vector3 ***&cellArray)
{
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            delete[] cellArray[i][j];
        }
        //   delete[] cellArray[i];
    }
    //delete[] cellArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Value the matrix field using finite volume method, put in the pointer
// of the MatrixField, value it, and return the pointer.
// Notice that the cell at corners should be absent in calculation.
Vector3 ***ValueCurlField(Vector3 ***curlArray_in, double ***ptrVolumeCellArray_in, GridsPoints *****ptrArray_in, int face_in, char field_in)
{
    //#pragma omp parallel for collapse(3)
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            for (int k = 0; k < fieldsGridsSize; k++)
            {
                if ((i == 0 && j == 0) ||
                    (i == 0 && j == fieldsGridsSize + 1) ||
                    (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1) ||
                    (i == fieldsGridsSize + 1 && j == 0))
                    continue;

                /*        Vector3 test = AreaVectorL( ptrArray_in, face_in, i, j, k);
                std::cout << test.x() << " " << test.y() << " " << test.z() << " " << test.norm();
                int pause;
                std::cin >> pause;
        */
                // for each cell, calculate sum(n X B( on face)) and devided by
                // Volume to get the curl B at the center of cell
                Vector3 temp = AreaVectorL(ptrArray_in, face_in, i, j, k).CrossProduct(FaceFieldVectorL(ptrArray_in, face_in, i, j, k, field_in));
                temp = temp.PlusProduct(
                    AreaVectorR(ptrArray_in, face_in, i, j, k).CrossProduct(FaceFieldVectorR(ptrArray_in, face_in, i, j, k, field_in)));
                temp = temp.PlusProduct(
                    AreaVectorT(ptrArray_in, face_in, i, j, k).CrossProduct(FaceFieldVectorT(ptrArray_in, face_in, i, j, k, field_in)));
                temp = temp.PlusProduct(
                    AreaVectorBot(ptrArray_in, face_in, i, j, k).CrossProduct(FaceFieldVectorBot(ptrArray_in, face_in, i, j, k, field_in)));
                temp = temp.PlusProduct(
                    AreaVectorF(ptrArray_in, face_in, i, j, k).CrossProduct(FaceFieldVectorF(ptrArray_in, face_in, i, j, k, field_in)));
                temp = temp.PlusProduct(
                    AreaVectorBack(ptrArray_in, face_in, i, j, k).CrossProduct(FaceFieldVectorBack(ptrArray_in, face_in, i, j, k, field_in)));
                double volumetemp = ptrVolumeCellArray_in[i][j][k];

                temp = temp.ScaleProduct(1.0 / volumetemp);
                curlArray_in[i][j][k].SetVector3(temp);
/*
                std::cout << face_in << field_in << i << j << k << " " << temp.norm() << std::endl;;
                int pause;
                if( temp.norm() !=0 ){
     //           std::cin >> pause;
                    }
 */           }
        }
    }
    return curlArray_in;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Value gradient field of Pe.
// gradientArray_in is in size of ( fsize+2 * fsize+2 * fsize) with vector3
// ptrVolumeCellArray is in size of ( fsize+2 * fsize+2 * fsize) with double
// Pe = n k T, in which n is the number density, k is the boltzmann constant, and T is the Te
Vector3 *****ValueGradient(GridsPoints *****ptrArray_in, Vector3 *****gradientArray_in, double ***ptrVolumeCellArray_in, char char_in)
{
#pragma omp parallel for collapse(4)
    // i, j, k are index of cells
    for (int face_in = 0; face_in < 6; face_in++)
    {
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 0; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                {
                    if ((i == 0 && j == 0) ||
                        (i == 0 && j == fieldsGridsSize + 1) ||
                        (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1) ||
                        (i == fieldsGridsSize + 1 && j == 0))
                        continue;

                    Vector3 temp = Vector3(0.0, 0.0, 0.0);

                    if (char_in == 'P') // electron pressure Pe
                    {
                        // for each cell, calculate sum of n(face vector) * densities and devided by
                        // Volume to get the gradient at the center of cell
                        temp = AreaVectorL(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNumberDensityL(ptrArray_in, face_in, i, j, k) * FaceTemperatureL(ptrArray_in, face_in, i, j, k) * boltzmann_k);
                        temp = temp.PlusProduct(
                            AreaVectorR(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNumberDensityR(ptrArray_in, face_in, i, j, k) * FaceTemperatureR(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                        temp = temp.PlusProduct(
                            AreaVectorT(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNumberDensityT(ptrArray_in, face_in, i, j, k) * FaceTemperatureT(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                        temp = temp.PlusProduct(
                            AreaVectorBot(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNumberDensityBot(ptrArray_in, face_in, i, j, k) * FaceTemperatureBot(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                        temp = temp.PlusProduct(
                            AreaVectorF(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNumberDensityF(ptrArray_in, face_in, i, j, k) * FaceTemperatureF(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                        temp = temp.PlusProduct(
                            AreaVectorBack(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNumberDensityBack(ptrArray_in, face_in, i, j, k) * FaceTemperatureBack(ptrArray_in, face_in, i, j, k) * boltzmann_k));
                        double volumetemp = ptrVolumeCellArray_in[i][j][k];
                        temp = temp.ScaleProduct(1.0 / volumetemp);

                        /*              std::cout << i << " " << j << " " << k << " " 
                <<" test gradPe " << temp.x() << " " << temp.y()<< " " <<temp.z() << std::endl << std::endl;

                if( temp.x() ==0)
                { 
                    std::cout << volumetemp << " " ;
                    std::cout << i << " " << j << " " << k << std::endl;
                    std::cout << " FaceNumberDensity " << FaceNumberDensityL(ptrArray_in, face_in, i, j, k) << " "
                    << FaceNumberDensityR(ptrArray_in, face_in, i, j, k) << " "
                    << FaceNumberDensityT(ptrArray_in, face_in, i, j, k) << " "
                    << FaceNumberDensityBot(ptrArray_in, face_in, i, j, k) << " "
                    << FaceNumberDensityF(ptrArray_in, face_in, i, j, k) << " "
                    << FaceNumberDensityBack(ptrArray_in, face_in, i, j, k) << std::endl;
                    std::cout << " --> " << temp.x() << " " << temp.y() << " " << temp.z() << std::endl;
                    int pause;
                    std::cin >> pause;
                }
            */
#pragma omp critical
                        gradientArray_in[face_in][i][j][k]->SetVector3(temp);
                    }
                    else if (char_in == 'B') // gradient |B|
                    {
                        temp = AreaVectorL(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNormBL(ptrArray_in, face_in, i, j, k));
                        temp = temp.PlusProduct(
                            AreaVectorR(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNormBR(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorT(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNormBT(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorBot(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNormBBot(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorF(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNormBF(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorBack(ptrArray_in, face_in, i, j, k).ScaleProduct(FaceNormBBack(ptrArray_in, face_in, i, j, k)));
                        double volumetemp = ptrVolumeCellArray_in[i][j][k];
                        temp = temp.ScaleProduct(1.0 / volumetemp);
#pragma omp critical
                        gradientArray_in[face_in][i][j][k]->SetVector3(temp);
                    }
                    else if (char_in == 'D')
                    {
                        temp = AreaVectorL(ptrArray_in, face_in, i, j, k).ScaleProduct(FacePotentialL(ptrArray_in, face_in, i, j, k));
                        temp = temp.PlusProduct(
                            AreaVectorR(ptrArray_in, face_in, i, j, k).ScaleProduct(FacePotentialR(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorT(ptrArray_in, face_in, i, j, k).ScaleProduct(FacePotentialT(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorBot(ptrArray_in, face_in, i, j, k).ScaleProduct(FacePotentialBot(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorF(ptrArray_in, face_in, i, j, k).ScaleProduct(FacePotentialF(ptrArray_in, face_in, i, j, k)));
                        temp = temp.PlusProduct(
                            AreaVectorBack(ptrArray_in, face_in, i, j, k).ScaleProduct(FacePotentialBack(ptrArray_in, face_in, i, j, k)));
                        double volumetemp = ptrVolumeCellArray_in[i][j][k];
                        temp = temp.ScaleProduct(1.0 / volumetemp);
#pragma omp critical
                        gradientArray_in[face_in][i][j][k]->SetVector3(temp);
                    }
                }
            }
        }
    }
    return gradientArray_in;
}

//************************************************************************
//************************************************************************
// FUNCTION
// UpdateVe3
void UpdateVe3(Vector3 ***curlField_in, GridsPoints *****ptrArray_in, int face_in)
{
    //#pragma omp parallel for collapse(3)
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 2; j++)
        {
            for (int k = 1; k < fieldsGridsSize; k++)
            {
                Vector3 curl_field = Vector3(0.0, 0.0, 0.0);

                if (i == 1 && j == 1)
                {
                    // gradPe at gridspoints
                    curl_field = curlField_in[1][0][k - 1].PlusProduct(
                        curlField_in[1][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[0][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][0][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][1][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[0][1][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == 1 && j == fieldsGridsSize + 1)
                {

                    curl_field = curlField_in[1][fieldsGridsSize + 1][k - 1].PlusProduct(
                        curlField_in[1][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[0][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][fieldsGridsSize + 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][fieldsGridsSize][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[0][fieldsGridsSize][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                {
                    curl_field = curlField_in[fieldsGridsSize][fieldsGridsSize + 1][k - 1].PlusProduct(
                        curlField_in[fieldsGridsSize][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize + 1][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][fieldsGridsSize + 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][fieldsGridsSize][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[fieldsGridsSize + 1][fieldsGridsSize][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == fieldsGridsSize + 1 && j == 1)
                {
                    curl_field = curlField_in[fieldsGridsSize][0][k - 1].PlusProduct(
                        curlField_in[fieldsGridsSize][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize + 1][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][0][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][1][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[fieldsGridsSize + 1][1][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else
                {
                    curl_field = curlField_in[i - 1][j - 1][k - 1].PlusProduct(
                        curlField_in[i][j - 1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i - 1][j][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i][j][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i - 1][j - 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i][j - 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i - 1][j][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[i][j][k])
                                     .ScaleProduct(1.0 / 8.0);
                }
                ptrArray_in[face_in][i][j][k]->updateve3(curl_field);
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION NOT used but is a good oppinion
// As in the updating curlField and gradientPe array, some variables are
// repeating calculating, it is suitable to put them in one function.
// Therefore, we need three matrix of curlB, curlE, and gradientPe.
// Assume they are curlB, curlE and gradPe, respectively.
void updateCellMatrix(Vector3 ****curlB_in, Vector3 ****curlE_in,
                      Vector3 ****gradPe_in, GridsPoints *****ptrArray_in, int face_in)
{
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            for (int k = 0; k < fieldsGridsSize + 2; k++)
            {
                if ((i == 0 && j == 0) ||
                    (i == 0 && j == fieldsGridsSize + 1) ||
                    (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1) ||
                    (i == fieldsGridsSize + 1 && j == 0))
                    continue;
                // Volume
                double volumetemp = CellVolume(ptrArray_in, face_in, i, j, k);
                // AreaVectors
                Vector3 nL = AreaVectorL(ptrArray_in, face_in, i, j, k);
                Vector3 nR = AreaVectorR(ptrArray_in, face_in, i, j, k);
                Vector3 nT = AreaVectorT(ptrArray_in, face_in, i, j, k);
                Vector3 nBot = AreaVectorBot(ptrArray_in, face_in, i, j, k);
                Vector3 nF = AreaVectorF(ptrArray_in, face_in, i, j, k);
                Vector3 nBack = AreaVectorBack(ptrArray_in, face_in, i, j, k);

                // for each cell, calculate sum of n(vector) * densities and devided by
                // Volume to get the gradient at the center of cell

                Vector3 tempGrad = nL.ScaleProduct(
                    FaceDensityL(ptrArray_in, face_in, i, j, k));
                tempGrad = tempGrad.PlusProduct(
                    nR.ScaleProduct(
                        FaceDensityR(ptrArray_in, face_in, i, j, k)));
                tempGrad = tempGrad.PlusProduct(
                    nT.ScaleProduct(
                        FaceDensityT(ptrArray_in, face_in, i, j, k)));
                tempGrad = tempGrad.PlusProduct(
                    nBot.ScaleProduct(
                        FaceDensityBot(ptrArray_in, face_in, i, j, k)));
                tempGrad = tempGrad.PlusProduct(
                    nF.ScaleProduct(
                        FaceDensityF(ptrArray_in, face_in, i, j, k)));
                tempGrad = tempGrad.PlusProduct(
                    nBack.ScaleProduct(
                        FaceDensityBack(ptrArray_in, face_in, i, j, k)));

                tempGrad = tempGrad.ScaleProduct(1 / volumetemp);
                gradPe_in[i][j][k]->SetVector3(tempGrad);

                Vector3 tempCurlB = nL.CrossProduct(
                    FaceFieldVectorL(ptrArray_in, face_in, i, j, k, 'B'));
                tempCurlB = tempCurlB.PlusProduct(
                    nR.CrossProduct(
                        FaceFieldVectorR(ptrArray_in, face_in, i, j, k, 'B')));
                tempCurlB = tempCurlB.PlusProduct(
                    nT.CrossProduct(
                        FaceFieldVectorT(ptrArray_in, face_in, i, j, k, 'B')));
                tempCurlB = tempCurlB.PlusProduct(
                    nBot.CrossProduct(
                        FaceFieldVectorBot(ptrArray_in, face_in, i, j, k, 'B')));
                tempCurlB = tempCurlB.PlusProduct(
                    nF.CrossProduct(
                        FaceFieldVectorF(ptrArray_in, face_in, i, j, k, 'B')));
                tempCurlB = tempCurlB.PlusProduct(
                    nBack.CrossProduct(
                        FaceFieldVectorBack(ptrArray_in, face_in, i, j, k, 'B')));

                tempCurlB = tempCurlB.ScaleProduct(1 / volumetemp);
                curlB_in[i][j][k]->SetVector3(tempCurlB);

                Vector3 tempCurlE = nL.CrossProduct(
                    FaceFieldVectorL(ptrArray_in, face_in, i, j, k, 'E'));
                tempCurlE = tempCurlE.PlusProduct(
                    nR.CrossProduct(
                        FaceFieldVectorR(ptrArray_in, face_in, i, j, k, 'E')));
                tempCurlE = tempCurlE.PlusProduct(
                    nT.CrossProduct(
                        FaceFieldVectorT(ptrArray_in, face_in, i, j, k, 'E')));
                tempCurlE = tempCurlE.PlusProduct(
                    nBot.CrossProduct(
                        FaceFieldVectorBot(ptrArray_in, face_in, i, j, k, 'E')));
                tempCurlE = tempCurlE.PlusProduct(
                    nF.CrossProduct(
                        FaceFieldVectorF(ptrArray_in, face_in, i, j, k, 'E')));
                tempCurlE = tempCurlE.PlusProduct(
                    nBack.CrossProduct(
                        FaceFieldVectorBack(ptrArray_in, face_in, i, j, k, 'E')));

                tempCurlE = tempCurlE.ScaleProduct(1 / volumetemp);
                curlE_in[i][j][k]->SetVector3(tempCurlE);
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Update E at grids for ve = vi ( no current)
// Update E at grids for ve ( with current)
// Used in the initialization function
void UpdateE3(Vector3 ***gradPe_in, GridsPoints *****ptrArray_in, int face_in)
{

#pragma omp parallel for collapse(3)

    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 2; j++)
        {
            for (int k = 1; k < fieldsGridsSize; k++)
            {
                Vector3 tempGradPe = Vector3(0.0, 0.0, 0.0);

                if (i == 1 && j == 1)
                {
                    // gradPe at gridspoints
                    tempGradPe = gradPe_in[1][0][k - 1].PlusProduct(
                        gradPe_in[1][1][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[0][1][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[1][0][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[1][1][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                                               gradPe_in[0][1][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == 1 && j == fieldsGridsSize + 1)
                {

                    tempGradPe = gradPe_in[1][fieldsGridsSize + 1][k - 1].PlusProduct(
                        gradPe_in[1][fieldsGridsSize][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[0][fieldsGridsSize][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[1][fieldsGridsSize + 1][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[1][fieldsGridsSize][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                                               gradPe_in[0][fieldsGridsSize][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                {
                    tempGradPe = gradPe_in[fieldsGridsSize][fieldsGridsSize + 1][k - 1].PlusProduct(
                        gradPe_in[fieldsGridsSize][fieldsGridsSize][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[fieldsGridsSize + 1][fieldsGridsSize][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[fieldsGridsSize][fieldsGridsSize + 1][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[fieldsGridsSize][fieldsGridsSize][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                                               gradPe_in[fieldsGridsSize + 1][fieldsGridsSize][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == fieldsGridsSize + 1 && j == 1)
                {
                    tempGradPe = gradPe_in[fieldsGridsSize][0][k - 1].PlusProduct(
                        gradPe_in[fieldsGridsSize][1][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[fieldsGridsSize + 1][1][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[fieldsGridsSize][0][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[fieldsGridsSize][1][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                                               gradPe_in[fieldsGridsSize + 1][1][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else
                {
                    tempGradPe = gradPe_in[i - 1][j - 1][k - 1].PlusProduct(
                        gradPe_in[i][j - 1][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[i - 1][j][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[i][j][k - 1]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[i - 1][j - 1][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[i][j - 1][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                        gradPe_in[i - 1][j][k]);
                    tempGradPe = tempGradPe.PlusProduct(
                                               gradPe_in[i][j][k])
                                     .ScaleProduct(1.0 / 8.0);
                }
                // update E

                ptrArray_in[face_in][i][j][k]->updateE(tempGradPe);

                if (k == 1)
                    ptrArray_in[face_in][i][j][0]->updateE(tempGradPe);
                if (k == fieldsGridsSize - 1)
                    ptrArray_in[face_in][i][j][fieldsGridsSize]->updateE(tempGradPe);
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION
// UpdateB3 vased on faraday's law

void UpdateB3(Vector3 ***curlField_in, GridsPoints *****ptrArray_in, int face_in)
{

#pragma omp parallel for collapse(3)

    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 2; j++)
        {
            for (int k = 1; k < fieldsGridsSize; k++)
            {
                Vector3 curl_field = Vector3(0.0, 0.0, 0.0);

                if (i == 1 && j == 1)
                {
                    // gradPe at gridspoints
                    curl_field = curlField_in[1][0][k - 1].PlusProduct(
                        curlField_in[1][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[0][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][0][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][1][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[0][1][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == 1 && j == fieldsGridsSize + 1)
                {

                    curl_field = curlField_in[1][fieldsGridsSize + 1][k - 1].PlusProduct(
                        curlField_in[1][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[0][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][fieldsGridsSize + 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[1][fieldsGridsSize][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[0][fieldsGridsSize][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                {
                    curl_field = curlField_in[fieldsGridsSize][fieldsGridsSize + 1][k - 1].PlusProduct(
                        curlField_in[fieldsGridsSize][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize + 1][fieldsGridsSize][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][fieldsGridsSize + 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][fieldsGridsSize][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[fieldsGridsSize + 1][fieldsGridsSize][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else if (i == fieldsGridsSize + 1 && j == 1)
                {
                    curl_field = curlField_in[fieldsGridsSize][0][k - 1].PlusProduct(
                        curlField_in[fieldsGridsSize][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize + 1][1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][0][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[fieldsGridsSize][1][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[fieldsGridsSize + 1][1][k])
                                     .ScaleProduct(1.0 / 6.0);
                }
                else
                {
                    curl_field = curlField_in[i - 1][j - 1][k - 1].PlusProduct(
                        curlField_in[i][j - 1][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i - 1][j][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i][j][k - 1]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i - 1][j - 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i][j - 1][k]);
                    curl_field = curl_field.PlusProduct(
                        curlField_in[i - 1][j][k]);
                    curl_field = curl_field.PlusProduct(
                                               curlField_in[i][j][k])
                                     .ScaleProduct(1.0 / 8.0);
                }
                // update B
                ptrArray_in[face_in][i][j][k]->updatedB(curl_field);
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the grad|B| on the gridspoints
// Calculate grad(Pe) on the gridspoints
void UpdateGrad(Vector3 *****gradVector, GridsPoints *****ptrArray_in, char char_in)
{
// i, j, k are index of girds points
#pragma omp parallel for collapse(4)
    for (int face_in = 0; face_in < totalFace; face_in++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    Vector3 tempGradVector = Vector3(0.0, 0.0, 0.0);
                    if (i == 1 && j == 1)
                    {
                        // gradPe at gridspoints
                        tempGradVector = gradVector[face_in][1][0][k - 1]->PlusProduct(
                            *gradVector[face_in][1][1][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][0][1][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][1][0][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][1][1][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                                                           *gradVector[face_in][0][1][k])
                                             .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == 1 && j == fieldsGridsSize + 1)
                    {

                        tempGradVector = gradVector[face_in][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            *gradVector[face_in][1][fieldsGridsSize][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][0][fieldsGridsSize][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][1][fieldsGridsSize + 1][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][1][fieldsGridsSize][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                                                           *gradVector[face_in][0][fieldsGridsSize][k])
                                             .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                    {
                        tempGradVector = gradVector[face_in][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            *gradVector[face_in][fieldsGridsSize][fieldsGridsSize][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][fieldsGridsSize + 1][fieldsGridsSize][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][fieldsGridsSize][fieldsGridsSize + 1][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][fieldsGridsSize][fieldsGridsSize][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                                                           *gradVector[face_in][fieldsGridsSize + 1][fieldsGridsSize][k])
                                             .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == 1)
                    {
                        tempGradVector = gradVector[face_in][fieldsGridsSize][0][k - 1]->PlusProduct(
                            *gradVector[face_in][fieldsGridsSize][1][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][fieldsGridsSize + 1][1][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][fieldsGridsSize][0][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][fieldsGridsSize][1][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                                                           *gradVector[face_in][fieldsGridsSize + 1][1][k])
                                             .ScaleProduct(1.0 / 6.0);
                    }
                    else
                    {
                        tempGradVector = gradVector[face_in][i - 1][j - 1][k - 1]->PlusProduct(
                            *gradVector[face_in][i][j - 1][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][i - 1][j][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][i][j][k - 1]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][i - 1][j - 1][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][i][j - 1][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                            *gradVector[face_in][i - 1][j][k]);
                        tempGradVector = tempGradVector.PlusProduct(
                                                           *gradVector[face_in][i][j][k])
                                             .ScaleProduct(1.0 / 8.0);
                    }
                    // update gradB3
                    //#pragma omp critical
                    if (char_in == 'B')
                    {
                        ptrArray_in[face_in][i][j][k]->SetGradNormB(tempGradVector);
                        if (k == 1)
                            ptrArray_in[face_in][i][j][0]->SetGradNormB(tempGradVector);
                        if (k == fieldsGridsSize - 1)
                            ptrArray_in[face_in][i][j][fieldsGridsSize]->SetGradNormB(tempGradVector);
                        //  std::cout << ptrArray_in[face_in][i][j][k]->GradB3().x() << " " << ptrArray_in[face_in][i][j][k]->GradB3().y() << " " << ptrArray_in[face_in][i][j][k]->GradB3().z() << std::endl;
                    }
                    else if (char_in == 'P' || char_in == 'D')
                    {
                        ptrArray_in[face_in][i][j][k]->SetGradPe(tempGradVector); // if gradPotential, the value is also stored in variable "gradPe" of CLASS gridsPoints 
                        if (k == 1)
                            ptrArray_in[face_in][i][j][0]->SetGradPe(tempGradVector);
                        if (k == fieldsGridsSize - 1)
                            ptrArray_in[face_in][i][j][fieldsGridsSize]->SetGradPe(tempGradVector);
                    }
                    else
                    {
                        std::cout << " no correct char put in \n";
                        exit(1);
                    }
                }
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the temprature of ions on the gridspoints
//void UpdateTempIons(GridsPoints *****ptrArray, GridsCells ****ptrArrayCells)
//{
//    for (int face_in = 0; face_in < totalFace; face_in++)
//    {
//        for (int i = 1; i < fieldsGridsSize + 2; i++)
//        {
//            for (int j = 1; j < fieldsGridsSize + 2; j++)
//            {
//                for (int k = 1; k < fieldsGridsSize; k++)
//                {
//                    double tempTH = 0.0;
//                    double tempTHe = 0.0;
//                    double tempTO = 0.0;
//                    if (i == 1 && j == 1)
//                    {
//                        tempTH = (ptrArrayCells[face_in][1][0][k - 1].TempH() +
//                                  ptrArrayCells[face_in][1][1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][0][1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][1][0][k].TempH() +
//                                  ptrArrayCells[face_in][1][1][k].TempH() +
//                                  ptrArrayCells[face_in][0][1][k].TempH()) /
//                                 6.0;
//                        tempTHe = (ptrArrayCells[face_in][1][0][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][1][1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][0][1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][1][0][k].TempHe() +
//                                   ptrArrayCells[face_in][1][1][k].TempHe() +
//                                   ptrArrayCells[face_in][0][1][k].TempHe()) /
//                                  6.0;
//                        tempTO = (ptrArrayCells[face_in][1][0][k - 1].TempO() +
//                                  ptrArrayCells[face_in][1][1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][0][1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][1][0][k].TempO() +
//                                  ptrArrayCells[face_in][1][1][k].TempO() +
//                                  ptrArrayCells[face_in][0][1][k].TempO()) /
//                                 6.0;
//                    }
//                    else if (i == 1 && j == fieldsGridsSize + 1)
//                    {
//                        tempTH = (ptrArrayCells[face_in][1][fieldsGridsSize + 1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][1][fieldsGridsSize][k - 1].TempH() +
//                                  ptrArrayCells[face_in][0][fieldsGridsSize][k - 1].TempH() +
//                                  ptrArrayCells[face_in][1][fieldsGridsSize + 1][k].TempH() +
//                                  ptrArrayCells[face_in][1][fieldsGridsSize][k].TempH() +
//                                  ptrArrayCells[face_in][0][fieldsGridsSize][k].TempH()) /
//                                 6.0;
//                        tempTHe = (ptrArrayCells[face_in][1][fieldsGridsSize + 1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][1][fieldsGridsSize][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][0][fieldsGridsSize][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][1][fieldsGridsSize + 1][k].TempHe() +
//                                   ptrArrayCells[face_in][1][fieldsGridsSize][k].TempHe() +
//                                   ptrArrayCells[face_in][0][fieldsGridsSize][k].TempHe()) /
//                                  6.0;
//                        tempTO = (ptrArrayCells[face_in][1][fieldsGridsSize + 1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][1][fieldsGridsSize][k - 1].TempO() +
//                                  ptrArrayCells[face_in][0][fieldsGridsSize][k - 1].TempO() +
//                                  ptrArrayCells[face_in][1][fieldsGridsSize + 1][k].TempO() +
//                                  ptrArrayCells[face_in][1][fieldsGridsSize][k].TempO() +
//                                  ptrArrayCells[face_in][0][fieldsGridsSize][k].TempO()) /
//                                 6.0;
//                    }
//                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
//                    {
//                        tempTH = (ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize + 1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize][k - 1].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][fieldsGridsSize][k - 1].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize + 1][k].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize][k].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][fieldsGridsSize][k].TempH()) /
//                                 6.0;
//                        tempTHe = (ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize + 1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize + 1][fieldsGridsSize][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize + 1][k].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize][k].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize + 1][fieldsGridsSize][k].TempHe()) /
//                                  6.0;
//                        tempTO = (ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize + 1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize][k - 1].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][fieldsGridsSize][k - 1].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize + 1][k].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][fieldsGridsSize][k].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][fieldsGridsSize][k].TempO()) /
//                                 6.0;
//                    }
//                    else if (i == fieldsGridsSize + 1 && j == 1)
//                    {
//                        tempTH = (ptrArrayCells[face_in][fieldsGridsSize][0][k - 1].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][0][k].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][1][k].TempH() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][1][k].TempH()) /
//                                 6.0;
//                        tempTHe = (ptrArrayCells[face_in][fieldsGridsSize][0][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize][1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize + 1][1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize][0][k].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize][1][k].TempHe() +
//                                   ptrArrayCells[face_in][fieldsGridsSize + 1][1][k].TempHe()) /
//                                  6.0;
//                        tempTO = (ptrArrayCells[face_in][fieldsGridsSize][0][k - 1].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][0][k].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize][1][k].TempO() +
//                                  ptrArrayCells[face_in][fieldsGridsSize + 1][1][k].TempO()) /
//                                 6.0;
//                    }
//                    else
//                    {
//                        tempTH = (ptrArrayCells[face_in][i - 1][j - 1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][i][j - 1][k - 1].TempH() +
//                                  ptrArrayCells[face_in][i - 1][j][k - 1].TempH() +
//                                  ptrArrayCells[face_in][i][j][k - 1].TempH() +
//                                  ptrArrayCells[face_in][i - 1][j - 1][k].TempH() +
//                                  ptrArrayCells[face_in][i][j - 1][k].TempH() +
//                                  ptrArrayCells[face_in][i - 1][j][k].TempH() +
//                                  ptrArrayCells[face_in][i][j][k].TempH()) /
//                                 8.0;
//                        tempTHe = (ptrArrayCells[face_in][i - 1][j - 1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][i][j - 1][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][i - 1][j][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][i][j][k - 1].TempHe() +
//                                   ptrArrayCells[face_in][i - 1][j - 1][k].TempHe() +
//                                   ptrArrayCells[face_in][i][j - 1][k].TempHe() +
//                                   ptrArrayCells[face_in][i - 1][j][k].TempHe() +
//                                   ptrArrayCells[face_in][i][j][k].TempHe()) /
//                                  8.0;
//                        tempTO = (ptrArrayCells[face_in][i - 1][j - 1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][i][j - 1][k - 1].TempO() +
//                                  ptrArrayCells[face_in][i - 1][j][k - 1].TempO() +
//                                  ptrArrayCells[face_in][i][j][k - 1].TempO() +
//                                  ptrArrayCells[face_in][i - 1][j - 1][k].TempO() +
//                                  ptrArrayCells[face_in][i][j - 1][k].TempO() +
//                                  ptrArrayCells[face_in][i - 1][j][k].TempO() +
//                                  ptrArrayCells[face_in][i][j][k].TempO()) /
//                                 8.0;
//                    }
//                    // update T of ions
//                    ptrArray[face_in][i][j][k]->SetTemperature(tempTH, tempTHe, tempTO);
//                    if (k == 1)
//                        ptrArray[face_in][i][j][0]->SetTemperature(tempTH, tempTHe, tempTO);
//                    if (k == fieldsGridsSize - 1)
//                        ptrArray[face_in][i][j][fieldsGridsSize]->SetTemperature(tempTH, tempTHe, tempTO);
//                }
//            }
//        }
//    }
//}

//************************************************************************
//************************************************************************
// FUNCTION
// Each calculation are on the girds.
// Calculate the grad|B| on the gridspoints
void UpdateGradPe(Vector3 *****gradPe, GridsPoints *****ptrArray)
{

#pragma omp parallel for collapse(3)
    for (int face_in = 0; face_in < totalFace; face_in++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize; k++)
                {
                    Vector3 tempGradPe = Vector3(0.0, 0.0, 0.0);

                    if (i == 1 && j == 1)
                    {
                        // gradPe at gridspoints
                        tempGradPe = gradPe[face_in][1][0][k - 1]->PlusProduct(
                            *gradPe[face_in][1][1][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][0][1][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][1][0][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][1][1][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                                                   *gradPe[face_in][0][1][k])
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == 1 && j == fieldsGridsSize + 1)
                    {

                        tempGradPe = gradPe[face_in][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            *gradPe[face_in][1][fieldsGridsSize][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][0][fieldsGridsSize][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][1][fieldsGridsSize + 1][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][1][fieldsGridsSize][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                                                   *gradPe[face_in][0][fieldsGridsSize][k])
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                    {
                        tempGradPe = gradPe[face_in][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            *gradPe[face_in][fieldsGridsSize][fieldsGridsSize][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][fieldsGridsSize + 1][fieldsGridsSize][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][fieldsGridsSize][fieldsGridsSize + 1][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][fieldsGridsSize][fieldsGridsSize][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                                                   *gradPe[face_in][fieldsGridsSize + 1][fieldsGridsSize][k])
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == 1)
                    {
                        tempGradPe = gradPe[face_in][fieldsGridsSize][0][k - 1]->PlusProduct(
                            *gradPe[face_in][fieldsGridsSize][1][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][fieldsGridsSize + 1][1][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][fieldsGridsSize][0][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][fieldsGridsSize][1][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                                                   *gradPe[face_in][fieldsGridsSize + 1][1][k])
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else
                    {
                        tempGradPe = gradPe[face_in][i - 1][j - 1][k - 1]->PlusProduct(
                            *gradPe[face_in][i][j - 1][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][i - 1][j][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][i][j][k - 1]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][i - 1][j - 1][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][i][j - 1][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                            *gradPe[face_in][i - 1][j][k]);
                        tempGradPe = tempGradPe.PlusProduct(
                                                   *gradPe[face_in][i][j][k])
                                         .ScaleProduct(1.0 / 8.0);
                    }
                    // update gradB3
                    ptrArray[face_in][i][j][k]->SetGradPe(tempGradPe);

                    if (k == 1)
                        ptrArray[face_in][i][j][0]->SetGradPe(ptrArray[face_in][i][j][k]->GradPe());
                    if (k == fieldsGridsSize - 1)
                        ptrArray[face_in][i][j][fieldsGridsSize]->SetGradPe(ptrArray[face_in][i][j][k]->GradPe());
                }
            }
        }
    }
}

//************************************************************************
// Function
// each calculation are on the gridpoints of main cell
// Ve and E
void UpdateVeGridsMain(GridsPoints *****ptrArray, Vector3 ******ptrEFace_dual, int type)
{
    // i, j, k are the index of main cell/ grids points
    // i, j range are 0-fsize+1
#pragma omp parallel for collapse(4)
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    Vector3 tempCurl = Vector3(0.0, 0.0, 0.0);
                    Vector3 ve = Vector3(0.0, 0.0, 0.0);
                    Vector3 tempE = Vector3(0.0, 0.0, 0.0);
                    
                    if (update_type == 1)   // include perturbation and waves
                    {
                        if (i == 1 && j == 1)
                        {
                            tempCurl = ptrEFace_dual[face][i][j][k][0]->PlusProduct(
                                *ptrEFace_dual[face][i][j - 1][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k - 1][2]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][2]);
                            tempCurl = tempCurl.ScaleProduct(0.2);
                        }
                        else if (i == 1 && j == fieldsGridsSize + 1)
                        {
                            tempCurl = ptrEFace_dual[face][i][j][k][0]->PlusProduct(
                                *ptrEFace_dual[face][i][j - 1][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k - 1][2]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][2]);
                            tempCurl = tempCurl.ScaleProduct(0.2);
                        }
                        else if (i == fieldsGridsSize + 1 && j == 1)
                        {
                            tempCurl = ptrEFace_dual[face][i - 1][j][k][0]->PlusProduct(
                                *ptrEFace_dual[face][i][j - 1][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k - 1][2]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][2]);
                            tempCurl = tempCurl.ScaleProduct(0.2);
                        }
                        else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                        {
                            tempCurl = ptrEFace_dual[face][i - 1][j][k][0]->PlusProduct(
                                *ptrEFace_dual[face][i][j - 1][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k - 1][2]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][2]);
                            tempCurl = tempCurl.ScaleProduct(0.2);
                        }
                        else
                        {
                            tempCurl = ptrEFace_dual[face][i - 1][j][k][0]->PlusProduct(
                                *ptrEFace_dual[face][i][j][k][0]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j - 1][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][1]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k - 1][2]);
                            tempCurl = tempCurl.PlusProduct(
                                *ptrEFace_dual[face][i][j][k][2]);
                            tempCurl = tempCurl.ScaleProduct(1.0 / 6.0);
                        }
                        //    if( tempCurl.norm() != 0.0)
                        //    std::cout << tempCurl.norm() << "  ===>  ";
                        //    Vector3 v3temp;
                        //    double Bnorm2;
                        if (ptrArray[face][i][j][k]->Density() > 0.0)   // normal situation
                        {
                            ve = ptrArray[face][i][j][k]->Vel3().MinusProduct(tempCurl.ScaleProduct(1.0 / qi0 / mu0 / ptrArray[face][i][j][k]->Density()));
                            tempE = ve.CrossProduct(ptrArray[face][i][j][k]->B3()).PlusProduct(
                                        ptrArray[face][i][j][k]->GradPe().ScaleProduct(1.0/qi0/ptrArray[face][i][j][k]->Density())).ScaleProduct(-1.0);
                            // test for no veXB term
                            // tempE = ptrArray[face][i][j][k]->GradPe().ScaleProduct(1.0/qi0/ptrArray[face][i][j][k]->Density()).ScaleProduct(-1.0);
                            // test for no grad Pe term
                            // tempE = ve.CrossProduct(ptrArray[face][i][j][k]->B3()).ScaleProduct(-1.0);
                        }
                        else    // in case of ZERO density
                        {
                            ve = ptrArray[face][i][j][k]->Vel3();
                            tempE = ve.CrossProduct(ptrArray[face][i][j][k]->B3()).ScaleProduct(-1.0);
                        }
                        //    if( ve.norm() != 0.0 && k == 31)
                        //    {
                        //        std::cout << ve.norm() << "  ===>  " << i << " " << j << " " << ptrArray[face][i][j][k]->Vel_e3().norm() << " \n";
                        //    }
                        ptrArray[face][i][j][k]->SetVel_e3(ve);
                        //    if( ve.norm() != 0.0 && k == 31)
                        //    std::cout<<  ptrArray[face][i][j][k]->Vel_e3().norm() << " \n";
                        ptrArray[face][i][j][k]->SetdE3(tempE);
                    }
                    else if (update_type == 0)  // cases don't include waves
                    {
#pragma omp critical
                        if (type == 0) // for case of gradPe to calculate parallel E, stored in variable "E" of CLASS gridspoints
                        {
                            //  update ve in CLASS gridspoints
                            ve = ptrArray[face][i][j][k]->Vel3();
                            ptrArray[face][i][j][k]->SetVel_e3(ve);
                            // E
                            if(ptrArray[face][i][j][k]->Density() > 0.0)    // normal situation
                            {
                                // calculate E related to gradPe and eXB in Ohm's Law
                                tempE = ve.CrossProduct(ptrArray[face][i][j][k]->B3()).PlusProduct(
                                        ptrArray[face][i][j][k]->GradPe().ScaleProduct(1.0/qi0/ptrArray[face][i][j][k]->Density())).ScaleProduct(-1.0);                        
                                // test for the term related to gradPe in Ohm's Law
                                //tempE = ptrArray[face][i][j][k]->GradPe().ScaleProduct(-1.0 / qi0 / ptrArray[face][i][j][k]->Density());
                            } else      // in case ZERO density
                            {
                                tempE = ve.CrossProduct(ptrArray[face][i][j][k]->B3()).ScaleProduct(-1.0);
                            }
                            // restrict tempE no greater than 1.0e-7
                            if (tempE.norm() > 1.0e-7)
                                tempE = tempE.NormalizedVector().ScaleProduct(1.0e-7);
                            // update E in CLASS gridspoints
                            ptrArray[face][i][j][k]->SetE3(tempE);
                            //  set edge layer
                            if (k == 1)
                            {
                                // set bottom layer k = 0, assume a 1.0e-7 magnitude of E
                                tempE = ptrArray[face][i][j][0]->Pos3().NormalizedVector().ScaleProduct(1.0e-7);
                                ptrArray[face][i][j][0]->SetE3(tempE);
                                ptrArray[face][i][j][0]->SetVel_e3(ve);
                            }
                            if (k == fieldsGridsSize - 1)
                            {
                                // set top layer, k = fieldsGridsSize
                                ptrArray[face][i][j][fieldsGridsSize]->SetE3(tempE);
                                ptrArray[face][i][j][fieldsGridsSize]->SetVel_e3(ve);
                            }
                        }
                        else if (type == 1) // add component of (grad Potential) for E, perpendicular E
                        {
                            tempE = ptrArray[face][i][j][k]->GradPe().ScaleProduct(-1.0); // variable "gradPe" is now grad potential instead
                            ptrArray[face][i][j][k]->SetE3(tempE.PlusProduct(ptrArray[face][i][j][k]->E3()));
                            // set edge layer
                            if (k == 1)
                                ptrArray[face][i][j][0]->SetE3(tempE.PlusProduct(ptrArray[face][i][j][k]->E3()));
                            if (k == fieldsGridsSize - 1)
                                ptrArray[face][i][j][fieldsGridsSize]->SetE3(tempE.PlusProduct(ptrArray[face][i][j][k]->E3()));
                        }
                    }
                }
            }
        }
    }
}

//************************************************************************
//  FOr velDist plot, temperature_pp&&pr, v_pp&&pr
void PrintOutHdf5Cells(GridsCells ****ptrArrayCells, int timeline)
{
    using namespace H5;
    char filename_H[80];
    char filename_He[80];
    char filename_O[80];

    char filen[80], filename[151];

    sprintf(filen, "moments_veldist%d.h5", timeline);
    strncpy(filename, outpdir, 150);
    strncat(filename, filen, 150);

    sprintf(filename_H, "velH_%d", timeline);
    sprintf(filename_He, "velHe_%d", timeline);
    sprintf(filename_O, "velO_%d", timeline);
    //
    H5std_string FILE_NAME(filename);
    H5std_string DATASET_NAME_H(filename_H);
    H5std_string DATASET_NAME_HE(filename_He);
    H5std_string DATASET_NAME_O(filename_O);
    //
    const int RANK = 6;
    // assign contingent memory
    double *data_mem = new double[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * velDistRange_para * velDistRange_mu];
    double ******array_data = new double *****[totalFace];
    for (int f = 0; f < totalFace; f++)
    {
        array_data[f] = new double ****[fieldsGridsSize];
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            array_data[f][i] = new double ***[fieldsGridsSize];
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                array_data[f][i][j] = new double **[fieldsGridsSize];
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    array_data[f][i][j][k] = new double *[velDistRange_para];
                    for (int r = 0; r < velDistRange_para; r++)
                    {
                        array_data[f][i][j][k][r] = data_mem + f * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * velDistRange_para * velDistRange_mu + i * fieldsGridsSize * fieldsGridsSize * velDistRange_para * velDistRange_mu + j * fieldsGridsSize * velDistRange_para * velDistRange_mu + k * velDistRange_para * velDistRange_mu + r * velDistRange_mu;
                    }
                }
            }
        }
    }
    //
    // set up data space
    Exception::dontPrint();
    hsize_t dim[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize, (hsize_t)fieldsGridsSize, 
                     (hsize_t)fieldsGridsSize, (hsize_t)velDistRange_para, (hsize_t)velDistRange_mu};
    DataSpace space_Dist(RANK, dim);
    //
    H5File file_Dist(FILE_NAME, H5F_ACC_TRUNC);
    // value and print
    //    vector<vector<double>> *ptr_a;
    //    //
    //    // H
    //    for (int f = 0; f < totalFace; f++)
    //    {
    //        for (int i = 0; i < fieldsGridsSize; i++)
    //        {
    //            for (int j = 0; j < fieldsGridsSize; j++)
    //            {
    //                for (int k = 0; k < fieldsGridsSize; k++)
    //                {
    //                    ptr_a = ptrArrayCells[f][i][j][k].VelDist_H();
    //                    for (int r = 0; r < velDistRange_para; r++)
    //                    {
    //                        for (int s = 0; s < velDistRange_mu; s++)
    //                        {
    //                            array_data[f][i][j][k][r][s] = (*ptr_a)[r][s];
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    //
    //    DataSet dataset_Dist_H;
    //    dataset_Dist_H = file_Dist.createDataSet(DATASET_NAME_H, PredType::NATIVE_DOUBLE, space_Dist);
    //    dataset_Dist_H.write(array_data[0][0][0][0][0], PredType::NATIVE_DOUBLE);
    //    dataset_Dist_H.close();
    //    //
    //    // He
    //    for (int f = 0; f < totalFace; f++)
    //    {
    //        for (int i = 0; i < fieldsGridsSize; i++)
    //        {
    //            for (int j = 0; j < fieldsGridsSize; j++)
    //            {
    //                for (int k = 0; k < fieldsGridsSize; k++)
    //                {
    //                    ptr_a = ptrArrayCells[f][i][j][k].VelDist_He();
    //                    for (int r = 0; r < velDistRange_para; r++)
    //                    {
    //                        for (int s = 0; s < velDistRange_mu; s++)
    //                        {
    //                            array_data[f][i][j][k][r][s] = (*ptr_a)[r][s];
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    //
    //    DataSet dataset_Dist_HE;
    //    dataset_Dist_HE = file_Dist.createDataSet(DATASET_NAME_HE, PredType::NATIVE_DOUBLE, space_Dist);
    //    dataset_Dist_HE.write(array_data[0][0][0][0][0], PredType::NATIVE_DOUBLE);
    //    dataset_Dist_HE.close();
    //    //
    //    // O
    //    for (int f = 0; f < totalFace; f++)
    //    {
    //        for (int i = 0; i < fieldsGridsSize; i++)
    //        {
    //            for (int j = 0; j < fieldsGridsSize; j++)
    //            {
    //                for (int k = 0; k < fieldsGridsSize; k++)
    //                {
    //                    ptr_a = ptrArrayCells[f][i][j][k].VelDist_O();
    //                    for (int r = 0; r < velDistRange_para; r++)
    //                    {
    //                        for (int s = 0; s < velDistRange_mu; s++)
    //                        {
    //                            array_data[f][i][j][k][r][s] = (*ptr_a)[r][s];
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    //
    //    DataSet dataset_Dist_O;
    //    dataset_Dist_O = file_Dist.createDataSet(DATASET_NAME_O, PredType::NATIVE_DOUBLE, space_Dist);
    //    dataset_Dist_O.write(array_data[0][0][0][0][0], PredType::NATIVE_DOUBLE);
    //    dataset_Dist_O.close();
    //    //
    //    file_Dist.close();
    //    //
    //    for (int f = 0; f < totalFace; f++)
    //    {
    //        for (int i = 0; i < fieldsGridsSize; i++)
    //        {
    //            for (int j = 0; j < fieldsGridsSize; j++)
    //            {
    //                for (int k = 0; k < fieldsGridsSize; k++)
    //                {
    //                    delete[] array_data[f][i][j][k];
    //                }
    //                delete[] array_data[f][i][j];
    //            }
    //            delete[] array_data[f][i];
    //        }
    //        delete[] array_data[f];
    //    }
    //    delete[] array_data;
    //    //
    //    delete[] data_mem;
    // temperature_pp&&pr and v_pp&&pr
    char filename_pp_pr[80];
    sprintf(filename_pp_pr, "pp_pr_%d", timeline);
    //H5std_string FILE_NAME("./MyData/DataForPlot.h5");
    H5std_string DATASET_NAME_PP_PR(filename_pp_pr);
    //
    const int RANK_PP_PR = 5;
    // apply for continue memory
    double *data_mem_pp_pr = new double[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * 12];
    double *****array_data_pp_pr = new double ****[totalFace];
    for (int f = 0; f < totalFace; f++)
    {
        array_data_pp_pr[f] = new double ***[fieldsGridsSize];
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            array_data_pp_pr[f][i] = new double **[fieldsGridsSize];
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                array_data_pp_pr[f][i][j] = new double *[fieldsGridsSize];
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    array_data_pp_pr[f][i][j][k] = data_mem_pp_pr + f * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * 12 + i * fieldsGridsSize * fieldsGridsSize * 12 + j * fieldsGridsSize * 12 + k * 12;
                }
            }
        }
    }
    //
    // set up data space
    Exception::dontPrint();
    hsize_t dim_pp_pr[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize, (hsize_t)fieldsGridsSize, 
                           (hsize_t)fieldsGridsSize, 12};
    DataSpace space_pp_pr(RANK_PP_PR, dim_pp_pr);
    //
    H5File file_pp_pr(FILE_NAME, H5F_ACC_RDWR);
    // value and print
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    //array_data_pp_pr[f][i][j][k][0] = ptrArrayCells[f][i][j][k].Temp_H_pp();
                    //array_data_pp_pr[f][i][j][k][1] = ptrArrayCells[f][i][j][k].Temp_H_pr();
                    //array_data_pp_pr[f][i][j][k][2] = ptrArrayCells[f][i][j][k].Temp_He_pp();
                    //array_data_pp_pr[f][i][j][k][3] = ptrArrayCells[f][i][j][k].Temp_He_pr();
                    //array_data_pp_pr[f][i][j][k][4] = ptrArrayCells[f][i][j][k].Temp_O_pp();
                    //array_data_pp_pr[f][i][j][k][5] = ptrArrayCells[f][i][j][k].Temp_O_pr();
                    //array_data_pp_pr[f][i][j][k][6] = ptrArrayCells[f][i][j][k].Vel_H_pp();
                    //array_data_pp_pr[f][i][j][k][7] = ptrArrayCells[f][i][j][k].Vel_H_pr();
                    //array_data_pp_pr[f][i][j][k][8] = ptrArrayCells[f][i][j][k].Vel_He_pp();
                    //array_data_pp_pr[f][i][j][k][9] = ptrArrayCells[f][i][j][k].Vel_He_pr();
                    //array_data_pp_pr[f][i][j][k][10] = ptrArrayCells[f][i][j][k].Vel_O_pp();
                    //array_data_pp_pr[f][i][j][k][11] = ptrArrayCells[f][i][j][k].Vel_O_pr();
                    //
                    //if( k ==0)
                    //std::cout << array_data_pp_pr[f][i][j][k][0] << " " <<  ptrArrayCells[f][i][j][k].Temp_H_pp() <<"\n";
                }
            }
        }
    }
    //
    DataSet dataset_pp_pr;
    dataset_pp_pr = file_pp_pr.createDataSet(DATASET_NAME_PP_PR, PredType::NATIVE_DOUBLE, space_pp_pr);
    dataset_pp_pr.write(array_data_pp_pr[0][0][0][0], PredType::NATIVE_DOUBLE);
    dataset_pp_pr.close();
    //
    file_pp_pr.close();
    //
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                delete[] array_data_pp_pr[f][i][j];
            }
            delete[] array_data_pp_pr[f][i];
        }
        delete[] array_data_pp_pr[f];
    }
    delete[] array_data_pp_pr;
    delete[] data_mem_pp_pr;
    //
}

//************************************************************************
//  FOr velDist plot    // have some issue
//void PrintOutHdf5Cells(GridsCells ****ptrArrayCells, int timeline, int h5FileCheck)
//{
//    using namespace H5;
//    char filename[80];
//    sprintf(filename, "ArrayOfCellsVelDist_%d", timeline);
//    //
//    H5std_string FILE_NAME("./Data/DataForPlot.h5");
//    H5std_string DATASET_NAME(filename);
//    //
//    const int RANK = 6;
//    // apply for continue memory
//    double *data_mem = new double[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * 6 * velDistRange_para];
//    double ******array_data = new double *****[totalFace];
//    for (int f = 0; f < totalFace; f++)
//    {
//        array_data[f] = new double ****[fieldsGridsSize];
//        for (int i = 0; i < fieldsGridsSize; i++)
//        {
//            array_data[f][i] = new double ***[fieldsGridsSize];
//            for (int j = 0; j < fieldsGridsSize; j++)
//            {
//                array_data[f][i][j] = new double **[fieldsGridsSize];
//                for (int k = 0; k < fieldsGridsSize; k++)
//                {
//                    array_data[f][i][j][k] = new double *[fieldsGridsSize];
//                    for (int m = 0; m < 6; m++)
//                    {
//                        array_data[f][i][j][k][m] = data_mem + f * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * 6 + i * fieldsGridsSize * fieldsGridsSize * 6 + j * fieldsGridsSize * 6 + m * 6;
//                        // array_data[f][i][j][k][m][n] is 'double'
//                    }
//                }
//            }
//        }
//    }
//    // value it
//    vector<double> *ptr_a;
//    for (int f = 0; f < totalFace; f++)
//    {
//        for (int i = 0; i < fieldsGridsSize; i++)
//        {
//            for (int j = 0; j < fieldsGridsSize; j++)
//            {
//                for (int k = 0; k < fieldsGridsSize; k++)
//                {
//                    // parallel H
//                    for (int n = 0; n < velDistRange_para; n++)
//                    {
//                        ptr_a = ptrArrayCells[f][i][j][k].VelDist_para_H();
//                        array_data[f][i][j][k][0][n] = (*ptr_a)[n];
//                    }
//                    // parallel He
//                    for (int n = 0; n < velDistRange_para; n++)
//                    {
//                        ptr_a = ptrArrayCells[f][i][j][k].VelDist_para_He();
//                        array_data[f][i][j][k][1][n] = (*ptr_a)[n];
//                    }
//                    // parallel O
//                    for (int n = 0; n < velDistRange_para; n++)
//                    {
//                        ptr_a = ptrArrayCells[f][i][j][k].VelDist_para_O();
//                        array_data[f][i][j][k][2][n] = (*ptr_a)[n];
//                    }
//                    // perp H
//                    for (int n = 0; n < velDistRange_mu; n++)
//                    {
//                        ptr_a = ptrArrayCells[f][i][j][k].VelDist_perp_H();
//                        array_data[f][i][j][k][3][n] = (*ptr_a)[n];
//                    }
//                    // perp He
//                    for (int n = 0; n < velDistRange_mu; n++)
//                    {
//                        ptr_a = ptrArrayCells[f][i][j][k].VelDist_perp_He();
//                        array_data[f][i][j][k][4][n] = (*ptr_a)[n];
//                    }
//                    // perp O
//                    for (int n = 0; n < velDistRange_mu; n++)
//                    {
//                        ptr_a = ptrArrayCells[f][i][j][k].VelDist_perp_O();
//                        array_data[f][i][j][k][5][n] = (*ptr_a)[n];
//                    }
//                }
//            }
//        }
//    }
//    // printout
//    Exception::dontPrint();
//    hsize_t dim[] = {totalFace, fieldsGridsSize, fieldsGridsSize, fieldsGridsSize, 6, velDistRange_mu};
//    DataSpace space_Dist(RANK, dim);
//    //
//    if (h5FileCheck == 1)
//    {
//        //    H5File* file = new H5File( FILE_NAME, H5F_ACC_RDWR);
//        H5File file_Dist(FILE_NAME, H5F_ACC_RDWR);
//        DataSet dataset_Dist;
//        dataset_Dist = file_Dist.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, space_Dist);
//        dataset_Dist.write(array_data[0][0][0][0][0], PredType::NATIVE_DOUBLE);
//        //
//        space_Dist.close();
//        dataset_Dist.close();
//    }
//    else
//    {
//        std::cout << " File for printout did not exist. -- Function(PrintOutHdf5Cells) \n";
//        exit(1);
//    }
//    //
//    delete data_mem;
//}

//************************************************************************
//************************************************************************
// FUNCTION
// Printout the gridpoints on the girds as hdf5 format
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
void PrintOutHdf5_const(GridsPoints *****ptrArray_in)
{
    using namespace H5;
    char filename[151];

    strncpy(filename, outpdir, 150);
    strncat(filename, "/arrayconstants.h5", 149);

    H5std_string FILE_NAME(filename);
    H5std_string DATASET_CONST_NAME("ArrayOfGrids_const");
    H5std_string MEMBERx("x");
    H5std_string MEMBERy("y");
    H5std_string MEMBERz("z");
    H5std_string MEMBER_pos3("pos3");
    H5std_string MEMBER_e3("e3");
    H5std_string MEMBER_b3("b3");
    H5std_string MEMBER_dB3("dB3");
    H5std_string MEMBER_gradB3("gradB3");
    H5std_string MEMBER_ve3("ve3");
    H5std_string MEMBER_v3("v3");
    H5std_string MEMBER_vH3("vH3");
    H5std_string MEMBER_vHe3("vHe3");
    H5std_string MEMBER_vO3("vO3");
    H5std_string MEMBER_density_H("densityH");
    H5std_string MEMBER_density_He("densityHe");
    H5std_string MEMBER_density_O("densityO");
    H5std_string MEMBER_density("density");
    H5std_string MEMBER_te("temperature");
    H5std_string MEMBER_potential("potential");

    const int RANK = 4;

    typedef struct //Vector3_h5
    {
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    typedef struct //GridsPoints_const_h5
    {
        Vector3_h5 pos3;
        Vector3_h5 b3;
        double temperature;

        Vector3_h5 e3;
        Vector3_h5 v3;

        double density_H;
        double density_He;
        double density_O;

        double potential;

    } GridsPoints_const_h5; // dont varies with timeline

    GridsPoints_const_h5 *data_mem_const = new GridsPoints_const_h5[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain)];
    GridsPoints_const_h5 ****array_data_const = new GridsPoints_const_h5 ***[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        array_data_const[face] = new GridsPoints_const_h5 **[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            array_data_const[face][i] = new GridsPoints_const_h5 *[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                //               array_data_const[face][i][j] = new GridsPoints_const_h5[1+fieldsGridsSize];
                array_data_const[face][i][j] = data_mem_const + face * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) +
                                               i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) +
                                               j * (1 + fieldsGridsSize * grid_domain);
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    array_data_const[face][i][j][k].pos3 = {ptrArray_in[face][i + 1][j + 1][k]->Pos3().x(), ptrArray_in[face][i + 1][j + 1][k]->Pos3().y(), ptrArray_in[face][i + 1][j + 1][k]->Pos3().z()};
                    array_data_const[face][i][j][k].b3 = {ptrArray_in[face][i + 1][j + 1][k]->B3_base().x(), ptrArray_in[face][i + 1][j + 1][k]->B3_base().y(), ptrArray_in[face][i + 1][j + 1][k]->B3_base().z()};
                    array_data_const[face][i][j][k].temperature = ptrArray_in[face][i + 1][j + 1][k]->Temperature();
                    array_data_const[face][i][j][k].e3 = {ptrArray_in[face][i + 1][j + 1][k]->E3().x(), ptrArray_in[face][i + 1][j + 1][k]->E3().y(), ptrArray_in[face][i + 1][j + 1][k]->E3().z()};
                    array_data_const[face][i][j][k].v3 = {ptrArray_in[face][i + 1][j + 1][k]->Vel3().x(), ptrArray_in[face][i + 1][j + 1][k]->Vel3().y(), ptrArray_in[face][i + 1][j + 1][k]->Vel3().z()};
                    array_data_const[face][i][j][k].density_H = ptrArray_in[face][i + 1][j + 1][k]->Density_H();
                    array_data_const[face][i][j][k].density_He = ptrArray_in[face][i + 1][j + 1][k]->Density_He();
                    array_data_const[face][i][j][k].density_O = ptrArray_in[face][i + 1][j + 1][k]->Density_O();
                    array_data_const[face][i][j][k].potential = ptrArray_in[face][i + 1][j + 1][k]->Potential();
                }
            }
        }
    }
    Exception::dontPrint();
    hsize_t dim[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize + 1, (hsize_t)fieldsGridsSize + 1, 
                    (hsize_t)fieldsGridsSize * grid_domain + 1};
    DataSpace space(RANK, dim);

    // vector elements
    CompType mtype_vector3(sizeof(Vector3_h5));
    mtype_vector3.insertMember(MEMBERx, HOFFSET(Vector3_h5, v_x), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERy, HOFFSET(Vector3_h5, v_y), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERz, HOFFSET(Vector3_h5, v_z), PredType::NATIVE_DOUBLE);

    // const group variables
    CompType mtype_grids_const(sizeof(GridsPoints_const_h5));
    mtype_grids_const.insertMember(MEMBER_pos3, HOFFSET(GridsPoints_const_h5, pos3), mtype_vector3);
    mtype_grids_const.insertMember(MEMBER_b3, HOFFSET(GridsPoints_const_h5, b3), mtype_vector3);
    mtype_grids_const.insertMember(MEMBER_te, HOFFSET(GridsPoints_const_h5, temperature), PredType::NATIVE_DOUBLE);
    mtype_grids_const.insertMember(MEMBER_e3, HOFFSET(GridsPoints_const_h5, e3), mtype_vector3);
    mtype_grids_const.insertMember(MEMBER_v3, HOFFSET(GridsPoints_const_h5, v3), mtype_vector3);
    mtype_grids_const.insertMember(MEMBER_density_H, HOFFSET(GridsPoints_const_h5, density_H), PredType::NATIVE_DOUBLE);
    mtype_grids_const.insertMember(MEMBER_density_He, HOFFSET(GridsPoints_const_h5, density_He), PredType::NATIVE_DOUBLE);
    mtype_grids_const.insertMember(MEMBER_density_O, HOFFSET(GridsPoints_const_h5, density_O), PredType::NATIVE_DOUBLE);
    mtype_grids_const.insertMember(MEMBER_potential, HOFFSET(GridsPoints_const_h5, potential), PredType::NATIVE_DOUBLE);

    //H5File *file = new H5File(FILE_NAME, H5F_ACC_TRUNC);
    //// print const values at timeline = 0
    //DataSet *dataset_const;
    //dataset_const = new DataSet(file->createDataSet(DATASET_CONST_NAME, mtype_grids_const, space));
    //dataset_const->write(array_data_const[0][0][0], mtype_grids_const);
    //delete dataset_const;
    //delete file;
    //
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    DataSet dataset_const;
    dataset_const = file.createDataSet(DATASET_CONST_NAME, mtype_grids_const, space);
    dataset_const.write(array_data_const[0][0][0], mtype_grids_const);
    dataset_const.close();
    file.close();

    //
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            delete[] array_data_const[face][i];
        }
        delete[] array_data_const[face];
    }
    delete[] array_data_const;
    delete[] data_mem_const; // notice where is the definition
}

//************************************************************************
//************************************************************************
// FUNCTION
// Printout the gridpoints on the girds as hdf5 format
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
void PrintOutHdf5(GridsPoints *****ptrArray_in,
              int timeline)
{
    using namespace H5;
    char filename[151], filen[80];

    sprintf(filen, "/arrayvariables%d.h5", timeline);
    strncpy(filename, outpdir, 150);
    strncat(filename, filen, 150);

    H5std_string FILE_NAME(filename);
    sprintf(filen, "ArrayofGrids_%d", timeline);
    H5std_string DATASET_NAME(filen);
    H5std_string MEMBERx("x");
    H5std_string MEMBERy("y");
    H5std_string MEMBERz("z");
    H5std_string MEMBER_pos3("pos3");
    H5std_string MEMBER_e3("e3");
    H5std_string MEMBER_b3("b3");
    H5std_string MEMBER_dB3("dB3");
    H5std_string MEMBER_gradB3("gradB3");
    H5std_string MEMBER_ve3("ve3");
    H5std_string MEMBER_v3("v3");
    H5std_string MEMBER_vH3("vH3");
    H5std_string MEMBER_vHe3("vHe3");
    H5std_string MEMBER_vO3("vO3");
    H5std_string MEMBER_density_H("densityH");
    H5std_string MEMBER_density_He("densityHe");
    H5std_string MEMBER_density_O("densityO");
    H5std_string MEMBER_density("density");
    H5std_string MEMBER_te("temperature");
    H5std_string MEMBER_potential("potential");

    const int RANK = 4;

    typedef struct //Vector3_h5
    {
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    typedef struct //GridsPoints_h5
    {
        Vector3_h5 e3;
        Vector3_h5 dB3;
        Vector3_h5 gradB3;

        Vector3_h5 ve3;
        Vector3_h5 v3;
        Vector3_h5 vH3;
        Vector3_h5 vHe3;
        Vector3_h5 vO3;

        double density;
        double density_H;
        double density_He;
        double density_O;
    } GridsPoints_h5; // varies with timeline

    // Apply for continus memory
    GridsPoints_h5 *data_mem = new GridsPoints_h5[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain)];
    GridsPoints_h5 ****array_data = new GridsPoints_h5 ***[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        array_data[face] = new GridsPoints_h5 **[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            array_data[face][i] = new GridsPoints_h5 *[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                //           array_data[face][i][j] = new GridsPoints_h5[1+fieldsGridsSize];
                array_data[face][i][j] = data_mem + face * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) +
                                         i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) +
                                         j * (1 + fieldsGridsSize * grid_domain);
                //
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    array_data[face][i][j][k].e3 = {ptrArray_in[face][i + 1][j + 1][k]->E3().x(), ptrArray_in[face][i + 1][j + 1][k]->E3().y(), ptrArray_in[face][i + 1][j + 1][k]->E3().z()};
                    //
                    array_data[face][i][j][k].dB3 = {ptrArray_in[face][i + 1][j + 1][k]->DB3().x(), ptrArray_in[face][i + 1][j + 1][k]->DB3().y(), ptrArray_in[face][i + 1][j + 1][k]->DB3().z()};
                    array_data[face][i][j][k].gradB3 = {ptrArray_in[face][i + 1][j + 1][k]->GradB3().x(), ptrArray_in[face][i + 1][j + 1][k]->GradB3().y(), ptrArray_in[face][i + 1][j + 1][k]->GradB3().z()};
                    //
                    array_data[face][i][j][k].ve3 = {ptrArray_in[face][i + 1][j + 1][k]->Vel_e3().x(), ptrArray_in[face][i + 1][j + 1][k]->Vel_e3().y(), ptrArray_in[face][i + 1][j + 1][k]->Vel_e3().z()};
                    //
                    array_data[face][i][j][k].v3 = {ptrArray_in[face][i + 1][j + 1][k]->Vel3().x(), ptrArray_in[face][i + 1][j + 1][k]->Vel3().y(), ptrArray_in[face][i + 1][j + 1][k]->Vel3().z()};
                    array_data[face][i][j][k].vH3 = {ptrArray_in[face][i + 1][j + 1][k]->VelH3().x(), ptrArray_in[face][i + 1][j + 1][k]->VelH3().y(), ptrArray_in[face][i + 1][j + 1][k]->VelH3().z()};
                    array_data[face][i][j][k].vHe3 = {ptrArray_in[face][i + 1][j + 1][k]->VelHe3().x(), ptrArray_in[face][i + 1][j + 1][k]->VelHe3().y(), ptrArray_in[face][i + 1][j + 1][k]->VelHe3().z()};
                    array_data[face][i][j][k].vO3 = {ptrArray_in[face][i + 1][j + 1][k]->VelO3().x(), ptrArray_in[face][i + 1][j + 1][k]->VelO3().y(), ptrArray_in[face][i + 1][j + 1][k]->VelO3().z()};
                    //
                    array_data[face][i][j][k].density = ptrArray_in[face][i + 1][j + 1][k]->Density();
                    array_data[face][i][j][k].density_H = ptrArray_in[face][i + 1][j + 1][k]->Density_H();
                    array_data[face][i][j][k].density_He = ptrArray_in[face][i + 1][j + 1][k]->Density_He();
                    array_data[face][i][j][k].density_O = ptrArray_in[face][i + 1][j + 1][k]->Density_O();
                    //    std::cout << array_data[face][i][j][k].b3.v_x << std::endl;
                }
            }
        }
    }

    Exception::dontPrint();
    hsize_t dim[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize + 1, (hsize_t)fieldsGridsSize + 1, 
                    (hsize_t)fieldsGridsSize * grid_domain + 1};
    DataSpace space(RANK, dim);
    // vector elements
    CompType mtype_vector3(sizeof(Vector3_h5));
    mtype_vector3.insertMember(MEMBERx, HOFFSET(Vector3_h5, v_x), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERy, HOFFSET(Vector3_h5, v_y), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERz, HOFFSET(Vector3_h5, v_z), PredType::NATIVE_DOUBLE);

    // non-const group variables
    CompType mtype_grids(sizeof(GridsPoints_h5));
    mtype_grids.insertMember(MEMBER_e3, HOFFSET(GridsPoints_h5, e3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_dB3, HOFFSET(GridsPoints_h5, dB3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_gradB3, HOFFSET(GridsPoints_h5, gradB3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_ve3, HOFFSET(GridsPoints_h5, ve3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_v3, HOFFSET(GridsPoints_h5, v3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_vH3, HOFFSET(GridsPoints_h5, vH3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_vHe3, HOFFSET(GridsPoints_h5, vHe3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_vO3, HOFFSET(GridsPoints_h5, vO3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_density, HOFFSET(GridsPoints_h5, density), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_density_H, HOFFSET(GridsPoints_h5, density_H), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_density_He, HOFFSET(GridsPoints_h5, density_He), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_density_O, HOFFSET(GridsPoints_h5, density_O), PredType::NATIVE_DOUBLE);

    //H5File *file = new H5File(FILE_NAME, H5F_ACC_RDWR);
    //DataSet *dataset;
    //dataset = new DataSet(file->createDataSet(DATASET_NAME, mtype_grids, space));
    //dataset->write(array_data[0][0][0], mtype_grids);
    //// dataset->write(&array_data[0][0][0][0], mtype_grids);
    //delete dataset;
    //delete file;
    H5File file(FILE_NAME, H5F_ACC_TRUNC);
    DataSet dataset;
    dataset = file.createDataSet(DATASET_NAME, mtype_grids, space);
    dataset.write(array_data[0][0][0], mtype_grids);
    dataset.close();
    file.close();

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            delete[] array_data[face][i];
        }
        delete[] array_data[face];
    }
    delete[] array_data;
    //
    delete[] data_mem;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Printout the gridpoints on the girds as hdf5 format
// Step 1: Generate a new matrix fulling with gridspoints class
// Step 2: Print out it as .h5
// parts only
void PrintOutHdf5_Particles(int timeline,
                            vector<Particles> &ptrParticlesList_H,
                            vector<Particles> &ptrParticlesList_He,
                            vector<Particles> &ptrParticlesList_O)
{
    using namespace H5;

    char filen[80], filename[151];

    sprintf(filen, "/particles%d.h5", timeline);
    strncpy(filename, outpdir, 150);
    strncat(filename, filen, 150);

    H5std_string FILE_NAME(filename);
    H5std_string DATASET_NAME_H("Particles_H");
    H5std_string DATASET_NAME_He("Particles_He");
    H5std_string DATASET_NAME_O("Particles_O");

    H5std_string MEMBERx("x");
    H5std_string MEMBERy("y");
    H5std_string MEMBERz("z");

    H5std_string MEMBER_posUint("posUint");
    H5std_string MEMBER_pos3("pos3");
    H5std_string MEMBER_v3("v3");

    H5std_string MEMBER_weight("weight");
    H5std_string MEMBER_mu("mu");

    typedef struct //Vector3_h5
    {
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    typedef struct //Particles_h5
    {
        unsigned long long posUint;
        Vector3_h5 pos3;
        Vector3_h5 v3;
        double weight;
        double mu;
    } Particles_h5;

    // Apply for continus memory
    int numParticle_H = ptrParticlesList_H.size();
    int numParticle_He = ptrParticlesList_He.size();
    int numParticle_O = ptrParticlesList_O.size();

    Particles_h5 *array_data_particles_H = new Particles_h5[numParticle_H];
    Particles_h5 *array_data_particles_He = new Particles_h5[numParticle_He];
    Particles_h5 *array_data_particles_O = new Particles_h5[numParticle_O];

    int RANK = 2;
    hsize_t dimsH[2];
    hsize_t dimsHe[2];
    hsize_t dimsO[2];
    dimsH[0] = numParticle_H;
    dimsH[1] = 1;

    dimsHe[0] = numParticle_He;
    dimsHe[1] = 1;

    dimsO[0] = numParticle_O;
    dimsO[1] = 1;

    for (int i = 0; i < numParticle_H; i++)
    {
        array_data_particles_H[i].posUint = ptrParticlesList_H[i].PosUint();
        array_data_particles_H[i].pos3 = {ptrParticlesList_H[i].PosParticles().x(), ptrParticlesList_H[i].PosParticles().y(), ptrParticlesList_H[i].PosParticles().z()};
        array_data_particles_H[i].v3 = {ptrParticlesList_H[i].VelParticles().x(), ptrParticlesList_H[i].VelParticles().y(), ptrParticlesList_H[i].VelParticles().z()};
        array_data_particles_H[i].weight = ptrParticlesList_H[i].WeightNi();
        array_data_particles_H[i].mu = ptrParticlesList_H[i].MagneticIvarient();
    }
    for (int i = 0; i < numParticle_He; i++)
    {
        array_data_particles_He[i].posUint = ptrParticlesList_He[i].PosUint();
        array_data_particles_He[i].pos3 = {ptrParticlesList_He[i].PosParticles().x(), ptrParticlesList_He[i].PosParticles().y(), ptrParticlesList_He[i].PosParticles().z()};
        array_data_particles_He[i].v3 = {ptrParticlesList_He[i].VelParticles().x(), ptrParticlesList_He[i].VelParticles().y(), ptrParticlesList_He[i].VelParticles().z()};
        array_data_particles_He[i].weight = ptrParticlesList_He[i].WeightNi();
        array_data_particles_He[i].mu = ptrParticlesList_He[i].MagneticIvarient();
    }
    for (int i = 0; i < numParticle_O; i++)
    {
        array_data_particles_O[i].posUint = ptrParticlesList_O[i].PosUint();
        array_data_particles_O[i].pos3 = {ptrParticlesList_O[i].PosParticles().x(), ptrParticlesList_O[i].PosParticles().y(), ptrParticlesList_O[i].PosParticles().z()};
        array_data_particles_O[i].v3 = {ptrParticlesList_O[i].VelParticles().x(), ptrParticlesList_O[i].VelParticles().y(), ptrParticlesList_O[i].VelParticles().z()};
        array_data_particles_O[i].weight = ptrParticlesList_O[i].WeightNi();
        array_data_particles_O[i].mu = ptrParticlesList_O[i].MagneticIvarient();
    }

    // Vector elements
    CompType mtype_vector3(sizeof(Vector3_h5));
    mtype_vector3.insertMember(MEMBERx, HOFFSET(Vector3_h5, v_x), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERy, HOFFSET(Vector3_h5, v_y), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERz, HOFFSET(Vector3_h5, v_z), PredType::NATIVE_DOUBLE);

    CompType mtype_particles(sizeof(Particles_h5));
    mtype_particles.insertMember(MEMBER_posUint, HOFFSET(Particles_h5, posUint), PredType::NATIVE_ULLONG);
    mtype_particles.insertMember(MEMBER_pos3, HOFFSET(Particles_h5, pos3), mtype_vector3);
    mtype_particles.insertMember(MEMBER_v3, HOFFSET(Particles_h5, v3), mtype_vector3);
    mtype_particles.insertMember(MEMBER_weight, HOFFSET(Particles_h5, weight), PredType::NATIVE_DOUBLE);
    mtype_particles.insertMember(MEMBER_mu, HOFFSET(Particles_h5, mu), PredType::NATIVE_DOUBLE);

    cout << " PrintParticles: H " << numParticle_H << " He " << numParticle_He << " O " << numParticle_O << endl;
    H5File file_particles(FILE_NAME, H5F_ACC_TRUNC);

    DataSpace dataspace_particles_H = DataSpace(RANK, dimsH);
    DataSet dataset_H;
    dataset_H = file_particles.createDataSet(DATASET_NAME_H, mtype_particles, dataspace_particles_H);
    dataset_H.write(array_data_particles_H, mtype_particles);
    dataspace_particles_H.close();
    dataset_H.close();

    DataSpace dataspace_particles_He = DataSpace(RANK, dimsHe);
    DataSet dataset_He;
    dataset_He = file_particles.createDataSet(DATASET_NAME_He, mtype_particles, dataspace_particles_He);
    dataset_He.write(array_data_particles_He, mtype_particles);
    dataspace_particles_He.close();
    dataset_He.close();

    DataSpace dataspace_particles_O = DataSpace(RANK, dimsO);
    DataSet dataset_O;
    dataset_O = file_particles.createDataSet(DATASET_NAME_O, mtype_particles, dataspace_particles_O);
    dataset_O.write(array_data_particles_O, mtype_particles);
    dataspace_particles_O.close();
    dataset_O.close();

    file_particles.close();

    delete array_data_particles_H;
    delete array_data_particles_He;
    delete array_data_particles_O;
}
// parts only
void PrintOutHdf5_Grids(int timeline,
                        GridsPoints *****ptrArray)
{
    using namespace H5;

    H5std_string FILE_NAME("./Data/GridsData.h5");
    H5std_string DATASET_NAME("Grids");

    H5std_string MEMBERx("x");
    H5std_string MEMBERy("y");
    H5std_string MEMBERz("z");

    H5std_string MEMBER_pos3("pos3");

    H5std_string MEMBER_e3("e3");

    H5std_string MEMBER_b3("b3");
    H5std_string MEMBER_dB3("dB3");

    H5std_string MEMBER_v3("v3");
    H5std_string MEMBER_vH3("vH3");
    H5std_string MEMBER_vHe3("vHe3");
    H5std_string MEMBER_vO3("vO3");

    H5std_string MEMBER_ve3("ve3");

    H5std_string MEMBER_gradB3("gradB3");
    H5std_string MEMBER_gradPe("gradPe");

    H5std_string MEMBER_density_H("densityH");
    H5std_string MEMBER_density_He("densityHe");
    H5std_string MEMBER_density_O("densityO");

    H5std_string MEMBER_te("temperature");
    H5std_string MEMBER_stopSign("stopSign");

    typedef struct //Vector3_h5
    {
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    typedef struct //GridsPoints_h5
    {
        Vector3_h5 pos3;

        Vector3_h5 e3;

        Vector3_h5 b3;
        Vector3_h5 dB3;

        Vector3_h5 vH3;
        Vector3_h5 vHe3;
        Vector3_h5 vO3;

        Vector3_h5 ve3;

        Vector3_h5 gradB3;
        Vector3_h5 gradPe;

        double density_H;
        double density_He;
        double density_O;

        double temperature;
        int stopSign;
    } GridsPoints_h5;

    // Apply for continus memory
    GridsPoints_h5 *data_mem = new GridsPoints_h5[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize)];
    GridsPoints_h5 ****array_data = new GridsPoints_h5 ***[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        array_data[face] = new GridsPoints_h5 **[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            array_data[face][i] = new GridsPoints_h5 *[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                array_data[face][i][j] = data_mem + face * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) +
                                         i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) +
                                         j * (1 + fieldsGridsSize);
                for (int k = 0; k < 1 + fieldsGridsSize; k++)
                {
                    array_data[face][i][j][k].pos3 = {ptrArray[face][i + 1][j + 1][k]->Pos3().x(), ptrArray[face][i + 1][j + 1][k]->Pos3().y(), ptrArray[face][i + 1][j + 1][k]->Pos3().z()};
                    array_data[face][i][j][k].e3 = {ptrArray[face][i + 1][j + 1][k]->E3().x(), ptrArray[face][i + 1][j + 1][k]->E3().y(), ptrArray[face][i + 1][j + 1][k]->E3().z()};
                    array_data[face][i][j][k].b3 = {ptrArray[face][i + 1][j + 1][k]->B3_base().x(), ptrArray[face][i + 1][j + 1][k]->B3_base().y(), ptrArray[face][i + 1][j + 1][k]->B3_base().z()};
                    array_data[face][i][j][k].dB3 = {ptrArray[face][i + 1][j + 1][k]->DB3().x(), ptrArray[face][i + 1][j + 1][k]->DB3().y(), ptrArray[face][i + 1][j + 1][k]->DB3().z()};

                    array_data[face][i][j][k].vH3 = {ptrArray[face][i + 1][j + 1][k]->VelH3().x(), ptrArray[face][i + 1][j + 1][k]->VelH3().y(), ptrArray[face][i + 1][j + 1][k]->VelH3().z()};
                    array_data[face][i][j][k].vHe3 = {ptrArray[face][i + 1][j + 1][k]->VelHe3().x(), ptrArray[face][i + 1][j + 1][k]->VelHe3().y(), ptrArray[face][i + 1][j + 1][k]->VelHe3().z()};
                    array_data[face][i][j][k].vO3 = {ptrArray[face][i + 1][j + 1][k]->VelO3().x(), ptrArray[face][i + 1][j + 1][k]->VelO3().y(), ptrArray[face][i + 1][j + 1][k]->VelO3().z()};
                    array_data[face][i][j][k].ve3 = {ptrArray[face][i + 1][j + 1][k]->Vel_e3().x(), ptrArray[face][i + 1][j + 1][k]->Vel_e3().y(), ptrArray[face][i + 1][j + 1][k]->Vel_e3().z()};

                    array_data[face][i][j][k].gradB3 = {ptrArray[face][i + 1][j + 1][k]->GradB3().x(), ptrArray[face][i + 1][j + 1][k]->GradB3().y(), ptrArray[face][i + 1][j + 1][k]->GradB3().z()};
                    array_data[face][i][j][k].gradPe = {ptrArray[face][i + 1][j + 1][k]->GradPe().x(), ptrArray[face][i + 1][j + 1][k]->GradPe().y(), ptrArray[face][i + 1][j + 1][k]->GradPe().z()};

                    array_data[face][i][j][k].density_H = ptrArray[face][i + 1][j + 1][k]->Density_H();
                    array_data[face][i][j][k].density_He = ptrArray[face][i + 1][j + 1][k]->Density_He();
                    array_data[face][i][j][k].density_O = ptrArray[face][i + 1][j + 1][k]->Density_O();

                    array_data[face][i][j][k].temperature = ptrArray[face][i + 1][j + 1][k]->Temperature();
                    array_data[face][i][j][k].stopSign = 0;
                }
            }
        }
    }

    Exception::dontPrint();

    int RANK = 4;
    hsize_t dim[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize + 1, (hsize_t)fieldsGridsSize + 1, 
                     (hsize_t)fieldsGridsSize + 1};
    DataSpace space(RANK, dim);

    // vector elements
    CompType mtype_vector3(sizeof(Vector3_h5));
    mtype_vector3.insertMember(MEMBERx, HOFFSET(Vector3_h5, v_x), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERy, HOFFSET(Vector3_h5, v_y), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERz, HOFFSET(Vector3_h5, v_z), PredType::NATIVE_DOUBLE);

    CompType mtype_grids(sizeof(GridsPoints_h5));
    mtype_grids.insertMember(MEMBER_pos3, HOFFSET(GridsPoints_h5, pos3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_e3, HOFFSET(GridsPoints_h5, e3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_b3, HOFFSET(GridsPoints_h5, b3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_dB3, HOFFSET(GridsPoints_h5, dB3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_vH3, HOFFSET(GridsPoints_h5, vH3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_vHe3, HOFFSET(GridsPoints_h5, vHe3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_vO3, HOFFSET(GridsPoints_h5, vO3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_ve3, HOFFSET(GridsPoints_h5, ve3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_gradB3, HOFFSET(GridsPoints_h5, gradB3), mtype_vector3);
    mtype_grids.insertMember(MEMBER_gradPe, HOFFSET(GridsPoints_h5, gradPe), mtype_vector3);
    mtype_grids.insertMember(MEMBER_density_H, HOFFSET(GridsPoints_h5, density_H), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_density_He, HOFFSET(GridsPoints_h5, density_He), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_density_O, HOFFSET(GridsPoints_h5, density_O), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_te, HOFFSET(GridsPoints_h5, temperature), PredType::NATIVE_DOUBLE);
    mtype_grids.insertMember(MEMBER_stopSign, HOFFSET(GridsPoints_h5, stopSign), PredType::NATIVE_INT);

    H5File file_Grids(FILE_NAME, H5F_ACC_TRUNC);

    DataSpace dataspace_Grids = DataSpace(RANK, dim);
    DataSet dataSetGrids;
    dataSetGrids = file_Grids.createDataSet(DATASET_NAME, mtype_grids, dataspace_Grids);
    dataSetGrids.write(array_data[0][0][0], mtype_grids);
    dataspace_Grids.close();
    dataSetGrids.close();

    delete data_mem;
}

//************************************************************************
//  Save grids and particle data for continuing run
void PrintOutHdf5_Particles_Grids(int timeline_in,
                                  GridsCells ****ptrArrayCells,
                                  GridsPoints *****ptrArray)
{

    using namespace H5;

    char filen[80], filename[151];

    sprintf(filen, "particles%d.h5", timeline_in);
    strncpy(filename, outpdir, 150);
    strncat(filename, filen, 150);

    H5std_string FILE_NAME(filename);

    H5std_string DATASET_NAME_H("Particles_H"); // pos, v, weight, mu
    H5std_string DATASET_NAME_He("Particles_He");
    H5std_string DATASET_NAME_O("Particles_O");
    H5std_string DATASET_NAME_H_Uint("Particles_H_Uint"); // Uint of pos
    H5std_string DATASET_NAME_He_Uint("Particles_He_Uint");
    H5std_string DATASET_NAME_O_Uint("Particles_O_Uint");

    H5std_string DATASET_NAME_Grids("Grids");
    H5std_string DATASET_NAME_Time("Time");

    int numParticle_H = 0;
    int numParticle_He = 0;
    int numParticle_O = 0;
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    numParticle_H += ptrArrayCells[f][i][j][k].Particles_H()->size();
                    numParticle_He += ptrArrayCells[f][i][j][k].Particles_He()->size();
                    numParticle_O += ptrArrayCells[f][i][j][k].Particles_O()->size();
                }
            }
        }
    }

    int par_size;
    if (numParticle_H >= numParticle_He)
    {
        par_size = numParticle_H;
    }
    else
    {
        par_size = numParticle_He;
    }
    if (par_size <= numParticle_O)
    {
        par_size = numParticle_O;
    }

    double *data_mem_particles = new double[8 * par_size];
    double **array_data_particles = new double *[par_size];
    unsigned long long *array_data_particles_Uint = new unsigned long long[par_size];

    int RANK = 2;
    hsize_t dimsH[2];
    hsize_t dimsHe[2];
    hsize_t dimsO[2];
    dimsH[0] = numParticle_H;
    dimsH[1] = 8;

    dimsHe[0] = numParticle_He;
    dimsHe[1] = 8;

    dimsO[0] = numParticle_O;
    dimsO[1] = 8;

    hsize_t dimsH_Uint[2];
    hsize_t dimsHe_Uint[2];
    hsize_t dimsO_Uint[2];

    dimsH_Uint[0] = numParticle_H;
    dimsH_Uint[1] = 1;

    dimsHe_Uint[0] = numParticle_He;
    dimsHe_Uint[1] = 1;

    dimsO_Uint[0] = numParticle_O;
    dimsO_Uint[1] = 1;

    H5File file_saved(FILE_NAME, H5F_ACC_TRUNC);

    int index = 0;
    vector<Particles> *ptrPar = NULL;
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    ptrPar = ptrArrayCells[f][i][j][k].Particles_H();
                    //
                    for (long unsigned int iter = 0; iter < ptrPar->size(); ++iter)
                    {
                        array_data_particles_Uint[index] = (*ptrPar)[iter].PosUint();
                        array_data_particles[index] = data_mem_particles + 8 * index;
                        array_data_particles[index][0] = (*ptrPar)[iter].PosParticles().x();
                        array_data_particles[index][1] = (*ptrPar)[iter].PosParticles().y();
                        array_data_particles[index][2] = (*ptrPar)[iter].PosParticles().z();
                        array_data_particles[index][3] = (*ptrPar)[iter].VelParticles().x();
                        array_data_particles[index][4] = (*ptrPar)[iter].VelParticles().y();
                        array_data_particles[index][5] = (*ptrPar)[iter].VelParticles().z();
                        array_data_particles[index][6] = (*ptrPar)[iter].WeightNi();
                        array_data_particles[index][7] = (*ptrPar)[iter].MagneticIvarient();
                        //
                        index += 1;
                    }
                }
            }
        }
    }
    //
    DataSpace dataspace_particles_H = DataSpace(RANK, dimsH);
    DataSet dataset_H;
    dataset_H = file_saved.createDataSet(DATASET_NAME_H, PredType::NATIVE_DOUBLE, dataspace_particles_H);
    dataset_H.write(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataspace_particles_H.close();
    dataset_H.close();
    //
    DataSpace dataspace_particles_H_Uint = DataSpace(RANK, dimsH_Uint);
    DataSet dataset_H_Uint;
    dataset_H_Uint = file_saved.createDataSet(DATASET_NAME_H_Uint, PredType::NATIVE_ULLONG, dataspace_particles_H_Uint);
    dataset_H_Uint.write(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    dataspace_particles_H_Uint.close();
    dataset_H_Uint.close();
    //
    index = 0;
    ptrPar = NULL;
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    ptrPar = ptrArrayCells[f][i][j][k].Particles_He();
                    //
                    for (long unsigned int iter = 0; iter < ptrPar->size(); ++iter)
                    {
                        array_data_particles_Uint[index] = (*ptrPar)[iter].PosUint();
                        array_data_particles[index] = data_mem_particles + 8 * index;
                        array_data_particles[index][0] = (*ptrPar)[iter].PosParticles().x();
                        array_data_particles[index][1] = (*ptrPar)[iter].PosParticles().y();
                        array_data_particles[index][2] = (*ptrPar)[iter].PosParticles().z();
                        array_data_particles[index][3] = (*ptrPar)[iter].VelParticles().x();
                        array_data_particles[index][4] = (*ptrPar)[iter].VelParticles().y();
                        array_data_particles[index][5] = (*ptrPar)[iter].VelParticles().z();
                        array_data_particles[index][6] = (*ptrPar)[iter].WeightNi();
                        array_data_particles[index][7] = (*ptrPar)[iter].MagneticIvarient();
                        //
                        index += 1;
                    }
                }
            }
        }
    }
    //
    DataSpace dataspace_particles_He = DataSpace(RANK, dimsHe);
    DataSet dataset_He;
    dataset_He = file_saved.createDataSet(DATASET_NAME_He, PredType::NATIVE_DOUBLE, dataspace_particles_He);
    dataset_He.write(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataspace_particles_He.close();
    dataset_He.close();
    //
    DataSpace dataspace_particles_He_Uint = DataSpace(RANK, dimsHe_Uint);
    DataSet dataset_He_Uint;
    dataset_He_Uint = file_saved.createDataSet(DATASET_NAME_He_Uint, PredType::NATIVE_ULLONG, dataspace_particles_He_Uint);
    dataset_He_Uint.write(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    dataspace_particles_He_Uint.close();
    dataset_He_Uint.close();
    //
    index = 0;
    ptrPar = NULL;
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    ptrPar = ptrArrayCells[f][i][j][k].Particles_O();
                    //
                    for (long unsigned int iter = 0; iter < ptrPar->size(); ++iter)
                    {
                        array_data_particles_Uint[index] = (*ptrPar)[iter].PosUint();
                        array_data_particles[index] = data_mem_particles + 8 * index;
                        array_data_particles[index][0] = (*ptrPar)[iter].PosParticles().x();
                        array_data_particles[index][1] = (*ptrPar)[iter].PosParticles().y();
                        array_data_particles[index][2] = (*ptrPar)[iter].PosParticles().z();
                        array_data_particles[index][3] = (*ptrPar)[iter].VelParticles().x();
                        array_data_particles[index][4] = (*ptrPar)[iter].VelParticles().y();
                        array_data_particles[index][5] = (*ptrPar)[iter].VelParticles().z();
                        array_data_particles[index][6] = (*ptrPar)[iter].WeightNi();
                        array_data_particles[index][7] = (*ptrPar)[iter].MagneticIvarient();
                        //
                        index += 1;
                    }
                }
            }
        }
    }
    //
    DataSpace dataspace_particles_O = DataSpace(RANK, dimsO);
    DataSet dataset_O;
    dataset_O = file_saved.createDataSet(DATASET_NAME_O, PredType::NATIVE_DOUBLE, dataspace_particles_O);
    dataset_O.write(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataspace_particles_O.close();
    dataset_O.close();

    DataSpace dataspace_particles_O_Uint = DataSpace(RANK, dimsO_Uint);
    DataSet dataset_O_Uint;
    dataset_O_Uint = file_saved.createDataSet(DATASET_NAME_O_Uint, PredType::NATIVE_ULLONG, dataspace_particles_O_Uint);
    dataset_O_Uint.write(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    dataspace_particles_O_Uint.close();
    dataset_O_Uint.close();

    cout << " PrintParticles: H " << numParticle_H << " He " << numParticle_He << " O " << numParticle_O << endl;

    delete[] data_mem_particles;
    delete[] array_data_particles;
    delete[] array_data_particles_Uint;

    // Apply for continus memory

    double *data_mem_grids = new double[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 57];
    double *****data_grids = new double ****[totalFace];

    for (int f = 0; f < totalFace; f++)
    {
        data_grids[f] = new double ***[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            data_grids[f][i] = new double **[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                data_grids[f][i][j] = new double *[1 + fieldsGridsSize * grid_domain];
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    data_grids[f][i][j][k] = data_mem_grids + f * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 57 +
                                             i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 57 +
                                             j * (1 + fieldsGridsSize * grid_domain) * 57 +
                                             k * 57;
                    // pos
                    data_grids[f][i][j][k][0] = ptrArray[f][i + 1][j + 1][k]->Pos3().x();
                    data_grids[f][i][j][k][1] = ptrArray[f][i + 1][j + 1][k]->Pos3().y();
                    data_grids[f][i][j][k][2] = ptrArray[f][i + 1][j + 1][k]->Pos3().z();
                    // e3
                    data_grids[f][i][j][k][3] = ptrArray[f][i + 1][j + 1][k]->E3().x();
                    data_grids[f][i][j][k][4] = ptrArray[f][i + 1][j + 1][k]->E3().y();
                    data_grids[f][i][j][k][5] = ptrArray[f][i + 1][j + 1][k]->E3().z();
                    // dE3
                    data_grids[f][i][j][k][6] = ptrArray[f][i + 1][j + 1][k]->E3().x();
                    data_grids[f][i][j][k][7] = ptrArray[f][i + 1][j + 1][k]->E3().y();
                    data_grids[f][i][j][k][8] = ptrArray[f][i + 1][j + 1][k]->E3().z();
                    // b3
                    data_grids[f][i][j][k][9] = ptrArray[f][i + 1][j + 1][k]->B3_base().x();
                    data_grids[f][i][j][k][10] = ptrArray[f][i + 1][j + 1][k]->B3_base().y();
                    data_grids[f][i][j][k][11] = ptrArray[f][i + 1][j + 1][k]->B3_base().z();
                    // dB3
                    data_grids[f][i][j][k][12] = ptrArray[f][i + 1][j + 1][k]->DB3().x();
                    data_grids[f][i][j][k][13] = ptrArray[f][i + 1][j + 1][k]->DB3().y();
                    data_grids[f][i][j][k][14] = ptrArray[f][i + 1][j + 1][k]->DB3().z();
                    // v3
                    data_grids[f][i][j][k][15] = ptrArray[f][i + 1][j + 1][k]->Vel3().x();
                    data_grids[f][i][j][k][16] = ptrArray[f][i + 1][j + 1][k]->Vel3().y();
                    data_grids[f][i][j][k][17] = ptrArray[f][i + 1][j + 1][k]->Vel3().z();
                    // vH3
                    data_grids[f][i][j][k][18] = ptrArray[f][i + 1][j + 1][k]->VelH3().x();
                    data_grids[f][i][j][k][19] = ptrArray[f][i + 1][j + 1][k]->VelH3().y();
                    data_grids[f][i][j][k][20] = ptrArray[f][i + 1][j + 1][k]->VelH3().z();
                    // vHe3
                    data_grids[f][i][j][k][21] = ptrArray[f][i + 1][j + 1][k]->VelHe3().x();
                    data_grids[f][i][j][k][22] = ptrArray[f][i + 1][j + 1][k]->VelHe3().y();
                    data_grids[f][i][j][k][23] = ptrArray[f][i + 1][j + 1][k]->VelHe3().z();
                    // vO3
                    data_grids[f][i][j][k][24] = ptrArray[f][i + 1][j + 1][k]->VelO3().x();
                    data_grids[f][i][j][k][25] = ptrArray[f][i + 1][j + 1][k]->VelO3().y();
                    data_grids[f][i][j][k][26] = ptrArray[f][i + 1][j + 1][k]->VelO3().z();
                    // vH3_cumu
                    data_grids[f][i][j][k][27] = ptrArray[f][i + 1][j + 1][k]->VelH3_cumu().x();
                    data_grids[f][i][j][k][28] = ptrArray[f][i + 1][j + 1][k]->VelH3_cumu().y();
                    data_grids[f][i][j][k][29] = ptrArray[f][i + 1][j + 1][k]->VelH3_cumu().z();
                    // vHe3_cumu
                    data_grids[f][i][j][k][30] = ptrArray[f][i + 1][j + 1][k]->VelHe3_cumu().x();
                    data_grids[f][i][j][k][31] = ptrArray[f][i + 1][j + 1][k]->VelHe3_cumu().y();
                    data_grids[f][i][j][k][32] = ptrArray[f][i + 1][j + 1][k]->VelHe3_cumu().z();
                    // vO3_cumu
                    data_grids[f][i][j][k][33] = ptrArray[f][i + 1][j + 1][k]->VelO3_cumu().x();
                    data_grids[f][i][j][k][34] = ptrArray[f][i + 1][j + 1][k]->VelO3_cumu().y();
                    data_grids[f][i][j][k][35] = ptrArray[f][i + 1][j + 1][k]->VelO3_cumu().z();
                    // ve3
                    data_grids[f][i][j][k][36] = ptrArray[f][i + 1][j + 1][k]->Vel_e3().x();
                    data_grids[f][i][j][k][37] = ptrArray[f][i + 1][j + 1][k]->Vel_e3().y();
                    data_grids[f][i][j][k][38] = ptrArray[f][i + 1][j + 1][k]->Vel_e3().z();
                    // gradB3
                    data_grids[f][i][j][k][39] = ptrArray[f][i + 1][j + 1][k]->GradB3().x();
                    data_grids[f][i][j][k][40] = ptrArray[f][i + 1][j + 1][k]->GradB3().y();
                    data_grids[f][i][j][k][41] = ptrArray[f][i + 1][j + 1][k]->GradB3().z();
                    // gradPe
                    data_grids[f][i][j][k][42] = ptrArray[f][i + 1][j + 1][k]->GradPe().x();
                    data_grids[f][i][j][k][43] = ptrArray[f][i + 1][j + 1][k]->GradPe().y();
                    data_grids[f][i][j][k][44] = ptrArray[f][i + 1][j + 1][k]->GradPe().z();
                    // density
                    data_grids[f][i][j][k][45] = ptrArray[f][i + 1][j + 1][k]->Density_H();
                    data_grids[f][i][j][k][46] = ptrArray[f][i + 1][j + 1][k]->Density_He();
                    data_grids[f][i][j][k][47] = ptrArray[f][i + 1][j + 1][k]->Density_O();
                    // density_cumu
                    data_grids[f][i][j][k][48] = ptrArray[f][i + 1][j + 1][k]->Density_H_cumu();
                    data_grids[f][i][j][k][49] = ptrArray[f][i + 1][j + 1][k]->Density_He_cumu();
                    data_grids[f][i][j][k][50] = ptrArray[f][i + 1][j + 1][k]->Density_O_cumu();
                    // temperature
                    data_grids[f][i][j][k][51] = ptrArray[f][i + 1][j + 1][k]->Temperature();
                    data_grids[f][i][j][k][52] = ptrArray[f][i + 1][j + 1][k]->Temperature_H();
                    data_grids[f][i][j][k][53] = ptrArray[f][i + 1][j + 1][k]->Temperature_He();
                    data_grids[f][i][j][k][54] = ptrArray[f][i + 1][j + 1][k]->Temperature_O();
                    // potential
                    data_grids[f][i][j][k][55] = ptrArray[f][i + 1][j + 1][k]->Potential();
                    //  psd_H
                    data_grids[f][i][j][k][56] = ptrArray[f][i + 1][j + 1][k]->PSD_H();
                }
            }
        }
    }
    int RANK_grids = 5;
    hsize_t dim_grids[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize + 1, (hsize_t)fieldsGridsSize + 1, 
                           (hsize_t)fieldsGridsSize * grid_domain + 1, 52};
    DataSpace dataspace_grids = DataSpace(RANK_grids, dim_grids);
    DataSet dataSetGrids;
    dataSetGrids = file_saved.createDataSet(DATASET_NAME_Grids, PredType::NATIVE_DOUBLE, dataspace_grids);
    dataSetGrids.write(data_grids[0][0][0][0], PredType::NATIVE_DOUBLE);
    dataspace_grids.close();
    dataSetGrids.close();
    //
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                delete[] data_grids[f][i][j];
            }
            delete[] data_grids[f][i];
        }
        delete[] data_grids[f];
    }
    delete[] data_grids;
    delete[] data_mem_grids;

    int data_time = timeline_in;
    if (continueParticles == 0)
        data_time = 0; // check if to be a initial start timaline

    int RANK_time = 1;        // rank
    hsize_t dim_time[] = {1}; // hsize dim
    DataSpace dataspace_time = DataSpace(RANK_time, dim_time);
    DataSet dataSetTime;
    dataSetTime = file_saved.createDataSet(DATASET_NAME_Time, PredType::NATIVE_INT, dataspace_time);
    dataSetTime.write(&data_time, PredType::NATIVE_INT);
    dataspace_time.close();
    dataSetTime.close();

    file_saved.close();
}
//************************************************************************
// Output data of particles and grids for continuing run: overloaded print function
void PrintOutHdf5_Particles_Grids(int timeline_in,
                                  vector<Particles> &ptrParticlesList_H,
                                  vector<Particles> &ptrParticlesList_He,
                                  vector<Particles> &ptrParticlesList_O,
                                  GridsPoints *****ptrArray)
{
    using namespace H5;

    char filen[80], filename[151];

    sprintf(filen, "particles_save%d.h5", timeline_in);
    strncpy(filename, outpdir, 150);
    strncat(filename, filen, 150);

    H5std_string FILE_NAME(filename);
    H5std_string DATASET_NAME_H("Particles_H"); // pos, v, weight, mu
    H5std_string DATASET_NAME_He("Particles_He");
    H5std_string DATASET_NAME_O("Particles_O");
    H5std_string DATASET_NAME_H_Uint("Particles_H_Uint"); // Uint of pos
    H5std_string DATASET_NAME_He_Uint("Particles_He_Uint");
    H5std_string DATASET_NAME_O_Uint("Particles_O_Uint");

    H5std_string DATASET_NAME_Grids("Grids");
    H5std_string DATASET_NAME_Time("Time");

    int numParticle_H = ptrParticlesList_H.size();
    int numParticle_He = ptrParticlesList_He.size();
    int numParticle_O = ptrParticlesList_O.size();
    int par_size;
    if (numParticle_H >= numParticle_He)
    {
        par_size = numParticle_H;
    }
    else
    {
        par_size = numParticle_He;
    }
    if (par_size <= numParticle_O)
    {
        par_size = numParticle_O;
    }

    double *data_mem_particles = new double[8 * par_size];
    double **array_data_particles = new double *[par_size];
    unsigned long long *array_data_particles_Uint = new unsigned long long[par_size];

    int RANK = 2;
    hsize_t dimsH[2];
    hsize_t dimsHe[2];
    hsize_t dimsO[2];
    dimsH[0] = numParticle_H;
    dimsH[1] = 8;

    dimsHe[0] = numParticle_He;
    dimsHe[1] = 8;

    dimsO[0] = numParticle_O;
    dimsO[1] = 8;

    hsize_t dimsH_Uint[2];
    hsize_t dimsHe_Uint[2];
    hsize_t dimsO_Uint[2];

    dimsH_Uint[0] = numParticle_H;
    dimsH_Uint[1] = 1;

    dimsHe_Uint[0] = numParticle_He;
    dimsHe_Uint[1] = 1;

    dimsO_Uint[0] = numParticle_O;
    dimsO_Uint[1] = 1;

    H5File file_saved(FILE_NAME, H5F_ACC_TRUNC);

#pragma omp parallel for
    for (int i = 0; i < numParticle_H; i++)
    {
        array_data_particles_Uint[i] = ptrParticlesList_H[i].PosUint();
        array_data_particles[i] = data_mem_particles + 8 * i;
        array_data_particles[i][0] = ptrParticlesList_H[i].PosParticles().x();
        array_data_particles[i][1] = ptrParticlesList_H[i].PosParticles().y();
        array_data_particles[i][2] = ptrParticlesList_H[i].PosParticles().z();
        array_data_particles[i][3] = ptrParticlesList_H[i].VelParticles().x();
        array_data_particles[i][4] = ptrParticlesList_H[i].VelParticles().y();
        array_data_particles[i][5] = ptrParticlesList_H[i].VelParticles().z();
        array_data_particles[i][6] = ptrParticlesList_H[i].WeightNi();
        array_data_particles[i][7] = ptrParticlesList_H[i].MagneticIvarient();
    }

    DataSpace dataspace_particles_H = DataSpace(RANK, dimsH);
    DataSet dataset_H;
    dataset_H = file_saved.createDataSet(DATASET_NAME_H, PredType::NATIVE_DOUBLE, dataspace_particles_H);
    dataset_H.write(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataspace_particles_H.close();
    dataset_H.close();

    DataSpace dataspace_particles_H_Uint = DataSpace(RANK, dimsH_Uint);
    DataSet dataset_H_Uint;
    dataset_H_Uint = file_saved.createDataSet(DATASET_NAME_H_Uint, PredType::NATIVE_ULLONG, dataspace_particles_H_Uint);
    dataset_H_Uint.write(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    dataspace_particles_H_Uint.close();
    dataset_H_Uint.close();

#pragma omp parallel for
    for (int i = 0; i < numParticle_He; i++)
    {
        array_data_particles_Uint[i] = ptrParticlesList_He[i].PosUint();
        array_data_particles[i] = data_mem_particles + 8 * i;
        array_data_particles[i][0] = ptrParticlesList_He[i].PosParticles().x();
        array_data_particles[i][1] = ptrParticlesList_He[i].PosParticles().y();
        array_data_particles[i][2] = ptrParticlesList_He[i].PosParticles().z();
        array_data_particles[i][3] = ptrParticlesList_He[i].VelParticles().x();
        array_data_particles[i][4] = ptrParticlesList_He[i].VelParticles().y();
        array_data_particles[i][5] = ptrParticlesList_He[i].VelParticles().z();
        array_data_particles[i][6] = ptrParticlesList_He[i].WeightNi();
        array_data_particles[i][7] = ptrParticlesList_He[i].MagneticIvarient();
    }
    DataSpace dataspace_particles_He = DataSpace(RANK, dimsHe);
    DataSet dataset_He;
    dataset_He = file_saved.createDataSet(DATASET_NAME_He, PredType::NATIVE_DOUBLE, dataspace_particles_He);
    dataset_He.write(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataspace_particles_He.close();
    dataset_He.close();

    DataSpace dataspace_particles_He_Uint = DataSpace(RANK, dimsHe_Uint);
    DataSet dataset_He_Uint;
    dataset_He_Uint = file_saved.createDataSet(DATASET_NAME_He_Uint, PredType::NATIVE_ULLONG, dataspace_particles_He_Uint);
    dataset_He_Uint.write(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    dataspace_particles_He_Uint.close();
    dataset_He_Uint.close();

#pragma omp parallel for
    for (int i = 0; i < numParticle_O; i++)
    {
        array_data_particles_Uint[i] = ptrParticlesList_O[i].PosUint();
        array_data_particles[i] = data_mem_particles + 8 * i;
        array_data_particles[i][0] = ptrParticlesList_O[i].PosParticles().x();
        array_data_particles[i][1] = ptrParticlesList_O[i].PosParticles().y();
        array_data_particles[i][2] = ptrParticlesList_O[i].PosParticles().z();
        array_data_particles[i][3] = ptrParticlesList_O[i].VelParticles().x();
        array_data_particles[i][4] = ptrParticlesList_O[i].VelParticles().y();
        array_data_particles[i][5] = ptrParticlesList_O[i].VelParticles().z();
        array_data_particles[i][6] = ptrParticlesList_O[i].WeightNi();
        array_data_particles[i][7] = ptrParticlesList_O[i].MagneticIvarient();
    }
    DataSpace dataspace_particles_O = DataSpace(RANK, dimsO);
    DataSet dataset_O;
    dataset_O = file_saved.createDataSet(DATASET_NAME_O, PredType::NATIVE_DOUBLE, dataspace_particles_O);
    dataset_O.write(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataspace_particles_O.close();
    dataset_O.close();

    DataSpace dataspace_particles_O_Uint = DataSpace(RANK, dimsO_Uint);
    DataSet dataset_O_Uint;
    dataset_O_Uint = file_saved.createDataSet(DATASET_NAME_O_Uint, PredType::NATIVE_ULLONG, dataspace_particles_O_Uint);
    dataset_O_Uint.write(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    dataspace_particles_O_Uint.close();
    dataset_O_Uint.close();

    cout << " PrintParticles: H " << numParticle_H << " He " << numParticle_He << " O " << numParticle_O << endl;

    delete data_mem_particles;
    delete array_data_particles;
    delete array_data_particles_Uint;

    // Apply for continus memory

    double *data_mem_grids = new double[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52];
    double *****data_grids = new double ****[totalFace];

    for (int f = 0; f < totalFace; f++)
    {
        data_grids[f] = new double ***[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            data_grids[f][i] = new double **[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                data_grids[f][i][j] = new double *[1 + fieldsGridsSize * grid_domain];
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    data_grids[f][i][j][k] = data_mem_grids + f * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             j * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             k * 57;
                    // pos
                    data_grids[f][i][j][k][0] = ptrArray[f][i + 1][j + 1][k]->Pos3().x();
                    data_grids[f][i][j][k][1] = ptrArray[f][i + 1][j + 1][k]->Pos3().y();
                    data_grids[f][i][j][k][2] = ptrArray[f][i + 1][j + 1][k]->Pos3().z();
                    // e3
                    data_grids[f][i][j][k][3] = ptrArray[f][i + 1][j + 1][k]->E3().x();
                    data_grids[f][i][j][k][4] = ptrArray[f][i + 1][j + 1][k]->E3().y();
                    data_grids[f][i][j][k][5] = ptrArray[f][i + 1][j + 1][k]->E3().z();
                    // dE3
                    data_grids[f][i][j][k][6] = ptrArray[f][i + 1][j + 1][k]->E3().x();
                    data_grids[f][i][j][k][7] = ptrArray[f][i + 1][j + 1][k]->E3().y();
                    data_grids[f][i][j][k][8] = ptrArray[f][i + 1][j + 1][k]->E3().z();
                    // b3
                    data_grids[f][i][j][k][9] = ptrArray[f][i + 1][j + 1][k]->B3_base().x();
                    data_grids[f][i][j][k][10] = ptrArray[f][i + 1][j + 1][k]->B3_base().y();
                    data_grids[f][i][j][k][11] = ptrArray[f][i + 1][j + 1][k]->B3_base().z();
                    // dB3
                    data_grids[f][i][j][k][12] = ptrArray[f][i + 1][j + 1][k]->DB3().x();
                    data_grids[f][i][j][k][13] = ptrArray[f][i + 1][j + 1][k]->DB3().y();
                    data_grids[f][i][j][k][14] = ptrArray[f][i + 1][j + 1][k]->DB3().z();
                    // v3
                    data_grids[f][i][j][k][15] = ptrArray[f][i + 1][j + 1][k]->Vel3().x();
                    data_grids[f][i][j][k][16] = ptrArray[f][i + 1][j + 1][k]->Vel3().y();
                    data_grids[f][i][j][k][17] = ptrArray[f][i + 1][j + 1][k]->Vel3().z();
                    // vH3
                    data_grids[f][i][j][k][18] = ptrArray[f][i + 1][j + 1][k]->VelH3().x();
                    data_grids[f][i][j][k][19] = ptrArray[f][i + 1][j + 1][k]->VelH3().y();
                    data_grids[f][i][j][k][20] = ptrArray[f][i + 1][j + 1][k]->VelH3().z();
                    // vHe3
                    data_grids[f][i][j][k][21] = ptrArray[f][i + 1][j + 1][k]->VelHe3().x();
                    data_grids[f][i][j][k][22] = ptrArray[f][i + 1][j + 1][k]->VelHe3().y();
                    data_grids[f][i][j][k][23] = ptrArray[f][i + 1][j + 1][k]->VelHe3().z();
                    // vO3
                    data_grids[f][i][j][k][24] = ptrArray[f][i + 1][j + 1][k]->VelO3().x();
                    data_grids[f][i][j][k][25] = ptrArray[f][i + 1][j + 1][k]->VelO3().y();
                    data_grids[f][i][j][k][26] = ptrArray[f][i + 1][j + 1][k]->VelO3().z();
                    // vH3_cumu
                    data_grids[f][i][j][k][27] = ptrArray[f][i + 1][j + 1][k]->VelH3_cumu().x();
                    data_grids[f][i][j][k][28] = ptrArray[f][i + 1][j + 1][k]->VelH3_cumu().y();
                    data_grids[f][i][j][k][29] = ptrArray[f][i + 1][j + 1][k]->VelH3_cumu().z();
                    // vHe3_cumu
                    data_grids[f][i][j][k][30] = ptrArray[f][i + 1][j + 1][k]->VelHe3_cumu().x();
                    data_grids[f][i][j][k][31] = ptrArray[f][i + 1][j + 1][k]->VelHe3_cumu().y();
                    data_grids[f][i][j][k][32] = ptrArray[f][i + 1][j + 1][k]->VelHe3_cumu().z();
                    // vO3_cumu
                    data_grids[f][i][j][k][33] = ptrArray[f][i + 1][j + 1][k]->VelO3_cumu().x();
                    data_grids[f][i][j][k][34] = ptrArray[f][i + 1][j + 1][k]->VelO3_cumu().y();
                    data_grids[f][i][j][k][35] = ptrArray[f][i + 1][j + 1][k]->VelO3_cumu().z();
                    // ve3
                    data_grids[f][i][j][k][36] = ptrArray[f][i + 1][j + 1][k]->Vel_e3().x();
                    data_grids[f][i][j][k][37] = ptrArray[f][i + 1][j + 1][k]->Vel_e3().y();
                    data_grids[f][i][j][k][38] = ptrArray[f][i + 1][j + 1][k]->Vel_e3().z();
                    // gradB3
                    data_grids[f][i][j][k][39] = ptrArray[f][i + 1][j + 1][k]->GradB3().x();
                    data_grids[f][i][j][k][40] = ptrArray[f][i + 1][j + 1][k]->GradB3().y();
                    data_grids[f][i][j][k][41] = ptrArray[f][i + 1][j + 1][k]->GradB3().z();
                    // gradPe
                    data_grids[f][i][j][k][42] = ptrArray[f][i + 1][j + 1][k]->GradPe().x();
                    data_grids[f][i][j][k][43] = ptrArray[f][i + 1][j + 1][k]->GradPe().y();
                    data_grids[f][i][j][k][44] = ptrArray[f][i + 1][j + 1][k]->GradPe().z();
                    // density
                    data_grids[f][i][j][k][45] = ptrArray[f][i + 1][j + 1][k]->Density_H();
                    data_grids[f][i][j][k][46] = ptrArray[f][i + 1][j + 1][k]->Density_He();
                    data_grids[f][i][j][k][47] = ptrArray[f][i + 1][j + 1][k]->Density_O();
                    // density_cumu
                    data_grids[f][i][j][k][48] = ptrArray[f][i + 1][j + 1][k]->Density_H_cumu();
                    data_grids[f][i][j][k][49] = ptrArray[f][i + 1][j + 1][k]->Density_He_cumu();
                    data_grids[f][i][j][k][50] = ptrArray[f][i + 1][j + 1][k]->Density_O_cumu();
                    // temperature
                    data_grids[f][i][j][k][51] = ptrArray[f][i + 1][j + 1][k]->Temperature();
                    data_grids[f][i][j][k][52] = ptrArray[f][i + 1][j + 1][k]->Temperature_H();
                    data_grids[f][i][j][k][53] = ptrArray[f][i + 1][j + 1][k]->Temperature_He();
                    data_grids[f][i][j][k][54] = ptrArray[f][i + 1][j + 1][k]->Temperature_O();
                    // potential
                    data_grids[f][i][j][k][55] = ptrArray[f][i + 1][j + 1][k]->Potential();
                    //  psd_H
                    data_grids[f][i][j][k][56] = ptrArray[f][i + 1][j + 1][k]->PSD_H();
                }
            }
        }
    }
    int RANK_grids = 5;
    hsize_t dim_grids[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize + 1, (hsize_t)fieldsGridsSize + 1, 
                           (hsize_t)fieldsGridsSize * grid_domain + 1, 52};
    DataSpace dataspace_grids = DataSpace(RANK_grids, dim_grids);
    DataSet dataSetGrids;
    dataSetGrids = file_saved.createDataSet(DATASET_NAME_Grids, PredType::NATIVE_DOUBLE, dataspace_grids);
    dataSetGrids.write(data_grids[0][0][0][0], PredType::NATIVE_DOUBLE);
    dataspace_grids.close();
    dataSetGrids.close();

    delete data_mem_grids;
    delete data_grids;

    int data_time = timeline_in;
    if (continueParticles == 0)
        data_time = 0; // check if to be a initial start timaline

    int RANK_time = 1;        // rank
    hsize_t dim_time[] = {1}; // hsize dim
    DataSpace dataspace_time = DataSpace(RANK_time, dim_time);
    DataSet dataSetTime;
    dataSetTime = file_saved.createDataSet(DATASET_NAME_Time, PredType::NATIVE_INT, dataspace_time);
    dataSetTime.write(&data_time, PredType::NATIVE_INT);
    dataspace_time.close();
    dataSetTime.close();

    file_saved.close();
}
//************************************************************************
//************************************************************************
// FUNCTION
// Printout the cell center divB as hdf5 format
// Step 1: Generate a new matrix fulling with  class
// Step 2: Print out it as .h5
void PrintOutHdf5_DivB(int timeline,
                       GridsPoints *****ptrArray,
                       Vector3 *****ptrBVectorFaceArray,
                       double ***ptrVolumeCellArray,
                       Vector3 *****ptrDivBVectorCellArray,
                       Vector3 *****ptrposVectorCellArray)
{

    int i, j, k, face_in;
    for (face_in = 0; face_in < totalFace; face_in++)
    {
        // I, J, K is the cell index
        for (int I = 1; I < fieldsGridsSize + 1; I++)
        {
            for (int J = 1; J < fieldsGridsSize + 1; J++)
            {
                for (int K = 0; K < fieldsGridsSize; K++)
                {
                    i = I - 1;
                    j = J - 1;
                    k = K;

                    double volumetemp = ptrVolumeCellArray[I][J][K];
                    if (volumetemp == 0)
                        continue;
                    double temp = AreaVectorL(ptrArray, face_in, I, J, K).DotProduct(ptrBVectorFaceArray[0][face_in][i][j][k]) +
                                 AreaVectorR(ptrArray, face_in, I, J, K).DotProduct(ptrBVectorFaceArray[0][face_in][i + 1][j][k]) +
                                 AreaVectorT(ptrArray, face_in, I, J, K).DotProduct(ptrBVectorFaceArray[1][face_in][i][j + 1][k]) +
                                 AreaVectorBot(ptrArray, face_in, I, J, K).DotProduct(ptrBVectorFaceArray[1][face_in][i][j][k]) +
                                 AreaVectorF(ptrArray, face_in, I, J, K).DotProduct(ptrBVectorFaceArray[2][face_in][i][j][k + 1]) +
                                 AreaVectorBack(ptrArray, face_in, I, J, K).DotProduct(ptrBVectorFaceArray[2][face_in][i][j][k]);

                    temp = temp / volumetemp;
                    *ptrDivBVectorCellArray[face_in][I][J][K] = Vector3(temp, 0.0, 0.0);

                    if (timeline == 0)
                    {
                        double posx = (ptrArray[face_in][I][J][K]->Pos3().x() +
                                      ptrArray[face_in][I + 1][J][K]->Pos3().x() +
                                      ptrArray[face_in][I][J + 1][K]->Pos3().x() +
                                      ptrArray[face_in][I + 1][J + 1][K]->Pos3().x() +
                                      ptrArray[face_in][I][J][K + 1]->Pos3().x() +
                                      ptrArray[face_in][I + 1][J][K + 1]->Pos3().x() +
                                      ptrArray[face_in][I][J + 1][K + 1]->Pos3().x() +
                                      ptrArray[face_in][I + 1][J + 1][K + 1]->Pos3().x()) /
                                     8.0;
                        double posy = (ptrArray[face_in][I][J][K]->Pos3().y() +
                                      ptrArray[face_in][I + 1][J][K]->Pos3().y() +
                                      ptrArray[face_in][I][J + 1][K]->Pos3().y() +
                                      ptrArray[face_in][I + 1][J + 1][K]->Pos3().y() +
                                      ptrArray[face_in][I][J][K + 1]->Pos3().y() +
                                      ptrArray[face_in][I + 1][J][K + 1]->Pos3().y() +
                                      ptrArray[face_in][I][J + 1][K + 1]->Pos3().y() +
                                      ptrArray[face_in][I + 1][J + 1][K + 1]->Pos3().y()) /
                                     8.0;

                        double posz = (ptrArray[face_in][I][J][K]->Pos3().z() +
                                      ptrArray[face_in][I + 1][J][K]->Pos3().z() +
                                      ptrArray[face_in][I][J + 1][K]->Pos3().z() +
                                      ptrArray[face_in][I + 1][J + 1][K]->Pos3().z() +
                                      ptrArray[face_in][I][J][K + 1]->Pos3().z() +
                                      ptrArray[face_in][I + 1][J][K + 1]->Pos3().z() +
                                      ptrArray[face_in][I][J + 1][K + 1]->Pos3().z() +
                                      ptrArray[face_in][I + 1][J + 1][K + 1]->Pos3().z()) /
                                     8.0;
                        Vector3 tempPos = Vector3(posx, posy, posz);
                        ptrposVectorCellArray[face_in][I][J][K]->SetVector3(tempPos);
                    }
                }
            }
        }
    }

    // print out divB
    using namespace H5;
    char filen[80], filename[151];

    sprintf(filen, "DivBinCells%d.h5", timeline);
    strncpy(filename, outpdir, 150);
    strncat(filename, filen, 150);

    H5std_string FILE_NAME(filename);
    H5std_string DATASET_NAME(filename);
    H5std_string DATASET_CONST_NAME("ArrayOfCellsPos_const");

    H5std_string MEMBERx("x");
    H5std_string MEMBERy("y");
    H5std_string MEMBERz("z");

    const int RANK = 4;

    typedef struct //Vector3_h5
    {
        double v_x;
        double v_y;
        double v_z;
    } Vector3_h5;

    Vector3_h5 *data_mem_pos = new Vector3_h5[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize];
    Vector3_h5 ****array_data_pos = new Vector3_h5 ***[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        array_data_pos[face] = new Vector3_h5 **[fieldsGridsSize];
        for (i = 0; i < fieldsGridsSize; i++)
        {
            array_data_pos[face][i] = new Vector3_h5 *[fieldsGridsSize];
            for (j = 0; j < fieldsGridsSize; j++)
            {
                array_data_pos[face][i][j] = data_mem_pos + face * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize + i * fieldsGridsSize * fieldsGridsSize + j * fieldsGridsSize;
                for (k = 0; k < fieldsGridsSize; k++)
                {
                    array_data_pos[face][i][j][k].v_x = ptrposVectorCellArray[face][i + 1][j + 1][k]->x();
                    array_data_pos[face][i][j][k].v_y = ptrposVectorCellArray[face][i + 1][j + 1][k]->y();
                    array_data_pos[face][i][j][k].v_z = ptrposVectorCellArray[face][i + 1][j + 1][k]->z();
                }
            }
        }
    }

    Vector3_h5 *data_mem_divB = new Vector3_h5[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize];
    Vector3_h5 ****array_data_divB = new Vector3_h5 ***[totalFace];

    for (int face = 0; face < totalFace; face++)
    {
        array_data_divB[face] = new Vector3_h5 **[fieldsGridsSize];
        for (i = 0; i < fieldsGridsSize; i++)
        {
            array_data_divB[face][i] = new Vector3_h5 *[fieldsGridsSize];
            for (j = 0; j < fieldsGridsSize; j++)
            {
                array_data_divB[face][i][j] = data_mem_divB + face * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize + i * fieldsGridsSize * fieldsGridsSize + j * fieldsGridsSize;
                for (k = 0; k < fieldsGridsSize; k++)
                {
                    array_data_divB[face][i][j][k].v_x = ptrDivBVectorCellArray[face][i + 1][j + 1][k]->x();
                    array_data_divB[face][i][j][k].v_y = ptrDivBVectorCellArray[face][i + 1][j + 1][k]->y();
                    array_data_divB[face][i][j][k].v_z = ptrDivBVectorCellArray[face][i + 1][j + 1][k]->z();
                }
            }
        }
    }
    Exception::dontPrint();

    hsize_t dim[] = {(hsize_t)totalFace, (hsize_t)fieldsGridsSize, (hsize_t)fieldsGridsSize, 
                     (hsize_t)fieldsGridsSize};
    DataSpace space(RANK, dim);

    CompType mtype_vector3(sizeof(Vector3_h5));
    mtype_vector3.insertMember(MEMBERx, HOFFSET(Vector3_h5, v_x), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERy, HOFFSET(Vector3_h5, v_y), PredType::NATIVE_DOUBLE);
    mtype_vector3.insertMember(MEMBERz, HOFFSET(Vector3_h5, v_z), PredType::NATIVE_DOUBLE);

    H5File *file = new H5File(FILE_NAME, H5F_ACC_RDWR);

    DataSet *dataset;
    dataset = new DataSet(file->createDataSet(DATASET_NAME, mtype_vector3, space));
    dataset->write(array_data_divB[0][0][0], mtype_vector3);
    delete dataset;
    delete data_mem_divB;

    // print our pos Cell array
    if (timeline == 0)
    {
        DataSet *dataset_pos;
        dataset_pos = new DataSet(file->createDataSet(DATASET_CONST_NAME, mtype_vector3, space));
        dataset_pos->write(array_data_pos[0][0][0], mtype_vector3);
        delete dataset_pos;
    }

    delete file;
    delete data_mem_pos;
}

int ReadSavedData(GridsPoints *****ptrArray,
                  GridsCells ****ptrArrayCells)
{
    using namespace H5;
    H5std_string FILE_NAME("./Data/DataSaved.h5");
    H5std_string DATASET_NAME_H("Particles_H"); // pos, v, weight, mu
    H5std_string DATASET_NAME_He("Particles_He");
    H5std_string DATASET_NAME_O("Particles_O");
    H5std_string DATASET_NAME_H_Uint("Particles_H_Uint"); // Uint of pos
    H5std_string DATASET_NAME_He_Uint("Particles_He_Uint");
    H5std_string DATASET_NAME_O_Uint("Particles_O_Uint");

    H5std_string DATASET_NAME_Grids("Grids");

    H5std_string DATASET_NAME_Time("Time");

    H5File file_saved(FILE_NAME, H5F_ACC_RDWR);

    hsize_t dims_read_H[2]={0,0};
    hsize_t dims_read_He[2]={0,0};
    hsize_t dims_read_O[2]={0,0};

    DataSet dataset_read_H;
    DataSet dataset_read_He;
    DataSet dataset_read_O;
    dataset_read_H = file_saved.openDataSet(DATASET_NAME_H);
    dataset_read_He = file_saved.openDataSet(DATASET_NAME_He);
    dataset_read_O = file_saved.openDataSet(DATASET_NAME_O);

    DataSpace dataspace_read_H = dataset_read_H.getSpace();
    DataSpace dataspace_read_He = dataset_read_He.getSpace();
    DataSpace dataspace_read_O = dataset_read_O.getSpace();

    dataspace_read_H.getSimpleExtentDims(dims_read_H, NULL);
    int numParticle_H = dims_read_H[0]; //ptrParticlesList_H.size();

    dataspace_read_He.getSimpleExtentDims(dims_read_He, NULL);
    int numParticle_He = dims_read_He[0]; //ptrParticlesList_H.size();

    dataspace_read_O.getSimpleExtentDims(dims_read_O, NULL);
    int numParticle_O = dims_read_O[0]; //ptrParticlesList_H.size();

    //hsize_t dims_read_H_Uint[2];
    //hsize_t dims_read_He_Uint[2];
    //hsize_t dims_read_O_Uint[2];

    DataSet dataset_read_H_Uint;
    DataSet dataset_read_He_Uint;
    DataSet dataset_read_O_Uint;
    dataset_read_H_Uint = file_saved.openDataSet(DATASET_NAME_H_Uint);
    dataset_read_He_Uint = file_saved.openDataSet(DATASET_NAME_He_Uint);
    dataset_read_O_Uint = file_saved.openDataSet(DATASET_NAME_O_Uint);
    DataSpace dataspace_read_H_Uint = dataset_read_H_Uint.getSpace();
    DataSpace dataspace_read_He_Uint = dataset_read_He_Uint.getSpace();
    DataSpace dataspace_read_O_Uint = dataset_read_O_Uint.getSpace();

    std::cout << " test \n"
              << numParticle_H << " " << numParticle_He << " " << numParticle_O << "\n";

    int par_size;

    if (numParticle_H >= numParticle_He)
    {
        par_size = numParticle_H;
    }
    else
    {
        par_size = numParticle_He;
    }
    if (par_size <= numParticle_O)
    {
        par_size = numParticle_O;
    }

    // apply continus memory
    double *data_mem_particles = new double[8 * par_size];
    double **array_data_particles = new double *[par_size];
    unsigned long long *array_data_particles_Uint = new unsigned long long[par_size];
    //    unsigned long long array_data_particles_Uint[par_size];

    for (int i = 0; i < numParticle_H; i++)
    {
        array_data_particles[i] = data_mem_particles + 8 * i;
    }

    // read particles H from HDF
    dataset_read_H.read(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataset_read_H_Uint.read(array_data_particles_Uint, PredType::NATIVE_ULLONG);
// insert into particles H
#pragma omp parallel for
    for (int i = 0; i < numParticle_H; i++)
    {
        if (array_data_particles_Uint[i] == 0)
            continue;
        else
        {
            Particles tempP;
            tempP.SetParticles(array_data_particles_Uint[i],
                               array_data_particles[i][0],
                               array_data_particles[i][1],
                               array_data_particles[i][2],
                               array_data_particles[i][3],
                               array_data_particles[i][4],
                               array_data_particles[i][5],
                               array_data_particles[i][6],
                               array_data_particles[i][7]);
            auto tempStr = tempP.InttoStrp1();
#pragma omp critical
            ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_H()->push_back(tempP);
        }
    }
    //
    for (int i = 0; i < numParticle_He; i++)
    {
        array_data_particles[i] = data_mem_particles + 8 * i;
    }
    // read particles He from HDF
    dataset_read_He.read(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataset_read_He_Uint.read(array_data_particles_Uint, PredType::NATIVE_ULLONG);
// insert into particles He
#pragma omp parallel for
    for (int i = 0; i < numParticle_He; i++)
    {
        if (array_data_particles_Uint[i] == 0)
            continue;
        else
        {
            Particles tempP;
            tempP.SetParticles(array_data_particles_Uint[i],
                               array_data_particles[i][0],
                               array_data_particles[i][1],
                               array_data_particles[i][2],
                               array_data_particles[i][3],
                               array_data_particles[i][4],
                               array_data_particles[i][5],
                               array_data_particles[i][6],
                               array_data_particles[i][7]);
            auto tempStr = tempP.InttoStrp1();
#pragma omp critical
            ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_He()->push_back(tempP);
        }
    }
    //
    for (int i = 0; i < numParticle_O; i++)
    {
        array_data_particles[i] = data_mem_particles + 8 * i;
    }
    //read particles O from HDF
    dataset_read_O.read(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataset_read_O_Uint.read(array_data_particles_Uint, PredType::NATIVE_ULLONG);
// insert into particles O
#pragma omp parallel for
    for (int i = 0; i < numParticle_O; i++)
    {
        if (array_data_particles_Uint[i] == 0)
            continue;
        else
        {
            Particles tempP;
            tempP.SetParticles(array_data_particles_Uint[i],
                               array_data_particles[i][0],
                               array_data_particles[i][1],
                               array_data_particles[i][2],
                               array_data_particles[i][3],
                               array_data_particles[i][4],
                               array_data_particles[i][5],
                               array_data_particles[i][6],
                               array_data_particles[i][7]);
            auto tempStr = tempP.InttoStrp1();
#pragma omp critical
            ptrArrayCells[tempStr.face][tempStr.ig][tempStr.jg][tempStr.kg].Particles_O()->push_back(tempP);
        }
    }
    //
    delete data_mem_particles;
    delete array_data_particles_Uint;
    delete array_data_particles;

    //*********************************************************************************
    // Grids reading
    //hsize_t dims_read_grids[2];
    DataSet dataset_read_grids;
    dataset_read_grids = file_saved.openDataSet(DATASET_NAME_Grids);
    DataSpace dataspace_read_grids = dataset_read_grids.getSpace();
    double *data_mem_grids = new double[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52];
    double *****data_grids = new double ****[totalFace];
    //
    for (int f = 0; f < totalFace; f++)
    {
        data_grids[f] = new double ***[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            data_grids[f][i] = new double **[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                data_grids[f][i][j] = new double *[1 + fieldsGridsSize * grid_domain];
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    data_grids[f][i][j][k] = data_mem_grids + f * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             j * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             k * 57;
                }
            }
        }
    }
    //
    dataset_read_grids.read(data_grids[0][0][0][0], PredType::NATIVE_DOUBLE);
    //
    dataset_read_grids.close();
    dataspace_read_grids.close();
    //
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    ptrArray[f][i + 1][j + 1][k]->ReadGridsPoints(data_grids[f][i][j][k][0], data_grids[f][i][j][k][1], data_grids[f][i][j][k][2],
                                                                  data_grids[f][i][j][k][3], data_grids[f][i][j][k][4], data_grids[f][i][j][k][5],
                                                                  data_grids[f][i][j][k][6], data_grids[f][i][j][k][7], data_grids[f][i][j][k][8],
                                                                  data_grids[f][i][j][k][9], data_grids[f][i][j][k][10], data_grids[f][i][j][k][11],
                                                                  data_grids[f][i][j][k][12], data_grids[f][i][j][k][13], data_grids[f][i][j][k][14],
                                                                  data_grids[f][i][j][k][15], data_grids[f][i][j][k][16], data_grids[f][i][j][k][17],
                                                                  data_grids[f][i][j][k][18], data_grids[f][i][j][k][19], data_grids[f][i][j][k][20],
                                                                  data_grids[f][i][j][k][21], data_grids[f][i][j][k][22], data_grids[f][i][j][k][23],
                                                                  data_grids[f][i][j][k][24], data_grids[f][i][j][k][25], data_grids[f][i][j][k][26],
                                                                  data_grids[f][i][j][k][27], data_grids[f][i][j][k][28], data_grids[f][i][j][k][29],
                                                                  data_grids[f][i][j][k][30], data_grids[f][i][j][k][31], data_grids[f][i][j][k][32],
                                                                  data_grids[f][i][j][k][33], data_grids[f][i][j][k][34], data_grids[f][i][j][k][35],
                                                                  data_grids[f][i][j][k][36], data_grids[f][i][j][k][37], data_grids[f][i][j][k][38],
                                                                  data_grids[f][i][j][k][39], data_grids[f][i][j][k][40], data_grids[f][i][j][k][41],
                                                                  data_grids[f][i][j][k][42], data_grids[f][i][j][k][43], data_grids[f][i][j][k][44],
                                                                  data_grids[f][i][j][k][45], data_grids[f][i][j][k][46], data_grids[f][i][j][k][47],
                                                                  data_grids[f][i][j][k][48], data_grids[f][i][j][k][49], data_grids[f][i][j][k][50],
                                                                  data_grids[f][i][j][k][51], data_grids[f][i][j][k][52], data_grids[f][i][j][k][53],
                                                                  data_grids[f][i][j][k][54], data_grids[f][i][j][k][55], data_grids[f][i][j][k][56],
                                                                  0);
                }
            }
        }
    }
    //
    delete data_grids;
    delete data_mem_grids;
    //*********************************************************************************
    // Time reading
    int data_time;
    DataSet dataset_read_time;
    dataset_read_time = file_saved.openDataSet(DATASET_NAME_Time);
    dataset_read_time.read(&data_time, PredType::NATIVE_DOUBLE);
    //
    dataset_read_time.close();
    dataset_read_time.close();
    //
    file_saved.close();
    //
    return data_time;
}

//************************************************************************
//************************************************************************
// Read particles and grids
int ReadSavedData(GridsPoints *****ptrArray,
                  vector<Particles> &ptrParticlesList_H,
                  vector<Particles> &ptrParticlesList_He,
                  vector<Particles> &ptrParticlesList_O,
                  vector<int> &ptrParticlesList_out_H,
                  vector<int> &ptrParticlesList_out_He,
                  vector<int> &ptrParticlesList_out_O)
{
    using namespace H5;
    H5std_string FILE_NAME("./Data/DataSaved.h5");
    H5std_string DATASET_NAME_H("Particles_H"); // pos, v, weight, mu
    H5std_string DATASET_NAME_He("Particles_He");
    H5std_string DATASET_NAME_O("Particles_O");
    H5std_string DATASET_NAME_H_Uint("Particles_H_Uint"); // Uint of pos
    H5std_string DATASET_NAME_He_Uint("Particles_He_Uint");
    H5std_string DATASET_NAME_O_Uint("Particles_O_Uint");

    H5std_string DATASET_NAME_Grids("Grids");

    H5std_string DATASET_NAME_Time("Time");

    H5File file_saved(FILE_NAME, H5F_ACC_RDWR);

    hsize_t dims_read_H[2]={0,0};
    hsize_t dims_read_He[2]={0,0};
    hsize_t dims_read_O[2]={0,0};

    DataSet dataset_read_H;
    DataSet dataset_read_He;
    DataSet dataset_read_O;
    dataset_read_H = file_saved.openDataSet(DATASET_NAME_H);
    dataset_read_He = file_saved.openDataSet(DATASET_NAME_He);
    dataset_read_O = file_saved.openDataSet(DATASET_NAME_O);

    DataSpace dataspace_read_H = dataset_read_H.getSpace();
    DataSpace dataspace_read_He = dataset_read_He.getSpace();
    DataSpace dataspace_read_O = dataset_read_O.getSpace();

    dataspace_read_H.getSimpleExtentDims(dims_read_H, NULL);
    int numParticle_H = dims_read_H[0]; //ptrParticlesList_H.size();

    dataspace_read_He.getSimpleExtentDims(dims_read_He, NULL);
    int numParticle_He = dims_read_He[0]; //ptrParticlesList_H.size();

    dataspace_read_O.getSimpleExtentDims(dims_read_O, NULL);
    int numParticle_O = dims_read_O[0]; //ptrParticlesList_H.size();

    //hsize_t dims_read_H_Uint[2];
    //hsize_t dims_read_He_Uint[2];
    //hsize_t dims_read_O_Uint[2];

    DataSet dataset_read_H_Uint;
    DataSet dataset_read_He_Uint;
    DataSet dataset_read_O_Uint;
    dataset_read_H_Uint = file_saved.openDataSet(DATASET_NAME_H_Uint);
    dataset_read_He_Uint = file_saved.openDataSet(DATASET_NAME_He_Uint);
    dataset_read_O_Uint = file_saved.openDataSet(DATASET_NAME_O_Uint);
    DataSpace dataspace_read_H_Uint = dataset_read_H_Uint.getSpace();
    DataSpace dataspace_read_He_Uint = dataset_read_He_Uint.getSpace();
    DataSpace dataspace_read_O_Uint = dataset_read_O_Uint.getSpace();

    std::cout << " test \n"
              << numParticle_H << " " << numParticle_He << " " << numParticle_O << "\n";

    int par_size;

    if (numParticle_H >= numParticle_He)
    {
        par_size = numParticle_H;
    }
    else
    {
        par_size = numParticle_He;
    }
    if (par_size <= numParticle_O)
    {
        par_size = numParticle_O;
    }

    // apply continus memory
    double *data_mem_particles = new double[8 * par_size];
    double **array_data_particles = new double *[par_size];
    unsigned long long *array_data_particles_Uint = new unsigned long long[par_size];
    //    unsigned long long array_data_particles_Uint[par_size];

    for (int i = 0; i < numParticle_H; i++)
    {
        array_data_particles[i] = data_mem_particles + 8 * i;
    }

    // read particles H from HDF
    dataset_read_H.read(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataset_read_H_Uint.read(array_data_particles_Uint, PredType::NATIVE_ULLONG);
    // insert into particles H

#pragma omp parallel for
    for (int i = 0; i < numParticle_H; i++)
    {
        if (array_data_particles_Uint[i] == 0)
            continue;
        else
        {
            Particles tempP;
            tempP.SetParticles(array_data_particles_Uint[i],
                               array_data_particles[i][0],
                               array_data_particles[i][1],
                               array_data_particles[i][2],
                               array_data_particles[i][3],
                               array_data_particles[i][4],
                               array_data_particles[i][5],
                               array_data_particles[i][6],
                               array_data_particles[i][7]);
#pragma omp critical
            ptrParticlesList_H.push_back(tempP);
        }
    }

    for (int i = 0; i < numParticle_He; i++)
    {
        array_data_particles[i] = data_mem_particles + 8 * i;
    }
    // read particles He from HDF
    dataset_read_He.read(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataset_read_He_Uint.read(array_data_particles_Uint, PredType::NATIVE_ULLONG);
// insert into particles He
#pragma omp parallel for
    for (int i = 0; i < numParticle_He; i++)
    {
        if (array_data_particles_Uint[i] == 0)
            continue;
        else
        {
            Particles tempP;
            tempP.SetParticles(array_data_particles_Uint[i],
                               array_data_particles[i][0],
                               array_data_particles[i][1],
                               array_data_particles[i][2],
                               array_data_particles[i][3],
                               array_data_particles[i][4],
                               array_data_particles[i][5],
                               array_data_particles[i][6],
                               array_data_particles[i][7]);
#pragma omp critical
            ptrParticlesList_He.push_back(tempP);
        }
    }

    for (int i = 0; i < numParticle_O; i++)
    {
        array_data_particles[i] = data_mem_particles + 8 * i;
    }
    //read particles O from HDF
    dataset_read_O.read(array_data_particles[0], PredType::NATIVE_DOUBLE);
    dataset_read_O_Uint.read(array_data_particles_Uint, PredType::NATIVE_ULLONG);
// insert into particles O
#pragma omp parallel for
    for (int i = 0; i < numParticle_O; i++)
    {
        if (array_data_particles_Uint[i] == 0)
            continue;
        else
        {
            Particles tempP;
            tempP.SetParticles(array_data_particles_Uint[i],
                               array_data_particles[i][0],
                               array_data_particles[i][1],
                               array_data_particles[i][2],
                               array_data_particles[i][3],
                               array_data_particles[i][4],
                               array_data_particles[i][5],
                               array_data_particles[i][6],
                               array_data_particles[i][7]);
#pragma omp critical
            ptrParticlesList_O.push_back(tempP);
        }
    }

    delete data_mem_particles;
    delete array_data_particles_Uint;
    delete array_data_particles;

    //*********************************************************************************
    // Grids reading
    //hsize_t dims_read_grids[2];

    DataSet dataset_read_grids;
    dataset_read_grids = file_saved.openDataSet(DATASET_NAME_Grids);
    DataSpace dataspace_read_grids = dataset_read_grids.getSpace();

    double *data_mem_grids = new double[totalFace * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52];
    double *****data_grids = new double ****[totalFace];

    for (int f = 0; f < totalFace; f++)
    {
        data_grids[f] = new double ***[1 + fieldsGridsSize];
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            data_grids[f][i] = new double **[1 + fieldsGridsSize];
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                data_grids[f][i][j] = new double *[1 + fieldsGridsSize * grid_domain];
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    data_grids[f][i][j][k] = data_mem_grids + f * (1 + fieldsGridsSize) * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             i * (1 + fieldsGridsSize) * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             j * (1 + fieldsGridsSize * grid_domain) * 52 +
                                             k * 57;
                }
            }
        }
    }

    dataset_read_grids.read(data_grids[0][0][0][0], PredType::NATIVE_DOUBLE);

    dataset_read_grids.close();
    dataspace_read_grids.close();

    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < 1 + fieldsGridsSize; i++)
        {
            for (int j = 0; j < 1 + fieldsGridsSize; j++)
            {
                for (int k = 0; k < 1 + fieldsGridsSize * grid_domain; k++)
                {
                    ptrArray[f][i + 1][j + 1][k]->ReadGridsPoints(data_grids[f][i][j][k][0], data_grids[f][i][j][k][1], data_grids[f][i][j][k][2],
                                                                  data_grids[f][i][j][k][3], data_grids[f][i][j][k][4], data_grids[f][i][j][k][5],
                                                                  data_grids[f][i][j][k][6], data_grids[f][i][j][k][7], data_grids[f][i][j][k][8],
                                                                  data_grids[f][i][j][k][9], data_grids[f][i][j][k][10], data_grids[f][i][j][k][11],
                                                                  data_grids[f][i][j][k][12], data_grids[f][i][j][k][13], data_grids[f][i][j][k][14],
                                                                  data_grids[f][i][j][k][15], data_grids[f][i][j][k][16], data_grids[f][i][j][k][17],
                                                                  data_grids[f][i][j][k][18], data_grids[f][i][j][k][19], data_grids[f][i][j][k][20],
                                                                  data_grids[f][i][j][k][21], data_grids[f][i][j][k][22], data_grids[f][i][j][k][23],
                                                                  data_grids[f][i][j][k][24], data_grids[f][i][j][k][25], data_grids[f][i][j][k][26],
                                                                  data_grids[f][i][j][k][27], data_grids[f][i][j][k][28], data_grids[f][i][j][k][29],
                                                                  data_grids[f][i][j][k][30], data_grids[f][i][j][k][31], data_grids[f][i][j][k][32],
                                                                  data_grids[f][i][j][k][33], data_grids[f][i][j][k][34], data_grids[f][i][j][k][35],
                                                                  data_grids[f][i][j][k][36], data_grids[f][i][j][k][37], data_grids[f][i][j][k][38],
                                                                  data_grids[f][i][j][k][39], data_grids[f][i][j][k][40], data_grids[f][i][j][k][41],
                                                                  data_grids[f][i][j][k][42], data_grids[f][i][j][k][43], data_grids[f][i][j][k][44],
                                                                  data_grids[f][i][j][k][45], data_grids[f][i][j][k][46], data_grids[f][i][j][k][47],
                                                                  data_grids[f][i][j][k][48], data_grids[f][i][j][k][49], data_grids[f][i][j][k][50],
                                                                  data_grids[f][i][j][k][51], data_grids[f][i][j][k][52], data_grids[f][i][j][k][53],
                                                                  data_grids[f][i][j][k][54], data_grids[f][i][j][k][55], data_grids[f][i][j][k][56],
                                                                  0);
                }
            }
        }
    }

    delete data_grids;
    delete data_mem_grids;

    //*********************************************************************************
    // Time reading
    int data_time;

    DataSet dataset_read_time;
    dataset_read_time = file_saved.openDataSet(DATASET_NAME_Time);

    dataset_read_time.read(&data_time, PredType::NATIVE_DOUBLE);

    dataset_read_time.close();
    dataset_read_time.close();

    file_saved.close();

    return data_time;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density
// at each grids points
// The volume means the sum volume of adjacent 6 or 8 cells
// Only consider the main domain so that the size is (fieldsGridsSize+1)*(fieldsGridsSize+1)*(fieldsGridsSize+1)
double ***VolumeGridsField(double ***ptrVolumeCellArray_in)
{
    static double *mem_VolumeGridsArray = new double[(fieldsGridsSize + 1) * (fieldsGridsSize + 1) * (fieldsGridsSize * grid_domain + 1)];
    double ***VolumeGridsArray = new double **[fieldsGridsSize + 1];

    for (int i = 0; i < fieldsGridsSize + 1; i++)
    {
        VolumeGridsArray[i] = new double *[fieldsGridsSize + 1];
        for (int j = 0; j < fieldsGridsSize + 1; j++)
        {
            //            VolumeGridsArray[i][j]= new double [fieldsGridsSize+1];
            VolumeGridsArray[i][j] = mem_VolumeGridsArray + i * (fieldsGridsSize + 1) * (fieldsGridsSize * grid_domain + 1) + j * (fieldsGridsSize * grid_domain + 1);
            for (int k = 0; k < fieldsGridsSize * grid_domain + 1; k++)
            {
                if (k == 0 || k == fieldsGridsSize * grid_domain)
                {
                    VolumeGridsArray[i][j][k] = 999.99; // should not be used
                }
                else
                {
                    if (i == 0 && j == 0)
                    {
                        VolumeGridsArray[i][j][k] = (3 * ptrVolumeCellArray_in[i + 1][j + 1][k] + 3 * ptrVolumeCellArray_in[i + 1][j + 1][k + 1]) / 6.0;
                    }

                    else if (i == 0 && j == fieldsGridsSize)
                    {
                        VolumeGridsArray[i][j][k] = (3 * ptrVolumeCellArray_in[i + 1][j][k] + 3 * ptrVolumeCellArray_in[i + 1][j][k + 1]) / 6.0;
                    }

                    else if (i == fieldsGridsSize && j == 0)
                    {
                        VolumeGridsArray[i][j][k] = (3 * ptrVolumeCellArray_in[i][j + 1][k] + 3 * ptrVolumeCellArray_in[i][j + 1][k + 1]) / 6.0;
                    }

                    else if (i == fieldsGridsSize && j == fieldsGridsSize)
                    {
                        VolumeGridsArray[i][j][k] = (3 * ptrVolumeCellArray_in[i][j][k] + 3 * ptrVolumeCellArray_in[i][j][k + 1]) / 6.0;
                    }

                    else if (i == 0 && j != 0 && j != fieldsGridsSize)
                    {
                        VolumeGridsArray[i][j][k] = (2 * ptrVolumeCellArray_in[i + 1][j + 1][k] + 2 * ptrVolumeCellArray_in[i + 1][j][k] + 2 * ptrVolumeCellArray_in[i + 1][j + 1][k + 1] + 2 * ptrVolumeCellArray_in[i + 1][j][k + 1] ) / 8.0;
                    }

                    else if (i == fieldsGridsSize && j != 0 && j != fieldsGridsSize)
                    {
                        VolumeGridsArray[i][j][k] = (2 * ptrVolumeCellArray_in[i][j + 1][k] + 2 * ptrVolumeCellArray_in[i][j][k] + 2 * ptrVolumeCellArray_in[i][j + 1][k + 1] + 2 * ptrVolumeCellArray_in[i][j][k + 1]) * 0.125;
                    }

                    else if (j == 0 && i != 0 && i != fieldsGridsSize)
                    {
                        VolumeGridsArray[i][j][k] = (2 * ptrVolumeCellArray_in[i + 1][j + 1][k] + 2 * ptrVolumeCellArray_in[i][j + 1][k] + 2 * ptrVolumeCellArray_in[i + 1][j + 1][k + 1] + 2 * ptrVolumeCellArray_in[i][j + 1][k]) * 0.125;
                    }

                    else if (j == fieldsGridsSize && i != 0 && i != fieldsGridsSize)
                    {
                        VolumeGridsArray[i][j][k] = (2 * ptrVolumeCellArray_in[i + 1][j][k] + 2 * ptrVolumeCellArray_in[i][j][k] + 2 * ptrVolumeCellArray_in[i + 1][j][k + 1] + 2 * ptrVolumeCellArray_in[i][j][k + 1]) * 0.125;
                    }

                    else
                    {
                        VolumeGridsArray[i][j][k] = (ptrVolumeCellArray_in[i + 1][j + 1][k] + ptrVolumeCellArray_in[i][j + 1][k] + ptrVolumeCellArray_in[i + 1][j][k] + ptrVolumeCellArray_in[i][j][k] + ptrVolumeCellArray_in[i + 1][j + 1][k + 1] + ptrVolumeCellArray_in[i][j + 1][k + 1] + ptrVolumeCellArray_in[i + 1][j][k + 1] + ptrVolumeCellArray_in[i][j][k + 1]) *0.125;
                    }
                }
            }
        }
    }

    return VolumeGridsArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume of cells needed to calculate density
// at each cell [fieldsize+2][fieldsize+2][fieldsize]
double ***VolumeCellsField(GridsPoints *****ptrArray_in)
{
    static double *mem_VolumeCellsArray = new double[(fieldsGridsSize + 2) * (fieldsGridsSize + 2) * (fieldsGridsSize * grid_domain)];
    double ***VolumeCellsArray = new double **[fieldsGridsSize + 2];
    for (int i = 0; i < fieldsGridsSize + 2; i++)
    {
        VolumeCellsArray[i] = new double *[fieldsGridsSize + 2];
        for (int j = 0; j < fieldsGridsSize + 2; j++)
        {
            //            VolumeCellsArray[i][j] = new double[fieldsGridsSize];
            VolumeCellsArray[i][j] = mem_VolumeCellsArray + i * (fieldsGridsSize + 2) * (fieldsGridsSize * grid_domain) + j * (fieldsGridsSize * grid_domain);
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                VolumeCellsArray[i][j][k] = CellVolume(ptrArray_in, 0, i, j, k);
            }
        }
    }
    return VolumeCellsArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Create a array to store volume weight for weighting calculation
// at each grids [fieldsize+1][fieldsize+1][fieldsize+1]
double ***VolumeWeightGridsField(GridsPoints *****ptrArray_in)
{
    static double *mem_VolumeWeightGridsArray = new double[(fieldsGridsSize + 1) * (fieldsGridsSize + 1) * (fieldsGridsSize * grid_domain + 1)];
    double ***VolumeWeightGridsArray = new double **[fieldsGridsSize + 1];
    for (int i = 0; i < fieldsGridsSize + 1; i++)
    {
        VolumeWeightGridsArray[i] = new double *[fieldsGridsSize + 1];
        for (int j = 0; j < fieldsGridsSize + 1; j++)
        {
            VolumeWeightGridsArray[i][j] = mem_VolumeWeightGridsArray + i * (fieldsGridsSize + 1) * (fieldsGridsSize * grid_domain + 1) + j * (fieldsGridsSize * grid_domain + 1);
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                if (i == 0 || j == 0)
                    continue;
                double x = ptrArray_in[0][i][j][k]->Pos3().x();
                double y = ptrArray_in[0][i][j][k]->Pos3().y();
                double z = ptrArray_in[0][i][j][k]->Pos3().z();
                VolumeWeightGridsArray[i][j][k] = WeightCalculation(x, y, z);
            }
        }
    }
    return VolumeWeightGridsArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Calculation the temprature
void Titheridge_Te(GridsPoints *****ptrArray_in, int DOY)
//  by Yuzhang Ma AUG 16,2019. JTU 3/4/2022: Revised SEC, DOY, Kp, and F107 specification
{
    //inputs
    double SEC = ihour*3600+imin*60+isec; // second of day
    //double DOY = 31;        // day of year

    //double Kp = 10;
    //double F107 = 100;

    //  consts
    const double rad2deg = 57.2957763671875;
    //const double pi = 3.1415927410125732;
    const double PLAT = 1.37846; //geographic latitude (rad) of Earth's north magnetic pole
    const double PLON = 5.04557; //geographic longitude (rad) of Earth's north magnetic pole
    const double h0 = 400.0e0, R0 = 1.0627825e0, R02 = 1.1295067e0;
    const double por = 0.2857142e0, Re = 6371.2e0, x1 = 1.0470869e0;
    const double a0 = 1.23e3, a1 = 2.2e3, a2 = -3.29e3, a3 = -2.6e-1, a4 = -6.8e-1;
    const double b0 = 4.65e0, b1 = -8.55e0, b2 = 4.14e0, b3 = -2.16e0, b4 = 1.45e0;
    const double aa0 = 9.85e2, aa1 = -9.63e2, aa2 = 1.125e3, aa3 = -0.6e0, aa4 = 0.10e0;
    const double bb0 = 7.56e-1, bb1 = -8.8e-1, bb2 = 2.9e-1, bb3 = -2.63e0, bb4 = 1.84e0;

    double GLONR = 15.0 * (12.0 - SEC / 3600.0) / rad2deg; //GLONR in rad

    if (GLONR < 0.0)
    {
        GLONR += 2.0 * pi; // GLONR = GLONR + 2*pi
    }
    if (GLONR > 2.0 * pi)
    {
        GLONR -= 2.0 * pi; // GLONR = GLONR - 2*pi
    }

    double MLON12 = atan(sin(PLAT) * tan(GLONR - PLON));

    if (MLON12 < 0.0)
    {
        MLON12 += 2.0 * pi; // MLON12 = MLON12 + 2*pi
    }

    //
    //  caution: Px_in py_in pz_in should be transformed to SM Coordinate before using.
    //

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 0; k < fieldsGridsSize * grid_domain + 1; k++)
                {

                    //check stopsign
                    if (ptrArray_in[face][i][j][k]->StopSign() == 1)
                        continue;

                    double px_in = ptrArray_in[face][i][j][k]->Pos3().x();
                    double py_in = ptrArray_in[face][i][j][k]->Pos3().y();
                    double pz_in = ptrArray_in[face][i][j][k]->Pos3().z();

                    //    std::cout << face << i<<j<< k<< " --> ";
                    //    std::cout << px_in << " " << py_in << " " << pz_in << " A--> ";

                    double rr = ptrArray_in[face][i][j][k]->Pos3().norm(); // in unit of m
                    double Lat = asin(pz_in / rr);
                    double phi = atan(py_in / px_in);
                    double L = rr / pow(cos(Lat), 2.0);
                    double heq = (L - 1.0e0) * Re;
                    double SL = 1.0e0 - x1 / L;
                    double BLATD = asin(sqrt(SL)) * rad2deg; // BLATD in degree geomagnetic latitude at 300km
                    double BLOND = (phi + MLON12) * rad2deg; // BLOND in degree // NaN for polar points

                    // std::cout << phi << " " << BLATD << " " << BLOND << " B==> ";

                    //GLAT calculate
                    double BLONR = BLOND / rad2deg;
                    double BLATR = BLATD / rad2deg;
                    double XM = cos(BLATR) * cos(BLONR);
                    double ZM = sin(BLATR);
                    //double XG = XM * sin(PLAT) + ZM * cos(PLAT);
                    double ZG = -XM * cos(PLAT) + ZM * sin(PLAT);
                    double GLATD = asin(ZG) * rad2deg; // GLATD in degree
                    double GLATR = GLATD / rad2deg;
                    //  day time T0 and G0  Eq.(19) in Titheridge, JGR1998
                    double T0d = (a0 + a1 * SL + a2 * SL * SL) / (1.0e0 + a3 * SL + a4 * SL * SL);
                    double G0d = (b0 + b1 * SL + b2 * SL * SL) / (1.0e0 + b3 * SL + b4 * SL * SL);

                    //  night time T0 and G0
                    double T0n = (aa0 + aa1 * SL + aa2 * SL * SL) / (1.0e0 + aa3 * SL + aa4 * SL * SL);
                    double G0n = (bb0 + bb1 * SL + bb2 * SL * SL) / (1.0e0 + bb3 * SL + bb4 * SL * SL);

                    //  change in solar activity
                    double delta_T0 = 3.4 * (F107 - 120.0);
                    G0d = G0d * T0d / (T0d + delta_T0);
                    T0d = T0d + delta_T0;
                    G0n = G0n * T0n / (T0n + delta_T0);
                    T0n = T0n + delta_T0;
                    //  change in magnetic activity
                    if (BLATD > 46.0e0)
                    {
                        double z = 0.135 * (BLATD - 46.0) / rad2deg;
                        double Dk = 37.0 + 1.33 * (BLATD - 46.0) - 37.0 * cos(z);
                        delta_T0 = Dk * (Kp - 1.5);
                        G0d = G0d * pow(T0d / (T0d + delta_T0), 2.5);
                        T0d = T0d + delta_T0;
                        G0n = G0n * pow(T0n / (T0n + delta_T0), 2.5);
                        T0n = T0n + delta_T0;
                    }

                    double G0T0d = G0d / T0d;
                    double G0T0n = G0n / T0n;

                    double alt = (rr / 1000 / Re - 1) * Re;
                    double alg = log(alt / h0);
                    double R2 = pow(1.0e0 + alt / Re, 2.0);
                    double Bh = 0.05e0 / (2.0e0 * L - R0) * (88.0e0 + alg * (10.5e0 - alg));

                    // std::cout << T0d << " " << Bh << " " << G0T0d << " " << G0T0n << " heq " << heq << " alt " << alt << " C==> ";
                    //  mean day values from Eq.(13) in Titheridge, 1998 JGR
                    double Tday0 = T0d * pow(1.0e0 + Bh * G0T0d * ((heq - h0) / R02 - (heq - alt) / R2), por);
                    //  mean night values from Eq.(13) in Titheridge, 1998 JGR
                    double Tnig0 = T0n * pow(1.0e0 + Bh * G0T0n * ((heq - h0) / R02 - (heq - alt) / R2), por);

                    //  std::cout << Tnig0 << " " << Tday0 << " D==> ";
                    //  local solar time (SAT) inputs: day of year,GLONR,local time
                    //  https://pvcdrom.pveducation.org/SUNLIGHT/SOLART.HTM
                    double B = 300 / 365 * (DOY - 81) / rad2deg; //B is degree; trans to rad
                    double EOT = 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B);
                    // get the time difference form GLONR
                    // https://blog.csdn.net/fct2001140269/article/details/86513925
                    int currentLon = GLONR * rad2deg;
                    if (currentLon > 180)
                    {
                        currentLon = -1 * (360 - currentLon);
                    }
                    int timeZone;
                    int shangValue = (int)(currentLon / 15);
                    double yushuValue = abs(currentLon % 15);
                    if (yushuValue <= 7.5)
                    {
                        timeZone = shangValue;
                    }
                    else
                    {
                        timeZone = shangValue + (currentLon > 0 ? 1 : -1);
                    }
                    double LSTM = 15 * timeZone;
                    double TC = 4 * (LSTM - GLONR * rad2deg) + EOT;
                    double SAT = SEC / 3600 + TC / 60;
                    if (SAT < 0)
                    {
                        SAT += 24.0;
                    }
                    if (SAT > 24.0)
                    {
                        SAT -= 24.0;
                    }
                    //DEC -- solar declination angle in radians
                    //solar declination angle inputs: day of year, output is in dgree, trans to rad)
                    //https://www.sciencedirect.com/topics/engineering/solar-declination
                    double DEC = 23.45 * sin(360 / 365 * (284 + DOY)) / rad2deg;
                    double S_LONG = -1.0e0 * tan(GLATR) * tan(DEC);
                    double s_r;
                    if (S_LONG >= -1.0e0 && S_LONG <= 1.0e0)
                    {
                        s_r = 12.0e0 - 3.82 * acos(S_LONG);
                    }
                    else
                    {
                        if (S_LONG <= -1.0e0)
                        {
                            s_r = 0.0e0;
                        }
                        else
                        {
                            s_r = 12.0e0;
                        }
                    }

                    //duration of the day-night transition
                    double Delta_tr = 1.2 + 0.5 * SL;
                    double Delta_ts = Delta_tr + 0.9;
                    //time of the center of the transition in hours after ground sunrise or sunset
                    double t_s = 0.15 * s_r;
                    if (t_s >= 1.0e0)
                    {
                        t_s = 1.0;
                    }
                    //double t_r = -0.5 * t_s;
                    double D;
                    if (SAT <= 12.0e0)
                    {
                        D = (s_r - 0.5 * t_s - SAT) / Delta_tr;
                    }
                    else
                    {
                        D = (SAT - t_s + s_r - 24.0e0) / Delta_ts;
                    }
                    phi = 1.0 / (1.0 + exp(3.2 * D));
                    double y = s_r - 12.5;
                    if (y <= -4.0)
                    {
                        y = -4.0;
                    }
                    double TDvar = 0.97 + 0.22 * pow((SAT - 12.5) / y, 2);
                    double Tday;
                    if (TDvar <= 1.2)
                    {
                        Tday = TDvar * Tday0;
                    }
                    else
                    {
                        Tday = 1.2 * Tday0;
                    }

                    if (SAT <= 12.0e0)
                    {
                        y = SAT;
                    }
                    else
                    {
                        y = SAT - 24.0;
                    }
                    double TNvar = 0.83 + 0.15 * exp(-0.2 * y * (1.4 - 0.8 * SL));
                    double Tnig;
                    if (TNvar <= 1.2)
                    {
                        Tnig = TNvar * Tnig0;
                    }
                    else
                    {
                        Tnig = 1.2 * Tnig0;
                    }

                    double Te = Tnig + (Tday - Tnig) * phi; // should be the avrage value of 2 Te points,it could be done out of the function

                    // std::cout << Te << std::endl;

                    ptrArray_in[face][i][j][k]->SetTemperature(Te);

                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(1);
                }
            }
        }
    }

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 0; k < fieldsGridsSize * grid_domain + 1; k++)
                {
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(0);
                }
            }
        }
    }
}


// //************************************************************************
// //************************************************************************
// // Function
// // Update PSD at gridpoints, Varney, 2016
// // Need s_parallel, alfvenic poynting flux
// // Update potential
// void PP_update_(GridsPoints *****ptrArray)
// {
//     double PI = 3.1415926535897;
//     double r_earth = radius + 100000.0; // 100 km earth shell altitude
//      // _top is values of
//     for (int f = 0; f < totalFace; f++)
//     {
//         for (int i = 1; i < fieldsGridsSize + 2; i++)
//         {
//             for (int j = 1; j < fieldsGridsSize + 2; j++)
//             {
//                 for (int k = 0; k < fieldsGridsSize * grid_domain + 1; k++)
//                 {
//                     double x = ptrArray[f][i][j][k]->Pos3().x();
//                     double y = ptrArray[f][i][j][k]->Pos3().y();
//                     double z = ptrArray[f][i][j][k]->Pos3().z();
//                     //
//                     double r_top = sqrt(x * x + y * y + z * z);
//                     double theta_top = acos(z / r_top);
//                     double phi_top;
//                     if (x != 0.0)
//                     {
//                         phi_top = atan(y / x);
//                         if (x < 0)
//                         {
//                             phi_top = phi_top + PI;
//                         }
//                     }
//                     else if (x == 0.0 && y > 0.0)
//                     {
//                         phi_top = PI / 2.0;
//                     }
//                     else if (x == 0.0 && y < 0.0)
//                     {
//                         phi_top = PI / 2.0 * (-1.0);
//                     }
//                     else if (x == 0.0 && y == 0.0)
//                     {
//                         phi_top = 0.0; // special points
//                     }
//                     // find related point on the earth shell in the north as (x_earth, y_earth, z_earth)
//                     double theta_earth = asin(sin(theta_top) * sqrt(r_earth / r_top));
//                     double phi_earth = phi_top;
//                     double x_earth = r_earth * sin(theta_earth) * cos(phi_earth);
//                     double y_earth = r_earth * sin(theta_earth) * sin(phi_earth);
//                     double z_earth = r_earth * cos(theta_earth);
//                     //
//                     if (x == 0.0 && y == 0.0)
//                     {
//                         x_earth = 0.0;
//                         y_earth = 0.0;
//                     }
//                     // calculate B
//                     double bmag_top = ptrArray[f][i][j][k]->B3_base().norm();
//                     double bmag_earth_x = 3.0 * dMoment * x_earth * z_earth / pow(r_earth, 5.0);
//                     double bmag_earth_y = 3.0 * dMoment * y_earth * z_earth / pow(r_earth, 5.0);
//                     double bmag_earth_z = dMoment * (3.0 * pow(z_earth, 2.0) - pow(r_earth, 2.0)) / pow(r_earth, 5.0);
//                     double bmag_earth = sqrt(bmag_earth_x * bmag_earth_x + bmag_earth_y * bmag_earth_y + bmag_earth_z * bmag_earth_z);
//                     //
//                     double s_parallel;
//                     // update s_parallel, poynting flux parallel at 100 km, gamma = 20.0
//                     double psd_w0 = 20.0 * s_parallel;
//                     //double w0 = 6.5; // Hz
//                     double psd_w = psd_w0 * pow(bmag_top / bmag_earth, -1.0 * alphaPSD);
//                     ptrArray[f][i][j][k]->SetPSD_H(psd_w);
//                     //
//                     double potential_;
//                     // update potential
//                     ptrArray[f][i][j][k]->SetPotential(potential_);
//                 }
//             }
//         }
//     }
// }

//************************************************************************
//************************************************************************
// Function
// initial the bot boundary for the  velocity of magnetic field line
// notice the k range
void SetRotationalVeBotBoundary(GridsPoints *****ptrArray_in,
                                int timeline_in)

{
    //double PI = 3.1415926535897;
    Vector3 EVector;      //manmaid E pointing out
    Vector3 externalE;    // from rotation increasing with time
    Vector3 original_vel; // rotating speed of electron
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                //    for( int k = 1; k < tempGridsCellLevel+1; k++)
                //if( update_type == 1)
                //{
                //    for( int k = 0; k < tempGridsCellLevelBot ; k++)
                //    {
                //        EVector= ptrArray_in[face][i][j][k]->Pos3().NormalizedVector().ScaleProduct(EBot_const);
                //        // Apply a external E because of rotation
                //        SetRotationalVel(ptrArray_in, face, i, j, k);
                //        original_vel = ptrArray_in[face][i][j][k]->Vel_e3();
                //    //            double x_time = (timeline_in*tstep-botBoundaryInitialTimeStart)/botBoundaryInitialTime;
                //    //            externalE = ptrArray_in[face][i][j][k]->B3_base().CrossProduct(original_vel.ScaleProduct(0.50* sin(PI*(x_time-0.50))+0.50));
                //        externalE = ptrArray_in[face][i][j][k]->B3_base().CrossProduct(original_vel);
                //    //        ptrArray_in[face][i][j][k]->SetE3( EVector.PlusProduct(externalE));
                //        ptrArray_in[face][i][j][k]->SetE3( EVector.PlusProduct(externalE));
                //    }
                //}else
                //{
                //    for( int k = 0; k < tempGridsCellLevelBot ; k++)
                //    {
                //        EVector= ptrArray_in[face][i][j][k]->Pos3().NormalizedVector().ScaleProduct(EBot_const);
                //        ptrArray_in[face][i][j][k]->SetE3( EVector);
                //    }
                //}

                for (int k = 0; k < tempGridsCellLevelBot; k++)
                {
                    EVector = ptrArray_in[face][i][j][k]->Pos3().NormalizedVector().ScaleProduct(EBot_const);
                    ptrArray_in[face][i][j][k]->SetE3(EVector);
                }
            }
        }
    }
}

void SetRotationalVeBotBoundary(GridsPoints *****ptrArray_in,
                                Vector3 *****ptrEVectorCellArray,
                                Vector3 *****ptrVeleVectorCellArray,
                                Vector3 *****ptrGradVectorCellArray,
                                int timeline_in)
{
    //double PI = 3.1415926535897;
    Vector3 externalE, original_vel;
    Vector3 EVector = Vector3(0.0, 0.0, 0.0);
    Vector3 VeleVector = Vector3(0.0, 0.0, 0.0);

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                //    for( int k = 1; k < tempGridsCellLevel+1; k++)
                for (int k = 1; k < coverGridsCellLevelBot + 1; k++)
                {
                    if (timeline_in == 0)
                    {

                        if (i == 1 && j == 1)
                        {
                            EVector = ptrEVectorCellArray[face][1][0][k - 1]->PlusProduct(
                                ptrEVectorCellArray[face][1][1][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][0][1][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][1][0][k]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][1][1][k]->V3());
                            EVector = EVector.PlusProduct(
                                                 ptrEVectorCellArray[face][0][1][k]->V3())
                                          .ScaleProduct(1.0 / 6.0);

                            VeleVector = ptrVeleVectorCellArray[face][1][0][k - 1]->PlusProduct(
                                ptrVeleVectorCellArray[face][1][1][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][0][1][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][1][0][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][1][1][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                                       ptrVeleVectorCellArray[face][0][1][k]->V3())
                                             .ScaleProduct(1.0 / 6.0);
                        }
                        else if (i == 1 && j == fieldsGridsSize + 1)
                        {

                            EVector = ptrEVectorCellArray[face][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                                ptrEVectorCellArray[face][1][fieldsGridsSize][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][0][fieldsGridsSize][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][1][fieldsGridsSize + 1][k]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][1][fieldsGridsSize][k]->V3());
                            EVector = EVector.PlusProduct(
                                                 ptrEVectorCellArray[face][0][fieldsGridsSize][k]->V3())
                                          .ScaleProduct(1.0 / 6.0);

                            VeleVector = ptrVeleVectorCellArray[face][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                                ptrVeleVectorCellArray[face][1][fieldsGridsSize][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][0][fieldsGridsSize][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][1][fieldsGridsSize + 1][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][1][fieldsGridsSize][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                                       ptrVeleVectorCellArray[face][0][fieldsGridsSize][k]->V3())
                                             .ScaleProduct(1.0 / 6.0);
                        }
                        else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                        {
                            EVector = ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k]->V3());
                            EVector = EVector.PlusProduct(
                                                 ptrEVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k]->V3())
                                          .ScaleProduct(1.0 / 6.0);

                            VeleVector = ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                                       ptrVeleVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k]->V3())
                                             .ScaleProduct(1.0 / 6.0);
                        }
                        else if (i == fieldsGridsSize + 1 && j == 1)
                        {
                            EVector = ptrEVectorCellArray[face][fieldsGridsSize][0][k - 1]->PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize][1][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize + 1][1][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize][0][k]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][fieldsGridsSize][1][k]->V3());
                            EVector = EVector.PlusProduct(
                                                 ptrEVectorCellArray[face][fieldsGridsSize + 1][1][k]->V3())
                                          .ScaleProduct(1.0 / 6.0);

                            VeleVector = ptrVeleVectorCellArray[face][fieldsGridsSize][0][k - 1]->PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize][1][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize + 1][1][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize][0][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][fieldsGridsSize][1][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                                       ptrVeleVectorCellArray[face][fieldsGridsSize + 1][1][k]->V3())
                                             .ScaleProduct(1.0 / 6.0);
                        }
                        else
                        {
                            EVector = ptrEVectorCellArray[face][i - 1][j - 1][k - 1]->PlusProduct(
                                ptrEVectorCellArray[face][i][j - 1][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][i - 1][j][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][i][j][k - 1]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][i - 1][j - 1][k]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][i][j - 1][k]->V3());
                            EVector = EVector.PlusProduct(
                                ptrEVectorCellArray[face][i - 1][j][k]->V3());
                            EVector = EVector.PlusProduct(
                                                 ptrEVectorCellArray[face][i][j][k]->V3())
                                          .ScaleProduct(1.0 / 8.0);

                            VeleVector = ptrVeleVectorCellArray[face][i - 1][j - 1][k - 1]->PlusProduct(
                                ptrVeleVectorCellArray[face][i][j - 1][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][i - 1][j][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][i][j][k - 1]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][i - 1][j - 1][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][i][j - 1][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                ptrVeleVectorCellArray[face][i - 1][j][k]->V3());
                            VeleVector = VeleVector.PlusProduct(
                                                       ptrVeleVectorCellArray[face][i][j][k]->V3())
                                             .ScaleProduct(1.0 / 8.0);
                        }

                        //       SetRotationalVel(ptrArray_in, face, i, j, k, timeline_in);

                        //          std::cout << EVector.x() << " " << EVector.y() << " " << EVector.z() << std::endl;
                        //          externalE = Vector3(0.0,0.0,0.0);
                        // update Es

                        /*            while (EVector.norm() < 1.0e-7)
                {
                    EVector = EVector.ScaleProduct( 10.0);
                }
    */
                        ptrArray_in[face][i][j][k]->SetE3(EVector);

                        //            ptrArray_in[face][i][j][k]->SetE3( externalE);

                        //                ptrArray_in[face][i][j][k]->SetVel_e3( VeleVector);

                        if (k == 1)
                        {
                            //        externalE = Vector3(0.0,0.0,0.0);
                            // update Es
                            ptrArray_in[face][i][j][0]->SetE3(EVector);
                            //                    ptrArray_in[face][i][j][0]->SetVel_e3( VeleVector);
                        }
                    }
                    else
                    {

                        // Apply a external E because of rotation
                        SetRotationalVel(ptrArray_in, face, i, j, k);
                        EVector = ptrArray_in[face][i][j][k]->E3();

                        original_vel = ptrArray_in[face][i][j][k]->Vel_e3();
                        double x_time = (timeline_in * tstep - botBoundaryInitialTimeStart) / botBoundaryInitialTime;
                        externalE = ptrArray_in[face][i][j][k]->B3_base().CrossProduct(original_vel.ScaleProduct(0.50 * sin(PI * (x_time - 0.50)) + 0.50));

                        ptrArray_in[face][i][j][0]->SetE3(EVector.PlusProduct(externalE));
                        if (k == 1)
                        {
                            EVector = ptrArray_in[face][i][j][1]->E3();
                            SetRotationalVel(ptrArray_in, face, i, j, 0);
                            original_vel = ptrArray_in[face][i][j][0]->Vel3();
                            externalE = ptrArray_in[face][i][j][0]->B3_base().CrossProduct(original_vel.ScaleProduct(0.50 * sin(PI * (x_time - 0.50)) + 0.50));

                            //            externalE = Vector3(0.0,0.0,0.0);
                            // update Es
                            ptrArray_in[face][i][j][0]->SetE3(EVector.PlusProduct(externalE));
                            //                        ptrArray_in[face][i][j][0]->SetVel_e3( VeleVector);
                        }
                    }
                }
            }
        }
    }

    // i, j, k are index of grids
}
/*
    for( int face = 0; face < totalFace; face++)
    {
        for( int i = 1; i < fieldsGridsSize+2; i++)
        {
            for( int j = 1; j < fieldsGridsSize+2; j++)
            {
                for( int k = 0; k < tempGridsCellLevel+1; k++)
                {
                // set stopSign  
                ptrArray_in[face][i][j][k]->SetStopSign(0);
                }
            }
        }
    } 
*/

//************************************************************************
//************************************************************************
// Function
// initial the top boundary for the  velocity of magnetic field line
void SetConvectionVelTopBoundary(GridsPoints *****ptrArray_in, int timeline_in)
{
    //double PI = 3.1415926535897;
    // input two const
    //double r0 = radius * cos(r0_latitude * PI / 180.0);
    //double c0 = radius * cos(c0_latitude * PI / 180.0);
    //double t0 = t0_convection;
    //double r_earth = radius;

    std::cout << " Set Top boundary " << std::endl;
    //non-circle
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                int k = fieldsGridsSize;

                if (ptrArray_in[face][i][j][k]->StopSign() == 1)
                    continue;
                SetConvectionVel(ptrArray_in, face, i, j, k);
                Vector3 original_vel = ptrArray_in[face][i][j][k]->Vel3();
                ptrArray_in[face][i][j][k]->SetVel_Boundary(original_vel.ScaleProduct(sin(PI / 2.0 * (timeline_in * tstep - topBoundaryInitialTimeStart) / topBoundaryInitialTime)));

                ptrArray_in[face][i][j][k]->SetGradPe(ptrArray_in[face][i][j][k - 1]->GradPe());
                ptrArray_in[face][i][j][k]->Density_H(ptrArray_in[face][i][j][k - 1]->Density_H());
                ptrArray_in[face][i][j][k]->Density_He(ptrArray_in[face][i][j][k - 1]->Density_He());
                ptrArray_in[face][i][j][k]->Density_O(ptrArray_in[face][i][j][k - 1]->Density_O());
                //double r = ptrArray_in[face][i][j][k]->Pos3().norm() / radius;
                //if( r > 0){
                //ptrArray_in[face][i][j][k]->Density_H( N0_H / r * ( 1.0 - tanh( r - 6.5)));
                //ptrArray_in[face][i][j][k]->Density_He( N0_He / r * ( 1.0 - tanh( r - 6.5)));
                //ptrArray_in[face][i][j][k]->Density_O( N0_O / r * ( 1.0 - tanh( r - 6.5)));
                //}
                // set E
                ptrArray_in[face][i][j][k]->updateE(ptrArray_in[face][i][j][k]->GradPe());

                // set stopSign
                ptrArray_in[face][i][j][k]->SetStopSign(1);
            }
        }
    }

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                int k = fieldsGridsSize;

                // set stopSign
                ptrArray_in[face][i][j][k]->SetStopSign(0);
            }
        }
    }
}

//************************************************************************
//************************************************************************
// Function
// Set velocity due to earth rotation
void SetRotationalVel(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{
    Vector3 tempPos = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3();
    Vector3 tempOmega = Vector3(0.0, 0.0, omega_earth);
    tempPos.Setz(0.0);
    ptrArray_in[face_in][i_in][j_in][k_in]->SetVe_Boundary(tempOmega.CrossProduct(tempPos));
}

//************************************************************************
//************************************************************************
// Function
// Set velocity due to earth rotation
void SetRotationalVel(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in, int timeline)
{
    Vector3 tempPos = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3();
    Vector3 tempOmega = Vector3(0.0, 0.0, omega_earth);
    tempPos.Setz(0.0);
    double x_time = (timeline * tstep - botBoundaryInitialTimeStart) / botBoundaryInitialTime;
    ptrArray_in[face_in][i_in][j_in][k_in]->SetVe_Boundary(tempOmega.CrossProduct(tempPos).ScaleProduct(0.50 * sin(3.1415926535898 * (x_time - 0.50)) + 0.50));
}
//************************************************************************
//************************************************************************
// Function
// Set velocity due to two convection cell patern
void SetConvectionVel(GridsPoints *****ptrArray_in, int face_in, int i_in, int j_in, int k_in)
{

    //double PI = 3.1415926535897;
    // input two const
    double r0 = radius * cos(r0_latitude * PI / 180.0);
    double c0 = radius * cos(c0_latitude * PI / 180.0);
    double t0 = t0_convection;
    double r_earth = radius;

    double x = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().x();
    double y = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().y();
    double z = ptrArray_in[face_in][i_in][j_in][k_in]->Pos3().z();

    // Set velocity
    double r_top = sqrt(x * x + y * y + z * z);
    double theta_top = acos(z / r_top);
    double phi_top;

    phi_top=atan2(y,x);
    /*if (x != 0.0)
    {
        phi_top = atan(y / x);
        if (x < 0)
        {
            phi_top = phi_top + PI;
        }
    }
    else if (x == 0.0 && y > 0.0)
    {
        phi_top = PI / 2.0;
    }
    else if (x == 0.0 && y < 0.0)
    {
        phi_top = PI / 2.0 * (-1.0);
    }
    else if (x == 0.0 && y == 0.0)
    {
        phi_top = 0.0; // special points
    }*/

    // Step 1
    // find related point on the earth shell in the north
    double theta_earth = asin(sin(theta_top) * sqrt(r_earth / r_top));
    double phi_earth = phi_top;
    double x_earth = r_earth * sin(theta_earth) * cos(phi_earth);
    double y_earth = r_earth * sin(theta_earth) * sin(phi_earth);
    //double z_earth = r_earth * cos(theta_earth);

    if (x == 0.0 && y == 0.0)
    {
        x_earth = 0.0;
        y_earth = 0.0;
    }
    //????????????????????????
    // Step 2
    // find the velocity of the point ( |x_earth|, y_earth) on the x-y plane
    double xx = x_earth;
    double yy = y_earth;
    double vx_earth = 0.0;
    double vy_earth = 0.0;
    double L; // distance, half of period for sin function

    // cout << xx << " " << yy << " " << r0 << endl;

    double x_prime, y_prime; // used for region 2
    if (y_earth < 0.0)
    {
        yy = -1.0 * yy;
    }
    if (xx <= r0 - r0 * yy / c0 && xx >= -1.0 * r0 + yy * r0 / c0) // region 1
    {
        //    cout << " test 1 " << endl;
        L = 2.0 * r0 * (1.0 - yy / c0);
        vx_earth = -1.0 * PI * L / t0 * sqrt(0.25 - xx * xx / L / L);
        vy_earth = 0.0;
    }
    else if (xx * xx + yy * yy <= r0 * r0) // region 2
    {
        //    cout << " test 2 " << endl;
        y_prime = (yy - r0 * r0 / c0 + sqrt((yy - r0 * r0 / c0) * (yy - r0 * r0 / c0) - (r0 * r0 / c0 / c0 - 1.0) * (r0 * r0 - xx * xx - yy * yy))) / (1.0 - r0 * r0 / c0 / c0);
        x_prime = r0 * (1.0 - y_prime / c0);

        L = 2.0 * x_prime; // x direction
        if (L != 0.0 && xx < x_prime - 1e-6)
        {
            vx_earth = PI * L / t0 * sqrt(0.25 - xx * xx / L / L);
        }
        else
        {
            vx_earth = 0.0;
        }

        L = x_prime; // y direction
        if (L != 0.0 && (yy - y_prime) / L < 1.0 - 1e-6)
        {
            vy_earth = PI * L / t0 * sqrt((yy - y_prime) / L * (1.0 - (yy - y_prime) / L));
        }
        else
        {
            vy_earth = 0.0;
        }

        if (xx > 0.0)
        {
            vy_earth = -1.0 * vy_earth;
        }
    }
    else // other places
    {
        //        cout << " test 3 " << endl;
        vx_earth = 0.0;
        vy_earth = 0.0;
    }
    if (y_earth < 0.0)
    {
        vy_earth *= -1.0;
    }

    /*
    // Rotate for a x-axis pointed to the sun
    double vx_earth_rotate = vy_earth;
    double vy_earth_rotate = vx_earth;
    vx_earth = vx_earth_rotate;
    vy_earth = vy_earth_rotate; 
*/
    // cout << vx_earth << " " << vy_earth << " >> " << endl;
    /*    if( k_in == 16 && i_in == fieldsGridsSize/2 +1 &&  j_in == 5 && face_in == 5) {
    cout << face_in << " " << i_in << " " << j_in << " " << k_in << " vel_earth " << vx_earth << " " << vy_earth <<
     " pos " << x << " " << y << " " << z << endl;}
 */
    // Step 3
    // find the realted velocity on the earth ( x_earth, y_earth, z_earth) or ( r_earth, theta_earth, phi_earth)
    // the velocity on the x and y direction is known as ( vx_earth, vy_earth)
    // Notice the polar points have phi_earth = 0 which is incorrect to use general form or vtheta_earth will be zero

    double vtheta_earth = (vx_earth * cos(phi_earth) + vy_earth * sin(phi_earth)) / cos(theta_earth);
    double vphi_earth = vy_earth * cos(phi_earth) - vx_earth * sin(phi_earth); // notice polar points
    if (theta_earth <= 1e-5)
    {
        vtheta_earth = vx_earth;
    }

    // cout << " vtheta_earth " << vtheta_earth << " vphi_earth " << vphi_earth ; // vx_earth & vy_earth are zero

    // Step 4
    // find the related velocity on the arbitrary shell as we want using the equation A21 and A22 of
    // Rasmussen et.al 1992
    double sinchi_top = 2.0 * cos(theta_top) / sqrt(1.0 + 3.0 * cos(theta_top) * cos(theta_top));
    double coschi_top = sin(theta_top) / sqrt(1.0 + 3.0 * cos(theta_top) * cos(theta_top));

    double vr_top, vtheta_top, vphi_top;
    if (theta_earth > 1e-5)
    {
        vr_top = r_top * coschi_top * coschi_top * vtheta_earth / r_earth * 2.0 * cos(theta_earth) / sin(theta_earth);
        vtheta_top = r_top * sinchi_top * coschi_top * vtheta_earth / r_earth * 2.0 * cos(theta_earth) / sin(theta_earth);
        vphi_top = vphi_earth / r_earth * r_top;
    }
    else
    {
        vr_top = 0.0;
        vtheta_top = r_top * vtheta_earth / r_earth;
        vphi_top = 0.0;
    }

    /*    
    if( k_in == 16 && i_in == fieldsGridsSize/2 +1 &&  j_in == 5 && face_in == 5) {
    cout << face_in << " " << i_in << " " << j_in << " " << k_in << " vr_top " << vr_top << " vtheta_top " << vtheta_top << " vphi_top " << vphi_top << " " <<
     " pos " << x << " " << y << " " << z << endl;}
*/
    //  cout << " vr_top " << vr_top << " vtheta_top " << vtheta_top << " vphi_top " << vphi_top << " >>>> " << endl; //

    double vx_top = vr_top * sin(theta_top) * cos(phi_top) +
                   vtheta_top * cos(theta_top) * cos(phi_top) -
                   vphi_top * sin(phi_top);
    double vy_top = vr_top * sin(theta_top) * sin(phi_top) +
                   vtheta_top * cos(theta_top) * sin(phi_top) +
                   vphi_top * cos(phi_top);
    double vz_top = vr_top * cos(theta_top) -
                   vtheta_top * sin(theta_top);

    if (theta_earth <= 1e-5)
    {
        vx_top = vtheta_top;
        vy_top = 0.0;
        vz_top = 0.0;
    }
    /*
//  cout << " r_top " << r_top << " theta_top " << theta_top << " >>> "; // vtheta_earth is zero
if( k_in == 16 && j_in == fieldsGridsSize / 2 && face_in == 5 )
{
    cout <<endl<< face_in << " " << i_in << " " << j_in << " " << k_in << endl;
    cout << vx_earth << " " << vy_earth << " theta " << theta_earth << " " << theta_top << " phi " << phi_earth << " " << phi_top
    << " phi_top " << phi_top << " vtheta " << vtheta_earth << " vr " << vr_top* sin(theta_top) << " vphi_top " << vphi_top
    << " v_top " << vx_top << " " << vy_top << " " << vz_top << endl;//
}
*/

    // cout << vx_top << " " << vy_top << " " << vz_top << endl;
    Vector3 temp = Vector3(vx_top, vy_top, vz_top);
    Vector3 original_vel = ptrArray_in[face_in][i_in][j_in][k_in]->Vel3();
    ptrArray_in[face_in][i_in][j_in][k_in]->SetVel_Boundary(original_vel.PlusProduct(temp));

    /*
    if( k_in == 16 && i_in == fieldsGridsSize/2 +1 &&  j_in == 5 && face_in == 5) {
    cout << face_in << " " << i_in << " " << j_in << " " << k_in << " vel " << vx_top << " " << vy_top << " " << vz_top << " " <<
     " pos " << x << " " << y << " " << z << endl;}
*/
}

//************************************************************************
//************************************************************************
// FUNCTION
// Calculate the gradient of normal of B
void GradBNorm(GridsPoints *****ptrArray_in)
{

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 0; k < fieldsGridsSize + 1; k++)
                {

                    //check stopsign
                    if (ptrArray_in[face][i][j][k]->StopSign() == 1)
                        continue;

                    ptrArray_in[face][i][j][k]->XYZtoGradBNorm();

                    // set stopSigns
                    ptrArray_in[face][i][j][k]->SetStopSign(1);
                }
            }
        }
    }

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

//************************************************************************
//************************************************************************
// FUNCTION
// Set zero for pho and v at each points
void ResetPhoVatGrids(GridsPoints *****ptrArray_in)
{
    // reset, clear the previous value: density and velocity
#pragma omp parallel for collapse(4)
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    ptrArray_in[face][i][j][k]->ResetParameters();
                }
            }
        }
    }
}

//************************************************************************
//************************************************************************
// FUNCTION
// finish culmulating pho and velocity and then average the density and velocity
// 1. add info from weightArray to gridsArray
// 2. divide by the volume to calculate the
void CalculatingAveragedPhoVatGrids(GridsPoints *****ptrArray_in,
                                    double ***ptrVolumeGridArray_in,
                                    int updateInfoPeriod_in,
                                    Vector3 ******ptrVelWeightGridsArray,
                                    double ******ptrMassWeightGridsArray)
{

    int dim = fieldsGridsSize / 2;
    for (int face = 0; face < totalFace; face++)
    {
        for (int pos = 0; pos < 8; pos++)
        {
#pragma omp parallel for collapse(3)
            for (int i = 0; i < dim + 1; i++)
            {
                for (int j = 0; j < dim + 1; j++)
                {
                    for (int k = 0; k < dim + 1; k++)
                    {
                        switch (pos)
                        {
                        case 0:
                            ptrArray_in[face][i + 1][j + 1][k]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 1:
                            ptrArray_in[face][i + 1 + dim][j + 1][k]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 2:
                            ptrArray_in[face][i + 1][j + 1 + dim][k]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 3:
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 4:
                            ptrArray_in[face][i + 1][j + 1][k + dim]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k + dim]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k + dim]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k + dim]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k + dim]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1][k + dim]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 5:
                            ptrArray_in[face][i + 1 + dim][j + 1][k + dim]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k + dim]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k + dim]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k + dim]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k + dim]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1][k + dim]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 6:
                            ptrArray_in[face][i + 1][j + 1 + dim][k + dim]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k + dim]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k + dim]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k + dim]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k + dim]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1][j + 1 + dim][k + dim]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;
                        case 7:
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k + dim]->PlusMassH3(ptrMassWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k + dim]->PlusMassHe3(ptrMassWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k + dim]->PlusMassO3(ptrMassWeightGridsArray[2][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k + dim]->PlusVelH3(ptrVelWeightGridsArray[0][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k + dim]->PlusVelHe3(ptrVelWeightGridsArray[1][face][pos][i][j][k]);
                            ptrArray_in[face][i + 1 + dim][j + 1 + dim][k + dim]->PlusVelO3(ptrVelWeightGridsArray[2][face][pos][i][j][k]);
                            break;

                        default:
                            std::cout << " calculate average pho and vel failed " << std::endl;
                            break;
                        }
                    }
                }
            }
            // omp end
        }
    }

// reset to zero
#pragma omp parallel for collapse(6)
    for (int ionType = 0; ionType < 3; ionType++)
    {
        for (int face = 0; face < totalFace; face++)
        {
            for (int pos = 0; pos < 8; pos++)
            {
                for (int i = 0; i < fieldsGridsSize / 2 + 1; i++)
                {
                    for (int j = 0; j < fieldsGridsSize / 2 + 1; j++)
                    {
                        for (int k = 0; k < fieldsGridsSize / 2 + 1; k++)
                        {
                            ptrMassWeightGridsArray[ionType][face][pos][i][j][k] = 0.0;
                            ptrVelWeightGridsArray[ionType][face][pos][i][j][k] = Vector3{0.0, 0.0, 0.0};
                        }
                    }
                }
            }
        }
    }

#pragma omp parallel for collapse(4)
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1 + tempGridsCellLevelBot - coverGridsCellLevelBot; k < fieldsGridsSize - tempGridsCellLevelTop + coverGridsCellLevelTop; k++)
                {
                    //check stopsign
                    if (ptrArray_in[face][i][j][k]->StopSign() == 1)
                        continue;

                    //double tempDensity;
                    // set volume
                    double volume = ptrVolumeGridArray_in[i - 1][j - 1][k]; // face of ptrArray is greater than that of ptrVolumeGridArray
                                                                           //                   std::cout << face << i << j << k << " " ;
                                                                           //                   std::cout << " volume " << volume << " density " << ptrArray_in[face][i][j][k]->Density() << " ==> " ;

                    ptrArray_in[face][i][j][k]->UpdateDueToWgt(ptrArray_in, volume, updateInfoPeriod_in);

                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(1);
                    //              std::cout << ptrArray_in[face][i][j][k]->Density() << std::endl;
                }
            }
        }
    }
// reset stopsign
#pragma omp parallel for collapse(4)
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

//************************************************************************
//************************************************************************
// FUNCTION
// finish culmulating pho and velocity and then average the density and velocity
// 1. add info from weightArray to gridsArray
// 2. divide by the volume to calculate the
void CalculatingAveragedPhoVatGrids(GridsPoints *****ptrArray_in,
                                    double ***ptrVolumeGridArray_in,
                                    int updateInfoPeriod_in)
{
    // i, j, k, index of grids
#pragma omp parallel for collapse(4)
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    //check stopsign
                    if (ptrArray_in[face][i][j][k]->StopSign() == 1)
                        continue;
                    //double tempDensity;
                    // set volume
                    double volume = ptrVolumeGridArray_in[i - 1][j - 1][k]; // face of ptrArray is greater than that of ptrVolumeGridArray
                                                                           //                   std::cout << face << i << j << k << " " ;
                                                                           //                   std::cout << " volume " << volume << " density " << ptrArray_in[face][i][j][k]->Density() << " ==> " ;
#pragma omp critical
                    ptrArray_in[face][i][j][k]->UpdateDueToWgt(ptrArray_in, volume, updateInfoPeriod_in);
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(1);
                    //              std::cout << ptrArray_in[face][i][j][k]->Density() << std::endl;
                }
            }
        }
    }
// reset stopsign
#pragma omp parallel for collapse(4)
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    // set stopSign
                    ptrArray_in[face][i][j][k]->SetStopSign(0);
                }
            }
        }
    }
}

//************************************************************************
//************************************************************************
// Create Cell centered field array for E for the type of Vector3
// The size of this array is [totalface * fsize+2 * fsize+2 * fsize*2]
Vector3 *****PtrVectorCellArray(Vector3 *mem_VectorCellArray)
{
    // define the pointer
    Vector3 *****ptrEArray = new Vector3 ****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrEArray[face] = new Vector3 ***[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrEArray[face][i] = new Vector3 **[fieldsGridsSize + 2];
            for (int j = 0; j < fieldsGridsSize + 2; j++)
            {
                ptrEArray[face][i][j] = new Vector3 *[fieldsGridsSize * grid_domain];
                for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                {
                    if (1 <= i && i < fieldsGridsSize + 1 && 1 <= j && j < fieldsGridsSize + 1)
                    {
                        ptrEArray[face][i][j][k] = mem_VectorCellArray + face * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain + (i - 1) * fieldsGridsSize * fieldsGridsSize * grid_domain + (j - 1) * fieldsGridsSize * grid_domain + k;
                    }
                    else
                    {
                        continue;
                    }
                }
            }
        }
    }
    // covered area
    // face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEArray[0][i][0][k] = ptrEArray[5][i][fieldsGridsSize][k]; // bot

            ptrEArray[0][fieldsGridsSize + 1][i][k] = ptrEArray[1][1][i][k]; // right

            ptrEArray[0][i][fieldsGridsSize + 1][k] = ptrEArray[2][i][1][k]; // top

            ptrEArray[0][0][i][k] = ptrEArray[4][fieldsGridsSize][i][k]; // left
        }
    }
    // face 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEArray[1][i][0][k] = ptrEArray[5][fieldsGridsSize][fieldsGridsSize + 1 - i][k]; // bot

            ptrEArray[1][fieldsGridsSize + 1][i][k] = ptrEArray[3][1][i][k]; // right

            ptrEArray[1][i][fieldsGridsSize + 1][k] = ptrEArray[2][fieldsGridsSize][i][k]; // top

            ptrEArray[1][0][i][k] = ptrEArray[0][fieldsGridsSize][i][k]; // left
        }
    }
    // face 2
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEArray[2][i][0][k] = ptrEArray[0][i][fieldsGridsSize][k]; // bot

            ptrEArray[2][fieldsGridsSize + 1][i][k] = ptrEArray[1][i][fieldsGridsSize][k]; // right

            ptrEArray[2][i][fieldsGridsSize + 1][k] = ptrEArray[3][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // top

            ptrEArray[2][0][i][k] = ptrEArray[4][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // left
        }
    }
    //face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEArray[3][i][0][k] = ptrEArray[5][fieldsGridsSize + 1 - i][1][k]; // bot

            ptrEArray[3][fieldsGridsSize + 1][i][k] = ptrEArray[4][1][i][k]; // right

            ptrEArray[3][i][fieldsGridsSize + 1][k] = ptrEArray[2][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // top

            ptrEArray[3][0][i][k] = ptrEArray[1][fieldsGridsSize][i][k]; // left
        }
    }
    // face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEArray[4][i][0][k] = ptrEArray[5][1][i][k]; // bot

            ptrEArray[4][fieldsGridsSize + 1][i][k] = ptrEArray[0][1][i][k]; // right

            ptrEArray[4][i][fieldsGridsSize + 1][k] = ptrEArray[2][1][fieldsGridsSize + 1 - i][k]; // top

            ptrEArray[4][0][i][k] = ptrEArray[3][fieldsGridsSize][i][k]; // left
        }
    }
    // face 5
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEArray[5][i][0][k] = ptrEArray[3][fieldsGridsSize + 1 - i][1][k]; // bot

            ptrEArray[5][fieldsGridsSize + 1][i][k] = ptrEArray[1][fieldsGridsSize + 1 - i][1][k]; // right

            ptrEArray[5][i][fieldsGridsSize + 1][k] = ptrEArray[0][i][1][k]; // top

            ptrEArray[5][0][i][k] = ptrEArray[4][i][1][k]; // left
        }
    }

    return ptrEArray;
}

//************************************************************************
//************************************************************************
// Prerun 1.6 // Create const array of B at center of each cells
// [totalface * fsize+2 * fsize+2 * fsize +2]
// Most code are the same with that of EVectorCellArray
Vector3 *****BVectorCellArray(GridsPoints *****ptrArray)
{
    // Apply space to store
    static Vector3 *mem_BVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize];
    // define the pointer
    Vector3 *****ptrBArray = new Vector3 ****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrBArray[face] = new Vector3 ***[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrBArray[face][i] = new Vector3 **[fieldsGridsSize + 2];
            for (int j = 0; j < fieldsGridsSize + 2; j++)
            {
                ptrBArray[face][i][j] = new Vector3 *[fieldsGridsSize];
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    if (1 <= i && i < fieldsGridsSize + 1 && 1 <= j && j < fieldsGridsSize + 1)
                    {
                        ptrBArray[face][i][j][k] = mem_BVectorCellArray + face * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize + (i - 1) * fieldsGridsSize * fieldsGridsSize + (j - 1) * fieldsGridsSize + k;
                    }
                    else
                    {
                        continue;
                    }
                }
            }
        }
    }
    // covered area
    // face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBArray[0][i][0][k] = ptrBArray[5][i][fieldsGridsSize][k]; // bot

            ptrBArray[0][fieldsGridsSize + 1][i][k] = ptrBArray[1][1][i][k]; // right

            ptrBArray[0][i][fieldsGridsSize + 1][k] = ptrBArray[2][i][1][k]; // top

            ptrBArray[0][0][i][k] = ptrBArray[4][fieldsGridsSize][i][k]; // left
        }
    }
    // face 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBArray[1][i][0][k] = ptrBArray[5][fieldsGridsSize][fieldsGridsSize + 1 - i][k]; // bot

            ptrBArray[1][fieldsGridsSize + 1][i][k] = ptrBArray[3][1][i][k]; // right

            ptrBArray[1][i][fieldsGridsSize + 1][k] = ptrBArray[2][fieldsGridsSize][i][k]; // top

            ptrBArray[1][0][i][k] = ptrBArray[0][fieldsGridsSize][i][k]; // left
        }
    }
    // face 2
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBArray[2][i][0][k] = ptrBArray[0][i][fieldsGridsSize][k]; // bot

            ptrBArray[2][fieldsGridsSize + 1][i][k] = ptrBArray[1][i][fieldsGridsSize][k]; // right

            ptrBArray[2][i][fieldsGridsSize + 1][k] = ptrBArray[3][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // top

            ptrBArray[2][0][i][k] = ptrBArray[4][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // left
        }
    }
    //face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBArray[3][i][0][k] = ptrBArray[5][fieldsGridsSize + 1 - i][1][k]; // bot

            ptrBArray[3][fieldsGridsSize + 1][i][k] = ptrBArray[4][1][i][k]; // right

            ptrBArray[3][i][fieldsGridsSize + 1][k] = ptrBArray[2][fieldsGridsSize + 1 - i][fieldsGridsSize][k]; // top

            ptrBArray[3][0][i][k] = ptrBArray[1][fieldsGridsSize][i][k]; // left
        }
    }
    // face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBArray[4][i][0][k] = ptrBArray[5][1][i][k]; // bot

            ptrBArray[4][fieldsGridsSize + 1][i][k] = ptrBArray[0][1][i][k]; // right

            ptrBArray[4][i][fieldsGridsSize + 1][k] = ptrBArray[2][1][fieldsGridsSize + 1 - i][k]; // top

            ptrBArray[4][0][i][k] = ptrBArray[3][fieldsGridsSize][i][k]; // left
        }
    }
    // face 5
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBArray[5][i][0][k] = ptrBArray[3][fieldsGridsSize + 1 - i][1][k]; // bot

            ptrBArray[5][fieldsGridsSize + 1][i][k] = ptrBArray[1][fieldsGridsSize + 1 - i][1][k]; // right

            ptrBArray[5][i][fieldsGridsSize + 1][k] = ptrBArray[0][i][1][k]; // top

            ptrBArray[5][0][i][k] = ptrBArray[4][i][1][k]; // left
        }
    }

    Vector3 tempB;
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    //Set the Vector at the center of cells
                }
            }
        }
    }

    return ptrBArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the B on face

// The size of this array is [face * (fsize+2) * (fsize+2) * (fsize) * direction
Vector3 ******BVectorFaceArray(Vector3 ******ptrBVectorFaceArray)
{
    ptrBVectorFaceArray = new Vector3 *****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrBVectorFaceArray[face] = new Vector3 ****[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrBVectorFaceArray[face][i] = new Vector3 ***[fieldsGridsSize + 2];
            {
                for (int j = 0; j < fieldsGridsSize + 2; j++)
                {
                    ptrBVectorFaceArray[face][i][j] = new Vector3 **[fieldsGridsSize * grid_domain];
                    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                    {
                        ptrBVectorFaceArray[face][i][j][k] = new Vector3 *[3];
                    }
                }
            }
        }
    }

    // ptr[face][i][j][k][dir]
    // i,j is cell index from 0-fsize+1; k is from 0-(fsize-1)
    // dir is the perticular face index, 0 is perpendicular to the first direction

    // vector initialized face by face
    // face 0 ( to us)
    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[0][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][1] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // right
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][0] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // face 1( on the right)
    // share with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[1][0][j][k][dir] = ptrBVectorFaceArray[0][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][1][j][k][0] = ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][0];
        }
    }
    // j = 1 rest faces
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][1][j][k][1] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[1][1][j][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[1][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][i][fieldsGridsSize + 1][k][1] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // right
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][0] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    //connect face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[1][1][j][k][1];
            ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[1][1][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][1][fieldsGridsSize + 1][k][1];
    }

    // face 3 ( on the back)
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[3][0][j][k][dir] = ptrBVectorFaceArray[1][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][1][j][k][0] = ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][0];
        }
    }
    // j = 1 rest faces
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][1][j][k][1] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[3][1][j][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[3][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][i][fieldsGridsSize + 1][k][1] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // right
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][0] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    //connect face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[3][1][j][k][1];
            ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[3][1][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[3][1][fieldsGridsSize + 1][k][1];
    }
    // face 4 ( on the left)
    // share with face 3
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[4][0][j][k][dir] = ptrBVectorFaceArray[3][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][1][j][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][0];
        }
    }
    // i = 1 rest faces
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][1][j][k][1] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[4][1][j][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[4][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][i][fieldsGridsSize + 1][k][1] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // right share with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[4][fieldsGridsSize + 1][j][k][dir] = ptrBVectorFaceArray[0][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[0][1][fieldsGridsSize + 1][k][1];
    }
    // connect with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[0][0][j][k][dir] = ptrBVectorFaceArray[4][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    // connect with face 3
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[4][1][j][k][1];
            ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[4][1][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][1][fieldsGridsSize + 1][k][1];
    }
    // face 2 ( on the top)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[2][i][0][k][dir] = ptrBVectorFaceArray[0][i][fieldsGridsSize][k][dir];
            }

            ptrBVectorFaceArray[2][i][1][k][1] = ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize][k][0];
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][fieldsGridsSize + 1][j][k][0] = ptrBVectorFaceArray[1][j][fieldsGridsSize + 1][k][1];
            ptrBVectorFaceArray[2][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[1][j][fieldsGridsSize][k][0];
            ptrBVectorFaceArray[2][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[1][j][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize][k][0];
    }
    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 2][fieldsGridsSize][k][0];
            ptrBVectorFaceArray[2][i][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][fieldsGridsSize + 1][k][1];
            ptrBVectorFaceArray[2][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[3][1][fieldsGridsSize][k][0];
    }
    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][0][j][k][0] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 1][fieldsGridsSize][k][1];
            ptrBVectorFaceArray[2][0][j][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 2][fieldsGridsSize][k][0];
            ptrBVectorFaceArray[2][0][j][k][2] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 1][fieldsGridsSize][k][2];

            ptrBVectorFaceArray[2][1][j][k][0] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 1][fieldsGridsSize + 1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][1][fieldsGridsSize][k][0];
    }
    // rest of j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][i][1][k][0] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[2][i][1][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // rest of i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][1][j][k][1] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[2][1][j][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // rest of i = 1 && j = 1
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][1][1][k][2] = new Vector3(0.0, 0.0, 0.0);
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[2][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // connect with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][i][1][k][0];
            ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][i][1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize + 1][1][k][0];
    }
    // connect with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][j][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize][j][k][1];
            ptrBVectorFaceArray[1][j][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][fieldsGridsSize][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }

    // connect with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize - i + 2][fieldsGridsSize][k][0];

            ptrBVectorFaceArray[3][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][fieldsGridsSize - i + 1][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][1][fieldsGridsSize][k][0];
    }

    // connect with face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][1][fieldsGridsSize - i + 2][k][1];
            ptrBVectorFaceArray[4][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][1][fieldsGridsSize - i + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][1][1][k][1];
    }
    // face 5 ( bot)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[5][i][fieldsGridsSize + 1][k][dir] = ptrBVectorFaceArray[0][i][1][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[0][fieldsGridsSize + 1][1][k][0];
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {

            ptrBVectorFaceArray[5][fieldsGridsSize + 1][j][k][0] = ptrBVectorFaceArray[1][1][fieldsGridsSize - j + 1][k][1];
            ptrBVectorFaceArray[5][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[1][1][fieldsGridsSize - j + 2][k][0];
            ptrBVectorFaceArray[5][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[1][1][fieldsGridsSize - j + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][1][1][k][0];
    }

    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[5][i][0][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 2][1][k][0];
            ptrBVectorFaceArray[5][i][0][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][2][k][1];
            ptrBVectorFaceArray[5][i][0][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][1][k][2];

            ptrBVectorFaceArray[5][i][1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[3][1][1][k][0];
    }

    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {

            ptrBVectorFaceArray[5][0][j][k][0] = ptrBVectorFaceArray[4][j][2][k][1];
            ptrBVectorFaceArray[5][0][j][k][1] = ptrBVectorFaceArray[4][j][1][k][0];
            ptrBVectorFaceArray[5][0][j][k][2] = ptrBVectorFaceArray[4][j][1][k][2];

            ptrBVectorFaceArray[5][1][j][k][0] = ptrBVectorFaceArray[4][j][1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize + 1][1][k][0];
    }
    // rest of j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[5][i][1][k][0] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[5][i][1][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // rest of i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[5][1][j][k][1] = new Vector3(0.0, 0.0, 0.0);
            ptrBVectorFaceArray[5][1][j][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // rest of i = 1 && j = 1
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][1][1][k][2] = new Vector3(0.0, 0.0, 0.0);
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[5][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // connect with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[0][i][0][k][dir] = ptrBVectorFaceArray[5][i][fieldsGridsSize][k][dir];
            }
        }
    }

    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize][k][0];
    }
    // connect with face 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][i][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 2][k][1];
            ptrBVectorFaceArray[1][i][0][k][1] = ptrBVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 1][k][0];
            ptrBVectorFaceArray[1][i][0][k][2] = ptrBVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize][1][k][1];
    }

    // connect with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][i][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize - i + 2][1][k][0];
            ptrBVectorFaceArray[3][i][0][k][1] = ptrBVectorFaceArray[5][fieldsGridsSize - i + 1][2][k][1];
            ptrBVectorFaceArray[3][i][0][k][2] = ptrBVectorFaceArray[5][fieldsGridsSize - i + 1][1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][1][1][k][0];
    }

    // connect with face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {

            ptrBVectorFaceArray[4][i][0][k][0] = ptrBVectorFaceArray[5][1][i][k][1];
            ptrBVectorFaceArray[4][i][0][k][1] = ptrBVectorFaceArray[5][2][i][k][0];
            ptrBVectorFaceArray[4][i][0][k][2] = ptrBVectorFaceArray[5][1][i][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][1][fieldsGridsSize + 1][k][1];
    }

    return ptrBVectorFaceArray;
}
//////////////////////////////////
/*

    for( int j = 1; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir =0; dir<3; dir++)
            {
                ptrBVectorFaceArray[3][0][j][k][dir] = ptrBVectorFaceArray[1][fieldsGridsSize][j][k][dir];
                ptrBVectorFaceArray[3][1][j][k][dir] = ptrBVectorFaceArray[1][fieldsGridsSize+1][j][k][dir];
            }
        }
    }
    for( int i = 2; i < fieldsGridsSize+2; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir<3 ; dir++)
                {
                   ptrBVectorFaceArray[3][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);
                }
            }
        }
    }
    // face 4 ( on the left)
    // share with face 3
    for( int j = 1; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir =0; dir<3; dir++)
            {
                ptrBVectorFaceArray[4][0][j][k][dir] = ptrBVectorFaceArray[3][fieldsGridsSize][j][k][dir];
                ptrBVectorFaceArray[4][1][j][k][dir] = ptrBVectorFaceArray[3][fieldsGridsSize+1][j][k][dir];
            }
        }
    }
    // share with face 0
    for( int j = 1; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir =0; dir<3; dir++)
            {
                ptrBVectorFaceArray[4][fieldsGridsSize+1][j][k][dir] = ptrBVectorFaceArray[0][1][j][k][dir];
            }
        }
    }
    for( int i = 2; i < fieldsGridsSize+1; i++)
    {
        for( int j = 1; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir<3 ; dir++)
                {
                   ptrBVectorFaceArray[4][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);
                }
            }
        }
    }
    // connect with face 0
    for( int j = 1; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir =0; dir<3; dir++)
            {
                ptrBVectorFaceArray[0][0][j][k][dir] = ptrBVectorFaceArray[4][fieldsGridsSize][j][k][dir];
            }
        }
    }
    // face 2( on the top)
    // share with face 0
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[2][i][0][k][dir] = ptrBVectorFaceArray[0][i][fieldsGridsSize][k][dir];
                ptrBVectorFaceArray[2][i][1][k][dir] = ptrBVectorFaceArray[0][i][fieldsGridsSize+1][k][dir];
            }
        }
    }
    // share with face 1
    for( int j = 2; j< fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
                ptrBVectorFaceArray[2][fieldsGridsSize+1][j][k][0] = ptrBVectorFaceArray[1][j][fieldsGridsSize+1][k][1];
                ptrBVectorFaceArray[2][fieldsGridsSize+1][j][k][1] = ptrBVectorFaceArray[1][j][fieldsGridsSize][k][0];
                ptrBVectorFaceArray[2][fieldsGridsSize+1][j][k][2] = ptrBVectorFaceArray[1][j][fieldsGridsSize][k][2];
                
                ptrBVectorFaceArray[2][fieldsGridsSize][j][k][0] = new Vector3( 0.0, 0.0, 0.0);
                ptrBVectorFaceArray[2][fieldsGridsSize][j][k][1] = ptrBVectorFaceArray[1][j][fieldsGridsSize+1][k][0];
                ptrBVectorFaceArray[2][fieldsGridsSize][j][k][2] = ptrBVectorFaceArray[1][j][fieldsGridsSize+1][k][2];
        }
    }
    // share with face 3 
    for( int i = 1; i < fieldsGridsSize; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
                ptrBVectorFaceArray[2][i][fieldsGridsSize+1][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize-i+2][fieldsGridsSize][k][0];
                ptrBVectorFaceArray[2][i][fieldsGridsSize+1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize-i+2][fieldsGridsSize][k][1];
                ptrBVectorFaceArray[2][i][fieldsGridsSize+1][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize-i+1][fieldsGridsSize][k][2];

                ptrBVectorFaceArray[2][i][fieldsGridsSize][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize-i+2][fieldsGridsSize+1][k][0];
                ptrBVectorFaceArray[2][i][fieldsGridsSize][k][1] = new Vector3( 0.0, 0.0, 0.0);
                ptrBVectorFaceArray[2][i][fieldsGridsSize][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize-i+1][fieldsGridsSize+1][k][2];

        }
    }
    // share with face 4
    for( int j = 1; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
                ptrBVectorFaceArray[2][0][j][k][0] = ptrBVectorFaceArray[4][fieldsGridsSize-j+1][fieldsGridsSize][k][1];
                ptrBVectorFaceArray[2][0][j][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize-j+2][fieldsGridsSize][k][0];
                ptrBVectorFaceArray[2][0][j][k][2] = ptrBVectorFaceArray[4][fieldsGridsSize-j+1][fieldsGridsSize][k][2];
                
                ptrBVectorFaceArray[2][1][j][k][0] = ptrBVectorFaceArray[4][fieldsGridsSize-j+1][fieldsGridsSize+1][k][1];
                ptrBVectorFaceArray[2][1][j][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize-j+2][fieldsGridsSize+1][k][0];
                ptrBVectorFaceArray[2][1][j][k][2] = ptrBVectorFaceArray[4][fieldsGridsSize-j+1][fieldsGridsSize+1][k][2];
        }
    }
    for( int i = 2; i < fieldsGridsSize; i++)
    {
        for( int j = 2; j < fieldsGridsSize; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir < 3; dir++)
                    {
                        ptrBVectorFaceArray[0][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);
                    }
            }
        }
    }
    // face 5
    // share with face 0
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[5][i][fieldsGridsSize+1][k][dir] = ptrBVectorFaceArray[0][i][1][k][dir];
            }
        }
    }

    // share with face 1
    for( int j = 1; j < fieldsGridsSize+1; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            
            ptrBVectorFaceArray[5][fieldsGridsSize+1][j][k][0] = ptrBVectorFaceArray[1][fieldsGridsSize-j+1][1][k][1];
            ptrBVectorFaceArray[5][fieldsGridsSize+1][j][k][1] = ptrBVectorFaceArray[1][fieldsGridsSize-j+2][1][k][0];
            ptrBVectorFaceArray[5][fieldsGridsSize+1][j][k][2] = ptrBVectorFaceArray[1][fieldsGridsSize-j+1][1][k][2];            
            
        }
    }
    // share with face 3
    for( int i = 1; i < fieldsGridsSize+2; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBVectorFaceArray[5][i][0][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize-i+2][1][k][0];
            ptrBVectorFaceArray[5][i][0][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize-i+1][2][k][1];
            ptrBVectorFaceArray[5][i][0][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize-i+1][1][k][2];
        }
    }

    // share with face 4
    for( int j = 0; j < fieldsGridsSize+1; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            
            ptrBVectorFaceArray[5][0][j][k][0] = ptrBVectorFaceArray[4][j][2][k][1];
            ptrBVectorFaceArray[5][0][j][k][1] = ptrBVectorFaceArray[4][j][1][k][0];
            ptrBVectorFaceArray[5][0][j][k][2] = ptrBVectorFaceArray[4][j][1][k][2];            
            
        }
    }
    //
    for( int i = 1; i < fieldsGridsSize+1; i++)
    {
        for( int j = 1; j < fieldsGridsSize+1; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir < 3; dir++)
                    {
                        ptrBVectorFaceArray[0][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);
                    }
            }
        }
    }
    // connect with face 0
    for( int i = 0; i< fieldsGridsSize+2; i++)
    {
        for( int k=0; k< fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[0][i][0][k][dir] = ptrBVectorFaceArray[5][i][fieldsGridsSize][k][dir];
            }
        }
    }
    // connect with face 1
    for( int j = 0; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            
            ptrBVectorFaceArray[1][fieldsGridsSize-j+1][0][k][1] = ptrBVectorFaceArray[5][fieldsGridsSize][j][k][0];
            ptrBVectorFaceArray[1][fieldsGridsSize-j+2][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize][j][k][1];
            ptrBVectorFaceArray[1][fieldsGridsSize-j+1][0][k][2] = ptrBVectorFaceArray[5][fieldsGridsSize][j][k][2];            
            
        }
    }
    // connect with face 3
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            ptrBVectorFaceArray[3][i][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize-i+2][1][k][0];
            ptrBVectorFaceArray[3][i][0][k][1] = ptrBVectorFaceArray[5][fieldsGridsSize-i+1][2][k][1];
            ptrBVectorFaceArray[3][i][0][k][2] = ptrBVectorFaceArray[5][fieldsGridsSize-i+1][1][k][2];
        }
    }
    // connect with face 4
    for( int j = 0; j < fieldsGridsSize+2; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            
            ptrBVectorFaceArray[4][j][0][k][1] = ptrBVectorFaceArray[5][2][j][k][0];
            ptrBVectorFaceArray[4][j][0][k][0] = ptrBVectorFaceArray[5][1][j][k][1];
            ptrBVectorFaceArray[4][j][0][k][2] = ptrBVectorFaceArray[5][1][j][k][2];            
            
        }
    }


    return ptrBVectorFaceArray;

}
*/

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the B on face

// The size of this array is [face * (fsize+2) * (fsize+2) * (fsize) * direction
int ******BVectorFaceArray(int ******ptrBVectorFaceArray)
{
    ptrBVectorFaceArray = new int *****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrBVectorFaceArray[face] = new int ****[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrBVectorFaceArray[face][i] = new int ***[fieldsGridsSize + 2];
            {
                for (int j = 0; j < fieldsGridsSize + 2; j++)
                {
                    ptrBVectorFaceArray[face][i][j] = new int **[fieldsGridsSize * grid_domain];
                    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                    {
                        ptrBVectorFaceArray[face][i][j][k] = new int *[3];
                    }
                }
            }
        }
    }

    // ptr[face][i][j][k][dir]
    // i,j is cell index from 0-fsize+1; k is from 0-(fsize-1)
    // dir is the perticular face index, 0 is perpendicular to the first direction

    // vector initialized face by face
    // face 0 ( to us)
    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[0][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][1] = new int();
        }
    }
    // right
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][0] = new int();
        }
    }
    // face 1( on the right)
    // share with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[1][0][j][k][dir] = ptrBVectorFaceArray[0][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][1][j][k][0] = ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][0];
        }
    }
    // j = 1 rest faces
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][1][j][k][1] = new int();
            ptrBVectorFaceArray[1][1][j][k][2] = new int();
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[1][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][i][fieldsGridsSize + 1][k][1] = new int();
        }
    }
    // right
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][0] = new int();
        }
    }
    //connect face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[1][1][j][k][1];
            ptrBVectorFaceArray[0][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[1][1][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][1][fieldsGridsSize + 1][k][1];
    }

    // face 3 ( on the back)
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[3][0][j][k][dir] = ptrBVectorFaceArray[1][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][1][j][k][0] = ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][0];
        }
    }
    // j = 1 rest faces
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][1][j][k][1] = new int();
            ptrBVectorFaceArray[3][1][j][k][2] = new int();
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[3][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][i][fieldsGridsSize + 1][k][1] = new int();
        }
    }
    // right
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][0] = new int();
        }
    }
    //connect face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[3][1][j][k][1];
            ptrBVectorFaceArray[1][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[3][1][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[3][1][fieldsGridsSize + 1][k][1];
    }
    // face 4 ( on the left)
    // share with face 3
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[4][0][j][k][dir] = ptrBVectorFaceArray[3][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][1][j][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][0];
        }
    }
    // i = 1 rest faces
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][1][j][k][1] = new int();
            ptrBVectorFaceArray[4][1][j][k][2] = new int();
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[4][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][i][fieldsGridsSize + 1][k][1] = new int();
        }
    }
    // right share with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[4][fieldsGridsSize + 1][j][k][dir] = ptrBVectorFaceArray[0][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[0][1][fieldsGridsSize + 1][k][1];
    }
    // connect with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[0][0][j][k][dir] = ptrBVectorFaceArray[4][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }
    // connect with face 3
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[4][1][j][k][1];
            ptrBVectorFaceArray[3][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[4][1][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][1][fieldsGridsSize + 1][k][1];
    }
    // face 2 ( on the top)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[2][i][0][k][dir] = ptrBVectorFaceArray[0][i][fieldsGridsSize][k][dir];
            }

            ptrBVectorFaceArray[2][i][1][k][1] = ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize][k][0];
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][fieldsGridsSize + 1][j][k][0] = ptrBVectorFaceArray[1][j][fieldsGridsSize + 1][k][1];
            ptrBVectorFaceArray[2][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[1][j][fieldsGridsSize][k][0];
            ptrBVectorFaceArray[2][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[1][j][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize][k][0];
    }
    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 2][fieldsGridsSize][k][0];
            ptrBVectorFaceArray[2][i][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][fieldsGridsSize + 1][k][1];
            ptrBVectorFaceArray[2][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[3][1][fieldsGridsSize][k][0];
    }
    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][0][j][k][0] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 1][fieldsGridsSize][k][1];
            ptrBVectorFaceArray[2][0][j][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 2][fieldsGridsSize][k][0];
            ptrBVectorFaceArray[2][0][j][k][2] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 1][fieldsGridsSize][k][2];

            ptrBVectorFaceArray[2][1][j][k][0] = ptrBVectorFaceArray[4][fieldsGridsSize - j + 1][fieldsGridsSize + 1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][1][fieldsGridsSize][k][0];
    }
    // rest of j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][i][1][k][0] = new int();
            ptrBVectorFaceArray[2][i][1][k][2] = new int();
        }
    }
    // rest of i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[2][1][j][k][1] = new int();
            ptrBVectorFaceArray[2][1][j][k][2] = new int();
        }
    }
    // rest of i = 1 && j = 1
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[2][1][1][k][2] = new int();
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[2][i][j][k][dir] = new int();
                }
            }
        }
    }
    // connect with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][i][1][k][0];
            ptrBVectorFaceArray[0][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][i][1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize + 1][1][k][0];
    }
    // connect with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][j][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize][j][k][1];
            ptrBVectorFaceArray[1][j][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][fieldsGridsSize][j][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize][fieldsGridsSize + 1][k][1];
    }

    // connect with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][fieldsGridsSize - i + 2][fieldsGridsSize][k][0];

            ptrBVectorFaceArray[3][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][fieldsGridsSize - i + 1][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][1][fieldsGridsSize][k][0];
    }

    // connect with face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[4][i][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][1][fieldsGridsSize - i + 2][k][1];
            ptrBVectorFaceArray[4][i][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray[2][1][fieldsGridsSize - i + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[2][1][1][k][1];
    }
    // face 5 ( bot)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[5][i][fieldsGridsSize + 1][k][dir] = ptrBVectorFaceArray[0][i][1][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrBVectorFaceArray[0][fieldsGridsSize + 1][1][k][0];
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {

            ptrBVectorFaceArray[5][fieldsGridsSize + 1][j][k][0] = ptrBVectorFaceArray[1][1][fieldsGridsSize - j + 1][k][1];
            ptrBVectorFaceArray[5][fieldsGridsSize + 1][j][k][1] = ptrBVectorFaceArray[1][1][fieldsGridsSize - j + 2][k][0];
            ptrBVectorFaceArray[5][fieldsGridsSize + 1][j][k][2] = ptrBVectorFaceArray[1][1][fieldsGridsSize - j + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[1][1][1][k][0];
    }

    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[5][i][0][k][0] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 2][1][k][0];
            ptrBVectorFaceArray[5][i][0][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][2][k][1];
            ptrBVectorFaceArray[5][i][0][k][2] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][1][k][2];

            ptrBVectorFaceArray[5][i][1][k][1] = ptrBVectorFaceArray[3][fieldsGridsSize - i + 1][1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[3][1][1][k][0];
    }

    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {

            ptrBVectorFaceArray[5][0][j][k][0] = ptrBVectorFaceArray[4][j][2][k][1];
            ptrBVectorFaceArray[5][0][j][k][1] = ptrBVectorFaceArray[4][j][1][k][0];
            ptrBVectorFaceArray[5][0][j][k][2] = ptrBVectorFaceArray[4][j][1][k][2];

            ptrBVectorFaceArray[5][1][j][k][0] = ptrBVectorFaceArray[4][j][1][k][1];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][0][fieldsGridsSize + 1][k][1] = ptrBVectorFaceArray[4][fieldsGridsSize + 1][1][k][0];
    }
    // rest of j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[5][i][1][k][0] = new int();
            ptrBVectorFaceArray[5][i][1][k][2] = new int();
        }
    }
    // rest of i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[5][1][j][k][1] = new int();
            ptrBVectorFaceArray[5][1][j][k][2] = new int();
        }
    }
    // rest of i = 1 && j = 1
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[5][1][1][k][2] = new int();
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrBVectorFaceArray[5][i][j][k][dir] = new int();
                }
            }
        }
    }
    // connect with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrBVectorFaceArray[0][i][0][k][dir] = ptrBVectorFaceArray[5][i][fieldsGridsSize][k][dir];
            }
        }
    }

    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[0][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize][k][0];
    }
    // connect with face 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[1][i][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 2][k][1];
            ptrBVectorFaceArray[1][i][0][k][1] = ptrBVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 1][k][0];
            ptrBVectorFaceArray[1][i][0][k][2] = ptrBVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[1][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize][1][k][1];
    }

    // connect with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrBVectorFaceArray[3][i][0][k][0] = ptrBVectorFaceArray[5][fieldsGridsSize - i + 2][1][k][0];
            ptrBVectorFaceArray[3][i][0][k][1] = ptrBVectorFaceArray[5][fieldsGridsSize - i + 1][2][k][1];
            ptrBVectorFaceArray[3][i][0][k][2] = ptrBVectorFaceArray[5][fieldsGridsSize - i + 1][1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[3][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][1][1][k][0];
    }

    // connect with face 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {

            ptrBVectorFaceArray[4][i][0][k][0] = ptrBVectorFaceArray[5][1][i][k][1];
            ptrBVectorFaceArray[4][i][0][k][1] = ptrBVectorFaceArray[5][2][i][k][0];
            ptrBVectorFaceArray[4][i][0][k][2] = ptrBVectorFaceArray[5][1][i][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrBVectorFaceArray[4][fieldsGridsSize + 1][0][k][0] = ptrBVectorFaceArray[5][1][fieldsGridsSize + 1][k][1];
    }

    return ptrBVectorFaceArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the E on face of dual cell
// The size of this array is [face * (fsize+2) * (fsize+2) * (fsize) * direction
// direction: 0-parallel i, 1-parallel j, 2-parallel k
Vector3 ******EVectorFaceArray(Vector3 ******ptrEVectorFaceArray)
{
    ptrEVectorFaceArray = new Vector3 *****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrEVectorFaceArray[face] = new Vector3 ****[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrEVectorFaceArray[face][i] = new Vector3 ***[fieldsGridsSize + 2];
            {
                for (int j = 0; j < fieldsGridsSize + 2; j++)
                {
                    ptrEVectorFaceArray[face][i][j] = new Vector3 **[fieldsGridsSize * grid_domain];
                    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                    {
                        ptrEVectorFaceArray[face][i][j][k] = new Vector3 *[3];
                    }
                }
            }
        }
    }
    // ptr[face][i][j][k][dir]
    // i, j is the cell index from 0-fsize+1, k is from 0-(fsize-1)

    // initialized face by face
    // face 0 ( to us)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[0][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][0] = new Vector3(0.0, 0.0, 0.0);
            ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }

    // face 1 ( on the right)
    // share with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[1][0][j][k][dir] = ptrEVectorFaceArray[0][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[1][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[1][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }

    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[1][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[1][i][fieldsGridsSize + 1][k][0] = new Vector3(0.0, 0.0, 0.0);
            ptrEVectorFaceArray[1][i][fieldsGridsSize + 1][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // connect with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[0][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[1][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][1][fieldsGridsSize + 1][k][2];
    }

    // face 3 ( on the back)
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[3][0][j][k][dir] = ptrEVectorFaceArray[1][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[3][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[3][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }

    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[3][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[3][i][fieldsGridsSize + 1][k][0] = new Vector3(0.0, 0.0, 0.0);
            ptrEVectorFaceArray[3][i][fieldsGridsSize + 1][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // connect with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[1][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[3][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][1][fieldsGridsSize + 1][k][2];
    }

    // face 4 ( on the left)
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[4][0][j][k][dir] = ptrEVectorFaceArray[3][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[4][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[4][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }

    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[4][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[4][i][fieldsGridsSize + 1][k][0] = new Vector3(0.0, 0.0, 0.0);
            ptrEVectorFaceArray[4][i][fieldsGridsSize + 1][k][2] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // connect 3
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[3][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[4][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[4][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[4][1][fieldsGridsSize + 1][k][2];
    }
    // share 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[4][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[0][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[0][1][fieldsGridsSize + 1][k][0];

        ptrEVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][1][fieldsGridsSize + 1][k][2];
    }

    // connect 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[0][0][j][k][dir] = ptrEVectorFaceArray[4][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }
    // face 2 ( on the top)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[2][i][0][k][dir] = ptrEVectorFaceArray[0][i][fieldsGridsSize][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize][k][2];
    }
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][i][1][k][0] = ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][i][1][k][2] = ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][2];
        }
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][fieldsGridsSize + 1][j][k][0] = ptrEVectorFaceArray[1][j][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][fieldsGridsSize + 1][j][k][1] = ptrEVectorFaceArray[1][j][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][fieldsGridsSize + 1][j][k][2] = ptrEVectorFaceArray[1][j][fieldsGridsSize + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2];
    }
    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][i][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 1][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][i][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][fieldsGridsSize + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[3][1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][1][fieldsGridsSize + 1][k][2];
    }
    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][0][j][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - j + 2][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][0][j][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize - j + 1][fieldsGridsSize][k][0];
            ptrEVectorFaceArray[2][0][j][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - j + 2][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][1][fieldsGridsSize][k][2];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][1][j][k][1] = ptrEVectorFaceArray[0][fieldsGridsSize - j + 1][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][1][j][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize - j + 2][fieldsGridsSize + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++) // not need
    {
        ptrEVectorFaceArray[2][1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][1][fieldsGridsSize + 1][k][2];
    }

    // rest j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][i][1][k][1] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // rest i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][1][j][k][0] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[2][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }
    // connect face 0
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][i][1][k][1];
        }
    }
    // connect face 1
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[1][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][fieldsGridsSize][i][k][0];
        }
    }
    // connect face 3
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[3][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][fieldsGridsSize - i + 2][fieldsGridsSize][k][1];
        }
    }
    // connect face 4
    for (int j = 1; j < fieldsGridsSize + 2; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[4][j][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][1][fieldsGridsSize - j + 2][k][0];
        }
    }
    // face 5 ( on the bot)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[5][i][fieldsGridsSize + 1][k][dir] = ptrEVectorFaceArray[0][i][1][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][1][k][1];
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][1][k][2];
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][fieldsGridsSize + 1][j][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize - j + 2][1][k][1];
            ptrEVectorFaceArray[5][fieldsGridsSize + 1][j][k][1] = ptrEVectorFaceArray[1][fieldsGridsSize - j + 1][1][k][0];
            ptrEVectorFaceArray[5][fieldsGridsSize + 1][j][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize - j + 2][1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][1][1][k][1];
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][1][1][k][2];
    }

    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][i][0][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 1][2][k][0];
            ptrEVectorFaceArray[5][i][0][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][1][k][1];
            ptrEVectorFaceArray[5][i][0][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[3][1][1][k][1];
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[3][1][2][k][2];
    }
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][i][1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 1][1][k][0];
            ptrEVectorFaceArray[5][i][1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][1][k][2];
        }
    }

    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][0][j][k][0] = ptrEVectorFaceArray[4][j][1][k][1];
            ptrEVectorFaceArray[5][0][j][k][1] = ptrEVectorFaceArray[4][j][2][k][0];
            ptrEVectorFaceArray[5][0][j][k][2] = ptrEVectorFaceArray[4][j][2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[4][fieldsGridsSize + 1][1][k][1];
        ptrEVectorFaceArray[5][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[4][fieldsGridsSize + 1][2][k][2];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][1][j][k][1] = ptrEVectorFaceArray[4][j][1][k][0];
            ptrEVectorFaceArray[5][1][j][k][2] = ptrEVectorFaceArray[4][j][1][k][2];
        }
    }
    // rest j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][i][1][k][1] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // rest i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][1][j][k][0] = new Vector3(0.0, 0.0, 0.0);
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[5][i][j][k][dir] = new Vector3(0.0, 0.0, 0.0);
                }
            }
        }
    }

    // connect 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[0][i][0][k][dir] = ptrEVectorFaceArray[5][i][fieldsGridsSize][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize][k][2];
    }
    // connect 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[1][i][0][k][0] = ptrEVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 1][k][1];
            ptrEVectorFaceArray[1][i][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 2][k][0];
            ptrEVectorFaceArray[1][i][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize][1][k][0];
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize][1][k][2];
    }
    // connect 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[3][i][0][k][0] = ptrEVectorFaceArray[5][fieldsGridsSize - i + 1][2][k][0];
            ptrEVectorFaceArray[3][i][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize - i + 2][1][k][1];
            ptrEVectorFaceArray[3][i][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize - i + 2][2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][1][1][k][1];
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][1][2][k][2];
    }
    // connect 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[4][i][0][k][0] = ptrEVectorFaceArray[5][2][i][k][1];
            ptrEVectorFaceArray[4][i][0][k][1] = ptrEVectorFaceArray[5][1][i][k][0];
            ptrEVectorFaceArray[4][i][0][k][2] = ptrEVectorFaceArray[5][2][i][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[4][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[4][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][2][fieldsGridsSize + 1][k][2];
    }

    return ptrEVectorFaceArray;
}
/*
    // face 1 ( on the right)
    // share with face 0
    for( int j =0; j< fieldsGridsSize+1; j++)
    {
        for( int k =0; k < fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir ++)
            {
                ptrEVectorFaceArray[1][0][j][k][dir] = ptrEVectorFaceArray[0][fieldsGridsSize][j][k][dir];
                ptrEVectorFaceArray[1][1][j][k][dir] = ptrEVectorFaceArray[0][fieldsGridsSize+1][j][k][dir];
            }
        }
    }
    for( int i = 2; i < fieldsGridsSize+2; i++)
    {
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[1][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);

                }
            }
        }
    }
    // face 3 ( on the back)
    // share with face 1
    for( int j =0; j< fieldsGridsSize+2; j++)
    {
        for( int k =0; k < fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir ++)
            {
                ptrEVectorFaceArray[3][0][j][k][dir] = ptrEVectorFaceArray[1][fieldsGridsSize][j][k][dir];
                ptrEVectorFaceArray[3][1][j][k][dir] = ptrEVectorFaceArray[1][fieldsGridsSize+1][j][k][dir];
            }
        }
    }
    for( int i = 2; i < fieldsGridsSize+2; i++)
    {
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[3][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);

                }
            }
        }
    }
    // face 4( on the left)
    // share with face 3
    // share with face 0
    for( int j =0; j< fieldsGridsSize+2; j++)
    {
        for( int k =0; k < fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir ++)
            {
                ptrEVectorFaceArray[4][0][j][k][dir] = ptrEVectorFaceArray[3][fieldsGridsSize][j][k][dir];
                ptrEVectorFaceArray[4][1][j][k][dir] = ptrEVectorFaceArray[3][fieldsGridsSize+1][j][k][dir];
                ptrEVectorFaceArray[4][fieldsGridsSize][j][k][dir] = ptrEVectorFaceArray[0][0][j][k][dir];
                ptrEVectorFaceArray[4][fieldsGridsSize+1][j][k][dir] = ptrEVectorFaceArray[0][1][j][k][dir];


            }
        }
    }
    for( int i = 2; i < fieldsGridsSize; i++)
    {
        for( int j = 0; j < fieldsGridsSize+2; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[4][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);

                }
            }
        }
    }
    // face 2 ( on the top)
    // share with face 0
    for( int i = 0; i < fieldsGridsSize+2; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            for( int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[2][i][0][k][dir] = ptrEVectorFaceArray[0][i][fieldsGridsSize][k][dir];
                ptrEVectorFaceArray[2][i][1][k][dir] = ptrEVectorFaceArray[0][i][fieldsGridsSize+1][k][dir];
            }
        }
    }
    // share with face 1
    for( int j = 1; j < fieldsGridsSize+1; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            ptrEVectorFaceArray[2][fieldsGridsSize+1][j][k][0] = ptrEVectorFaceArray[1][j][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][fieldsGridsSize+1][j][k][1] = ptrEVectorFaceArray[1][j][fieldsGridsSize-1][k][0];
            ptrEVectorFaceArray[2][fieldsGridsSize+1][j][k][2] = ptrEVectorFaceArray[1][j][fieldsGridsSize-1][k][2];

            ptrEVectorFaceArray[2][fieldsGridsSize][j][k][0] = ptrEVectorFaceArray[1][j][fieldsGridsSize+1][k][1];
            ptrEVectorFaceArray[2][fieldsGridsSize][j][k][1] = ptrEVectorFaceArray[1][j][fieldsGridsSize][k][0];
            ptrEVectorFaceArray[2][fieldsGridsSize][j][k][2] = ptrEVectorFaceArray[1][j][fieldsGridsSize][k][2];

        }
    }
    // sharewith face 4
    for ( int j = 1 ; j < fieldsGridsSize+1; j++)
    {
        for( int k = 0; k< fieldsGridsSize; k++)
        {
            
            ptrEVectorFaceArray[2][0][j][k][0] = ptrEVectorFaceArray[4][fieldsGridsSize+1-j][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][0][j][k][1] = ptrEVectorFaceArray[4][fieldsGridsSize+1-j+1][fieldsGridsSize][k][0];
            ptrEVectorFaceArray[2][0][j][k][2] = ptrEVectorFaceArray[4][fieldsGridsSize+1-j][fieldsGridsSize][k][2];

            ptrEVectorFaceArray[2][1][j][k][0] = ptrEVectorFaceArray[4][fieldsGridsSize+1-j][fieldsGridsSize+1][k][1];
            ptrEVectorFaceArray[2][1][j][k][1] = ptrEVectorFaceArray[4][fieldsGridsSize+1-j+1][fieldsGridsSize+1][k][0];
            ptrEVectorFaceArray[2][1][j][k][2] = ptrEVectorFaceArray[4][fieldsGridsSize+1-j][fieldsGridsSize+1][k][2];

        }
    }
    // share with face 3
    for( int i = 1; i < fieldsGridsSize; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            ptrEVectorFaceArray[2][i][fieldsGridsSize][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize+1-j][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][i][fieldsGridsSize][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize+1-j-1][fieldsGridsSize+1][k][0];
            ptrEVectorFaceArray[2][i][fieldsGridsSize][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize+1-j][fieldsGridsSize][k][2];

            ptrEVectorFaceArray[2][i][fieldsGridsSize+1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize+1-j][fieldsGridsSize-1][k][1];
            ptrEVectorFaceArray[2][i][fieldsGridsSize+1][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize+1-j+1][fieldsGridsSize-1][k][0];
            ptrEVectorFaceArray[2][i][fieldsGridsSize+1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize+1-j][fieldsGridsSize-1][k][2];

        }
    }
    for( int i = 2; i < fieldsGridsSize; i++)
    {
        for( int j = 2; j< fieldsGridsSize; j++)
        {
            for( int k = 0; k< fieldsGridsSize; k++)
            {
                for( int dir = 0; dir< 3; dir++)
                {
                    ptrEVectorFaceArray[2][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);

                }

            }
        }
    }
    // face 5 ( on the bot)
    // share with face 0
    for( int i = 0; i < fieldsGridsSize+1; i++)
    {
         for( int k = 0; k< fieldsGridsSize; k++)
         {
             for( int dir = 0; dir < 3; dir++)
             {
                 ptrEVectorFaceArray[5][i][fieldsGridsSize][k][dir] = ptrEVectorFaceArray[0][i][0][k][dir];
             }
         }
    }
    // share with face 1
    for( int j = 0; j< fieldsGridsSize; j++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            ptrEVectorFaceArray[5][fieldsGridsSize][j][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize-j][0][k][1];
            ptrEVectorFaceArray[5][fieldsGridsSize][j][k][1] = ptrEVectorFaceArray[1][fieldsGridsSize-j+1][0][k][0];
            ptrEVectorFaceArray[5][fieldsGridsSize][j][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize-j][0][k][2];
        }
    }
    // share with face 4
    for ( int j = 0 ; j < fieldsGridsSize; j++)
    {
        for( int k = 0; k< fieldsGridsSize; k++)
        {
            ptrEVectorFaceArray[5][0][j][k][0] = ptrEVectorFaceArray[4][j][1][k][1];
            ptrEVectorFaceArray[5][0][j][k][1] = ptrEVectorFaceArray[4][j][0][k][0];
            ptrEVectorFaceArray[5][0][j][k][2] = ptrEVectorFaceArray[4][j][0][k][2];
        }
    }
    // share with face 3
    for( int i = 1; i < fieldsGridsSize; i++)
    {
        for( int k = 0; k < fieldsGridsSize; k++)
        {
            ptrEVectorFaceArray[5][i][0][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize-i+1][0][k][0];
            ptrEVectorFaceArray[5][i][0][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize-i][1][k][1];
            ptrEVectorFaceArray[5][i][0][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize-i][0][k][2];
        }
    }
    for( int i = 1; i < fieldsGridsSize; i++)
    {
        for( int j = 1; j < fieldsGridsSize; j++)
        {
            for( int k = 0; k < fieldsGridsSize; k++)
            {
                for( int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[5][i][j][k][dir] = new Vector3( 0.0, 0.0, 0.0);
                }
            }
        }
    }

    return ptrEVectorFaceArray;
    */

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the E on face of dual cell
// The size of this array is [face * (fsize+2) * (fsize+2) * (fsize) * direction
// direction: 0-parallel i, 1-parallel j, 2-parallel k
int ******EVectorFaceArray(int ******ptrEVectorFaceArray)
{
    ptrEVectorFaceArray = new int *****[totalFace];
    for (int face = 0; face < totalFace; face++)
    {
        ptrEVectorFaceArray[face] = new int ****[fieldsGridsSize + 2];
        for (int i = 0; i < fieldsGridsSize + 2; i++)
        {
            ptrEVectorFaceArray[face][i] = new int ***[fieldsGridsSize + 2];
            {
                for (int j = 0; j < fieldsGridsSize + 2; j++)
                {
                    ptrEVectorFaceArray[face][i][j] = new int **[fieldsGridsSize * grid_domain];
                    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                    {
                        ptrEVectorFaceArray[face][i][j][k] = new int *[3];
                    }
                }
            }
        }
    }
    // ptr[face][i][j][k][dir]
    // i, j is the cell index from 0-fsize+1, k is from 0-(fsize-1)

    // initialized face by face
    // face 0 ( to us)
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[0][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][0] = new int();
            ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][2] = new int();
        }
    }

    // face 1 ( on the right)
    // share with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[1][0][j][k][dir] = ptrEVectorFaceArray[0][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[1][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[1][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }

    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[1][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[1][i][fieldsGridsSize + 1][k][0] = new int();
            ptrEVectorFaceArray[1][i][fieldsGridsSize + 1][k][2] = new int();
        }
    }
    // connect with face 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[0][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[1][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][1][fieldsGridsSize + 1][k][2];
    }

    // face 3 ( on the back)
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[3][0][j][k][dir] = ptrEVectorFaceArray[1][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[3][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[3][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }

    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[3][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[3][i][fieldsGridsSize + 1][k][0] = new int();
            ptrEVectorFaceArray[3][i][fieldsGridsSize + 1][k][2] = new int();
        }
    }
    // connect with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[1][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[3][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][1][fieldsGridsSize + 1][k][2];
    }

    // face 4 ( on the left)
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[4][0][j][k][dir] = ptrEVectorFaceArray[3][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[4][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[4][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }

    // main
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 1; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[4][i][j][k][dir] = new int();
                }
            }
        }
    }
    // top
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[4][i][fieldsGridsSize + 1][k][0] = new int();
            ptrEVectorFaceArray[4][i][fieldsGridsSize + 1][k][2] = new int();
        }
    }
    // connect 3
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[3][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[4][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[4][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[4][1][fieldsGridsSize + 1][k][2];
    }
    // share 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[4][fieldsGridsSize + 1][j][k][dir] = ptrEVectorFaceArray[0][1][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[0][1][fieldsGridsSize + 1][k][0];

        ptrEVectorFaceArray[4][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][1][fieldsGridsSize + 1][k][2];
    }

    // connect 0
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[0][0][j][k][dir] = ptrEVectorFaceArray[4][fieldsGridsSize][j][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize][fieldsGridsSize + 1][k][2];
    }
    // face 2 ( on the top)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[2][i][0][k][dir] = ptrEVectorFaceArray[0][i][fieldsGridsSize][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][fieldsGridsSize][k][2];
    }
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][i][1][k][0] = ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][i][1][k][2] = ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][2];
        }
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][fieldsGridsSize + 1][j][k][0] = ptrEVectorFaceArray[1][j][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][fieldsGridsSize + 1][j][k][1] = ptrEVectorFaceArray[1][j][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][fieldsGridsSize + 1][j][k][2] = ptrEVectorFaceArray[1][j][fieldsGridsSize + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2];
    }
    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][i][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 1][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][i][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][fieldsGridsSize + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[3][1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][1][fieldsGridsSize + 1][k][2];
    }
    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][0][j][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - j + 2][fieldsGridsSize][k][1];
            ptrEVectorFaceArray[2][0][j][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize - j + 1][fieldsGridsSize][k][0];
            ptrEVectorFaceArray[2][0][j][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - j + 2][fieldsGridsSize][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[2][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[3][1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[2][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[3][1][fieldsGridsSize][k][2];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][1][j][k][1] = ptrEVectorFaceArray[0][fieldsGridsSize - j + 1][fieldsGridsSize + 1][k][0];
            ptrEVectorFaceArray[2][1][j][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize - j + 2][fieldsGridsSize + 1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++) // not need
    {
        ptrEVectorFaceArray[2][1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][1][fieldsGridsSize + 1][k][2];
    }

    // rest j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][i][1][k][1] = new int();
        }
    }
    // rest i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[2][1][j][k][0] = new int();
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[2][i][j][k][dir] = new int();
                }
            }
        }
    }
    // connect face 0
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[0][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][i][1][k][1];
        }
    }
    // connect face 1
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[1][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][fieldsGridsSize][i][k][0];
        }
    }
    // connect face 3
    for (int i = 1; i < fieldsGridsSize + 2; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[3][i][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][fieldsGridsSize - i + 2][fieldsGridsSize][k][1];
        }
    }
    // connect face 4
    for (int j = 1; j < fieldsGridsSize + 2; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[4][j][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[2][1][fieldsGridsSize - j + 2][k][0];
        }
    }
    // face 5 ( on the bot)
    // share with face 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[5][i][fieldsGridsSize + 1][k][dir] = ptrEVectorFaceArray[0][i][1][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][1] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][1][k][1];
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[0][fieldsGridsSize + 1][1][k][2];
    }
    // share with face 1
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][fieldsGridsSize + 1][j][k][0] = ptrEVectorFaceArray[1][fieldsGridsSize - j + 2][1][k][1];
            ptrEVectorFaceArray[5][fieldsGridsSize + 1][j][k][1] = ptrEVectorFaceArray[1][fieldsGridsSize - j + 1][1][k][0];
            ptrEVectorFaceArray[5][fieldsGridsSize + 1][j][k][2] = ptrEVectorFaceArray[1][fieldsGridsSize - j + 2][1][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[1][1][1][k][1];
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[1][1][1][k][2];
    }

    // share with face 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][i][0][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 1][2][k][0];
            ptrEVectorFaceArray[5][i][0][k][1] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][1][k][1];
            ptrEVectorFaceArray[5][i][0][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[3][1][1][k][1];
        ptrEVectorFaceArray[5][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[3][1][2][k][2];
    }
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][i][1][k][0] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 1][1][k][0];
            ptrEVectorFaceArray[5][i][1][k][2] = ptrEVectorFaceArray[3][fieldsGridsSize - i + 2][1][k][2];
        }
    }

    // share with face 4
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][0][j][k][0] = ptrEVectorFaceArray[4][j][1][k][1];
            ptrEVectorFaceArray[5][0][j][k][1] = ptrEVectorFaceArray[4][j][2][k][0];
            ptrEVectorFaceArray[5][0][j][k][2] = ptrEVectorFaceArray[4][j][2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[5][0][fieldsGridsSize + 1][k][0] = ptrEVectorFaceArray[4][fieldsGridsSize + 1][1][k][1];
        ptrEVectorFaceArray[5][0][fieldsGridsSize + 1][k][2] = ptrEVectorFaceArray[4][fieldsGridsSize + 1][2][k][2];
    }
    for (int j = 1; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][1][j][k][1] = ptrEVectorFaceArray[4][j][1][k][0];
            ptrEVectorFaceArray[5][1][j][k][2] = ptrEVectorFaceArray[4][j][1][k][2];
        }
    }
    // rest j = 1
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][i][1][k][1] = new int();
        }
    }
    // rest i = 1
    for (int j = 2; j < fieldsGridsSize + 1; j++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[5][1][j][k][0] = new int();
        }
    }
    // main
    for (int i = 2; i < fieldsGridsSize + 1; i++)
    {
        for (int j = 2; j < fieldsGridsSize + 1; j++)
        {
            for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    ptrEVectorFaceArray[5][i][j][k][dir] = new int();
                }
            }
        }
    }

    // connect 0
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                ptrEVectorFaceArray[0][i][0][k][dir] = ptrEVectorFaceArray[5][i][fieldsGridsSize][k][dir];
            }
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize][k][1];
        ptrEVectorFaceArray[0][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize + 1][fieldsGridsSize][k][2];
    }
    // connect 1
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[1][i][0][k][0] = ptrEVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 1][k][1];
            ptrEVectorFaceArray[1][i][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 2][k][0];
            ptrEVectorFaceArray[1][i][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize][fieldsGridsSize - i + 2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize][1][k][0];
        ptrEVectorFaceArray[1][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize][1][k][2];
    }
    // connect 3
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[3][i][0][k][0] = ptrEVectorFaceArray[5][fieldsGridsSize - i + 1][2][k][0];
            ptrEVectorFaceArray[3][i][0][k][1] = ptrEVectorFaceArray[5][fieldsGridsSize - i + 2][1][k][1];
            ptrEVectorFaceArray[3][i][0][k][2] = ptrEVectorFaceArray[5][fieldsGridsSize - i + 2][2][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][1][1][k][1];
        ptrEVectorFaceArray[3][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][1][2][k][2];
    }
    // connect 4
    for (int i = 1; i < fieldsGridsSize + 1; i++)
    {
        for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
        {
            ptrEVectorFaceArray[4][i][0][k][0] = ptrEVectorFaceArray[5][2][i][k][1];
            ptrEVectorFaceArray[4][i][0][k][1] = ptrEVectorFaceArray[5][1][i][k][0];
            ptrEVectorFaceArray[4][i][0][k][2] = ptrEVectorFaceArray[5][2][i][k][2];
        }
    }
    for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
    {
        ptrEVectorFaceArray[4][fieldsGridsSize + 1][0][k][1] = ptrEVectorFaceArray[5][1][fieldsGridsSize + 1][k][0];
        ptrEVectorFaceArray[4][fieldsGridsSize + 1][0][k][2] = ptrEVectorFaceArray[5][2][fieldsGridsSize + 1][k][2];
    }

    return ptrEVectorFaceArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Set up a vector array to store the B on face

// The size of this array is [direction * face * (fsize+1) * (fsize+1) * (fsize+1)]
// The size of this array is [direction * face * (fsize+2) * (fsize+2) * (fsize+1)] X

Vector3 *****BVectorFaceArray(GridsPoints *****ptrArray_in)
{
    static Vector3 *mem_BVectorFaceArray = new Vector3[3 * totalFace * (fieldsGridsSize + 1) * (fieldsGridsSize + 1) * (fieldsGridsSize + 1)];
    Vector3 *****ptrBFaceArray = new Vector3 ****[3];
    for (int direction = 0; direction < 3; direction++)
    {
        ptrBFaceArray[direction] = new Vector3 ***[totalFace];
        for (int face = 0; face < totalFace; face++)
        {
            ptrBFaceArray[direction][face] = new Vector3 **[fieldsGridsSize + 1];
            for (int i = 0; i < fieldsGridsSize + 1; i++)
            {
                ptrBFaceArray[direction][face][i] = new Vector3 *[fieldsGridsSize + 1];
                for (int j = 0; j < fieldsGridsSize + 1; j++)
                {
                    //                    ptrBFaceArray[direction][face][i][j] = new Vector3[fieldsGridsSize+1];
                    ptrBFaceArray[direction][face][i][j] = mem_BVectorFaceArray + direction * totalFace * (fieldsGridsSize + 1) * (fieldsGridsSize + 1) * (fieldsGridsSize + 1) + face * (fieldsGridsSize + 1) * (fieldsGridsSize + 1) * (fieldsGridsSize + 1) + i * (fieldsGridsSize + 1) * (fieldsGridsSize + 1) + j * (fieldsGridsSize + 1);
                    for (int k = 0; k < fieldsGridsSize + 1; k++)
                    {
                        int I = i + 1;
                        int J = j + 1;
                        int K = k;
                        Vector3 tempB, tempPos;
                        if (direction == 0) // face perpendicular to i direction
                        {
                            if (k == fieldsGridsSize || j == fieldsGridsSize)
                                continue;
                            tempPos = ptrArray_in[face][I][J][K]->Pos3().PlusProduct(ptrArray_in[face][I][J + 1][K]->Pos3());
                            tempPos = tempPos.PlusProduct(ptrArray_in[face][I][J][K + 1]->Pos3()).PlusProduct(ptrArray_in[face][I][J + 1][K + 1]->Pos3());
                            tempPos = tempPos.ScaleProduct(0.25);
                        }
                        else if (direction == 1) // face perpendicular to j direction
                        {
                            if (k == fieldsGridsSize || i == fieldsGridsSize)
                                continue;
                            tempPos = ptrArray_in[face][I][J][K]->Pos3().PlusProduct(ptrArray_in[face][I + 1][J][K]->Pos3());
                            tempPos = tempPos.PlusProduct(ptrArray_in[face][I][J][K + 1]->Pos3()).PlusProduct(ptrArray_in[face][I + 1][J][K + 1]->Pos3());
                            tempPos = tempPos.ScaleProduct(0.25);
                        }
                        else if (direction == 2) // face perpendicular to k direction
                        {
                            if (i == fieldsGridsSize || j == fieldsGridsSize)
                                continue;
                            tempPos = ptrArray_in[face][I][J][K]->Pos3().PlusProduct(ptrArray_in[face][I + 1][J][K]->Pos3());
                            tempPos = tempPos.PlusProduct(ptrArray_in[face][I][J + 1][K]->Pos3()).PlusProduct(ptrArray_in[face][I + 1][J + 1][K]->Pos3());
                            tempPos = tempPos.ScaleProduct(0.25);
                        }

                        double r = sqrt(pow(tempPos.x(), 2.0) + pow(tempPos.y(), 2.0) + pow(tempPos.z(), 2.0));
                        tempB.Setx(3 * dMoment * tempPos.x() * tempPos.z() / pow(r, 5.0));
                        tempB.Sety(3 * dMoment * tempPos.y() * tempPos.z() / pow(r, 5.0));
                        tempB.Setz(dMoment * (3 * pow(tempPos.z(), 2.0) - pow(r, 2.0)) / pow(r, 5.0));
                        //    ptrBFaceArray[direction][face][i][j][k] = tempB;
                        // dB initialization
                        // BFace only contain dB
                        ptrBFaceArray[direction][face][i][j][k] = {0.0, 0.0, 0.0};
                    }
                }
            }
        }
    }

    return ptrBFaceArray;
}

//************************************************************************
//************************************************************************
// FUNCTION
// Update the B on each face from BVectorFaceArray
// Calculate dB/dt over adjacent faces [madsen1995]
void BVectorFaceArrayUpdate(GridsPoints *****ptrArray,
                            Vector3 ******ptrBVectorFaceArray_main,
                            Vector3 ******ptrBVectorFaceArray_main_backup,
                            int ******ptrBCheckFaceArray_main)
{
    // i direction
    // face i j k
    // adjacent grid points:ptaArray[face][i+1][j+1][k]
    //                      ptaArray[face][i+1][j+2][k]
    //                      ptaArray[face][i+1][j+1][k+1]
    //                      ptaArray[face][i+1][j+2][k+1]
    // adjacent grid points outside:ptaArray[face][i][j+1][k]
    //                              ptaArray[face][i][j+2][k]
    //                              ptaArray[face][i][j+1][k+1]
    //                              ptaArray[face][i][j+2][k+1]
    // adjacent grid points inside: ptaArray[face][i+2][j+1][k]
    //                              ptaArray[face][i+2][j+2][k]
    //                              ptaArray[face][i+2][j+1][k+1]
    //                              ptaArray[face][i+2][j+2][k+1]
    // 20 length with E must be indentified
    //
    //  Using gridspoints E to average length E to simplify the calculation
    // matrix for calculating the face B need:
    // 3 face vector, 3 circle intergration with length E
    //

    // index are for the [fieldgridsize+2][fieldgridsize+2][fieldgridsize+2]
    Vector3 dBOnFace, sumtempB;
    Vector3 AFace, BFace, CFace;
    double weightB, sumWeightB;
    double AInteger, BInteger, CInteger;

    for (int direction = 0; direction < 3; direction++)
    {
        for (int face = 0; face < totalFace; face++)
        {
            // index of main cell No., total range is from [0, fsize+1], [0, ksize-1]
            for (int i = 1; i < fieldsGridsSize + 2; i++)
            {
                for (int j = 1; j < fieldsGridsSize + 2; j++)
                {
                    for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                    {
                        if (direction == 0) // perpendicular to i direction
                        {
                            if (j == fieldsGridsSize + 1)
                                continue;

                            sumtempB = Vector3(0.0, 0.0, 0.0);
                            //face vector
                            AFace = AreaVectorT(ptrArray, face, i - 1, j, k);
                            BFace = AreaVectorF(ptrArray, face, i - 1, j, k);
                            CFace = AreaVectorR(ptrArray, face, i - 1, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationT(ptrArray, face, i - 1, j, k);
                            BInteger = EIntegrationF(ptrArray, face, i - 1, j, k);
                            CInteger = EIntegrationR(ptrArray, face, i - 1, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = dBOnFace.ScaleProduct(weightB);
                            sumWeightB = weightB;

                            AFace = AreaVectorF(ptrArray, face, i - 1, j, k);
                            BFace = AreaVectorBot(ptrArray, face, i - 1, j, k);
                            CFace = AreaVectorR(ptrArray, face, i - 1, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationF(ptrArray, face, i - 1, j, k);
                            BInteger = EIntegrationBot(ptrArray, face, i - 1, j, k);
                            CInteger = EIntegrationR(ptrArray, face, i - 1, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBot(ptrArray, face, i - 1, j, k);
                            BFace = AreaVectorBack(ptrArray, face, i - 1, j, k);
                            CFace = AreaVectorR(ptrArray, face, i - 1, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBot(ptrArray, face, i - 1, j, k);
                            BInteger = EIntegrationBack(ptrArray, face, i - 1, j, k);
                            CInteger = EIntegrationR(ptrArray, face, i - 1, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBack(ptrArray, face, i - 1, j, k);
                            BFace = AreaVectorT(ptrArray, face, i - 1, j, k);
                            CFace = AreaVectorR(ptrArray, face, i - 1, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBack(ptrArray, face, i - 1, j, k);
                            BInteger = EIntegrationT(ptrArray, face, i - 1, j, k);
                            CInteger = EIntegrationR(ptrArray, face, i - 1, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorT(ptrArray, face, i, j, k);
                            BFace = AreaVectorBack(ptrArray, face, i, j, k);
                            CFace = AreaVectorL(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationT(ptrArray, face, i, j, k);
                            BInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            CInteger = EIntegrationL(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBack(ptrArray, face, i, j, k);
                            BFace = AreaVectorBot(ptrArray, face, i, j, k);
                            CFace = AreaVectorL(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            BInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            CInteger = EIntegrationL(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBot(ptrArray, face, i, j, k);
                            BFace = AreaVectorF(ptrArray, face, i, j, k);
                            CFace = AreaVectorL(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            BInteger = EIntegrationF(ptrArray, face, i, j, k);
                            CInteger = EIntegrationL(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorF(ptrArray, face, i, j, k);
                            BFace = AreaVectorT(ptrArray, face, i, j, k);
                            CFace = AreaVectorL(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationF(ptrArray, face, i, j, k);
                            BInteger = EIntegrationT(ptrArray, face, i, j, k);
                            CInteger = EIntegrationL(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            dBOnFace = sumtempB.ScaleProduct(1.0 / sumWeightB);

                            *ptrBVectorFaceArray_main[face][i][j][k][direction] =
                                ptrBVectorFaceArray_main[face][i][j][k][direction]->PlusProduct(dBOnFace.ScaleProduct(tstep / div_max));

                            if (*ptrBCheckFaceArray_main[face][i][j][k][direction] == 0)
                            {
                                *ptrBVectorFaceArray_main_backup[face][i][j][k][direction] =
                                    ptrBVectorFaceArray_main_backup[face][i][j][k][direction]->PlusProduct(
                                        *ptrBVectorFaceArray_main[face][i][j][k][direction]);
                                *ptrBCheckFaceArray_main[face][i][j][k][direction] = 1;
                            }
                        }
                        else if (direction == 1) // perpendicular to j direction
                        {
                            if (i == fieldsGridsSize + 1)
                                continue;

                            sumtempB = Vector3(0.0, 0.0, 0.0);
                            AFace = AreaVectorR(ptrArray, face, i, j - 1, k);
                            BFace = AreaVectorBack(ptrArray, face, i, j - 1, k);
                            CFace = AreaVectorT(ptrArray, face, i, j - 1, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationR(ptrArray, face, i, j - 1, k);
                            BInteger = EIntegrationBack(ptrArray, face, i, j - 1, k);
                            CInteger = EIntegrationT(ptrArray, face, i, j - 1, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = dBOnFace.ScaleProduct(weightB);
                            sumWeightB = weightB;

                            AFace = AreaVectorBack(ptrArray, face, i, j - 1, k);
                            BFace = AreaVectorL(ptrArray, face, i, j - 1, k);
                            CFace = AreaVectorT(ptrArray, face, i, j - 1, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBack(ptrArray, face, i, j - 1, k);
                            BInteger = EIntegrationL(ptrArray, face, i, j - 1, k);
                            CInteger = EIntegrationT(ptrArray, face, i, j - 1, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorL(ptrArray, face, i, j - 1, k);
                            BFace = AreaVectorF(ptrArray, face, i, j - 1, k);
                            CFace = AreaVectorT(ptrArray, face, i, j - 1, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationL(ptrArray, face, i, j - 1, k);
                            BInteger = EIntegrationF(ptrArray, face, i, j - 1, k);
                            CInteger = EIntegrationT(ptrArray, face, i, j - 1, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorF(ptrArray, face, i, j - 1, k);
                            BFace = AreaVectorR(ptrArray, face, i, j - 1, k);
                            CFace = AreaVectorT(ptrArray, face, i, j - 1, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationF(ptrArray, face, i, j - 1, k);
                            BInteger = EIntegrationR(ptrArray, face, i, j - 1, k);
                            CInteger = EIntegrationT(ptrArray, face, i, j - 1, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorR(ptrArray, face, i, j, k);
                            BFace = AreaVectorF(ptrArray, face, i, j, k);
                            CFace = AreaVectorBot(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationR(ptrArray, face, i, j, k);
                            BInteger = EIntegrationF(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorF(ptrArray, face, i, j, k);
                            BFace = AreaVectorL(ptrArray, face, i, j, k);
                            CFace = AreaVectorBot(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationF(ptrArray, face, i, j, k);
                            BInteger = EIntegrationL(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorL(ptrArray, face, i, j, k);
                            BFace = AreaVectorBack(ptrArray, face, i, j, k);
                            CFace = AreaVectorBot(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationL(ptrArray, face, i, j, k);
                            BInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBack(ptrArray, face, i, j, k);
                            BFace = AreaVectorR(ptrArray, face, i, j, k);
                            CFace = AreaVectorBot(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            BInteger = EIntegrationR(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            dBOnFace = sumtempB.ScaleProduct(1.0 / sumWeightB);
                            *ptrBVectorFaceArray_main[face][i][j][k][direction] =
                                ptrBVectorFaceArray_main[face][i][j][k][direction]->PlusProduct(dBOnFace.ScaleProduct(tstep / div_max));

                            if (*ptrBCheckFaceArray_main[face][i][j][k][direction] == 0)
                            {
                                *ptrBVectorFaceArray_main_backup[face][i][j][k][direction] =
                                    ptrBVectorFaceArray_main_backup[face][i][j][k][direction]->PlusProduct(
                                        *ptrBVectorFaceArray_main[face][i][j][k][direction]);
                                *ptrBCheckFaceArray_main[face][i][j][k][direction] = 1;
                            }
                        }
                        else
                        {
                            if (i == fieldsGridsSize + 1 || j == fieldsGridsSize + 1)
                                continue;

                            sumtempB = Vector3(0.0, 0.0, 0.0);
                            AFace = AreaVectorT(ptrArray, face, i, j, k - 1);
                            BFace = AreaVectorL(ptrArray, face, i, j, k - 1);
                            CFace = AreaVectorF(ptrArray, face, i, j, k - 1);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationT(ptrArray, face, i, j, k - 1);
                            BInteger = EIntegrationL(ptrArray, face, i, j, k - 1);
                            CInteger = EIntegrationF(ptrArray, face, i, j, k - 1);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB = weightB;

                            AFace = AreaVectorL(ptrArray, face, i, j, k - 1);
                            BFace = AreaVectorBot(ptrArray, face, i, j, k - 1);
                            CFace = AreaVectorF(ptrArray, face, i, j, k - 1);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationL(ptrArray, face, i, j, k - 1);
                            BInteger = EIntegrationBot(ptrArray, face, i, j, k - 1);
                            CInteger = EIntegrationF(ptrArray, face, i, j, k - 1);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBot(ptrArray, face, i, j, k - 1);
                            BFace = AreaVectorR(ptrArray, face, i, j, k - 1);
                            CFace = AreaVectorF(ptrArray, face, i, j, k - 1);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBot(ptrArray, face, i, j, k - 1);
                            BInteger = EIntegrationR(ptrArray, face, i, j, k - 1);
                            CInteger = EIntegrationF(ptrArray, face, i, j, k - 1);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorR(ptrArray, face, i, j, k - 1);
                            BFace = AreaVectorT(ptrArray, face, i, j, k - 1);
                            CFace = AreaVectorF(ptrArray, face, i, j, k - 1);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationR(ptrArray, face, i, j, k - 1);
                            BInteger = EIntegrationT(ptrArray, face, i, j, k - 1);
                            CInteger = EIntegrationF(ptrArray, face, i, j, k - 1);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorT(ptrArray, face, i, j, k);
                            BFace = AreaVectorR(ptrArray, face, i, j, k);
                            CFace = AreaVectorBack(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationT(ptrArray, face, i, j, k);
                            BInteger = EIntegrationR(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorR(ptrArray, face, i, j, k);
                            BFace = AreaVectorBot(ptrArray, face, i, j, k);
                            CFace = AreaVectorBack(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationR(ptrArray, face, i, j, k);
                            BInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorBot(ptrArray, face, i, j, k);
                            BFace = AreaVectorL(ptrArray, face, i, j, k);
                            CFace = AreaVectorBack(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationBot(ptrArray, face, i, j, k);
                            BInteger = EIntegrationL(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            AFace = AreaVectorL(ptrArray, face, i, j, k);
                            BFace = AreaVectorT(ptrArray, face, i, j, k);
                            CFace = AreaVectorBack(ptrArray, face, i, j, k);
                            weightB = AFace.CrossProduct(BFace).DotProduct(CFace);
                            AInteger = EIntegrationL(ptrArray, face, i, j, k);
                            BInteger = EIntegrationT(ptrArray, face, i, j, k);
                            CInteger = EIntegrationBack(ptrArray, face, i, j, k);
                            dBOnFace.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                            sumtempB = sumtempB.PlusProduct(dBOnFace.ScaleProduct(weightB));
                            sumWeightB += weightB;

                            dBOnFace = sumtempB.ScaleProduct(1.0 / sumWeightB);

                            *ptrBVectorFaceArray_main[face][i][j][k][direction] =
                                ptrBVectorFaceArray_main[face][i][j][k][direction]->PlusProduct(dBOnFace.ScaleProduct(tstep / div_max));

                            if (*ptrBCheckFaceArray_main[face][i][j][k][direction] == 0)
                            {
                                *ptrBVectorFaceArray_main_backup[face][i][j][k][direction] =
                                    ptrBVectorFaceArray_main_backup[face][i][j][k][direction]->PlusProduct(
                                        *ptrBVectorFaceArray_main[face][i][j][k][direction]);
                                *ptrBCheckFaceArray_main[face][i][j][k][direction] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    for (int dir = 0; dir < 3; dir++)
                    {
                        if (dir == 0) // perpendicular to i direction
                        {
                            if (j == fieldsGridsSize + 1)
                                continue;
                            *ptrBCheckFaceArray_main[face][i][j][k][dir] = 0;
                        }
                        else if (dir == 1)
                        {
                            if (i == fieldsGridsSize + 1)
                                continue;
                            *ptrBCheckFaceArray_main[face][i][j][k][dir] = 0;
                        }
                        else if (dir == 2)
                        {
                            if (i == fieldsGridsSize + 1 || j == fieldsGridsSize + 1)
                                continue;
                            *ptrBCheckFaceArray_main[face][i][j][k][dir] = 0;
                        }
                    }
                }
            }
        }
    }
}
//    double Xinteger, sumX;
//    sumX =  EIntegrationT( ptrArray_in, 0, 1, 1, 4) + EIntegrationF( ptrArray_in, 0, 1, 1, 4) +
//            EIntegrationR( ptrArray_in, 0, 1, 1, 4) + EIntegrationBot( ptrArray_in, 0, 1, 1, 4) +
//            EIntegrationBack( ptrArray_in, 0, 1, 1, 4) + EIntegrationL( ptrArray_in, 0, 1, 1, 4);
//    std::cout << " sumX " << sumX << " " <<  EIntegrationT( ptrArray_in, 0, 1, 1, 4) <<
//                                     " " <<  EIntegrationF( ptrArray_in, 0, 1, 1, 4) <<
//                                     " " <<  EIntegrationR( ptrArray_in, 0, 1, 1, 4) <<
//                                     " " <<  EIntegrationBot( ptrArray_in, 0, 1, 1, 4) <<
//                                     " " <<  EIntegrationBack( ptrArray_in, 0, 1, 1, 4) <<
//                                     " " <<  EIntegrationL( ptrArray_in, 0, 1, 1, 4) << std::endl;
//
//    double XVector;
//    XVector =   AreaVectorF( ptrArray_in, 0, 1, 1, 4).NormalizedVector().DotProduct(ptrBFaceArray_in[2][0][1][1][5]) +
//                AreaVectorBack( ptrArray_in, 0, 1, 1, 4).NormalizedVector().DotProduct(ptrBFaceArray_in[2][0][1][1][4]) +
//                AreaVectorT( ptrArray_in, 0, 1, 1, 4).NormalizedVector().DotProduct(ptrBFaceArray_in[1][0][1][2][4]) +
//                AreaVectorBot( ptrArray_in, 0, 1, 1, 4).NormalizedVector().DotProduct(ptrBFaceArray_in[1][0][1][1][4]) +
//                AreaVectorR( ptrArray_in, 0, 1, 1, 4).NormalizedVector().DotProduct(ptrBFaceArray_in[0][0][2][1][4]) +
//                AreaVectorL( ptrArray_in, 0, 1, 1, 4).NormalizedVector().DotProduct(ptrBFaceArray_in[0][0][1][1][4]);
//    std::cout << " XVector " << XVector << std::endl;

// ***************************************************************************************
// Calculate curl B at the center of cells

Vector3 *****CurlBCellArray(GridsPoints *****ptrArray_in,
                            Vector3 *****ptrVectorCellArray,
                            Vector3 *****ptrBVectorFaceArray,
                            double ***ptrVolumeCellArray)
{
    int i, j, k, face_in;
    for (face_in = 0; face_in < totalFace; face_in++)
    {
        // I, J, K is the cell index
        for (int I = 1; I < fieldsGridsSize + 1; I++)
        {
            for (int J = 1; J < fieldsGridsSize + 1; J++)
            {
                for (int K = 0; K < fieldsGridsSize; K++)
                {
                    i = I - 1;
                    j = J - 1;
                    k = K;

                    double volumetemp = ptrVolumeCellArray[I][J][K];
                    if (volumetemp == 0)
                        continue;
                    Vector3 temp = AreaVectorL(ptrArray_in, face_in, I, J, K).CrossProduct(ptrBVectorFaceArray[0][face_in][i][j][k]);

                    temp = temp.PlusProduct(
                        AreaVectorR(ptrArray_in, face_in, I, J, K).CrossProduct(ptrBVectorFaceArray[0][face_in][i + 1][j][k]));

                    temp = temp.PlusProduct(
                        AreaVectorT(ptrArray_in, face_in, I, J, K).CrossProduct(ptrBVectorFaceArray[1][face_in][i][j + 1][k]));

                    temp = temp.PlusProduct(
                        AreaVectorBot(ptrArray_in, face_in, I, J, K).CrossProduct(ptrBVectorFaceArray[1][face_in][i][j][k]));

                    temp = temp.PlusProduct(
                        AreaVectorF(ptrArray_in, face_in, I, J, K).CrossProduct(ptrBVectorFaceArray[2][face_in][i][j][k + 1]));

                    temp = temp.PlusProduct(
                        AreaVectorBack(ptrArray_in, face_in, I, J, K).CrossProduct(ptrBVectorFaceArray[2][face_in][i][j][k]));

                    temp = temp.ScaleProduct(1.0 / volumetemp);
                    ptrVectorCellArray[face_in][I][J][K]->SetVector3(temp);

                    /*    if( temp.x()!=0){
                std::cout << " curlBCell " << std::endl;
                std::cout << i << " " << j << " " << k << std::endl;
                std::cout << " volume " << volumetemp << std::endl;
                std::cout << " curlB " << temp.x() << " " << temp.y() << " " << temp.z() << std::endl;
                int pause;
                std::cin >> pause;
                } 
            */
                }
            }
        }
    }
    return ptrVectorCellArray;
}

// *********************************************************************
// Update E at the grids
void UpdateECellArray(GridsPoints *****ptrArray,
                      Vector3 *****ptrEVectorCellArray,
                      Vector3 *****ptrGradVectorCellArray)
{
    Vector3 tempP;
    double density;
    // i, j, k are index of cells
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 1; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 1; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    density = (ptrArray[face][i][j][k]->Density() +
                               ptrArray[face][i + 1][j][k]->Density() +
                               ptrArray[face][i][j + 1][k]->Density() +
                               ptrArray[face][i + 1][j + 1][k]->Density() +
                               ptrArray[face][i][j][k + 1]->Density() +
                               ptrArray[face][i + 1][j][k + 1]->Density() +
                               ptrArray[face][i][j + 1][k + 1]->Density() +
                               ptrArray[face][i + 1][j + 1][k + 1]->Density()) /
                              8.0;

                    if (density != 0)
                    {
                        tempP = ptrGradVectorCellArray[face][i][j][k]->ScaleProduct(1.0 / qi0 / density);
                    }
                    else
                    {
                        tempP = Vector3{0.0, 0.0, 0.0};
                    }

                    ptrEVectorCellArray[face][i][j][k]->SetVector3(tempP.ScaleProduct(-1.0));
                }
            }
        }
    }
    // update Vel on grids
    // i, j , k are index of grids
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = tempGridsCellLevelBot - coverGridsCellLevelBot + 1; k < fieldsGridsSize - tempGridsCellLevelTop + coverGridsCellLevelTop; k++)
                //    for( int k = 1; k < fieldsGridsSize; k++)
                {
                    Vector3 EVector = Vector3(0.0, 0.0, 0.0);
                    //Vector3 VeleVector = Vector3(0.0, 0.0, 0.0);
                    if (i == 1 && j == 1)
                    {
                        EVector = ptrEVectorCellArray[face][1][0][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][1][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][0][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][0][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][1][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][0][1][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == 1 && j == fieldsGridsSize + 1)
                    {

                        EVector = ptrEVectorCellArray[face][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][1][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][0][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][fieldsGridsSize + 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][fieldsGridsSize][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][0][fieldsGridsSize][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                    {
                        EVector = ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == 1)
                    {
                        EVector = ptrEVectorCellArray[face][fieldsGridsSize][0][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize + 1][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][0][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][1][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][fieldsGridsSize + 1][1][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);
                    }
                    else
                    {
                        EVector = ptrEVectorCellArray[face][i - 1][j - 1][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][i][j - 1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i - 1][j][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i][j][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i - 1][j - 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i][j - 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i - 1][j][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][i][j][k]->V3())
                                      .ScaleProduct(0.1250);
                    }
                    // update E

                    ptrArray[face][i][j][k]->SetE3(EVector);
                    /*              
                    EVector = ptrArray_in[face][i][j][k]->Pos3().NormalizedVector().ScaleProduct(EBot_const);
                    // Apply a external E because of rotation
                    SetRotationalVel(ptrArray_in, face, i, j, k);

                    original_vel = ptrArray_in[face][i][j][k]->Vel_e3();
                    double x_time = (timeline_in*tstep-botBoundaryInitialTimeStart)/botBoundaryInitialTime;
                    externalE = ptrArray_in[face][i][j][k]->B3_base().CrossProduct(original_vel.ScaleProduct(0.50* sin(PI*(x_time-0.50))+0.50));

                    ptrArray_in[face][i][j][k]->SetE3( EVector.PlusProduct(externalE));
*/
                }
            }
        }
    }
}
// *********************************************************************
// Update E at the centers of cells
// with the (curl B), (ve), (grad Pe) at the centers of cells, index from 1 to fsize
// Notice the size of ptrEVectorCellArray is [6][fsize+2][fsize+2][fsize]
void UpdateECellArray(GridsPoints *****ptrArray,
                      Vector3 *****ptrEVectorCellArray,
                      Vector3 *****ptrVeleVectorCellArray,
                      Vector3 *****curlBCellArray,
                      Vector3 *****ptrGradVectorCellArray)
{
    Vector3 Vele, Veli;
    double density;
    Vector3 tempP, tempB;

    // i, j, k are index of cells
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 1; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 1; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    Veli = Vector3(ptrArray[face][i][j][k]->Vel3().x() +
                                       ptrArray[face][i + 1][j][k]->Vel3().x() +
                                       ptrArray[face][i][j + 1][k]->Vel3().x() +
                                       ptrArray[face][i + 1][j + 1][k]->Vel3().x() +
                                       ptrArray[face][i][j][k + 1]->Vel3().x() +
                                       ptrArray[face][i + 1][j][k + 1]->Vel3().x() +
                                       ptrArray[face][i][j + 1][k + 1]->Vel3().x() +
                                       ptrArray[face][i + 1][j + 1][k + 1]->Vel3().x(),
                                   ptrArray[face][i][j][k]->Vel3().y() +
                                       ptrArray[face][i + 1][j][k]->Vel3().y() +
                                       ptrArray[face][i][j + 1][k]->Vel3().y() +
                                       ptrArray[face][i + 1][j + 1][k]->Vel3().y() +
                                       ptrArray[face][i][j][k + 1]->Vel3().y() +
                                       ptrArray[face][i + 1][j][k + 1]->Vel3().y() +
                                       ptrArray[face][i][j + 1][k + 1]->Vel3().y() +
                                       ptrArray[face][i + 1][j + 1][k + 1]->Vel3().y(),
                                   ptrArray[face][i][j][k]->Vel3().z() +
                                       ptrArray[face][i + 1][j][k]->Vel3().z() +
                                       ptrArray[face][i][j + 1][k]->Vel3().z() +
                                       ptrArray[face][i + 1][j + 1][k]->Vel3().z() +
                                       ptrArray[face][i][j][k + 1]->Vel3().z() +
                                       ptrArray[face][i + 1][j][k + 1]->Vel3().z() +
                                       ptrArray[face][i][j + 1][k + 1]->Vel3().z() +
                                       ptrArray[face][i + 1][j + 1][k + 1]->Vel3().z())
                               .ScaleProduct(0.125);

                    density = (ptrArray[face][i][j][k]->Density() +
                               ptrArray[face][i + 1][j][k]->Density() +
                               ptrArray[face][i][j + 1][k]->Density() +
                               ptrArray[face][i + 1][j + 1][k]->Density() +
                               ptrArray[face][i][j][k + 1]->Density() +
                               ptrArray[face][i + 1][j][k + 1]->Density() +
                               ptrArray[face][i][j + 1][k + 1]->Density() +
                               ptrArray[face][i + 1][j + 1][k + 1]->Density()) /
                              8.0;

                    tempB = Vector3(ptrArray[face][i][j][k]->B3_base().x() +
                                        ptrArray[face][i + 1][j][k]->B3_base().x() +
                                        ptrArray[face][i][j + 1][k]->B3_base().x() +
                                        ptrArray[face][i + 1][j + 1][k]->B3_base().x() +
                                        ptrArray[face][i][j][k + 1]->B3_base().x() +
                                        ptrArray[face][i + 1][j][k + 1]->B3_base().x() +
                                        ptrArray[face][i][j + 1][k + 1]->B3_base().x() +
                                        ptrArray[face][i + 1][j + 1][k + 1]->B3_base().x(),
                                    ptrArray[face][i][j][k]->B3_base().y() +
                                        ptrArray[face][i + 1][j][k]->B3_base().y() +
                                        ptrArray[face][i][j + 1][k]->B3_base().y() +
                                        ptrArray[face][i + 1][j + 1][k]->B3_base().y() +
                                        ptrArray[face][i][j][k + 1]->B3_base().y() +
                                        ptrArray[face][i + 1][j][k + 1]->B3_base().y() +
                                        ptrArray[face][i][j + 1][k + 1]->B3_base().y() +
                                        ptrArray[face][i + 1][j + 1][k + 1]->B3_base().y(),
                                    ptrArray[face][i][j][k]->B3_base().z() +
                                        ptrArray[face][i + 1][j][k]->B3_base().z() +
                                        ptrArray[face][i][j + 1][k]->B3_base().z() +
                                        ptrArray[face][i + 1][j + 1][k]->B3_base().z() +
                                        ptrArray[face][i][j][k + 1]->B3_base().z() +
                                        ptrArray[face][i + 1][j][k + 1]->B3_base().z() +
                                        ptrArray[face][i][j + 1][k + 1]->B3_base().z() +
                                        ptrArray[face][i + 1][j + 1][k + 1]->B3_base().z())
                                .ScaleProduct(0.125);

                    if (density != 0)
                    {
                        Vele = Veli.MinusProduct(curlBCellArray[face][i][j][k]->ScaleProduct(1.0 / (mu0 * qi0 * density)));
                        //test for no vi effects
                        //                Vele =  curlBCellArray[face][i][j][k]->ScaleProduct(-1.0 / (mu0 * qi0 *density));

                        /*              if( Vele.x() !=0)
                {        
                std::cout << face << " " << i << " " << j << " " << k << std::endl;
                std::cout << " curlBCellArray " << curlBCellArray[i][j][k].x() << std::endl;
                std::cout << " curl " << curlBCellArray[i][j][k].x() << " " << curlBCellArray[i][j][k].y() << " " << curlBCellArray[i][j][k].z() << " -- ";
                std::cout << " Vi " << Veli.x() << " " << Veli.y() << " " << Veli.z() << " -- " ;
                std::cout << " Ve " << Vele.x() << " " << Vele.y() << " " << Vele.z() << std::endl;
                int pause ;
                std::cin >> pause;
                }
 */
                        tempP = ptrGradVectorCellArray[face][i][j][k]->ScaleProduct(1.0 / qi0 / density);
                    }
                    else
                    {
                        Vele = Vector3{0.0, 0.0, 0.0};
                        tempP = Vector3{0.0, 0.0, 0.0};
                    }

                    ptrVeleVectorCellArray[face][i][j][k]->SetVector3(Vele);

                    //              ptrEVectorCellArray[face][i][j][k]->SetVector3( tempB.CrossProduct(Vele));

                    //  std::cout << " veli" <<  Veli.x() << " " << Veli.y() << " " << Veli.z() << std::endl;
                    // test for no vi effects
                    //                ptrEVectorCellArray[face][i][j][k]->SetVector3( tempP.ScaleProduct(-1.0));
                    ptrEVectorCellArray[face][i][j][k]->SetVector3(tempB.CrossProduct(Vele).PlusProduct(tempP).ScaleProduct(-1.0));

                    /*
                if( tempP.x() == 0)
                {
                    std::cout << i << " " << j << " " << k << std::endl;
                    std::cout <<" tempP in gradPe " <<  tempP.x() << " " << tempP.y() << " " << tempP.z() << std::endl;
                }
                */
                }
            }
        }
    }
    // update Vel on grids
    // i, j , k are index of grids
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = tempGridsCellLevelBot - coverGridsCellLevelBot + 1; k < fieldsGridsSize - tempGridsCellLevelTop + coverGridsCellLevelTop; k++)
                //    for( int k = 1; k < fieldsGridsSize; k++)
                {
                    Vector3 EVector = Vector3(0.0, 0.0, 0.0);
                    Vector3 VeleVector = Vector3(0.0, 0.0, 0.0);
                    if (i == 1 && j == 1)
                    {
                        EVector = ptrEVectorCellArray[face][1][0][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][1][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][0][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][0][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][1][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][0][1][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);

                        VeleVector = ptrVeleVectorCellArray[face][1][0][k - 1]->PlusProduct(
                            ptrVeleVectorCellArray[face][1][1][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][0][1][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][1][0][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][1][1][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                                                   ptrVeleVectorCellArray[face][0][1][k]->V3())
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == 1 && j == fieldsGridsSize + 1)
                    {

                        EVector = ptrEVectorCellArray[face][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][1][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][0][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][fieldsGridsSize + 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][1][fieldsGridsSize][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][0][fieldsGridsSize][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);

                        VeleVector = ptrVeleVectorCellArray[face][1][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            ptrVeleVectorCellArray[face][1][fieldsGridsSize][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][0][fieldsGridsSize][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][1][fieldsGridsSize + 1][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][1][fieldsGridsSize][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                                                   ptrVeleVectorCellArray[face][0][fieldsGridsSize][k]->V3())
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                    {
                        EVector = ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);

                        VeleVector = ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k - 1]->PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize + 1][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize][fieldsGridsSize][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                                                   ptrVeleVectorCellArray[face][fieldsGridsSize + 1][fieldsGridsSize][k]->V3())
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == 1)
                    {
                        EVector = ptrEVectorCellArray[face][fieldsGridsSize][0][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize + 1][1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][0][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][fieldsGridsSize][1][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][fieldsGridsSize + 1][1][k]->V3())
                                      .ScaleProduct(1.0 / 6.0);

                        VeleVector = ptrVeleVectorCellArray[face][fieldsGridsSize][0][k - 1]->PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize][1][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize + 1][1][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize][0][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][fieldsGridsSize][1][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                                                   ptrVeleVectorCellArray[face][fieldsGridsSize + 1][1][k]->V3())
                                         .ScaleProduct(1.0 / 6.0);
                    }
                    else
                    {
                        EVector = ptrEVectorCellArray[face][i - 1][j - 1][k - 1]->PlusProduct(
                            ptrEVectorCellArray[face][i][j - 1][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i - 1][j][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i][j][k - 1]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i - 1][j - 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i][j - 1][k]->V3());
                        EVector = EVector.PlusProduct(
                            ptrEVectorCellArray[face][i - 1][j][k]->V3());
                        EVector = EVector.PlusProduct(
                                             ptrEVectorCellArray[face][i][j][k]->V3())
                                      .ScaleProduct(0.1250);

                        VeleVector = ptrVeleVectorCellArray[face][i - 1][j - 1][k - 1]->PlusProduct(
                            ptrVeleVectorCellArray[face][i][j - 1][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][i - 1][j][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][i][j][k - 1]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][i - 1][j - 1][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][i][j - 1][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                            ptrVeleVectorCellArray[face][i - 1][j][k]->V3());
                        VeleVector = VeleVector.PlusProduct(
                                                   ptrVeleVectorCellArray[face][i][j][k]->V3())
                                         .ScaleProduct(0.1250);
                    }
                    // update E

                    ptrArray[face][i][j][k]->SetE3(EVector);
                    ptrArray[face][i][j][k]->SetVel_e3(VeleVector);
                    //    std::cout << " VeleVector"<< VeleVector.x() << " " << VeleVector.y() << " " << VeleVector.z() << std::endl;

                    //    std::cout << i << " " << j << " " << k << " " << ptrVeleVectorCellArray[i][j][k].x() << " " <<ptrVeleVectorCellArray[i][j][k].y() << " " << ptrVeleVectorCellArray[i][j][k].z() << std::endl;
                    //        std::cout << ptrArray[face][i][j][k]->Vel_e3().x() << " " << ptrArray[face][i][j][k]->Vel_e3().y() << " " << ptrArray[face][i][j][k]->Vel_e3().z() << std::endl;
                    //        std::cout << " EVector"<< EVector.x() << " " << EVector.y() << " " << EVector.z() << std::endl;
                    //        std::cout << ptrArray[face][i][j][k]->E3().x() << " " <<ptrArray[face][i][j][k]->E3().y() << " " <<ptrArray[face][i][j][k]->E3().z() << std::endl << std::endl;
                }
            }
        }
    }
}
//*******************************************************
// update E on each grids
// Using Vector3****** ptrEFace_dual_backup
void EVectorGridsArrayUpdate(GridsPoints *****ptrArray,
                             Vector3 ******ptrEFace_dual)
{
    Vector3 tempCurl = Vector3(0.0, 0.0, 0.0);
    Vector3 ve = Vector3(0.0, 0.0, 0.0);
    // i, j, k are the index of main cell
    // i, j range are 0-fsize+1
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    if (i == 1 && j == 1)
                    {
                        tempCurl = ptrEFace_dual[face][i][j][k][0]->PlusProduct(
                            *ptrEFace_dual[face][i][j - 1][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k - 1][2]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][2]);

                        tempCurl = tempCurl.ScaleProduct(0.20);
                    }
                    else if (i == 1 && j == fieldsGridsSize + 1)
                    {
                        tempCurl = ptrEFace_dual[face][i][j][k][0]->PlusProduct(
                            *ptrEFace_dual[face][i][j - 1][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k - 1][2]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][2]);

                        tempCurl = tempCurl.ScaleProduct(0.20);
                    }
                    else if (i == fieldsGridsSize + 1 && j == 1)
                    {
                        tempCurl = ptrEFace_dual[face][i - 1][j][k][0]->PlusProduct(
                            *ptrEFace_dual[face][i][j - 1][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k - 1][2]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][2]);

                        tempCurl = tempCurl.ScaleProduct(0.20);
                    }
                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                    {
                        tempCurl = ptrEFace_dual[face][i - 1][j][k][0]->PlusProduct(
                            *ptrEFace_dual[face][i][j - 1][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k - 1][2]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][2]);

                        tempCurl = tempCurl.ScaleProduct(0.20);
                    }
                    else
                    {
                        tempCurl = ptrEFace_dual[face][i - 1][j][k][0]->PlusProduct(
                            *ptrEFace_dual[face][i][j][k][0]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j - 1][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][1]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k - 1][2]);
                        tempCurl = tempCurl.PlusProduct(
                            *ptrEFace_dual[face][i][j][k][2]);

                        tempCurl = tempCurl.ScaleProduct(1.0 / 6.0);
                    }

                    ve = ptrArray[face][i][j][k]->Vel3().MinusProduct(tempCurl.ScaleProduct(1.0 / div_max / qi0 / mu0 / ptrArray[face][i][j][k]->Density()));

                    //        tempCurl = ve.CrossProduct(ptrArray[face][i][j][k]->B3()).PlusProduct(ptrArray[face][i][j][k]->GradPe().ScaleProduct(1.0/qi0/ptrArray[face][i][j][k]->Density())).ScaleProduct(-1.0);
                    tempCurl = ve.CrossProduct(ptrArray[face][i][j][k]->B3());

                    ptrArray[face][i][j][k]->SetE3(tempCurl);
                }
            }
        }
    }

    // reset backup for main array
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    for (int dir = 0; dir < 3; dir++)
                    {
                        if (dir == 0) // perpendicular to i direction
                        {
                            if (j == fieldsGridsSize + 1)
                                continue;
                            *ptrEFace_dual[face][i][j][k][dir] = Vector3{0.0, 0.0, 0.0};
                        }
                        else if (dir == 1)
                        {
                            if (i == fieldsGridsSize + 1)
                                continue;
                            *ptrEFace_dual[face][i][j][k][dir] = Vector3{0.0, 0.0, 0.0};
                        }
                        else if (dir == 2)
                        {
                            if (i == fieldsGridsSize + 1 || j == fieldsGridsSize + 1)
                                continue;
                            *ptrEFace_dual[face][i][j][k][dir] = Vector3{0.0, 0.0, 0.0};
                        }
                    }
                }
            }
        }
    }
}

// ******************************************************
// Update the dB on each grids
// Using  Vector3****** ptrBVectorFaceArray_backup
void BVectorGridsArrayUpdate(GridsPoints *****ptrArray,
                             Vector3 ******ptrBVectorFaceArray)
{
    Vector3 Btemp = Vector3(0.0, 0.0, 0.0);
    for (int face = 0; face < totalFace; face++)
    {
        // main cell index, the total range is from 0 - fsize+1, [0, fsize-1]
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {

                    if (i == 1 && (j == 1))
                    {
                        Btemp = ptrBVectorFaceArray[face][i][j][k][0]->PlusProduct(
                                                                         *ptrBVectorFaceArray[face][i][j - 1][k][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k][2])
                                    .ScaleProduct(1.0 / 9.0);
                    }
                    else if (i == 1 && j == fieldsGridsSize + 1)
                    {
                        Btemp = ptrBVectorFaceArray[face][i][j][k][0]->PlusProduct(
                                                                         *ptrBVectorFaceArray[face][i][j - 1][k][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j - 1][k][2])
                                    .ScaleProduct(1.0 / 9.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == 1)
                    {
                        Btemp = ptrBVectorFaceArray[face][i][j][k][0]->PlusProduct(
                                                                         *ptrBVectorFaceArray[face][i][j - 1][k][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k - 1][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j - 1][k][2])
                                    .ScaleProduct(1.0 / 9.0);
                    }
                    else if (i == fieldsGridsSize + 1 && j == fieldsGridsSize + 1)
                    {
                        Btemp = ptrBVectorFaceArray[face][i][j][k][0]->PlusProduct(
                                                                         *ptrBVectorFaceArray[face][i][j - 1][k][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k - 1][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j - 1][k][2])
                                    .ScaleProduct(1.0 / 9.0);
                    }
                    else
                    {
                        Btemp = ptrBVectorFaceArray[face][i][j][k][0]->PlusProduct(
                                                                         *ptrBVectorFaceArray[face][i][j - 1][k][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k - 1][0])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k - 1][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k - 1][1])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i][j - 1][k][2])
                                    .PlusProduct(
                                        *ptrBVectorFaceArray[face][i - 1][j - 1][k][2])
                                    .ScaleProduct(1.0 / 12.0);
                    }

                    double average = 1.0 / div_max;
                    ptrArray[face][i][j][k]->SetdB3(Btemp.ScaleProduct(average));
                }
            }
        }
    }
    // reset backup for main array
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 1; i < fieldsGridsSize + 2; i++)
        {
            for (int j = 1; j < fieldsGridsSize + 2; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain; k++)
                {
                    for (int dir = 0; dir < 3; dir++)
                    {
                        if (dir == 0) // perpendicular to i direction
                        {
                            if (j == fieldsGridsSize + 1)
                                continue;
                            *ptrBVectorFaceArray[face][i][j][k][dir] = Vector3{0.0, 0.0, 0.0};
                        }
                        else if (dir == 1)
                        {
                            if (i == fieldsGridsSize + 1)
                                continue;
                            *ptrBVectorFaceArray[face][i][j][k][dir] = Vector3{0.0, 0.0, 0.0};
                        }
                        else if (dir == 2)
                        {
                            if (i == fieldsGridsSize + 1 || j == fieldsGridsSize + 1)
                                continue;
                            *ptrBVectorFaceArray[face][i][j][k][dir] = Vector3{0.0, 0.0, 0.0};
                        }
                    }
                }
            }
        }
    }
}

// calculate curl B on face of E dual cell ( current j)
void EVectorFaceArrayUpdate(Vector3 *****ptrArray_dual,
                            Vector3 ******ptrBVectorFaceArray_main,
                            Vector3 ******ptrEVectorFaceArray_dual,
                            Vector3 ******ptrEVectorFaceArray_dual_backup,
                            int ******ptrECheckFaceArray_dual)
{
    Vector3 curldB, sumtemp;
    Vector3 AFace, BFace, CFace;
    //double weight = 0.0;
    //double sumWeight = 0.0;
    //double AInteger = 0.0;
    //double BInteger = 0.0;
    //double CInteger = 0.0;

    for (int dir = 0; dir < 3; dir++)
    {
        for (int face = 0; face < totalFace; face++)
        {
            if (dir == 0) // i direction
            {
                for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                {
                    ptrArray_dual[face][0][0][k] = ptrArray_dual[face][0][1][k];
                    ptrArray_dual[face][0][fieldsGridsSize + 1][k] = ptrArray_dual[face][0][fieldsGridsSize][k];
                    ptrArray_dual[face][fieldsGridsSize + 1][0][k] = ptrArray_dual[face][fieldsGridsSize + 1][1][k];
                    ptrArray_dual[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k] = ptrArray_dual[face][fieldsGridsSize + 1][fieldsGridsSize][k];

                    ptrBVectorFaceArray_main[face][0][0][k][2] = ptrBVectorFaceArray_main[face][0][1][k][2];
                    ptrBVectorFaceArray_main[face][0][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray_main[face][0][fieldsGridsSize][k][2];
                    ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][0][k][2] = ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][1][k][2];
                    ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][fieldsGridsSize][k][2];
                }
            }
            else if (dir == 1) // j direction
            {
                for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                {
                    ptrArray_dual[face][0][0][k] = ptrArray_dual[face][1][0][k];
                    ptrArray_dual[face][0][fieldsGridsSize + 1][k] = ptrArray_dual[face][1][fieldsGridsSize + 1][k];
                    ptrArray_dual[face][fieldsGridsSize + 1][0][k] = ptrArray_dual[face][fieldsGridsSize][0][k];
                    ptrArray_dual[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k] = ptrArray_dual[face][fieldsGridsSize][fieldsGridsSize + 1][k];

                    ptrBVectorFaceArray_main[face][0][0][k][2] = ptrBVectorFaceArray_main[face][1][0][k][2];
                    ptrBVectorFaceArray_main[face][0][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray_main[face][1][fieldsGridsSize + 1][k][2];
                    ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][0][k][2] = ptrBVectorFaceArray_main[face][fieldsGridsSize][0][k][2];
                    ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray_main[face][fieldsGridsSize][fieldsGridsSize + 1][k][2];
                }
            }
            else if (dir == 2) // left two corner have no left face, while right two corner have no right face
            {
                for (int k = 0; k < fieldsGridsSize * grid_domain; k++)
                {
                    ptrArray_dual[face][0][0][k] = ptrArray_dual[face][0][1][k];
                    ptrArray_dual[face][0][fieldsGridsSize + 1][k] = ptrArray_dual[face][0][fieldsGridsSize][k];
                    ptrArray_dual[face][fieldsGridsSize + 1][0][k] = ptrArray_dual[face][fieldsGridsSize + 1][1][k];
                    ptrArray_dual[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k] = ptrArray_dual[face][fieldsGridsSize + 1][fieldsGridsSize][k];

                    ptrBVectorFaceArray_main[face][0][0][k][2] = ptrBVectorFaceArray_main[face][0][1][k][2];
                    ptrBVectorFaceArray_main[face][0][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray_main[face][0][fieldsGridsSize][k][2];
                    ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][0][k][2] = ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][1][k][2];
                    ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][fieldsGridsSize + 1][k][2] = ptrBVectorFaceArray_main[face][fieldsGridsSize + 1][fieldsGridsSize][k][2];
                }
            }
            // index of dual cell, k is from 1 - (fSize-1)
            //
            for (int i = 0; i < fieldsGridsSize + 1; i++)
            {
                for (int j = 0; j < fieldsGridsSize + 1; j++)
                {
                    for (int k = 1; k < fieldsGridsSize * grid_domain; k++) // notice k starts at 1
                    {

                        if (dir == 0 && i > 0 && k < fieldsGridsSize * grid_domain - 1) // i direction
                        {
                            //               std::cout << " d0 " << face << " " << i << " " << j << " " << k << std::endl;
                            *ptrEVectorFaceArray_dual[face][i][j + 1][k + 1][0] = CurldBOnDualCellFace(ptrArray_dual,
                                                                                                       ptrBVectorFaceArray_main,
                                                                                                       0,
                                                                                                       face,
                                                                                                       i,
                                                                                                       j,
                                                                                                       k);
                            //    if( k == fieldsGridsSize*grid_domain-1)
                            //    std::cout << " E face check " << ptrEVectorFaceArray_dual[face][i][j+1][k+1][0]->norm() << "\n";

                            if (*ptrECheckFaceArray_dual[face][i][j + 1][k + 1][0] == 0)
                            {
                                *ptrEVectorFaceArray_dual_backup[face][i][j + 1][k + 1][0] = ptrEVectorFaceArray_dual_backup[face][i][j + 1][k + 1][0]->PlusProduct(
                                    *ptrEVectorFaceArray_dual[face][i][j + 1][k + 1][0]);
                                *ptrECheckFaceArray_dual[face][i][j + 1][k + 1][0] = 1;
                            }
                        }
                        else if (dir == 1 && j > 0 && k < fieldsGridsSize * grid_domain - 1) // j direction
                        {
                            *ptrEVectorFaceArray_dual[face][i + 1][j][k + 1][1] = CurldBOnDualCellFace(ptrArray_dual,
                                                                                                       ptrBVectorFaceArray_main,
                                                                                                       1,
                                                                                                       face,
                                                                                                       i,
                                                                                                       j,
                                                                                                       k);
                            if (*ptrECheckFaceArray_dual[face][i + 1][j][k + 1][1] == 0)
                            {
                                *ptrEVectorFaceArray_dual_backup[face][i + 1][j][k + 1][1] =
                                    ptrEVectorFaceArray_dual_backup[face][i + 1][j][k + 1][1]->PlusProduct(
                                        *ptrEVectorFaceArray_dual[face][i + 1][j][k + 1][1]);
                                *ptrECheckFaceArray_dual[face][i + 1][j][k + 1][1] = 1;
                            }
                        }
                        else if (dir == 2) // k direction
                        {
                            *ptrEVectorFaceArray_dual[face][i + 1][j + 1][k][2] = CurldBOnDualCellFace(ptrArray_dual,
                                                                                                       ptrBVectorFaceArray_main,
                                                                                                       2,
                                                                                                       face,
                                                                                                       i,
                                                                                                       j,
                                                                                                       k);
                            if (*ptrECheckFaceArray_dual[face][i + 1][j + 1][k][2] == 0)
                            {
                                *ptrEVectorFaceArray_dual_backup[face][i + 1][j + 1][k][2] =
                                    ptrEVectorFaceArray_dual_backup[face][i + 1][j + 1][k][2]->PlusProduct(
                                        *ptrEVectorFaceArray_dual[face][i + 1][j + 1][k][2]);
                                *ptrECheckFaceArray_dual[face][i + 1][j + 1][k][2] = 1;
                            }
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
            }
        }
    }

    // reset check 0 for dual array
    for (int face = 0; face < totalFace; face++)
    {
        for (int i = 0; i < fieldsGridsSize + 1; i++)
        {
            for (int j = 0; j < fieldsGridsSize + 1; j++)
            {
                for (int k = 1; k < fieldsGridsSize * grid_domain - 1; k++)
                {
                    for (int dir = 0; dir < 3; dir++)
                    {
                        if (dir == 0 && i > 0) // i direction
                        {
                            *ptrECheckFaceArray_dual[face][i][j + 1][k + 1][0] = 0;
                        }
                        else if (dir == 1 && j > 0)
                        {
                            *ptrECheckFaceArray_dual[face][i + 1][j][k + 1][1] = 0;
                        }
                        else if (dir == 2)
                        {
                            *ptrECheckFaceArray_dual[face][i + 1][j + 1][k][2] = 0;
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
            }
        }
    }
}
// calculate curl dB
// on the length
// to be consistant with the function
// i j k are the dual cell index
// (i j k 0) in dual-cell is (i, j+1, k+1, 0) in length-index
// (i j k 1) in dual-cell is (i+1, j, k+1, 1) in length-index
// (i j k 2) in dual-cell is (i+1, j+1, k, 2) in length-index

// (i, i-1), (j, j-1), (k, k-1)
// range: i = 0 - fsize+1; j = 0 - fsize+1; k = 0 - fsize-2
// value range: i = 1 - fsize; j = 1 - fsize; k = 1 - fsize-2
Vector3 CurldBOnDualCellFace(Vector3 *****ptrArray_dual,
                             Vector3 ******ptrBVectorFaceArray_main,
                             int dir,
                             int face,
                             int i,
                             int j,
                             int k)
{
    Vector3 curldB, sumtemp;
    Vector3 AFace, BFace, CFace;
    double weight = 0.0;
    double sumWeight = 0.0;
    double AInteger = 0.0;
    double BInteger = 0.0;
    double CInteger = 0.0;

    if (dir == 0)
    {
        // face vector
        AFace = AreaVectorT(ptrArray_dual, face, i - 1, j, k);
        BFace = AreaVectorF(ptrArray_dual, face, i - 1, j, k);
        CFace = AreaVectorR(ptrArray_dual, face, i - 1, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        BInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        CInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight = weight;

        AFace = AreaVectorF(ptrArray_dual, face, i - 1, j, k);
        BFace = AreaVectorBot(ptrArray_dual, face, i - 1, j, k);
        CFace = AreaVectorR(ptrArray_dual, face, i - 1, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        CInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        AFace = AreaVectorBot(ptrArray_dual, face, i - 1, j, k);
        BFace = AreaVectorBack(ptrArray_dual, face, i - 1, j, k);
        CFace = AreaVectorR(ptrArray_dual, face, i - 1, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        BInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        CInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        AFace = AreaVectorBack(ptrArray_dual, face, i - 1, j, k);
        BFace = AreaVectorT(ptrArray_dual, face, i - 1, j, k);
        CFace = AreaVectorR(ptrArray_dual, face, i - 1, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        CInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i - 1, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        AFace = AreaVectorT(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorL(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        AFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorL(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        AFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorF(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorL(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        AFace = AreaVectorF(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorT(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorL(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight += weight;

        curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
    }
    else if (dir == 1)
    {
        AFace = AreaVectorR(ptrArray_dual, face, i, j - 1, k);
        BFace = AreaVectorBack(ptrArray_dual, face, i, j - 1, k);
        CFace = AreaVectorT(ptrArray_dual, face, i, j - 1, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        BInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        CInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = curldB.ScaleProduct(weight);
        sumWeight = weight;

        AFace = AreaVectorBack(ptrArray_dual, face, i, j - 1, k);
        BFace = AreaVectorL(ptrArray_dual, face, i, j - 1, k);
        CFace = AreaVectorT(ptrArray_dual, face, i, j - 1, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        CInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;

        AFace = AreaVectorL(ptrArray_dual, face, i, j - 1, k);
        BFace = AreaVectorF(ptrArray_dual, face, i, j - 1, k);
        CFace = AreaVectorT(ptrArray_dual, face, i, j - 1, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        BInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        CInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;

        AFace = AreaVectorF(ptrArray_dual, face, i, j - 1, k);
        BFace = AreaVectorR(ptrArray_dual, face, i, j - 1, k);
        CFace = AreaVectorT(ptrArray_dual, face, i, j - 1, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        CInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j - 1, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;

        AFace = AreaVectorR(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorF(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;

        AFace = AreaVectorF(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorL(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;

        AFace = AreaVectorL(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;

        AFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
        BFace = AreaVectorR(ptrArray_dual, face, i, j, k);
        CFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
        weight = AFace.CrossProduct(BFace).DotProduct(CFace);
        AInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        CInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
        curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
        sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
        sumWeight += weight;
        curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
    }
    else if (dir == 2)
    {
        if (k == fieldsGridsSize - 1)
        {
            if (i == 0 && (j == fieldsGridsSize || j == 0)) // left two corner, dual cell have no left face
            {
                AFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = curldB.ScaleProduct(weight);
                sumWeight = weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
            }
            else if (i == fieldsGridsSize && (j == fieldsGridsSize || j == 0)) // right two corner, dual cell have no right face
            {
                AFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = curldB.ScaleProduct(weight);
                sumWeight = weight;

                AFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
            }
            else
            {

                AFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = curldB.ScaleProduct(weight);
                sumWeight = weight;

                AFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
            }
        }
        else // != fieldsGridsSize
        {
            if (i == 0 && (j == fieldsGridsSize || j == 0)) // left two corner, dual cell have no left face
            {
                AFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = curldB.ScaleProduct(weight);
                sumWeight = weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorT(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorR(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorR(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
            }
            else if (i == fieldsGridsSize && (j == fieldsGridsSize || j == 0)) // right two corner, dual cell have no right face
            {
                AFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = curldB.ScaleProduct(weight);
                sumWeight = weight;

                AFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorT(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorL(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorL(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
            }
            else
            {

                AFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = curldB.ScaleProduct(weight);
                sumWeight = weight;

                AFace = AreaVectorL(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorR(ptrArray_dual, face, i, j, k - 1);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k - 1);
                CFace = AreaVectorF(ptrArray_dual, face, i, j, k - 1);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                CInteger = dBIntegrationF(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k - 1);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorT(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorR(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorR(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationR(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorBot(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorL(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationBot(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                AFace = AreaVectorL(ptrArray_dual, face, i, j, k);
                BFace = AreaVectorT(ptrArray_dual, face, i, j, k);
                CFace = AreaVectorBack(ptrArray_dual, face, i, j, k);
                weight = AFace.CrossProduct(BFace).DotProduct(CFace);
                AInteger = dBIntegrationL(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                BInteger = dBIntegrationT(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                CInteger = dBIntegrationBack(ptrBVectorFaceArray_main, ptrArray_dual, face, i, j, k);
                curldB.FaceBSolver(AFace, BFace, CFace, AInteger, BInteger, CInteger);
                sumtemp = sumtemp.PlusProduct(curldB.ScaleProduct(weight));
                sumWeight += weight;

                curldB = sumtemp.ScaleProduct(1.0 / sumWeight);
            }
        }
    }

    return curldB;
}
