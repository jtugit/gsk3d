#ifndef _MODULE_H_
#define _MODULE_H_
#include <iostream>
#include "parameters.h"
#include <memory>
#include "particles.h"
#include <vector>
#include "fieldsgrids.h"
#include "module_0.h"
#include "module_1.h"
#include <cmath>
#include <limits>
#include <bitset>
#include <omp.h>
#include "module_base.h"
#include "gridscells.h"
#include <fstream>
#include <vector>
#include "geopack.h"
using std::cout;
using std::endl;
using std::string;
using std::vector;
//************************************************************************
//************************************************************************
// FUNCTION
// Main general control funtion
//************************************************************************
//************************************************************************
void ProcessFunc(double, double, double, double, double, double, double, string);

//************************************************************************
//************************************************************************
// FUNCTION
// test
//
// test vector3 trans to uint_64
inline void TestTransV3_uint64()
{
    uint_64 intPos = UniDisInCell(4,15,3,0);
    Vector3 tempV3 = Uint64ToVector3(intPos);
    Vector3 vel = { 1000.0, 100.0, 10.0};
    Particles tempPar = Particles(intPos, tempV3, vel, 0.0, 0.0);
    std::cout << "a= "<< tempPar.PosUint() << " \n" ;
    tempPar.UpdateUint_64();
    std::cout << "b= "<< tempPar.PosUint() << " \n" ;
    tempPar.UpdateUint_64();
    std::cout << "c= "<< tempPar.PosUint() << " \n" ;

    std::cout << " TransTest end\n";
    std::cin.get();

}
// test filename coding
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
inline void TestFileName(GridsPoints *****ptrArrayGrids, int timeline)
{
    //double pi = acos(-1.0);
    string iriDataPath = "./ExternalData/INPUT/IRI.txt";
    string tiegcmDataPath = "./ExternalData/INPUT/TIEGCM.txt";
    std::ifstream iriFile, tiegcmFile;
    //
    iriFile.open(iriDataPath, std::ios::in);
    tiegcmFile.open(tiegcmDataPath, std::ios::in);
    if (!iriFile || !tiegcmFile)
    {
        std::cout << " unable to open data file";
        exit(1);
    }
    else
    {
        std::cout << " Open data sucessfully \n";
    }
    // apply continue memory to store IRI data
    // iri: Ti, H+, He+, O+
    vector<vector<vector<double>>> *iriData = new vector<vector<vector<double>>>(201, vector<vector<double>>(201, vector<double>(4)));
    int rowIri = 201;
    int columnIri = 201;
    //char bufferIri[2048];
    //double tempIri[4];
    double iri[17];
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
                   &tempIri[0],            // Ti
                   &tempIri[1],            //NO+
                   &tempIri[2],            // NH+
                   &tempIri[3]);           // NHe+
            (*iriData)[i][j][0] = tempIri[0]; // Ti
            (*iriData)[i][j][1] = tempIri[2]; // NH+
            (*iriData)[i][j][2] = tempIri[3]; // NHe+
            (*iriData)[i][j][3] = tempIri[1]; // NO+*/
        }
    }
    // TIEGCM
    vector<vector<vector<double>>> *tiegcmData = new vector<vector<vector<double>>>(38, vector<vector<double>>(73, vector<double>(4)));
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
    // (xyz) in GM -> (long,lat) in GEO
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
    for( int f = 0; f < totalFace; f++)
    {
        for( int i = 1; i < fieldsGridsSize +2; i++)
        {
            for( int j = 1; j < fieldsGridsSize +2; j++)
            {
                for(k = 0; k < fieldsGridsSize; k++)
                {
                    // get (xmag, ymag, zmag)
                    x = ptrArrayGrids[f][i][j][k]->Pos3().x();
                    y = ptrArrayGrids[f][i][j][k]->Pos3().y();
                    z = ptrArrayGrids[f][i][j][k]->Pos3().z();
                    // trans into (xmag, ymag, zmag) on earth at 500km  // why 750 km does not work
                    DipoleRelatedOnEarthPolar( x, y, z, xmag, ymag, zmag, 700000.0);
                    // trans into (xgeo, ygeo, zgeo)
                    GEOMAG_08(xgeo, ygeo, zgeo, xmag, ymag, zmag, -1);
                    // trans into (thet, phis) in GEO
                    SPHCAR_08(ra, thet, phis, xgeo, ygeo, zgeo, -1);
                    // bottom layer, for only one layer profile, we set k =0 and 1 are the same density
                    if( k <= tempGridsCellLevelBot)
                    {
                        // locate the index in IRI and TIEGCM array by "thet" and "phis", and do a rough interpolation
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
                        ti = (*iriData)[a][b][0] * w1 + (*iriData)[a][b + 1][0] *w3+ (*iriData)[a + 1][b][0] * w2+ (*iriData)[a + 1][b + 1][0] *w4;
                        nH = (*iriData)[a][b][1] * w1+ (*iriData)[a][b + 1][1] *w3+ (*iriData)[a + 1][b][1] * w2+ (*iriData)[a + 1][b + 1][1] *w4;
                        nHe = (*iriData)[a][b][2] * w1+ (*iriData)[a][b + 1][2] *w3+ (*iriData)[a + 1][b][2] * w2+ (*iriData)[a + 1][b + 1][2] *w4;
                        nO = (*iriData)[a][b][3] * w1+ (*iriData)[a][b + 1][3] *w3+ (*iriData)[a + 1][b][3] * w2+ (*iriData)[a + 1][b + 1][3] *w4;
                        
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
                        phi = (*tiegcmData)[a][b][0] * w1+ (*tiegcmData)[a][b + 1][0] *w3+ (*tiegcmData)[a + 1][b][0] * w2+ (*tiegcmData)[a + 1][b + 1][0] *w4;
                        //vi_lon = (*tiegcmData)[a][b][1] * w1+ (*tiegcmData)[a][b + 1][1] *w3+ (*tiegcmData)[a + 1][b][1] * w2+ (*tiegcmData)[a + 1][b + 1][1]*w4;
                        vi_lat = (*tiegcmData)[a][b][2] * w1+ (*tiegcmData)[a][b + 1][2] *w3+ (*tiegcmData)[a + 1][b][2] * w2+ (*tiegcmData)[a + 1][b + 1][2] *w4;
                        vi_vert = (*tiegcmData)[a][b][3]* w1 + (*tiegcmData)[a][b + 1][3]*w3 + (*tiegcmData)[a + 1][b][3] * w2+ (*tiegcmData)[a + 1][b + 1][3] *w4;
                        // trans vi in (long, lat, vertical) into (vxGEO, vyGEO, vzGEO)
                        BCARSPRR_08(thet, phis, vi_vert, vi_lat, vi_lat, vxGEO, vyGEO, vzGEO);
                        // trans velocity (vxGEO, vyGEO, vzGEO) into (vxMAG, vyMAG, vzMAG)
                        GEOMAG_08(vxGEO, vyGEO, vzGEO, vxMAG, vyMAG, vzMAG, 1);
                        // value the related values
                        ptrArrayGrids[f][i][j][k]->SetTemperatureIons( ti);  // Ti
                        ptrArrayGrids[f][i][j][k]->Density_H(nH);   // nH
                        ptrArrayGrids[f][i][j][k]->Density_H(nHe);   // nHe
                        ptrArrayGrids[f][i][j][k]->Density_O(nO);   // nO
                        ptrArrayGrids[f][i][j][k]->SetPotential(phi);   // PHI
                        ptrArrayGrids[f][i][j][k]->SetVel3(Vector3(vxMAG, vyMAG, vzMAG));   // vi
                        //
                    } else
                    {
                        // TIEGCM: phi, only!,  - 0.5 for a is because the data type of TIEGCMaDouble = thet * 180.0 / pi / 5.0 + 0.5;
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
                        phi = (*tiegcmData)[a][b][0] * w1+ (*tiegcmData)[a][b + 1][0] *w3+ (*tiegcmData)[a + 1][b][0] * w2+ (*tiegcmData)[a + 1][b + 1][0] *w4;
                        ptrArrayGrids[f][i][j][k]->SetPotential(phi);   // Potential
                        //
                    }
                    //
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
// Test bit coding
inline void TestBitCodine()
{
    uint_64 a = 0;
    a = ~0;
    std::cout << " checka " << a << "\n";
    std::cout << std::bitset<64>(a) << " " << std::bitset<64>(~a) << std::endl;
    //uint_64 b = 11529215046068469760;
    //std::cout << " checkb " << b << "\n";
    //std::cout << std::bitset<64>(b) << " " << std::bitset<64>(~b) << std::endl;
}
// TestSplitOrCombine
inline void TestSplitOrCombie(GridsPoints *****ptrArray,
                              GridsCells ****ptrArrayCells)
{
    int fc = 0;
    int ic = 0;
    int jc = 0;
    int kc = 0;
    int parNum = 187;
    double random = 1.0;
    GridsCells ptrCell = ptrArrayCells[fc][ic][jc][kc];
    // 100 H+, He+, O+
    uint_64 intPos;
    Vector3 tempP_Hos, vel;
    double mu_simu, Ni_simu;
    //Vector3 b3Cell = ptrArrayCells[fc][ic][jc][kc].B3Cell().NormalizedVector();
    for (int i = 0; i < parNum; i++)
    {
        intPos = UniDisInCell(fc, ic, jc, kc);
        tempP_Hos = Uint64ToVector3(intPos);
        MaxwellDisV(ptrArray, intPos, mi0_H, 0.0, vel, mu_simu);
        Ni_simu = 1.0e+25 + 1.0e+25 * (dRand() - 0.5) * random;
        Particles testP = Particles(intPos, tempP_Hos, vel, Ni_simu, mu_simu);
        ptrCell.Particles_H()->push_back(testP);
    }
    //
    std::cout << ptrArrayCells[fc][ic][jc][kc].Particles_H()->size() << "\n";

    ptrArrayCells[fc][ic][jc][kc].SplitOrCombine(0);
    std::cout << ptrArrayCells[fc][ic][jc][kc].Particles_H()->size() << "\n";
}
// test check 
inline void CheckOutParticles(GridsCells ****ptrArrayCells)
{
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    auto ptrParticles = ptrArrayCells[f][i][j][k].Particles_H();
                    for (auto iter = ptrParticles->begin(); iter < ptrParticles->end(); ++iter)
                    {
                        if(!iter->AliveParticle())
                        std::cout << " PH " << f << " " << i << " " << j << " " << k << " " <<ptrParticles->size() << "\n";
                    }
                }
            }
        }
    }
    //
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    auto ptrParticles = ptrArrayCells[f][i][j][k].Particles_He();
                    for (auto iter = ptrParticles->begin(); iter < ptrParticles->end(); ++iter)
                    {
                        if(!iter->AliveParticle())
                        std::cout << " PHe "<< f << " " << i << " " << j << " " << k << " " <<ptrParticles->size() << "\n";
                    }
                }
            }
        }
    }
    //
    for (int f = 0; f < totalFace; f++)
    {
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            for (int j = 0; j < fieldsGridsSize; j++)
            {
                for (int k = 0; k < fieldsGridsSize; k++)
                {
                    auto ptrParticles = ptrArrayCells[f][i][j][k].Particles_O();
                    for (auto iter = ptrParticles->begin(); iter < ptrParticles->end(); ++iter)
                    {
                        if(!iter->AliveParticle())
                        std::cout << " PO "<< f << " " << i << " " << j << " " << k << " " <<ptrParticles->size() << "\n";
                    }
                }
            }
        }
    }
}
// test random
inline void TestSpeciesRandom(GridsPoints *****ptrArray,
                              GridsCells ****ptrArrayCells)
{
    int fc = 0;
    int ic = 0;
    int jc = 0;
    int kc = 0;
    int parNum = 18;
    double random = 1.0;
    GridsCells ptrCell = ptrArrayCells[fc][ic][jc][kc];
    // 100 H+
    uint_64 intPos;
    Vector3 tempP_Hos, vel;
    double mu_simu, Ni_simu;
    //Vector3 b3Cell = ptrArrayCells[fc][ic][jc][kc].B3Cell().NormalizedVector();
    for (int i = 0; i < parNum; i++)
    {
        intPos = UniDisInCell(fc, ic, jc, kc);
        tempP_Hos = Uint64ToVector3(intPos);
        MaxwellDisV(ptrArray, intPos, mi0_H, 0.0, vel, mu_simu);
        Ni_simu = 1.0e+25 + 1.0e+25 * (dRand() - 0.5) * random;
        Particles testP = Particles(intPos, tempP_Hos, vel, Ni_simu, mu_simu);
        ptrCell.Particles_H()->push_back(testP);
    }
    int sizeP = ptrCell.Particles_H()->size();
    for( int i = 0; i < sizeP; i++)
    {
        std::cout << (*(ptrCell.Particles_H()))[i].WeightNi() << "\n";
    }
    std::cout << "size "<<ptrArrayCells[fc][ic][jc][kc].Particles_H()->size() << "\n";
    //
    (*(ptrCell.Particles_H()))[10].SetOutParticles();   // setout one of them
    //
    ptrArrayCells[fc][ic][jc][kc].SpeciesRandom(0);
    //
    sizeP = ptrCell.Particles_H()->size();
    for( int i = 0; i < sizeP; i++)
    {
        std::cout << (*(ptrCell.Particles_H()))[i].WeightNi() << "\n";
    }
    std::cout << "size "<<ptrArrayCells[fc][ic][jc][kc].Particles_H()->size() << "\n";
    //
}
// TestDataAccess
inline void TestDataAccess(GridsPoints *****ptrArrayGrids,
                           GridsCells ****ptrArrayCells)
{
    std::cout << " ********** Accessing test Data ************* \n";
    vector<vector<double>> *ptr_b;
    ptr_b = ptrArrayCells[0][0][0][0].VelDist_H();
    
    std::cout << " ---- \n";
    std::cout << ptr_b << " \n";
    std::cout << (*ptr_b)[0][0] << " \n";
    std::cout << (*ptr_b)[1][1] << " \n";
    auto iter = ptr_b->begin();
    //    auto *iter_ptr = *ptr_b->begin();
    std::cout << *(((iter + 1))->begin() + 1) << "\n";
    //    std::cout << *iter_ptr->begin() << "\n";
}
// E calculation test
// test E when only a few particles in domain
inline void TestTemperatureCalculation(GridsPoints *****ptrArray,
                                       GridsCells ****ptrArrayCells,
                                       double ***ptrVolumeCellArray,
                                       double ***ptrVolumeGridArray,
                                       Vector3 *****ptrGradVectorCellArray)
{
    Vector3 ******ptrEVectorFaceArray_dual = NULL;
    //
    PrintOutHdf5_const(ptrArray);
    // put one particles in a cell
    int fc = 0;
    int ic = 8;
    int jc = 8;
    int kc = 8;
    uint_64 intPos1 = UniDisInCell(fc, ic+1, jc+1, kc);
    // test
    uint_64 face_test = 0, ig_test = 8, jg_test = 8, kg_test = 8;
            //iw_test = 0, jw_test = 0, kw_test = 0;
    ig_test = ig_test<<(particlesGridsLevel-fieldsGridsLevel);
    jg_test = jg_test<<(particlesGridsLevel-fieldsGridsLevel);
    kg_test = kg_test<<(particlesGridsLevel-fieldsGridsLevel);
    intPos1 = 0;
    intPos1 = face_test<<61;
    for (int i = 0; i < particlesGridsLevel; i++)
    {
        intPos1 += (((ig_test >> (particlesGridsLevel - 1 - i)) & 1) << (60 - i * 3)) 
                    + (((jg_test >> (particlesGridsLevel - 1 - i)) & 1) << (60 - 1 - i * 3)) 
                    + (((kg_test >> (particlesGridsLevel - 1 - i)) & 1) << (60 - 2 - i * 3));
    }
    //
    Vector3 tempVector1 = Uint64ToVector3(intPos1);
    Vector3 vp_simu1= {0.0,0.0,0.0};
    double mu_simu1;
    MaxwellDisV(ptrArray, intPos1, mi0_H, 0.0, vp_simu1, mu_simu1);
    double Ni_simu1 = 1.0e+26;
    Particles testP1 = Particles(intPos1, tempVector1, vp_simu1, Ni_simu1, mu_simu1);
    struct structg structG = testP1.InttoStrp1();
    std::cout<< " structG_initial " <<structG.face << " " << structG.ig << " " << structG.jg << " " << structG.kg << "\n";
        
    ptrArrayCells[fc][ic][jc][kc].Particles_H()->push_back(testP1);
    //
        IterateParticlesInCells(ptrArray,
                                ptrArrayCells,
                                0);
    //
        CalculatingAveragedPhoVatGrids( ptrArray,
                                        ptrVolumeGridArray,
                                        updateInfoPeriod);
     // dE- parallel, from grad Pe
     // grad at center of cell
     ptrGradVectorCellArray = ValueGradient(ptrArray,
                                            ptrGradVectorCellArray,
                                            ptrVolumeCellArray,
                                            'P');
    // grad at vertex of cells
    UpdateGrad(ptrGradVectorCellArray,
               ptrArray,
               'P');
    // Calculate ve and E at vertex of main cell
    UpdateVeGridsMain(ptrArray,
                      ptrEVectorFaceArray_dual,
                      0);
    //
    Vector3 E_test = ptrArray[fc][ic+1][jc+1][kc]->E3();   
    Vector3 gradPe = ptrArray[fc][ic+1][jc+1][kc]->GradPe();
    double density = ptrArray[fc][ic+1][jc+1][kc]->Density_H();    
    double density_cumu = ptrArray[fc][ic+1][jc+1][kc]->Density_H_cumu();
    std::cout << "E_test " << E_test.norm() << "\n"
                << "gradPe " << gradPe.norm() << "\n"
                << "density " << density << " " << density_cumu << "\n";
    for( auto iter = ptrArrayCells[fc][ic][jc][kc].Particles_H()->begin(); iter != ptrArrayCells[fc][ic][jc][kc].Particles_H()->end(); ++iter)
    {
        structG = iter->InttoStrp1();
        std::cout<< " structG " <<structG.face << " " << structG.ig << " " << structG.jg << " " << structG.kg << "\n";
        //
    }    

    PrintOutHdf5(ptrArray, 1);

    std::cout << " E calculation Test end\n";
    std::cin.get();
}
// particles moving test
inline void TestParticles(GridsPoints *****ptrArray_in)
{
    std::cout << " ********** performing test particle ************ \n "
          <<" set tstep = 0.2s, V0 = (0.0,0.0,0.0)\n";
    //
    //PrintOutHdf5(ptrArray_in, 0, 0);
    //PrintOutHdf5(ptrArray_in, 1, 1);
    int fc1 = 0;
    int ic1 = 8;
    int jc1 = 15;
    int kc1 = 0;

    uint_64 intPos1 = UniDisInCell(fc1, ic1+1, jc1+1, kc1);
    Vector3 tempVector1 = Uint64ToVector3(intPos1);
    Vector3 vp_simu1= {0.0,0.0,0.0};
    double mu_simu1;
    MaxwellDisV(ptrArray_in, intPos1, mi0_H, 0.0, vp_simu1, mu_simu1);
    double Ni_simu1 = 1.0e+20;
    vp_simu1= {0.0,0.0,0.0};
    mu_simu1 *= 1.0;
    Particles testP1 = Particles(intPos1, tempVector1, vp_simu1, Ni_simu1, mu_simu1);
    //
    int fc2 = 2;
    int ic2 = 8;
    int jc2 = 4;
    int kc2 = 0;
    //
    uint_64 intPos2 = UniDisInCell(fc2, ic2+1, jc2+1, kc2);
    Vector3 tempVector2 = Uint64ToVector3(intPos2);
    Vector3 vp_simu2= {0.0,0.0,0.0};
    double mu_simu2;
    MaxwellDisV(ptrArray_in, intPos2, mi0_H, 0.0, vp_simu2, mu_simu2);
    double Ni_simu2 = 1.0e+20;
    vp_simu2= {0.0,0.0,0.0};
    // show the starting point
    //testP2 = testP1;
    Particles testP2 = Particles(intPos2, tempVector2, vp_simu2, Ni_simu2, mu_simu2);
    //
    
    int fc3 = 2;
    int ic3 = 8;
    int jc3 = 7;
    int kc3 = 0;
    //
    uint_64 intPos3 = UniDisInCell(fc3, ic3+1, jc3+1, kc3);
    Vector3 tempVector3 = Uint64ToVector3(intPos3);
    Vector3 vp_simu3= {0.0,0.0,0.0};
    double mu_simu3;
    MaxwellDisV(ptrArray_in, intPos3, mi0_H, 0.0, vp_simu3, mu_simu3);
    double Ni_simu3 = 1.0e+20;
    vp_simu3= {0.0,0.0,0.0};
    // show the starting point
    //testP2 = testP1;
    Particles testP3 = Particles(intPos3, tempVector3, vp_simu3, Ni_simu3, mu_simu3);
    //
    std::cout <<"mu "<< mu_simu1 << " " <<  mu_simu2 << " " << mu_simu3 << "\n";
    //
    std::ofstream outP("./MyData/PosOut.txt");

    int check1 = 0;
    int check2 = 0;
    int check3 = 0;

    // test! particles moving
    for (int timeline = 0; timeline < 800000; timeline++)
    {
        struct structg strg_in1 = testP1.InttoStrp1();
        struct structPar str_Par1;
        StructPar( strg_in1, str_Par1);
        if (check1 == 0)
            check1 = testP1.BorisMethod(str_Par1, ptrArray_in, mi0_H, 0);
        double velPar = testP1.VelParticles().norm();
        //
        struct structg strg_in2 = testP2.InttoStrp1();
        struct structPar str_Par2;
        StructPar( strg_in2, str_Par2);
        if (check2 == 0)
            check2 = testP2.BorisMethod(str_Par2, ptrArray_in, mi0_H, 0);
            //
        struct structg strg_in3 = testP3.InttoStrp1();
        struct structPar str_Par3;
        StructPar( strg_in3, str_Par3);
        if (check3 == 0)
            check3 = testP3.BorisMethod(str_Par3, ptrArray_in, mi0_H, 0);
        // remove drift velocity
    //    Vector3 tempb;
    //    Vector3 b3_pa;
    //    Vector3 tempb1 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 1][strg_in1.kg]->B3();
    //    Vector3 tempb2 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 1][strg_in1.kg]->B3();
    //    Vector3 tempb3 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 2][strg_in1.kg]->B3();
    //    Vector3 tempb4 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 2][strg_in1.kg]->B3();
    //    Vector3 tempb5 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 1][strg_in1.kg + 1]->B3();
    //    Vector3 tempb6 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 1][strg_in1.kg + 1]->B3();
    //    Vector3 tempb7 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 2][strg_in1.kg + 1]->B3();
    //    Vector3 tempb8 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 2][strg_in1.kg + 1]->B3();
    //    double w1 = str_Par.w000; //= 1.0 - 1.0 * (strg_in1.iw + 0.5) * (strg_in1.jw + 0.5) * (strg_in1.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w2 = str_Par.w100; //= 1.0 - 1.0 * (cellSize1 - strg_in1.iw - 0.5) * (strg_in1.jw + 0.5) * (strg_in1.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w3 = str_Par.w110; //= 1.0 - 1.0 * (cellSize1 - strg_in1.iw - 0.5) * (cellSize1 - strg_in1.jw - 0.5) * (strg_in1.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w4 = str_Par.w010; //= 1.0 - 1.0 * (strg_in1.iw + 0.5) * (cellSize1 - strg_in1.jw - 0.5) * (strg_in1.kw + 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w5 = str_Par.w001; //= 1.0 - 1.0 * (strg_in1.iw + 0.5) * (strg_in1.jw + 0.5) * (cellSize1 - strg_in1.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w6 = str_Par.w101; //= 1.0 - 1.0 * (cellSize1 - strg_in1.iw - 0.5) * (strg_in1.jw + 0.5) * (cellSize1 - strg_in1.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w7 = str_Par.w111; //= 1.0 - 1.0 * (cellSize1 - strg_in1.iw - 0.5) * (cellSize1 - strg_in1.jw - 0.5) * (cellSize1 - strg_in1.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    double w8 = str_Par.w011; //= 1.0 - 1.0 * (strg_in1.iw + 0.5) * (cellSize1 - strg_in1.jw - 0.5) * (cellSize1 - strg_in1.kw - 0.5) / cellSize1 / cellSize1 / cellSize1;
    //    tempb.Setx(tempb1.x() * w1 + tempb2.x() * w2 + tempb3.x() * w3 + tempb4.x() * w4 + tempb5.x() * w5 + tempb6.x() * w6 + tempb7.x() * w7 + tempb8.x() * w8);
    //    tempb.Sety(tempb1.y() * w1 + tempb2.y() * w2 + tempb3.y() * w3 + tempb4.y() * w4 + tempb5.y() * w5 + tempb6.y() * w6 + tempb7.y() * w7 + tempb8.y() * w8);
    //    tempb.Setz(tempb1.z() * w1 + tempb2.z() * w2 + tempb3.z() * w3 + tempb4.z() * w4 + tempb5.z() * w5 + tempb6.z() * w6 + tempb7.z() * w7 + tempb8.z() * w8);
    //    b3_pa = Vector3(tempb).NormalizedVector();
    //    Vector3 tempV = b3_pa.ScaleProduct(testP1.VelParticles().DotProduct(b3_pa));
    //    testP1.SetVelocity(tempV);
        //

        //    if( check2 ==0)
        //    check2 = testP2.BorisMethod( &tempStr2, ptrArray, mi0_O, 0);
        //    std::cout << testP1.VelParticles().x() << " " << testP1.VelParticles().y()
        //              << " " << testP1.VelParticles().z() << " =>> ";
        //    std::cout << testP1.PosParticles().x() << " " << testP1.PosParticles().y()
        //              << " " << testP1.PosParticles().z() << std::endl;
        //    if( check1 != 0 )
        //    std::cout << " " << timeline << " ";
        //   if( check2 != 0 )
        //   {
        //   std::cout << " " <<  timeline << " ";
        //   std::cin.get();
        //   }

        if (timeline % 2 == 0)
        {

            //    std::cout<< testP2.VelParticles().norm() << std::endl;
            outP << testP1.PosParticles().x() << " " << testP1.PosParticles().y() << " " << testP1.PosParticles().z() << " "
                 << testP2.PosParticles().x() << " " << testP2.PosParticles().y() << " " << testP2.PosParticles().z() << " " 
                 << testP3.PosParticles().x() << " " << testP3.PosParticles().y() << " " << testP3.PosParticles().z() << " " 
                 << velPar <<"\n";
        }
    }
    outP.close();
    //
    std::cout << " test Particles moving end\n";
    std::cin.get();
}
//*******************************************************************************
// collision test
inline void TestColisionTwoParticles(GridsPoints *****ptrArray,
                                     GridsCells ****ptrArrayCells)
{
    int fc = 0;
    int ic = 10;
    int jc = 10;
    int kc = 10;
    double random = 1.0;
    //GridsCells ptrCell = ptrArrayCells[fc][ic][jc][kc];
    // H+, He+, O+
    //uint_64 intPos;
    Vector3 tempP_Hos, vel;
    //double mu_simu, Ni_simu;
    //
    double mi_a = mi0_H;
    double mi_b = mi0_He;
    //
    //Vector3 b3Cell = ptrArrayCells[fc][ic][jc][kc].B3Cell().NormalizedVector();
    // H
    uint_64 intPos1 = UniDisInCell(fc, ic, jc, kc);
    Vector3 tempVector1 = Uint64ToVector3(intPos1);
    Vector3 vp_simu1;
    double mu_simu1;
    MaxwellDisV(ptrArray, intPos1, mi_a, 0.0, vp_simu1, mu_simu1);
    double Ni_simu1 = 1.0e+20 + 1.0e+20 * (dRand() - 0.5) * 2 * random;
    Particles testP1 = Particles(intPos1, tempVector1, vp_simu1, Ni_simu1, mu_simu1);
    // He
    uint_64 intPos2 = UniDisInCell(fc, ic, jc, kc);
    Vector3 tempVector2 = Uint64ToVector3(intPos2);
    Vector3 vp_simu2;
    double mu_simu2;
    MaxwellDisV(ptrArray, intPos2, mi_b, 0.0, vp_simu2, mu_simu2);
    double Ni_simu2 = 1.0e+20 + 1.0e+20 * (dRand() - 0.5) * 2 * random;
    Particles testP2 = Particles(intPos2, tempVector2, vp_simu2, Ni_simu2, mu_simu2);
    //
    double lnAab = 10.0;
    //
    int len = 200;
    double test_mome[len] = {0.0};

    for (int ii = 0; ii < len; ii++)
    {
        //coll_vel_chng(  0, 0,
        //                0.000001, 1.0e+8,
        //                lnAab,
        //                ptrArray,
        //                testP1, testP2);
        coll_vel_chng(0, 1,
                      1, 1.0e+8,
                      lnAab,
                      ptrArray,
                      testP1, testP2);

        Vector3 ba = testP1.B3atParticles(ptrArray);
        Vector3 bb = testP2.B3atParticles(ptrArray);
        std::cout << testP1.VelCollParticles(ba, mi_a).norm2() * mi_a * testP1.WeightNi() + testP2.VelCollParticles(bb, mi_b).norm2() * mi_b * testP2.WeightNi() << " \n";
        test_mome[ii] = testP1.VelCollParticles(ba, mi_a).norm2() * mi_a * testP1.WeightNi() + testP2.VelCollParticles(bb, mi_b).norm2() * mi_b * testP2.WeightNi();
    }
    //
    std::ofstream outP("./MyData/test_mome.txt");
    for (int ii = 0; ii < len; ii++)
    {
        outP << test_mome[ii] << "\n";
    }
    outP.close();
}
// Test maxwell distribution
inline void TestMaxwellDistri()
{
    double mkb[3] = { 1.2123442e-4, 4.813845268e-4,1.9242837e-3};
    // 1-D
    std::cout << " 1-D \n";
    int len = 10000;
    vector<double> ptrVel(len);
    double culVel = 0.0;
    double culT = 0.0;
    for( int i = 0; i < len; i++)
    {
        ptrVel[i] =  MaxwellDisV(ikT, 0.0, mi0_H);
        std::cout << ptrVel[i] << "\n";
    }
    // calculate average
    for( int i = 0 ; i < len; i++)
    {
        culVel+=ptrVel[i];
    }
    culVel = culVel / len;
    std::cout << " averageVel= " << culVel << "\n";
    // T
    for( int i = 0; i< len; i++)
    {
        culT += ptrVel[i] * ptrVel[i];
    }
    culT = (culT / len - culVel * culVel) * mkb[0];
    std::cout << " T= " << culT << "\n";
    //
    // 3-D
    std::cout << " 3-D \n";
    vector<Vector3> ptrVel3(len);
    Vector3 culVel3 = { 0.0, 0.0, 0.0};
    double culT3 = 0.0;
    for( int i = 0; i < len; i++)
    {
        ptrVel3[i] =  Vector3(MaxwellDisV(ikT, 0.0, mi0_H),MaxwellDisV(ikT, 0.0, mi0_H),MaxwellDisV(ikT, 0.0, mi0_H));
        std::cout << ptrVel3[i].x() << " " << ptrVel3[i].y() << " " << ptrVel3[i].z() << "\n";
    }
    // calculate average
    for( int i = 0 ; i < len; i++)
    {
        culVel3= ptrVel3[i].PlusProduct(culVel3);
    }
    culVel3 = culVel3.ScaleProduct(1.0/len);
    // T3
    for( int i = 0; i< len; i++)
    {
        culT3 += ptrVel3[i].norm2();
    }
    culT3 = (culT3 / len - culVel * culVel) * mkb[0] / 3.0;
    std::cout << " averageVel3= " << culVel3.norm() << "\n";
    std::cout << " T3= " << culT3 << "\n";


    //
    std::cout << " MaxwellTest end \n";
    std::cin.get();
}
// test collision in Cells
inline void TestCollisionInCell(GridsPoints *****ptrArray,
                                GridsCells ****ptrArrayCells)
{
    int fc = 0;
    int ic = 8;
    int jc = 8;
    int kc = 12;
    int parNum = 1000;
    double random = 1.0;
        double mkb[3] = { 1.2123442e-4, 4.813845268e-4,1.9242837e-3};
    GridsCells ptrCell = ptrArrayCells[fc][ic][jc][kc];
    // 100 H+, He+, O+
    uint_64 intPos;
    Vector3 tempP_Hos, vel;
    double mu_simu, Ni_simu;
    //Vector3 b3Cell = ptrArrayCells[fc][ic][jc][kc].B3Cell().NormalizedVector();
    
    //double Te, ne, TH, THe, TO;
    Vector3 vH_para, vHe_para, vO_para;
    double T_H_initial= 3000.0, T_He_initial = 1000.0, T_O_initial = 0.0;   // + 1000K
    double lnAab;
    int f = fc;
    int i = ic;
    int j = jc;
    int k = kc;

    for (i = 0; i < parNum; i++)
    {
        intPos = UniDisInCell(fc, ic, jc, kc);
        tempP_Hos = Uint64ToVector3(intPos);
        MaxwellDisV(ptrArray, intPos, mi0_H, T_H_initial, vel, mu_simu);
        Ni_simu = 1.0e+25 + 1.0e+25 * (dRand() - 0.5) * random;
        // test for pp & pr temperature
        //vel = Vector3(0.0,0.0,0.0);
        Particles testP = Particles(intPos, tempP_Hos, vel, Ni_simu, mu_simu);
        // test for vel dist
        Vector3 b3 = testP.B3atParticles(ptrArray);
        if( i < parNum/2)
            testP.SetVelocity( b3.NormalizedVector().ScaleProduct( 25000.0 + testP.VelParticles().DotProduct(b3.NormalizedVector())* (1.0)));
        else
            testP.SetVelocity( b3.NormalizedVector().ScaleProduct( -25000.0 + testP.VelParticles().DotProduct(b3.NormalizedVector())* (1.0)));
        //
        ptrCell.Particles_H()->push_back(testP);
    }
    for (i = 0; i < parNum; i++)
    {
        intPos = UniDisInCell(fc, ic, jc, kc);
        tempP_Hos = Uint64ToVector3(intPos);
        MaxwellDisV(ptrArray, intPos, mi0_He, T_He_initial, vel, mu_simu);
        Ni_simu = 1.0e+25 + 1.0e+25 * (dRand() - 0.5) * random;
        Particles testP = Particles(intPos, tempP_Hos, vel, Ni_simu, mu_simu);
        ptrCell.Particles_He()->push_back(testP);
    }
    for (i = 0; i < parNum; i++)
    {
        intPos = UniDisInCell(fc, ic, jc, kc);
        tempP_Hos = Uint64ToVector3(intPos);
        MaxwellDisV(ptrArray, intPos, mi0_O, T_O_initial, vel, mu_simu);
        Ni_simu = 1.0e+25 + 1.0e+25 * (dRand() - 0.5) * random;
        Particles testP = Particles(intPos, tempP_Hos, vel, Ni_simu, mu_simu);
        ptrCell.Particles_O()->push_back(testP);
    }
    // calculate temperature Test for H
    Vector3 b3 = {0.0, 0.0, 0.0};
    Vector3 culVel_H = {0.0, 0.0, 0.0}, culVel_He = {0.0, 0.0, 0.0}, culVel_O = {0.0, 0.0, 0.0};
    double temperature_H = 0.0, temperature_He = 0.0, temperature_O = 0.0;
    double totalWeight_H = 0.0, totalWeight_He = 0.0, totalWeight_O = 0.0;
    // H
    auto pp_H = ptrArrayCells[f][i][j][k].Particles_H();
    // average vel
    for (auto iter = pp_H->begin(); iter != pp_H->end(); ++iter)
    {
        b3 = iter->B3atParticles(ptrArray);
        culVel_H = culVel_H.PlusProduct( iter->VelCollParticles(b3, mi0_H).ScaleProduct( iter->WeightNi()));
        totalWeight_H += iter->WeightNi();
    }
    culVel_H = culVel_H.ScaleProduct( 1.0/totalWeight_H);
    //
    for (auto iter = pp_H->begin(); iter != pp_H->end(); ++iter)
    {
        b3 = iter->B3atParticles(ptrArray);
        temperature_H += iter->VelCollParticles(b3, mi0_H).norm2()  * iter->WeightNi();
        //temperature_H += iter->VelCollParticles(b3, mi0_H).MinusProduct(culVel_H).norm2() * iter->WeightNi();
    }
    temperature_H = (temperature_H / totalWeight_H - culVel_H.norm2() )* mkb[0] / 3.0;
    std::cout << " T_H= " << temperature_H << " culVel_H " << culVel_H.norm() << "\n";
    // test collision
    // speciesRandom
    //ptrArrayCells[f][i][j][k].SpeciesRandom(0);
    //ptrArrayCells[f][i][j][k].SpeciesRandom(1);
    //ptrArrayCells[f][i][j][k].SpeciesRandom(2);
    //// LnA, need species (a,b), Te, ne, Ta, Tb, va_parallel, vb_parallel
    ////
    //Te = ptrArray[f][i+1][j+1][k]->Temperature() + ptrArray[f][i+2][j+1][k]->Temperature()
    //    +ptrArray[f][i+1][j+2][k]->Temperature() + ptrArray[f][i+2][j+2][k]->Temperature()
    //    +ptrArray[f][i+1][j+1][k+1]->Temperature() + ptrArray[f][i+2][j+1][k+1]->Temperature()
    //    +ptrArray[f][i+1][j+2][k+1]->Temperature() + ptrArray[f][i+2][j+2][k+1]->Temperature();
    //Te = Te / 8.0;
    //ne = ptrArray[f][i+1][j+1][k]->Density() + ptrArray[f][i+2][j+1][k]->Density()
    //    +ptrArray[f][i+1][j+2][k]->Density() + ptrArray[f][i+2][j+2][k]->Density()
    //    +ptrArray[f][i+1][j+1][k+1]->Density() + ptrArray[f][i+2][j+1][k+1]->Density()
    //    +ptrArray[f][i+1][j+2][k+1]->Density() + ptrArray[f][i+2][j+2][k+1]->Density();
    //ne = ne / 8.0;
    ////
    //ne = 1.0e+10;
    ////  velocity in cell and Temperature
    //vH_para=Vector3( ptrArray[f][i+1][j+1][k]->VelH3().x() + ptrArray[f][i+2][j+1][k]->VelH3().x()
    //                +ptrArray[f][i+1][j+2][k]->VelH3().x() + ptrArray[f][i+2][j+2][k]->VelH3().x()
    //                +ptrArray[f][i+1][j+1][k+1]->VelH3().x() + ptrArray[f][i+2][j+1][k+1]->VelH3().x()
    //                +ptrArray[f][i+1][j+2][k+1]->VelH3().x() + ptrArray[f][i+2][j+2][k+1]->VelH3().x(),
    //                ptrArray[f][i+1][j+1][k]->VelH3().y() + ptrArray[f][i+2][j+1][k]->VelH3().y()
    //                +ptrArray[f][i+1][j+2][k]->VelH3().y() + ptrArray[f][i+2][j+2][k]->VelH3().y()
    //                +ptrArray[f][i+1][j+1][k+1]->VelH3().y() + ptrArray[f][i+2][j+1][k+1]->VelH3().y()
    //                +ptrArray[f][i+1][j+2][k+1]->VelH3().y() + ptrArray[f][i+2][j+2][k+1]->VelH3().y(),
    //                ptrArray[f][i+1][j+1][k]->VelH3().z() + ptrArray[f][i+2][j+1][k]->VelH3().z()
    //                +ptrArray[f][i+1][j+2][k]->VelH3().z() + ptrArray[f][i+2][j+2][k]->VelH3().z()
    //                +ptrArray[f][i+1][j+1][k+1]->VelH3().z() + ptrArray[f][i+2][j+1][k+1]->VelH3().z()
    //                +ptrArray[f][i+1][j+2][k+1]->VelH3().z() + ptrArray[f][i+2][j+2][k+1]->VelH3().z()
    //                );
    //vH_para = vH_para.ScaleProduct(0.125);
    //TH = ptrArrayCells[f][i][j][k].TempColl(0, vH_para);
    //vHe_para=Vector3( ptrArray[f][i+1][j+1][k]->VelHe3().x() + ptrArray[f][i+2][j+1][k]->VelHe3().x()
    //                 +ptrArray[f][i+1][j+2][k]->VelHe3().x() + ptrArray[f][i+2][j+2][k]->VelHe3().x()
    //                 +ptrArray[f][i+1][j+1][k+1]->VelHe3().x() + ptrArray[f][i+2][j+1][k+1]->VelHe3().x()
    //                 +ptrArray[f][i+1][j+2][k+1]->VelHe3().x() + ptrArray[f][i+2][j+2][k+1]->VelHe3().x(),
    //                 ptrArray[f][i+1][j+1][k]->VelHe3().y() + ptrArray[f][i+2][j+1][k]->VelHe3().y()
    //                 +ptrArray[f][i+1][j+2][k]->VelHe3().y() + ptrArray[f][i+2][j+2][k]->VelHe3().y()
    //                 +ptrArray[f][i+1][j+1][k+1]->VelHe3().y() + ptrArray[f][i+2][j+1][k+1]->VelHe3().y()
    //                 +ptrArray[f][i+1][j+2][k+1]->VelHe3().y() + ptrArray[f][i+2][j+2][k+1]->VelHe3().y(),
    //                 ptrArray[f][i+1][j+1][k]->VelHe3().z() + ptrArray[f][i+2][j+1][k]->VelHe3().z()
    //                 +ptrArray[f][i+1][j+2][k]->VelHe3().z() + ptrArray[f][i+2][j+2][k]->VelHe3().z()
    //                 +ptrArray[f][i+1][j+1][k+1]->VelHe3().z() + ptrArray[f][i+2][j+1][k+1]->VelHe3().z()
    //                 +ptrArray[f][i+1][j+2][k+1]->VelHe3().z() + ptrArray[f][i+2][j+2][k+1]->VelHe3().z()
    //                 );
    //vHe_para = vHe_para.ScaleProduct(0.125);
    //THe = ptrArrayCells[f][i][j][k].TempColl(1, vHe_para);
    //vO_para=Vector3( ptrArray[f][i+1][j+1][k]->VelO3().x() + ptrArray[f][i+2][j+1][k]->VelO3().x()
    //                +ptrArray[f][i+1][j+2][k]->VelO3().x() + ptrArray[f][i+2][j+2][k]->VelO3().x()
    //                +ptrArray[f][i+1][j+1][k+1]->VelO3().x() + ptrArray[f][i+2][j+1][k+1]->VelO3().x()
    //                +ptrArray[f][i+1][j+2][k+1]->VelO3().x() + ptrArray[f][i+2][j+2][k+1]->VelO3().x(),
    //                ptrArray[f][i+1][j+1][k]->VelO3().y() + ptrArray[f][i+2][j+1][k]->VelO3().y()
    //                +ptrArray[f][i+1][j+2][k]->VelO3().y() + ptrArray[f][i+2][j+2][k]->VelO3().y()
    //                +ptrArray[f][i+1][j+1][k+1]->VelO3().y() + ptrArray[f][i+2][j+1][k+1]->VelO3().y()
    //                +ptrArray[f][i+1][j+2][k+1]->VelO3().y() + ptrArray[f][i+2][j+2][k+1]->VelO3().y(),
    //                ptrArray[f][i+1][j+1][k]->VelO3().z() + ptrArray[f][i+2][j+1][k]->VelO3().z()
    //                +ptrArray[f][i+1][j+2][k]->VelO3().z() + ptrArray[f][i+2][j+2][k]->VelO3().z()
    //                +ptrArray[f][i+1][j+1][k+1]->VelO3().z() + ptrArray[f][i+2][j+1][k+1]->VelO3().z()
    //                +ptrArray[f][i+1][j+2][k+1]->VelO3().z() + ptrArray[f][i+2][j+2][k+1]->VelO3().z()
    //                );
    //vO_para = vO_para.ScaleProduct(0.125);
    //TO = ptrArrayCells[f][i][j][k].TempColl(2, vO_para);
    ////  parallel velocity
    //vH_para = b3Cell.ScaleProduct( vH_para.DotProduct( b3Cell));
    //vHe_para= b3Cell.ScaleProduct( vHe_para.DotProduct( b3Cell));
    //vO_para = b3Cell.ScaleProduct( vO_para.DotProduct( b3Cell));
    //
    int len = 1000;
    double test_mome[len];
    double test_energy_H[len], test_energy_He[len], test_energy_O[len];
    
    double test_temp_H[len],test_temp_He[len],test_temp_O[len];
    double test_temp_pp[len], test_temp_pr[len];
    double test_vel_H[len], test_vel_He[len], test_vel_O[len];
    double test_vel_H_pr[len],test_vel_H_pp[len],test_vel_He_pr[len],test_vel_He_pp[len],test_vel_O_pr[len],test_vel_O_pp[len];
    double velDist_H_output[len/20][velDistRange_para][velDistRange_mu];
    lnAab = 10.0;
    //
    double energy = 0.0; 
    //double mu;
    //double cab = (qi0 * qi0 / e_const / mi0_H) * (qi0 * qi0 / e_const / mi0_H)/ 4.0 / 3.1415926535;
    //std::cout << " aa " << cab / 0.5 / 0.5 << " ab " << cab / ( 16.0/25.0) << " ac " << cab / (16.0*16.0/17.0*17.0) << "\n"
    //        << " bb " << cab / ( 4.0) << " bc " << cab/( 16.0*16.0/5.0*5.0) << " cc " << cab / 64.0 << " \n";
    //
    for (int ii = 0; ii < len; ii++)
    {
        energy = 0.0;
        b3 = {0.0, 0.0, 0.0};
        culVel_H = {0.0, 0.0, 0.0}, culVel_He = {0.0, 0.0, 0.0}, culVel_O = {0.0, 0.0, 0.0};
        temperature_H = 0.0, temperature_He = 0.0, temperature_O = 0.0;
        totalWeight_H = 0.0, totalWeight_He = 0.0, totalWeight_O = 0.0;
        // speciesRandom
        ptrArrayCells[f][i][j][k].SpeciesRandom(0);
        ptrArrayCells[f][i][j][k].SpeciesRandom(1);
        ptrArrayCells[f][i][j][k].SpeciesRandom(2);
        // print all particles
        //std::cout << " H \n";
        //for( auto iter = ptrArrayCells[f][i][j][k].Particles_H()->begin(); iter != ptrArrayCells[f][i][j][k].Particles_H()->end(); ++iter)
        //{
        //    std::cout << iter->WeightNi() << " " << iter->VelParticles().norm() << " " << iter->PosParticles().norm() << " \n";
        //}
        //std::cout << " He \n";
        //for( auto iter = ptrArrayCells[f][i][j][k].Particles_He()->begin(); iter != ptrArrayCells[f][i][j][k].Particles_He()->end(); ++iter)
        //{
        //    std::cout << iter->WeightNi() << " " << iter->VelParticles().norm() << " " << iter->PosParticles().norm() << " \n";
        //}
        //std::cout << " timestep " << ii << "end \n";
        //std::cin.get();
        //
        //
        //  H
        // average vel
        for (auto iter = pp_H->begin(); iter != pp_H->end(); ++iter)
        {
            b3 = iter->B3atParticles(ptrArray);
            energy += iter->VelCollParticles(b3, mi0_H).norm2() * mi0_H * iter->WeightNi();
            culVel_H = culVel_H.PlusProduct( iter->VelCollParticles(b3, mi0_H).ScaleProduct( iter->WeightNi()));
            totalWeight_H += iter->WeightNi();
        }
        culVel_H = culVel_H.ScaleProduct( 1.0/totalWeight_H);
        //
        for (auto iter = pp_H->begin(); iter != pp_H->end(); ++iter)
        {
            b3 = iter->B3atParticles(ptrArray);
            temperature_H += iter->VelCollParticles(b3, mi0_H).norm2()  * iter->WeightNi();
            //temperature_H += iter->VelCollParticles(b3, mi0_H).MinusProduct(culVel_H).norm2() * iter->WeightNi();
        }
        temperature_H = (temperature_H / totalWeight_H - culVel_H.norm2() )* mkb[0] / 3.0;
     //   test_energy_H[ii] = energy;
        //temperature_H = temperature_H / totalWeight_H * mkb[0];
        //  He
        energy = 0.0;
        auto pp_He = ptrArrayCells[f][i][j][k].Particles_He();
        // average vel
        for (auto iter = pp_He->begin(); iter != pp_He->end(); ++iter)
        {
            b3 = iter->B3atParticles(ptrArray);
            energy += iter->VelCollParticles(b3, mi0_He).norm2() * mi0_He * iter->WeightNi();
            culVel_He = culVel_He.PlusProduct( iter->VelCollParticles(b3, mi0_He).ScaleProduct( iter->WeightNi()));
            totalWeight_He += iter->WeightNi();
        }
        culVel_He = culVel_He.ScaleProduct( 1.0/totalWeight_He);
        //
        for (auto iter = pp_He->begin(); iter != pp_He->end(); ++iter)
        {
            b3 = iter->B3atParticles(ptrArray);
            temperature_He += iter->VelCollParticles(b3, mi0_He).norm2()  * iter->WeightNi();
            //temperature_He += iter->VelCollParticles(b3, mi0_He).MinusProduct(culVel_He).norm2()* iter->WeightNi();
        }
        temperature_He = (temperature_He / totalWeight_He - culVel_He.norm2() )* mkb[1]/ 3.0;
    //    test_energy_He[ii] = energy;
        //temperature_He = temperature_He / totalWeight_He * mkb[1];
        // O
        energy = 0.0;
        auto pp_O = ptrArrayCells[f][i][j][k].Particles_O();
        // average vel
        for (auto iter = pp_O->begin(); iter != pp_O->end(); ++iter)
        {
            b3 = iter->B3atParticles(ptrArray);
            energy += iter->VelCollParticles(b3, mi0_O).norm2() * mi0_O * iter->WeightNi();
            culVel_O = culVel_O.PlusProduct( iter->VelCollParticles(b3, mi0_O).ScaleProduct( iter->WeightNi()));
            totalWeight_O += iter->WeightNi();
        }
        culVel_O = culVel_O.ScaleProduct( 1.0/totalWeight_O);
        //
        for (auto iter = pp_O->begin(); iter != pp_O->end(); ++iter)
        {
            b3 = iter->B3atParticles(ptrArray);
            temperature_O += iter->VelCollParticles(b3, mi0_O).norm2()  * iter->WeightNi();
            //temperature_O += iter->VelCollParticles(b3, mi0_O).MinusProduct(culVel_O).norm2()* iter->WeightNi();
        }
        
        temperature_O = (temperature_O / totalWeight_O - culVel_O.norm2() )* mkb[2]/ 3.0;
    //    test_energy_O[ii] = energy;
        //temperature_O = temperature_O / totalWeight_O * mkb[2];
        //
        test_mome[ii] = ptrArrayCells[fc][ic][jc][kc].TestCollEnergy(ptrArray, test_energy_H[ii], test_energy_He[ii], test_energy_O[ii]);
        //std::cout << test_mome[ii] << " " << test_energy_H[ii] << " " << test_energy_He[ii] << " " << test_energy_O[ii] << "\n";
        //std::cin.get();
        //
        test_temp_H[ii] = temperature_H;
        test_temp_He[ii] = temperature_He;
        test_temp_O[ii] = temperature_O;
        //*********************************************************************************************
        // Test para and perp temperature && vel distribution velArray[vel_pr][vel_pp] && <v_pp^2>, <v_pr>
        ptrArrayCells[f][i][j][k].velDistArray(0, ptrArray);
        ptrArrayCells[f][i][j][k].velDistArray(1, ptrArray);
        ptrArrayCells[f][i][j][k].velDistArray(2, ptrArray);
        ////  <v_pp^2> and <v_pr>
        //test_vel_H_pr[ii] = sqrt(ptrArrayCells[f][i][j][k].Vel_H_pr());
        //test_vel_H_pp[ii] = sqrt(ptrArrayCells[f][i][j][k].Vel_H_pp());
        //test_vel_He_pr[ii] = sqrt(ptrArrayCells[f][i][j][k].Vel_He_pr());
        //test_vel_He_pp[ii] = sqrt(ptrArrayCells[f][i][j][k].Vel_He_pp());
        //test_vel_O_pr[ii] = sqrt(ptrArrayCells[f][i][j][k].Vel_O_pr());
        //test_vel_O_pp[ii] = sqrt(ptrArrayCells[f][i][j][k].Vel_O_pp());
        ////
        //test_temp_pp[ii] = ptrArrayCells[f][i][j][k].Temp_H_pp();
        //test_temp_pr[ii] = ptrArrayCells[f][i][j][k].Temp_H_pr();
        // Average velocity
        test_vel_H[ii] = culVel_H.norm(); 
        test_vel_He[ii] = culVel_He.norm();
        test_vel_O[ii] = culVel_O.norm();
        // vel distribution H
        auto velDist_H = ptrArrayCells[f][i][j][k].VelDist_H();
        //auto velDist_He= ptrArrayCells[f][i][j][k].VelDist_He();
        //auto velDist_O = ptrArrayCells[f][i][j][k].VelDist_O();
        if( ii % 20 ==0)
        {
            for(i = 0; i < velDistRange_para; i++)
            {
                for(j = 0; j < velDistRange_mu; j++)
                {
                    velDist_H_output[ii/20][i][j] = (*velDist_H)[(uint64_t)i][(uint64_t)j];
                //    std::cout << velDist_H_output[ii%20][i][j] <<"\n";
                }
            }
            //
            //for( int i = 0; i < velDistRange_para; i++)
            //{
            //    for( int j = 0; j < velDistRange_mu; j++)
            //    {
            //    //    std::cout <<  (*velDist_H)[i][j];
            //        if( j == velDistRange_mu-1)
            //        {
            //            std::cout  << "\n";
            //        }
            //        else 
            //        {
            //            std::cout  << " ";
            //        }
            //    }
            //}
            //
        //    std::cout << " ************************\n ";
            //for( int i = 0; i < velDistRange_para; i++)
            //{
            //    for( int j = 0; j < velDistRange_mu; j++)
            //    {
            //    //    std::cout <<  velDist_H_output[ii/20][i][j] ;
            //        if( j == velDistRange_mu-1)
            //        {
            //            std::cout  << "\n";
            //        }
            //        else 
            //        {
            //            std::cout  << " ";
            //        }
            //    }
            //}
        }
        // reset vel dist
        ptrArrayCells[f][i][j][k].ResetVelDistArray(1);
        // collision 
        //*****************************************************************************
        //  H - H
        //lnAab = LnAab(0, 0, Te, ne, TH, TH, vH_para, vH_para);
        ptrArrayCells[f][i][j][k].self_coll(0, lnAab, ptrArray);
        //////////// He - He
        ////////////lnAab = LnAab(1, 1, Te, ne, THe, THe, vHe_para, vHe_para);
    //    ptrArrayCells[f][i][j][k].self_coll(1, lnAab, ptrArray);
    //    ////////////// O - O
    //    //////////lnAab = LnAab(2, 2, Te, ne, TO, TO, vO_para, vO_para);
    //    ptrArrayCells[f][i][j][k].self_coll(2, lnAab, ptrArray);
    //    //// H - He
    //    //lnAab = LnAab(0, 1, Te, ne, TH, THe, vH_para, vHe_para);
    //    ptrArrayCells[f][i][j][k].mutual_coll(0, 1, lnAab, ptrArray);
    //    //// H - O
    //    //////lnAab = LnAab(0, 2, Te, ne, TH, TO, vH_para, vO_para);
    //    ptrArrayCells[f][i][j][k].mutual_coll(0, 2, lnAab, ptrArray);
    //    //////////// He - O
    //    ////////lnAab = LnAab(1, 2, Te, ne, THe, TO, vHe_para, vO_para);
    //    ptrArrayCells[f][i][j][k].mutual_coll(1, 2, lnAab, ptrArray);
        //
    }
    //
    std::ofstream outP("./MyData/test_mome.txt");
    for (int ii = 0; ii < len; ii++)
    {
        outP << test_energy_H[ii] << " " << test_energy_He[ii] << " " << test_energy_O[ii] << " " << test_mome[ii] << "\n";
    }
    outP.close();
    // 
    std::ofstream outTemperature( "./MyData/test_temperature.txt");
    for (int ii = 0; ii < len; ii++)
    {
        outTemperature << test_temp_H[ii] << " " <<test_temp_He[ii]<< " " << test_temp_O[ii] << "\n";
    }
    outTemperature.close();
    //
    std::ofstream outTemperature_pp_pr( "./MyData/test_temperature_pp_pr.txt");
    for (int ii = 0; ii < len; ii++)
    {
        outTemperature_pp_pr <<  test_temp_pp[ii] << " " << test_temp_pr[ii]<< "\n";
    }
    outTemperature_pp_pr.close();
    // average velocity
    std::ofstream outAverageVel( "./MyData/test_averageVel.txt");
    for (int ii = 0; ii < len; ii++)
    {
        outAverageVel <<  test_vel_H[ii] << " " << test_vel_He[ii]<< " " << test_vel_O[ii] <<"\n";
    }
    outAverageVel.close();
    // velocity distri
    std::ofstream outVelDist( "./MyData/test_velDist.txt");
    for (int ii = 0; ii < len/20; ii++)
    {
        for(i = 0; i < velDistRange_para; i++)
        {
            for(j = 0; j < velDistRange_mu; j++)
            {
                outVelDist <<  velDist_H_output[ii][i][j] ;
            //    std::cout << velDist_H_output[ii][i][j] <<"\n";
                if( j == velDistRange_mu-1)
                {
                    outVelDist << "\n";
                }
                else 
                {
                    outVelDist << " ";
                }
            }
        }
    }
    outVelDist.close();
    // <v_pp^2> and <v_pr>
    std::ofstream outVel_pp_pr( "./MyData/test_vel_pp_pr.txt");
    for (int ii = 0; ii < len; ii++)
    {
        outVel_pp_pr <<  test_vel_H_pr[ii] << " " << test_vel_H_pp[ii] << " " << test_vel_He_pr[ii] << " "
                        << test_vel_He_pp[ii] << " " << test_vel_O_pr[ii] << " " << test_vel_O_pp[ii]<<"\n";
    }
    outVel_pp_pr.close();




    std::cout << " TestCollisionInCells end \n";
    std::cin.get();
}

#endif