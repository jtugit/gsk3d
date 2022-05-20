#include <iostream>
#include <list>
#include <vector>
#include <memory>
#include <string>
#include "parameters.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module.h"
#include "module_0.h"
#include "module_1.h"
#include <cmath>
#include "H5Cpp.h"
#include <bitset>
#include <omp.h>
#include <fstream>
#include "module_base.h"
#include "gridscells.h"
#include "geopack.h"

using std::cout;
using std::endl;
using std::vector;

void update_timedate(double);
foot3 ****calcFootprints(GridsPoints *****ptrArray,int IYEAR,int IDAY,int IHOUR,int MIN,int ISEC,
    float pdyn,float dst,float Byimf,float Bzimf);

//************************************************************************
//************************************************************************
// Test FUNCTION // Put at end of module.cpp
// This test fuction only apply for fixed global E and B, particle moving would not
// affect the E and B at gridspoints.
// This function only apply for testing the moving particles.
// Assume we have a list of Particles
//************************************************************************
//************************************************************************
void ProcessFunc(double VGSEX, double VGSEY, double VGSEZ, double pdyn, double dst,
    double Byimf, double Bzimf, string workdir)
{
    std::cout << "\n Function starts..." << std::endl;

    // set random seed, separately for each individual thread
    srand((unsigned)time(0));
    for( int i = 0; i < total_thread_num; ++i)
    {
        seed_thread[i] = (uint32_t)rand();
    }

    // Prerun 1.0 // Create Grids, including B and Pos. And then Velocity (corotation) and N (exponential )
    // Prerun 1.1 // And then E (electron momentum equation).
    //std::cout << " Function starts. Create Grid system" << std::endl;
    int iday = dayno(iyr,mon,date);
    RECALC_08(iyr,iday,ihour,imin,isec,VGSEX,VGSEY,VGSEZ);
    // trans the rotating angle speed
    GEOMAG_08(omega_earth_wx,omega_earth_wy,omega_earth_wz,omega_wx,omega_wy,omega_wz,1);

    // generate
    std::cout << " Create Grid points..." << std::endl;

    GridsPoints *****ptrArray = GridsCreation();
    //    PrintOutHdf5( ptrArray, 0, 0);
    //    PrintOutHdf5( ptrArray, 1, 1);
    std::cout << " Create Grids Cells..." << std::endl;
    // Allocate memory for moments array and partcile array
    GridsCells *mem_GridsCells = new GridsCells[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    //
    GridsCells ****ptrArrayCells = GridsCellsCreation(mem_GridsCells,
                                                      ptrArray);
    //
    std::cout << " Create Dual cell..." << std::endl;
    Vector3 *****ptrArray_dual = GridsCreation_dual(ptrArray);
    //    GridsPoints***** ptrArray_bot = GridsCreation( ptrArray_bot, tempGridsCellLevel);
    //    GridsPoints***** ptrArray_top = GridsCreation( ptrArray_top, tempGridsCellLevel);
    //std::cout << " Calculate electron temperature using Titheridge model..." << std::endl;
    //
    //Titheridge_Te(ptrArray, iday); // initial Temprature of electron
                             //    PP_update(ptrArray);     // update PSD at grid points
    // Prerun 1.2 // Create Cell centered field array for nesseary calculation for one face of six
    // The size is [fsize+2][fsize+2][fsize+2]
    /*    Vector3*** ptrVectorCellArray = NULL;
    Vector3*** ptrVeleVectorCellArray = NULL;
    Vector3*** ptrGradVectorCellArray= NULL;
    VectorCellField( ptrVectorCellArray);   // curl B or grad |B|
    VectorCellField_Vel( ptrVeleVectorCellArray);    // vel e
    VectorCellField_Grad( ptrGradVectorCellArray);  // grad Pe 
*/
    /*    Vector3*** ptrVectorCellArray = VectorCellField();  
    Vector3*** ptrVelVectorCellArray = VectorCellField_Vel();
    Vector3*** ptrGradVectorCellArray= VectorCellField_Grad();
*/
    // Prerun 1.3 // Create grids field array of volum for one face of six
    // The size is [fsize+2][fsize+2][fsize+2]
    std::cout << " Create array of cell volumes" << std::endl;
    double ***ptrVolumeCellArray = VolumeCellsField(ptrArray);
    std::cout << " Create array of staggered cell volumes centered with grids" << std::endl;
    // The size is [fsize+1][fsize+1][fsize+1]
    double ***ptrVolumeGridArray = VolumeGridsField(ptrVolumeCellArray);
    // THe size is [fieldsize+1][fieldsize+1][fieldsize+1]
    //std::cout << " Create array of volume weight " << std::endl;
    //double ***ptrVolumeWeightGridArray = VolumeWeightGridsField(ptrArray);
    // Create const array of B at center of each cells
    // [totalface * fsize+2 * fsize+2 * fsize +2]
    //    Vector3***** ptrBVectorCellArray = BVectorCellArray( ptrArray);
    // Presun 1.4 // Create Cell centered field array for E
    // [totalface * fsize+2 * fsize+2 * fsize]
    // Apply space to store
    //**//std::cout << " Create array of E at cell center " << std::endl;
    //**//static Vector3 *mem_EVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    //**//Vector3 *****ptrEVectorCellArray = PtrVectorCellArray(mem_EVectorCellArray);
    //**//// Create Cell centered field array for curl B
    //**//std::cout << " Create array of dB at cell center " << std::endl;
    //**//static Vector3 *mem_CurlBVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    //**//Vector3 *****ptrCurlBVectorCellArray = PtrVectorCellArray(mem_CurlBVectorCellArray);
    // Create Cell centered field array for grad Pe or grad |B|, total index range is totalface*fsize+2*fsize+2*fsize
    std::cout << " Create array of gradients of Pe and |B| " << std::endl;
    static Vector3 *mem_gradVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    Vector3 *****ptrGradVectorCellArray = PtrVectorCellArray(mem_gradVectorCellArray);
    //**//// Create Cell centered field array for Vele
    //**//std::cout << " Create array of Ve at cell center " << std::endl;
    //**//static Vector3 *mem_VeleVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    //**//Vector3 *****ptrVeleVectorCellArray = PtrVectorCellArray(mem_VeleVectorCellArray);
    // Create Cell center field array for dual cell in which it is the grids points of dual cell
    //
    // Create grids array for openmp weighting calculation  // ?
    //**//std::cout << " Create volume weight array " << std::endl;
    //**//static Vector3 *mem_VelWeightGridsArray = new Vector3[3 * totalFace * 8 * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1)];
    //**//static double *mem_MassWeightGridsArray = new double[3 * totalFace * 8 * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1) * (fieldsGridsSize / 2 + 1)];
    //**//Vector3 ******ptrVelWeightGridsArray = PtrVelWeightGridsArray(mem_VelWeightGridsArray);
    //**//double ******ptrMassWeightGridsArray = PtrMassWeightGridsArray(mem_MassWeightGridsArray);
    //**////
    //**//static Vector3 *mem_DivBScalarVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    //**//Vector3 *****ptrDivBVectorCellArray = PtrVectorCellArray(mem_DivBScalarVectorCellArray);
    //**//static Vector3 *mem_posVectorCellArray = new Vector3[totalFace * fieldsGridsSize * fieldsGridsSize * fieldsGridsSize * grid_domain];
    //**//Vector3 *****ptrposVectorCellArray = PtrVectorCellArray(mem_posVectorCellArray);
    // Presun 1.5 // Create Face centered field array for dB, perturbation of magnetic field
    std::cout << " Create B vector of main cells and E vector of dual cells " << std::endl;
    Vector3 ******ptrBVectorFaceArray_main = NULL;
    ptrBVectorFaceArray_main = BVectorFaceArray(ptrBVectorFaceArray_main);
    Vector3 ******ptrEVectorFaceArray_dual = NULL;
    ptrEVectorFaceArray_dual = EVectorFaceArray(ptrEVectorFaceArray_dual);
    // culmulate array
    std::cout << " Create backups array of B vector of main cells and E vector of dual cells " << std::endl;
    Vector3 ******ptrBVectorFaceArray_main_backup = NULL;
    ptrBVectorFaceArray_main_backup = BVectorFaceArray(ptrBVectorFaceArray_main_backup);
    Vector3 ******ptrEVectorFaceArray_dual_backup = NULL;
    ptrEVectorFaceArray_dual_backup = EVectorFaceArray(ptrEVectorFaceArray_dual_backup);
    // indicator
    int ******ptrBCheckFaceArray_main = NULL;
    ptrBCheckFaceArray_main = BVectorFaceArray(ptrBCheckFaceArray_main);
    int ******ptrECheckFaceArray_dual = NULL;
    ptrECheckFaceArray_dual = EVectorFaceArray(ptrECheckFaceArray_dual); // pause
    // [direction * face * (fsize+1) * (fsize+1) * (fsize+1)]
    //    Vector3***** ptrBVectorFaceArray = BVectorFaceArray( ptrArray);
    // Initialize condition for no current case
    // move_type: 0-particles at rest, 1-particles moving with the earth
    //    if( move_type ==1)
    //    SetInitialCondition( ptrArray, ptrVectorCellArray, ptrVolumeCellArray);
    // Initiallize temp_bot and temp_top
    //    InitializeTempGrids( ptrArray, ptrArray_bot, ptrArray_top, tempGridsCellLevel);
    // Prerun 1.4 // Create particles list, initialize the velocity and position of each particles
    //cout << " Create particles list of main domain" << endl;
    //
    //vector<Particles> ptrParticlesList_H;
    //vector<Particles> ptrParticlesList_He;
    //vector<Particles> ptrParticlesList_O;
    //vector<int> ptrParticlesList_out_H;
    //vector<int> ptrParticlesList_out_He;
    //vector<int> ptrParticlesList_out_O;
    ////ptrParticlesList_H.reserve(500000000);
    ////ptrParticlesList_He.reserve(500000000);
    ////ptrParticlesList_O.reserve(500000000);
    ////ptrParticlesList_out_H.reserve(50000000);
    ////ptrParticlesList_out_He.reserve(50000000);
    ////ptrParticlesList_out_O.reserve(50000000);
    ////
    //vector<Particles> ptrParticlesListTemp_H;
    //vector<Particles> ptrParticlesListTemp_He;
    //vector<Particles> ptrParticlesListTemp_O;
    //vector<Particles> ptrParticlesListTempBackup_H;
    //vector<Particles> ptrParticlesListTempBackup_He;
    //vector<Particles> ptrParticlesListTempBackup_O;
    //
    //*********************************test*********************************
    
    //test part
    //TestTransV3_uint64();
    //TestParticles(ptrArray);
    //TestTemperatureCalculation(ptrArray,
    //                           ptrArrayCells,
    //                           ptrVolumeCellArray,
    //                           ptrVolumeGridArray,
    //                           ptrGradVectorCellArray);
    //TestMaxwellDistri();
    //TestColisionTwoParticles(ptrArray, ptrArrayCells);
    //TestCollisionInCell( ptrArray, ptrArrayCells);
    //TestDataAccess( ptrArray, ptrArrayCells);
    //TestSplitOrCombie(ptrArray, ptrArrayCells);
    //TestSpeciesRandom(ptrArray,ptrArrayCells);
    //PrintOutHdf5_Particles_Grids(0,
    //                             ptrArrayCells,
    //                             ptrArray);
    //TestBitCodine();
    // test data
    //char datafile[] = "./ExternalData/";
    //TestFileName(ptrArray, 0);
    //std::cout << " Endtest pause\n " ;
    //std::cout << " \n end test \n start program \n " << std::endl;
    // test double of particles
    // Particles testPar;
    // testPar.SetPosPar(Vector3(4109354.75, -4107348.75, -3780686.0));
    // struct structg testStr;
    // struct structPar testStrPar;
    // //int check = testPar.UpdateUint_64_test();
    // testStr = testPar.InttoStrp1();
    // StructPar(testStr, testStrPar); 
    // for( int i = 0 ; i < 10; ++i)
    // {
        
    // std::cout << fixed << setprecision(16) << " test1 " << " " 
    //             << testPar.PosParticles().x() << " "
    //             << testPar.PosParticles().y() << " "
    //             << testPar.PosParticles().z() << " " <<
    //        testStr.face << " " << testStr.face << " " << testStr.ig << " " << testStr.ig << " " << 
    //        testStr.jg << " " << testStr.jg << " " << testStr.kg << " " << testStr.kg << " " << "\n"
    //        <<  std::bitset<64>(testPar.PosUint()) << "\n";
    // //
    // //check = testPar.BorisMethod(testStrPar, ptrArray, mi0_H,0);
    // //check = testPar.UpdateUint_64_test();
    // std::cout << fixed << setprecision(16) << " test2 " << " " 
    //             << testPar.PosParticles().x() << " "
    //             << testPar.PosParticles().y() << " "
    //             << testPar.PosParticles().z() << " " <<
    //        testStr.face << " " << testStr.face << " " << testStr.ig << " " << testStr.ig << " " << 
    //        testStr.jg << " " << testStr.jg << " " << testStr.kg << " " << testStr.kg << " " << "\n"
    //        <<  std::bitset<64>(testPar.PosUint()) << "\n";
    // }

    // std::cout << " test end \n";         
    //std::cin.get();
    //*********************************************************************************
    int timestart = 0;
    //
    // mirror force
    ptrGradVectorCellArray = ValueGradient(ptrArray,
                                           ptrGradVectorCellArray,
                                           ptrVolumeCellArray,
                                           'B');
    UpdateGrad(ptrGradVectorCellArray,
               ptrArray,
               'B');
    // input data
    vector<vector<vector<double>>> *iriData = new vector<vector<vector<double>>>(201, vector<vector<double>>(201, vector<double>(4)));
    vector<vector<vector<double>>> *tiegcmData = new vector<vector<vector<double>>>(38, vector<vector<double>>(73, vector<double>(4)));
    InputDataFileFromExternalData(ptrArray,timestart,workdir,iriData,tiegcmData);
    UpdateDatafileFromExternalDataForBottomBoundary(ptrArray,timestart,workdir,iriData,tiegcmData);

    //JTU 3/18/2022. trace field lines to 100 km
    cout << " Calculate footprints of grids" << endl;
    foot3 ****footArray = calcFootprints(ptrArray,iyr,iday,ihour,imin,isec,(float)pdyn,(float)dst,
        (float)Byimf,(float)Bzimf);

    // Initialize particles
    if (continueParticles == 0)
    {
        cout << " Initialize Particles\n";
        ParticlesInCellsInitial(ptrArray,
                                ptrArrayCells,
                                0);
    }
    else if (continueParticles == 1)
    {
        // read particles array
        cout << " Reload Particles\n";
        timestart = ReadSavedData(ptrArray,
                                  ptrArrayCells);
    }
    else
    {
        std::cout << " Incorrect value of continueParticles \n";
        std::cin.get();
    }
    cout << " Loading particles completed.\n";

    /*********************************************************************/
    /*----------------- Time advance ---------------------------------- */
    /*********************************************************************/
    // Run 2.0
    for (int timeline = timestart; timeline <= numberTimeStep; ++timeline)
    {
        std::cout << std::endl
                  << " Time step: " << timeline << std::endl;
        // iterate main particles
        double startTime = omp_get_wtime();
        // const info print out
        if (timeline == 0)
        {
            std::cout << " Output grid system..." << std::endl;
            PrintOutHdf5_const(ptrArray);

            std::cout << " Output moments on grids at time step: " << timeline << std::endl; // grids print
            PrintOutHdf5(ptrArray,timeline);
            std::cout << " Output moments in cells at time step: " << timeline << std::endl; // cells print
            PrintOutHdf5Cells(ptrArrayCells, timeline);
        }
        //******************************************************
        // move particles, reconstruct data structure in each cell
        MoveParticles(ptrArray,
                      ptrArrayCells,
                      0);
        MoveParticles(ptrArray,
                      ptrArrayCells,
                      1);
        MoveParticles(ptrArray,
                      ptrArrayCells,
                      2);
        
        double endTime1 = omp_get_wtime();
        cout << "Time elapsed for moving particles: " << (endTime1 - startTime) << endl;

                ////*****************************************************
                //// test collisions      
                //double Te, ne, TH, THe, TO;
                //Vector3 vH_para, vHe_para, vO_para;
                //double lnAab;
                //ptrArrayCells[4][5][5][5].SplitOrCombine(0);
                ////ptrArrayCells[f][i][j][k].SplitOrCombine(1);
                ////ptrArrayCells[f][i][j][k].SplitOrCombine(2);
                //// assume
                //lnAab = 10.0;
                //ptrArrayCells[4][5][5][1].self_coll(0, lnAab, ptrArray);
                //ptrArrayCells[4][5][5][1].mutual_coll(0, 1, lnAab, ptrArray);
                //std::cout << " pause \n";
                //std::cin.get();
                //continue;
                //// *************************************************
        
        //******************************************************
        //  (1) splitOrCombine && (2) collisions && (3) wave-particles interaction 
        OperationsInCells(ptrArray,
                          ptrArrayCells,
                          timeline);
        //
        double endTime2 = omp_get_wtime();
        cout << "Time elapsed for evaluating Coulomb collisions: " <<  (endTime2 - endTime1) << endl; 
        //******************************************************      
        // assign particles to (1) grids && (2) velocity distribution array
        AssignParticles(ptrArray,
                        ptrArrayCells,
                        timeline,
                        0);
        AssignParticles(ptrArray,
                        ptrArrayCells,
                        timeline,
                        1);
        AssignParticles(ptrArray,
                        ptrArrayCells,
                        timeline,
                        2);
        //
        double endTime3 = omp_get_wtime();
        cout << "Time elapsed for calculating particles moments: "<< (endTime3 - endTime2) << endl;   

        // JTu 3/5/2022: update time date (will add updating of Kp, F107, solar wind velocity)
        double secs = (double)isec + tstep;
        update_timedate(secs);

        iday = dayno(iyr,mon,date);
        RECALC_08(iyr,iday,ihour,imin,isec,VGSEX,VGSEY,VGSEZ);

        // trans the rotating angle speed
        GEOMAG_08(omega_earth_wx,omega_earth_wy,omega_earth_wz,omega_wx,omega_wy,omega_wz,1);

        // update electron temperature
        Titheridge_Te(ptrArray, iday); // Temprature of electron

        //******************************************************
        //average rho, v, update grids info B, E & reset pho, v
        if(timeline % updateInfoPeriod == 0)
        {
            // average pho and v during previous "updateInfoPeriod" time steps accumulation
            CalculatingAveragedPhoVatGrids(ptrArray,
                                           ptrVolumeGridArray,
                                           updateInfoPeriod);
            // average velDist during previous "updateVelDist" time steps accumulation
            AverVelDistInCells(ptrArrayCells,
                               updateVelDist);
            if (update_type == 0) //
            {
                // "ptrEVectorFaceArray_dual" is not used in this section
                // but need a variable named by it
                Vector3 ******ptrEVectorFaceArray_dual = NULL;
                // *****************************************************************************************
                // update E, contains two components: <1> component calculated from Ohm's Law, consist 
                // of eXB term and gradPe term. This component is marked as E_parallel <2> component 
                // calculated from gradPotential whose potential is mapping along field lines from that 
                // in low altitude(like external ionosphere datafile). This component is marked as E_perpendicular
                // <1>
                // E_parallel, related to grad Pe
                // calculate gradPe at center of cell, stored in "ptrGradVectorCellArray"
                ptrGradVectorCellArray = ValueGradient(ptrArray,
                                                       ptrGradVectorCellArray,
                                                       ptrVolumeCellArray,
                                                       'P');
                // average gradPe to grids, stored as "gradPe" at CLASS gridsPoints
                UpdateGrad(ptrGradVectorCellArray,
                           ptrArray,
                           'P');
                // Calculate ve and E at vertex of main cell
                // in this section, ve = vi, E is calculated from Ohm's Law
                UpdateVeGridsMain(ptrArray,
                                  NULL,
                                  0);
                // <2> 
                // E_perpendicular, from grad potential
                // calculate gradPotential at center of cell, stored in "ptrGradVectorcellArray"
                ptrGradVectorCellArray = ValueGradient(ptrArray,
                                                       ptrGradVectorCellArray,
                                                       ptrVolumeCellArray,
                                                       'D');
                // average gradPotential to grids, stored as "gradPe" in CLASS gridsPoints
                // because no variable "gradPotential" in CLASS gridsPoints
                UpdateGrad(ptrGradVectorCellArray,
                           ptrArray,
                           'D');
                // Add the component of E, calculated from grad Potential
                UpdateVeGridsMain(ptrArray,
                                  NULL,
                                  1);
                // Tempratuer H, He, O  // currently not used
                // UpdateTempIons(ptrArray, ptrArrayCells);
            }
            else if (update_type == 1) //
            {
                for (int div = 0; div < div_max; div++)
                {
                    // 1. calculate curl B on face of E dual cell ( current j)
                    // 2. interpolate current j into gridpoint of main cell
                    // 3. calculate ve at gridpoints
                    // 4. calculate gradient of Pe at center of main cell
                    // 5. interpolate grad Pe into gridpoint of main cell
                    // 6. calculate E at gridpoint of main cell
                    // 7. calculate dB on face of main cell
                    // 8. interpolate dB into gridpoint of main cell
                    //
                    // update curldB on dual cell face
                    //*******************************************************
                    // Update the dB by calculating curl E
                    // Update the dE by calculating curl dB
                    BVectorFaceArrayUpdate(ptrArray,
                                           ptrBVectorFaceArray_main,
                                           ptrBVectorFaceArray_main_backup,
                                           ptrBCheckFaceArray_main);
                    EVectorFaceArrayUpdate(ptrArray_dual,
                                           ptrBVectorFaceArray_main,
                                           ptrEVectorFaceArray_dual,
                                           ptrEVectorFaceArray_dual_backup,
                                           ptrECheckFaceArray_dual);
                    // update ve and E at gridpoints of main cell after calculating grad Pe
                    // Calculate the gradient of Pe at center of cells and interpolate at grids of main cell
                    // ( fsize+2 * fsize+2 * fsize)
                    ptrGradVectorCellArray = ValueGradient(ptrArray,
                                                           ptrGradVectorCellArray,
                                                           ptrVolumeCellArray,
                                                           'P');
                    UpdateGrad(ptrGradVectorCellArray,
                               ptrArray,
                               'P');
                    // Calculate ve at grid points of main cell, interpolating curl E from lengs
                    UpdateVeGridsMain(ptrArray,
                                      ptrEVectorFaceArray_dual,
                                      0);
                    // Calculate the dB and E at grids points
                    if (div == div_max - 1)
                    {
                        BVectorGridsArrayUpdate(ptrArray,
                                                ptrBVectorFaceArray_main_backup);
                        EVectorGridsArrayUpdate(ptrArray,
                                                ptrEVectorFaceArray_dual_backup);
                    }
                    // Update gradient norm B
                    if (timeline == 0)
                    {
                        ptrGradVectorCellArray = ValueGradient(ptrArray,
                                                               ptrGradVectorCellArray,
                                                               ptrVolumeCellArray,
                                                               'B');
                        UpdateGrad(ptrGradVectorCellArray,
                                   ptrArray,
                                   'B');
                    }
                }
            }
        }
        double endTime4 = omp_get_wtime();
        cout << "Time elapsed for updating e-field & electron velocity: " << (endTime4 - endTime3) << endl;       

        // update datafiles of external model
        if (timeline % updataDataFiles == 0 && timeline != 0)
        {
            UpdateDatafileFromExternalDataForBottomBoundary(ptrArray,
                                                            timeline,
                                                            workdir,
                                                            iriData,
                                                            tiegcmData);
            cout << " Updated bottom boundary conditions\n";
            // YH 4/1, update potential at each gridpoints based on external model
            UpdatePotentialAtGrids(ptrArray, footArray);
        }

        // printout
        if (timeline % printTimePeriod == 0 && timeline!=0)
        {
            std::cout << " Output moments on grids at time step :" << timeline << std::endl; // grids print
            PrintOutHdf5(ptrArray,timeline);
            std::cout << " Output moments in cells at time step: " << timeline << std::endl; // cells print
            PrintOutHdf5Cells(ptrArrayCells, timeline);
        }
        //
        if (timeline % printTimePeriodParticles == 0 && timeline != 0) // for continue run
        {                                             //&& timeline!=0
            std::cout << " Output particle properties at time step: " << timeline << std::endl;
            PrintOutHdf5_Particles_Grids(timeline,
                                         ptrArrayCells,
                                         ptrArray);
        }
        //  reset particles at bottom layer
        if (timeline % resetBottomLayerParticles == 0 && timeline != 0)
        {
            ParticlesInCellsInitial(ptrArray,
                                    ptrArrayCells,
                                    1);
        }
        // reset pho and v
        if (timeline % updateInfoPeriod == 0 )
            ResetPhoVatGrids(ptrArray);
        // reset VelDist
        if (timeline % updateVelDist == 0 )
            ResetVelDistInCells(ptrArrayCells);
        //
        double endTime = omp_get_wtime();
        cout << "Time elapsed for output & reset: " << (endTime - endTime4) << endl;
    }
    //
    delete mem_gradVectorCellArray;
    //delete mem_EVectorCellArray;
    //delete mem_CurlBVectorCellArray;
    //delete mem_VeleVectorCellArray;
    //delete mem_VelWeightGridsArray;
    //delete mem_MassWeightGridsArray;
    //delete mem_DivBScalarVectorCellArray;
    //delete mem_posVectorCellArray;
    //delete ptrArray;
    //ptrParticlesList_H.clear();
    //ptrParticlesList_He.clear();
    //ptrParticlesList_O.clear();
    //ptrParticlesList_out_H.clear();
    //ptrParticlesList_out_He.clear();
    //ptrParticlesList_out_O.clear();
    //ptrParticlesList_H.shrink_to_fit();
    //ptrParticlesList_He.shrink_to_fit();
    //ptrParticlesList_O.shrink_to_fit();
    //ptrParticlesList_out_H.shrink_to_fit();
    //ptrParticlesList_out_He.shrink_to_fit();
    //ptrParticlesList_out_O.shrink_to_fit();
    ////
    //ptrParticlesListTemp_H.clear();
    //ptrParticlesListTemp_He.clear();
    //ptrParticlesListTemp_O.clear();
    //ptrParticlesListTemp_H.shrink_to_fit();
    //ptrParticlesListTemp_He.shrink_to_fit();
    //ptrParticlesListTemp_O.shrink_to_fit();
}

//************************************************************************
//************************************************************************
// FUNCTION
// Process control function, basicly it has two parts and repeat between them.
// 1 the particles parts that is mainly for moving the particles, and 2 the
// grid parts that is mainly for updating E B or other general variables on
// grid nodes.
//
// Procedures:
// Assume E and B on gridspoints are know before the first loop
// Particles Part I: Only need to go through the particles list once, for
// each particle class:
// 1. Get local E and B, details in particles.cpp, update locations
// and velocity of each particles, and return a structp of
// (ig, jg, kg, iw, jw, kw, vp) for each particles.
// 2. For each strucp of each particles, calculate related information of
// density and vi on each gridspoints. Old density and vi can be covered
// directly.
// Gridspoints Part II:
// Follwed Part I, assume density and vi are updated at gridspoints.
// Matrix of curl E and B need to be applied in this part to store curl
// E and B.
// The information update face by face.
// The calculation of curl E and B would be demonstrated elsewhere.
// 1. Ampere's Law: Calculate curl B at each gridspoints, and then calculate
// ve using Ampere's Law.
// 2. Electron's momentum equation: Calculate gradient of Pe, and then
// calculate new E at each gridspoints.
// 3. Faraday's Law: Calculate curl E, and then update B.
//************************************************************************
//************************************************************************
