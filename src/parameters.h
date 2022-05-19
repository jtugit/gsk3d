#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_
#include <iostream>
#include <cmath>
#include "H5Cpp.h"
using uint_64 = unsigned long long;
const double pi = acos(-1.0);
const double PI = acos(-1.0);
//************************************************************************
// random control
extern int total_thread_num;
extern unsigned int *seed_thread;
extern int openmpCores;
//************************************************************************
// earth physics parameters control

// Parameters constants of the earth
// Length of face in 6
const double length = 10007543.398; // 1/4 of qeuator perimeter of the Earth; 
// radius of earth
const double radius = 6371000.0;
// dipole moment for earth, T*m^3
const double dMoment = -8.0e+15;

// charge (C)
const double qi0 = 1.602e-19;
// mass of ion in unit of (kg) for Hydrogen
const double mi0_H = 1.66e-27;
// mass of ion in unit of kg for Helium
const double mi0_He = 6.64e-27;
// mass of ion in unit of kg for Oxygen
const double mi0_O = 2.656e-26;
// Vacuum permeability mu0
const double mu0 = 1.256637e-6;
// Const 1eV =
const double energy1eV = 1.60217662e-19;
// c light
const double lightSpeed = 3.0e+8;
// ion kT k =  1.38 x 10−23 J·K−1, T = 1000 K
const double ikT = 1.38e-20;
// Boltzmann_k constant in unit J·K−1
const double boltzmann_k = 1.38e-23;
// permittivity in free space
const double e_const = 8.85418782e-12;
// g (m / s2)
const double gravity = 9.8;
// angular velocity of Earth ( rad/s)
const double omega_earth = 7.292e-5;

// const for normalized case
const double normalized_N = 1.0e+10;
//************************************************************************
//************************************************************************
// grids size control

// These two levels are between 1 and 20, and particlesgridslevel is greater than the fieldgridslevel.
// The radialgridslevel is calculated from the fieldgridslevel to make the fieldgrids-cell similar to
// a cubic. Why do we want a cubic similar cell? The reason is that, the first, we want to use Z-curve
// to represent the location, therefore we need the same number grids on each three coordinates; the second,
// it is easy to calculate the weighting of each particle on the gridspoints, otherwise, with a cuboid cell,
// it may be hard to judge which particles in the cell will have weighting on the grids points.
// fieldGridsLevel max shouble better be 9
// particlesGridsLevel max should be 10 greater than fieldGridsLevel

extern int fieldsGridsLevel; // level of grids
                                //    const int particlesCellLevel = 10;       // better as 10
extern int particlesGridsLevel;  // max = 20;8
extern int cellSize1;
//    const int cellSize1 = 1 << (particlesGridsLevel - fieldsGridsLevel) <<  (particlesGridsLevel - fieldsGridsLevel) <<  (particlesGridsLevel - fieldsGridsLevel);
extern int fieldsGridsSize;
extern int particlesGridsSize;

const int tempGridsCellLevelBot = 1;  
extern int coverGridsCellLevelBot; //number of overlapped layers  at top boundary

const int tempGridsCellLevelTop = 1;
extern int coverGridsCellLevelTop;  //number of overlapped layers at bottom boundary

extern int cellBitlength;      //2^cellBitlength = mini cell number in one grid cell
extern int blankBitlength;  //unused mini cell number

// Determine the bottom L: LMin; Calculate the top L: LMax
const int totalFace = 6;
extern double AltitudeMin; // unit: m altitude of lower boundary
extern double LMin;

extern int grid_domain; // or 2 (max 2) number of domains in k-direction (radial direction) 
                           // can be double 

extern double LMax, LMid, LMin_maindomain, LMax_maindomain, LMin_centraldomain, LMax_centraldomain;

//*************************** domain 1
extern double ratioK;
extern double logRatio;

extern double LMax1;
extern double LMin_maindomain1;
extern double LMax_maindomain1;
extern double LMin_centraldomain1;
extern double LMax_centraldomain1;

//*************************** domain 2
const double grid_N1 = pow(2.0, 18) * 1.1;
extern double const_sinh1;
const double grid_N2 = pow(2.0, 20);
extern double const_sinh2;

extern double LMid2;
extern double LMin_maindomain2;
extern double LMin_centraldomain2;
extern double LMax2;
extern double LMax_maindomain2;
extern double LMax_centraldomain2;

//************************************************************************
//  For coordinates trans coefficients
extern double AAP[106], G[106], H[106], REC[106];
// solar wind in GST coordinate
//extern double VGSEX, VGSEY, VGSEZ;
// earth's rotation speed in GM coordinate
extern double omega_wx, omega_wy, omega_wz;
extern double omega_earth_wx, omega_earth_wy, omega_earth_wz;
//************************************************************************
//  For DATE
extern int iyr, mon, date, ihour, imin, isec;

//************************************************************************
// Geophysical parameters
extern double Kp, Ap, F107, F107A;

// run time control

// Simulation parameters 
extern double tstep;            //(unit) s
extern int numberTimeStep;   //  length of running
extern int updateInfoPeriod; // how many steps to recalculate E and other variable at gridspoints
                                   // A, base
extern int updateVelDist;    // how many time steps to calculate velocity distribution in each cell
                                    // B, multiple of A
extern int printTimePeriod;  // how many time steps to output state variables on grids
                                    //60000; // C, multiple of B
extern int printTimePeriodParticles;  // how many time steps to output particles and grids infor
                                            // for continue run
                                            // D, multiple of C
extern int printTimeVelDist_notused;    // print out vel distribution in cells not used
extern int printDivBTimePeriod;   // calculate DIV of dB for checking non-div
extern int div_max;                     // how many times E & dB are updated in one paricle time step
const int resetBottomLayerParticles = 5; // reset particle in bottom layer
extern int updataDataFiles;     // update bottom data files

extern int collision_perPeriod; // collision per %d step
const double isotropy_correction = 10.0;    // increse efficience of collision
//************************************************************************
//************************************************************************
// For printout

extern int h5FileCheck; // 0: create a new h5, 1: open exist h5
extern char outpdir[151];

//************************************************************************
//************************************************************************
// Initial particles control
const double mu_MaxwellDis = 2.0;
const double sigma_MaxwellDis = 0.15;
const double neutralizeRate = 0.5;       // rate of neutralization 
const double neutralizeRateCover = 0.00; // in the cover domain
const double updateParticlesRate = 0.0001;
// number density at base level ( / m^3) for H
const double N0_H = 100000000000.0; // may not used
// number density at base level for He  // may not used
const double N0_He = 100000000000.0; // may not used
// number density at base level for O   // may not used
const double N0_O = 100000000000.0; // may not used
// initial particle numbers per cell ( count) ****XXX   no greater than particlesNumMax
extern int iniParticleNumberPerCell; //initial number of partilces in a cell H+, He+, O+
// temp cell particle number per cell ( count)  // old
const int tempParticleNumberPerCellH = 50;
const int tempParticleNumberPerCellHe = 50;
const int tempParticleNumberPerCellO = 200;
//// particles number per cell    ***xxx initial
//const int particleNumberPerCell_H = 10;
//const int particleNumberPerCell_He = 10;
//const int particleNumberPerCell_O = 10;
// particles number range per cell  ****XXX important
extern int particlesNumMin; //60      minimum number of paricles. If less, split particles
extern int particlesNumMax;    //120 max num of particles. If larger, combine some particles
// refilling, 0-initialize all domain; 1- only initialize only in bottom layer
extern int refilling;
// Initialize particles: 0-create a new particles array; 1-read from file
extern int continueParticles;
// particles velocity distribution parameters
extern int velDistRange_para;//51; number of velocity cells in parallel direction
extern int velDistRange_mu; //51; number of velocity cells in perpendicular direction
// speed range  in (m/s)
//  perp in [-velPerpMax,velPerpMax]; para in [-velParaMax, +velParaMax]  [60000]
extern double velPerpMax;
extern double velParaMax;
extern double velSpace; //dV

//************************************************************************
//************************************************************************
// For some control
// Coordinate system
// 0- lab sys
// 1- rotating sys with earth
const int coordinate_rotate = 1;
// 0- no current 1- with current
extern int update_type;
// 0- test particle 1- not test
const int move_type = 0;   // should always be 0

// work only with update_type = 1
// 0- zero velocity
//    increased inner boundary velocity to normal
//    normal density for all
// 1- normal velocity for all ( not applied)
//    increased outer boundary velocity
//    normal density for all

// 0- no convection bot boundary
// 1- with convection bot boundary
extern int initial_bot_type;
// 0- no convection top boundary
// 1- with convection top boundary
extern int initial_top_type;

//time (s) for transition from 0 to full convection
const double botBoundaryInitialTimeStart = 0.0;
const double topBoundaryInitialTimeStart = 0.0;
const double botBoundaryInitialTime = 60.0; // 600
const double topBoundaryInitialTime = 60.0;

//************************************************************************
//************************************************************************
// For top boundary initialization, in degree
// c0_latitude > r0_latitude
const double r0_latitude = 65.0;
const double c0_latitude = 75.0;
const double t0_convection = 3600.0;
//************************************************************************
//************************************************************************
// For bot boundary initialization, in degree

const double EBot_const = 1.0e-6; //5.0e-7; e-filed boundary condition
                                  //**************************************
                                  //**************************************
                                  //    const double rho_max = 8.0e-17;
                                  //    const double rho_min = 2.0e-17;

const double ne_density_max = 1.0e+12;
const double ne_density_min = 1.0e+11;

const double ratioHight = 2000000; // in m

const double ratioH2000 = 0.56;
const double ratioHe2000 = 0.41; // not used
const double ratioO2000 = 0.3;

const double ratioH_bot = 0.04;
const double ratioHe_bot = 0.16;
const double ratioO_bot = 0.8;

const double ratioH_top = 0.98;
const double ratioHe_top = 0.0195;
const double ratioO_top = 0.0005;

const double initialDensityRate = 0.001;

//************************************************************************
//************************************************************************
const double collisionCurrent = 0.10;
const double diffusionCollisionCoefficient = 10.0; // for system stable
extern int collision_particles_control;      // 0 - off; 1- on

//************ wave-particle interaction**********************
extern double alphaPSD; // = 1.1 to 1.7; //power index of BBLF wave power spectrum
extern double P_0;      //1.0e-11 - 1.0e-9 sectrum power at H+ gyro-frequency 

//****************************Boris method force control******************
//
const int centrifugal_force = 1;
const int coriolis_force = 1;
const int gravity_force = 1;
const int mirror_force = 1;
const int electric_force = 1;
//*************************************************************
const double ss[400] = {0.01e0, 0.02e0, 0.03e0, 0.04e0, 0.05e0, 0.06e0, 0.07e0, 0.08e0, 0.09e0, 0.10e0,
                        0.11e0, 0.12e0, 0.13e0, 0.14e0, 0.15e0, 0.16e0, 0.17e0, 0.18e0, 0.19e0, 0.20e0,
                        0.21e0, 0.22e0, 0.23e0, 0.24e0, 0.25e0, 0.26e0, 0.27e0, 0.28e0, 0.29e0, 0.30e0,
                        0.31e0, 0.32e0, 0.33e0, 0.34e0, 0.35e0, 0.36e0, 0.37e0, 0.38e0, 0.39e0, 0.40e0,
                        0.41e0, 0.42e0, 0.43e0, 0.44e0, 0.45e0, 0.46e0, 0.47e0, 0.48e0, 0.49e0, 0.50e0,
                        0.51e0, 0.52e0, 0.53e0, 0.54e0, 0.55e0, 0.56e0, 0.57e0, 0.58e0, 0.59e0, 0.60e0,
                        0.61e0, 0.62e0, 0.63e0, 0.64e0, 0.65e0, 0.66e0, 0.67e0, 0.68e0, 0.69e0, 0.70e0,
                        0.71e0, 0.72e0, 0.73e0, 0.74e0, 0.75e0, 0.76e0, 0.77e0, 0.78e0, 0.79e0, 0.80e0,
                        0.81e0, 0.82e0, 0.83e0, 0.84e0, 0.85e0, 0.86e0, 0.87e0, 0.88e0, 0.89e0, 0.90e0,
                        0.91e0, 0.92e0, 0.93e0, 0.94e0, 0.95e0, 0.96e0, 0.97e0, 0.98e0, 0.99e0, 1.00e0,
                        1.01e0, 1.02e0, 1.03e0, 1.04e0, 1.05e0, 1.06e0, 1.07e0, 1.08e0, 1.09e0, 1.10e0,
                        1.11e0, 1.12e0, 1.13e0, 1.14e0, 1.15e0, 1.16e0, 1.17e0, 1.18e0, 1.19e0, 1.20e0,
                        1.21e0, 1.22e0, 1.23e0, 1.24e0, 1.25e0, 1.26e0, 1.27e0, 1.28e0, 1.29e0, 1.30e0,
                        1.31e0, 1.32e0, 1.33e0, 1.34e0, 1.35e0, 1.36e0, 1.37e0, 1.38e0, 1.39e0, 1.40e0,
                        1.41e0, 1.42e0, 1.43e0, 1.44e0, 1.45e0, 1.46e0, 1.47e0, 1.48e0, 1.49e0, 1.50e0,
                        1.51e0, 1.52e0, 1.53e0, 1.54e0, 1.55e0, 1.56e0, 1.57e0, 1.58e0, 1.59e0, 1.60e0,
                        1.61e0, 1.62e0, 1.63e0, 1.64e0, 1.65e0, 1.66e0, 1.67e0, 1.68e0, 1.69e0, 1.70e0,
                        1.71e0, 1.72e0, 1.73e0, 1.74e0, 1.75e0, 1.76e0, 1.77e0, 1.78e0, 1.79e0, 1.80e0,
                        1.81e0, 1.82e0, 1.83e0, 1.84e0, 1.85e0, 1.86e0, 1.87e0, 1.88e0, 1.89e0, 1.90e0,
                        1.91e0, 1.92e0, 1.93e0, 1.94e0, 1.95e0, 1.96e0, 1.97e0, 1.98e0, 1.99e0, 2.00e0,
                        2.01e0, 2.02e0, 2.03e0, 2.04e0, 2.05e0, 2.06e0, 2.07e0, 2.08e0, 2.09e0, 2.10e0,
                        2.11e0, 2.12e0, 2.13e0, 2.14e0, 2.15e0, 2.16e0, 2.17e0, 2.18e0, 2.19e0, 2.20e0,
                        2.21e0, 2.22e0, 2.23e0, 2.24e0, 2.25e0, 2.26e0, 2.27e0, 2.28e0, 2.29e0, 2.30e0,
                        2.31e0, 2.32e0, 2.33e0, 2.34e0, 2.35e0, 2.36e0, 2.37e0, 2.38e0, 2.39e0, 2.40e0,
                        2.41e0, 2.42e0, 2.43e0, 2.44e0, 2.45e0, 2.46e0, 2.47e0, 2.48e0, 2.49e0, 2.50e0,
                        2.51e0, 2.52e0, 2.53e0, 2.54e0, 2.55e0, 2.56e0, 2.57e0, 2.58e0, 2.59e0, 2.60e0,
                        2.61e0, 2.62e0, 2.63e0, 2.64e0, 2.65e0, 2.66e0, 2.67e0, 2.68e0, 2.69e0, 2.70e0,
                        2.71e0, 2.72e0, 2.73e0, 2.74e0, 2.75e0, 2.76e0, 2.77e0, 2.78e0, 2.79e0, 2.80e0,
                        2.81e0, 2.82e0, 2.83e0, 2.84e0, 2.85e0, 2.86e0, 2.87e0, 2.88e0, 2.89e0, 2.90e0,
                        2.91e0, 2.92e0, 2.93e0, 2.94e0, 2.95e0, 2.96e0, 2.97e0, 2.98e0, 2.99e0, 3.00e0,
                        3.01e0, 3.02e0, 3.03e0, 3.04e0, 3.05e0, 3.06e0, 3.07e0, 3.08e0, 3.09e0, 3.10e0,
                        3.11e0, 3.12e0, 3.13e0, 3.14e0, 3.15e0, 3.16e0, 3.17e0, 3.18e0, 3.19e0, 3.20e0,
                        3.21e0, 3.22e0, 3.23e0, 3.24e0, 3.25e0, 3.26e0, 3.27e0, 3.28e0, 3.29e0, 3.30e0,
                        3.31e0, 3.32e0, 3.33e0, 3.34e0, 3.35e0, 3.36e0, 3.37e0, 3.38e0, 3.39e0, 3.40e0,
                        3.41e0, 3.42e0, 3.43e0, 3.44e0, 3.45e0, 3.46e0, 3.47e0, 3.48e0, 3.49e0, 3.50e0,
                        3.51e0, 3.52e0, 3.53e0, 3.54e0, 3.55e0, 3.56e0, 3.57e0, 3.58e0, 3.59e0, 3.60e0,
                        3.61e0, 3.62e0, 3.63e0, 3.64e0, 3.65e0, 3.66e0, 3.67e0, 3.68e0, 3.69e0, 3.70e0,
                        3.71e0, 3.72e0, 3.73e0, 3.74e0, 3.75e0, 3.76e0, 3.77e0, 3.78e0, 3.79e0, 3.80e0,
                        3.81e0, 3.82e0, 3.83e0, 3.84e0, 3.85e0, 3.86e0, 3.87e0, 3.88e0, 3.89e0, 3.90e0,
                        3.91e0, 3.92e0, 3.93e0, 3.94e0, 3.95e0, 3.96e0, 3.97e0, 3.98e0, 3.99e0, 4.00e0};
const double AA[400] = {0.1005008e+03, 0.5050167e+02, 0.3383583e+02, 0.2550333e+02, 0.2050417e+02,
                        0.1717167e+02, 0.1479155e+02, 0.1300667e+02, 0.1161861e+02, 0.1050833e+02,
                        0.9600073e+01, 0.8843328e+01, 0.8203128e+01, 0.7654494e+01, 0.7179102e+01,
                        0.6763205e+01, 0.6396285e+01, 0.6070154e+01, 0.5778342e+01, 0.5515671e+01,
                        0.5277941e+01, 0.5061716e+01, 0.4864155e+01, 0.4682890e+01, 0.4515930e+01,
                        0.4361594e+01, 0.4218447e+01, 0.4085265e+01, 0.3960991e+01, 0.3844712e+01,
                        0.3735635e+01, 0.3633067e+01, 0.3536400e+01, 0.3445102e+01, 0.3358701e+01,
                        0.3276778e+01, 0.3198964e+01, 0.3124926e+01, 0.3054370e+01, 0.2987030e+01,
                        0.2922667e+01, 0.2861066e+01, 0.2802032e+01, 0.2745390e+01, 0.2690979e+01,
                        0.2638654e+01, 0.2588282e+01, 0.2539742e+01, 0.2492922e+01, 0.2447721e+01,
                        0.2404044e+01, 0.2361806e+01, 0.2320927e+01, 0.2281333e+01, 0.2242956e+01,
                        0.2205734e+01, 0.2169606e+01, 0.2134520e+01, 0.2100425e+01, 0.2067272e+01,
                        0.2035019e+01, 0.2003623e+01, 0.1973047e+01, 0.1943254e+01, 0.1914210e+01,
                        0.1885884e+01, 0.1858245e+01, 0.1831266e+01, 0.1804921e+01, 0.1779183e+01,
                        0.1754030e+01, 0.1729439e+01, 0.1705389e+01, 0.1681861e+01, 0.1658835e+01,
                        0.1636293e+01, 0.1614220e+01, 0.1592597e+01, 0.1571412e+01, 0.1550647e+01,
                        0.1530291e+01, 0.1510329e+01, 0.1490750e+01, 0.1471540e+01, 0.1452690e+01,
                        0.1434187e+01, 0.1416021e+01, 0.1398183e+01, 0.1380663e+01, 0.1363452e+01,
                        0.1346540e+01, 0.1329921e+01, 0.1313585e+01, 0.1297525e+01, 0.1281733e+01,
                        0.1266203e+01, 0.1250927e+01, 0.1235900e+01, 0.1221114e+01, 0.1206564e+01,
                        0.1192244e+01, 0.1178148e+01, 0.1164271e+01, 0.1150607e+01, 0.1137152e+01,
                        0.1123901e+01, 0.1110849e+01, 0.1097991e+01, 0.1085324e+01, 0.1072842e+01,
                        0.1060542e+01, 0.1048420e+01, 0.1036471e+01, 0.1024694e+01, 0.1013082e+01,
                        0.1001635e+01, 0.9903465e+00, 0.9792152e+00, 0.9682374e+00, 0.9574099e+00,
                        0.9467299e+00, 0.9361945e+00, 0.9258008e+00, 0.9155462e+00, 0.9054280e+00,
                        0.8954436e+00, 0.8855906e+00, 0.8758665e+00, 0.8662689e+00, 0.8567956e+00,
                        0.8474444e+00, 0.8382130e+00, 0.8290994e+00, 0.8201015e+00, 0.8112173e+00,
                        0.8024448e+00, 0.7937821e+00, 0.7852274e+00, 0.7767789e+00, 0.7684348e+00,
                        0.7601934e+00, 0.7520529e+00, 0.7440118e+00, 0.7360685e+00, 0.7282214e+00,
                        0.7204689e+00, 0.7128097e+00, 0.7052421e+00, 0.6977648e+00, 0.6903765e+00,
                        0.6830757e+00, 0.6758611e+00, 0.6687314e+00, 0.6616853e+00, 0.6547217e+00,
                        0.6478392e+00, 0.6410367e+00, 0.6343130e+00, 0.6276670e+00, 0.6210976e+00,
                        0.6146036e+00, 0.6081839e+00, 0.6018376e+00, 0.5955636e+00, 0.5893609e+00,
                        0.5832285e+00, 0.5771653e+00, 0.5711706e+00, 0.5652432e+00, 0.5593824e+00,
                        0.5535872e+00, 0.5478567e+00, 0.5421901e+00, 0.5365864e+00, 0.5310450e+00,
                        0.5255649e+00, 0.5201453e+00, 0.5147855e+00, 0.5094847e+00, 0.5042421e+00,
                        0.4990569e+00, 0.4939285e+00, 0.4888560e+00, 0.4838389e+00, 0.4788763e+00,
                        0.4739677e+00, 0.4691122e+00, 0.4643093e+00, 0.4595584e+00, 0.4548586e+00,
                        0.4502095e+00, 0.4456105e+00, 0.4410607e+00, 0.4365598e+00, 0.4321071e+00,
                        0.4277020e+00, 0.4233439e+00, 0.4190322e+00, 0.4147665e+00, 0.4105461e+00,
                        0.4063706e+00, 0.4022393e+00, 0.3981518e+00, 0.3941076e+00, 0.3901060e+00,
                        0.3861467e+00, 0.3822292e+00, 0.3783529e+00, 0.3745173e+00, 0.3707221e+00,
                        0.3669667e+00, 0.3632507e+00, 0.3595736e+00, 0.3559349e+00, 0.3523343e+00,
                        0.3487713e+00, 0.3452455e+00, 0.3417564e+00, 0.3383037e+00, 0.3348868e+00,
                        0.3315055e+00, 0.3281594e+00, 0.3248479e+00, 0.3215708e+00, 0.3183276e+00,
                        0.3151180e+00, 0.3119416e+00, 0.3087980e+00, 0.3056869e+00, 0.3026079e+00,
                        0.2995607e+00, 0.2965448e+00, 0.2935601e+00, 0.2906060e+00, 0.2876823e+00,
                        0.2847887e+00, 0.2819248e+00, 0.2790903e+00, 0.2762849e+00, 0.2735082e+00,
                        0.2707600e+00, 0.2680400e+00, 0.2653477e+00, 0.2626831e+00, 0.2600456e+00,
                        0.2574351e+00, 0.2548513e+00, 0.2522938e+00, 0.2497625e+00, 0.2472569e+00,
                        0.2447769e+00, 0.2423222e+00, 0.2398924e+00, 0.2374874e+00, 0.2351068e+00,
                        0.2327505e+00, 0.2304181e+00, 0.2281094e+00, 0.2258242e+00, 0.2235621e+00,
                        0.2213230e+00, 0.2191067e+00, 0.2169128e+00, 0.2147411e+00, 0.2125914e+00,
                        0.2104636e+00, 0.2083572e+00, 0.2062722e+00, 0.2042083e+00, 0.2021652e+00,
                        0.2001429e+00, 0.1981409e+00, 0.1961592e+00, 0.1941975e+00, 0.1922557e+00,
                        0.1903334e+00, 0.1884305e+00, 0.1865469e+00, 0.1846822e+00, 0.1828364e+00,
                        0.1810091e+00, 0.1792003e+00, 0.1774097e+00, 0.1756371e+00, 0.1738824e+00,
                        0.1721454e+00, 0.1704259e+00, 0.1687236e+00, 0.1670385e+00, 0.1653704e+00,
                        0.1637190e+00, 0.1620842e+00, 0.1604659e+00, 0.1588638e+00, 0.1572779e+00,
                        0.1557078e+00, 0.1541536e+00, 0.1526150e+00, 0.1510918e+00, 0.1495839e+00,
                        0.1480911e+00, 0.1466133e+00, 0.1451504e+00, 0.1437021e+00, 0.1422684e+00,
                        0.1408490e+00, 0.1394439e+00, 0.1380529e+00, 0.1366758e+00, 0.1353125e+00,
                        0.1339629e+00, 0.1326268e+00, 0.1313041e+00, 0.1299946e+00, 0.1286983e+00,
                        0.1274149e+00, 0.1261444e+00, 0.1248866e+00, 0.1236414e+00, 0.1224087e+00,
                        0.1211883e+00, 0.1199802e+00, 0.1187841e+00, 0.1176000e+00, 0.1164277e+00,
                        0.1152672e+00, 0.1141182e+00, 0.1129808e+00, 0.1118547e+00, 0.1107399e+00,
                        0.1096363e+00, 0.1085437e+00, 0.1074620e+00, 0.1063911e+00, 0.1053309e+00,
                        0.1042813e+00, 0.1032422e+00, 0.1022135e+00, 0.1011950e+00, 0.1001868e+00,
                        0.9918860e-01, 0.9820038e-01, 0.9722203e-01, 0.9625345e-01, 0.9529455e-01,
                        0.9434522e-01, 0.9340537e-01, 0.9247491e-01, 0.9155374e-01, 0.9064176e-01,
                        0.8973888e-01, 0.8884502e-01, 0.8796008e-01, 0.8708398e-01, 0.8621661e-01,
                        0.8535791e-01, 0.8450777e-01, 0.8366611e-01, 0.8283286e-01, 0.8200791e-01,
                        0.8119120e-01, 0.8038263e-01, 0.7958213e-01, 0.7878962e-01, 0.7800501e-01,
                        0.7722823e-01, 0.7645919e-01, 0.7569783e-01, 0.7494405e-01, 0.7419780e-01,
                        0.7345898e-01, 0.7272754e-01, 0.7200338e-01, 0.7128645e-01, 0.7057666e-01,
                        0.6987395e-01, 0.6917825e-01, 0.6848948e-01, 0.6780758e-01, 0.6713248e-01,
                        0.6646410e-01, 0.6580239e-01, 0.6514727e-01, 0.6449868e-01, 0.6385656e-01,
                        0.6322084e-01, 0.6259145e-01, 0.6196833e-01, 0.6135143e-01, 0.6074067e-01,
                        0.6013600e-01, 0.5953735e-01, 0.5894467e-01, 0.5835789e-01, 0.5777696e-01,
                        0.5720182e-01, 0.5663240e-01, 0.5606867e-01, 0.5551054e-01, 0.5495798e-01};

const double cab[3][3] = {{0.97048,0.379094,0.000947734},
                          {0.379094,0.060655,0.000947734},
                          {0.000947734,0.000947734,0.00379094}};

#endif