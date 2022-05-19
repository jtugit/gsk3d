    #include <string>
    #include "H5Cpp.h"

    int total_thread_num;
    unsigned int *seed_thread;
    int openmpCores;

    int fieldsGridsLevel;
    int cellSize1;

    int fieldsGridsSize;
    int particlesGridsLevel;
    int particlesGridsSize;

    int coverGridsCellLevelBot;
    int coverGridsCellLevelTop;

    int cellBitlength;      //2^cellBitlength = mini cell number in one grid cell
    int blankBitlength;  //unused mini cell number

    double AltitudeMin;

    double LMin;

    int grid_domain;

    int iyr, mon, date, ihour, imin, isec;

    double Kp, Ap, F107, F107A;

    double tstep;

    int numberTimeStep;   //  length of running
    int updateInfoPeriod; // how many steps to recalculate E and other variable at gridspoints
                                   // A, base
    int updateVelDist;    // how many time steps to calculate velocity distribution in each cell
                                    // B, multiple of A
    int printTimePeriod;  // how many time steps to output state variables on grids
                                    //60000; // C, multiple of B
    int printTimePeriodParticles;  // how many time steps to output particles and grids infor
                                            // for continue run
                                            // D, multiple of C
    int printTimeVelDist_notused;    // print out vel distribution in cells not used
    int printDivBTimePeriod;   // calculate DIV of dB for checking non-div
    int div_max;                     // how many times E & dB are updated in one paricle time step
    int updataDataFiles;     // update bottom data files
    int collision_perPeriod; // collision per %d step

    int iniParticleNumberPerCell; //initial number of partilces in a cell H+, He+, O+

    int particlesNumMin; //60      minimum number of paricles. If less, split particles
    int particlesNumMax;    //120 max num of particles. If larger, combine some particles
// refilling, 0-initialize all domain; 1- only initialize only in bottom layer
    int refilling;
// Initialize particles: 0-create a new particles array; 1-read from file
    int continueParticles;
// particles velocity distribution parameters
    int velDistRange_para;//51; number of velocity cells in parallel direction
    int velDistRange_mu; //51; number of velocity cells in perpendicular direction
// speed range  in (m/s)
//  perp in [-velPerpMax,velPerpMax]; para in [-velParaMax, +velParaMax]  [60000]
    double velPerpMax;
    double velParaMax;

    double velSpace;

    int update_type;

    int initial_bot_type;
    // 0- no convection top boundary 1- with convection top boundary
    int initial_top_type;

    int collision_particles_control;      // 0 - off; 1- on

    double alphaPSD;
    double P_0;

//----------------------------------------------------------------------------------------------
    double LMax, LMid, LMin_maindomain, LMax_maindomain, LMin_centraldomain, LMax_centraldomain;

//domain 1--------------------------------------------------
    double ratioK;
    double logRatio;

    double LMax1;
    double LMin_maindomain1;
    double LMax_maindomain1;
    double LMin_centraldomain1;
    double LMax_centraldomain1;

//domain 2--------------------------------------------------
    double const_sinh1;
    double const_sinh2;

    double LMid2;
    double LMin_maindomain2;
    double LMin_centraldomain2;
    double LMax2;
    double LMax_maindomain2;
    double LMax_centraldomain2;

    double AAP[106], G[106], H[106], REC[106];

    double omega_earth_wx = 0.0;
    double omega_earth_wy = 0.0;
    double omega_earth_wz = omega_earth;

    double omega_wx = 0.0;
    double omega_wy = 0.0;
    double omega_wz = 0.0;
    
    int h5FileCheck = 0;
    char outpdir[151];

  //  uint_64 fieldsGridsSize = 1U << (fieldsGridsLevel);
   // uint_64 particlesGridsSize = 1U << particlesGridsLevel;
   // uint_64 radialGridsSize = 1U << (fieldsGridsLevel); // maybe not use

