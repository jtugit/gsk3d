/******************************************************************************
 *  input_param
 *   Input various parameters
 *
 * Jiannan Tu 1/24/2022
 ******************************************************************************/
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>

using namespace std;

#include "parameters.h"

int input_param(string workdir, string &ParticleFileName)
{
    string fname;

    fname = workdir + "/gsk3din.dat";

    fstream infstr(fname, fstream::in);
    if(!infstr) {
        cout << "Can't open file: " << fname << endl;
        cout << "Exit from function input_param" << endl;
        return -1;
    }

    infstr >> continueParticles;
    infstr.ignore(200, '\n');
    infstr >> total_thread_num;
    infstr.ignore(200, '\n');
    infstr >> iyr;
    infstr.ignore(200, '\n');
    infstr >> mon;
    infstr.ignore(200, '\n');
    infstr >> date;
    infstr.ignore(200, '\n');
    infstr >> ihour;
    infstr.ignore(200, '\n');
    infstr >> imin;
    infstr.ignore(200, '\n');
    infstr >> isec;
    infstr.ignore(200, '\n');
    infstr >> F107;
    infstr.ignore(200, '\n');
    infstr >> F107A;
    infstr.ignore(200, '\n');
    infstr >> Kp;
    infstr.ignore(200, '\n');
    infstr >> Ap;
    infstr.ignore(200, '\n');
    infstr >> fieldsGridsLevel;
    infstr.ignore(200, '\n');
    infstr >> particlesGridsLevel;
    infstr.ignore(200, '\n');
    infstr >> coverGridsCellLevelBot;
    infstr.ignore(200, '\n');
    infstr >> coverGridsCellLevelTop;
    infstr.ignore(200, '\n');
    infstr >> AltitudeMin;
    infstr.ignore(200, '\n');
    infstr >> grid_domain;
    infstr.ignore(200, '\n');
    infstr >> tstep;
    infstr.ignore(200, '\n');
    infstr >> numberTimeStep;
    infstr.ignore(200, '\n');
    infstr >> updateInfoPeriod;
    infstr.ignore(200, '\n');
    infstr >> updateVelDist;
    infstr.ignore(200, '\n');
    infstr >> printTimePeriod;
    infstr.ignore(200, '\n');
    infstr >> printTimePeriodParticles;
    infstr.ignore(200, '\n');
    infstr >> printTimeVelDist_notused;
    infstr.ignore(200, '\n');
    infstr >> printDivBTimePeriod;
    infstr.ignore(200, '\n');
    infstr >> div_max;
    infstr.ignore(200, '\n');
    infstr >> updataDataFiles;
    infstr.ignore(200, '\n');
    infstr >> collision_perPeriod;
    infstr.ignore(200, '\n');
    infstr >> h5FileCheck;
    infstr.ignore(200, '\n');
    infstr >> iniParticleNumberPerCell;
    infstr.ignore(200, '\n');
    infstr >> particlesNumMin;
    infstr.ignore(200, '\n');
    infstr >> particlesNumMax;
    infstr.ignore(200, '\n');
    infstr >> refilling;
    infstr.ignore(200, '\n');
    infstr >> velDistRange_para;
    infstr.ignore(200, '\n');
    infstr >> velDistRange_mu;
    infstr.ignore(200, '\n');
    infstr >> velPerpMax;
    infstr.ignore(200, '\n');
    infstr >> velParaMax;
    infstr.ignore(200, '\n');
    infstr >> update_type;
    infstr.ignore(200, '\n');
    infstr >> initial_bot_type;
    infstr.ignore(200, '\n');
    infstr >> initial_top_type;
    infstr.ignore(200, '\n');
    infstr >> collision_particles_control;
    infstr.ignore(200, '\n');
    infstr >> alphaPSD;
    infstr.ignore(200, '\n');
    infstr >> P_0;
    infstr.ignore(200, '\n');
    infstr >> outpdir;
    infstr.ignore(200, '\n');
    infstr >> ParticleFileName;

    infstr.close();

    ParticleFileName = workdir + "/ExternalData/" + ParticleFileName;

    cellSize1 = 1 << (particlesGridsLevel - fieldsGridsLevel);
    fieldsGridsSize = 1 << fieldsGridsLevel;
    particlesGridsSize = 1 << particlesGridsLevel;

    cellBitlength = 61 - fieldsGridsLevel * 3;      //2^cellBitlength = mini cell number in one grid cell
    blankBitlength = 61 - particlesGridsLevel * 3;  //unused mini cell number

    LMin = 1 + AltitudeMin / radius;

    velSpace = velPerpMax / ( velDistRange_mu - 1); //dV

    ratioK = 1.0 + (length / particlesGridsSize) / radius;
    logRatio = log10(1.0 + (length / particlesGridsSize) / radius);

    LMax1 = LMin * pow(ratioK, particlesGridsSize);
    LMin_maindomain1 = LMin * pow(ratioK, (tempGridsCellLevelBot - coverGridsCellLevelBot) * cellSize1);
    LMax_maindomain1 = LMin * pow(ratioK, (fieldsGridsSize - tempGridsCellLevelTop + coverGridsCellLevelTop) * cellSize1);
    LMin_centraldomain1 = LMin * pow(ratioK, tempGridsCellLevelBot *cellSize1);
    LMax_centraldomain1 = LMin * pow(ratioK, (fieldsGridsSize - tempGridsCellLevelTop) * cellSize1);

    const_sinh1 = sinh(particlesGridsSize / grid_N1);
    const_sinh2 = sinh(particlesGridsSize / grid_N2);

    LMid2 = LMin + 0.5;
    LMin_maindomain2 = LMin + (LMid2 - LMin) * sinh((tempGridsCellLevelBot - coverGridsCellLevelBot) * cellSize1 / grid_N1) / const_sinh1;
    LMin_centraldomain2 = LMin + (LMid2 - LMin) * sinh(tempGridsCellLevelBot * cellSize1 / grid_N1) / const_sinh1;
    LMax2 = LMin + 2.5;
    LMax_maindomain2 = LMid2 + (LMax2 - LMid2) * sinh((fieldsGridsSize - tempGridsCellLevelTop + coverGridsCellLevelTop) * cellSize1 / grid_N2) / const_sinh2;
    LMax_centraldomain2 = LMid2 + (LMax2 - LMid2) * sinh((fieldsGridsSize - tempGridsCellLevelTop) * cellSize1 / grid_N2) / const_sinh2;

    seed_thread = new unsigned int [total_thread_num];

    return 0;
}

//********* output run parameters
int output_param(string ParticleFileName)
{
    string fname;

    fname = string(outpdir) + "/gsk3drun_param.log";

    fstream outfstr(fname, fstream::out);
    if(!outfstr) {
        cout << "Can't open file " << fname << endl;
        return -1;
    }

    outfstr << "Values of input parameters:" << endl; 
    outfstr << continueParticles << "   -- continueParticles: new or continue run" << endl;
    outfstr << total_thread_num << "   -- openmpCores: number of openmp cores" << endl;
    outfstr << iyr << "    -- iyr: year" << endl;
    outfstr << mon << "    -- mon: month" << endl;
    outfstr << date << "    -- date: date" <<endl;
    outfstr << ihour << "    -- ihour: UT hour" <<endl;
    outfstr << imin << "    -- imin: UT minutes" <<endl;
    outfstr << isec << "    -- isec: UT secods" << endl;
    outfstr << fieldsGridsLevel << "    -- fieldsGridsLevel" << endl ;
    outfstr << particlesGridsLevel << "   -- particlesGridsLevel" << endl;
    outfstr << coverGridsCellLevelBot << "    -- coverGridsCellLevelBot" << endl;
    outfstr << coverGridsCellLevelTop << "    -- coverGridsCellLevelTop" << endl;
    outfstr << AltitudeMin << "    -- AltitudeMin" << endl;
    outfstr << grid_domain << "    -- grid_domain" << endl;
    outfstr << tstep << "    -- tstep: time step (s)" << endl;
    outfstr << numberTimeStep << "    -- numberTimeStep" << endl;
    outfstr << updateInfoPeriod << "    -- updateInfoPeriod" << endl;
    outfstr << updateVelDist << "    -- updateVelDist" << endl;
    outfstr << printTimePeriod << "    -- printTimePeriod" << endl;
    outfstr << printTimePeriodParticles << "    -- printTimePeriodParticles" << endl;
    outfstr << printTimeVelDist_notused << "    -- printTimeVelDist_notused" << endl;
    outfstr << printDivBTimePeriod << "    -- printDivBTimePeriod" << endl;
    outfstr << div_max << "    -- div_max" << endl;
    outfstr << updataDataFiles << "    -- updataDataFiles" << endl;
    outfstr << collision_perPeriod << "    -- collision_perPeriod" << endl;
    outfstr << h5FileCheck << "    -- h5FileCheck" << endl;
    outfstr << iniParticleNumberPerCell << "    -- iniParticleNumberPerCell" << endl;
    outfstr << particlesNumMin << "    -- particlesNumMin" << endl;
    outfstr << particlesNumMax << "    -- particlesNumMax" << endl;
    outfstr << refilling << "    -- refilling" << endl;
    outfstr << velDistRange_para << "    -- velDistRange_para" << endl;
    outfstr << velDistRange_mu << "    -- velDistRange_mu" << endl;
    outfstr << velPerpMax << "    -- velPerpMax" << endl;
    outfstr << velParaMax << "    -- velParaMax" << endl;
    outfstr << initial_bot_type << "    -- initial_bot_type" << endl;
    outfstr << initial_top_type << "    -- initial_top_type" << endl;
    outfstr << collision_particles_control << "    -- collision_particles_control" << endl;
    outfstr << alphaPSD << "    -- power index of BBLF wave spectrum" << endl;
    outfstr << P_0 << "    -- wave spectrum power at H+ gyro-frequency" << endl;
    if (continueParticles)
        outfstr << ParticleFileName << "    -- Continue run from this particle data file" << endl;

    outfstr.close();

    return 0;
}
