#include <iostream>
#include <cmath>
#include "parameters.h"
#include "parameters_def.h"
#include "fieldsgrids.h"
#include "particles.h"
#include "vector3.h"
#include "module.h"
#include "vector"
#include <omp.h>
#include <unistd.h>

using namespace std;

int input_param(string workdi, string &ParticleFileName);
int output_param(string ParticleFileName);

int main()
{
    std::cout << " Welcome to 3-D Generalized Semi-Kinetic (GSK3D) Model " << std::endl;
    extern double LMax, LMid, LMin_maindomain, LMax_maindomain, LMin_centraldomain, LMax_centraldomain;
    double VGSEX = -400.0;
    double VGSEY = 0.0;
    double VGSEZ = 0.0;
    double pdyn=0.5, dst=-30.0, Byimf=0.0, Bzimf=-2.0;
    
    char   buff[150];
    string workdir, ParticleFileName;

    if (getcwd(buff, sizeof(buff)) == NULL) {
        perror("getcwd() error");
        exit(EXIT_FAILURE);
    }
    size_t pos = string(buff).find("/src");
    workdir = string(buff).substr(0,pos);

    //get current time
    time_t start_t=time(NULL);

    if (input_param(workdir, ParticleFileName) < 0) exit(EXIT_FAILURE);
    if (output_param(ParticleFileName) < 0) exit(EXIT_FAILURE);

    if( grid_domain == 1)
    {   
        LMax = LMax1;
        LMid = LMax1;
        LMin_maindomain = LMin_maindomain1;
        LMax_maindomain = LMax_maindomain1;
        LMin_centraldomain = LMin_centraldomain1;
        LMax_centraldomain = LMax_centraldomain1;
    } else if( grid_domain == 2 )
    {
        LMax = LMax2;
        LMid = LMid2;
        LMin_maindomain = LMin_maindomain2;
        LMax_maindomain = LMax_maindomain2;
        LMin_centraldomain = LMin_centraldomain2;
        LMax_centraldomain = LMax_centraldomain2;
    } 
    else
    {
        std::cout << " grid_domain error \n";
        std::cin.get();
    }

    std::cout << " The simulation is for the domain from L= <" << LMin << "> to <" << LMax <<  "> " <<std::endl;
    std::cout << " The simulation is for the main-domain from L= <" << LMin_maindomain << "> to <" << LMax_maindomain <<  "> " <<std::endl;
    std::cout << " The simulation is for the central-domain from L= <" << LMin_centraldomain << "> to <" << LMax_centraldomain <<  "> " <<std::endl;

    std::cout << " GridsSize of cubic on each face of total six faces is " << fieldsGridsSize << std::endl;
    std::cout << " Time step is <" << tstep << "> second" << endl;
    std::cout << " fields update every <" << updateInfoPeriod << "> timestep" << endl;
    std::cout << " output moments & fields every <" << printTimePeriod << "> timesteps" << std::endl;

    if (continueParticles == 0) {
        std::cout << " In each cell, the initial number of particles (H+, He+, O+) is <" << iniParticleNumberPerCell << ">" << std::endl;
        if( move_type == 0) {
            std::cout << " Initialized particles are at rest, " << std::endl;
        }else {
            std::cout << " Initiailized prticles are rotation with Earth, " << std::endl;
        }
    }
    else cout << " Continuous from a previous run" << endl;
    
    if( update_type ==0) {
        std::cout << " No currents. Wave-particle interactions neglected." << std::endl;
    }else {
        std::cout << " Current included. Wave-particle interactions considered." << std::endl;
    }

    if( initial_bot_type == 1) {
        std::cout << " Rotation with Earth at bottom boundary applied. Applied rotation E fields will " <<
        "propagate to entire domain. " << std::endl;
    }
    if( initial_top_type == 1) {
        std::cout << " Two-cell convention pattern on poles map to outside/top boundary, and applied "
        << "E fields will propagate to entire domain. " << std::endl;
    }

    std::cout << " Number of threads: " << total_thread_num << endl;
    std::cout << " Total Number of  time steps to run: " << numberTimeStep << endl;

    // t= 0. put in initial condition and create particles initial condition
    // t = n. 
    // Go through the particles list, 
    // 1. calculate the local E and B
    // 2. calculate the new velocity using Boris' Method
    // 3. update the location of the particles
    // 4. update the weighting of N and V on grids points
    // t = n+1 , repeat t =n

    ProcessFunc(VGSEX, VGSEY, VGSEZ, pdyn, dst, Byimf, Bzimf, workdir);

    time_t end_t=time(NULL);
    cout<<"Program completed! Wall time used: " <<difftime(end_t, start_t)<<" s"<<endl;

    return 0;
}