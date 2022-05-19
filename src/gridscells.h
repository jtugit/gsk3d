#ifndef _FIELDSCELLS_H_
#define _FIELDSCELLS_H_
/// Field grids information
#include <iostream>
#include <vector>
#include "parameters.h"
#include "vector3.h"
#include "particles.h"
#include "module_base.h"
#include "module_2.h"
#include <algorithm>
using std::max;
using std::vector;

class GridsCells
{
public:
    friend class Particles;
    friend class Vector3;

    //
    inline void InitialParticlesCells()
    {
        particles_H = new vector<Particles>();
        particles_H->reserve(particlesNumMax * 2);
        particles_He = new vector<Particles>();
        particles_He->reserve(particlesNumMax * 2);
        particles_O = new vector<Particles>();
        particles_O->reserve(particlesNumMax * 2);
    }
    //
    inline void InitialParVelDistArray()
    {
        velDist_H = new vector<vector<double> >(velDistRange_para, vector<double>(velDistRange_mu));
        //    velDist_H ->reserve(velDistRange);
        velDist_He = new vector<vector<double> >(velDistRange_para, vector<double>(velDistRange_mu));
        velDist_O = new vector<vector<double> >(velDistRange_para, vector<double>(velDistRange_mu));
    }
    inline void InitialRandomIndexPar()
    {
        par_random = new double[particlesNumMax];
    }
    //
    inline int ParCapacity()
    {
        return particles_H->capacity();
    }
    inline int ParSize()
    {
        return particles_H->size();
    }
    inline void SetParticlesCells(vector<Particles> *par_H,
                                  vector<Particles> *par_He,
                                  vector<Particles> *par_O)
    {
        particles_H = par_H;
        particles_He = par_He;
        particles_O = par_O;
    }
    // //
    // inline void SetParRandom(vector<double> parR)
    // {
    //     par_random = &parR;
    // }
    //
    inline void SetVolume(double vol)
    {
        volume = vol;
    }
    //
    inline void SetB3(Vector3 other)
    {
        b3Cell = other;
    }
    //
    inline Vector3 B3Cell()
    {
        return b3Cell;
    }
    //
    inline double Volume()
    {
        return volume;
    }
    inline vector<Particles> *Particles_H()
    {
        return particles_H;
    }
    //
    inline void Particles_H_pushback(Particles tempPar)
    {
        particles_H->push_back(tempPar);
    }
    //
    inline vector<Particles> *Particles_He()
    {
        return particles_He;
    }
    //
    inline vector<Particles> *Particles_O()
    {
        return particles_O;
    }
    //
    inline vector<vector<double>> *VelDist_H()
    {
        return velDist_H;
    }
    //
    inline vector<vector<double>> *VelDist_He()
    {
        return velDist_He;
    }
    //
    inline vector<vector<double>> *VelDist_O()
    {
        return velDist_O;
    }
    //***********************************************************************
    inline void TestPar()
    {
        std::cout << " mu " << (*particles_H)[0].MagneticIvarient() << " "
                  << " vel " << (*particles_H)[0].VelParticles().norm();
    }
    //***********************************************************************
    inline double TestCollMoment(GridsPoints *****ptrArray)
    {
        Vector3 sumMoment = {0.0, 0.0, 0.0};
        Vector3 b3Par = {0.0, 0.0, 0.0};
        for (auto iter = particles_H->begin(); iter != particles_H->end(); ++iter)
        {
            b3Par = iter->B3atParticles(ptrArray);
            sumMoment = sumMoment.PlusProduct(iter->VelCollParticles(b3Par, mi0_H).ScaleProduct(mi0_H * iter->WeightNi()));
        }
        for (auto iter = particles_He->begin(); iter != particles_He->end(); ++iter)
        {
            b3Par = iter->B3atParticles(ptrArray);
            sumMoment = sumMoment.PlusProduct(iter->VelCollParticles(b3Par, mi0_He).ScaleProduct(mi0_He * iter->WeightNi()));
        }
        for (auto iter = particles_O->begin(); iter != particles_O->end(); ++iter)
        {
            b3Par = iter->B3atParticles(ptrArray);
            sumMoment = sumMoment.PlusProduct(iter->VelCollParticles(b3Par, mi0_O).ScaleProduct(mi0_O * iter->WeightNi()));
        }
        //    std::cout<< " vel " << sumMoment.DotProduct(b3Cell) ;
        return sumMoment.DotProduct(b3Cell);
    }

    //***********************************************************************
    inline double TestCollEnergy(GridsPoints *****ptrArray, double &a, double &b, double &c)
    {
        Vector3 b3Par = {0.0, 0.0, 0.0};
        double sumEnergy = 0.0;
        for (auto iter = particles_H->begin(); iter != particles_H->end(); ++iter)
        {
            b3Par = iter->B3atParticles(ptrArray);
            sumEnergy = sumEnergy + iter->VelCollParticles(b3Par, mi0_H).norm2() * mi0_H * iter->WeightNi();
        }
        a = sumEnergy;
        for (auto iter = particles_He->begin(); iter != particles_He->end(); ++iter)
        {
            b3Par = iter->B3atParticles(ptrArray);
            sumEnergy = sumEnergy + iter->VelCollParticles(b3Par, mi0_He).norm2() * mi0_He * iter->WeightNi();
        }
        b = sumEnergy - a;
        for (auto iter = particles_O->begin(); iter != particles_O->end(); ++iter)
        {
            b3Par = iter->B3atParticles(ptrArray);
            sumEnergy = sumEnergy + iter->VelCollParticles(b3Par, mi0_O).norm2() * mi0_O * iter->WeightNi();
        }
        c = sumEnergy - a - b;
        //    std::cout<< " vel " << sumMoment.DotProduct(b3Cell) ;
        return sumEnergy;
    }

    //***********************************************************************
    // Removing slot of out-cell particles
    inline void SpeciesRearrange(int a)
    {
        Particles tempP;
        vector<Particles> *par_s = NULL;
        if (a == 0)
        {
            par_s = particles_H;
        }
        else if (a == 1)
        {
            par_s = particles_He;
        }
        else if (a == 2)
        {
            par_s = particles_O;
        }
        // consider the first particles may be out,
        // initial set up slotPosition[] are all zero,
        // so the index of begin() is set as 1 instead of 0
        //
        if(par_s->size() != 0)
        {
            int slotPosition[particlesNumMax * 2] = {0};
            int s = 0;
            // iter from begining to find the slots
            for (auto iter = par_s->begin(); iter != par_s->end(); ++iter)
            {
                if(iter->AliveParticle() == true)
                    continue;
                slotPosition[s] = iter - par_s->begin() + 1;
                s = s + 1;
            }
            // insert the slots
            for (int i = 0; i < particlesNumMax * 2; ++i)
            {
                //
                if (slotPosition[i] == 0)
                    break;
                // slot position: par_s->begin() + slotPosition[i] - 1
                auto iter = par_s->end() - 1 ;
                while (iter->AliveParticle() == false)
                {
                    --iter;
                }
                int checkPos = par_s->begin() + slotPosition[i] - 1 - iter;
                if(checkPos < 0)    
                {
                    *(par_s->begin() + slotPosition[i] - 1) = *iter;
                    iter->SetOutParticles();
                }
            }
            // erease: s totally
            for (int i = 0; i < s; ++i)
            {
                par_s->pop_back();
            }
        }
    }

    //***********************************************************************
    // Rearrange particles for species, removing slot of out-cell particles, and
    // randomize array of particles
    inline void SpeciesRandom(int a)
    {
        Particles tempP;
        vector<Particles> *par_s = NULL;
        if (a == 0)
        {
            par_s = particles_H;
        }
        else if (a == 1)
        {
            par_s = particles_He;
        }
        else if (a == 2)
        {
            par_s = particles_O;
        }
        //
        if (par_s->size() != 0)
        {
        // test   // iter from end to begin
        // test   vector<Particles>::iterator bw_step = par_s->end() - 1;
        // test   //  Move blank element to the end of vector
        // test   for (auto iter = par_s->begin(); iter != par_s->end(); ++iter)
        // test   {
        // test       if (iter == bw_step)
        // test           break;
        // test       //  Find the first Blank element
        // test       if (iter->AliveParticle())
        // test           continue;
        // test       //  Find the last NON-BLANK element
        // test       while (!bw_step->AliveParticle() && bw_step != par_s->begin())
        // test       {
        // test           bw_step = bw_step - 1;
        // test       }
        // test       //  swap element
        // test       if (!iter->AliveParticle() && iter < bw_step)
        // test       {
        // test           *iter = Particles(*bw_step);
        // test           bw_step->SetOutParticles();
        // test           bw_step = bw_step - 1;
        // test       }
        // test   }
        // test   // erase BLANK element at the end
        // test   int num = 0;
        // test   for (auto iter = par_s->rbegin(); iter != par_s->rend(); ++iter)
        // test   {
        // test       if (!iter->AliveParticle())
        // test           num = num + 1;
        // test       else
        // test           break;
        // test   }
        // test   for (int i = 0; i < num; ++i)
        // test   {
        // test       par_s->pop_back();
        // test   }

            //  old version of erase Blank element at the end
            //bw_step = par_s->end() - 1;
            //while (bw_step->AliveParticle() == 0)
            //{
            //    if(bw_step - par_s->begin() != 0) // takes huge time at this step
            //        bw_step = bw_step - 1;
            //    par_s->pop_back();
            //}
            // Randomize element
            for (auto iter = par_s->begin(); iter != par_s->end(); ++iter)
            {
                //    int fw_step = static_cast<int>(std::floor((par_s->end() - iter) * dRand()));
                int fw_step = floor((par_s->end() - iter) * dRand());
                //int fw_step = rand() % (par_s->end() - iter);
                if (fw_step != 0 && fw_step != par_s->end() - iter)
                {
                    //
                    //tempP = Particles(*iter);
                    //*iter = Particles(*(iter + fw_step));
                    //*(iter + fw_step) = Particles(tempP);
                    // using swap
                    std::swap(*iter, *(iter+fw_step));
                }
            }
        }
    }
    //  trans guiding center to absolute velocity
    inline void TransAbsoluteVelocity(  int a,
                                        vector<Vector3>& b3atvertex)
 {
        double Pi = 3.14159265359;
        vector<Particles> *ptrP;
        //int sizeP;
        //Particles *par[10];
        //Particles *tempP = NULL;
        //int numP = 0;
        //
        if (a == 0)
        {
            ptrP = particles_H;
        }
        else if (a == 1)
        {
            ptrP = particles_He;
        }
        else //if (a == 2)
        {
            ptrP = particles_O;
        }
        //
        double mi0;
        if( a == 0)
        {   
            mi0 = mi0_H;
        }else if( a == 1)
        {
            mi0 = mi0_He;
        } else //if( a == 2)
        {
            mi0 = mi0_O;
        }
        //
        for( auto iter = ptrP->begin(); iter != ptrP->end(); ++iter)
        {
            if(iter->AliveParticle() == false)
            {
                std::cout << " TransAbsoluteVelocity error \n";
                exit(1);
            }
            //
            struct structg tempStr;
            struct structPar tempStrPar;
            tempStr = iter->InttoStrp1();
            StructPar(tempStr, tempStrPar);
            // find the local b3 
            Vector3 b3cell;
            b3cell.Setx(b3atvertex[0].x() * tempStrPar.w000 + b3atvertex[1].x() * tempStrPar.w100 + b3atvertex[2].x() * tempStrPar.w110 + b3atvertex[3].x() * tempStrPar.w010 + b3atvertex[4].x() * tempStrPar.w001 + b3atvertex[5].x() * tempStrPar.w101 + b3atvertex[6].x() * tempStrPar.w111 + b3atvertex[7].x() * tempStrPar.w011);
            b3cell.Sety(b3atvertex[0].y() * tempStrPar.w000 + b3atvertex[1].y() * tempStrPar.w100 + b3atvertex[2].y() * tempStrPar.w110 + b3atvertex[3].y() * tempStrPar.w010 + b3atvertex[4].y() * tempStrPar.w001 + b3atvertex[5].y() * tempStrPar.w101 + b3atvertex[6].y() * tempStrPar.w111 + b3atvertex[7].y() * tempStrPar.w011);
            b3cell.Setz(b3atvertex[0].z() * tempStrPar.w000 + b3atvertex[1].z() * tempStrPar.w100 + b3atvertex[2].z() * tempStrPar.w110 + b3atvertex[3].z() * tempStrPar.w010 + b3atvertex[4].z() * tempStrPar.w001 + b3atvertex[5].z() * tempStrPar.w101 + b3atvertex[6].z() * tempStrPar.w111 + b3atvertex[7].z() * tempStrPar.w011);
            // find two eigen vector perpendicular to b3cell
            Vector3 velEigen1, velEigen2;
            if (b3cell.x() == 0.0 && b3cell.y() == 0.0)
            {
                velEigen1 = Vector3(1.0, 0.0, 0.0);
                velEigen2 = Vector3(0.0, 1.0, 0.0);
            }
            else
            {
                // set velEigen1: parallel to xOy plane
                velEigen1 = Vector3(1.0, -b3cell.x() / b3cell.y(), 0.0);
                velEigen1 = velEigen1.NormalizedVector();
                // Calculate velEigen2
                velEigen2 = velEigen1.CrossProduct(b3cell).NormalizedVector();
            }
            //
            double angle = dRand() * Pi * 2.0;
            double vel_mu = sqrt(iter->MagneticIvarient() * b3cell.norm() * 2.0 / mi0);
            Vector3 vel_perp = velEigen1.ScaleProduct(vel_mu * cos(angle)).PlusProduct(velEigen2.ScaleProduct(vel_mu * sin(angle)));
            Vector3 v_abs = vel_perp.PlusProduct(iter->VelParticles());
            //
            iter->SetVelocity(v_abs);
        }
    }
    
    //************************************************************************
    //************************************************************************
    // wave-particles interactions, lin (1992), varney (2016)
    // int a : species identifier
    // double P_0 : PSD base level (can be adjusted)
    // double timestep : timestep which is equal to tstep
    // vector<Vector3> b3atvertex : b3 at eight corner vertex of the cell
    inline void WaveParticlesInteraction(  int a, double timestep,
                                    vector<Vector3> &b3atvertex)
    {
        vector<Particles> *ptrP;
        //
        if (a == 0)
        {
            ptrP = particles_H;
        }
        else if (a == 1)
        {
            ptrP = particles_He;
        }
        else if (a == 2)
        {
            ptrP = particles_O;
        }
        //
        double mi0 = 0.0;
        if( a == 0)
        {   
            mi0 = mi0_H;
        }else if( a == 1)
        {
            mi0 = mi0_He;
        } else if( a == 2)
        {
            mi0 = mi0_O;
        }
        //
        for( auto iter = ptrP->begin(); iter != ptrP->end(); ++iter)
        {
            if(iter->AliveParticle() == false)
            {
                std::cout << " TransGuidingCenterVelocity error \n";
                exit(1);
            }
            //
            struct structg tempStr;
            struct structPar tempStrPar;
            tempStr = iter->InttoStrp1();
            StructPar(tempStr, tempStrPar);
            // find the local b3 
            Vector3 b3cell;
            b3cell.Setx(b3atvertex[0].x() * tempStrPar.w000 + b3atvertex[1].x() * tempStrPar.w100 + b3atvertex[2].x() * tempStrPar.w110 + b3atvertex[3].x() * tempStrPar.w010 + b3atvertex[4].x() * tempStrPar.w001 + b3atvertex[5].x() * tempStrPar.w101 + b3atvertex[6].x() * tempStrPar.w111 + b3atvertex[7].x() * tempStrPar.w011);
            b3cell.Sety(b3atvertex[0].y() * tempStrPar.w000 + b3atvertex[1].y() * tempStrPar.w100 + b3atvertex[2].y() * tempStrPar.w110 + b3atvertex[3].y() * tempStrPar.w010 + b3atvertex[4].y() * tempStrPar.w001 + b3atvertex[5].y() * tempStrPar.w101 + b3atvertex[6].y() * tempStrPar.w111 + b3atvertex[7].y() * tempStrPar.w011);
            b3cell.Setz(b3atvertex[0].z() * tempStrPar.w000 + b3atvertex[1].z() * tempStrPar.w100 + b3atvertex[2].z() * tempStrPar.w110 + b3atvertex[3].z() * tempStrPar.w010 + b3atvertex[4].z() * tempStrPar.w001 + b3atvertex[5].z() * tempStrPar.w101 + b3atvertex[6].z() * tempStrPar.w111 + b3atvertex[7].z() * tempStrPar.w011);
            // calculate W_perp using W_perp = mu * |b3|
            double W_perp = iter->MagneticIvarient() * b3cell.norm();
            // calculate cos(latitude) at position of particle
            double cosLatitude = sqrt( iter->PosParticles().x() * iter->PosParticles().x() + iter->PosParticles().y() * iter->PosParticles().y()) / iter->PosParticles().norm();
            // calcualte P_h
            double P_h = P_0 * pow(cosLatitude, 8.0) * exp( - W_perp * W_perp / 50.0 / 50.0);
            // adjust P_h for He and O because of different gyro-frequency, alphaPSD is a parameter =1.7 mentioned in varney(2016)
            if (a == 1)
            {
                P_h *= pow(4.0, alphaPSD);
            }
            else if (a == 2)
            {
                P_h *= pow(16.0, alphaPSD);
            }
            // calculate sigma_perp : sigma_perp = q^2 / 2m * P_H * tstep
            double sigma_perp = qi0 * qi0 / 2.0 / mi0 * P_h * timestep;
            // calculate dW_perp from a Gaussian probability distribution : P(dW) = (1 / 2Pi / sigma_perp) * exp( - dW^2 / 2sigma_perp^2) // should be sqrt(2pi)?
            double dW_perp = MaxwellDis(sigma_perp, 0.0);
            // calculate W_perp(new) 
            double cos_theta = cos(dRand() * 3.14159265358979323846);
            double W1_perp = W_perp + abs(dW_perp) + sqrt(W_perp) * sqrt(abs(dW_perp)) * 2.0 * cos_theta;
            // calculate new mu- invariant magnetic moment
            iter->SetMu( W1_perp / b3cell.norm());
        }
    }
    //
    //  trans guiding center to absolute velocity
    inline void TransGuidingCenterVelocity(    int a,
                                                vector<Vector3>& b3atvertex)
 {
        vector<Particles> *ptrP;
        //
        if (a == 0)
        {
            ptrP = particles_H;
        }
        else if (a == 1)
        {
            ptrP = particles_He;
        }
        else //if (a == 2)
        {
            ptrP = particles_O;
        }
        //
        double mi0;
        if( a == 0)
        {   
            mi0 = mi0_H;
        }else if( a == 1)
        {
            mi0 = mi0_He;
        } else //if( a == 2)
        {
            mi0 = mi0_O;
        }
        //
        for( auto iter = ptrP->begin(); iter != ptrP->end(); ++iter)
        {
            if(iter->AliveParticle() == false)
            {
                std::cout << " TransGuidingCenterVelocity error \n";
                exit(1);
            }
            //
            struct structg tempStr;
            struct structPar tempStrPar;
            tempStr = iter->InttoStrp1();
            StructPar(tempStr, tempStrPar);
            // find the local b3 
            Vector3 b3cell;
            b3cell.Setx(b3atvertex[0].x() * tempStrPar.w000 + b3atvertex[1].x() * tempStrPar.w100 + b3atvertex[2].x() * tempStrPar.w110 + b3atvertex[3].x() * tempStrPar.w010 + b3atvertex[4].x() * tempStrPar.w001 + b3atvertex[5].x() * tempStrPar.w101 + b3atvertex[6].x() * tempStrPar.w111 + b3atvertex[7].x() * tempStrPar.w011);
            b3cell.Sety(b3atvertex[0].y() * tempStrPar.w000 + b3atvertex[1].y() * tempStrPar.w100 + b3atvertex[2].y() * tempStrPar.w110 + b3atvertex[3].y() * tempStrPar.w010 + b3atvertex[4].y() * tempStrPar.w001 + b3atvertex[5].y() * tempStrPar.w101 + b3atvertex[6].y() * tempStrPar.w111 + b3atvertex[7].y() * tempStrPar.w011);
            b3cell.Setz(b3atvertex[0].z() * tempStrPar.w000 + b3atvertex[1].z() * tempStrPar.w100 + b3atvertex[2].z() * tempStrPar.w110 + b3atvertex[3].z() * tempStrPar.w010 + b3atvertex[4].z() * tempStrPar.w001 + b3atvertex[5].z() * tempStrPar.w101 + b3atvertex[6].z() * tempStrPar.w111 + b3atvertex[7].z() * tempStrPar.w011);
            //
            Vector3 b_unit = b3cell.NormalizedVector();
            Vector3 v_guidingCenter = b_unit.ScaleProduct(iter->VelParticles().DotProduct(b_unit));
            Vector3 vel_perp = iter->VelParticles().MinusProduct(v_guidingCenter);
            double mu = vel_perp.norm2() * 0.5 * mi0 / b3cell.norm();
            //
            iter->UpdateTempPar(v_guidingCenter, mu); 
        }
    }
    //***********************************************************************
    //  Need the size of each particles array be in range of (Nmin, NMax),
    //  if the size is small, it needs split, the biggest particle split into a few more
    //  if the size is big, it need combine, the smallest a few particles combine into one
    inline void SplitOrCombine(int a)
    {
        vector<Particles> *ptrP;
        int sizeP;
        Particles *par[10];
        Particles *tempP = NULL;
        int numP = 0;
        //
        if (a == 0)
        {
            ptrP = particles_H;
        }
        else if (a == 1)
        {
            ptrP = particles_He;
        }
        else if (a == 2)
        {
            ptrP = particles_O;
        }
        //
        sizeP = ptrP->size();
        // 1: combine situation
        while (sizeP > particlesNumMax)
        {
            tempP = &(*ptrP)[0];
            // initialize par[10]- the first 10 element and tempP- max in par[10]
            for (int ii = 0; ii < 10; ii++)
            {
                par[ii] = &(*ptrP)[ii];
                if (par[ii]->WeightNi() > tempP->WeightNi())
                {
                    tempP = par[ii];
                    numP = ii;
                }
            }
            // find out the min 10 elements in ptrP
            for (int i = 10; i < sizeP; i++)
            {
                if ((*ptrP)[i].WeightNi() < tempP->WeightNi())
                {
                    par[numP] = &(*ptrP)[i];
                    // re-identify tempP- max in par[10]
                    tempP = par[0];
                    for (int ii = 1; ii < 10; ii++)
                    {
                        if (par[ii]->WeightNi() > tempP->WeightNi())
                        {
                            tempP = par[ii];
                            numP = ii;
                        }
                    }
                }
            }
            //  combine 10 into 1
            double pos[3] = {0.0, 0.0, 0.0};
            double vel[3] = {0.0, 0.0, 0.0};
            double weight = 0.0;
            double mu = 0.0;
            for (int i = 0; i < 10; i++)
            {
                weight += par[i]->WeightNi();
                pos[0] += par[i]->WeightNi() * par[i]->PosParticles().x();
                pos[1] += par[i]->WeightNi() * par[i]->PosParticles().y();
                pos[2] += par[i]->WeightNi() * par[i]->PosParticles().z();
                vel[0] += par[i]->WeightNi() * par[i]->VelParticles().x();
                vel[1] += par[i]->WeightNi() * par[i]->VelParticles().y();
                vel[2] += par[i]->WeightNi() * par[i]->VelParticles().z();
                mu += par[i]->WeightNi() * par[i]->MagneticIvarient();
                par[i]->SetOutParticles();
            }
            Vector3 posV3 = Vector3(pos[0] / weight, pos[1] / weight, pos[2] / weight);
            Vector3 velV3 = Vector3(vel[0] / weight, vel[1] / weight, vel[2] / weight);
            mu = mu / weight;
            Particles newP = Particles(static_cast<uint_64>(0), posV3, velV3, weight, mu);
            newP.CalculateUint_64();
            ptrP->push_back(newP);
            //  check empty slot in array
            SpeciesRearrange(a);
            sizeP = ptrP->size();
            //
        }
        //  split
        int size_small = 0;
        while (sizeP < particlesNumMin && sizeP != 0)
        {
            //  In case sizeP is smaller than 10
            if (sizeP < 10)
                size_small = sizeP;
            else
            {
                size_small = 10;
            }
            //  tempP is the smallest
            tempP = &(*ptrP)[0];
            // initialize par[10]- the first 10 element and tempP- min in par[10]
            for (int ii = 0; ii < size_small; ii++)
            {
                par[ii] = &(*ptrP)[ii];
                if (par[ii]->WeightNi() < tempP->WeightNi())
                {
                    tempP = par[ii];
                    numP = ii;
                }
            }
            // find out the max 10 elements in ptrP
            for (int i = size_small; i < sizeP; i++)
            {
                if ((*ptrP)[i].WeightNi() > tempP->WeightNi())
                {
                    par[numP] = &(*ptrP)[i];
                    // re-identify tempP- min in par[10]
                    tempP = par[0];
                    for (int ii = 1; ii < size_small; ii++)
                    {
                        if (par[ii]->WeightNi() < tempP->WeightNi())
                        {
                            tempP = par[ii];
                            numP = ii;
                        }
                    }
                }
            }
            // split: 1 into 2
            // posP, mu remain same; vp, weight change
            double ran;
            Vector3 pos_temp;
            Vector3 vel_temp;
            double weight_temp;
            double mu_temp;
            for (int i = 0; i < size_small; i++)
            {
                ran = dRand();
                pos_temp = (*ptrP)[i].PosParticles();
                vel_temp = (*ptrP)[i].VelParticles();
                weight_temp = (*ptrP)[i].WeightNi();
                (*ptrP)[i].SetVelocity(vel_temp.ScaleProduct(1.0 + 0.05 * ran));
                (*ptrP)[i].SetWeightNi(weight_temp / 2.0);
                mu_temp = (*ptrP)[i].MagneticIvarient();
                Particles newP = Particles(static_cast<uint_64>(0), pos_temp, vel_temp.ScaleProduct(1.0 - 0.05 * ran), weight_temp / 2.0, mu_temp);
                newP.CalculateUint_64();
                ptrP->push_back(newP);
            }
            sizeP = ptrP->size();
        }
        // random particles in heap
        SpeciesRandom(a);
    }
    //***********************************************************************
    // Calculate post-collision particles for two species with different
    // weight of different particles.
    inline void mutual_coll(int a,
                            int b,
                            double lnA)
    {
        //double pi2 = 6.283185307f; // 2pi
        //
        if (a == b)
            std::cout << " mutual collision different species " << std::endl;
        //
        vector<Particles> *par_a = NULL;
        if (a == 0)
        {
            par_a = particles_H;
        }
        else if (a == 1)
        {
            par_a = particles_He;
        }
        else //if (a == 2)
        {
            par_a = particles_O;
        }
        vector<Particles> *par_b = NULL;
        if (b == 0)
        {
            par_b = particles_H;
        }
        else if (b == 1)
        {
            par_b = particles_He;
        }
        else //if (b == 2)
        {
            par_b = particles_O;
        }
        // a is longer array
        int inta = par_a->size();
        int intb = par_b->size();
        vector<Particles> *par_aa = NULL, *par_bb = NULL;
        //
        if (inta >= intb)
        {
            par_aa = par_a;
            par_bb = par_b;
        }
        else
        {
            par_aa = par_b;
            par_bb = par_a;
        }
        inta = par_aa->size();
        intb = par_bb->size();
        if (inta == 0 || intb == 0)
            return;
        //
        for (int i = 0; i < inta; i++)
            par_random[i] = dRand();
        // ddt
        double na = 0.0;
        double nb = 0.0;
        double nab = 0.0;
        int leng=0;
        double wa = 0.0;
        double wb = 0.0;
        double ddt = 0.0;
        for (auto iter = par_aa->begin(); iter < par_aa->end(); ++iter)
        {
            na += iter->WeightNi();
            leng = iter - par_aa->begin();
            wa = iter->WeightNi();
            wb = (*par_bb)[floor(intb * par_random[leng])].WeightNi();
            // nab += wa * wb / max(wa, wb);
            nab += wa  / max(wa, wb) * wb;
        }
        for (auto iter = par_bb->begin(); iter < par_bb->end(); ++iter)
        {
            nb += iter->WeightNi();
        }
        nb = nb / volume;
        ddt = tstep * na / nab;
        //
        for (auto iter = par_aa->begin(); iter != par_aa->end(); ++iter)
        {
            // particle1 = *iter, particles2 = par_bb[(int)(intb * ptrPar_random[leng])]
            coll_vel_chng(a, b, ddt, nb, lnA,
                          *iter,
                          (*par_bb)[floor(intb * par_random[leng])]);
        }
    }
    //***********************************************************************
    // Collisions between ions of the same species
    inline void self_coll(int a,
                          double lnA)
    {
        //double pi2 = 6.283185307f;
        //
        vector<Particles> *par_a = NULL;
        if (a == 0)
        {
            par_a = particles_H;
        }
        else if (a == 1)
        {
            par_a = particles_He;
        }
        else if (a == 2)
        {
            par_a = particles_O;
        }
        //
        int inta = par_a->size();
        if (inta == 0)
            return;
        int inthf = inta / 2;
        // ddt
        double na = 0.0;
        double nab = 0.0;
        double wa = 0.0;
        double wb = 0.0;
        double ddt = 0.0;
        for (auto iter = par_a->begin(); iter != par_a->end(); ++iter)
        {
            na += iter->WeightNi();
        }
        // iter the first half particles
        for (auto iter = par_a->begin(); iter < par_a->begin() + inthf; ++iter)
        {
            wa = iter->WeightNi();
            wb = (iter + inthf)->WeightNi();
            //nab += wa * wb / max(wa, wb);
            nab += wa  / max(wa, wb) * wb;
        }
        if (inta % 2 != 0)
        {
            wa = par_a->begin()->WeightNi();
            wb = (par_a->end() - 1)->WeightNi();
            nab += wa  / max(wa, wb) * wb;
        }
        nab = nab * 2.0;
        ddt = tstep * na / nab;
        //
        na = na / volume; // number density of species a
        //
        for (auto iter = par_a->begin(); iter != par_a->begin() + inthf; ++iter)
        {
            coll_vel_chng(a, a, ddt, na, lnA,
                          *iter,
                          *(iter + inthf));
        }
        if (inta % 2 != 0 && inta != 1)
        {
            coll_vel_chng(a, a, ddt, na, lnA,
                          *(par_a->end() - 1),
                          *(par_a->begin()));
        }
    }

    //***********************************************************************
    // Calculate post-collision particles for two species with different
    // weight of different particles.
    inline void mutual_coll(int a,
                            int b,
                            double lnA,
                          GridsPoints ***** ptrArrayGrids)
    {
        //double pi2 = 6.283185307f; // 2pi
        //
        if (a == b)
            std::cout << " mutual collision different species " << std::endl;
        //
        vector<Particles> *par_a = NULL;
        if (a == 0)
        {
            par_a = particles_H;
        }
        else if (a == 1)
        {
            par_a = particles_He;
        }
        else if (a == 2)
        {
            par_a = particles_O;
        }
        vector<Particles> *par_b = NULL;
        if (b == 0)
        {
            par_b = particles_H;
        }
        else if (b == 1)
        {
            par_b = particles_He;
        }
        else if (b == 2)
        {
            par_b = particles_O;
        }
        // a is longer array
        int inta = par_a->size();
        int intb = par_b->size();
        vector<Particles> *par_aa = NULL, *par_bb = NULL;
        //
        if (inta >= intb)
        {
            par_aa = par_a;
            par_bb = par_b;
        }
        else
        {
            par_aa = par_b;
            par_bb = par_a;
        }
        inta = par_aa->size();
        intb = par_bb->size();
        if (inta == 0 || intb == 0)
            return;
        //
        for (int i = 0; i < inta; i++)
            par_random[i] = dRand();
        // ddt
        double na = 0.0;
        double nb = 0.0;
        double nab = 0.0;
        int leng;
        double wa = 0.0;
        double wb = 0.0;
        double ddt = 0.0;
        for (auto iter = par_aa->begin(); iter < par_aa->end(); ++iter)
        {
            na += iter->WeightNi();
            leng = iter - par_aa->begin();
            wa = iter->WeightNi();
            wb = (*par_bb)[floor(intb * par_random[leng])].WeightNi();
            // nab += wa * wb / max(wa, wb);
            nab += wa  / max(wa, wb) * wb;
        }
        for (auto iter = par_bb->begin(); iter < par_bb->end(); ++iter)
        {
            nb += iter->WeightNi();
        }
        nb = nb / volume;
        ddt = tstep * na / nab;
        //
        for (auto iter = par_aa->begin(); iter != par_aa->end(); ++iter)
        {
            // particle1 = *iter, particles2 = par_bb[(int)(intb * par_random[leng])]
            coll_vel_chng(a, b, ddt, nb, lnA,
                          ptrArrayGrids,
                          *iter,
                          (*par_bb)[floor(intb * par_random[leng])]);
        }
    }
    // Collisions between ions of the same species
    inline void self_coll(int a,
                          double lnA,
                          GridsPoints ***** ptrArrayGrids)
    {
        //double pi2 = 6.283185307f;
        //
        vector<Particles> *par_a = NULL;
        if (a == 0)
        {
            par_a = particles_H;
        }
        else if (a == 1)
        {
            par_a = particles_He;
        }
        else if (a == 2)
        {
            par_a = particles_O;
        }
        //
        int inta = par_a->size();
        if (inta == 0)
            return;
        int inthf = inta / 2;
        // ddt
        double na = 0.0;
        double nab = 0.0;
        double wa = 0.0;
        double wb = 0.0;
        double ddt = 0.0;
        for (auto iter = par_a->begin(); iter != par_a->end(); ++iter)
        {
            na += iter->WeightNi();
        }
        // iter the first half particles
        for (auto iter = par_a->begin(); iter < par_a->begin() + inthf; ++iter)
        {
            wa = iter->WeightNi();
            wb = (iter + inthf)->WeightNi();
            //nab += wa * wb / max(wa, wb);
            nab += wa  / max(wa, wb) * wb;
        }
        if (inta % 2 != 0)
        {
            wa = par_a->begin()->WeightNi();
            wb = (par_a->end() - 1)->WeightNi();
            nab += wa  / max(wa, wb) * wb;
        }
        nab = nab * 2.0;
        ddt = tstep * na / nab;
        //
        na = na / volume; // number density of species a
        //
        for (auto iter = par_a->begin(); iter != par_a->begin() + inthf; ++iter)
        {
            coll_vel_chng(a, a, ddt, na, lnA,
                          ptrArrayGrids,
                          *iter,
                          *(iter + inthf));
        }
        if (inta % 2 != 0 && inta != 1)
        {
            coll_vel_chng(a, a, ddt, na, lnA,
                          ptrArrayGrids,
                          *(par_a->end() - 1),
                          *(par_a->begin()));
        }
    }
    //***********************************************************************
    // Calculate Temperature of species a and return it
    // Assume that temperature is isotropic
    inline double TempColl(int a, Vector3 vel)
    {
        //
        //double Re = 6371200.0;
        double temp = 0.0;
        //double Re2 = Re * Re;
        //double Re3 = Re * Re * Re;
        //vector<Particles> *par_a = NULL;
        //double mkb[3] = { 1.2123442e-4, 4.813845268e-4,1.9242837e-3};
        //if (a == 0)
        //{
        //    par_a = particles_H;
        //}
        //else if (a == 1)
        //{
        //    par_a = particles_He;
        //}
        //else if (a == 2)
        //{
        //    par_a = particles_O;
        //}
        ////
        //int den = par_a->size();
        ////
        //for (auto iter = par_a->begin(); iter != par_a->end(); ++iter)
        //{
        //    temp += iter->VelParticles().MinusProduct(vel).norm2();
        //}
        ////
        //if (den > 0)
        //    temp = temp / den * mkb[a];
        //
        return temp;
    }

    //***********************************************************************
    // Velocity distribution calculation
    // Calculate particles weight of species at velocity range[-60km/s, 60km/s]
    // 121 intervals
    // for array[vel_pr][vel_pp]
    inline void velDistArray(int a, GridsPoints *****ptrArray)
    {
        vector<Particles> *par_a = NULL;
        vector<vector<double> > *velDistPtr = NULL;
        double mi0;
        //double mkb[3] = {1.2123442e-4F, 4.813845268e-4F, 1.9242837e-3F};
        //double pi = acos(-1.0);

        //
        Vector3 b3Par, b_unit;
        double vel_mu; // vel_perp
        double vel_para;
        Vector3 v_para;
        if (a == 0)
        {
            par_a = particles_H;
            mi0 = mi0_H;
            velDistPtr = velDist_H;
        }
        else if (a == 1)
        {
            par_a = particles_He;
            mi0 = mi0_He;
            velDistPtr = velDist_He;
        }
        else //if (a == 2)
        {
            par_a = particles_O;
            mi0 = mi0_O;
            velDistPtr = velDist_O;
        }
        // index in array[array_a][array_b]
        int array_a, array_b;
        //double tempArray;
        //
        double totalWeight = 0.0;
        double cul_vel_pr = 0.0; //, cul_vel_pp = 0.0;
        double vel2_pr = 0.0, vel2_pp = 0.0;
        b3Par = b3Cell;
        b_unit = b3Par.NormalizedVector();
        //
        for (auto iter = par_a->begin(); iter != par_a->end(); ++iter)
        {
            //weight
            totalWeight += iter->WeightNi();
            //
            //b3Par = iter->B3atParticles(ptrArray);
            //b_unit = b3Par.NormalizedVector();
            //
            // para velocity from guiding center motion
            vel_para = iter->VelParticles().DotProduct(b_unit);
            cul_vel_pr += vel_para * iter->WeightNi();
            vel2_pr += vel_para * vel_para * iter->WeightNi();
            //
            if (vel_para >= -velParaMax && vel_para <= velParaMax)
            {
                array_a = std::floor((vel_para + velParaMax) / 2.0 / velParaMax * velDistRange_para);
            }
            else if (vel_para < -velParaMax)
            {
                array_a = 0;
            }
            else //if (vel_para > velParaMax)
            {
                array_a = velDistRange_para - 1;
            }
            //  perp velocity from mu;
            vel_mu = sqrt(iter->MagneticIvarient() * b3Par.norm() * 2.0 / mi0);
            //cul_vel_pp = 0.0; // mu direction should be zero
            vel2_pp += vel_mu * vel_mu * iter->WeightNi();
            // velocity distribution for only one direction in perpendicular plane
            vel_mu = vel_mu * cos(pi * dRand());
            if (vel_mu >= -velPerpMax && vel_mu <= velPerpMax)
            {
                array_b = std::floor((vel_mu + velPerpMax) / 2.0 / velPerpMax * velDistRange_mu);
            }
            else if (vel_mu < -velPerpMax)
            {
                array_b = 0;
            }
            else //if (vel_mu > velPerpMax)
            {
                array_b = velDistRange_mu - 1;
            }
            //
            // culmulate weight of VelArray[][]
            //std::cout << "array index " << array_a << " " << array_b << "\n";
            if (array_a >= 0 && array_a < velDistRange_para && array_b >= 0 && array_b < velDistRange_mu)
                (*velDistPtr)[array_a][array_b] += iter->WeightNi();
            //
        }

        //for (auto iter = par_a->begin(); iter != par_a->end(); ++iter)
        //{
        //    b3Par = iter->B3atParticles(ptrArray);
        //    b_unit = b3Par.NormalizedVector();
        //    // para velocity from guiding center motion
        //    vel_para = iter->VelParticles().DotProduct(b_unit);
        //    if (vel_para > -velParaMax && vel_para < velParaMax)
        //    {
        //        array_a = static_cast<int>(std::floor((vel_para + velParaMax) / 2.0 / velParaMax * velDistRange_para));
        //    }
        //    else if (vel_para < -velParaMax)
        //    {
        //        array_a = 0;
        //    }
        //    else if (vel_para > velParaMax)
        //    {
        //        array_a = velDistRange_para - 1;
        //    }
        //    //  perp velocity from mu;
        //    vel_mu = sqrt(iter->MagneticIvarient() * b3Par.norm() * 2.0 / mi0);
        //    array_b = (int)(vel_mu / velPerpMax * velDistRange_mu);
        //    if (vel_mu < velPerpMax)
        //        array_b = static_cast<int>(std::floor(vel_mu / velPerpMax * velDistRange_mu));
        //    else
        //    {
        //        array_b = velDistRange_mu - 1;
        //    }
        //    //
        //    (*velDistPtr)[array_a][array_b] += iter->WeightNi();
        //}
        //
    }
    // //***********************************************************************
    // // Velocity distribution calculation
    // // Calculate particles weight of species at velocity range[-50km/s, 50km/s]
    // // 100 intervals
    // // For two array[], only for parallel or perpendicular
    // inline void VelDistCalculation(int a, GridsPoints *****ptrArray)
    // {
    //     vector<Particles> *par_a = NULL;
    //     vector<double> *velDist_para = NULL;
    //     vector<double> *velDist_mu = NULL;
    //     double mi0;
    //     //
    //     Vector3 b3Par, b_unit;
    //     double vel_mu; // vel_perp
    //     double vel_para;
    //     Vector3 v_para;
    //     if (a == 0)
    //     {
    //         par_a = particles_H;
    //         mi0 = mi0_H;
    //     }
    //     else if (a == 1)
    //     {
    //         par_a = particles_He;
    //         mi0 = mi0_He;
    //     }
    //     else if (a == 2)
    //     {
    //         par_a = particles_O;
    //         mi0 = mi0_O;
    //     }
    //     //
    //     for (auto iter = par_a->begin(); iter != par_a->end(); ++iter)
    //     {
    //         b3Par = iter->B3atParticles(ptrArray);
    //         b_unit = b3Par.NormalizedVector();
    //         //  perp velocity from mu;
    //         vel_mu = sqrt(iter->MagneticIvarient() * b3Par.norm() * 2.0 / mi0);
    //         if (vel_mu < velPerpMax)
    //             (*velDist_mu)[(int)(vel_mu / velPerpMax * velDistRange_mu)] += iter->WeightNi();
    //         else
    //         {
    //             (*velDist_mu)[velDistRange_mu - 1] += iter->WeightNi();
    //         }
    //         // para velocity from guiding center motion
    //         vel_para = iter->VelParticles().DotProduct(b_unit);
    //         if (vel_para > -velParaMax && vel_para < velParaMax)
    //             (*velDist_para)[(int)((vel_para + velParaMax) / 2.0 / velParaMax * velDistRange_para)];
    //         else if (vel_para < -velParaMax)
    //         {
    //             (*velDist_para)[0] += iter->WeightNi();
    //         }
    //         else if (vel_para > velParaMax)
    //         {
    //             (*velDist_para)[velDistRange_para - 1] += iter->WeightNi();
    //         }
    //     }
    //     //
    // }
    //****************************************************************
    // calculate average of particles for print
    // a - species
    inline void AverVelDisArray(int average_steps, int a)
    {
        vector<vector<double>> *par_a = NULL;
        if (a == 0)
        {
            par_a = velDist_H;
        }
        else if (a == 1)
        {
            par_a = velDist_He;
        }
        else if (a == 2)
        {
            par_a = velDist_O;
        }
        for (int i = 0; i < velDistRange_para; i++)
        {
            for (int j = 0; j < velDistRange_mu; j++)
            {
                (*par_a)[i][j] = (*par_a)[i][j] / average_steps / volume / velSpace;
            }
        }
    }
    // Reset velDist
    inline void ResetVelDistArray(int mode)
    {
        for (int i = 0; i < velDistRange_para; i++)
        {
            for (int j = 0; j < velDistRange_mu; j++)
            {
                (*velDist_H)[i][j] = 0.0;
            }
        }
        for (int i = 0; i < velDistRange_para; i++)
        {
            for (int j = 0; j < velDistRange_mu; j++)
            {
                (*velDist_He)[i][j] = 0.0;
            }
        }
        for (int i = 0; i < velDistRange_para; i++)
        {
            for (int j = 0; j < velDistRange_mu; j++)
            {
                (*velDist_O)[i][j] = 0.0;
            }
        }
    }
    // Constructors
    GridsCells();
    GridsCells(double gradPe_x, double gradPe_y, double gradPe_z,
               double gradPo_x, double gradPo_y, double gradPo_z,
               double b3_x, double b3_y, double b3_z,
               vector<Particles> part_H,
               vector<Particles> part_He,
               vector<Particles> part_O,
               vector<vector<double>> velDist_HH,
               vector<vector<double>> velDist_HEHE,
               vector<vector<double>> velDist_OO,
               double* par_r,
               double vol);
    //
private:
    Vector3 gradPe;
    Vector3 gradPotential;
    Vector3 b3Cell; // maybe not used
    //
    vector<Particles> *particles_H;
    vector<Particles> *particles_He;
    vector<Particles> *particles_O;
    //
    vector<vector<double>> *velDist_H;
    vector<vector<double>> *velDist_He;
    vector<vector<double>> *velDist_O;
    //
    double *par_random; // used in mutual-collision
    double volume;
    //// temperature
    //double temp_coll_H_pr;
    //double temp_coll_H_pp;
    //double temp_coll_He_pr;
    //double temp_coll_He_pp;
    //double temp_coll_O_pr;
    //double temp_coll_O_pp;
    //// <|v_pp|> and <v_pr>
    //double v_pp_H;
    //double v_pr_H;
    //double v_pp_He;
    //double v_pr_He;
    //double v_pp_O;
    //double v_pr_O;
};
#endif