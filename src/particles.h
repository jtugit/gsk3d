#ifndef _PARTICLES_H_
#define _PARTICLES_H_
/// Particle information
#include <iostream>
#include "mathutil.h"
#include "parameters.h"
#include "vector3.h"
#include "structdef.h"
#include "fieldsgrids.h"
#include <bitset>
#include "module_base.h"
#include <vector>

class Particles
{
public:
    friend class GridsPoints;
    friend class Vector3;

    //************************************************************************
    //************************************************************************
    // Transfore uint_64 posuInt to F I J K of the grids, from which we can latey determin the total
    // 8 grids-points of the cell in which the particles is.

    //************************************************************************
    //************************************************************************
    inline structg InttoStrp1()
    {
        struct structg strg = {0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
        strg.face = (int)(posUint >> 61);
        //std::cout << posUint << " "<< strg.face <<" "<< strg.ig<<" "<< strg.jg<<" "<< strg.kg << std::endl;
        for (int i = 0; i < fieldsGridsLevel; i++)
        {
            strg.ig = (strg.ig << 1) + ((int)(posUint >> (60 - i * 3)) & 1);
            strg.jg = (strg.jg << 1) + ((int)(posUint >> (59 - i * 3)) & 1);
            strg.kg = (strg.kg << 1) + ((int)(posUint >> (58 - i * 3)) & 1);
        }
        // adjust to domain 1/2
        strg.kg += (int)(posUint & 1 << fieldsGridsLevel);
        //std::cout << posUint << " "<< strg.face <<" "<< strg.ig<<" "<< strg.jg<<" "<< strg.kg << std::endl;

        for (int i = fieldsGridsLevel; i < particlesGridsLevel; i++)
        {
            strg.iw = (strg.iw << 1) + ((int)(posUint >> (60 - i * 3)) & 1);
            strg.jw = (strg.jw << 1) + ((int)(posUint >> (59 - i * 3)) & 1);
            strg.kw = (strg.kw << 1) + ((int)(posUint >> (58 - i * 3)) & 1);
        }
        strg.vx = vp.x();
        strg.vy = vp.y();
        strg.vz = vp.z();
        //    strg.mass = mass1;

        /*    std::cout << std::bitset<64>(posUint) << " "<< strg.face <<" "<< strg.ig<<" "<< strg.jg<<" "<< strg.kg << " "
              << strg.iw << " " << strg.jw << " " << strg.kw << std::endl;
    int pause ;
    std::cin >> pause;
*/
        return strg;
    }
    //************************************************************************
    //************************************************************************
    // Transfor uint_64 posuInt to i j k of the particle in the cell, which can help to determin the
    // weighting of particle density and velocity on the grids-points of the cell in which the particles
    // is.
    //************************************************************************
    //************************************************************************
    inline structg InttoStrp2()
    {
        struct structg strg;
        strg.face = (int)(posUint >> 61);
        for (int i = 0; i < fieldsGridsSize; i++)
        {
            strg.ig = (strg.ig << 1) + ((int)(posUint >> (60 - i * 3)) & 1);
            strg.jg = (strg.jg << 1) + ((int)(posUint >> (60 - 1 - i * 3)) & 1);
            strg.kg = (strg.kg << 1) + ((int)(posUint >> (60 - 2 - i * 3)) & 1);
        }
        for (int i = fieldsGridsSize + 1; i < particlesGridsSize; i++)
        {
            strg.iw = (strg.iw << 1) + ((int)(posUint >> (60 - i * 3) & 1));
            strg.jw = (strg.jw << 1) + ((int)(posUint >> (60 - 1 - i * 3) & 1));
            strg.kw = (strg.kw << 1) + ((int)(posUint >> (60 - 2 - i * 3) & 1));
        }
        strg.vx = vp.x();
        strg.vy = vp.y();
        strg.vz = vp.z();
        //    strg.mass = mass2;
        return strg;
    }

    //************************************************************************
    //************************************************************************
    // Calculate Latitude of magnetic field lines in which the particles are along
    //************************************************************************
    //************************************************************************
    inline int LatMagParticles()
    {
        double x = posP.x();
        double y = posP.y();
        double z = posP.z();
        //double PI = 3.1415926535897;

        double r_top = sqrt(x * x + y * y + z * z); // r of particles
        double theta_top = acos(z / r_top);         // theta of particles

        double theta_earth = asin(sin(theta_top) * sqrt(radius / r_top)); // theta of mapping point on earth
                                                                          //    std::cout << theta_earth * 180.0 / PI << "  ";
        if (theta_earth < PI / 2.0 - c0_latitude * PI / 180.0)            // in the open field area
            return 1;
        else
        {
            return 0;
        }
    }

    // // calculate wave-particle interaction, and renew mu
    // // Get interpolated PSD from 8 vertex
    // inline void WaveParticle(struct structg *strg_in, GridsPoints *****ptrArray_in, double mi0_in)
    // {
    //     // sigma
    //     double w000 = 1.0 * (cellSize1 - strg_in->iw - 0.5) * (cellSize1 - strg_in->jw - 0.5) * (cellSize1 - strg_in->kw - 0.5);
    //     double w100 = 1.0 * (0.5 + strg_in->iw) * (cellSize1 - strg_in->jw - 0.5) * (cellSize1 - strg_in->kw - 0.5);
    //     double w010 = 1.0 * (cellSize1 - strg_in->iw - 0.5) * (1 + strg_in->jw - 0.5) * (cellSize1 - strg_in->kw - 0.5);
    //     double w110 = 1.0 * (0.5 + strg_in->iw) * (0.5 + strg_in->jw) * (cellSize1 - strg_in->kw - 0.5);
    //     double w001 = 1.0 * (cellSize1 - strg_in->iw - 0.5) * (cellSize1 - strg_in->jw - 0.5) * (0.5 + strg_in->kw);
    //     double w101 = 1.0 * (0.5 + strg_in->iw) * (cellSize1 - strg_in->jw - 0.5) * (0.5 + strg_in->kw);
    //     double w011 = 1.0 * (cellSize1 - strg_in->iw - 0.5) * (0.5 + strg_in->jw) * (0.5 + strg_in->kw);
    //     double w111 = 1.0 * (0.5 + strg_in->iw) * (0.5 + strg_in->jw) * (0.5 + strg_in->kw);
    //     double totalWeight = w000 + w100 + w010 + w110 + w001 + w101 + w011 + w111;
    //     w000 = w000 / totalWeight;
    //     w100 = w100 / totalWeight;
    //     w010 = w010 / totalWeight;
    //     w110 = w110 / totalWeight;
    //     w001 = w001 / totalWeight;
    //     w101 = w101 / totalWeight;
    //     w011 = w011 / totalWeight;
    //     w111 = w111 / totalWeight;

    //     double temp1 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 1][strg_in->kg]->PSD_H();
    //     double temp2 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 1][strg_in->kg]->PSD_H();
    //     double temp3 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 2][strg_in->kg]->PSD_H();
    //     double temp4 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 2][strg_in->kg]->PSD_H();
    //     double temp5 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 1][strg_in->kg + 1]->PSD_H();
    //     double temp6 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 1][strg_in->kg + 1]->PSD_H();
    //     double temp7 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 2][strg_in->kg + 1]->PSD_H();
    //     double temp8 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 2][strg_in->kg + 1]->PSD_H();

    //     double psd = temp1 * w000 + temp2 * w100 + temp3 * w110 + temp4 * w010 + temp5 * w001 + temp6 * w101 + temp7 * w111 + temp8 * w011;

    //     Vector3 tempb1 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 1][strg_in->kg]->B3();
    //     Vector3 tempb2 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 1][strg_in->kg]->B3();
    //     Vector3 tempb3 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 2][strg_in->kg]->B3();
    //     Vector3 tempb4 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 2][strg_in->kg]->B3();
    //     Vector3 tempb5 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 1][strg_in->kg + 1]->B3();
    //     Vector3 tempb6 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 1][strg_in->kg + 1]->B3();
    //     Vector3 tempb7 = ptrArray_in[strg_in->face][strg_in->ig + 2][strg_in->jg + 2][strg_in->kg + 1]->B3();
    //     Vector3 tempb8 = ptrArray_in[strg_in->face][strg_in->ig + 1][strg_in->jg + 2][strg_in->kg + 1]->B3();

    //     Vector3 tempb;
    //     tempb.Setx(tempb1.x() * w000 + tempb2.x() * w100 + tempb3.x() * w110 + tempb4.x() * w010 + tempb5.x() * w001 + tempb6.x() * w101 + tempb7.x() * w111 + tempb8.x() * w011);

    //     tempb.Sety(tempb1.y() * w000 + tempb2.y() * w100 + tempb3.y() * w110 + tempb4.y() * w010 + tempb5.y() * w001 + tempb6.y() * w101 + tempb7.y() * w111 + tempb8.y() * w011);

    //     tempb.Setz(tempb1.z() * w000 + tempb2.z() * w100 + tempb3.z() * w110 + tempb4.z() * w010 + tempb5.z() * w001 + tempb6.z() * w101 + tempb7.z() * w111 + tempb8.z() * w011);

    //     if (mi0_in == mi0_He)
    //     {
    //         psd = psd * pow(4.0, alphaPSD);
    //     }
    //     else if (mi0_in == mi0_O)
    //     {
    //         psd = psd * pow(16.0, alphaPSD);
    //     }

    //     double sigma = qi0 * qi0 / 2.0 / mi0_H * psd * tstep / div_max;
    //     if (mi0_in == mi0_He)
    //     {
    //         sigma = sigma / 4.0;
    //     }
    //     else if (mi0_in == mi0_O)
    //     {
    //         sigma = sigma / 16.0;
    //     }

    //     // calculate dW_perp
    //     double temp = dRand();
    //     double dW_perp = erfinv(2.0 * temp - 1.0) * sqrt(2.0) * sigma;
    //     if (dW_perp < 0.0)
    //         dW_perp = -1.0 * dW_perp;

    //     // update mu
    //     double W_perp = mu * tempb.norm();
    //     double temprand = (double)rand() / (RAND_MAX);
    //     double Wnew = W_perp + dW_perp + 2.0 * sqrt(W_perp) * sqrt(dW_perp) * cos(2.0 * 3.14159265358979323846 * temprand);
    //     mu = Wnew / tempb.norm();
    // }


// YH 4/1 
// return potential at postition of unreal particle, method is similar to the beginning part
// of BorisMethod() below
    double PotentialAtUnrealPar(struct structPar &strg_in, GridsPoints *****ptrArray_in);
    //************************************************************************
    //************************************************************************
    // Calculate the new vp by Boris' Method, update vector3 vp
    // And return a int "0" means in the main domain
    // "1" means out of the main domain
    //
    //************************************************************************
    //************************************************************************
    int BorisMethod(struct structPar &strg_in, GridsPoints *****ptrArray_in, double mi0_in, int maindomain);

    // int BorisMethod(struct structPar &strg_in, std::vector<GridsPoints> gridsArray, int type_a, int maindomain);
    
    //************************************************************************
    //************************************************************************
    // Calculate the new position of the particles in uint_64
    // And return a int "0" means in the main domain
    // "1" means out of the main domain
    //
    int UpdateUint_64();
    // int UpdateUint_64_temp();
    
    int UpdateUint_64_test();
    // A new position and calculate the new Uint_64
    inline void CalculateUint_64()
    {
        uint_64 face = 0, ip = 0, jp = 0, kp = 0;
        double px, py, pz;
        double temp[2];
        px = posP.x();
        py = posP.y();
        pz = posP.z();
        double L = sqrt(posP.x() * posP.x() + posP.y() * posP.y() + posP.z() * posP.z()) / radius;
        //
        if (grid_domain == 1)
        {
            kp = static_cast<uint_64>(floor(log10(L / LMin) / logRatio));
        }
        else if (grid_domain == 2)
        {
            double cc;
            if (L < LMid)
            {
                cc = const_sinh1 * (L - LMin) / (LMid - LMin);
                kp = static_cast<uint_64>(floor(grid_N1 * log(cc + sqrt(cc * cc + 1.0))));
            }
            else if (L >= LMid)
            {
                cc = const_sinh2 * (L - LMid) / (LMax - LMid);
                kp = static_cast<uint_64>(floor(grid_N2 * log(cc + sqrt(cc * cc + 1.0))));
                kp = kp + (2 ^ 20);
            }
        }
        else
        {
            std::cout << " grid_domain error \n";
            std::cin.get();
        }
        //
        face = Getface(px, py, pz);
        switch (face)
        {
        case 0:
            temp[0] = py / px;
            temp[1] = pz / px;
            break;
        case 1:
            temp[0] = -px / py;
            temp[1] = pz / py;
            break;
        case 2:
            temp[0] = py / pz;
            temp[1] = -px / pz;
            break;
        case 3:
            temp[0] = py / px;
            temp[1] = -pz / px;
            break;
        case 4:
            temp[0] = -px / py;
            temp[1] = -pz / py;
            break;
        default:
            temp[0] = -py / pz;
            temp[1] = -px / pz;
            break;
        }
        // 4.3 UVtoST, note that 0<ST<1
        for (int i = 0; i <= 1; i++)
        {
            if (temp[i] >= 0)
                temp[i] = 0.5 * std::sqrt(1 + 3 * temp[i]);
            else
                temp[i] = 1 - 0.5 * std::sqrt(1 - 3 * temp[i]);
        }
        // 4.4 STtoIpJp
        ip = static_cast<unsigned int>(floor(temp[0] * particlesGridsSize));
        jp = static_cast<unsigned int>(floor(temp[1] * particlesGridsSize));
        // 5. F ip jp kp to Uint_64
        posUint = face << 61;
        for (int i = 0; i < particlesGridsLevel; i++)
        {
            posUint += (((ip >> (particlesGridsLevel - 1 - i)) & 1) << (60 - i * 3))
                    + (((jp >> (particlesGridsLevel - 1 - i)) & 1) << (60 - 1 - i * 3))
                    + (((kp >> (particlesGridsLevel - 1 - i)) & 1) << (60 - 2 - i * 3));
        }
        if( L >= LMid) {
            std::cout << " L > LMid, L = " << L << "\n" << " ip " << ip << " jp " << jp << " kp " << kp <<" posUint " << posUint;
            exit(0);
            }
        if (L >= LMid && grid_domain ==2)
            posUint = posUint + 1;
    }
    //
    inline void SetOutParticles()
    {
        posUint = ~static_cast<uint_64>(0);
    }

    inline uint_64 PosUint()
    {
        return posUint;
    }
    //
    //inline uint_64 OppoPosUint()
    //{
    //    uint_64 temp = ~posUint;
    //    return temp;
    //}
    //
    inline bool AliveParticle()
    {
        if( ~posUint == 0)
            return false;
        else
            return true;
    }

    inline Vector3 PosParticles()
    {
        return posP;
    }

    inline Vector3 VelParticles()
    {
        return vp;
    }
    // return b3 at location of particles
    inline Vector3 B3atParticles(GridsPoints *****ptrArray_in)
    {
        struct structg strg_in1 = InttoStrp1();
        Vector3 tempb;
        Vector3 tempb1 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 1][strg_in1.kg]->B3();
        Vector3 tempb2 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 1][strg_in1.kg]->B3();
        Vector3 tempb3 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 2][strg_in1.kg]->B3();
        Vector3 tempb4 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 2][strg_in1.kg]->B3();
        Vector3 tempb5 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 1][strg_in1.kg + 1]->B3();
        Vector3 tempb6 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 1][strg_in1.kg + 1]->B3();
        Vector3 tempb7 = ptrArray_in[strg_in1.face][strg_in1.ig + 2][strg_in1.jg + 2][strg_in1.kg + 1]->B3();
        Vector3 tempb8 = ptrArray_in[strg_in1.face][strg_in1.ig + 1][strg_in1.jg + 2][strg_in1.kg + 1]->B3();

        double dcellSize1 = (double)cellSize1;
        double cellSize3 = (double)cellSize1;
        double strg_in1iw = (double)strg_in1.iw;
        double strg_in1jw = (double)strg_in1.jw;
        double strg_in1kw = (double)strg_in1.kw;

        double w1 = 1.0F - 1.0F * (strg_in1iw + 0.5F) * (strg_in1jw + 0.5F) * (strg_in1kw + 0.5F) / cellSize3;
        double w2 = 1.0F - 1.0F * (dcellSize1 - strg_in1iw - 0.5F) * (strg_in1jw + 0.5F) * (strg_in1kw + 0.5F) / cellSize3;
        double w3 = 1.0F - 1.0F * (dcellSize1 - strg_in1iw - 0.5F) * (cellSize1 - strg_in1jw - 0.5F) * (strg_in1kw + 0.5F) / cellSize3;
        double w4 = 1.0F - 1.0F * (strg_in1iw + 0.5F) * (dcellSize1 - strg_in1jw - 0.5F) * (strg_in1kw + 0.5F) / cellSize3;
        double w5 = 1.0F - 1.0F * (strg_in1iw + 0.5F) * (strg_in1jw + 0.5F) * (dcellSize1 - strg_in1kw - 0.5F) / cellSize3;
        double w6 = 1.0F - 1.0F * (cellSize1 - strg_in1iw - 0.5F) * (strg_in1jw + 0.5F) * (dcellSize1 - strg_in1kw - 0.5F) / cellSize3;
        double w7 = 1.0F - 1.0F * (cellSize1 - strg_in1iw - 0.5F) * (dcellSize1 - strg_in1jw - 0.5F) * (cellSize1 - strg_in1kw - 0.5F) / cellSize3;
        double w8 = 1.0F - 1.0F * (strg_in1iw + 0.5F) * (dcellSize1 - strg_in1jw - 0.5F) * (dcellSize1 - strg_in1kw - 0.5F) / cellSize3;
        tempb.Setx(tempb1.x() * w1 + tempb2.x() * w2 + tempb3.x() * w3 + tempb4.x() * w4 + tempb5.x() * w5 + tempb6.x() * w6 + tempb7.x() * w7 + tempb8.x() * w8);
        tempb.Sety(tempb1.y() * w1 + tempb2.y() * w2 + tempb3.y() * w3 + tempb4.y() * w4 + tempb5.y() * w5 + tempb6.y() * w6 + tempb7.y() * w7 + tempb8.y() * w8);
        tempb.Setz(tempb1.z() * w1 + tempb2.z() * w2 + tempb3.z() * w3 + tempb4.z() * w4 + tempb5.z() * w5 + tempb6.z() * w6 + tempb7.z() * w7 + tempb8.z() * w8);
        return tempb;
    }
    // return abs vel including mu information
    // used for calculating collisions
    // b3cell is magnetic field at particles' location
    inline Vector3 VelCollParticles(Vector3 b3cell, double mi0)
    {
        Vector3 velAbs;
        Vector3 velEigen1, velEigen2;
        //double lxy = sqrt(b3cell.x() * b3cell.x() + b3cell.y() * b3cell.y());
        //double r = b3cell.norm();
        double Pi = 3.14159265359F;
        Vector3 b_unit = b3cell.NormalizedVector();
        // find two eigen vector perpendicular to b3cell
        if (b3cell.x() == 0.0F && b3cell.y() == 0.0F)
        {
            velEigen1 = Vector3(1.0F, 0.0F, 0.0F);
            velEigen2 = Vector3(0.0F, 1.0F, 0.0F);
        }
        else
        {
            // set velEigen1: parallel to xOy plane
            velEigen1 = Vector3(1.0F, -b3cell.x() / b3cell.y(), 0.0F);
            velEigen1 = velEigen1.NormalizedVector();
            // Calculate velEigen2
            velEigen2 = velEigen1.CrossProduct(b3cell).NormalizedVector();
        }
        //
        double angle = dRand() * Pi * 2.0F;
        double vel_mu = sqrt(mu * b3cell.norm() * 2.0F / mi0);
        Vector3 vel_perp = velEigen1.ScaleProduct(vel_mu * cos(angle)).PlusProduct(velEigen2.ScaleProduct(vel_mu * sin(angle)));
        Vector3 v_para = b_unit.ScaleProduct(vp.DotProduct(b_unit));
        //
        return vel_perp.PlusProduct(v_para);
    }
    // drift velocity
    inline Vector3 VelDriftParticles(Vector3 b3cell)
    {
        Vector3 b_unit = b3cell.NormalizedVector();
        Vector3 v_para = b_unit.ScaleProduct(vp.DotProduct(b_unit));
        return vp.MinusProduct(v_para);
    }

    inline double WeightNi()
    {
        return weightNi;
    }

    inline void SetWeightNi(double weight_new)
    {
        weightNi = weight_new;
    }
    inline void ResetWeightNi() // reset weight for testing openmp
    {
        weightNi = 0.0;
    }

    inline double MagneticIvarient()
    {
        return mu;
    }

    inline void SetPosPar( Vector3 testPosPar)
    {
        posP = testPosPar;
        // std::cout << " warning SetPosPar applied \n" ;
    }
    inline void SetParticles(const Particles &other_par)
    {
        posUint = other_par.posUint;
        posP = other_par.posP;
        vp = other_par.vp;
        weightNi = other_par.weightNi;
        mu = other_par.mu;
    }

    inline void SetParticles(unsigned long long &posUint_in,
                             double &posPx_in,
                             double &posPy_in,
                             double &posPz_in,
                             double &vx_in,
                             double &vy_in,
                             double &vz_in,
                             double &weight_in,
                             double &mu_in)
    {
        posUint = posUint_in;
        posP = Vector3(posPx_in, posPy_in, posPz_in);
        vp = Vector3(vx_in, vy_in, vz_in);
        weightNi = weight_in;
        mu = mu_in;
    }

    inline void UpdateTempPar(Vector3 &Vel, double &mu_simu)
    {
        vp = Vel;
        mu = mu_simu;
    }
    //
    inline void SetVelocity(Vector3 Vel)
    {
        vp = Vel;
    }
    inline void SetMu( double _mu)
    {
        mu = _mu;
    }
    // set velocity and mu
    inline void SetMuVelocity(Vector3 velAbs, Vector3 velDrift, Vector3 b3cell, double mi0)
    {
        // velDrift should be perpendicular to b3cell
        // velAbs = vel_para + vel_gyro + velDrift
        // vel_perp = vel_gyro + velDrift
        vp = b3cell.NormalizedVector().ScaleProduct(velAbs.DotProduct(b3cell.NormalizedVector())).PlusProduct(velDrift);
        Vector3 vel_perp = velAbs.MinusProduct(vp);
        mu = vel_perp.norm2() * 0.5 * mi0 / b3cell.norm();

        //Vector3 vel_para = b3cell.NormalizedVector().ScaleProduct( velAbs.DotProduct(b3cell.NormalizedVector()));
        //Vector3 vel_perp = velAbs.MinusProduct(vel_para);
        //Vector3 vel_gyro = vel_perp.MinusProduct( velDrift);
        //mu = vel_gyro.norm2() * mi0 * 0.5 / b3cell.norm();
        //vp = vel_para.PlusProduct(velDrift);
    }
    //////////////////////////////////Constructor//////////////////////////////
    Particles(const Particles &other_par);
    Particles(uint_64 posUint_in,
              Vector3 posP_in,
              Vector3 vp_in,
              double weightNi_in,
              double mu_in);
    // FUNCTION //Default Constructor
    Particles();

private:
    uint_64 posUint; // single unsigned int for position
    Vector3 posP;
    Vector3 vp; // velocity of particles guiding center, including drift velocity and parallal vel

    double weightNi; // weight for number of real particles
    double mu;       // magnetic moment/ adiabatic invarient
    //    double px; double py; double pz;
    //    double vx; double vy; double vz;
    //    int face; int ri; int rj; int rk; //face, (i,j) in fieldsgrids, radial
};
#endif