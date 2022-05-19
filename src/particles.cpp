#include <iostream>
#include "mathutil.h"
#include "parameters.h"
#include "particles.h"
#include "vector3.h"
#include "structdef.h"
#include <bitset>
#include "module_base.h"
#include <vector>
// private:
// uint_64 posUint;
// Vector3 vp;

    Particles::Particles(const Particles &other_par)
    {
        posUint = other_par.posUint;
        posP = other_par.posP;
        vp = other_par.vp;
        weightNi = other_par.weightNi;
        mu = other_par.mu;
    }
    Particles::Particles(uint_64 posUint_in,
              Vector3 posP_in,
              Vector3 vp_in,
              double weightNi_in,
              double mu_in)
    {
        posUint = posUint_in;
        posP = posP_in;
        vp = vp_in;
        weightNi = weightNi_in;
        mu = mu_in;
    }

    // FUNCTION //Default Constructor
    Particles::Particles()
    {
        posUint = 0;
        posP = Vector3(0.0, 0.0, 0.0);
        vp = Vector3(0.0, 0.0, 0.0);
        weightNi = 0.0;
        mu = 0.0;
    }
//************************************************************************
//************************************************************************
// FUNCTION // Update uint_64 IN
// And return a int "0" means in the main domain
// "1" means out of the main domain
//************************************************************************
//************************************************************************
int Particles::UpdateUint_64()
{
    uint_64 face = 0, ip = 0, jp = 0, kp = 0;
    double px, py, pz;
    double temp[2];
    Vector3 posP_saved = posP;
    px = posP.x();
    py = posP.y();
    pz = posP.z();
    // update double x y z, vp is velocity
    px += vp.x() * tstep;
    py += vp.y() * tstep;
    pz += vp.z() * tstep;
    posP = Vector3(px, py, pz);
    // 4. transfor to face ip kp jp
    // 4.1 radial kp
    double L = sqrt(posP.x() * posP.x() + posP.y() * posP.y() + posP.z() * posP.z()) / radius;
    // check if in our simulation domain
    if(L > LMax)
    {
        if (LatMagParticles() == 0) // in closed field lines, assume the particles trans into the other hemisphere
        {
            posP = Vector3(posP_saved.x(), posP_saved.y(), -posP_saved.z());
            vp = Vector3(-vp.x(), -vp.y(), vp.z());
            return 0;
        } else
        { 
            // set out particles
            //std::cout << " return1 type1 " << L << " " << LMax << " " << LMax1 << "\n";
            return 1;
        }
    }else if(L < LMin)
    {
        return 1;
    }else
    {
        if(grid_domain == 1)
        {
            kp = static_cast<uint_64>(floor(log10(L / LMin) / logRatio));
        }
        else if(grid_domain == 2)
        {
            double cc;
            if(L < LMid)
            {
                cc = const_sinh1 * (L - LMin) / (LMid - LMin);
                kp = static_cast<uint_64>(floor(grid_N1 * log(cc + sqrt(cc * cc + 1.0))));
            }
            else if(L >= LMid)
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
        face = Getface(posP.x(), posP.y(), posP.z());
        switch (face)
        {
        case 0:
            temp[0] = posP.y() / posP.x();
            temp[1] = posP.z() / posP.x();
            break;
        case 1:
            temp[0] = -posP.x() / posP.y();
            temp[1] = posP.z() / posP.y();
            break;
        case 2:
            temp[0] = posP.y() / posP.z();
            temp[1] = -posP.x() / posP.z();
            break;
        case 3:
            temp[0] = posP.y() / posP.x();
            temp[1] = -posP.z() / posP.x();
            break;
        case 4:
            temp[0] = -posP.x() / posP.y();
            temp[1] = -posP.z() / posP.y();
            break;
        default:
            temp[0] = -posP.y() / posP.z();
            temp[1] = -posP.x() / posP.z();
            break;
        }
        // 4.3 UVtoST, note that -1<UV<1, 0<ST<1
        // tan(theta) = temp[i], theta range [-pi/4,pi/4]
        temp[0] = atan(temp[0]) * 2.0 / pi + 0.5;
        temp[1] = atan(temp[1]) * 2.0 / pi + 0.5;
        // // faster access version
        // for(int i = 0; i <= 1; i++)
        // {
        //     if (temp[i] >= 0)
        //         temp[i] = 0.5 * sqrt(1.0 + 3.0 * temp[i]);
        //     else
        //         temp[i] = 1.0 - 0.5 * sqrt(1.0 - 3.0 * temp[i]);
        // }
        // 4.4 STtoIpJp
        ip = static_cast<uint_64>(floor(temp[0] * particlesGridsSize));
        jp = static_cast<uint_64>(floor(temp[1] * particlesGridsSize));
        // 5. F ip jp kp to Uint_64
        posUint = static_cast<uint_64>(0);
        posUint = face << 61;
        for(int i = 0; i < particlesGridsLevel; i++)
        {
            posUint += (((ip >> (particlesGridsLevel - 1 - i)) & 1) << (60 - i * 3))
                    + (((jp >> (particlesGridsLevel - 1 - i)) & 1) << (59 - i * 3))
                    + (((kp >> (particlesGridsLevel - 1 - i)) & 1) << (58 - i * 3));
        }
        if(L >= LMid && grid_domain == 2)
        {
            posUint = posUint + 1;
            std::cout << " grid_domain = 2 \n";
            exit(2);
        }
        //
        return 0;
    }
}

//**************************test******************************************
int Particles::UpdateUint_64_test()
{
    uint_64 face = 0, ip = 0, jp = 0, kp = 0;
    double px, py, pz;
    double temp[2];
    Vector3 posP_saved = posP;
    px = posP.x();
    py = posP.y();
    pz = posP.z();
    // update double x y z, vp is velocity
    px += vp.x() * tstep;
    py += vp.y() * tstep;
    pz += vp.z() * tstep;
    // test for fixed particles  ?** posP = Vector3(px, py, pz);
    // 4. transfor to face ip kp jp
    // 4.1 radial kp
    double L = sqrt(posP.x() * posP.x() + posP.y() * posP.y() + posP.z() * posP.z()) / radius;
    // check if in our simulation domain
    if(L > LMax)
    {
        if (LatMagParticles() == 0) // in closed field lines, assume the particles trans into the other hemisphere
        {
            posP = Vector3(posP_saved.x(), posP_saved.y(), -posP_saved.z());
            vp = Vector3(-vp.x(), -vp.y(), vp.z());
            return 0;
        } else
        { 
            // set out particles
            //std::cout << " return1 type1 " << L << " " << LMax << " " << LMax1 << "\n";
            return 1;
        }
    }else if(L < LMin)
    {
        //std::cout << " return1 type2\n";
        return 1;
    }else
    {
        if(grid_domain == 1)
        {
            kp = static_cast<uint_64>(floor(log10(L / LMin) / logRatio));
        }
        else if(grid_domain == 2)
        {
            double cc;
            if(L < LMid)
            {
                cc = const_sinh1 * (L - LMin) / (LMid - LMin);
                kp = static_cast<uint_64>(floor(grid_N1 * log(cc + sqrt(cc * cc + 1.0))));
            }
            else if(L >= LMid)
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
        face = Getface(posP.x(), posP.y(), posP.z());
        switch (face)
        {
        case 0:
            temp[0] = posP.y() / posP.x();
            temp[1] = posP.z() / posP.x();
            break;
        case 1:
            temp[0] = -posP.x() / posP.y();
            temp[1] = posP.z() / posP.y();
            break;
        case 2:
            temp[0] = posP.y() / posP.z();
            temp[1] = -posP.x() / posP.z();
            break;
        case 3:
            temp[0] = posP.y() / posP.x();
            temp[1] = -posP.z() / posP.x();
            break;
        case 4:
            temp[0] = -posP.x() / posP.y();
            temp[1] = -posP.z() / posP.y();
            break;
        default:
            temp[0] = -posP.y() / posP.z();
            temp[1] = -posP.x() / posP.z();
            break;
        }
        // 4.3 UVtoST, note that -1<UV<1, 0<ST<1
        // tan(theta) = temp[i], theta range [-pi/4,pi/4]
        temp[0] = atan(temp[0]) * 2.0 / pi + 0.5;
        temp[1] = atan(temp[1]) * 2.0 / pi + 0.5;
        // // faster access version
        // for(int i = 0; i <= 1; i++)
        // {
        //     if (temp[i] >= 0)
        //         temp[i] = 0.5 * sqrt(1.0 + 3.0 * temp[i]);
        //     else
        //         temp[i] = 1.0 - 0.5 * sqrt(1.0 - 3.0 * temp[i]);
        // }
        // 4.4 STtoIpJp
        ip = static_cast<uint_64>(floor(temp[0] * particlesGridsSize));
        jp = static_cast<uint_64>(floor(temp[1] * particlesGridsSize));
        // 5. F ip jp kp to Uint_64
        posUint = static_cast<uint_64>(0);
        posUint = face << 61;
        for(int i = 0; i < particlesGridsLevel; i++)
        {
            posUint += (((ip >> (particlesGridsLevel - 1 - i)) & 1) << (60 - i * 3))
                    + (((jp >> (particlesGridsLevel - 1 - i)) & 1) << (59 - i * 3))
                    + (((kp >> (particlesGridsLevel - 1 - i)) & 1) << (58 - i * 3));
        }
        if(L >= LMid && grid_domain == 2)
        {
            posUint = posUint + 1;
            std::cout << " grid_domain = 2 \n";
            exit(2);
        }
        //
        return 0;
    }
}

//************************************************************************
//************************************************************************
// FUNCTION // Update uint_64 IN
// And return a int "0" means in the main domain
// // "1" means out of the main domain
// int Particles::UpdateUint_64_temp()
// {
//     uint_64 face = 0, ip = 0, jp = 0, kp = 0;
//     double px, py, pz;
//     double temp[2];
//     int check = 0;

//     px = posP.x();
//     py = posP.y();
//     pz = posP.z();
//     // 3. update double x y z
//     px += vp.x() * tstep;
//     py += vp.y() * tstep;
//     pz += vp.z() * tstep;

//     posP = Vector3(px, py, pz);

//     // 4. transfor to face ip kp jp
//     // 4.1 radial kp
//     double L = sqrt(px * px + py * py + pz * pz) / radius;
//     //    double L = sqrt( posP.x() * posP.x() + posP.y() * posP.y() + posP.z() * posP.z())/radius;

//     // check if in the main domain
//     if (L > LMax || L < LMin)
//     {
//         // in this case, set out of whole domain
//         //        if(L > LMax_centraldomain) std::cout << "p";
//         //        if(L < LMin_centraldomain) std::cout << "q";
//         check = 2;
//         return check;
//     }
//     else
//     {
//         double tempRandom = (double)rand() / (RAND_MAX);
//         if ((L > LMax_centraldomain || L < LMin_centraldomain) && tempRandom < neutralizeRateCover)
//         {
//             // out the central domain, check random number to set out the particle
//             //            if(L > LMax_centraldomain) std::cout << "p";
//             //            if(L < LMin_centraldomain) std::cout << "q";
//             check = 2;
//             return check;
//         }
//         else
//         {
//             //        kp = static_cast<uint_64>( floor( log10(L / LMin)/logRatio *cellSize1));
//             if (grid_domain == 1)
//                 kp = static_cast<uint_64>(floor(log10(L / LMin) / logRatio));
//             else if (grid_domain == 2)
//             {
//                 double cc;
//                 if (L < LMid)
//                 {
//                     cc = const_sinh1 * (L - LMin) / (LMid - LMin);
//                     kp = static_cast<uint_64>(floor(grid_N1 * log(cc + sqrt(cc * cc + 1.0))));
//                 }
//                 else if (L >= LMid)
//                 {
//                     cc = const_sinh2 * (L - LMid) / (LMax - LMid);
//                     kp = static_cast<uint_64>(floor(grid_N2 * log(cc + sqrt(cc * cc + 1.0))));
//                     kp = kp + (2 ^ 20);
//                 }
//             }
//             else
//             {

//                 std::cout << " grid_domain error \n";
//                 std::cin.get();
//             }

//             // 4.2 XYZtoUV, note that -1<UV<1
//             face = Getface(px, py, pz);
//             switch (face)
//             {
//             case 0:
//                 temp[0] = py / px;
//                 temp[1] = pz / px;
//                 break;
//             case 1:
//                 temp[0] = -px / py;
//                 temp[1] = pz / py;
//                 break;
//             case 2:
//                 temp[0] = py / pz;
//                 temp[1] = -px / pz;
//                 break;
//             case 3:
//                 temp[0] = py / px;
//                 temp[1] = -pz / px;
//                 break;
//             case 4:
//                 temp[0] = -px / py;
//                 temp[1] = -pz / py;
//                 break;
//             default:
//                 temp[0] = -py / pz;
//                 temp[1] = -px / pz;
//                 break;
//             }
//             // 4.3 UVtoST, note that 0<ST<1
//             for (int i = 0; i <= 1; i++)
//             {
//                 if (temp[i] >= 0)
//                     temp[i] = 0.5 * std::sqrt(1 + 3 * temp[i]);
//                 else
//                     temp[i] = 1 - 0.5 * std::sqrt(1 - 3 * temp[i]);
//             }
//             // 4.4 STtoIpJp
//             ip = static_cast<unsigned int>(floor(temp[0] * particlesGridsSize));
//             jp = static_cast<unsigned int>(floor(temp[1] * particlesGridsSize));

//             // 5. F ip jp kp to Uint_64
//             posUint = face << 61;
//             for (int i = 0; i < particlesGridsLevel; i++)
//             {
//                 posUint += (((ip >> (particlesGridsLevel - 1 - i)) & 1) << (60 - i * 3)) 
//                         + (((jp >> (particlesGridsLevel - 1 - i)) & 1) << (60 - 1 - i * 3))
//                         + (((kp >> (particlesGridsLevel - 1 - i)) & 1) << (60 - 2 - i * 3));
//             }
//             if (L >= LMid)
//                 posUint = posUint + 1;
//             /*        std::cout << " next " << std::endl;
//     std::cout << std::bitset<64>(face) << " " << face << std::endl;
//     std::cout << std::bitset<64>(ip) << " " << ip << std::endl;
//     std::cout << std::bitset<64>(jp) << " " << jp << std::endl;
//     std::cout << std::bitset<64>(kp) << " " << kp << std::endl << std::endl;
//     int pause;
//     std::cin >>pause;
// */
//             if (L > LMax_centraldomain || L < LMin_centraldomain)
//             {
//                 // still in the temp domain
//                 check = 1;
//             }
//             else if (L < LMax_centraldomain && L > LMin_centraldomain)
//             {
//                 // in the central domain
//                 check = 0;

//                 //           std::cout << " check =0 " << vp.x() << " " << vp.y() << " " << vp.z() << std::endl;
//                 //        if(L > LMax_centraldomain) std::cout << "p";
//                 //        if(L < LMin_centraldomain) std::cout << "q";
//             }

//             return check;
//         }
//     }
// }

//************************************************************************
//************************************************************************
// FUNCTION // Boris' Method update velocity
// And return a int "0" means in the main domain
// // "1" means out of the main domain
// //  tried a pre-stored grids for faster accessing, but not successful
// //************************************************************************
// //************************************************************************
// int Particles::BorisMethod(struct structPar &strg_in, std::vector<GridsPoints> gridsPointsAtVertex, int type_a, int maindomain)
// {
//     double mi0_in = 0.0;
//     if( type_a == 0)
//     {
//         mi0_in = mi0_H;
//     }else if (type_a ==1)
//     {
//         mi0_in = mi0_He;
//     }else if(type_a ==2)
//     {
//         mi0_in = mi0_O;
//     }else
//     {
//         std::cout << " BorisMethod error \n";
//     }
    
//     // 1. pre-set some variables
//     // 1.1 qtm = q * dt / m /2
//     double gama = sqrt(1.0 - vp.norm2() / lightSpeed / lightSpeed);
//     mi0_in = mi0_in / gama;
//     double qtm = qi0 * tstep / mi0_in / 2.0;

//     // 1.2 local B and E
//     //    std::cout << std::endl;
//     //    std::cout << strg_in.face << strg_in.ig << strg_in.jg << strg_in.kg << std::endl;
    
//     Vector3 temp1 = gridsPointsAtVertex[0].B3();
//     Vector3 temp2 = gridsPointsAtVertex[1].B3();
//     Vector3 temp3 = gridsPointsAtVertex[2].B3();
//     Vector3 temp4 = gridsPointsAtVertex[3].B3();
//     Vector3 temp5 = gridsPointsAtVertex[4].B3();
//     Vector3 temp6 = gridsPointsAtVertex[5].B3();
//     Vector3 temp7 = gridsPointsAtVertex[6].B3();
//     Vector3 temp8 = gridsPointsAtVertex[7].B3();

//     Vector3 tempb;
//     tempb.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

//     tempb.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

//     tempb.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

//     //                // test test ***********
//     //                double r = sqrt(pow(posP.x(),2.0) + pow(posP.y(),2.0) + pow(posP.z(),2.0));
//     //                tempb.Setx( 3 * dMoment * posP.x() * posP.z() / pow(r,5.0));
//     //                tempb.Sety( 3 * dMoment * posP.y() * posP.z() / pow(r,5.0));
//     //                tempb.Setz( dMoment * (3 * pow(posP.z(),2.0) - pow(r,2.0)) / pow(r,5.0));

//     //    std::cout << " tempb " << tempb.x() << " " << temp1.x() << " ===>>> " ;
//     temp1 = gridsPointsAtVertex[0].E3();
//     temp2 = gridsPointsAtVertex[1].E3();    
//     temp3 = gridsPointsAtVertex[2].E3();
//     temp4 = gridsPointsAtVertex[3].E3();
//     temp5 = gridsPointsAtVertex[4].E3();    
//     temp6 = gridsPointsAtVertex[5].E3();    
//     temp7 = gridsPointsAtVertex[6].E3();    
//     temp8 = gridsPointsAtVertex[7].E3();    
    
//     Vector3 tempe;
//     tempe.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

//     tempe.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

//     tempe.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

//     // test test ***********

//     //                Vector3 tempOmega1 = Vector3( 0.0, 0.0, omega_earth);
//     //                Vector3 temptemppos = posP;
//     //                temptemppos.Setz(0.0);
//     //                tempe = tempb.CrossProduct(tempOmega1.CrossProduct(temptemppos));

//     //    std::cout << ptrArray_in[strg_in.face][strg_in.ig+1][strg_in.jg+1][strg_in.kg]->E3().x() << " ==>> ";
//     temp1 = gridsPointsAtVertex[0].GradB3();
//     temp2 = gridsPointsAtVertex[1].GradB3();    
//     temp3 = gridsPointsAtVertex[2].GradB3();
//     temp4 = gridsPointsAtVertex[3].GradB3();
//     temp5 = gridsPointsAtVertex[4].GradB3();    
//     temp6 = gridsPointsAtVertex[5].GradB3();    
//     temp7 = gridsPointsAtVertex[6].GradB3();    
//     temp8 = gridsPointsAtVertex[7].GradB3();    

//     // 2 Revised E ( including gradB at the location of unique particles)
//     Vector3 tempGradB;
//     tempGradB.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);
//     tempGradB.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);
//     tempGradB.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

//     temp1 = gridsPointsAtVertex[0].Pos3();
//     temp2 = gridsPointsAtVertex[1].Pos3();    
//     temp3 = gridsPointsAtVertex[2].Pos3();
//     temp4 = gridsPointsAtVertex[3].Pos3();
//     temp5 = gridsPointsAtVertex[4].Pos3();    
//     temp6 = gridsPointsAtVertex[5].Pos3();    
//     temp7 = gridsPointsAtVertex[6].Pos3();    
//     temp8 = gridsPointsAtVertex[7].Pos3();    
//     // 1 Revised E ( including gravity which is not accurate to the loation of unique particles)
//     Vector3 tempPos;

//     tempPos.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

//     tempPos.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

//     tempPos.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

//     // ve
//     temp1 = gridsPointsAtVertex[0].Vel_e3();
//     temp2 = gridsPointsAtVertex[1].Vel_e3();    
//     temp3 = gridsPointsAtVertex[2].Vel_e3();
//     temp4 = gridsPointsAtVertex[3].Vel_e3();
//     temp5 = gridsPointsAtVertex[4].Vel_e3();    
//     temp6 = gridsPointsAtVertex[5].Vel_e3();    
//     temp7 = gridsPointsAtVertex[6].Vel_e3();    
//     temp8 = gridsPointsAtVertex[7].Vel_e3();   
//     Vector3 tempVe3;

//     tempVe3.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

//     tempVe3.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

//     tempVe3.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

//     // vi
//     temp1 = gridsPointsAtVertex[0].Vel3();
//     temp2 = gridsPointsAtVertex[1].Vel3();    
//     temp3 = gridsPointsAtVertex[2].Vel3();
//     temp4 = gridsPointsAtVertex[3].Vel3();
//     temp5 = gridsPointsAtVertex[4].Vel3();    
//     temp6 = gridsPointsAtVertex[5].Vel3();    
//     temp7 = gridsPointsAtVertex[6].Vel3();    
//     temp8 = gridsPointsAtVertex[7].Vel3();  
//     Vector3 tempVi3;

//     tempVi3.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

//     tempVi3.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

//     tempVi3.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

//     // ambipolar diffusion ? what does it mean?
//     //Vector3 diffusionVector;
//     //double L = sqrt(posP.x() * posP.x() + posP.y() * posP.y() + posP.z() * posP.z()) / radius;
//     //double diffusionFactor = exp(LMin - 0.0001 - L) * diffusionCollisionCoefficient;
//     //diffusionVector = tempVe3.MinusProduct(tempVi3).ScaleProduct(diffusionFactor);
//     //diffusionVector = diffusionVector.ScaleProduct(mi0_in / qi0);
//     // Definition centrifugal && coriolis
//     Vector3 centrifugal = Vector3(0.0, 0.0, 0.0);
//     Vector3 coriolis = Vector3(0.0, 0.0, 0.0);
//     Vector3 tempGravity = Vector3(0.0, 0.0, 0.0);
//     Vector3 tempOmega = Vector3(omega_wx, omega_wy, omega_wz);
//     // *(- mu /m * m / e) = - mu / e , where mu = m * v_perp^2 / 2 / B
//     if (mirror_force == 1)
//         tempGradB = tempGradB.ScaleProduct(-1.0 * mu / qi0); // grad B: mirror force effects
//     // Centrifugal force: need to figure out the accleration of rotation first
//     // a = w x ( w x r ), then the centrifugal force (a) is: a = - w x ( w x r )
//     // Coriolis force is negligable to centrifugal force
//     if (centrifugal_force == 1)
//         centrifugal = tempOmega.CrossProduct(tempPos).CrossProduct(tempOmega); // acceleration of centrifugal: (w X r) X w
//     // Coriolis force
//     if (coriolis_force == 1)
//         coriolis = vp.CrossProduct(tempOmega).ScaleProduct(2.0); // acceleration of  coriolis: -2(w X v)
//     // new gravity
//     if (gravity_force == 1)
//         tempGravity = tempPos.NormalizedVector().ScaleProduct(-1.0 * gravity * radius * radius / tempPos.norm2());
//     //std::cout << " test force magnitude " << " tempE " << tempe.DotProduct(tempb.NormalizedVector()) << " gravity " << tempGravity.DotProduct(tempb.NormalizedVector()) <<
//     //            " centrifugal " << centrifugal.DotProduct(tempb.NormalizedVector()) << " coriolis " << coriolis.DotProduct(tempb.NormalizedVector()) << 
//     //            " gradB " << tempGradB.DotProduct(tempb.NormalizedVector()) << " mu force  " << tempGradB.DotProduct(tempb.NormalizedVector())* mu / mi0_in<<  " \n";
//     //std::cin.get();
//     // put together
//     if (coordinate_rotate == 1)
//         tempGravity = tempGravity.PlusProduct(centrifugal).PlusProduct(coriolis);
//     //
//     tempGravity = tempGravity.ScaleProduct(mi0_in / qi0); // Gravity effects into reduced E
//     // revised E including gravity and gradB, and centrifugal and coriolis force if initilized particles
//     if (electric_force == 0)
//         tempe = Vector3(0.0, 0.0, 0.0);
//     tempe = tempe.PlusProduct(tempGravity).PlusProduct(tempGradB); // .PlusProduct( diffusionVector);
//     // test only Grad B
//     tempe = tempGradB;

//     //   tempe = tempe.MinusProduct(tempb.NormalizedVector().ScaleProduct( tempe.DotProduct(tempb.NormalizedVector()) )).PlusProduct(tempGradB);
//     //   tempe = tempe.MinusProduct(tempb.NormalizedVector().ScaleProduct( tempe.DotProduct(tempb.NormalizedVector()) ));
//     // test particles only under mirror force
//     //if( move_type ==0)
//     //tempe = tempe.MinusProduct(tempb.NormalizedVector().ScaleProduct( tempe.DotProduct(tempb.NormalizedVector()) )).PlusProduct(tempGradB); // bounce motion
//     //    tempe = tempGradB;
//     // test gravity
//     //    tempe = tempPos;
//     // test mirror force effects
//     //    tempe = tempGradB;
//     //
//     // 1.3 Vector3 t ; Vector3 s

//     // test *************************************
//     //tempb = { 0.0,0.0,0.0};
//     //vp = { 0.0,0.0,0.0};
//     // 
//     Vector3 t = tempb.ScaleProduct(qtm);
//     Vector3 s = t.ScaleProduct(2.0 / (1.0 + t.norm2()));
//     // 2. Boris equations
//     // 2.1 equation I: v1 = v + E * qtm
//     Vector3 v1 = vp.PlusProduct(tempe.ScaleProduct(qtm));
//     // 2.2 equation II: v2 = v1 + v1 X t
//     Vector3 v2 = v1.PlusProduct(v1.CrossProduct(t));
//     // 2.3 equation III: v3 = v1 + v2 X s
//     Vector3 v3 = v1.PlusProduct(v2.CrossProduct(s));
//     // 2.4 equation IV: v = v3 + E * qtm
//     vp = v3.PlusProduct(tempe.ScaleProduct(qtm));
//     // 3. update the postion
//     if (maindomain == 0)
//     {
//         int check = UpdateUint_64();
//         return check;
//     }
//     else
//     {
//         std::cout << " maindomain != 0 \n";
//         exit(0);
//         return UpdateUint_64_temp();
//     }
// }

// YH 4/1 
// return potential at postition of unreal particle, method is similar to the beginning part
// of BorisMethod() below
double Particles::PotentialAtUnrealPar(struct structPar &strg_in, GridsPoints *****ptrArray_in)
{
    double temp1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Potential();
    double temp2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->Potential();
    double temp3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->Potential();
    double temp4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->Potential();
    double temp5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->Potential();
    double temp6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->Potential();
    double temp7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->Potential();
    double temp8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->Potential();

    double tempPotential = temp1 * strg_in.w000 + temp2 * strg_in.w100 + temp3 * strg_in.w110 + temp4 * strg_in.w010 + temp5 * strg_in.w001 + temp6 * strg_in.w101 + temp7 * strg_in.w111 + temp8 * strg_in.w011;
    return tempPotential;
}
//
int Particles::BorisMethod(struct structPar &strg_in, GridsPoints *****ptrArray_in, double mi0_in, int maindomain)
{
    // 1. pre-set some variables
    // 1.1 qtm = q * dt / m /2

    double gama = sqrt(1.0 - vp.norm2() / lightSpeed / lightSpeed);
    mi0_in = mi0_in / gama;
    double qtm = qi0 * tstep / mi0_in / 2.0;

    // 1.2 local B and E
    //    std::cout << std::endl;
    //    std::cout << strg_in.face << strg_in.ig << strg_in.jg << strg_in.kg << std::endl;

    Vector3 temp1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 temp2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->B3();
    Vector3 temp3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 temp4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->B3();
    Vector3 temp5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 temp6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->B3();
    Vector3 temp7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->B3();
    Vector3 temp8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->B3();

    Vector3 tempb;
    tempb.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

    tempb.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

    tempb.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

    //                // test test ***********
    //                double r = sqrt(pow(posP.x(),2.0) + pow(posP.y(),2.0) + pow(posP.z(),2.0));
    //                tempb.Setx( 3 * dMoment * posP.x() * posP.z() / pow(r,5.0));
    //                tempb.Sety( 3 * dMoment * posP.y() * posP.z() / pow(r,5.0));
    //                tempb.Setz( dMoment * (3 * pow(posP.z(),2.0) - pow(r,2.0)) / pow(r,5.0));

    //    std::cout << " tempb " << tempb.x() << " " << temp1.x() << " ===>>> " ;

    temp1 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->E3());
    temp2 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->E3());
    temp3 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->E3());
    temp4 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->E3());
    temp5 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->E3());
    temp6 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->E3());
    temp7 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->E3());
    temp8 = Vector3(ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->E3());

    Vector3 tempe;
    tempe.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

    tempe.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

    tempe.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

    // test test ***********

    //                Vector3 tempOmega1 = Vector3( 0.0, 0.0, omega_earth);
    //                Vector3 temptemppos = posP;
    //                temptemppos.Setz(0.0);
    //                tempe = tempb.CrossProduct(tempOmega1.CrossProduct(temptemppos));

    //    std::cout << ptrArray_in[strg_in.face][strg_in.ig+1][strg_in.jg+1][strg_in.kg]->E3().x() << " ==>> ";

    temp1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->GradB3();
    temp2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->GradB3();
    temp3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->GradB3();
    temp4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->GradB3();
    temp5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->GradB3();
    temp6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->GradB3();
    temp7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->GradB3();
    temp8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->GradB3();

    // 2 Revised E ( including gradB at the location of unique particles)
    Vector3 tempGradB;
    tempGradB.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);
    tempGradB.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);
    tempGradB.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

    temp1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Pos3();
    temp2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->Pos3();
    temp3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->Pos3();
    temp4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->Pos3();
    temp5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->Pos3();
    temp6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->Pos3();
    temp7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->Pos3();
    temp8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->Pos3();

    // 1 Revised E ( including gravity which is not accurate to the loation of unique particles)
    Vector3 tempPos;

    tempPos.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

    tempPos.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

    tempPos.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

    // ve
    temp1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Vel_e3();
    temp2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->Vel_e3();
    temp3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->Vel_e3();
    temp4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->Vel_e3();
    temp5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->Vel_e3();
    temp6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->Vel_e3();
    temp7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->Vel_e3();
    temp8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->Vel_e3();

    Vector3 tempVe3;

    tempVe3.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

    tempVe3.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

    tempVe3.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

    // vi
    temp1 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg]->Vel3();
    temp2 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg]->Vel3();
    temp3 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg]->Vel3();
    temp4 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg]->Vel3();
    temp5 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 1][strg_in.kg + 1]->Vel3();
    temp6 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 1][strg_in.kg + 1]->Vel3();
    temp7 = ptrArray_in[strg_in.face][strg_in.ig + 2][strg_in.jg + 2][strg_in.kg + 1]->Vel3();
    temp8 = ptrArray_in[strg_in.face][strg_in.ig + 1][strg_in.jg + 2][strg_in.kg + 1]->Vel3();

    Vector3 tempVi3;

    tempVi3.Setx(temp1.x() * strg_in.w000 + temp2.x() * strg_in.w100 + temp3.x() * strg_in.w110 + temp4.x() * strg_in.w010 + temp5.x() * strg_in.w001 + temp6.x() * strg_in.w101 + temp7.x() * strg_in.w111 + temp8.x() * strg_in.w011);

    tempVi3.Sety(temp1.y() * strg_in.w000 + temp2.y() * strg_in.w100 + temp3.y() * strg_in.w110 + temp4.y() * strg_in.w010 + temp5.y() * strg_in.w001 + temp6.y() * strg_in.w101 + temp7.y() * strg_in.w111 + temp8.y() * strg_in.w011);

    tempVi3.Setz(temp1.z() * strg_in.w000 + temp2.z() * strg_in.w100 + temp3.z() * strg_in.w110 + temp4.z() * strg_in.w010 + temp5.z() * strg_in.w001 + temp6.z() * strg_in.w101 + temp7.z() * strg_in.w111 + temp8.z() * strg_in.w011);

    // ambipolar diffusion ? what does it mean?
    //Vector3 diffusionVector;
    //double L = sqrt(posP.x() * posP.x() + posP.y() * posP.y() + posP.z() * posP.z()) / radius;
    //double diffusionFactor = exp(LMin - 0.0001 - L) * diffusionCollisionCoefficient;
    //diffusionVector = tempVe3.MinusProduct(tempVi3).ScaleProduct(diffusionFactor);
    //diffusionVector = diffusionVector.ScaleProduct(mi0_in / qi0);
    // Definition centrifugal && coriolis
    Vector3 centrifugal = Vector3(0.0, 0.0, 0.0);
    Vector3 coriolis = Vector3(0.0, 0.0, 0.0);
    Vector3 tempGravity = Vector3(0.0, 0.0, 0.0);
    Vector3 tempOmega = Vector3(omega_wx, omega_wy, omega_wz);
    // *(- mu /m * m / e) = - mu / e , where mu = m * v_perp^2 / 2 / B
    if (mirror_force == 1)
        tempGradB = tempGradB.ScaleProduct(-1.0 * mu / qi0); // grad B: mirror force effects
    // Centrifugal force: need to figure out the accleration of rotation first
    // a = w x ( w x r ), then the centrifugal force (a) is: a = - w x ( w x r )
    // Coriolis force is negligable to centrifugal force
    if (centrifugal_force == 1)
        centrifugal = tempOmega.CrossProduct(tempPos).CrossProduct(tempOmega); // acceleration of centrifugal: (w X r) X w
    // Coriolis force
    if (coriolis_force == 1)
        coriolis = vp.CrossProduct(tempOmega).ScaleProduct(2.0); // acceleration of  coriolis: -2(w X v)
    // new gravity
    if (gravity_force == 1)
        tempGravity = tempPos.NormalizedVector().ScaleProduct(-1.0 * gravity * radius * radius / tempPos.norm2());
    //std::cout << " test force magnitude " << " tempE " << tempe.DotProduct(tempb.NormalizedVector()) << " gravity " << tempGravity.DotProduct(tempb.NormalizedVector()) <<
    //            " centrifugal " << centrifugal.DotProduct(tempb.NormalizedVector()) << " coriolis " << coriolis.DotProduct(tempb.NormalizedVector()) << 
    //            " gradB " << tempGradB.DotProduct(tempb.NormalizedVector()) << " mu force  " << tempGradB.DotProduct(tempb.NormalizedVector())* mu / mi0_in<<  " \n";
    //std::cin.get();
    // put together
    if (coordinate_rotate == 1)
        tempGravity = tempGravity.PlusProduct(centrifugal).PlusProduct(coriolis);
    //
    tempGravity = tempGravity.ScaleProduct(mi0_in / qi0); // Gravity effects into reduced E
    // revised E including gravity and gradB, and centrifugal and coriolis force if initilized particles
    if (electric_force == 0)
        tempe = Vector3(0.0, 0.0, 0.0);
    tempe = tempe.PlusProduct(tempGravity).PlusProduct(tempGradB); // .PlusProduct( diffusionVector);
    // test only Grad B
    // tempe = tempGradB;

    //   tempe = tempe.MinusProduct(tempb.NormalizedVector().ScaleProduct( tempe.DotProduct(tempb.NormalizedVector()) )).PlusProduct(tempGradB);
    //   tempe = tempe.MinusProduct(tempb.NormalizedVector().ScaleProduct( tempe.DotProduct(tempb.NormalizedVector()) ));
    // test particles only under mirror force
    //if( move_type ==0)
    //tempe = tempe.MinusProduct(tempb.NormalizedVector().ScaleProduct( tempe.DotProduct(tempb.NormalizedVector()) )).PlusProduct(tempGradB); // bounce motion
    //    tempe = tempGradB;
    // test gravity
    //    tempe = tempPos;
    // test mirror force effects
    //    tempe = tempGradB;
    //
    // 1.3 Vector3 t ; Vector3 s
    // 
    Vector3 t = tempb.ScaleProduct(qtm);
    Vector3 s = t.ScaleProduct(2.0 / (1.0 + t.norm2()));
    // 2. Boris equations
    // 2.1 equation I: v1 = v + E * qtm
    Vector3 v1 = vp.PlusProduct(tempe.ScaleProduct(qtm));
    // 2.2 equation II: v2 = v1 + v1 X t
    Vector3 v2 = v1.PlusProduct(v1.CrossProduct(t));
    // 2.3 equation III: v3 = v1 + v2 X s
    Vector3 v3 = v1.PlusProduct(v2.CrossProduct(s));
    // 2.4 equation IV: v = v3 + E * qtm
    vp = v3.PlusProduct(tempe.ScaleProduct(qtm));
    // 3. update the postion
    if (maindomain == 0)
    {
        return UpdateUint_64();
    }
    else
    {
        return UpdateUint_64_test();
    }
}
